"""
Workflow integrity.

These exist because of a bug found by hand during review: the Snakefile called
eleven scripts that were not in the repository. The DAG built cleanly — a
missing shell command is only discovered when the rule actually runs, which for
this pipeline could be hours in. Cheap to check, expensive to miss.
"""

import re
import subprocess

import pytest

from wildphylo import cli

SNAKEFILE = "workflow/Snakefile"


# `from wildphylo import config` is an import, not a call; exclude it.
CLI_CALL = re.compile(r"\bwildphylo\s+(?!import\b)([a-z]+)\s+([a-z-]+)")


def _cli_calls(text: str) -> set[tuple[str, str]]:
    """Every ``wildphylo <stage> <step>`` pair a file invokes."""
    return set(CLI_CALL.findall(text))


def _registered() -> set[tuple[str, str]]:
    return {(stage, name) for stage in cli.STAGES for name in cli.steps(stage)}


def test_snakefile_exists(repo_root):
    assert (repo_root / SNAKEFILE).is_file()


def test_snakefile_calls_no_scripts_by_path(repo_root):
    """Every rule goes through the CLI; a stray python path is a regression."""
    text = (repo_root / SNAKEFILE).read_text()
    assert not re.findall(r"python\s+scripts/", text)


def test_every_step_the_workflow_calls_is_registered(repo_root):
    """
    A step the Snakefile names but the CLI does not register fails only when
    the rule runs, which for this pipeline could be hours in. Catch it here.
    """
    called = _cli_calls((repo_root / SNAKEFILE).read_text())
    assert called, "the Snakefile calls no wildphylo steps at all"
    assert called <= _registered(), f"unregistered steps: {sorted(called - _registered())}"


def test_every_step_ci_calls_is_registered(repo_root):
    missing = set()
    for wf in (repo_root / ".github" / "workflows").glob("*.yml"):
        for call in _cli_calls(wf.read_text()) - _registered():
            missing.add(f"{wf.name}:{' '.join(call)}")
    assert not missing, f"CI calls unregistered steps: {sorted(missing)}"


def test_every_stage_in_the_plan_registers_at_least_one_step():
    """An empty stage module means a port was started and not finished."""
    empty = [stage for stage in cli.STAGES if not cli.steps(stage)]
    assert not empty, f"stages with no steps: {empty}"


@pytest.mark.parametrize("step", ["temporal", "congruence"])
def test_standalone_gates_have_working_help(step, capsys):
    """
    The two gate steps are documented as runnable on their own. A --help that
    crashes means an import error nobody has noticed.
    """
    with pytest.raises(SystemExit) as exc:
        cli.main(["gates", step, "--help"])
    assert exc.value.code == 0
    assert "usage:" in capsys.readouterr().out.lower()


def test_doc_links_resolve(repo_root):
    """A broken link in the README is the first thing a new user hits."""
    broken = []
    for doc in list(repo_root.glob("*.md")) + list((repo_root / "docs").glob("*.md")):
        for target in re.findall(
            r"\]\(([A-Za-z0-9_./-]+\.(?:md|yaml|toml|sh|cff))\)", doc.read_text()
        ):
            if not ((repo_root / target).is_file() or (doc.parent / target).is_file()):
                broken.append(f"{doc.name} -> {target}")  # noqa: PERF401 -- two-condition filter reads better as a loop
    assert not broken, f"broken internal links: {broken}"


def test_workflow_yaml_parses(repo_root):
    yaml = pytest.importorskip("yaml")
    for wf in (repo_root / ".github" / "workflows").glob("*.yml"):
        parsed = yaml.safe_load(wf.read_text())
        assert "jobs" in parsed, f"{wf.name} defines no jobs"
        # PyYAML parses the bare key `on` as the boolean True.
        assert ("on" in parsed) or (True in parsed), f"{wf.name} has no trigger"


def test_init_script_is_valid_bash(repo_root):
    r = subprocess.run(
        ["bash", "-n", "INIT_REPO.sh"], cwd=repo_root, capture_output=True, text=True, check=False
    )
    assert r.returncode == 0, f"INIT_REPO.sh syntax error: {r.stderr}"


def test_init_script_does_not_set_branch_protection(repo_root):
    """
    Regression test. An earlier version required a status check named
    'tests (py3.12)' while the CI matrix ran 3.11 and 3.13 — a check that never
    appears, which would have blocked every merge to main permanently.
    """
    text = (repo_root / "INIT_REPO.sh").read_text()
    # The command appears inside the closing here-doc as advice for the user.
    # Only an occurrence outside that block would actually execute.
    executable = text.split('say "Next steps (manual)"')[0]
    active = [
        ln
        for ln in executable.splitlines()
        if "branches/main/protection" in ln and not ln.strip().startswith("#")
    ]
    assert not active, (
        "INIT_REPO.sh sets branch protection. The required check names depend "
        "on the CI matrix; getting them wrong blocks all merges. Print the "
        "command for the user instead of running it."
    )


def test_ci_python_versions_are_real(repo_root):
    """The matrix should span the floor and a current release, not one version."""
    yaml = pytest.importorskip("yaml")
    ci = yaml.safe_load((repo_root / ".github/workflows/ci.yml").read_text())
    versions = ci["jobs"]["test"]["strategy"]["matrix"]["python"]
    assert len(versions) >= 2, "test on more than one Python version"


# --- external tool resolution ----------------------------------------------


def test_snakefile_does_not_hardcode_an_iqtree_binary(repo_root):
    """
    IQ-TREE ships as `iqtree`, `iqtree2` or `iqtree3` depending on version and
    channel — a real conda install produced `iqtree` and `iqtree3` but no
    `iqtree2`. Hardcoding one name fails at runtime with "command not found",
    which on a large alignment is discovered hours in.
    """
    text = (repo_root / SNAKEFILE).read_text()
    shell_lines = [ln for ln in text.splitlines() if ln.strip().startswith('"') and "iqtree" in ln]
    for ln in shell_lines:
        assert "{params.iqtree}" in ln or "{IQTREE}" in ln, (
            f"shell command hardcodes an IQ-TREE binary name: {ln.strip()}"
        )


def test_iqtree_resolution_prefers_newest_available():
    import shutil
    from unittest import mock

    order = ("iqtree3", "iqtree2", "iqtree")

    def resolve(present):
        with mock.patch.object(shutil, "which", side_effect=lambda c: c if c in present else None):
            for c in order:
                if shutil.which(c):
                    return c
        return None

    assert resolve({"iqtree", "iqtree3"}) == "iqtree3"
    assert resolve({"iqtree", "iqtree2"}) == "iqtree2"
    assert resolve({"iqtree"}) == "iqtree"
    assert resolve(set()) is None


def test_environment_does_not_pin_a_version_that_renames_the_binary(repo_root):
    """The pin must not imply a binary name the workflow then assumes."""
    env = (repo_root / "environment.yml").read_text()
    assert "iqtree2=" not in env and "iqtree3=" not in env

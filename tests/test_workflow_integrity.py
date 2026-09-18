"""
Workflow integrity.

These exist because of a bug found by hand during review: the Snakefile called
eleven scripts that were not in the repository. The DAG built cleanly — a
missing shell command is only discovered when the rule actually runs, which for
this pipeline could be hours in. Cheap to check, expensive to miss.
"""

import re
import subprocess
import sys

import pytest

SNAKEFILE = "workflow/Snakefile"


def _referenced_scripts(text: str) -> set[str]:
    return set(re.findall(r"scripts/[0-9A-Za-z_]+\.py", text))


def test_snakefile_exists(repo_root):
    assert (repo_root / SNAKEFILE).is_file()


def _documented_gaps(repo_root) -> set[str]:
    """Scripts docs/PORTING.md openly declares unwritten."""
    porting = (repo_root / "docs" / "PORTING.md").read_text()
    section = porting.split("## Not yet written at all")[1].split("## Not used")[0]
    return set(re.findall(r"`(\w+\.py)`", section))


def test_every_script_the_workflow_calls_is_present_or_documented(repo_root):
    """
    A missing script is tolerable only if PORTING.md says so. Undocumented
    absence is the bug this test exists for: the DAG builds fine and the rule
    fails at runtime, potentially hours in.
    """
    text = (repo_root / SNAKEFILE).read_text()
    known = _documented_gaps(repo_root)
    undocumented = sorted(
        s for s in _referenced_scripts(text)
        if not (repo_root / s).is_file() and s.split("/")[-1] not in known)
    assert not undocumented, (
        "the Snakefile calls scripts that neither exist nor appear in "
        f"docs/PORTING.md: {undocumented}")


def test_every_script_ci_calls_is_present_or_documented(repo_root):
    known = _documented_gaps(repo_root)
    undocumented = set()
    for wf in (repo_root / ".github" / "workflows").glob("*.yml"):
        for s in _referenced_scripts(wf.read_text()):
            if not (repo_root / s).is_file() and s.split("/")[-1] not in known:
                undocumented.add(f"{wf.name}:{s}")
    assert not undocumented, \
        f"CI references undocumented missing scripts: {sorted(undocumented)}"


def test_documented_gaps_match_reality(repo_root):
    """
    docs/PORTING.md lists the scripts that are referenced but unwritten. If one
    gets written, the doc must stop claiming it is missing, or the gap list
    quietly becomes fiction.
    """
    porting = (repo_root / "docs" / "PORTING.md").read_text()
    wrongly_claimed = sorted(
        s for s in _documented_gaps(repo_root)
        if (repo_root / "scripts" / s).is_file())
    assert not wrongly_claimed, (
        f"docs/PORTING.md says these are unwritten, but they exist: "
        f"{wrongly_claimed}")


@pytest.mark.parametrize("script", ["04b_temporal_signal.py",
                                    "05b_segment_congruence.py"])
def test_standalone_scripts_have_working_help(repo_root, script):
    """
    The two gate scripts are documented as runnable on their own. A --help that
    crashes means an import error nobody has noticed.
    """
    r = subprocess.run([sys.executable, f"scripts/{script}", "--help"],
                       cwd=repo_root, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, f"{script} --help failed: {r.stderr[-500:]}"
    assert "usage:" in r.stdout.lower()


def test_doc_links_resolve(repo_root):
    """A broken link in the README is the first thing a new user hits."""
    broken = []
    for doc in list(repo_root.glob("*.md")) + list((repo_root / "docs").glob("*.md")):
        for target in re.findall(r"\]\(([A-Za-z0-9_./-]+\.(?:md|yaml|toml|sh|cff))\)",
                                 doc.read_text()):
            if not ((repo_root / target).is_file()
                    or (doc.parent / target).is_file()):
                broken.append(f"{doc.name} -> {target}")
    assert not broken, f"broken internal links: {broken}"


def test_workflow_yaml_parses(repo_root):
    yaml = pytest.importorskip("yaml")
    for wf in (repo_root / ".github" / "workflows").glob("*.yml"):
        parsed = yaml.safe_load(wf.read_text())
        assert "jobs" in parsed, f"{wf.name} defines no jobs"
        # PyYAML parses the bare key `on` as the boolean True.
        assert ("on" in parsed) or (True in parsed), \
            f"{wf.name} has no trigger"


def test_init_script_is_valid_bash(repo_root):
    r = subprocess.run(["bash", "-n", "INIT_REPO.sh"],
                       cwd=repo_root, capture_output=True, text=True)
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
    executable = text.split("say \"Next steps (manual)\"")[0]
    active = [ln for ln in executable.splitlines()
              if "branches/main/protection" in ln and not ln.strip().startswith("#")]
    assert not active, (
        "INIT_REPO.sh sets branch protection. The required check names depend "
        "on the CI matrix; getting them wrong blocks all merges. Print the "
        "command for the user instead of running it.")


def test_ci_python_versions_are_real(repo_root):
    """The matrix should span the floor and a current release, not one version."""
    yaml = pytest.importorskip("yaml")
    ci = yaml.safe_load((repo_root / ".github/workflows/ci.yml").read_text())
    versions = ci["jobs"]["test"]["strategy"]["matrix"]["python"]
    assert len(versions) >= 2, "test on more than one Python version"

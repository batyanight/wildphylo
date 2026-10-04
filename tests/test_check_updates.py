"""
Tests for the scheduled-rebuild decision.

Entrez is never contacted: the counting function is stubbed. A test suite that
needed the network would be skipped in CI and would then be verifying nothing
on exactly the code path CI runs.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location(
        "chk13", ROOT / "scripts" / "13_check_updates.py")
    m = importlib.util.module_from_spec(spec)
    # Register before exec: @dataclass resolves cls.__module__ through
    # sys.modules while building __init__, and a path-loaded module is not
    # there unless it is put there.
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


@pytest.fixture
def cfg():
    return {"id": "demo",
            "nextstrain": {"update": {"min_new_accessions": 5,
                                      "max_new_accessions": 500}}}


# --- the decision table -----------------------------------------------------

def test_no_previous_build_proceeds(mod, cfg):
    d = mod.decide(cfg, "demo", current=100, previous=None, prev_build=None)
    assert d.proceed and d.delta is None


def test_no_change_skips(mod, cfg):
    d = mod.decide(cfg, "demo", 100, 100, "20260101")
    assert not d.proceed and d.delta == 0


def test_below_the_minimum_skips(mod, cfg):
    d = mod.decide(cfg, "demo", 103, 100, "20260101")
    assert not d.proceed
    assert "min_new_accessions" in d.reason


def test_within_range_proceeds(mod, cfg):
    d = mod.decide(cfg, "demo", 110, 100, "20260101")
    assert d.proceed and d.delta == 10 and d.level == "note"


def test_above_the_maximum_proceeds_but_is_flagged(mod, cfg):
    """
    A bulk submission is when a rebuild is most worth doing AND most worth
    looking at. Skipping it would be the wrong response to the most
    interesting case.
    """
    d = mod.decide(cfg, "demo", 900, 100, "20260101")
    assert d.proceed and d.level == "warn"


def test_a_shrinking_dataset_does_not_rebuild(mod, cfg):
    """
    GenBank withdraws records very rarely. Rebuilding here would replace a
    larger dataset with a smaller one and the Auspice build would silently
    lose sequences.
    """
    d = mod.decide(cfg, "demo", 90, 100, "20260101")
    assert not d.proceed and d.level == "error"
    assert "LOST" in d.reason


def test_force_overrides_everything_including_a_shrinking_dataset(mod, cfg):
    assert mod.decide(cfg, "demo", 100, 100, "20260101", force=True).proceed
    d = mod.decide(cfg, "demo", 90, 100, "20260101", force=True)
    assert d.proceed and d.level == "error", "forced, but still reported as an error"


# --- reading the previous build ---------------------------------------------

def test_previous_build_reads_the_newest_dated_directory(mod, tmp_path):
    for date, n in (("20260101", 3), ("20260301", 7), ("20260201", 5)):
        acc = tmp_path / "demo" / date / "raw" / "demo.acc"
        acc.parent.mkdir(parents=True)
        acc.write_text("\n".join(f"ACC{i}" for i in range(n)) + "\n")
    name, count = mod.previous_build(tmp_path, "demo")
    assert (name, count) == ("20260301", 7)


def test_a_build_with_no_accession_file_is_skipped_for_the_one_below_it(mod, tmp_path):
    """
    An interrupted build leaves a directory with no .acc. Treating that as zero
    would report every existing record as new.
    """
    good = tmp_path / "demo" / "20260101" / "raw" / "demo.acc"
    good.parent.mkdir(parents=True)
    good.write_text("A\nB\nC\n")
    (tmp_path / "demo" / "20260301" / "raw").mkdir(parents=True)
    name, count = mod.previous_build(tmp_path, "demo")
    assert (name, count) == ("20260101", 3)


def test_missing_builds_directory_is_not_an_error(mod, tmp_path):
    assert mod.previous_build(tmp_path, "nothing-here") == (None, None)


def test_blank_lines_in_the_accession_file_are_not_counted(mod, tmp_path):
    acc = tmp_path / "demo" / "20260101" / "raw" / "demo.acc"
    acc.parent.mkdir(parents=True)
    acc.write_text("A\n\nB\n\n\n")
    assert mod.previous_build(tmp_path, "demo")[1] == 2


# --- config discovery -------------------------------------------------------

def test_template_config_is_skipped(mod, repo_root):
    got = mod.load_configs(repo_root / "config" / "pathogen", None)
    assert all(not p.name.startswith("_") for p, _ in got)
    assert len(got) >= 2


def test_only_restricts_to_one_pathogen(mod, repo_root):
    got = mod.load_configs(repo_root / "config" / "pathogen", "cdv-1200")
    assert [c["id"] for _, c in got] == ["cdv-1200"]


def test_an_unknown_only_is_an_error_not_an_empty_rebuild(mod, repo_root):
    with pytest.raises(SystemExit):
        mod.load_configs(repo_root / "config" / "pathogen", "does-not-exist")


def test_an_invalid_config_stops_the_check_rather_than_being_skipped(mod, tmp_path):
    """
    Skipping it would make the scheduled rebuild silently narrow: the pathogen
    just stops being rebuilt and nothing says so.
    """
    d = tmp_path / "pathogen"
    d.mkdir()
    (d / "broken.yaml").write_text(yaml.safe_dump({"id": "broken"}))
    with pytest.raises(SystemExit):
        mod.load_configs(d, None)


def test_the_query_comes_from_01_not_a_copy(mod, repo_root):
    """
    If the count used a reimplemented query, a drift between the two would show
    up as a delta that looks like a real change in GenBank.
    """
    fetch = mod._fetch_module()
    assert hasattr(fetch, "build_query")
    cfg = yaml.safe_load(
        (repo_root / "config" / "pathogen" / "cdv-1200.yaml").read_text())
    q = fetch.build_query(cfg)
    assert f"txid{cfg['fetch']['taxid']}" in q


# --- CLI, with Entrez stubbed ----------------------------------------------

def _run_cli(repo_root, tmp_path, count, extra_env=None):
    """
    Run the CLI with count_records monkeypatched at the source level, via a
    tiny driver script. Done in a subprocess so the $GITHUB_OUTPUT handling and
    the exit code are exercised as CI exercises them.
    """
    driver = tmp_path / "driver.py"
    driver.write_text(f"""
import importlib.util, sys
from pathlib import Path
ROOT = Path({str(repo_root)!r})
spec = importlib.util.spec_from_file_location("chk13", ROOT / "scripts" / "13_check_updates.py")
m = importlib.util.module_from_spec(spec)
sys.modules["chk13"] = m
spec.loader.exec_module(m)
m.count_records = lambda q: {count}
sys.argv = ["13_check_updates.py", "--config-dir", str(ROOT / "config" / "pathogen"),
            "--builds-dir", {str(tmp_path / 'builds')!r},
            "--only", "cdv-1200", "--github-output",
            "--out-json", {str(tmp_path / 'out.json')!r}]
sys.exit(m.main())
""")
    env = dict(os.environ)
    env["GITHUB_OUTPUT"] = str(tmp_path / "gh_output")
    env.update(extra_env or {})
    (tmp_path / "gh_output").write_text("")
    r = subprocess.run([sys.executable, str(driver)], cwd=repo_root,
                       capture_output=True, text=True, timeout=120, env=env)
    return r, (tmp_path / "gh_output").read_text(), tmp_path / "out.json"


def _stage_previous(tmp_path, n, pathogen="cdv-1200", date="20260101"):
    acc = tmp_path / "builds" / pathogen / date / "raw" / f"{pathogen}.acc"
    acc.parent.mkdir(parents=True, exist_ok=True)
    acc.write_text("\n".join(f"ACC{i}" for i in range(n)) + "\n")


def test_cli_writes_github_outputs_when_there_is_work(repo_root, tmp_path):
    _stage_previous(tmp_path, 100)
    r, gh, out = _run_cli(repo_root, tmp_path, count=200)
    assert r.returncode == 0, r.stderr[-1500:]
    assert "proceed=true" in gh
    assert json.loads(gh.split("pathogens=")[1].strip()) == ["cdv-1200"]
    assert json.loads(out.read_text())["decisions"][0]["delta"] == 100


def test_cli_reports_false_when_nothing_changed(repo_root, tmp_path):
    _stage_previous(tmp_path, 200)
    r, gh, out = _run_cli(repo_root, tmp_path, count=200)
    assert r.returncode == 0
    assert "proceed=false" in gh
    assert json.loads(gh.split("pathogens=")[1].strip()) == []


def test_cli_force_env_overrides_no_change(repo_root, tmp_path):
    _stage_previous(tmp_path, 200)
    r, gh, _ = _run_cli(repo_root, tmp_path, count=200, extra_env={"FORCE": "true"})
    assert "proceed=true" in gh


def test_cli_help_works():
    r = subprocess.run([sys.executable, "scripts/13_check_updates.py", "--help"],
                       cwd=ROOT, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0
    assert "usage:" in r.stdout.lower()

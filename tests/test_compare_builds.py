"""
Tests for the build-to-build comparison gate.

The framing being tested is that a large change between builds is a CURATION
REGRESSION until shown otherwise. New data arriving monthly moves a TMRCA by
months; a dataset that shrinks has almost always lost records to a changed
query. A gate that treated these as findings would publish the pipeline's own
bugs as results.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from wildphylo import cli
from wildphylo import config as cfgmod
from wildphylo.layout import BuildLayout
from wildphylo.stages import gates

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def mod():
    return gates


def make_build(root, date, accs=None, temporal=None, tmrca=None, gates=None):
    b = root / date
    if accs is not None:
        (b / "raw").mkdir(parents=True, exist_ok=True)
        (b / "raw" / "x.acc").write_text("\n".join(accs) + "\n")
    if temporal is not None:
        (b / "temporal").mkdir(parents=True, exist_ok=True)
        (b / "temporal" / "H_temporal_signal.json").write_text(json.dumps(temporal))
    if tmrca is not None:
        (b / "auspice").mkdir(parents=True, exist_ok=True)
        (b / "auspice" / "x.json").write_text(
            json.dumps({"tree": {"name": "ROOT", "node_attrs": {"num_date": {"value": tmrca}}}})
        )
    for rel, doc in (gates or {}).items():
        p = b / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(doc))
    b.mkdir(parents=True, exist_ok=True)
    return b


def auspice(tmrca):
    return {"tree": {"name": "ROOT", "node_attrs": {"num_date": {"value": tmrca}}}}


# --- the quiet-success case -------------------------------------------------


def test_a_normal_monthly_rebuild_is_clean(mod, tmp_path):
    prev = make_build(
        tmp_path,
        "20260901",
        accs=[f"A{i}" for i in range(100)],
        temporal={"verdict": "pass", "r_squared": 0.42, "slope": 7.1e-4},
    )
    curr = make_build(
        tmp_path,
        "20261001",
        accs=[f"A{i}" for i in range(108)],
        temporal={"verdict": "pass", "r_squared": 0.43, "slope": 7.2e-4},
    )
    c = mod.compare_builds(
        prev, curr, "x", current_auspice=auspice(1965.0), previous_auspice=auspice(1964.6)
    )
    assert c.verdict == "ok"


# --- the regressions --------------------------------------------------------


def test_losing_sequences_is_an_error_and_names_them(mod, tmp_path):
    prev = make_build(tmp_path, "20260901", accs=[f"A{i}" for i in range(100)])
    curr = make_build(tmp_path, "20261001", accs=[f"A{i}" for i in range(80)])
    c = mod.compare_builds(prev, curr, "x")
    assert c.verdict == "regression"
    msg = " ".join(f.message for f in c.findings if f.level == "error")
    assert "A80" in msg, "the lost accessions must be named, not just counted"


def test_sequence_loss_can_be_downgraded(mod, tmp_path):
    prev = make_build(tmp_path, "20260901", accs=["A", "B", "C"])
    curr = make_build(tmp_path, "20261001", accs=["A", "B"])
    c = mod.compare_builds(prev, curr, "x", alert_on_loss=False)
    assert c.verdict == "review"


def test_a_stable_count_hiding_a_composition_change_is_caught(mod, tmp_path):
    """
    The failure a count-only check cannot see: twenty records swapped out for
    twenty others. That is what a broken query produces.
    """
    prev = make_build(tmp_path, "20260901", accs=[f"A{i}" for i in range(100)])
    curr = make_build(
        tmp_path, "20261001", accs=[f"A{i}" for i in range(80)] + [f"B{i}" for i in range(20)]
    )
    c = mod.compare_builds(prev, curr, "x", alert_on_loss=False)
    msgs = " ".join(f.message for f in c.findings)
    assert "swapped out" in msgs
    assert "count is unchanged" in msgs


def test_a_big_tmrca_shift_is_an_error(mod, tmp_path):
    prev = make_build(tmp_path, "20260901", accs=["A"])
    curr = make_build(tmp_path, "20261001", accs=["A"])
    c = mod.compare_builds(
        prev,
        curr,
        "x",
        max_tmrca_shift=5.0,
        previous_auspice=auspice(1965.0),
        current_auspice=auspice(1925.0),
    )
    assert c.verdict == "regression"
    assert any("TMRCA" in f.message and f.level == "error" for f in c.findings)


def test_a_small_tmrca_shift_is_only_a_note(mod, tmp_path):
    prev = make_build(tmp_path, "20260901", accs=["A"])
    curr = make_build(tmp_path, "20261001", accs=["A"])
    c = mod.compare_builds(
        prev,
        curr,
        "x",
        max_tmrca_shift=5.0,
        previous_auspice=auspice(1965.0),
        current_auspice=auspice(1963.2),
    )
    assert c.verdict == "ok"


def test_a_gate_verdict_dropping_is_a_regression(mod, tmp_path):
    prev = make_build(
        tmp_path,
        "20260901",
        accs=["A"],
        temporal={"verdict": "pass", "r_squared": 0.4, "slope": 1e-4},
    )
    curr = make_build(
        tmp_path,
        "20261001",
        accs=["A"],
        temporal={"verdict": "fail", "r_squared": 0.4, "slope": 1e-4},
    )
    c = mod.compare_builds(prev, curr, "x")
    assert c.verdict == "regression"
    assert any("regressed" in f.message for f in c.findings)


def test_a_gate_verdict_improving_is_not_a_regression(mod, tmp_path):
    prev = make_build(
        tmp_path,
        "20260901",
        accs=["A"],
        temporal={"verdict": "weak", "r_squared": 0.4, "slope": 1e-4},
    )
    curr = make_build(
        tmp_path,
        "20261001",
        accs=["A"],
        temporal={"verdict": "pass", "r_squared": 0.4, "slope": 1e-4},
    )
    c = mod.compare_builds(prev, curr, "x")
    assert c.verdict == "ok"
    assert any("improved" in f.message for f in c.findings)


def test_per_clade_gates_are_compared_separately(mod, tmp_path):
    """
    One clade's trait-signal failing must not be hidden by another's passing.
    """
    g_prev = {
        "clades/america2/trait_signal.json": {"verdict": "pass"},
        "clades/asia1/trait_signal.json": {"verdict": "pass"},
    }
    g_curr = {
        "clades/america2/trait_signal.json": {"verdict": "fail"},
        "clades/asia1/trait_signal.json": {"verdict": "pass"},
    }
    prev = make_build(tmp_path, "20260901", accs=["A"], gates=g_prev)
    curr = make_build(tmp_path, "20261001", accs=["A"], gates=g_curr)
    c = mod.compare_builds(prev, curr, "x")
    assert c.verdict == "regression"
    assert any("america2" in f.message for f in c.findings)


def test_relative_change_keeps_its_sign(mod, tmp_path):
    """
    Regression. R^2 falling and R^2 rising are not the same event, and a
    magnitude printed with a leading + reads as the opposite of what happened.
    """
    prev = make_build(
        tmp_path,
        "20260901",
        accs=["A"],
        temporal={"verdict": "pass", "r_squared": 0.42, "slope": 1e-4},
    )
    curr = make_build(
        tmp_path,
        "20261001",
        accs=["A"],
        temporal={"verdict": "pass", "r_squared": 0.11, "slope": 1e-4},
    )
    c = mod.compare_builds(prev, curr, "x")
    r2 = next(f.message for f in c.findings if "R^2" in f.message)
    assert "-74%" in r2


# --- the cases that must NOT be errors ---------------------------------------


def test_a_first_build_is_not_a_failure(tmp_path):
    """A new pathogen has nothing to compare against. That must not block it."""
    make_build(tmp_path / "p", "20261001", accs=["A"])
    out = tmp_path / "r.md"
    rc = cli.main(["gates", "compare", "--builds-dir", str(tmp_path / "p"), "--out", str(out)])
    assert rc == 0
    assert "nothing to compare" in out.read_text()


def test_a_missing_builds_directory_is_not_a_failure(tmp_path):
    rc = cli.main(["gates", "compare", "--builds-dir", str(tmp_path / "nope")])
    assert rc == 0


def test_a_tree_with_no_date_is_not_compared_as_zero(mod):
    """
    A tree that was not time-scaled has no TMRCA. Reading that as 0.0 would
    produce a two-thousand-year shift out of nothing.
    """
    assert mod.auspice_tmrca({"tree": {"name": "ROOT", "node_attrs": {}}}) is None
    assert mod.auspice_tmrca({}) is None
    assert mod.auspice_tmrca(None) is None


def test_auspice_tree_as_a_list_is_handled(mod):
    """Some Augur versions wrap the root in a single-element list."""
    doc = {"tree": [{"name": "ROOT", "node_attrs": {"num_date": {"value": 1965.0}}}]}
    assert mod.auspice_tmrca(doc) == pytest.approx(1965.0)


# --- both call signatures ---------------------------------------------------


def test_the_ci_invocation_works_without_an_auspice_json(tmp_path):
    """
    CI stops at the temporal gate, so there is no Auspice build to pass. The
    script must still compare what exists.
    """
    root = tmp_path / "cdv"
    make_build(
        root,
        "20260901",
        accs=["A", "B", "C"],
        temporal={"verdict": "pass", "r_squared": 0.4, "slope": 1e-4},
    )
    make_build(
        root,
        "20261001",
        accs=["A", "B"],
        temporal={"verdict": "pass", "r_squared": 0.4, "slope": 1e-4},
    )
    out = tmp_path / "r.md"
    rc = cli.main(
        ["gates", "compare", "--pathogen", "cdv", "--builds-dir", str(root), "--out", str(out)]
    )
    assert rc == 1, "a lost sequence should fail the step"
    assert "cdv" in out.read_text()


def test_the_snakemake_invocation_works_with_one(tmp_path):
    root = tmp_path / "cdv"
    make_build(root, "20260901", accs=["A"], tmrca=1965.0)
    make_build(root, "20261001", accs=["A"])
    cur = tmp_path / "current.json"
    cur.write_text(json.dumps(auspice(1964.8)))
    out = tmp_path / "r.md"
    rc = cli.main(
        [
            "gates",
            "compare",
            "--current",
            str(cur),
            "--builds-dir",
            str(root),
            "--max-tmrca-shift",
            "5.0",
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    assert "TMRCA" in out.read_text()


def test_cli_help_works(capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main(["gates", "compare", "--help"])
    rc = exc.value.code
    assert rc == 0
    assert "usage:" in capsys.readouterr().out.lower()


# --- --build-dir ------------------------------------------------------------

CDV_CONFIG = ROOT / "tests" / "fixtures" / "configs" / "cdv.yaml"


def _two_builds(tmp_path):
    """Two dated builds of the cdv fixture config; the newer one has its Auspice dataset."""
    cfg = cfgmod.load(CDV_CONFIG)
    root = tmp_path / "builds" / cfg["id"]
    make_build(root, "20260901", accs=["A"], tmrca=1965.0)
    make_build(root, "20261001", accs=["A"])
    layout = BuildLayout.at(root / "20261001", cfg)
    layout.auspice_dataset(cfg["nextstrain"]["build_name"]).parent.mkdir(parents=True)
    layout.auspice_dataset(cfg["nextstrain"]["build_name"]).write_text(json.dumps(auspice(1964.8)))
    return layout


def test_build_dir_compares_against_the_sibling_builds_and_writes_into_the_build(tmp_path):
    layout = _two_builds(tmp_path)
    rc = cli.main(
        ["gates", "compare", "--config", str(CDV_CONFIG), "--build-dir", str(layout.root)]
    )
    assert rc == 0
    assert "TMRCA" in layout.build_comparison().read_text(), "the build's own dataset was used"
    assert json.loads(layout.build_comparison_json().read_text())


def test_build_dir_without_an_auspice_dataset_still_compares(tmp_path):
    layout = _two_builds(tmp_path)
    for path in layout.root.joinpath("auspice").iterdir():
        path.unlink()
    rc = cli.main(
        ["gates", "compare", "--config", str(CDV_CONFIG), "--build-dir", str(layout.root)]
    )
    assert rc == 0
    assert layout.build_comparison().is_file()


def test_an_explicit_flag_beats_the_build_dir_default(tmp_path):
    layout = _two_builds(tmp_path)
    elsewhere = tmp_path / "elsewhere.md"
    argv = ["gates", "compare", "--config", str(CDV_CONFIG), "--build-dir", str(layout.root)]
    assert cli.main([*argv, "--out", str(elsewhere)]) == 0
    assert elsewhere.is_file()
    assert not layout.build_comparison().exists()
    assert layout.build_comparison_json().is_file()


def test_the_sequence_loss_policy_follows_the_config_unless_a_flag_overrides_it(tmp_path):
    cfg = yaml.safe_load(CDV_CONFIG.read_text())
    cfg["nextstrain"]["update"]["alert_on_sequence_loss"] = False
    config = tmp_path / "cdv.yaml"
    config.write_text(yaml.safe_dump(cfg))
    root = tmp_path / "builds" / cfg["id"]
    make_build(root, "20260901", accs=["A", "B", "C"])
    make_build(root, "20261001", accs=["A"])
    argv = ["gates", "compare", "--config", str(config), "--build-dir", str(root / "20261001")]
    assert cli.main(argv) == 0, "alert_on_sequence_loss: false in the config"
    assert cli.main([*argv, "--alert-on-sequence-loss"]) == 1

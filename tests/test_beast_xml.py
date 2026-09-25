"""
Tests for 07_make_beast_xml.py and lib/beastxml.py.

The original baked one model into a single f-string — HKY+G4, relaxed lognormal
clock, constant-size coalescent — and defaulted the clock rate to 7.46e-4, a CDV
H-gene value, whatever pathogen was running. Every config declaration was
ignored. These tests exist to keep that from coming back.
"""

import json
import random
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from lib.beastxml import (  # noqa: E402
    resolve_clock_rate, tree_prior_blocks, clock_blocks, tip_date_blocks,
)

MEASURED = 5.12e-4


@pytest.fixture(scope="module")
def aln(tmp_path_factory):
    d = tmp_path_factory.mktemp("beast")
    rng = random.Random(5)
    rows = []
    for i in range(30):
        date = 1995 + rng.random() * 28
        host = "procyonid" if i % 3 else "wild_canid"
        seq = "".join(rng.choice("ACGT") for _ in range(600))
        rows.append(f">ACC{i:05d}|{host}|{date:.3f}\n{seq}\n")
    p = d / "aln.fasta"
    p.write_text("".join(rows))
    return p


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    p = tmp_path_factory.mktemp("rep") / "ts.json"
    p.write_text(json.dumps({"regressions": [
        {"slope": MEASURED, "r_squared": 0.8, "x_intercept": 1976.0}]}))
    return p


def cfg_variant(tmp_path, **beast):
    c = yaml.safe_load((ROOT / "config/pathogen/cdv.yaml").read_text())
    c["beast"].update(beast)
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(c, sort_keys=False))
    return p


def run(aln, out, *extra):
    return subprocess.run(
        [sys.executable, "scripts/07_make_beast_xml.py", "--aln", str(aln),
         "--out-xml", str(out), "--out-traits", str(out.with_suffix(".traits.txt")),
         "--chain", "100000", *map(str, extra)],
        capture_output=True, text=True, cwd=ROOT)


# --- clock rate provenance --------------------------------------------------

def test_config_rate_wins_when_declared():
    cfg = {"beast": {"clock_rate_prior": 7.46e-4,
                     "clock_rate_source": "Panzera et al."}}
    rate, note = resolve_clock_rate(cfg, MEASURED)
    assert rate == pytest.approx(7.46e-4)
    assert "Panzera" in note


def test_measured_slope_used_when_config_has_none():
    """
    The default that matters for a new pathogen: centre the prior on THIS
    dataset rather than inheriting a rate from another virus.
    """
    rate, note = resolve_clock_rate({"beast": {"clock_rate_prior": None}}, MEASURED)
    assert rate == pytest.approx(MEASURED)
    assert "measured" in note


def test_per_locus_rate_overrides_the_pathogen_rate():
    cfg = {"beast": {"clock_rate_prior": 7.46e-4, "clock_rate_source": "x"}}
    locus = {"clock_rate_prior": 4.26e-4, "clock_rate_source": "Seg-10"}
    rate, note = resolve_clock_rate(cfg, None, locus)
    assert rate == pytest.approx(4.26e-4)
    assert "Seg-10" in note


def test_no_rate_anywhere_is_an_error():
    """Silently inventing a rate is the failure mode this prevents."""
    with pytest.raises(ValueError, match="no clock rate"):
        resolve_clock_rate({"beast": {"clock_rate_prior": None}}, None)


def test_time_scaled_slope_is_rejected(aln, tmp_path):
    """
    A slope of ~1.0 comes from a time-scaled tree and is a units check, not a
    rate. Using it would centre the clock on 1.0 subs/site/year.
    """
    rep = tmp_path / "t.json"
    rep.write_text(json.dumps({"regressions": [{"slope": 1.002}]}))
    cfg = cfg_variant(tmp_path, clock_rate_prior=None, clock_rate_source=None)
    r = run(aln, tmp_path / "x.xml", "--config", cfg, "--temporal-report", rep)
    assert "not a substitution rate" in r.stdout
    assert r.returncode == 1          # nothing left to centre the prior on


# --- the safety gate --------------------------------------------------------

def test_prior_disagreeing_with_data_refuses_to_write(aln, report, tmp_path):
    """
    A prior 49x off the measured slope is what a rate copied from the wrong
    pathogen looks like. Writing the XML anyway would burn days of compute on
    an answer determined by the prior.
    """
    cfg = cfg_variant(tmp_path, clock_rate_prior=2.5e-2)
    out = tmp_path / "bad.xml"
    r = run(aln, out, "--config", cfg, "--temporal-report", report)
    assert r.returncode == 1
    assert not out.exists(), "an XML was written despite the gate failing"


def test_agreeing_prior_passes(aln, report, tmp_path):
    cfg = cfg_variant(tmp_path, clock_rate_prior=6.0e-4)
    out = tmp_path / "ok.xml"
    r = run(aln, out, "--config", cfg, "--temporal-report", report)
    assert r.returncode == 0, r.stderr[-800:]
    assert out.is_file()


# --- model blocks actually change the XML -----------------------------------

@pytest.mark.parametrize("prior,marker", [
    ("coalescent_constant", "CoalescentConstant"),
    ("coalescent_skyline", "BayesianSkyline"),
    ("birth_death", "BirthDeath"),
])
def test_tree_prior_is_honoured(aln, tmp_path, prior, marker):
    """
    A constant-size coalescent assumes a stable, randomly sampled population.
    The CDV clade holding the 1994 Serengeti epidemic has 29 of 82 tips from
    one outbreak year — the wrong model, and the config has to be able to say so.
    """
    cfg = cfg_variant(tmp_path, tree_prior=prior)
    out = tmp_path / f"{prior}.xml"
    r = run(aln, out, "--config", cfg)
    assert r.returncode == 0, r.stderr[-800:]
    text = out.read_text()
    assert marker in text
    ET.fromstring(text)


@pytest.mark.parametrize("model,marker", [
    ("relaxed_lognormal", "UCRelaxedClockModel"),
    ("strict", "StrictClockModel"),
])
def test_clock_model_is_honoured(aln, tmp_path, model, marker):
    cfg = cfg_variant(tmp_path, clock_model=model)
    out = tmp_path / f"{model}.xml"
    assert run(aln, out, "--config", cfg).returncode == 0
    assert marker in out.read_text()


def test_every_tree_prior_produces_well_formed_xml(aln, tmp_path):
    for prior in ("coalescent_constant", "coalescent_skyline", "birth_death"):
        cfg = cfg_variant(tmp_path, tree_prior=prior)
        out = tmp_path / f"wf_{prior}.xml"
        assert run(aln, out, "--config", cfg).returncode == 0
        ET.parse(out)


def test_unknown_tree_prior_is_rejected():
    with pytest.raises(ValueError, match="unknown tree prior"):
        tree_prior_blocks("magic", "p", 30)


def test_clock_block_reports_which_model_to_reference():
    assert clock_blocks("strict", "p", 1e-3, 30)["ref"].startswith("StrictClock")
    assert clock_blocks("relaxed_lognormal", "p", 1e-3, 30)["ref"].startswith("RelaxedClock")


def test_relaxed_clock_rate_categories_scale_with_tips():
    """One category per branch: 2n-2 for an unrooted tree of n tips."""
    assert 'dimension="58"' in clock_blocks("relaxed_lognormal", "p", 1e-3, 30)["state"]


def test_skyline_groups_scale_down_for_small_datasets():
    """5 groups on 30 tips over-parameterises; the block has to adapt."""
    small = tree_prior_blocks("coalescent_skyline", "p", 30)["state"]
    large = tree_prior_blocks("coalescent_skyline", "p", 400)["state"]
    assert 'dimension="2"' in small
    assert 'dimension="5"' in large


# --- sampled tip dates ------------------------------------------------------

DATES = {"A|x|2001.5": 2001.5, "B|x|2010.5": 2010.5, "C|x|2020.0": 2020.0}


def test_interval_policy_adds_a_tip_date_operator():
    blocks = tip_date_blocks("interval", "p", {"A|x|2001.5": 0.5}, DATES)
    assert "TipDatesRandomWalker" in blocks["operators"]
    assert blocks["n_sampled"] == 1


def test_midpoint_policy_adds_nothing():
    blocks = tip_date_blocks("midpoint", "p", {"A|x|2001.5": 0.5}, DATES)
    assert blocks["operators"] == ""


def test_each_sampled_tip_is_bounded_by_its_own_uncertainty():
    """
    THE regression. The first version emitted a bare TipDatesRandomWalker with
    a fixed window and no bounds, so a tip's date could wander without limit —
    worse than a midpoint, which is at least wrong in a known, bounded way.
    BEAST gave it away by tuning toward a 24-year proposal window for tips
    whose real uncertainty is half a year.
    """
    blocks = tip_date_blocks("interval", "p", {"A|x|2001.5": 0.5}, DATES)
    assert "MRCAPrior" in blocks["priors"]
    import re
    lo, hi = map(float, re.search(r'lower="([\d.]+)" upper="([\d.]+)"',
                                  blocks["priors"]).groups())
    # height measured back from the most recent sample (2020.0)
    assert hi - lo == pytest.approx(1.0, abs=1e-4)      # +/- 0.5 yr
    assert lo == pytest.approx(2020.0 - 2002.0, abs=1e-4)


def test_bounds_scale_with_the_uncertainty():
    tight = tip_date_blocks("interval", "p", {"A|x|2001.5": 0.04}, DATES)
    wide = tip_date_blocks("interval", "p", {"A|x|2001.5": 0.5}, DATES)
    import re
    def span(b):
        lo, hi = map(float, re.search(r'lower="([\d.]+)" upper="([\d.]+)"',
                                      b["priors"]).groups())
        return hi - lo
    assert span(tight) < span(wide)
    assert span(tight) == pytest.approx(0.08, abs=1e-4)


def test_no_tip_height_can_go_negative():
    """The most recent tip cannot sit above itself."""
    blocks = tip_date_blocks("interval", "p", {"C|x|2020.0": 0.5}, DATES)
    import re
    lo = float(re.search(r'lower="([\d.]+)"', blocks["priors"]).group(1))
    assert lo >= 0.0


def test_proposal_window_is_not_left_to_the_tuner():
    """
    A window far larger than the bound just wastes steps being rejected at the
    boundary. It should track the uncertainty.
    """
    blocks = tip_date_blocks("interval", "p",
                             {"A|x|2001.5": 0.5, "B|x|2010.5": 0.5}, DATES)
    import re
    w = float(re.search(r'windowSize="([\d.]+)"', blocks["operators"]).group(1))
    assert w <= 0.5


def test_tip_dates_are_logged_so_convergence_can_be_checked():
    """
    Sampled tip dates mix slowly. If they are not in the trace, nobody can
    check their ESS and the warning in check_convergence has nothing to read.
    """
    blocks = tip_date_blocks("interval", "p", {"A|x|2001.5": 0.5}, DATES)
    assert "tipdate." in blocks["log"]


def test_no_dates_supplied_means_no_sampling():
    """Without dates the bounds cannot be computed; silently sampling unbounded
    is the bug this replaced."""
    blocks = tip_date_blocks("interval", "p", {"A|x|2001.5": 0.5}, None)
    assert blocks["n_sampled"] == 0
    assert blocks["operators"] == ""


def test_interval_without_uncertainty_data_is_reported(aln, tmp_path):
    """
    Silently falling back to midpoints would report an HPD narrower than the
    data support, which is the whole thing DR-001 exists to prevent.
    """
    cfg = cfg_variant(tmp_path)          # cdv.yaml already sets interval
    r = run(aln, tmp_path / "i.xml", "--config", cfg)
    assert "no per-tip uncertainty" in r.stdout


# --- date randomisation -----------------------------------------------------

def test_single_replicate_mode_writes_one_file(aln, tmp_path):
    out = tmp_path / "drt3.xml"
    cfg = cfg_variant(tmp_path)
    r = run(aln, out, "--config", cfg, "--replicate", 3)
    assert r.returncode == 0, r.stderr[-800:]
    text = out.read_text()
    assert "DATE-RANDOMISED REPLICATE 3" in text
    ET.fromstring(text)


def test_replicates_differ_from_each_other(aln, tmp_path):
    cfg = cfg_variant(tmp_path)
    run(aln, tmp_path / "r1.xml", "--config", cfg, "--replicate", 1)
    run(aln, tmp_path / "r2.xml", "--config", cfg, "--replicate", 2)
    assert (tmp_path / "r1.xml").read_text() != (tmp_path / "r2.xml").read_text()


def test_randomised_replicate_does_not_sample_tip_dates(aln, tmp_path):
    """
    The DRT asks whether the REAL dates carry signal. Sampling shuffled dates
    within their windows would blur the very thing being tested.
    """
    cfg = cfg_variant(tmp_path)
    out = tmp_path / "drt1.xml"
    run(aln, out, "--config", cfg, "--replicate", 1)
    assert "TipDatesRandomWalker" not in out.read_text()


# --- provenance in the output -----------------------------------------------

def test_xml_header_records_where_the_rate_came_from(aln, report, tmp_path):
    """
    An uncited rate in a generated file is how an assumption becomes a result.
    The header has to say which it was.
    """
    cfg = cfg_variant(tmp_path, clock_rate_prior=None, clock_rate_source=None)
    out = tmp_path / "prov.xml"
    assert run(aln, out, "--config", cfg,
               "--temporal-report", report).returncode == 0
    header = out.read_text()[:1200]
    assert "rate provenance" in header
    assert "measured root-to-tip" in header

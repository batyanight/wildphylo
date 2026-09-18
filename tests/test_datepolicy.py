"""
Tests for the midpoint-vs-interval decision rule.

The rule must give different answers on different data. A recommender that
always says 'interval' is just the literature default wearing a lab coat.
"""

import math
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lib.datepolicy import recommend, rayleigh  # noqa: E402


def dataset(n_year=100, n_precise=40, season_centre=None, season_conc=0.0, seed=5):
    """Build (decimal_years, precisions). season_conc 0 = uniform within year."""
    rng = random.Random(seed)
    years, precs = [], []
    for _ in range(n_year):
        years.append(1990 + rng.randint(0, 30) + 0.5)
        precs.append("year")
    for _ in range(n_precise):
        if season_centre is None or season_conc == 0:
            w = rng.random()
        else:
            w = (season_centre + rng.gauss(0, max(1e-6, (1 - season_conc) / 4))) % 1
        years.append(1990 + rng.randint(0, 30) + w)
        precs.append("day")
    return years, precs


# --- the circular test ------------------------------------------------------

def test_rayleigh_uniform_is_not_significant():
    rng = random.Random(1)
    R, _Z, p = rayleigh([rng.random() for _ in range(200)])
    assert R < 0.2 and p > 0.05


def test_rayleigh_detects_tight_season():
    rng = random.Random(1)
    R, _Z, p = rayleigh([(0.8 + rng.gauss(0, 0.04)) % 1 for _ in range(60)])
    assert R > 0.8 and p < 0.001


def test_rayleigh_too_few_points_is_not_significant():
    """Refuse to call seasonality off five samples."""
    _R, _Z, p = rayleigh([0.8, 0.81, 0.82, 0.79, 0.80])
    assert p == 1.0


def test_rayleigh_bimodal_season_does_not_fire():
    """
    Two peaks six months apart is strongly non-uniform but does NOT bias the
    midpoint, because the modes cancel. The test is about bias, not seasonality,
    so a null here is correct behaviour rather than a miss.
    """
    rng = random.Random(2)
    w = [(0.1 + rng.gauss(0, 0.05)) % 1 for _ in range(40)] + \
        [(0.6 + rng.gauss(0, 0.05)) % 1 for _ in range(40)]
    R, _Z, _p = rayleigh(w)
    assert R < 0.3


# --- the two failure modes, separately -------------------------------------

def test_seasonal_bias_forces_interval():
    years, precs = dataset(season_centre=0.85, season_conc=0.9)
    r = recommend(years, precs, alignment_length=1800, clock_rate=8.5e-4)
    assert r.policy == "interval"
    assert r.confidence == "strong"
    assert any("seasonal" in x for x in r.reasons)


def test_slow_clock_forces_midpoint_even_with_many_imprecise_tips():
    """
    A slowly-evolving locus cannot resolve within-year timing. Interval
    sampling there adds an unidentifiable parameter. This is the case the
    'always use interval' habit gets wrong.
    """
    years, precs = dataset()
    r = recommend(years, precs, alignment_length=1000, clock_rate=1e-6)
    assert r.policy == "midpoint"
    assert r.resolvable_subs < 0.3
    assert any("not identifiable" in x or "unidentifiable" in x for x in r.reasons)


def test_fast_clock_and_many_imprecise_tips_gives_interval():
    years, precs = dataset()
    r = recommend(years, precs, alignment_length=1800, clock_rate=8.5e-4)
    assert r.policy == "interval"
    assert r.resolvable_subs >= 1.0


def test_few_imprecise_tips_makes_the_choice_moot():
    years, precs = dataset(n_year=3, n_precise=97)
    r = recommend(years, precs, alignment_length=1800, clock_rate=8.5e-4)
    assert r.policy == "midpoint"
    assert r.confidence == "strong"


def test_all_precise_dates_needs_no_policy():
    years = [2000 + i * 0.01 for i in range(50)]
    r = recommend(years, ["day"] * 50, alignment_length=1800, clock_rate=8.5e-4)
    assert r.policy == "midpoint" and r.confidence == "strong"
    assert r.n_imprecise == 0


# --- honesty about what was not tested -------------------------------------

def test_missing_clock_rate_is_flagged_not_guessed():
    years, precs = dataset()
    r = recommend(years, precs, alignment_length=1800, clock_rate=None)
    assert r.confidence == "marginal"
    assert any("no clock rate" in c for c in r.caveats)
    assert r.resolvable_subs is None


def test_too_few_precise_tips_flags_bias_test_as_unassessed():
    years, precs = dataset(n_year=100, n_precise=3)
    r = recommend(years, precs, alignment_length=1800, clock_rate=8.5e-4)
    assert any("too few" in c for c in r.caveats)


def test_missing_at_random_caveat_always_present_when_bias_tested():
    """The bias test's key assumption must never be silent."""
    years, precs = dataset()
    r = recommend(years, precs, alignment_length=1800, clock_rate=8.5e-4)
    assert any("missing at random" in c for c in r.caveats)


def test_recommendation_is_json_serialisable():
    import json
    years, precs = dataset()
    r = recommend(years, precs, alignment_length=1800, clock_rate=8.5e-4)
    json.dumps(r.to_dict())


def test_empty_dataset_does_not_crash():
    r = recommend([], [], alignment_length=1800, clock_rate=8.5e-4)
    assert r.n_total == 0


# --- regression against the real CDV numbers -------------------------------

def test_cdv_clade3_profile_recommends_interval():
    """
    The actual clade-3 profile: 126/162 year-only, 1824 nt, 8.54e-4 subs/site/yr,
    bimodal season so no significant directional bias. The false-precision test
    is what should decide it, not the bias test.
    """
    rng = random.Random(9)
    years = [1992 + rng.random() * 31 + 0.5 for _ in range(126)]
    precs = ["year"] * 126
    for _ in range(36):
        w = rng.choice([rng.gauss(0.1, 0.08), rng.gauss(0.95, 0.08)]) % 1
        years.append(1992 + rng.random() * 31 + w)
        precs.append("day")
    r = recommend(years, precs, alignment_length=1824, clock_rate=8.54e-4)
    assert r.policy == "interval"
    assert r.resolvable_subs == pytest.approx(1.56, abs=0.05)
    assert r.seasonality_p is not None and r.seasonality_p > 0.05

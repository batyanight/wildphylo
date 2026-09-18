"""
Tests for the temporal-signal gate.

The point of these is that the gate must *fail* on data that cannot support a
clock. A gate that always passes is worse than no gate, because it launders
a guess into a result.
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lib.temporal import (  # noqa: E402
    assess, root_to_tip, permutation_test, sampling_diagnostics,
)


def clocklike(n=60, rate=8.5e-4, start=1990, span=30, noise=0.02, seed=1):
    rng = random.Random(seed)
    dates = [start + rng.random() * span for _ in range(n)]
    d0 = min(dates)
    dist = [(d - d0) * rate * (1 + rng.gauss(0, noise)) + 0.001 for d in dates]
    return dates, dist


def test_recovers_a_known_rate():
    """Sanity: on simulated clocklike data the slope must return the input rate."""
    dates, dist = clocklike(rate=8.5e-4)
    r = root_to_tip(dates, dist)
    assert r.slope == pytest.approx(8.5e-4, rel=0.10)
    assert r.r_squared > 0.95


def test_clocklike_data_passes():
    dates, dist = clocklike()
    rep = assess(dates, dist, n_perm=200)
    assert rep.verdict == "pass"


def test_random_distances_fail():
    """No clock signal must produce 'fail', not a confident bad answer."""
    rng = random.Random(7)
    dates = [1990 + rng.random() * 30 for _ in range(60)]
    dist = [rng.random() * 0.05 for _ in range(60)]
    rep = assess(dates, dist, n_perm=500)
    assert rep.verdict == "fail"


def test_inverted_time_axis_is_caught():
    """
    This repo's own history contains a run with the time axis reversed
    (date-backward vs date-forward) that survived 10^8 MCMC states. A negative
    root-to-tip slope is the cheapest possible detector for it.
    """
    dates, dist = clocklike()
    rep = assess(dates, list(reversed(dist)), n_perm=200)
    assert rep.verdict == "fail"
    assert any("non-positive" in n or "not distinguishable" in n for n in rep.notes)


def test_short_sampling_span_fails():
    dates, dist = clocklike(start=2020, span=1)
    rep = assess(dates, dist, n_perm=200, min_span=5.0)
    assert rep.verdict == "fail"
    assert any("span" in n for n in rep.notes)


def test_permutation_test_is_deterministic():
    dates, dist = clocklike()
    a = permutation_test(dates, dist, n_perm=200, seed=42)
    b = permutation_test(dates, dist, n_perm=200, seed=42)
    assert a == b


def test_permutation_p_is_small_for_real_signal():
    dates, dist = clocklike()
    assert permutation_test(dates, dist, n_perm=500) < 0.05


# --- sampling design diagnostics -------------------------------------------

def test_temporal_concentration_flagged():
    """
    The CDV dataset has 90% of tips in 2010-2019 but a TMRCA reported at 1976.
    The gate should say so out loud rather than leaving it to a footnote.
    """
    rng = random.Random(3)
    dates = [2010 + rng.random() * 9 for _ in range(90)] + \
            [1992 + rng.random() * 15 for _ in range(10)]
    d0 = min(dates)
    dist = [(d - d0) * 8.5e-4 + 0.001 for d in dates]
    rep = assess(dates, dist, n_perm=200)
    assert rep.concentration > 0.75
    assert any("10-year window" in n for n in rep.notes)
    assert rep.verdict in ("weak", "fail")


def test_year_only_dates_flagged():
    dates, dist = clocklike(n=100)
    precisions = ["year"] * 80 + ["day"] * 20
    rep = assess(dates, dist, precisions=precisions, n_perm=200)
    assert rep.frac_year_only == pytest.approx(0.8)
    assert any("year-only" in n for n in rep.notes)


def test_sparse_series_flagged():
    dates = [1990.5, 1991.5, 1992.5] + [2020 + i * 0.1 for i in range(30)]
    d0 = min(dates)
    dist = [(d - d0) * 8.5e-4 + 0.001 for d in dates]
    rep = assess(dates, dist, n_perm=200)
    assert rep.occupancy < 0.3
    assert any("sparse" in n or "gappy" in n for n in rep.notes)


def test_diagnostics_on_empty_input():
    assert sampling_diagnostics([]) == {}


def test_report_is_json_serialisable():
    import json
    dates, dist = clocklike()
    rep = assess(dates, dist, n_perm=100)
    json.dumps(rep.to_dict())


# --- time-scaled tree detection --------------------------------------------

def test_time_scaled_tree_recognised():
    """On a BEAST MCC tree branch lengths are years, so the slope must be +1."""
    rng = random.Random(11)
    dates = [1992 + rng.random() * 31 for _ in range(162)]
    root = 1976.0
    dist = [d - root for d in dates]          # exactly elapsed time
    rep = assess(dates, dist, n_perm=200)
    assert any("time-scaled tree" in n for n in rep.notes)


def test_inverted_time_scaled_tree_is_fatal():
    """
    The date-backward bug: tip heights run backwards. Slope is exactly -1.
    This must be caught and must be fatal, not a warning.
    """
    rng = random.Random(11)
    dates = [1992 + rng.random() * 31 for _ in range(162)]
    latest = max(dates)
    dist = [latest - d for d in dates]        # reversed axis
    rep = assess(dates, dist, n_perm=200)
    assert rep.verdict == "fail"
    assert any("TIME AXIS INVERTED" in n for n in rep.notes)

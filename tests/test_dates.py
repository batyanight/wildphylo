"""Tests for collection-date parsing."""

import datetime as dt
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lib.dates import parse_collection_date, to_decimal_year, PRECISION_RANK  # noqa: E402


# --- regression tests for bugs found in the original parser -----------------

@pytest.mark.parametrize("bad", [
    "2013-13-45",   # month 13, day 45
    "2013-02-30",   # 30 February
    "2013-00-00",   # month zero
    "0000-01-01",   # year zero
    "2013-02-29",   # not a leap year
])
def test_malformed_dates_do_not_raise(bad):
    """The original parser raised ValueError here and killed the whole run."""
    result = parse_collection_date(bad)
    assert not result.ok
    assert result.reason, "a failure must carry a machine-readable reason"


@pytest.mark.parametrize("year", ["3013", "1013", "0500", "9999"])
def test_implausible_years_rejected(year):
    """A typo'd year used to parse silently and shift a calibration by centuries."""
    result = parse_collection_date(year)
    assert not result.ok
    assert "year_out_of_range" in result.reason


def test_plausibility_window_is_configurable():
    assert parse_collection_date("1850", min_year=1800).ok
    assert not parse_collection_date("1850").ok


def test_future_years_rejected_but_current_year_allowed():
    this_year = dt.date.today().year
    assert parse_collection_date(str(this_year)).ok
    assert not parse_collection_date(str(this_year + 5)).ok


# --- correctness ------------------------------------------------------------

def test_iso_day_precision():
    r = parse_collection_date("2013-08-14")
    assert r.precision == "day"
    assert r.iso == "2013-08-14"
    assert r.decimal_year == pytest.approx(2013.616, abs=1e-3)


def test_genbank_style_day():
    a = parse_collection_date("14-Aug-2013")
    b = parse_collection_date("2013-08-14")
    assert a.decimal_year == pytest.approx(b.decimal_year)
    assert a.iso == b.iso


def test_year_only_uses_midpoint_not_january():
    """
    Anchoring year-only dates at January biases every such tip half a year
    into the past. With 126/162 year-only tips in the CDV set, that is a
    systematic pull on the TMRCA, not rounding noise.
    """
    assert parse_collection_date("2013").decimal_year == pytest.approx(2013.5)


def test_leap_year_handled():
    assert parse_collection_date("2012-02-29").ok
    assert not parse_collection_date("2013-02-29").ok


def test_range_midpoint_and_uncertainty():
    r = parse_collection_date("2011/2013")
    assert r.precision == "range"
    assert r.decimal_year == pytest.approx(2012.5)
    assert r.uncertainty_years == pytest.approx(1.0)


def test_explicit_unknowns_are_labelled_not_guessed():
    for s in ["unknown", "not collected", "N/A", "missing"]:
        r = parse_collection_date(s)
        assert not r.ok
        assert r.reason == "explicitly_unknown"


def test_empty_input():
    assert not parse_collection_date("").ok
    assert not parse_collection_date(None).ok


def test_uncertainty_scales_with_precision():
    day = parse_collection_date("2013-08-14")
    month = parse_collection_date("2013-08")
    year = parse_collection_date("2013")
    assert day.uncertainty_years < month.uncertainty_years < year.uncertainty_years


def test_precision_rank_orders_correctly():
    assert PRECISION_RANK["day"] < PRECISION_RANK["month"] < PRECISION_RANK["year"]


def test_decimal_year_monotonic_within_year():
    vals = [to_decimal_year(2013, m, 15) for m in range(1, 13)]
    assert vals == sorted(vals)


def test_december_month_midpoint_does_not_overflow():
    """December needs a year rollover to find the month end; easy off-by-one."""
    r = parse_collection_date("2013-12")
    assert 2013.9 < r.decimal_year < 2014.0

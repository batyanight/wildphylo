"""
dates.py — collection-date parsing for tip-dated phylogenetics.

Replaces the inline parser in 02_curate_metadata.py. Three changes that matter:

  1. Malformed dates return a reason, they do not raise. GenBank contains
     "2013-13-45" and "0000-01-01". The old parser crashed the run on these.
  2. Years are bounds-checked against a plausibility window. A typo'd "3013"
     used to parse silently and would have dragged a clock calibration by a
     millennium.
  3. Precision is reported explicitly so downstream code can weight or exclude
     imprecise tips instead of treating year-midpoints as if they were days.

Nothing here is taxon-specific.
"""

from __future__ import annotations

import datetime as dt
import re
from dataclasses import dataclass

MONTHS = {m: i for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun",
     "jul", "aug", "sep", "oct", "nov", "dec"], start=1)}

# Precision ordering, best first. Used by subsampling to prefer precise tips.
PRECISION_RANK = {"day": 0, "month": 1, "range": 2, "year": 3, "none": 4}


@dataclass(frozen=True)
class ParsedDate:
    """Result of parsing one collection-date string."""
    raw: str
    iso: str                 # normalised ISO form, or "" if unparseable
    decimal_year: float | None
    precision: str           # day | month | year | range | none
    uncertainty_years: float  # half-width of the interval the date could fall in
    reason: str = ""         # why it failed, empty on success

    @property
    def ok(self) -> bool:
        return self.decimal_year is not None


def _fail(raw: str, reason: str) -> ParsedDate:
    return ParsedDate(raw, "", None, "none", 0.0, reason)


def to_decimal_year(year: int, month: int | None, day: int | None) -> float:
    """
    Decimal year, using the midpoint of the known interval when month or day
    are missing. Midpoint is the honest choice: a year-only date of 2013 is
    2013.5, not 2013.0, and using 2013.0 biases every year-only tip half a year
    into the past.
    """
    start = dt.date(year, 1, 1)
    days_in_year = (dt.date(year + 1, 1, 1) - start).days
    if month is None:
        return year + 0.5
    if day is None:
        first = dt.date(year, month, 1)
        nxt = dt.date(year + 1, 1, 1) if month == 12 else dt.date(year, month + 1, 1)
        mid = first + (nxt - first) / 2
        return year + (mid - start).days / days_in_year
    return year + (dt.date(year, month, day) - start).days / days_in_year


def _uncertainty(precision: str, year: int, month: int | None) -> float:
    """Half-width, in years, of the interval a date of this precision covers."""
    if precision == "day":
        return 0.5 / 365.25
    if precision == "month":
        days = 31 if month is None else (
            dt.date(year + (month == 12), (month % 12) + 1, 1) - dt.date(year, month, 1)).days
        return days / 2 / 365.25
    if precision == "year":
        return 0.5
    return 0.0


def parse_collection_date(
    raw: str,
    *,
    min_year: int = 1900,
    max_year: int | None = None,
) -> ParsedDate:
    """
    Parse a GenBank-style collection date.

    Accepts ISO (2013-08-14, 2013-08), GenBank (14-Aug-2013, Aug-2013), bare
    years, and slash-separated ranges. Never raises on malformed input.

    min_year / max_year bound what is considered plausible. max_year defaults
    to the current year plus one, which catches typos without rejecting
    legitimately recent samples. Set min_year from your pathogen config: a
    sequence dated before the sequencing era is a data-entry error, not a
    discovery.
    """
    if max_year is None:
        max_year = dt.date.today().year + 1

    if raw is None or not str(raw).strip():
        return _fail(raw or "", "empty")
    s = str(raw).strip()

    low = s.lower()
    if low in {"unknown", "not collected", "missing", "n/a", "na", "none", "null"}:
        return _fail(s, "explicitly_unknown")

    def check_year(y: int) -> str:
        if not (min_year <= y <= max_year):
            return f"year_out_of_range({y}; expected {min_year}-{max_year})"
        return ""

    # Range: "2011/2013" or "2011-01-01/2011-12-31". Midpoint, with the
    # half-span recorded as uncertainty.
    if "/" in s:
        parts = [p.strip() for p in s.split("/") if p.strip()]
        sub = [parse_collection_date(p, min_year=min_year, max_year=max_year)
               for p in parts]
        vals = [p.decimal_year for p in sub if p.ok]
        if not vals:
            bad = next((p.reason for p in sub if p.reason), "unparseable_range")
            return _fail(s, bad)
        lo, hi = min(vals), max(vals)
        return ParsedDate(s, s, (lo + hi) / 2, "range", (hi - lo) / 2)

    # ISO full date
    m = re.fullmatch(r"(\d{4})-(\d{1,2})-(\d{1,2})", s)
    if m:
        y, mo, d = map(int, m.groups())
        if (r := check_year(y)):
            return _fail(s, r)
        try:
            dec = to_decimal_year(y, mo, d)
        except ValueError as e:
            return _fail(s, f"invalid_calendar_date({e})")
        return ParsedDate(s, f"{y:04d}-{mo:02d}-{d:02d}", dec, "day",
                          _uncertainty("day", y, mo))

    # ISO year-month
    m = re.fullmatch(r"(\d{4})-(\d{1,2})", s)
    if m:
        y, mo = map(int, m.groups())
        if (r := check_year(y)):
            return _fail(s, r)
        if not 1 <= mo <= 12:
            return _fail(s, f"invalid_month({mo})")
        return ParsedDate(s, f"{y:04d}-{mo:02d}", to_decimal_year(y, mo, None),
                          "month", _uncertainty("month", y, mo))

    # GenBank day-month-year
    m = re.fullmatch(r"(\d{1,2})-([A-Za-z]{3})-(\d{4})", s)
    if m:
        d, mon, y = int(m.group(1)), m.group(2).lower(), int(m.group(3))
        if mon not in MONTHS:
            return _fail(s, f"unknown_month_name({m.group(2)})")
        if (r := check_year(y)):
            return _fail(s, r)
        mo = MONTHS[mon]
        try:
            dec = to_decimal_year(y, mo, d)
        except ValueError as e:
            return _fail(s, f"invalid_calendar_date({e})")
        return ParsedDate(s, f"{y:04d}-{mo:02d}-{d:02d}", dec, "day",
                          _uncertainty("day", y, mo))

    # GenBank month-year
    m = re.fullmatch(r"([A-Za-z]{3})-(\d{4})", s)
    if m:
        mon, y = m.group(1).lower(), int(m.group(2))
        if mon not in MONTHS:
            return _fail(s, f"unknown_month_name({m.group(1)})")
        if (r := check_year(y)):
            return _fail(s, r)
        mo = MONTHS[mon]
        return ParsedDate(s, f"{y:04d}-{mo:02d}", to_decimal_year(y, mo, None),
                          "month", _uncertainty("month", y, mo))

    # Bare year
    m = re.fullmatch(r"(\d{4})", s)
    if m:
        y = int(m.group(1))
        if (r := check_year(y)):
            return _fail(s, r)
        return ParsedDate(s, f"{y:04d}", to_decimal_year(y, None, None),
                          "year", 0.5)

    return _fail(s, "unrecognised_format")

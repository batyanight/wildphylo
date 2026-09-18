"""
temporal.py — does the data actually support tip-dating?

This replaces the manual TempEst step. Nothing in the original pipeline asked
whether a molecular clock was estimable before writing a clock prior into a
BEAST XML and spending days of compute on it.

Three things are computed:

  1. Root-to-tip regression on a rooted ML tree. Slope is a crude clock-rate
     estimate; R^2 and the correlation say whether divergence tracks time at
     all. The regression is run under both midpoint and best-fitting (residual
     -minimising) roots, because the two disagreeing is itself informative.

  2. A date-shuffling permutation test. Tip dates are permuted N times and the
     observed slope is compared to the null distribution. This is the cheap
     screen that should run *before* the expensive BEAST date-randomisation
     test, not instead of it -- they test different things (this one tests the
     point estimate, the DRT tests the posterior).

  3. Sampling-design diagnostics: temporal span, date precision mix, and how
     concentrated the sampling is in time. A dataset with 90% of its tips in
     one decade can produce a confident-looking TMRCA centuries deep; the
     number that matters is how much of the span is actually populated.

Outputs a dict that the workflow turns into a gate: if the clock is not
supported, the Bayesian step should not run unsupervised.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, asdict, field


@dataclass
class RegressionResult:
    rooting: str
    n: int
    slope: float              # substitutions/site/year
    intercept: float
    r_squared: float
    correlation: float
    x_intercept: float | None  # implied TMRCA / root date


@dataclass
class TemporalReport:
    n_tips: int
    span_years: float
    earliest: float
    latest: float
    precision_counts: dict = field(default_factory=dict)
    frac_year_only: float = 0.0
    occupancy: float = 0.0       # fraction of year-bins in the span with >=1 tip
    concentration: float = 0.0   # largest fraction of tips in any 10-year window
    regressions: list = field(default_factory=list)
    permutation_p: float | None = None
    permutation_n: int = 0
    verdict: str = "unknown"
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["regressions"] = [asdict(r) if not isinstance(r, dict) else r
                            for r in self.regressions]
        return d


def _linreg(xs: list[float], ys: list[float]) -> tuple[float, float, float, float]:
    n = len(xs)
    if n < 3:
        return 0.0, 0.0, 0.0, 0.0
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0 or syy == 0:
        return 0.0, my, 0.0, 0.0
    slope = sxy / sxx
    intercept = my - slope * mx
    r = sxy / math.sqrt(sxx * syy)
    return slope, intercept, r * r, r


def root_to_tip(dates: list[float], distances: list[float],
                rooting: str = "given") -> RegressionResult:
    """Regress root-to-tip genetic distance on sampling date."""
    slope, intercept, r2, r = _linreg(dates, distances)
    x_int = (-intercept / slope) if slope > 0 else None
    return RegressionResult(rooting, len(dates), slope, intercept, r2, r, x_int)


def permutation_test(dates: list[float], distances: list[float],
                     n_perm: int = 1000, seed: int = 12345) -> float:
    """
    One-sided p: how often does a random date assignment produce a slope at
    least as steep as the observed one? A p above ~0.05 means the apparent
    clock signal is indistinguishable from noise, and any TMRCA estimated
    from these dates is an artefact of the tree shape, not of time.
    """
    obs, *_ = _linreg(dates, distances)
    rng = random.Random(seed)
    shuffled = list(dates)
    hits = 0
    for _ in range(n_perm):
        rng.shuffle(shuffled)
        s, *_ = _linreg(shuffled, distances)
        if s >= obs:
            hits += 1
    return (hits + 1) / (n_perm + 1)


def sampling_diagnostics(dates: list[float],
                         precisions: list[str] | None = None) -> dict:
    """
    Describe the temporal design of the dataset. These are the numbers that
    decide whether a deep TMRCA is inference or extrapolation.
    """
    if not dates:
        return {}
    lo, hi = min(dates), max(dates)
    span = hi - lo
    bins = {int(d) for d in dates}
    n_years = max(1, int(hi) - int(lo) + 1)
    occupancy = len(bins) / n_years

    # Densest 10-year window.
    concentration = 0.0
    if span > 0:
        sd = sorted(dates)
        j = 0
        for i, start in enumerate(sd):
            while j < len(sd) and sd[j] < start + 10:
                j += 1
            concentration = max(concentration, (j - i) / len(dates))
    else:
        concentration = 1.0

    out = {
        "earliest": lo, "latest": hi, "span_years": span,
        "occupancy": occupancy, "concentration": concentration,
    }
    if precisions:
        from collections import Counter
        c = Counter(precisions)
        out["precision_counts"] = dict(c)
        out["frac_year_only"] = (c.get("year", 0) + c.get("range", 0)) / len(precisions)
    return out


def assess(dates: list[float], distances: list[float],
           precisions: list[str] | None = None,
           n_perm: int = 1000,
           min_r2: float = 0.10,
           min_span: float = 5.0,
           seed: int = 12345) -> TemporalReport:
    """
    Run the whole assessment and return a verdict of pass / weak / fail.

    The thresholds are deliberately not strict statistical tests -- there is no
    universally right R^2 for a molecular clock. They are tripwires. 'weak'
    means proceed but do not trust a point estimate; 'fail' means the Bayesian
    step should not run unattended.
    """
    diag = sampling_diagnostics(dates, precisions)
    rep = TemporalReport(
        n_tips=len(dates),
        span_years=diag.get("span_years", 0.0),
        earliest=diag.get("earliest", 0.0),
        latest=diag.get("latest", 0.0),
        precision_counts=diag.get("precision_counts", {}),
        frac_year_only=diag.get("frac_year_only", 0.0),
        occupancy=diag.get("occupancy", 0.0),
        concentration=diag.get("concentration", 0.0),
    )

    reg = root_to_tip(dates, distances, "given")
    rep.regressions = [reg]
    rep.permutation_p = permutation_test(dates, distances, n_perm, seed)
    rep.permutation_n = n_perm

    notes = []

    # Time-scaled tree detection. On a tree whose branch lengths are already in
    # years (a BEAST MCC tree), root-to-tip distance IS elapsed time, so the
    # slope must be +1.0. A slope near -1.0 means the time axis is reversed --
    # the date-backward/date-forward error that survived 10^8 MCMC states in
    # this project's own history. Nothing else produces a slope of -1.
    if abs(abs(reg.slope) - 1.0) < 0.05 and reg.r_squared > 0.80:
        if reg.slope < 0:
            notes.append(
                "TIME AXIS INVERTED: this looks like a time-scaled tree, and "
                "the root-to-tip slope is -1.0. Tip heights run backwards "
                "relative to sampling dates. Check date-forward vs "
                "date-backward before interpreting anything downstream")
            rep.verdict = "fail"
            rep.notes = notes
            return rep
        notes.append(
            "slope is +1.0: this is a time-scaled tree, so the slope confirms "
            "the time axis is correct but is not a substitution-rate estimate. "
            "For a rate, run this on the ML tree in substitutions per site")

    if rep.span_years < min_span:
        notes.append(
            f"sampling span is {rep.span_years:.1f} years; below {min_span} "
            "there is rarely enough accumulated change to calibrate a clock")
    if reg.slope <= 0:
        notes.append(
            "root-to-tip slope is non-positive: divergence does not increase "
            "with sampling date, which rules out a simple clock on this rooting")
    if reg.r_squared < min_r2:
        notes.append(
            f"R^2 = {reg.r_squared:.3f} is below {min_r2}; the clock signal is weak")
    if rep.permutation_p is not None and rep.permutation_p > 0.05:
        notes.append(
            f"date-shuffling p = {rep.permutation_p:.3f}: the observed slope is "
            "not distinguishable from randomly assigned dates")
    if rep.concentration > 0.75:
        notes.append(
            f"{rep.concentration:.0%} of tips fall within a single 10-year window; "
            "any TMRCA far outside the sampled period is extrapolation, and its "
            "credible interval will understate that")
    if rep.frac_year_only > 0.5:
        notes.append(
            f"{rep.frac_year_only:.0%} of dates are year-only or ranges; consider "
            "sampling tip dates as intervals in BEAST rather than fixing them at "
            "the year midpoint")
    if rep.occupancy < 0.3 and rep.span_years >= min_span:
        notes.append(
            f"only {rep.occupancy:.0%} of years in the span contain a sample; "
            "the series is sparse and gappy rather than continuous")

    fatal = (reg.slope <= 0
             or (rep.permutation_p is not None and rep.permutation_p > 0.05)
             or rep.span_years < min_span)
    if fatal:
        rep.verdict = "fail"
    elif reg.r_squared < min_r2 or rep.concentration > 0.75:
        rep.verdict = "weak"
    else:
        rep.verdict = "pass"
    rep.notes = notes
    return rep

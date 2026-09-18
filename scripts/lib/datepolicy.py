"""
datepolicy.py — choose between midpoint and interval tip dating, per dataset.

The question
------------
Most GenBank collection dates are year-only. Two ways to handle them:

  midpoint   fix the tip at the middle of its uncertainty window (2013 -> 2013.5)
  interval   let BEAST sample the tip date within [2013.0, 2014.0)

The literature default is "use interval if you can", which is correct but not
actionable -- interval costs chain length and mixing, and on many datasets it
changes nothing. This module decides from the data instead of from habit.

Two independent failure modes
-----------------------------
Midpoint can be wrong in two unrelated ways, and they need separate tests
because a dataset can have one without the other.

  1. BIAS. Midpoint assumes true collection dates are uniform within the year.
     Wildlife sampling is rarely uniform: carcass recovery, outbreak response
     and field-season effort all cluster seasonally. If real dates concentrate
     in autumn, every year-only tip is placed months too early, and the error
     is directional -- it does not average out across tips, it accumulates into
     the rate and TMRCA.

     Testable, because the subset WITH day-precision dates is a sample of the
     same seasonal process. Rayleigh test for circular uniformity on their
     within-year positions.

     Caveat, and it is a real one: that subset is only representative if date
     precision is missing at random. If one lab submitted all the precise dates
     and sampled one season, the test measures that lab, not the pathogen.
     `n_precise` and the source breakdown are reported so this can be judged.

  2. FALSE PRECISION. Even with zero bias, midpoint asserts a date the data
     did not supply. Whether that matters depends on whether the alignment
     could have resolved it:

         lambda = clock_rate * alignment_length   (subs per genome per year)
         resolvable = lambda * 2u                 (subs across the window)

     If resolvable << 1, no amount of sequence data distinguishes January from
     December -- the tip date is unidentifiable and interval sampling just adds
     a poorly-mixing parameter for nothing. Use midpoint.

     If resolvable >~ 1, the alignment carries roughly a substitution's worth
     of signal across the window. Midpoint throws that away AND reports a rate
     HPD tighter than the data support. Use interval.

Pitfalls of each choice, stated so they can be cited in methods
---------------------------------------------------------------
midpoint  - understates uncertainty in rate and TMRCA; HPDs are too narrow
          - directionally biased if sampling is seasonal
          - but: cheap, stable, and comparable with most published analyses
interval  - wider, more honest HPDs
          - more parameters; needs longer chains and ESS must be checked on the
            sampled tip dates, not just on the rate
          - tip-date operators interact badly with strict clocks on shallow
            trees; mixing can collapse
          - not comparable with a previously published midpoint analysis
            without rerunning both
exclude   - only sane when imprecise tips are a small minority AND dropping
            them does not gut a host group or a time period. Almost never the
            right answer for GenBank-derived wildlife datasets, where the
            imprecise tips are usually the old ones carrying the temporal span
"""

from __future__ import annotations

import math
from dataclasses import dataclass, asdict, field


@dataclass
class DatePolicyRecommendation:
    policy: str                       # midpoint | interval | exclude
    confidence: str                   # strong | moderate | marginal
    n_total: int = 0
    n_imprecise: int = 0
    frac_imprecise: float = 0.0
    mean_uncertainty_years: float = 0.0
    n_precise: int = 0
    seasonality_R: float | None = None
    seasonality_p: float | None = None
    mean_within_year: float | None = None
    implied_shift_years: float | None = None
    subs_per_genome_year: float | None = None
    resolvable_subs: float | None = None
    reasons: list = field(default_factory=list)
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def rayleigh(within_year: list[float]) -> tuple[float, float, float]:
    """
    Circular uniformity test on within-year positions (each in [0, 1)).
    Returns (R, Z, p). R near 0 is uniform; R near 1 is tight clustering.

    Note this is a test for UNIMODAL clustering. A bimodal season -- two
    outbreak peaks six months apart -- gives a low R and a non-significant p
    while still being strongly non-uniform. It will not, however, bias the
    midpoint much, because the two modes cancel. That is the behaviour we
    want: this test is asking about bias, not about seasonality per se.
    """
    n = len(within_year)
    if n < 8:
        return 0.0, 0.0, 1.0
    ang = [2 * math.pi * w for w in within_year]
    C = sum(math.cos(a) for a in ang) / n
    S = sum(math.sin(a) for a in ang) / n
    R = math.hypot(C, S)
    Z = n * R * R
    # Standard second-order approximation to the Rayleigh p-value.
    p = math.exp(-Z) * (1 + (2 * Z - Z * Z) / (4 * n))
    return R, Z, min(1.0, max(0.0, p))


def recommend(
    decimal_years: list[float],
    precisions: list[str],
    *,
    alignment_length: int,
    clock_rate: float | None,
    resolvable_floor: float = 0.3,
    resolvable_ceiling: float = 1.0,
    seasonality_alpha: float = 0.05,
    max_shift_years: float = 0.08,
    exclude_threshold: float = 0.10,
) -> DatePolicyRecommendation:
    """
    Recommend a date policy for this dataset.

    clock_rate may be None early in a run; pass the root-to-tip slope once the
    ML tree exists. Without it only the bias test can be run, and the
    recommendation is downgraded to 'marginal'.
    """
    n = len(decimal_years)
    imprecise_kinds = {"year", "range", "unknown"}
    imprecise = [p in imprecise_kinds for p in precisions]
    n_imp = sum(imprecise)
    frac = n_imp / n if n else 0.0

    # Uncertainty half-width by precision class.
    halfwidth = {"day": 0.5 / 365.25, "month": 0.5 / 12, "year": 0.5,
                 "range": 0.5, "unknown": 0.5}
    mean_u = (sum(halfwidth.get(p, 0.5) for p in precisions) / n) if n else 0.0

    rec = DatePolicyRecommendation(
        policy="midpoint", confidence="marginal",
        n_total=n, n_imprecise=n_imp, frac_imprecise=frac,
        mean_uncertainty_years=mean_u,
    )

    if n == 0:
        rec.reasons.append("no dated tips")
        return rec

    if n_imp == 0:
        rec.policy = "midpoint"
        rec.confidence = "strong"
        rec.reasons.append(
            "every tip has day or month precision; there is nothing for "
            "interval sampling to estimate")
        return rec

    # --- test 1: bias -------------------------------------------------------
    precise_w = [y % 1 for y, p in zip(decimal_years, precisions)
                 if p not in imprecise_kinds]
    rec.n_precise = len(precise_w)
    seasonal = False
    if len(precise_w) >= 8:
        R, _Z, p = rayleigh(precise_w)
        mean_w = sum(precise_w) / len(precise_w)
        rec.seasonality_R = R
        rec.seasonality_p = p
        rec.mean_within_year = mean_w
        rec.implied_shift_years = mean_w - 0.5
        seasonal = (p < seasonality_alpha
                    and abs(mean_w - 0.5) > max_shift_years)
        if seasonal:
            rec.reasons.append(
                f"precise-dated tips cluster seasonally (Rayleigh p={p:.3g}) "
                f"with mean within-year position {mean_w:.3f}; midpoint would "
                f"shift every year-only tip by {mean_w - 0.5:+.3f} yr in a "
                "consistent direction")
        else:
            rec.reasons.append(
                f"no directional seasonal bias detected (Rayleigh p={p:.3g}, "
                f"mean within-year position {mean_w:.3f}, implied shift "
                f"{mean_w - 0.5:+.3f} yr)")
        rec.caveats.append(
            f"the bias test uses the {len(precise_w)} tips that carry precise "
            "dates, and assumes date precision is missing at random. If those "
            "came disproportionately from one submitter or one field season, "
            "the test describes that subset rather than the sampling process")
    else:
        rec.caveats.append(
            f"only {len(precise_w)} tips carry precise dates; too few to test "
            "for seasonal bias, so that failure mode is unassessed")

    # --- test 2: false precision -------------------------------------------
    if clock_rate and alignment_length:
        lam = clock_rate * alignment_length
        resolvable = lam * 2 * 0.5      # window is 1 year for year-only dates
        rec.subs_per_genome_year = lam
        rec.resolvable_subs = resolvable
        if resolvable < resolvable_floor:
            rec.reasons.append(
                f"the alignment accumulates {lam:.2f} substitutions per genome "
                f"per year, so a one-year window spans only {resolvable:.2f} "
                "substitutions. Tip dates are not identifiable at this "
                "resolution and interval sampling would add a parameter the "
                "data cannot inform")
        elif resolvable >= resolvable_ceiling:
            rec.reasons.append(
                f"the alignment accumulates {lam:.2f} substitutions per genome "
                f"per year, so a one-year window spans {resolvable:.2f} "
                "substitutions. The data can resolve within-window timing, and "
                "fixing the midpoint discards that information while reporting "
                "an HPD narrower than the data support")
        else:
            rec.reasons.append(
                f"a one-year window spans {resolvable:.2f} substitutions -- "
                "borderline. Interval sampling is defensible but will mostly "
                "widen intervals rather than move estimates")
    else:
        rec.caveats.append(
            "no clock rate supplied, so the false-precision test was not run. "
            "Re-run this after the ML tree exists, passing the root-to-tip "
            "slope")

    # --- decide -------------------------------------------------------------
    resolvable = rec.resolvable_subs

    if frac < exclude_threshold and not seasonal:
        rec.policy = "midpoint"
        rec.confidence = "strong"
        rec.reasons.append(
            f"only {frac:.0%} of tips are imprecise; the choice has little "
            "leverage on the result either way")
        return rec

    if seasonal:
        rec.policy = "interval"
        rec.confidence = "strong"
        rec.reasons.append(
            "seasonal bias is the deciding factor: a directional error does "
            "not average out over tips, it propagates into the rate")
        return rec

    if resolvable is None:
        rec.policy = "midpoint"
        rec.confidence = "marginal"
        rec.reasons.append(
            "defaulting to midpoint pending a clock-rate estimate; this is a "
            "provisional choice, not a finding")
        return rec

    if resolvable >= resolvable_ceiling and frac >= exclude_threshold:
        rec.policy = "interval"
        rec.confidence = "moderate" if frac < 0.5 else "strong"
        rec.reasons.append(
            f"{frac:.0%} of tips are imprecise and the window is resolvable; "
            "midpoint would report false precision across a majority of the "
            "dataset")
        return rec

    if resolvable < resolvable_floor:
        rec.policy = "midpoint"
        rec.confidence = "strong"
        rec.reasons.append(
            "tip dates are unidentifiable at this substitution rate; interval "
            "sampling would cost chain length for no gain")
        return rec

    rec.policy = "interval"
    rec.confidence = "marginal"
    rec.reasons.append(
        "borderline on both tests; interval is the conservative choice because "
        "its failure mode is intervals that are too wide, which is visible, "
        "whereas midpoint's failure mode is intervals that are too narrow, "
        "which is not")
    return rec

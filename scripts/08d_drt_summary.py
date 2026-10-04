#!/usr/bin/env python3
"""
08d_drt_summary.py — the date-randomisation test, summarised and judged.

    python scripts/08d_drt_summary.py --real combined.log \\
        --reps drt/rep*/x_drt.log \\
        --out-json drt_summary.json --out-plot drt.png

What the test is for
--------------------
Root-to-tip regression (04b) asks whether genetic distance increases with
sampling date. It is a screen, and it can be passed by a dataset with no clock
signal at all: with enough tips, a trivial correlation is significant, and
population structure alone produces one when older samples happen to belong to
a diverged lineage. The CDV H tree passed by permutation at R^2 = 0.024.

The date-randomisation test asks the question properly. Shuffle the sampling
dates among the tips, so the tree and the sequences are untouched and only the
date-to-tip assignment is destroyed, then re-estimate the clock rate. If the
real rate is indistinguishable from rates estimated on randomised dates, the
rate is being driven by tree shape rather than by time, and every date in the
analysis — TMRCA, node ages, the whole timescale — is unidentifiable.

The criterion
-------------
Two are in use and they disagree, so both are reported.

  STRICT (Duan et al. 2016; Murray et al. 2016)
      the real rate's 95% HPD must not overlap the HPD of ANY replicate.

  PERMISSIVE (Ramsden et al. 2008, as originally applied)
      the real point estimate must fall outside the range of the replicate
      point estimates.

The verdict uses the strict form, because the permissive one is passed
routinely by datasets that go on to produce unstable TMRCAs: comparing a point
estimate against point estimates ignores that the replicate posteriors are
often very wide, and a real rate sitting just outside a cloud of overlapping
intervals is not distinguishable from them.

The pitfall this script exists to catch
---------------------------------------
Replicates are normally run with a shorter chain than the real analysis — in
this repo's configs, 20M against 100M. A replicate that has not converged has a
WIDE posterior, and a wide posterior overlaps everything. So poor mixing in the
replicates makes the test fail, and the failure looks exactly like absent
temporal signal.

That is the wrong way round: the test is supposed to be hard to pass, not hard
to pass for reasons unrelated to the data. Every replicate's ESS is therefore
computed and reported, and a verdict resting on low-ESS replicates is
downgraded to `inconclusive` rather than being allowed to read as `fail`.

The opposite failure is quieter and is also checked: if the randomisation did
not actually randomise — the same dates reassigned to the same tips, or a
replicate seed reused — the replicate rates cluster tightly on the real one and
the test passes trivially. Identical replicate estimates are flagged.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.beastlog import (                                      # noqa: E402
    read_beast_log, ess, hpd, find_clock_parameter,
)


def rate_summary(path: Path, burnin: float, clock_param: str | None) -> dict:
    """Clock-rate posterior from one log: mean, HPD, ESS, retained samples."""
    header, data = read_beast_log(path)
    cut = int(round(burnin * len(data)))
    kept = data[cut:]
    if len(kept) < 10:
        raise ValueError(
            f"{path}: {len(kept)} samples after a {burnin:.0%} burn-in; "
            "too short to summarise")
    name = find_clock_parameter(header, clock_param)
    col = kept[:, header.index(name)]
    lo, hi = hpd(col)
    return {"path": str(path), "param": name,
            "mean": float(col.mean()), "median": float(np.median(col)),
            "hpd_low": lo, "hpd_high": hi, "ess": round(ess(col), 1),
            "n_kept": int(len(kept)), "n_total": int(len(data))}


def intervals_overlap(a: dict, b: dict) -> bool:
    return a["hpd_low"] <= b["hpd_high"] and b["hpd_low"] <= a["hpd_high"]


def assess(real: dict, reps: list[dict], min_ess: int = 100) -> dict:
    """Compare the real clock rate against the randomised replicates."""
    overlapping = [r for r in reps if intervals_overlap(real, r)]
    rep_means = [r["mean"] for r in reps]

    # The permissive criterion compares the real estimate against the ENVELOPE
    # of the replicate posteriors, not against the range of their point
    # estimates. The range of point estimates shrinks toward the standard
    # error as replicates are added -- with 20 replicates it is far narrower
    # than any single HPD -- so testing against it makes almost any real rate
    # look "outside", and a genuine failure gets rescued by the weaker test.
    lo_m = min(r["hpd_low"] for r in reps)
    hi_m = max(r["hpd_high"] for r in reps)
    outside_range = not (lo_m <= real["mean"] <= hi_m)

    low_ess = [r for r in reps if r["ess"] < min_ess]
    # A replicate posterior wider than the real one is the signature of poor
    # mixing rather than of a flat likelihood, and it is what makes an
    # unconverged replicate overlap everything.
    real_width = real["hpd_high"] - real["hpd_low"]
    wide = [r for r in reps
            if (r["hpd_high"] - r["hpd_low"]) > 3 * real_width]

    reasons: list[str] = []
    if not overlapping:
        verdict = "pass"
        reasons.append(
            f"the real rate's HPD overlaps none of the {len(reps)} replicates: "
            "the clock signal is not reproducible from randomised dates, which "
            "is what temporal signal means.")
    elif outside_range:
        verdict = "weak"
        reasons.append(
            f"{len(overlapping)} of {len(reps)} replicate HPDs overlap the "
            "real one, but the real point estimate falls outside the range of "
            "replicate estimates. This passes the permissive criterion and "
            "fails the strict one; dates estimated here carry wider real "
            "uncertainty than their HPDs show.")
    else:
        verdict = "fail"
        reasons.append(
            f"the real rate is inside the spread of randomised ones "
            f"({len(overlapping)} of {len(reps)} HPDs overlap). The rate is "
            "being driven by tree shape rather than by sampling dates, so the "
            "TMRCA and every node age derived from it are not identifiable.")

    # An unconverged replicate overlaps everything, which turns poor mixing
    # into an apparent absence of temporal signal. That must not read as fail.
    if low_ess and verdict in ("fail", "weak"):
        verdict = "inconclusive"
        reasons.append(
            f"{len(low_ess)} replicate(s) have clock-rate ESS below {min_ess}. "
            "A replicate that has not mixed has a wide posterior that overlaps "
            "everything, so this result cannot be distinguished from poor "
            "sampling. Lengthen the replicate chains "
            "(beast.date_randomisation.chain_length) and rerun before reading "
            "this as absent temporal signal.")
    elif low_ess:
        reasons.append(
            f"{len(low_ess)} replicate(s) have ESS below {min_ess}. The "
            "verdict stands — poor mixing widens replicate HPDs, which makes "
            "a pass harder rather than easier — but the replicates are "
            "under-sampled.")
    if wide:
        reasons.append(
            f"{len(wide)} replicate(s) have an HPD more than 3x the width of "
            "the real one, which usually means the chain never settled rather "
            "than that the randomised likelihood is flat.")

    # The reverse failure: randomisation that did not randomise -- the same
    # permutation reused, or the same seed, so the replicates are one replicate
    # counted many times and the test passes trivially.
    #
    # The test cannot be a fixed threshold on the spread of replicate means.
    # Independent replicates ALSO have tightly clustered means: the spread of a
    # mean goes as SD/sqrt(ESS), so a fixed cutoff measures sample size rather
    # than independence and fires on exactly the well-sampled runs it should
    # ignore.
    #
    # What separates the two is that genuinely different permutations give
    # different point estimates, so the between-replicate spread should EXCEED
    # the Monte Carlo error of a single replicate's mean. When it does not, the
    # replicates are statistically one draw.
    spread = float(np.std(rep_means)) if len(rep_means) > 1 else 0.0
    mc_errors = [((r["hpd_high"] - r["hpd_low"]) / 3.92) / np.sqrt(max(r["ess"], 1))
                 for r in reps]
    expected = float(np.median(mc_errors)) if mc_errors else 0.0
    independent = spread > expected if expected > 0 else True
    if len(reps) > 1 and not independent:
        reasons.append(
            f"the replicate estimates vary by less between replicates "
            f"({spread:.2e}) than the Monte Carlo error within one of them "
            f"({expected:.2e}). Different date permutations give different "
            "point estimates, so this is the signature of replicates sharing a "
            "permutation or a seed rather than of a tight null. A test whose "
            "replicates are one replicate passes trivially.")
        if verdict == "pass":
            verdict = "inconclusive"

    return {"verdict": verdict, "reasons": reasons,
            "n_replicates": len(reps),
            "n_overlapping": len(overlapping),
            "overlapping": [r["path"] for r in overlapping],
            "real_outside_replicate_range": bool(outside_range),
            "replicate_envelope": [lo_m, hi_m],
            "replicate_mean_spread": spread,
            "replicate_mc_error": expected,
            "low_ess_replicates": [r["path"] for r in low_ess],
            "min_ess": min_ess}


def plot(real: dict, reps: list[dict], out: Path) -> bool:
    """Real HPD against the replicate HPDs. Returns False if matplotlib is absent."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    fig, ax = plt.subplots(figsize=(7, 0.3 * len(reps) + 2.2))
    for i, r in enumerate(reps):
        ax.plot([r["hpd_low"], r["hpd_high"]], [i, i], lw=2, color="0.55",
                solid_capstyle="butt")
        ax.plot([r["mean"]], [i], "o", ms=3, color="0.3")
    y = len(reps) + 0.6
    ax.plot([real["hpd_low"], real["hpd_high"]], [y, y], lw=3, color="#b2182b",
            solid_capstyle="butt")
    ax.plot([real["mean"]], [y], "o", ms=5, color="#b2182b")
    ax.axvspan(real["hpd_low"], real["hpd_high"], color="#b2182b", alpha=0.08)
    ax.set_yticks([y] + list(range(len(reps))))
    ax.set_yticklabels(["real"] + [f"rep {i + 1}" for i in range(len(reps))],
                       fontsize=7)
    ax.set_xlabel(f"{real['param']}  (subs/site/year)")
    ax.set_title("Date-randomisation test: real clock rate vs randomised dates",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return True


def report(real: dict, reps: list[dict], res: dict) -> str:
    W = 78
    L = ["=" * W, "DATE-RANDOMISATION TEST", "=" * W,
         f"clock parameter : {real['param']}",
         f"replicates      : {len(reps)}",
         "",
         f"  {'':<14}{'mean':>12}{'95% HPD':>24}{'ESS':>8}",
         "  " + "-" * (W - 4),
         f"  {'REAL':<14}{real['mean']:>12.3e}"
         f"{real['hpd_low']:>12.3e}{real['hpd_high']:>12.3e}{real['ess']:>8.0f}",
         ""]
    for i, r in enumerate(reps, 1):
        flag = ""
        if intervals_overlap(real, r):
            flag = "  <- overlaps real"
        if r["ess"] < res["min_ess"]:
            flag += "  <- low ESS"
        L.append(f"  {'rep ' + str(i):<14}{r['mean']:>12.3e}"
                 f"{r['hpd_low']:>12.3e}{r['hpd_high']:>12.3e}"
                 f"{r['ess']:>8.0f}{flag}")
    L += ["",
          f"  strict     : {res['n_overlapping']} of {res['n_replicates']} "
          f"replicate HPDs overlap the real HPD  "
          f"({'PASS' if res['n_overlapping'] == 0 else 'FAIL'})",
          f"  permissive : the real estimate is "
          f"{'OUTSIDE' if res['real_outside_replicate_range'] else 'INSIDE'} "
          f"the replicate envelope "
          f"[{res['replicate_envelope'][0]:.3e}, "
          f"{res['replicate_envelope'][1]:.3e}]",
          "",
          f"VERDICT: {res['verdict'].upper()}", ""]
    for r in res["reasons"]:
        L.append("  - " + r)
    L += ["", "=" * W,
          "A pass means the timescale is identifiable from the data, not that",
          "it is correct. It says nothing about whether the clock model, the",
          "tree prior or the subsample are appropriate.",
          "=" * W]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--real", required=True, type=Path,
                    help="Combined post-burn-in log from the real analysis "
                         "(what 08b_convergence writes).")
    ap.add_argument("--reps", nargs="+", required=True, type=Path,
                    help="One log per date-randomised replicate.")
    ap.add_argument("--burnin", type=float, default=0.10,
                    help="Applied to the replicates. The combined real log is "
                         "already post-burn-in, so 0 is used for it unless "
                         "--real-burnin says otherwise.")
    ap.add_argument("--real-burnin", type=float, default=0.0)
    ap.add_argument("--clock-param", default=None,
                    help="Column holding the clock rate. Auto-detected from "
                         "the usual BEAST names when not given.")
    ap.add_argument("--min-ess", type=int, default=100,
                    help="Replicate ESS below this makes the verdict "
                         "inconclusive rather than a failure.")
    ap.add_argument("--out-json", type=Path)
    ap.add_argument("--out-plot", type=Path)
    ap.add_argument("--out-md", type=Path)
    ap.add_argument("--exit-nonzero-on-fail", action="store_true")
    a = ap.parse_args()

    missing = [p for p in [a.real, *a.reps] if not p.is_file()]
    if missing:
        print(f"ERROR: no such log file: {missing[0]}", file=sys.stderr)
        return 1
    if len(a.reps) < 5:
        print(f"WARNING: only {len(a.reps)} replicates. Ten is a usual "
              "minimum; with fewer, the replicate range is too poorly "
              "estimated for the permissive criterion to mean much.",
              file=sys.stderr)

    try:
        real = rate_summary(a.real, a.real_burnin, a.clock_param)
        reps = [rate_summary(p, a.burnin, a.clock_param or real["param"])
                for p in a.reps]
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    res = assess(real, reps, min_ess=a.min_ess)
    text = report(real, reps, res)
    print(text)

    if a.out_plot:
        if plot(real, reps, a.out_plot):
            print(f"\nwrote {a.out_plot}")
        else:
            print("\nmatplotlib not available; skipped the plot",
                  file=sys.stderr)
    if a.out_json:
        a.out_json.parent.mkdir(parents=True, exist_ok=True)
        a.out_json.write_text(json.dumps(
            {**res, "real": real, "replicates": reps}, indent=2))
        print(f"wrote {a.out_json}")
    if a.out_md:
        a.out_md.parent.mkdir(parents=True, exist_ok=True)
        a.out_md.write_text("```\n" + text + "\n```\n")
        print(f"wrote {a.out_md}")

    if res["verdict"] == "fail" and a.exit_nonzero_on_fail:
        print("\nexiting non-zero: the timescale is not identifiable.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

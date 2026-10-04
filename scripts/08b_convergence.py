#!/usr/bin/env python3
"""
08b_convergence.py — ESS and between-chain agreement, checked automatically.

    python scripts/08b_convergence.py --logs seed12345/x.log seed54321/x.log \\
        --burnin 0.10 --min-ess 200 \\
        --out-json convergence.json --out-combined combined.log

Replaces `loganalyser -b 10` run by hand, and adds the check loganalyser does
not do: whether the chains agree WITH EACH OTHER. A single chain can look
immaculate — every ESS in the thousands, trace flat as a ruler — while sitting
in a region of parameter space the other chain never visits. ESS measures how
well a chain sampled where it went, not whether it went to the right place.

What is computed
----------------
ESS   per parameter, by Tracer's algorithm, so the numbers match what you would
      read off Tracer or loganalyser rather than being a second opinion. That
      matters: a gate that disagrees with the tool the user checks it against
      gets overridden, and then it is not a gate.

OVL   the overlapping coefficient between chains, per parameter: the area
      shared by their marginal posteriors, 1.0 for identical and 0.0 for
      disjoint. Chosen over comparing HPD intervals because two bimodal
      posteriors can have near-identical 95% intervals while placing their mass
      in different places. Over Gelman-Rubin because R-hat assumes roughly
      normal marginals, which clock rates and tree heights are not, and because
      R-hat on two chains is noisy.

The combined log
----------------
Post-burn-in samples from every chain, pooled, with the Sample column
renumbered so downstream tools read it as one chain. Written here rather than
by logcombiner because the pooling must happen only AFTER the chains are shown
to agree — pooling first would average a converged chain with a stuck one and
produce a posterior that is an artefact of the pooling.

Exit code is non-zero when the gate fails, so a workflow stops rather than
handing an unconverged chain to TreeAnnotator.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config import load, ConfigError            # noqa: E402
from lib.preflight import check_convergence         # noqa: E402

# The log reader, ESS, HPD and overlap live in lib/beastlog.py, shared with
# 08d_drt_summary. One parser and one HPD, so a convergence verdict and a
# date-randomisation verdict describe the same chain.
from lib.beastlog import (                               # noqa: E402
    MAX_LAG, NON_PARAMETERS, read_beast_log, ess, hpd, overlap,
)


def analyse(paths: list[Path], burnin: float, min_ess: int) -> dict:
    chains = []
    for p in paths:
        header, data = read_beast_log(p)
        cut = int(round(burnin * len(data)))
        chains.append({"path": str(p), "header": header,
                       "data": data[cut:], "n_total": len(data),
                       "n_kept": len(data) - cut})
        if len(data) - cut < 10:
            raise ValueError(
                f"{p}: only {len(data) - cut} samples survive a "
                f"{burnin:.0%} burn-in; the chain is too short to assess")

    # Every chain must log the same parameters. A mismatch means the XMLs
    # differ, which is a different analysis, not a different seed -- exactly
    # the mistake of editing one seed's XML and forgetting the other.
    ref = [h for h in chains[0]["header"] if h not in NON_PARAMETERS]
    for c in chains[1:]:
        other = [h for h in c["header"] if h not in NON_PARAMETERS]
        if other != ref:
            only_a = sorted(set(ref) - set(other))
            only_b = sorted(set(other) - set(ref))
            raise ValueError(
                f"chains log different parameters, so they are not replicates "
                f"of one analysis. Only in {chains[0]['path']}: {only_a}; "
                f"only in {c['path']}: {only_b}")

    cols = {name: i for i, name in enumerate(chains[0]["header"])}

    # ESS is reported on the POOLED post-burn-in samples, which is what
    # loganalyser on a combined log gives, plus per chain so a single bad
    # chain is visible rather than being rescued by the other.
    ess_pooled: dict[str, float] = {}
    ess_per_chain: dict[str, list[float]] = {}
    hpds: dict[str, list[float]] = {}
    for name in ref:
        j = cols[name]
        per = [ess(c["data"][:, j]) for c in chains]
        ess_per_chain[name] = [round(e, 1) for e in per]
        # The pooled ESS of independent chains is their sum. Concatenating and
        # running the autocorrelation estimator over the join would read the
        # discontinuity as signal and understate it.
        ess_pooled[name] = round(float(sum(per)), 1)
        pooled = np.concatenate([c["data"][:, j] for c in chains])
        lo, hi = hpd(pooled)
        hpds[name] = [round(float(pooled.mean()), 6), round(lo, 6), round(hi, 6)]

    ovl: dict[str, float] = {}
    if len(chains) > 1:
        for name in ref:
            j = cols[name]
            pairwise = [overlap(chains[a]["data"][:, j], chains[b]["data"][:, j])
                        for a in range(len(chains))
                        for b in range(a + 1, len(chains))]
            ovl[name] = round(min(pairwise), 4)

    return {"chains": chains, "params": ref, "cols": cols,
            "ess": ess_pooled, "ess_per_chain": ess_per_chain,
            "overlap": ovl, "hpd": hpds, "min_ess": min_ess,
            "burnin": burnin}


def write_combined(res: dict, out: Path) -> int:
    """Pool post-burn-in samples with a renumbered Sample column."""
    chains = res["chains"]
    header = chains[0]["header"]
    sample_col = 0 if header[0] in NON_PARAMETERS else None
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out.open("w") as fh:
        fh.write("# pooled post-burn-in samples from: "
                 + ", ".join(c["path"] for c in chains) + "\n")
        fh.write("# burn-in " + f"{res['burnin']:.0%}"
                 + " applied within each chain before pooling\n")
        fh.write("\t".join(header) + "\n")
        for c in chains:
            for row in c["data"]:
                vals = list(row)
                if sample_col is not None:
                    vals[sample_col] = n
                fh.write("\t".join(
                    (f"{int(v)}" if i == sample_col else f"{v:.6g}")
                    for i, v in enumerate(vals)) + "\n")
                n += 1
    return n


def report(res: dict, checks) -> str:
    W = 78
    L = ["=" * W, "CONVERGENCE", "=" * W,
         f"chains        : {len(res['chains'])}",
         f"burn-in       : {res['burnin']:.0%} within each chain"]
    for c in res["chains"]:
        L.append(f"  {c['path']}: {c['n_kept']} of {c['n_total']} samples kept")
    L += ["",
          f"  {'parameter':<28}{'ESS':>10}{'per chain':>20}{'OVL':>8}",
          "  " + "-" * (W - 4)]
    worst = sorted(res["params"], key=lambda p: res["ess"][p])
    for name in worst:
        per = "/".join(f"{e:.0f}" for e in res["ess_per_chain"][name])
        o = res["overlap"].get(name)
        flag = ""
        if res["ess"][name] < res["min_ess"]:
            flag = "  <- ESS"
        elif o is not None and o < 0.5:
            flag = "  <- chains disagree"
        L.append(f"  {name[:28]:<28}{res['ess'][name]:>10.0f}{per:>20}"
                 f"{('' if o is None else f'{o:.2f}'):>8}{flag}")
    L += ["",
          "  ESS is the pooled figure (chains summed). OVL is the smallest",
          "  pairwise overlap of the chains' marginal posteriors: 1.0 means",
          "  they agree, below 0.5 means at least one chain is somewhere the",
          "  others never went.",
          "", checks.render(), "=" * W]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--logs", nargs="+", required=True,
                    help="BEAST .log files, one per chain.")
    ap.add_argument("--burnin", type=float, default=0.10,
                    help="Fraction discarded from the head of EACH chain. "
                         "Applied per chain, not to the pooled set.")
    ap.add_argument("--min-ess", type=int, default=200)
    ap.add_argument("--config", type=Path,
                    help="Pathogen config, for min_ess and the date-policy "
                         "check in preflight.check_convergence.")
    ap.add_argument("--out-json", type=Path)
    ap.add_argument("--out-combined", type=Path)
    ap.add_argument("--out-md", type=Path)
    ap.add_argument("--allow-fail", action="store_true",
                    help="Report and exit 0 even on failure. For looking at a "
                         "chain that is still running.")
    a = ap.parse_args()

    paths = [Path(p) for p in a.logs]
    missing = [p for p in paths if not p.is_file()]
    if missing:
        print(f"ERROR: no such log file: {missing[0]}", file=sys.stderr)
        return 1

    cfg: dict = {"beast": {"min_ess": a.min_ess}}
    if a.config:
        try:
            cfg = load(a.config)
        except ConfigError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1
        cfg.setdefault("beast", {})["min_ess"] = (
            cfg.get("beast", {}).get("min_ess") or a.min_ess)

    try:
        res = analyse(paths, a.burnin, cfg["beast"]["min_ess"])
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    checks = check_convergence(cfg, res["ess"], len(res["chains"]),
                               res["overlap"] or None)
    text = report(res, checks)
    print(text)

    if a.out_combined:
        n = write_combined(res, a.out_combined)
        print(f"\nwrote {a.out_combined}  ({n} pooled samples)")

    if a.out_json:
        a.out_json.parent.mkdir(parents=True, exist_ok=True)
        a.out_json.write_text(json.dumps({
            "verdict": "pass" if checks.ok else "fail",
            "n_chains": len(res["chains"]),
            "burnin": res["burnin"],
            "min_ess": res["min_ess"],
            "logs": [c["path"] for c in res["chains"]],
            "samples_kept": [c["n_kept"] for c in res["chains"]],
            "ess": res["ess"],
            "ess_per_chain": res["ess_per_chain"],
            "overlap": res["overlap"],
            "hpd": res["hpd"],
            "checks": checks.to_dict(),
        }, indent=2))
        print(f"wrote {a.out_json}")

    if a.out_md:
        a.out_md.parent.mkdir(parents=True, exist_ok=True)
        a.out_md.write_text("```\n" + text + "\n```\n")
        print(f"wrote {a.out_md}")

    if not checks.ok and not a.allow_fail:
        print("\nexiting non-zero: the chains have not converged.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

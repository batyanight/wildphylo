#!/usr/bin/env python3
"""
04b_temporal_signal.py — automated root-to-tip regression and clock gate.

Replaces the manual TempEst step. Run before committing compute to BEAST.

    python scripts/04b_temporal_signal.py --tree tree.treefile \\
        --out-json temporal_signal.json --out-md temporal_signal.md

Tip labels are expected as ACCESSION|group|decimal_year, or dates can be
supplied via --metadata.

Exit code is 0 on pass/weak and 1 on fail, so it works as a shell gate as well
as a Snakemake checkpoint.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.temporal import assess  # noqa: E402


def root_to_tip_distances(tree):
    """Distance from root to every tip, in substitutions per site."""
    return {t.name: tree.distance(tree.root, t) for t in tree.get_terminals()}


def date_from_label(label: str) -> float | None:
    try:
        return float(label.rsplit("|", 1)[-1])
    except (ValueError, IndexError):
        return None


def midpoint_root(tree):
    try:
        tree.root_at_midpoint()
    except Exception:
        pass
    return tree


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--tree-format", default="newick")
    ap.add_argument("--metadata", help="TSV with columns name/strain and decimal_year")
    ap.add_argument("--permutations", type=int, default=1000)
    ap.add_argument("--min-r2", type=float, default=0.10)
    ap.add_argument("--min-span", type=float, default=5.0)
    ap.add_argument("--no-midpoint-root", action="store_true")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-md")
    # A FAIL verdict is a result, not a crash. Exiting non-zero made Snakemake
    # treat the rule as failed and DELETE the report -- which is the single
    # most valuable artefact precisely when the gate fails. The workflow reads
    # the verdict from the JSON and decides what to schedule; it does not need
    # the exit code. Standalone shell use can opt back in.
    ap.add_argument("--exit-nonzero-on-fail", action="store_true",
                    help="Return exit code 1 when the verdict is 'fail'. Off by "
                         "default so the report always survives.")
    ap.add_argument("--out-plot")
    a = ap.parse_args()

    from Bio import Phylo
    tree = Phylo.read(a.tree, a.tree_format)
    if not a.no_midpoint_root:
        tree = midpoint_root(tree)

    dists = root_to_tip_distances(tree)

    dates, precisions = {}, {}
    if a.metadata:
        import pandas as pd
        df = pd.read_csv(a.metadata, sep="\t")
        namecol = next((c for c in ("name", "strain", "accession", "taxon")
                        if c in df.columns), df.columns[0])
        for _, r in df.iterrows():
            dates[str(r[namecol])] = r.get("decimal_year")
            precisions[str(r[namecol])] = r.get("date_precision", "unknown")

    xs, ys, ps = [], [], []
    unresolved = []
    for name, d in dists.items():
        dec = dates.get(name)
        if dec is None:
            dec = date_from_label(name)
        if dec is None:
            unresolved.append(name)
            continue
        xs.append(float(dec)); ys.append(d)
        ps.append(precisions.get(name, "unknown"))

    if len(xs) < 3:
        print(f"ERROR: only {len(xs)} tips carry a usable date; cannot regress.",
              file=sys.stderr)
        return 1
    if unresolved:
        print(f"NOTE: {len(unresolved)} tips had no parseable date and were "
              f"excluded from the regression (e.g. {unresolved[:3]})")

    rep = assess(xs, ys, precisions=ps, n_perm=a.permutations,
                 min_r2=a.min_r2, min_span=a.min_span)

    Path(a.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_json).write_text(json.dumps(rep.to_dict(), indent=2))

    r = rep.regressions[0]
    lines = [
        "# Temporal signal assessment", "",
        f"**Verdict: {rep.verdict.upper()}**", "",
        "| Metric | Value |", "|---|---|",
        f"| Tips | {rep.n_tips} |",
        f"| Sampling span | {rep.earliest:.2f} – {rep.latest:.2f} ({rep.span_years:.1f} yr) |",
        f"| Root-to-tip slope | {r.slope:.3e} subs/site/yr |",
        f"| R² | {r.r_squared:.4f} |",
        f"| Correlation | {r.correlation:.4f} |",
        f"| Implied root date | {r.x_intercept:.1f} |" if r.x_intercept else "| Implied root date | undefined |",
        f"| Date-shuffling p | {rep.permutation_p:.4f} ({rep.permutation_n} perms) |",
        f"| Year-bin occupancy | {rep.occupancy:.0%} |",
        f"| Densest 10-yr window | {rep.concentration:.0%} of tips |",
        f"| Year-only or range dates | {rep.frac_year_only:.0%} |",
        "",
    ]
    if rep.notes:
        lines += ["## Cautions", ""] + [f"- {n}" for n in rep.notes] + [""]
    if a.out_md:
        Path(a.out_md).write_text("\n".join(lines))
    print("\n".join(lines))

    if a.out_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.scatter(xs, ys, s=18, alpha=0.7, edgecolor="none")
            lo, hi = min(xs), max(xs)
            ax.plot([lo, hi], [r.intercept + r.slope * lo,
                               r.intercept + r.slope * hi], lw=1.5)
            ax.set_xlabel("Sampling date (decimal year)")
            ax.set_ylabel("Root-to-tip divergence (subs/site)")
            ax.set_title(f"Root-to-tip regression — {rep.verdict.upper()}  "
                         f"(R²={r.r_squared:.3f}, p={rep.permutation_p:.3f})")
            fig.tight_layout()
            Path(a.out_plot).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(a.out_plot, dpi=150)
        except Exception as e:                        # noqa: BLE001
            print(f"NOTE: plot skipped ({e})")

    if rep.verdict == "fail":
        print("\nVerdict is FAIL. The report above has still been written — read "
              "it rather than re-running. A failed gate is a finding about the "
              "data, not a broken step.")
        if a.exit_nonzero_on_fail:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

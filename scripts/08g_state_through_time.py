#!/usr/bin/env python3
"""
08g_state_through_time.py — how many lineages carry each host state, year by year.

Why this exists
---------------
08f_jump_history.py counts transitions and dates them. It cannot tell you whether
a state that dominates the deep tree is a root-era phenomenon or a persistent
feature of the reconstruction, because a count has no denominator: a state can
emit many transitions early simply because that is where the tree is narrow and
every emission is visible.

This script supplies the denominator. It slices the posterior trees at fixed
time points and counts, at each slice, how many lineages are painted with each
state — the discrete-trait version of a lineages-through-time plot. Read it as:

  * a state ENRICHED relative to its tip share only in the deep tree, decaying
    to roughly its tip share once more lineages exist, is being driven by the
    root, not by ongoing presence in the tree
  * a state that holds enrichment above its tip share across decades is a
    structural feature of the topology, which a randomization test then has to
    rule in or out as signal

Enrichment, not raw share, is the quantity to read. A state holding 20% of
lineages means nothing on its own; 20% against a 10% tip share is a 2x
enrichment, and against a 40% tip share is a depletion.

Slices where the tree holds only one or two lineages cannot support a
proportion at all — "40% of 2 lineages" is one lineage — so slices below
--min-lineages are excluded from both the table and the verdict.

Method
------
Each branch spans [parent_year, child_year] and has a state at each end. For a
slice falling in the older half of a branch the parent's state is used, for the
younger half the child's. That is the standard midpoint convention and it is an
approximation: a branch whose endpoints differ made its transition at some
unknown point along its length, and the midpoint is the expectation under a
uniform prior, not the truth. Branches with identical endpoints — the large
majority — are unaffected.

Counts are averaged over the retained posterior trees, so a "mean lineages" of
3.4 means that at that moment the average tree had 3.4 lineages in that state.

Performance
-----------
Cost scales with trees x branches x slices. --thin 10 (the default) keeps every
10th tree, which for a 10,000-tree run leaves ~900 per seed and is far more than
enough: adjacent logged trees are highly autocorrelated and add almost nothing.
Raise it to --thin 1 if you want the full file and have a minute to spare.

Usage
-----
    python scripts/08g_state_through_time.py \
        --trees builds/.../seed*/america2_host_group.trees \
        --burnin 0.10 --thin 10 \
        --out builds/.../state_through_time.tsv

    # focus the verdict lines on one state
    --focus felid
"""

import argparse
import glob as globmod
import sys
from collections import Counter, defaultdict
from pathlib import Path

# The parser lives in lib/traittrees.py, shared with 08f_jump_history. It used
# to be loaded out of that script by importlib, because a module whose name
# starts with a digit cannot be imported normally -- which broke the moment the
# scripts were renamed out of the 12/13 slots. Shared code does not belong in a
# file whose name encodes its position in a pipeline.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.traittrees import (                                    # noqa: E402
    iter_trees, node_depths, label_date, smart_open,
    detect_tag, annotation_value, NUMBER_RE, TREE_LINE_RE,
)

try:
    import numpy as np
except ImportError:
    raise SystemExit(
        "numpy is required for this script.\n"
        "  conda activate cdv-phylo   (numpy is already in environment.yml)\n"
        "  or: pip install numpy"
    )


def main():
    ap = argparse.ArgumentParser(
        description="Lineages carrying each discrete state, through time.",
    )
    ap.add_argument("--trees", nargs="+", required=True)
    ap.add_argument("--burnin", type=float, default=0.10,
                    help="fraction discarded from the START of each file (default 0.10)")
    ap.add_argument("--thin", type=int, default=10,
                    help="keep every Nth post-burnin tree (default 10)")
    ap.add_argument("--tag", default=None,
                    help="annotation key holding the state (default: auto-detect)")
    ap.add_argument("--code-map", default=None,
                    help='BEAST codeMap string, if the trees store numeric codes')
    ap.add_argument("--most-recent", type=float, default=None,
                    help="decimal year of the most recent tip (default: from tip labels)")
    ap.add_argument("--min-lineages", type=float, default=3.0,
                    help="suppress slices holding fewer than this many lineages; "
                         "proportions there are not interpretable (default 3)")
    ap.add_argument("--step", type=float, default=1.0,
                    help="slice interval in years (default 1.0)")
    ap.add_argument("--focus", default=None,
                    help="state to summarise in the verdict lines, e.g. felid")
    ap.add_argument("--out", type=Path, default=None,
                    help="write the year x state table here")
    args = ap.parse_args()

    code_lookup = {}
    if args.code_map:
        for chunk in args.code_map.split(","):
            if "=" not in chunk:
                continue
            name, codes = chunk.split("=", 1)
            codes = codes.strip().split()
            if len(codes) == 1 and codes[0].isdigit():
                code_lookup[codes[0]] = name.strip()

    def decode(st):
        return None if st is None else code_lookup.get(st, st)

    files = []
    for pattern in args.trees:
        hits = ([Path(h) for h in sorted(globmod.glob(pattern))]
                if any(c in pattern for c in "*?[") else [Path(pattern)])
        if not hits:
            sys.exit(f"no file matched: {pattern}")
        files.extend(hits)

    # ---- pass 1: parse retained trees into flat branch lists ---------------- #
    per_tree = []          # [[(y_parent, y_child, state_parent, state_child), ...], ...]
    root_years = []
    tip_census = None
    resolved_tag = args.tag
    kept = 0

    for path in files:
        n_trees = 0
        with smart_open(path) as fh:
            for line in fh:
                if TREE_LINE_RE.match(line.strip()):
                    n_trees += 1
        if n_trees == 0:
            sys.exit(f"{path}: no trees found")
        skip = int(round(args.burnin * n_trees))

        for tree_i, nodes, tag in iter_trees(path, resolved_tag):
            resolved_tag = tag
            if tree_i < skip or (tree_i - skip) % args.thin:
                continue
            depth = node_depths(nodes)
            tip_idx = [i for i, nd in enumerate(nodes) if nd["is_tip"]]
            max_depth = max(depth[i] for i in tip_idx)
            if args.most_recent is not None:
                most_recent = args.most_recent
            else:
                dates = [label_date(nodes[i]["label"]) for i in tip_idx]
                dates = [d for d in dates if d is not None]
                if not dates:
                    sys.exit("Tip labels carry no date and --most-recent was not given.")
                most_recent = max(dates)

            def year(i):
                return most_recent - (max_depth - depth[i])

            if tip_census is None:
                tip_census = Counter()
                for i in tip_idx:
                    tip_census[decode(nodes[i]["state"]) or "?"] += 1

            branches = []
            for i, nd in enumerate(nodes):
                p = nd["parent"]
                if p is None:
                    root_years.append(year(i))
                    continue
                branches.append(
                    (year(p), year(i), decode(nodes[p]["state"]), decode(nd["state"]))
                )
            if not branches:
                continue
            per_tree.append(branches)
            kept += 1

    if not per_tree:
        sys.exit("No trees retained — check --burnin and --thin.")

    # ---- pass 2: accumulate lineage counts on a fixed grid ------------------ #
    y_lo = min(min(b[0] for b in br) for br in per_tree)
    y_hi = max(max(b[1] for b in br) for br in per_tree)
    grid = np.arange(np.floor(y_lo), np.ceil(y_hi) + args.step, args.step)

    states = sorted({s for br in per_tree for b in br for s in (b[2], b[3]) if s})
    counts = {s: np.zeros(len(grid)) for s in states}

    for branches in per_tree:
        for y0, y1, s0, s1 in branches:
            if y1 <= y0:
                continue
            mid = (y0 + y1) / 2.0
            # older half carries the parent state, younger half the child state
            for lo, hi, st in ((y0, mid, s0), (mid, y1, s1)):
                if st is None or hi <= lo:
                    continue
                i0 = int(np.searchsorted(grid, lo, side="left"))
                i1 = int(np.searchsorted(grid, hi, side="right"))
                if i1 > i0:
                    counts[st][i0:i1] += 1.0

    for s in states:
        counts[s] /= len(per_tree)
    total = np.sum([counts[s] for s in states], axis=0)

    # ---- output ------------------------------------------------------------- #
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("year\tstate\tmean_lineages\tproportion\n")
            for k, yr in enumerate(grid):
                if total[k] <= 0:
                    continue
                for s in states:
                    fh.write(f"{yr:.1f}\t{s}\t{counts[s][k]:.4f}\t"
                             f"{counts[s][k] / total[k]:.4f}\n")

    W = 78
    print("=" * W)
    print("STATE THROUGH TIME  (mean lineages per posterior tree)")
    print("=" * W)
    print(f"state annotation key : {resolved_tag}")
    print(f"trees retained       : {kept}  (burnin {args.burnin:.0%}, thin {args.thin})")
    print(f"root year            : mean {np.mean(root_years):.1f}, "
          f"95% range {np.percentile(root_years, 2.5):.1f}-"
          f"{np.percentile(root_years, 97.5):.1f}")
    if tip_census:
        n = sum(tip_census.values())
        print("tips                 : " +
              ", ".join(f"{s} {c} ({c / n:.0%})" for s, c in tip_census.most_common()))

    tip_total = sum(tip_census.values()) if tip_census else 0
    tip_share = {s: (tip_census.get(s, 0) / tip_total if tip_total else float("nan"))
                 for s in states}

    print("\n  " + "year".ljust(7) + "lin".rjust(6) + "   " +
          "".join(s[:8].rjust(10) for s in states))
    suppressed = 0
    for k, yr in enumerate(grid):
        if abs(yr - round(yr)) > 1e-9 or int(round(yr)) % 5:
            continue
        if total[k] < args.min_lineages:
            if total[k] > 0:
                suppressed += 1
            continue
        row = "".join(f"{counts[s][k] / total[k]:>9.0%} " for s in states)
        print(f"  {yr:<7.0f}{total[k]:>6.1f}   {row}")
    print("\n  Columns are the share of lineages in each state at that year.")
    print("  'lin' is the mean number of lineages alive in the tree then.")
    if suppressed:
        print(f"  {suppressed} slice(s) hidden: fewer than {args.min_lineages:.0f} "
              "lineages, so a percentage there describes one or two branches.")

    print("\n  ENRICHMENT (share of lineages / share of tips; 1.0 = as expected)")
    print("  " + "year".ljust(7) + "lin".rjust(6) + "   " +
          "".join(s[:8].rjust(10) for s in states))
    for k, yr in enumerate(grid):
        if abs(yr - round(yr)) > 1e-9 or int(round(yr)) % 5:
            continue
        if total[k] < args.min_lineages:
            continue
        row = ""
        for s in states:
            e = (counts[s][k] / total[k]) / tip_share[s] if tip_share[s] else float("nan")
            row += f"{e:>9.1f} "
        print(f"  {yr:<7.0f}{total[k]:>6.1f}   {row}")
    print("  A state well above 1.0 only in the deep tree, falling to ~1.0 once")
    print("  lineages accumulate, is anchored by its oldest tips rather than")
    print("  persistent. Sustained values above 1.0 are the case for structure.")

    if args.focus:
        f = args.focus
        if f not in states:
            sys.exit(f"--focus {f} not among states: {', '.join(states)}")
        share = counts[f] / np.where(total > 0, total, np.nan)
        valid = [k for k in range(len(grid)) if total[k] >= args.min_lineages]
        if not valid:
            print(f"\n[focus] no slice holds {args.min_lineages:.0f}+ lineages.")
        else:
            ts = tip_share[f]
            print(f"\n--- FOCUS: {f} " + "-" * (W - 12 - len(f)))
            print(f"  tip share: {ts:.0%}  ({tip_census[f]} of {tip_total} tips)")
            print(f"  interpretable window: {grid[valid[0]]:.0f}-{grid[valid[-1]]:.0f} "
                  f"({len(valid)} slices with >= {args.min_lineages:.0f} lineages)")

            # thirds of the interpretable window, weighted by lineages present
            third = max(1, len(valid) // 3)
            blocks = [("earliest", valid[:third]),
                      ("middle", valid[third:2 * third]),
                      ("latest", valid[2 * third:])]
            enr = {}
            print(f"\n  {'window':<10}{'years':>14}{'mean lin':>10}"
                  f"{'share':>9}{'enrichment':>12}")
            for name, ks in blocks:
                if not ks:
                    continue
                w = np.array([total[k] for k in ks])
                sh = float(np.average([share[k] for k in ks], weights=w))
                enr[name] = sh / ts if ts else float("nan")
                print(f"  {name:<10}{f'{grid[ks[0]]:.0f}-{grid[ks[-1]]:.0f}':>14}"
                      f"{float(np.mean(w)):>10.1f}{sh:>9.0%}{enr[name]:>12.1f}x")

            e0, e1 = enr.get("earliest", float("nan")), enr.get("latest", float("nan"))
            print()
            if e0 >= 1.5 and e1 < 1.25:
                print(f"  READ: {e0:.1f}x enrichment in the deep tree decaying to "
                      f"{e1:.1f}x. This state's")
                print("  prominence is carried by its oldest tips. Outbound transition")
                print("  counts from it are a counting consequence of that anchoring,")
                print("  not a rate. Do not report them as transition frequencies.")
            elif e1 >= 1.25:
                print(f"  READ: enrichment persists ({e0:.1f}x early, {e1:.1f}x late).")
                print("  Not explained by deep-node anchoring. Whether it is signal or")
                print("  label noise is now a question for tip-label randomization.")
            else:
                print(f"  READ: no meaningful enrichment at any depth "
                      f"({e0:.1f}x early, {e1:.1f}x late).")
                print("  The state tracks its tip frequency; there is nothing here to")
                print("  explain as transmission structure.")

    print("\n" + "=" * W)
    print("Branch states come from BEAST's sampled reconstruction; this script")
    print("summarises them, it does not re-estimate them. A state can be")
    print("persistent here and still be noise — only a tip-label randomization")
    print("test can establish that the trait carries signal at all.")
    print("=" * W)


if __name__ == "__main__":
    main()

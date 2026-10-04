#!/usr/bin/env python3
"""
08f_jump_history.py — count host-state transitions in BEAST2 DTA trait trees and
record WHEN each one happened.

Why this exists
---------------
Counting transitions alone (X -> Y happened N times) cannot distinguish a real
transmission pattern from a sampling artefact. If the oldest tips in a dataset
all carry one host label, ancestral reconstruction pulls the deep nodes toward
that label, and every lineage descending from those deep nodes is scored as an
outbound transition from it. The count is identical either way; only the
*timing* of the transitions separates the two cases:

  - jumps concentrated near the root, around the oldest sample dates
        -> consistent with deep-node anchoring by the oldest tips
  - jumps spread through the tree, including recent internal nodes
        -> consistent with ongoing transmission between those hosts

This script parses the `.host.trees` file written by beastclassic's
TreeWithTraitLogger (AncestralStateTreeLikelihood), walks every branch of every
posterior tree, and writes one row per state change with the calendar year of
the parent node, the child node, and the branch midpoint.

It is a branch-level counter, not a full Markov-jump history. A branch whose two
endpoints differ is scored as ONE transition even if the underlying continuous
-time process made three. Counts are therefore a lower bound and slightly biased
against long branches. For exact jump histories you need BEAST's complete
history logger; this is the diagnostic you can run on output you already have.

What it reports
---------------
  * tip-date consistency check (parser sanity: reconstructed tip dates vs the
    dates embedded in the tip labels)
  * a TIP CENSUS: how many tips carry each state and over what date range.
    This is the comparison that makes the root state interpretable — a state
    held by few, old tips that nonetheless wins the root is anchoring.
  * posterior distribution of the ROOT state
  * transitions ranked by frequency, with median and IQR of jump year
  * for each transition, the share of its jumps older than a cutoff year and
    the share falling on TERMINAL branches. A jump on a terminal branch
    describes one tip's own label and carries no information about onward
    transmission; a transition that is mostly terminal is not a corridor.
  * optional focus block comparing one transition against all others, with a
    by-decade breakdown that exposes bimodal (mixture) timing

Usage
-----
    python scripts/08f_jump_history.py \
        --trees beast/clade3/run100M_fixed/seed*/clade3_5state.host.trees \
        --burnin 0.10 \
        --out beast/clade3/run100M_fixed/jumps.tsv

    # if the annotation stores numeric codes instead of names, pass the codeMap
    # string straight out of the BEAST XML (userDataType codeMap="...")
    --code-map "domestic_dog=0,mustelid=1,procyonid=2,wild_canid=3,wild_felid=4"

    # ask directly about one transition
    --focus wild_felid:procyonid

Tip labels are expected in the format written by 07_make_beast_xml.py:
    ACCESSION|host_group|decimal_year
The most recent tip date is taken from those labels unless --most-recent is given.
"""

import argparse
import glob as globmod
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

# --------------------------------------------------------------------------- #
# small numeric helpers (kept stdlib-only so this runs anywhere the repo runs)
# --------------------------------------------------------------------------- #


def quantile(sorted_vals, q):
    """Linear-interpolated quantile of an already-sorted list."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def mean(vals):
    return sum(vals) / len(vals) if vals else float("nan")


def sd(vals):
    if len(vals) < 2:
        return float("nan")
    m = mean(vals)
    return (sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5


# --------------------------------------------------------------------------- #
# The parser lives in lib/traittrees.py, shared with 08g_state_through_time.
# A transition count and a lineage share computed by two different parsers
# could not both be true of the same posterior, so there is one parser.
# --------------------------------------------------------------------------- #

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.traittrees import (                                    # noqa: E402
    NUMBER_RE, TREE_LINE_RE, TRANSLATE_ROW_RE,
    smart_open, read_bracket, annotation_value, NewickParser,
    strip_leading_comment, detect_tag, iter_trees,
    node_depths, label_date,
)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #


def main():
    ap = argparse.ArgumentParser(
        description="Height-stratified host-transition counter for BEAST2 DTA trees.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--trees", nargs="+", required=True,
                    help="one or more .host.trees files (.gz ok); pass every seed")
    ap.add_argument("--burnin", type=float, default=0.10,
                    help="fraction of trees discarded from the START of each file (default 0.10)")
    ap.add_argument("--tag", default=None,
                    help="annotation key holding the state (default: auto-detect; "
                         "this repo's XML uses 'location')")
    ap.add_argument("--code-map", default=None,
                    help='BEAST codeMap string, e.g. "domestic_dog=0,mustelid=1,..." '
                         "— only needed if the trees store numeric codes")
    ap.add_argument("--most-recent", type=float, default=None,
                    help="decimal year of the most recent tip (default: read from tip labels)")
    ap.add_argument("--min-prob", type=float, default=None,
                    help="drop jumps where either endpoint's state probability is below this "
                         "(only works if the logger wrote <tag>.prob)")
    ap.add_argument("--cutoff", type=float, default=None,
                    help="year separating 'deep' from 'shallow' jumps "
                         "(default: median of all jump midpoints)")
    ap.add_argument("--focus", default=None,
                    help="transition to interrogate, as FROM:TO, e.g. wild_felid:procyonid")
    ap.add_argument("--out", type=Path, default=None,
                    help="write the per-jump TSV here (default: no file, summary only)")
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

    def decode(state):
        if state is None:
            return None
        return code_lookup.get(state, state)

    files = []
    for pattern in args.trees:
        if any(ch in pattern for ch in "*?["):
            hits = [Path(h) for h in sorted(globmod.glob(pattern))]
        else:
            hits = [Path(pattern)]
        if not hits:
            sys.exit(f"no file matched: {pattern}")
        files.extend(hits)

    jumps = []              # see TSV header below for column order
    tip_census = None       # [(state, decimal_year)] taken from the first tree
    per_tree_counts = []    # jumps per retained tree
    root_states = Counter()
    trees_seen = defaultdict(int)
    trees_used = 0
    date_errors = []
    missing_states = [0]
    resolved_tag = args.tag

    for path in files:
        # first pass only to know how many trees exist, so burn-in is a fraction
        n_trees = 0
        with smart_open(path) as fh:
            for line in fh:
                if TREE_LINE_RE.match(line.strip()):
                    n_trees += 1
        if n_trees == 0:
            sys.exit(f"{path}: no trees found — is this a tree file?")
        skip = int(round(args.burnin * n_trees))
        trees_seen[str(path)] = n_trees

        for tree_i, nodes, tag in iter_trees(path, resolved_tag):
            resolved_tag = tag
            if tree_i < skip:
                continue
            trees_used += 1

            depth = node_depths(nodes)
            tip_idx = [i for i, nd in enumerate(nodes) if nd["is_tip"]]
            max_depth = max(depth[i] for i in tip_idx)

            if args.most_recent is not None:
                most_recent = args.most_recent
            else:
                dates = [label_date(nodes[i]["label"]) for i in tip_idx]
                dates = [d for d in dates if d is not None]
                if not dates:
                    sys.exit(
                        "Tip labels carry no parseable date and --most-recent was not given."
                    )
                most_recent = max(dates)

            def year(i):
                return most_recent - (max_depth - depth[i])

            # parser sanity check: do reconstructed tip dates match the labels?
            if trees_used <= 5:
                for i in tip_idx:
                    d = label_date(nodes[i]["label"])
                    if d is not None:
                        date_errors.append(abs(year(i) - d))

            if tip_census is None:
                tip_census = []
                for i in tip_idx:
                    st = decode(nodes[i]["state"])
                    if st is None:  # fall back to the host field of the label
                        parts = (nodes[i]["label"] or "").split("|")
                        st = parts[1] if len(parts) > 2 else "?"
                    tip_census.append((st, label_date(nodes[i]["label"])))

            # tolerate a unary stem above the true root (some loggers emit one)
            r, stem = 0, set()
            while nodes[r]["state"] is None and nodes[r].get("nchild") == 1:
                stem.add(r)
                r = next(i for i, nd in enumerate(nodes) if nd["parent"] == r)
            root_states[decode(nodes[r]["state"])] += 1
            missing_states[0] += sum(
                1 for i, nd in enumerate(nodes)
                if nd["state"] is None and i != r and i not in stem
            )

            n_here = 0
            for i, nd in enumerate(nodes):
                p = nd["parent"]
                if p is None:
                    continue
                s_from, s_to = decode(nodes[p]["state"]), decode(nd["state"])
                if s_from is None or s_to is None or s_from == s_to:
                    continue
                if args.min_prob is not None:
                    pp, cp = nodes[p]["prob"], nd["prob"]
                    if pp is not None and pp < args.min_prob:
                        continue
                    if cp is not None and cp < args.min_prob:
                        continue
                y_p, y_c = year(p), year(i)
                jumps.append(
                    (path.name, tree_i, s_from, s_to, y_p, y_c, (y_p + y_c) / 2.0,
                     nd["length"], nodes[p]["prob"], nd["prob"], nd["is_tip"])
                )
                n_here += 1
            per_tree_counts.append(n_here)

    if not jumps:
        sys.exit("No transitions found. Check --tag / --code-map against the tree file.")

    # ---------------------------------------------------------------- output #
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("file\ttree\tfrom\tto\tparent_year\tchild_year\tmid_year\t"
                     "branch_length\tparent_prob\tchild_prob\tchild_is_tip\n")
            for row in jumps:
                fh.write(
                    "{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.6f}\t{}\t{}\t{}\n".format(
                        row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7],
                        "" if row[8] is None else f"{row[8]:.4f}",
                        "" if row[9] is None else f"{row[9]:.4f}",
                        int(row[10]),
                    )
                )

    W = 78
    print("=" * W)
    print("HOST-TRANSITION HISTORY  (branch-level, height-stratified)")
    print("=" * W)
    print(f"state annotation key : {resolved_tag}")
    for f, n in trees_seen.items():
        print(f"  {f}: {n} trees, {int(round(args.burnin * n))} discarded as burn-in")
    print(f"trees analysed       : {trees_used}")
    print(f"transitions recorded : {len(jumps)}")
    print(f"per tree             : mean {mean(per_tree_counts):.2f}, sd {sd(per_tree_counts):.2f}")
    if missing_states[0]:
        print(f"WARNING              : {missing_states[0]} nodes carried no "
              f"'{resolved_tag}' state and were skipped — transition counts are "
              "incomplete. Check --tag.")
    if date_errors:
        print(f"tip-date check       : max |reconstructed - label| = {max(date_errors):.4f} yr "
              f"({'OK' if max(date_errors) < 0.01 else 'MISMATCH — heights are suspect'})")

    census = defaultdict(list)
    for st, yr in (tip_census or []):
        census[st].append(yr)
    n_tips = sum(len(v) for v in census.values())

    print("\n--- TIP CENSUS " + "-" * (W - 15))
    print(f"  {'state':<16}{'tips':>6}{'% tips':>8}{'oldest':>10}{'newest':>10}")
    oldest_overall, oldest_states = None, []
    for st, yrs in sorted(census.items(), key=lambda kv: -len(kv[1])):
        ys = [y for y in yrs if y is not None]
        lo = min(ys) if ys else float("nan")
        hi = max(ys) if ys else float("nan")
        print(f"  {str(st):<16}{len(yrs):>6}{len(yrs) / n_tips:>8.1%}{lo:>10.1f}{hi:>10.1f}")
        if ys and (oldest_overall is None or lo < oldest_overall):
            oldest_overall = lo
    if oldest_overall is not None:
        for st, yrs in census.items():
            ys = [y for y in yrs if y is not None]
            if ys and min(ys) <= oldest_overall + 1.0:
                oldest_states.append(st)

    print("\n--- ROOT STATE (posterior) " + "-" * (W - 27))
    total_roots = sum(root_states.values())
    top_state, top_n = root_states.most_common(1)[0]
    for state, n in root_states.most_common():
        print(f"  {str(state):<16} {n / total_roots:6.1%}  ({n})")
    tip_share = len(census.get(top_state, [])) / n_tips if n_tips else float("nan")
    print(f"\n  {top_state} holds the root in {top_n / total_roots:.0%} of trees "
          f"while carrying {tip_share:.0%} of the tips.")
    if oldest_states:
        print(f"  Oldest tips (within 1 yr of {oldest_overall:.1f}): "
              f"{', '.join(sorted(oldest_states))}")
    if top_state in oldest_states and tip_share < 0.25:
        print("  READ: a minority state that also anchors the oldest tips is winning")
        print("  the root. Outbound transitions from it are inflated by that anchoring;")
        print("  treat its rank in the table as an upper bound, not an estimate.")

    mids = sorted(j[6] for j in jumps)
    cutoff = args.cutoff if args.cutoff is not None else quantile(mids, 0.5)

    by_type = defaultdict(list)
    term_type = defaultdict(int)
    for j in jumps:
        by_type[(j[2], j[3])].append(j[6])
        if j[10]:
            term_type[(j[2], j[3])] += 1

    print(f"\n--- TRANSITIONS, ranked  (cutoff for 'deep' = {cutoff:.1f}) " + "-" * 12)
    print(f"  {'transition':<30}{'n':>7}{'%':>7}{'median':>8}{'IQR':>16}"
          f"{'%deep':>7}{'%term':>7}")
    for (a, b), yrs in sorted(by_type.items(), key=lambda kv: -len(kv[1])):
        ys = sorted(yrs)
        pct = len(ys) / len(jumps)
        deep = sum(1 for y in ys if y < cutoff) / len(ys)
        term = term_type[(a, b)] / len(ys)
        print(f"  {a + ' -> ' + b:<30}{len(ys):>7}{pct:>7.1%}"
              f"{quantile(ys, 0.5):>8.1f}"
              f"{f'{quantile(ys, 0.25):.1f}-{quantile(ys, 0.75):.1f}':>16}"
              f"{deep:>7.0%}{term:>7.0%}")
    print("  %deep = share older than the cutoff. %term = share landing on a")
    print("  terminal branch, i.e. explaining one tip's label rather than any")
    print("  onward transmission.")

    if args.focus:
        try:
            fa, fb = args.focus.split(":")
        except ValueError:
            sys.exit("--focus must look like FROM:TO, e.g. wild_felid:procyonid")
        ys = sorted(by_type.get((fa, fb), []))
        focus_rows = [j for j in jumps if (j[2], j[3]) == (fa, fb)]
        if not ys:
            print(f"\n[focus] no {fa} -> {fb} transitions recorded.")
        else:
            others = sorted(y for (a, b), v in by_type.items() if (a, b) != (fa, fb) for y in v)
            print(f"\n--- FOCUS: {fa} -> {fb} " + "-" * (W - 13 - len(fa) - len(fb)))
            print(f"  n = {len(ys)}  ({len(ys) / len(jumps):.1%} of all transitions)")
            print(f"  jump year   median {quantile(ys, 0.5):.1f}, "
                  f"IQR {quantile(ys, 0.25):.1f}-{quantile(ys, 0.75):.1f}, "
                  f"range {ys[0]:.1f}-{ys[-1]:.1f}")
            print(f"  all others  median {quantile(others, 0.5):.1f}, "
                  f"IQR {quantile(others, 0.25):.1f}-{quantile(others, 0.75):.1f}")
            oldest_decade = sum(1 for y in ys if y < ys[0] + 10) / len(ys)
            print(f"  share within 10 yr of its own oldest jump: {oldest_decade:.0%}")
            n_term = sum(1 for j in focus_rows if j[10])
            print(f"  on terminal branches: {n_term / len(focus_rows):.0%}")

            print("\n  by decade (count, and share on terminal branches):")
            lo_d = int(min(ys) // 10 * 10)
            hi_d = int(max(ys) // 10 * 10)
            for d in range(lo_d, hi_d + 10, 10):
                rows = [j for j in focus_rows if d <= j[6] < d + 10]
                if not rows:
                    continue
                bar = "#" * max(1, round(40 * len(rows) / len(focus_rows)))
                t = sum(1 for j in rows if j[10]) / len(rows)
                print(f"    {d}s {len(rows):>7} ({len(rows) / len(focus_rows):>4.0%}) "
                      f"term {t:>4.0%}  {bar}")
            print("\n  Two separated humps mean two processes are being pooled: deep")
            print("  jumps from root anchoring, and shallow ones from something else.")
            print("  If the shallow hump is mostly terminal, it is tip labels, not spread.")

    print("\n" + "=" * W)
    print("Branch-level counts are a LOWER BOUND: a branch with two endpoints in")
    print("different states is scored once regardless of how many jumps the")
    print("continuous-time process actually made along it.")
    print("=" * W)


if __name__ == "__main__":
    main()

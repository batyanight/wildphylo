#!/usr/bin/env python3
"""
05b_segment_congruence.py — compare segment topologies before treating them as
one history.

Run after the per-segment ML trees exist and before subsampling. For an
unsegmented pathogen this step does not run at all.

    python scripts/05b_segment_congruence.py \\
        --trees tree/seg2_ml.treefile tree/seg6_ml.treefile tree/seg10_ml.treefile \\
        --names seg2 seg6 seg10 \\
        --out-json congruence.json --out-md congruence.md \\
        --out-tsv reassortant_candidates.tsv

Exit code is 0 when concatenation is legitimate or the verdict is merely mixed,
and 1 when segments are decoupled AND the config asked to concatenate anyway.
Being decoupled is not itself an error -- for BTV it is the expected result.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.segments import assess_congruence  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", nargs="+", required=True)
    ap.add_argument("--names", nargs="+",
                    help="Segment names, in the same order as --trees. "
                         "Defaults to the tree filenames.")
    ap.add_argument("--tree-format", default="newick")
    ap.add_argument("--collapse-support", type=float, default=70.0,
                    help="Ignore splits below this support before comparing. "
                         "Set to a negative number to disable.")
    ap.add_argument("--congruent-max-rf", type=float, default=0.15)
    ap.add_argument("--decoupled-min-rf", type=float, default=0.40)
    ap.add_argument("--neighbour-k", type=int, default=5)
    ap.add_argument("--reassortant-threshold", type=float, default=0.6)
    ap.add_argument("--require-concatenation", action="store_true",
                    help="Fail if the segments cannot legitimately be joined.")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-md")
    ap.add_argument("--out-tsv", help="Reassortant candidates, one per row.")
    a = ap.parse_args()

    names = a.names or [Path(t).stem for t in a.trees]
    if len(names) != len(a.trees):
        print(f"ERROR: {len(a.trees)} trees but {len(names)} names",
              file=sys.stderr)
        return 2
    if len(a.trees) < 2:
        print("ERROR: at least two segment trees are needed; congruence "
              "cannot be assessed from one", file=sys.stderr)
        return 2

    from Bio import Phylo
    trees = {}
    for name, path in zip(names, a.trees):
        try:
            trees[name] = Phylo.read(path, a.tree_format)
        except Exception as e:                       # noqa: BLE001
            print(f"ERROR: could not read {path}: {e}", file=sys.stderr)
            return 2

    support = None if a.collapse_support < 0 else a.collapse_support
    res = assess_congruence(
        trees,
        collapse_support=support,
        congruent_max_rf=a.congruent_max_rf,
        decoupled_min_rf=a.decoupled_min_rf,
        neighbour_k=a.neighbour_k,
        reassortant_threshold=a.reassortant_threshold,
    )

    Path(a.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_json).write_text(json.dumps(res.to_dict(), indent=2))

    lines = [
        "# Segment congruence", "",
        f"**Verdict: {res.verdict.upper()}** — "
        f"concatenation {'permitted' if res.concatenate else 'refused'}", "",
        "| Metric | Value |", "|---|---|",
        f"| Segments | {', '.join(res.segments)} |",
        f"| Taxa present in every segment | {res.n_shared_taxa} |",
        f"| Mean normalised RF | {res.mean_rf:.4f} |",
        f"| Max normalised RF | {res.max_rf:.4f} |",
        f"| Reassortant candidates | {len(res.reassortant_candidates)} |",
        "",
    ]
    if res.pairwise_rf:
        lines += ["## Pairwise topological distance", "",
                  "| Pair | normalised RF |", "|---|---|"]
        # Pair keys contain "|", which would break the Markdown table.
        lines += [f"| {k.replace('|', ' vs ')} | {v:.4f} |"
                  for k, v in sorted(res.pairwise_rf.items())]
        lines.append("")
    if res.reassortant_candidates:
        lines += ["## Reassortant candidates", "",
                  "Each sits among different relatives depending on the "
                  "segment. In a discrete-trait host analysis such a taxon "
                  "produces a host transition that reflects reassortment "
                  "rather than transmission.", "",
                  "| Taxon | Incongruence |", "|---|---|"]
        lines += [f"| {c['taxon']} | {c['incongruence']:.3f} |"
                  for c in res.reassortant_candidates[:25]]
        if len(res.reassortant_candidates) > 25:
            lines.append(f"| … and {len(res.reassortant_candidates) - 25} more | |")
        lines.append("")
    if res.notes:
        lines += ["## Notes", ""] + [f"- {n}" for n in res.notes] + [""]

    # Write every output BEFORE printing. If stdout is a pipe that closes
    # early (a `| head` in a terminal, a truncated CI log), a print that comes
    # first takes the process down with SIGPIPE and the files are never
    # written -- silently, since the exit status looks fine.
    if a.out_md:
        Path(a.out_md).write_text("\n".join(lines))
    if a.out_tsv:
        Path(a.out_tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(a.out_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["taxon", "incongruence"])
            for c in res.reassortant_candidates:
                w.writerow([c["taxon"], c["incongruence"]])

    try:
        print("\n".join(lines))
    except BrokenPipeError:
        pass

    if a.require_concatenation and not res.concatenate:
        print("\nERROR: concatenation was required but the segments do not "
              "support it.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

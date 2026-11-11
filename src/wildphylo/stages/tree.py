"""tree stage: cut the ML tree into clades and test each for temporal signal.

``wildphylo tree clades`` exists because a global tree of one gene across every
lineage of a pathogen is usually NOT a tip-dating dataset. Root-to-tip
regression on the full CDV H tree returned R^2 = 0.024 with an implied root of
1656: significant by permutation only because 1,126 tips make a trivial
correlation significant, and nonsensical because the tree is saturated (total
length 16.5 subs/site, with IQ-TREE warning about long pairwise distances).

Clock signal lives WITHIN lineages, below saturation. This step finds those
clades from the tree itself rather than from a reference panel, because a
reference panel that covers only part of the diversity leaves most of the data
unclassified and gives confidently wrong calls for the rest.

How clades are cut: descend from the root and take each maximal subtree whose
internal diversity (maximum tip-to-tip patristic distance) falls below
``--max-diversity``. That threshold is the operational definition of "below
saturation": groups tight enough that distance still tracks time approximately
linearly. Subtrees smaller than ``--min-size`` are reported but not carried
forward; a clade with a handful of tips cannot support a clock estimate however
clean it looks.

This is a data-driven grouping, not a taxonomy. Name the clades afterwards by
seeing which reference sequences land in which group (``--references``). Naming
them before you know which groups exist is how a partial reference panel
produces a confident mislabelling.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any

from wildphylo.cli import step
from wildphylo.temporal import assess

if TYPE_CHECKING:
    import argparse


def _say(text: str = "") -> None:
    """Write one line of the step's report to stdout."""
    print(text, file=sys.stdout)  # noqa: T201 -- the report is the output


def tip_date(name: str) -> float | None:
    """The decimal year in the last ``|``-separated field of a tip label, if numeric."""
    try:
        return float(name.rsplit("|", 1)[1])
    except (ValueError, IndexError):
        return None


def clade_support(clade: Any) -> tuple[float | None, float | None]:
    """Return (support of the clade's own stem, minimum support among its internal nodes).

    A clade cut purely on diversity can be an artefact of the tree search. The
    stem value says whether the group itself is supported; the internal minimum
    says whether its structure is. A clade held together by a single weak node
    is one to treat carefully before committing days of BEAST to it.

    Returns (None, None) when the tree carries no support values: a tree built
    without bootstraps is not a tree with zero support, and conflating the two
    would flag every clade on such a tree as suspect.
    """

    def conf(node: Any) -> float | None:
        c = getattr(node, "confidence", None)
        return None if c is None else float(c)

    stem = conf(clade)
    inner = [v for v in (conf(n) for n in clade.get_nonterminals()) if v is not None]
    return stem, (min(inner) if inner else None)


def max_pairwise(tree: Any, clade: Any) -> float:
    """Maximum tip-to-tip distance inside a clade.

    Computed as the two largest depths below the clade's own root, which is the
    diameter for an additive tree and avoids the O(n^2) all-pairs cost that
    makes this unusable on a few thousand tips.
    """
    depths = tree.depths(unit_branch_lengths=False)
    base = depths.get(clade, 0.0)
    below = [depths[t] - base for t in clade.get_terminals() if t in depths]
    if len(below) < 2:
        return 0.0
    below.sort(reverse=True)
    return float(below[0] + below[1])


def cut(tree: Any, max_diversity: float, min_size: int) -> tuple[list[Any], list[Any]]:
    """Return (kept, small): maximal subtrees whose diameter is under the threshold."""
    kept: list[Any] = []
    small: list[Any] = []
    stack = [tree.root]
    while stack:
        node = stack.pop()
        tips = node.get_terminals()
        if len(tips) < 2:
            small.append(node)
            continue
        if max_pairwise(tree, node) <= max_diversity:
            (kept if len(tips) >= min_size else small).append(node)
            continue
        if not node.clades:
            small.append(node)
            continue
        stack.extend(node.clades)
    return kept, small


def _load_references(path: str) -> dict[str, str]:
    """Read ``accession -> lineage`` for the nucleotide rows of a lineage reference table."""
    refs: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 4 and parts[3].strip() == "nuc":
            refs[parts[0].strip()] = parts[1].strip()
    _say(f"{len(refs)} usable nucleotide references loaded")
    return refs


def _prune_short_tips(tree: Any, qc_path: str, min_ungapped: int) -> None:
    """Drop tips with less than ``min_ungapped`` ungapped nt according to the alignment QC table.

    Ragged coverage makes root-to-tip measure sequence length rather than
    elapsed time, so this has to happen before the cut, not after.
    """
    import pandas as pd

    qc = pd.read_csv(qc_path, sep="\t")
    keep = set(qc.loc[qc["ungapped_len"] >= min_ungapped, "name"])
    drop = [t for t in tree.get_terminals() if t.name not in keep]
    for t in drop:
        tree.prune(t)
    _say(f"pruned {len(drop)} tips below {min_ungapped} nt; {len(tree.get_terminals())} remain")


def _assess_clades(
    tree: Any, kept: list[Any], refs: dict[str, str], args: argparse.Namespace
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Fit root-to-tip within each kept clade; return the summary rows and per-clade reports."""
    rows: list[dict[str, Any]] = []
    reports: dict[str, Any] = {}
    for i, clade in enumerate(sorted(kept, key=lambda c: -len(c.get_terminals())), 1):
        names = [t.name for t in clade.get_terminals()]
        xs: list[float] = []
        ys: list[float] = []
        for t in clade.get_terminals():
            d = tip_date(t.name)
            if d is not None:
                xs.append(d)
                ys.append(tree.distance(clade, t))
        if len(xs) < 3:
            continue
        rep = assess(xs, ys, n_perm=args.permutations)
        reg = rep.regressions[0]
        hits = sorted({refs[n.split("|")[0]] for n in names if n.split("|")[0] in refs})
        stem, min_inner = clade_support(clade)
        row = {
            "clade": f"clade_{i:02d}",
            "n_tips": len(names),
            "n_dated": len(xs),
            "diversity": round(max_pairwise(tree, clade), 4),
            "span_years": round(rep.span_years, 1),
            "earliest": round(rep.earliest, 2),
            "latest": round(rep.latest, 2),
            "slope": f"{reg.slope:.3e}",
            "r_squared": round(reg.r_squared, 4),
            "perm_p": round(rep.permutation_p or 0.0, 4),
            "implied_root": (round(reg.x_intercept, 1) if reg.x_intercept else ""),
            "stem_support": "" if stem is None else round(stem, 1),
            "min_internal_support": "" if min_inner is None else round(min_inner, 1),
            "verdict": rep.verdict,
            "reference_lineages": ";".join(hits),
        }
        rows.append(row)
        reports[row["clade"]] = rep.to_dict()
        if args.out_dir:
            tips_dir = Path(args.out_dir)
            tips_dir.mkdir(parents=True, exist_ok=True)
            (tips_dir / f"{row['clade']}_tips.txt").write_text("\n".join(names) + "\n")
    return rows, reports


def _write_clade_table(rows: list[dict[str, Any]], out_tsv: str) -> None:
    """Write the clade summary TSV (empty file when no clade qualified)."""
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with Path(out_tsv).open("w", newline="") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)


def _report_clades(rows: list[dict[str, Any]], min_support: float | None) -> None:
    """Print the clade table, support warnings and the candidates for Bayesian analysis."""
    have_support = any(r["stem_support"] != "" for r in rows)
    hdr = (
        f"{'clade':<10}{'tips':>6}{'div':>8}{'span':>7}{'slope':>12}"
        f"{'R2':>8}{'p':>8}{'root':>8}"
        + (f"{'supp':>7}" if have_support else "")
        + f"  {'verdict':<8} refs"
    )
    _say(hdr)
    _say("-" * len(hdr))
    for r in rows:
        _say(
            f"{r['clade']:<10}{r['n_tips']:>6}{r['diversity']:>8.3f}"
            f"{r['span_years']:>7.0f}{r['slope']:>12}{r['r_squared']:>8.3f}"
            f"{r['perm_p']:>8.3f}{r['implied_root']!s:>8}"
            + (f"{r['stem_support']!s:>7}" if have_support else "")
            + f"  {r['verdict']:<8} {r['reference_lineages']}"
        )

    if have_support and min_support is not None:
        weak_stem = [
            r for r in rows if r["stem_support"] != "" and r["stem_support"] < min_support
        ]
        if weak_stem:
            _say(f"\nClades with stem support below {min_support:g}:")
            for r in weak_stem:
                _say(
                    f"  {r['clade']}: {r['stem_support']} — the group itself "
                    "may be an artefact of the tree search, so a clock fitted "
                    "to it is fitted to a group that may not exist"
                )
    elif not have_support:
        _say(
            "\nNOTE: the tree carries no support values, so clade reliability "
            "cannot be assessed. Build with -B 1000 if a clade is going to "
            "carry a Bayesian analysis."
        )

    good = [r for r in rows if r["verdict"] in ("pass", "weak")]
    _say(f"\n{len(good)} of {len(rows)} clades show usable temporal signal.")
    if good:
        _say("Candidates for Bayesian analysis, largest first:")
        for r in sorted(good, key=lambda x: -x["n_tips"])[:5]:
            lineages = f"[{r['reference_lineages']}]" if r["reference_lineages"] else "(unnamed)"
            _say(
                f"  {r['clade']}: {r['n_tips']} tips, {r['span_years']:.0f} yr, "
                f"R2={r['r_squared']:.3f}, root~{r['implied_root']} {lineages}"
            )
    else:
        _say(
            "No clade passed. Try a smaller --max-diversity: the groups are "
            "probably still spanning too much divergence to be clocklike."
        )
    _say(
        "\nThese groups come from the tree, not from a taxonomy. Name them by "
        "which references fall where, and say so in the methods."
    )


def add_arguments_clades(parser: argparse.ArgumentParser) -> None:
    """``wildphylo tree clades --tree TREE --out-tsv TSV [--max-diversity D] [...]``."""
    parser.add_argument("--tree", required=True)
    parser.add_argument("--tree-format", default="newick")
    parser.add_argument("--qc", help="alignment QC TSV, to filter by coverage first")
    parser.add_argument(
        "--min-ungapped",
        type=int,
        default=0,
        help="Drop tips with less than this much ungapped sequence. "
        "Ragged coverage makes root-to-tip measure sequence "
        "length rather than elapsed time.",
    )
    parser.add_argument(
        "--max-diversity",
        type=float,
        default=0.10,
        help="Maximum within-clade tip-to-tip distance (subs/site).",
    )
    parser.add_argument("--min-size", type=int, default=30)
    parser.add_argument(
        "--min-support",
        type=float,
        default=None,
        help="Warn about clades whose stem support is below this "
        "(UFBoot: 95 is the usual threshold). Ignored on a "
        "tree with no support values.",
    )
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--references", help="lineage_references.tsv, to label clades afterwards")
    parser.add_argument("--no-midpoint-root", action="store_true")
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-json")
    parser.add_argument("--out-dir", help="Write one tip list per clade here.")


@step("tree", "clades", "Cut the ML tree into clades and test each for temporal signal")
def clades(args: argparse.Namespace) -> int:
    """Cut clades from the ML tree, fit root-to-tip in each, write the summary table."""
    from Bio import Phylo

    # Bio.Phylo ships neither an explicit export of `read` nor annotations for it.
    tree = Phylo.read(args.tree, args.tree_format)  # type: ignore[attr-defined,no-untyped-call]

    if args.qc and args.min_ungapped > 0:
        _prune_short_tips(tree, args.qc, args.min_ungapped)

    if not args.no_midpoint_root:
        try:
            tree.root_at_midpoint()
        except Exception as exc:  # noqa: BLE001 -- Bio raises bare Exceptions; keep the given root
            _say(f"NOTE: midpoint rooting failed ({exc}); using the given root")

    refs = _load_references(args.references) if args.references else {}

    kept, small = cut(tree, args.max_diversity, args.min_size)
    _say(
        f"\n{len(kept)} clades at >= {args.min_size} tips, {len(small)} smaller groups set aside\n"
    )

    rows, reports = _assess_clades(tree, kept, refs, args)
    _write_clade_table(rows, args.out_tsv)
    if args.out_json:
        Path(args.out_json).write_text(json.dumps(reports, indent=2))
    _report_clades(rows, args.min_support)
    return 0

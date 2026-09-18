#!/usr/bin/env python3
"""
05_tree_summary.py — summarise an ML tree: major clades and their host composition.

Answers two questions you can't easily eyeball in FigTree once you have more than
a few dozen tips:

  1. What is the host composition of the dataset overall?
  2. How do wild carnivore sequences distribute across the major clades?

The tree is midpoint-rooted, then repeatedly split at its deepest divisions until
there are roughly `--n-clades` groups. Those groups approximate the major lineages.
This is a structural summary, NOT formal lineage assignment — for that you need
reference sequences with published lineage labels in the alignment.

Usage
-----
    python scripts/05_tree_summary.py --tree data/processed/H_ml.treefile

    # more or fewer groups
    python scripts/05_tree_summary.py --tree ... --n-clades 12

    # write a per-tip table with clade membership
    python scripts/05_tree_summary.py --tree ... --write-table

Tip labels are expected in the format written by 03_align_and_tree.py:
    ACCESSION|host_group|decimal_year
"""

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    from Bio import Phylo
except ImportError as exc:
    raise ImportError(
        "Biopython is required.\n"
        "  terminal: pip install biopython\n"
        "  notebook: %pip install biopython   (then restart the kernel)"
    ) from exc

WILD = {"wild_felid", "wild_canid", "mustelid", "procyonid",
        "pinniped", "ursid", "ailurid", "viverrid"}


def parse_label(name: str):
    """ACCESSION|host_group|decimal_year -> (accession, host_group, year|None)"""
    if not name:
        return ("", "unknown", None)
    parts = str(name).split("|")
    if len(parts) < 3:
        return (parts[0], "unparsed", None)
    try:
        year = float(parts[-1])
    except ValueError:
        year = None
    return (parts[0], parts[1], year)


def split_into_clades(tree, target: int, min_size: int):
    """Repeatedly split the largest group at its own deepest division."""
    groups = [tree.root]
    while len(groups) < target:
        splittable = [g for g in groups
                      if len(g.get_terminals()) >= max(min_size * 2, 4) and g.clades]
        if not splittable:
            break
        biggest = max(splittable, key=lambda c: len(c.get_terminals()))
        children = [c for c in biggest.clades if len(c.get_terminals()) >= 1]
        if len(children) < 2:
            break
        groups.remove(biggest)
        groups.extend(children)
    return sorted(groups, key=lambda c: -len(c.get_terminals()))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True, type=Path)
    ap.add_argument("--n-clades", type=int, default=8)
    ap.add_argument("--min-size", type=int, default=3)
    ap.add_argument("--no-midpoint-root", action="store_true",
                    help="Skip midpoint rooting (use if the tree is already rooted)")
    ap.add_argument("--write-table", action="store_true",
                    help="Write a per-tip table with clade membership")
    ap.add_argument("--outdir", type=Path, default=Path("data/processed"))
    args = ap.parse_args()

    if not args.tree.is_file():
        print(f"Tree not found: {args.tree}", file=sys.stderr)
        parent = args.tree.parent
        if parent.is_dir():
            found = sorted(p.name for p in parent.glob("*.tree*"))
            print(f"Tree files in {parent}: {found or '(none)'}", file=sys.stderr)
        print("\nBuild the tree first:\n"
              "    python scripts/03_align_and_tree.py "
              "--fasta data/processed/sequences_H.fasta "
              "--metadata data/processed/metadata_clean.tsv", file=sys.stderr)
        return 1

    tree = Phylo.read(str(args.tree), "newick")
    if not args.no_midpoint_root:
        try:
            tree.root_at_midpoint()
        except Exception as exc:                       # noqa: BLE001
            print(f"(midpoint rooting failed: {exc}; using the tree as-is)\n")

    tips = tree.get_terminals()
    parsed = [parse_label(t.name) for t in tips]

    # ---- overall composition -------------------------------------------
    print("=" * 64)
    print(f"  {len(tips)} tips in {args.tree.name}")
    print("=" * 64)

    hosts = Counter(h for _, h, _ in parsed)
    years = [y for _, _, y in parsed if y is not None]
    print("\nHOST COMPOSITION")
    for host, n in hosts.most_common():
        marker = "  <- wild" if host in WILD else ""
        print(f"  {host:<18} {n:>5}   {n/len(tips):>5.1%}{marker}")
    n_wild = sum(n for h, n in hosts.items() if h in WILD)
    print(f"\n  wild carnivores total: {n_wild} ({n_wild/len(tips):.1%})")
    if years:
        print(f"  sampling years: {min(years):.1f} - {max(years):.1f}")
    if "unparsed" in hosts:
        print(f"\n  WARNING: {hosts['unparsed']} tips didn't parse as "
              "ACCESSION|host_group|year — check the tip labels.")

    # ---- clade structure -------------------------------------------------
    groups = split_into_clades(tree, args.n_clades, args.min_size)
    print("\n" + "=" * 64)
    print(f"  MAJOR CLADES ({len(groups)} groups)")
    print("  Structural split of the tree, not published lineage assignment.")
    print("=" * 64)

    tip_to_clade = {}
    rows = []
    for i, clade in enumerate(groups, start=1):
        members = [parse_label(t.name) for t in clade.get_terminals()]
        if len(members) < args.min_size:
            label = f"clade_{i}_small"
        else:
            label = f"clade_{i}"
        comp = Counter(h for _, h, _ in members)
        yrs = [y for _, _, y in members if y is not None]
        wild_here = sum(n for h, n in comp.items() if h in WILD)
        support = getattr(clade, "confidence", None)

        print(f"\n  {label}  ({len(members)} tips"
              + (f", support {support}" if support is not None else "") + ")")
        if yrs:
            print(f"    years: {min(yrs):.1f} - {max(yrs):.1f}  (span {max(yrs)-min(yrs):.1f})")
        for host, n in comp.most_common():
            marker = "  <-" if host in WILD else ""
            print(f"    {host:<18} {n:>4}{marker}")
        if wild_here:
            print(f"    >> {wild_here} wild carnivore tips in this clade")

        for acc, host, year in members:
            tip_to_clade[acc] = label
            rows.append({"accession": acc, "host_group": host,
                         "decimal_year": year, "clade": label})

    # ---- where the wild sequences are ------------------------------------
    print("\n" + "=" * 64)
    print("  WILD CARNIVORE TIPS BY CLADE")
    print("=" * 64)
    by_clade = defaultdict(Counter)
    for r in rows:
        if r["host_group"] in WILD:
            by_clade[r["clade"]][r["host_group"]] += 1
    if not by_clade:
        print("\n  No wild carnivore tips found. Check host_group values.")
    for clade in sorted(by_clade, key=lambda c: -sum(by_clade[c].values())):
        total = sum(by_clade[clade].values())
        detail = ", ".join(f"{h} {n}" for h, n in by_clade[clade].most_common())
        print(f"  {clade:<16} {total:>4} wild   ({detail})")

    print("\n  The clade with the most wild carnivore tips and a decent time span")
    print("  is the natural candidate for a focused within-lineage analysis.")

    if args.write_table:
        try:
            import pandas as pd
            args.outdir.mkdir(parents=True, exist_ok=True)
            out = args.outdir / (args.tree.stem + "_clades.tsv")
            pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
            print(f"\nper-tip clade table -> {out}")
            print("Use it to pull a subset alignment for the focused analysis.")
        except ImportError:
            print("\n(pandas not installed — skipping the table)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

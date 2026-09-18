#!/usr/bin/env python3
"""
06_subsample.py — reduce an alignment to a BEAST-tractable size, or pull one clade.

Two modes:

  GLOBAL   stratified subsample across the whole dataset, preserving clade,
           host-group and temporal structure
  CLADE    extract all tips belonging to one clade, optionally subsampled

Why stratify: a random subsample of a dataset that is 52% domestic dog will be 52%
domestic dog, discarding most of the rare host groups your host-jump analysis depends
on. Stratifying by clade x host x time bin keeps the informative tips and drops the
redundant ones.

Sampling within a stratum is deliberate, not random:
  1. sequences with day-precision dates before year-only dates
  2. longer ungapped sequences before shorter ones
  3. one representative per identical-sequence set

Usage
-----
    # global stratified subsample to ~400 tips
    python scripts/06_subsample.py --mode global --target 400 \\
        --aln data/processed/H_aligned.fasta \\
        --clades data/processed/H_ml_clades.tsv \\
        --metadata data/processed/metadata_clean.tsv

    # one clade, everything in it
    python scripts/06_subsample.py --mode clade --clade clade_3 \\
        --aln data/processed/H_aligned.fasta \\
        --clades data/processed/H_ml_clades.tsv \\
        --metadata data/processed/metadata_clean.tsv

    # one clade, capped
    python scripts/06_subsample.py --mode clade --clade clade_2 --target 250 ...

Outputs a FASTA, a metadata table for the subset, and a record of what was kept and
dropped. Use --seed to make the selection reproducible (it is by default).
"""

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    import pandas as pd
except ImportError as exc:
    raise ImportError(
        "pandas is required.\n"
        "  terminal: pip install pandas\n"
        "  notebook: %pip install pandas   (then restart the kernel)"
    ) from exc

WILD = {"wild_felid", "wild_canid", "mustelid", "procyonid",
        "pinniped", "ursid", "ailurid", "viverrid"}
GAP = set("-.")


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, chunks = {}, None, []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def write_fasta(path: Path, seqs: dict[str, str], width: int = 60) -> None:
    with path.open("w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def strip_all_gap_columns(seqs: dict[str, str]) -> tuple[dict[str, str], int]:
    """Drop columns that are gaps in every remaining sequence.

    A subset of an alignment keeps the full column count of the original, so any
    column that only ever held bases belonging to dropped sequences survives as
    an all-gap column. These are not harmless. Anything that treats an alignment
    column span as a coordinate — a CDS annotation, a codon partition, a site
    index quoted in a figure — is thrown off by them, and the offset is invisible
    because the columns look like ordinary gaps.
    """
    if not seqs:
        return seqs, 0
    L = len(next(iter(seqs.values())))
    keep = [i for i in range(L) if any(s[i] not in GAP for s in seqs.values())]
    if len(keep) == L:
        return seqs, 0
    return {n: "".join(s[i] for i in keep) for n, s in seqs.items()}, L - len(keep)


def parse_label(label: str):
    parts = str(label).split("|")
    if len(parts) < 3:
        return (parts[0], "unparsed", None)
    try:
        year = float(parts[-1])
    except ValueError:
        year = None
    return (parts[0], parts[1], year)


def build_records(seqs, clades_df, meta_df):
    """One row per aligned sequence, with clade, host, year, quality fields."""
    clade_of = dict(zip(clades_df["accession"], clades_df["clade"])) if clades_df is not None else {}
    prec_of, len_of = {}, {}
    if meta_df is not None:
        base = meta_df["accession"].astype(str).str.split(".").str[0]
        prec_of = dict(zip(base, meta_df.get("date_precision", pd.Series(["year"] * len(meta_df)))))
        len_of = dict(zip(base, meta_df.get("H_length", pd.Series([0] * len(meta_df)))))

    rows = []
    for label, seq in seqs.items():
        acc, host, year = parse_label(label)
        ungapped = "".join(c for c in seq if c not in GAP)
        rows.append({
            "label": label,
            "accession": acc,
            "host_group": host,
            "decimal_year": year,
            "clade": clade_of.get(acc, "unassigned"),
            "date_precision": prec_of.get(acc, "unknown"),
            "ungapped_len": len(ungapped),
            "ambig": sum(1 for c in ungapped.upper() if c in "NX?"),
            "seq_key": ungapped.upper(),
        })
    return pd.DataFrame(rows)


PRECISION_RANK = {"day": 0, "month": 1, "range": 2, "year": 3, "unknown": 4, "none": 5}


def rank_within_stratum(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["prec_rank"] = d["date_precision"].map(lambda p: PRECISION_RANK.get(p, 4))
    d["ambig_frac"] = d["ambig"] / d["ungapped_len"].clip(lower=1)
    return d.sort_values(["prec_rank", "ambig_frac", "ungapped_len"],
                         ascending=[True, True, False])


def drop_duplicates(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Keep one representative per identical sequence, preferring the best-dated."""
    ranked = rank_within_stratum(df)
    kept = ranked.drop_duplicates(subset="seq_key", keep="first")
    return kept, len(df) - len(kept)


def stratified_sample(df: pd.DataFrame, target: int, year_bin: int, seed: int,
                      min_per_stratum: int = 1) -> pd.DataFrame:
    """
    Allocate the target across clade x host x time-bin strata.

    Every non-empty stratum gets at least min_per_stratum, then the remaining budget
    is distributed proportional to the square root of stratum size — which keeps large
    strata dominant without letting them crowd out the rare host groups entirely.
    """
    d = df.copy()
    d["year_bin"] = d["decimal_year"].apply(
        lambda y: int(y // year_bin) * year_bin if pd.notna(y) else -1)
    d["stratum"] = d["clade"] + "|" + d["host_group"] + "|" + d["year_bin"].astype(str)

    strata = {k: g for k, g in d.groupby("stratum")}
    if len(strata) >= target:
        # more strata than budget: one per stratum, prioritising wild-host strata
        picks = []
        order = sorted(strata.items(),
                       key=lambda kv: (kv[1]["host_group"].iloc[0] not in WILD, -len(kv[1])))
        for _, g in order[:target]:
            picks.append(rank_within_stratum(g).head(1))
        return pd.concat(picks)

    allocation = {k: min_per_stratum for k in strata}
    remaining = target - sum(allocation.values())
    if remaining > 0:
        weights = {k: max(len(g) - min_per_stratum, 0) ** 0.5 for k, g in strata.items()}
        total_w = sum(weights.values())
        if total_w > 0:
            for k, w in weights.items():
                allocation[k] += int(remaining * w / total_w)

    # hand out any rounding leftovers to the largest under-sampled strata
    used = sum(min(allocation[k], len(g)) for k, g in strata.items())
    leftover = target - used
    if leftover > 0:
        for k, g in sorted(strata.items(), key=lambda kv: -len(kv[1])):
            room = len(g) - allocation[k]
            if room > 0:
                take = min(room, leftover)
                allocation[k] += take
                leftover -= take
            if leftover <= 0:
                break

    picks = [rank_within_stratum(g).head(min(allocation[k], len(g)))
             for k, g in strata.items()]
    return pd.concat(picks)


def summarise(df: pd.DataFrame, title: str) -> None:
    print(f"\n{title}: {len(df)} sequences")
    if df.empty:
        return
    hosts = Counter(df["host_group"])
    n_wild = sum(n for h, n in hosts.items() if h in WILD)
    for host, n in hosts.most_common():
        print(f"    {host:<16} {n:>5}{'  <- wild' if host in WILD else ''}")
    print(f"    {'wild total':<16} {n_wild:>5}  ({n_wild/len(df):.1%})")
    yrs = df["decimal_year"].dropna()
    if len(yrs):
        print(f"    years {yrs.min():.1f} - {yrs.max():.1f}")
    clades = Counter(df["clade"])
    if len(clades) > 1:
        print("    clades: " + ", ".join(f"{c}={n}" for c, n in clades.most_common(8)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["global", "clade"], required=True)
    ap.add_argument("--aln", required=True, type=Path)
    ap.add_argument("--clades", type=Path, help="H_ml_clades.tsv from 05_tree_summary.py")
    ap.add_argument("--metadata", type=Path)
    ap.add_argument("--clade", help="Clade name, required for --mode clade")
    ap.add_argument("--target", type=int, default=None,
                    help="Target sequence count (global default 400; clade default = keep all)")
    ap.add_argument("--year-bin", type=int, default=5, help="Time bin width in years")
    ap.add_argument("--seed", type=int, default=20260819)
    ap.add_argument("--keep-duplicates", action="store_true")
    ap.add_argument("--keep-all-gap-columns", action="store_true",
                    help="Keep columns that are all-gap in the subset. Off by default: "
                         "such columns shift every downstream column coordinate and are "
                         "invisible when you look at the alignment.")
    ap.add_argument("--outdir", type=Path, default=Path("data/processed"))
    ap.add_argument("--prefix", default=None)
    args = ap.parse_args()

    if args.mode == "clade" and not args.clade:
        print("--mode clade requires --clade (e.g. --clade clade_3)", file=sys.stderr)
        return 1
    if not args.aln.is_file():
        print(f"Alignment not found: {args.aln}", file=sys.stderr)
        return 1

    seqs = read_fasta(args.aln)
    clades_df = pd.read_csv(args.clades, sep="\t") if args.clades and args.clades.is_file() else None
    meta_df = pd.read_csv(args.metadata, sep="\t") if args.metadata and args.metadata.is_file() else None
    if clades_df is None and args.mode == "clade":
        print("--mode clade needs --clades (run 05_tree_summary.py --write-table)", file=sys.stderr)
        return 1

    df = build_records(seqs, clades_df, meta_df)
    print(f"Loaded {len(df)} aligned sequences")
    if (df["clade"] == "unassigned").any():
        n = int((df["clade"] == "unassigned").sum())
        print(f"  note: {n} sequences had no clade assignment "
              "(present in the alignment but not the tree)")

    # ---- select pool ----------------------------------------------------
    if args.mode == "clade":
        pool = df[df["clade"] == args.clade].copy()
        if pool.empty:
            print(f"No sequences in {args.clade}. Available: "
                  f"{sorted(df['clade'].unique())}", file=sys.stderr)
            return 1
        title = args.clade
    else:
        pool = df.copy()
        title = "global"

    summarise(pool, f"Pool ({title})")

    # ---- deduplicate ------------------------------------------------------
    if not args.keep_duplicates:
        pool, n_dropped = drop_duplicates(pool)
        if n_dropped:
            print(f"\nRemoved {n_dropped} exact duplicate sequences "
                  "(kept the best-dated representative of each)")

    # ---- subsample --------------------------------------------------------
    target = args.target
    if args.mode == "global" and target is None:
        target = 400
    if target is not None and len(pool) > target:
        pool = stratified_sample(pool, target, args.year_bin, args.seed)
        method = f"stratified subsample to {target}"
    else:
        method = "all sequences retained"
    print(f"\n{method}")

    summarise(pool, "FINAL SUBSET")

    # sanity checks the analysis actually depends on
    hosts = Counter(pool["host_group"])
    thin = {h: n for h, n in hosts.items() if h in WILD and n < 5}
    if thin:
        print(f"\n  WARNING: wild host groups with <5 sequences: {thin}")
        print("  Discrete trait analysis will be unreliable for these. Consider merging")
        print("  them into a broader category, or raising --target.")
    yrs = pool["decimal_year"].dropna()
    if len(yrs) and (yrs.max() - yrs.min()) < 10:
        print(f"\n  WARNING: temporal span is only {yrs.max()-yrs.min():.1f} years — "
              "weak clock signal likely.")

    # ---- write -------------------------------------------------------------
    args.outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or (f"H_{args.clade}" if args.mode == "clade" else f"H_global{len(pool)}")

    keep_labels = list(pool["label"])
    subset = {lbl: seqs[lbl] for lbl in keep_labels}
    if not args.keep_all_gap_columns:
        subset, n_stripped = strip_all_gap_columns(subset)
        if n_stripped:
            L_before = len(next(iter(seqs.values())))
            print(f"\nStripped {n_stripped} all-gap columns "
                  f"({L_before} -> {L_before - n_stripped}); they held only bases from "
                  "sequences this subset dropped")
    out_fasta = args.outdir / f"{prefix}.fasta"
    write_fasta(out_fasta, subset)

    out_meta = args.outdir / f"{prefix}_metadata.tsv"
    pool.drop(columns=["seq_key"]).to_csv(out_meta, sep="\t", index=False)

    out_dates = args.outdir / f"{prefix}_dates.tsv"
    with out_dates.open("w") as fh:
        fh.write("taxon\tdate\n")
        for lbl in keep_labels:
            fh.write(f"{lbl}\t{lbl.rsplit('|', 1)[1]}\n")

    print(f"\nwrote {out_fasta}")
    print(f"wrote {out_meta}")
    print(f"wrote {out_dates}")
    print("\nNOTE: this FASTA is already aligned and already has BEAST-style tip labels.")
    print("Do NOT feed it back to 03_align_and_tree.py — that script expects unaligned")
    print("sequences with bare accessions and would drop everything. Run IQ-TREE directly:")
    print(f"\n    iqtree2 -s {out_fasta} -m MFP -B 1000 --alrt 1000 -T AUTO \\")
    print(f"        --prefix {args.outdir / prefix}_ml -redo")
    print("\nAll-gap columns have already been stripped (--keep-all-gap-columns to keep them).")
    print("Re-check TempEst — a subset with better temporal balance often has")
    print("markedly stronger clock signal than the full tree did.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
04_alignment_qc.py — inspect an alignment before spending compute on a tree.

This is the "look at your alignment" step. On a desktop you'd open it in AliView
or Jalview; on a hosted kernel you probably can't, so this reports the same things
numerically and as a coverage plot.

What it checks
--------------
  * per-sequence gap and ambiguity fractions, with outliers flagged
  * per-column coverage, and suggested trim points for ragged ends
  * identity of each sequence to the consensus (catches misaligned or wrong-gene entries)
  * in-frame stop codons (catches frameshifts and wrong reading frames)
  * exact duplicate sequences
  * all-gap columns, which shift every downstream column coordinate
  * tip labels against the metadata table, so a stale host group cannot survive
  * host-group composition of what survives

Usage
-----
    python scripts/04_alignment_qc.py --aln data/processed/H_aligned.fasta \\
                                      --metadata data/processed/metadata_clean.tsv

    # after deciding trim points from the coverage plot
    python scripts/04_alignment_qc.py --aln ... --trim-to 42 1801 --write-trimmed

Nothing here decides anything for you. It tells you which sequences are worth
looking at by eye and which columns are mostly gaps.
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

try:
    import pandas as pd
except ImportError as exc:
    raise ImportError(
        "pandas is required.\n"
        "  terminal: pip install pandas\n"
        "  notebook: %pip install pandas   (then restart the kernel)"
    ) from exc

AMBIG = set("NnXx?")
GAP = set("-.")
STOPS = {"TAA", "TAG", "TGA"}


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


def consensus(seqs: dict[str, str]) -> str:
    if not seqs:
        return ""
    length = len(next(iter(seqs.values())))
    out = []
    for i in range(length):
        col = Counter(s[i].upper() for s in seqs.values()
                      if s[i] not in GAP and s[i] not in AMBIG)
        out.append(col.most_common(1)[0][0] if col else "-")
    return "".join(out)


def identity_to(seq: str, ref: str) -> float:
    """Fraction of comparable positions that match. Gaps in either are skipped."""
    same = comp = 0
    for a, b in zip(seq.upper(), ref):
        if a in GAP or b in GAP or a in AMBIG:
            continue
        comp += 1
        same += (a == b)
    return same / comp if comp else 0.0


def count_stops(seq: str, frame: int = 0) -> int:
    ungapped = "".join(c for c in seq.upper() if c not in GAP)[frame:]
    codons = [ungapped[i:i + 3] for i in range(0, len(ungapped) - 2, 3)]
    return sum(1 for c in codons[:-1] if c in STOPS)   # ignore a terminal stop


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True, type=Path)
    ap.add_argument("--metadata", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=Path("data/processed"))
    ap.add_argument("--gap-threshold", type=float, default=0.5,
                    help="Flag sequences more than this fraction gaps (default %(default)s)")
    ap.add_argument("--identity-threshold", type=float, default=0.80,
                    help="Flag sequences below this identity to consensus (default %(default)s)")
    ap.add_argument("--coverage-threshold", type=float, default=0.5,
                    help="Columns with less than this coverage count as ragged (default %(default)s)")
    ap.add_argument("--trim-to", nargs=2, type=int, metavar=("START", "END"),
                    help="Trim to these 1-based column positions (inclusive)")
    ap.add_argument("--write-trimmed", action="store_true")
    ap.add_argument("--strip-all-gap", action="store_true",
                    help="Remove columns that are gaps in every sequence")
    ap.add_argument("--write-stripped", action="store_true",
                    help="Save the result of --strip-all-gap")
    ap.add_argument("--plot", type=Path, default=None,
                    help="Write a coverage plot PNG (default: alongside the alignment)")
    args = ap.parse_args()

    if not args.aln.is_file():
        print(f"Alignment not found: {args.aln}\n", file=sys.stderr)
        parent = args.aln.parent
        if parent.is_dir():
            found = sorted(p.name for p in parent.glob("*.fasta"))
            print(f"FASTA files in {parent}: {found or '(none)'}", file=sys.stderr)
        else:
            print(f"Directory {parent} doesn't exist either.", file=sys.stderr)
        print("\nRun the alignment step first:\n"
              "    python scripts/03_align_and_tree.py \\\n"
              "        --fasta data/processed/sequences_H.fasta \\\n"
              "        --metadata data/processed/metadata_clean.tsv \\\n"
              "        --stop-after align", file=sys.stderr)
        return 1

    seqs = read_fasta(args.aln)
    if not seqs:
        print(f"No sequences read from {args.aln}", file=sys.stderr)
        return 1

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        print(f"ERROR: sequences have different lengths {sorted(lengths)[:5]} — "
              "this file is not aligned.", file=sys.stderr)
        return 1
    L = lengths.pop()
    N = len(seqs)
    print(f"{N} sequences, {L} alignment columns\n")

    # ---- per-column coverage -------------------------------------------
    cov = [sum(1 for s in seqs.values() if s[i] not in GAP) / N for i in range(L)]
    good = [i for i, c in enumerate(cov) if c >= args.coverage_threshold]
    if good:
        first, last = good[0], good[-1]
        print(f"Columns at >= {args.coverage_threshold:.0%} coverage: "
              f"{first + 1} to {last + 1} (1-based)")
        if first > 0 or last < L - 1:
            print(f"  ragged: {first} leading and {L - 1 - last} trailing columns are sparse")
            print(f"  suggested:  --trim-to {first + 1} {last + 1} --write-trimmed")
    print(f"Mean column coverage: {sum(cov) / L:.1%}")

    # ---- all-gap columns -------------------------------------------------
    # A column that is a gap in every sequence carries no data but still
    # occupies a position. Anything that treats a column span as a coordinate —
    # a CDS annotation, a codon partition, a site index in a figure caption — is
    # silently shifted by it. These usually appear when an alignment is subset
    # without re-stripping, so the columns held bases from sequences now gone.
    empty = [i for i, c in enumerate(cov) if c == 0]
    if empty:
        runs, start, prev = [], empty[0], empty[0]
        for i in empty[1:]:
            if i != prev + 1:
                runs.append((start, prev)); start = i
            prev = i
        runs.append((start, prev))
        interior = [(a, b) for a, b in runs if a > (good[0] if good else 0)
                    and b < (good[-1] if good else L - 1)]
        print(f"\n{len(empty)} all-gap columns in {len(runs)} run(s)"
              + (f", {sum(b - a + 1 for a, b in interior)} of them interior"
                 if interior else ""))
        for a, b in runs[:8]:
            print(f"    columns {a + 1}-{b + 1}  ({b - a + 1} nt)")
        if len(runs) > 8:
            print(f"    ... and {len(runs) - 8} more")
        print("  These shift every column coordinate after them. Remove with")
        print("  --strip-all-gap --write-stripped, or re-run 06_subsample.py, which")
        print("  now strips them when it writes a subset.")
    else:
        print("No all-gap columns.")
    print()

    # ---- per-sequence stats --------------------------------------------
    cons = consensus(seqs)
    rows = []
    for name, s in seqs.items():
        ung = "".join(c for c in s if c not in GAP)
        rows.append({
            "name": name,
            "ungapped_len": len(ung),
            "gap_frac": sum(1 for c in s if c in GAP) / L,
            "ambig_frac": (sum(1 for c in ung if c in AMBIG) / len(ung)) if ung else 1.0,
            "identity_to_consensus": identity_to(s, cons),
        })
    df = pd.DataFrame(rows)

    # ---- flags ----------------------------------------------------------
    df["flag"] = ""
    def flag(mask, label):
        df.loc[mask, "flag"] = (df.loc[mask, "flag"] + ";" + label).str.lstrip(";")
        return int(mask.sum())

    n_gap = flag(df["gap_frac"] > args.gap_threshold, "mostly_gaps")
    n_amb = flag(df["ambig_frac"] > 0.05, "ambiguous>5%")
    n_id  = flag(df["identity_to_consensus"] < args.identity_threshold, "low_identity")
    # Determine the alignment's reading frame from the whole dataset: the frame
    # with the fewest total internal stops. Then flag sequences that still have
    # stops in THAT frame — those are the genuinely suspect ones (frameshift,
    # wrong strand, pseudogene, sequencing error). Testing each sequence against
    # its own best frame is meaningless, because frames 1 and 2 of any real
    # coding sequence are full of stops.
    frame_totals = {f: sum(count_stops(s, f) for s in seqs.values()) for f in (0, 1, 2)}
    aln_frame = min(frame_totals, key=frame_totals.get)
    df["stops_in_aln_frame"] = [count_stops(seqs[n], aln_frame) for n in df["name"]]

    n_st = flag(df["stops_in_aln_frame"] > 0, "internal_stop")

    print("Flagged sequences")
    print(f"  mostly gaps (>{args.gap_threshold:.0%})        : {n_gap}")
    print(f"  >5% ambiguous bases                : {n_amb}")
    print(f"  identity to consensus <{args.identity_threshold:.0%}      : {n_id}")
    print(f"  internal stops in frame {aln_frame}          : {n_st}")
    print()
    print(f"  Alignment reading frame inferred as {aln_frame} "
          f"(total internal stops per frame: {frame_totals}).")
    print("  If frame 0 isn't the winner, the alignment starts mid-codon — normal for")
    print("  GenBank records, but use this offset when setting up codon partitions.")
    print("  Sequences with stops in the inferred frame are worth checking: frameshift,")
    print("  wrong strand, or sequencing error.")
    print()

    flagged = df[df["flag"] != ""].sort_values("identity_to_consensus")
    if not flagged.empty:
        print("Worth looking at by eye (worst first):")
        print(flagged[["name", "ungapped_len", "gap_frac",
                       "identity_to_consensus", "flag"]].head(20).to_string(index=False))
        print()

    # ---- duplicates ------------------------------------------------------
    by_seq: dict[str, list[str]] = {}
    for name, s in seqs.items():
        by_seq.setdefault("".join(c for c in s.upper() if c not in GAP), []).append(name)
    dups = {k: v for k, v in by_seq.items() if len(v) > 1}
    if dups:
        print(f"{len(dups)} sets of identical sequences ({sum(len(v) for v in dups.values())} seqs):")
        for names in list(dups.values())[:5]:
            print("   ", ", ".join(names[:4]) + (" ..." if len(names) > 4 else ""))
        print("  Identical sequences are legitimate data, not errors — but they")
        print("  contribute no phylogenetic signal and slow BEAST down. Consider")
        print("  keeping one per host/year/country combination.\n")

    # ---- labels vs metadata ----------------------------------------------
    # Tip labels encode host group and date, so they are a copy of the metadata
    # that can go stale independently of it. A label that disagrees with the
    # table is a real error: it will silently propagate into the tree, the BEAST
    # run and every figure, and nothing downstream re-derives it.
    if args.metadata and args.metadata.is_file():
        meta_df = pd.read_csv(args.metadata, sep="\t", dtype=str)
        if "accession" in meta_df.columns:
            by_acc = {str(a).split(".")[0]: r
                      for a, r in zip(meta_df["accession"], meta_df.to_dict("records"))}
            unknown, host_bad, date_bad = [], [], []
            for name in seqs:
                parts = name.split("|")
                if len(parts) < 3:
                    continue
                acc, host, year = parts[0].split(".")[0], parts[1], parts[-1]
                row = by_acc.get(acc)
                if row is None:
                    unknown.append(name); continue
                m_host = (row.get("host_group") or "").strip()
                if m_host and m_host != host:
                    host_bad.append((name, host, m_host))
                m_year = (row.get("decimal_year") or "").strip()
                try:
                    if m_year and abs(float(m_year) - float(year)) > 0.01:
                        date_bad.append((name, year, m_year))
                except ValueError:
                    pass

            print("Labels vs metadata")
            print(f"  accession not in metadata          : {len(unknown)}")
            print(f"  host group disagrees with metadata : {len(host_bad)}")
            print(f"  date disagrees with metadata       : {len(date_bad)}")
            for name, in_label, in_meta in (host_bad + date_bad)[:10]:
                print(f"    {name}")
                print(f"      label says {in_label!r}, metadata says {in_meta!r}")
            for name in unknown[:5]:
                print(f"    {name}  (no metadata row)")
            if host_bad or date_bad:
                print("  Fix the label, not the metadata, unless the table is the stale one.")
                print("  Relabelling changes tip names, so anything built from these labels")
                print("  — trees, BEAST XML, exported datasets — has to be regenerated.")
            print()

    # ---- host composition ------------------------------------------------
    if args.metadata and args.metadata.is_file():
        groups = Counter(n.split("|")[1] for n in seqs if "|" in n)
        if groups:
            print("Host composition of the alignment:")
            for g, c in groups.most_common():
                print(f"    {g:<16} {c:>5}")
            thin = [g for g, c in groups.items() if c < 5]
            if thin:
                print(f"  WARNING: fewer than 5 sequences in: {', '.join(thin)}")
            print()

    # ---- outputs ---------------------------------------------------------
    args.outdir.mkdir(parents=True, exist_ok=True)
    qc_path = args.outdir / (args.aln.stem + "_qc.tsv")
    df.sort_values("identity_to_consensus").to_csv(qc_path, sep="\t", index=False)
    print(f"per-sequence QC table -> {qc_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(10, 5.5))
        ax[0].fill_between(range(1, L + 1), cov, color="#4B3A79", lw=0)
        ax[0].axhline(args.coverage_threshold, color="#C0697F", ls="--", lw=1)
        ax[0].set_title("Per-column coverage (fraction of sequences with a base)")
        ax[0].set_xlabel("alignment column"); ax[0].set_ylim(0, 1)
        ax[1].hist(df["identity_to_consensus"], bins=40, color="#6E5CA0")
        ax[1].axvline(args.identity_threshold, color="#C0697F", ls="--", lw=1)
        ax[1].set_title("Identity to consensus")
        plt.tight_layout()
        png = args.plot or args.outdir / (args.aln.stem + "_qc.png")
        plt.savefig(png, dpi=130)
        print(f"coverage plot          -> {png}")
    except ImportError:
        print("(matplotlib not installed — skipping the plot)")

    # ---- optional strip --------------------------------------------------
    if args.strip_all_gap:
        keep = [i for i in range(L) if cov[i] > 0]
        if len(keep) == L:
            print("\nNo all-gap columns to strip.")
        else:
            stripped = {n: "".join(s[i] for i in keep) for n, s in seqs.items()}
            print(f"\nStripped {L - len(keep)} all-gap columns ({L} -> {len(keep)})")
            if args.write_stripped:
                out = args.outdir / (args.aln.stem + "_stripped.fasta")
                write_fasta(out, stripped)
                print(f"wrote {out}")
            else:
                print("(add --write-stripped to save it)")

    # ---- optional trim ---------------------------------------------------
    if args.trim_to:
        s0, s1 = args.trim_to
        trimmed = {n: s[s0 - 1:s1] for n, s in seqs.items()}
        empty = [n for n, s in trimmed.items() if set(s) <= GAP]
        if empty:
            print(f"\n{len(empty)} sequences become all-gap after trimming and were dropped:")
            print("   ", ", ".join(empty[:6]))
            trimmed = {n: s for n, s in trimmed.items() if n not in empty}
        print(f"\nTrimmed to columns {s0}-{s1}: {len(trimmed)} sequences, {s1 - s0 + 1} columns")
        if args.write_trimmed:
            out = args.outdir / (args.aln.stem + "_trimmed.fasta")
            write_fasta(out, trimmed)
            print(f"wrote {out}")
        else:
            print("(add --write-trimmed to save it)")

    print("\nNext: if the alignment looks sound, build the tree:")
    print("    python scripts/03_align_and_tree.py --fasta data/processed/sequences_H.fasta \\")
    print("        --metadata data/processed/metadata_clean.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

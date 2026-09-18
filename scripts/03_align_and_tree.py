#!/usr/bin/env python3
"""
03_align_and_tree.py — alignment, tip relabelling, ML phylogeny, TempEst prep.

Runs the Phase 3 pipeline end to end:
  1. Relabel FASTA headers with host group and decimal year (BEAST/TempEst need this)
  2. Align with MAFFT
  3. Optionally trim poorly aligned columns with trimAl
  4. Build an ML tree in IQ-TREE with ModelFinder + UFBoot + SH-aLRT
  5. Write a TempEst-ready date file for the root-to-tip regression

Usage
-----
    python scripts/03_align_and_tree.py \
        --fasta data/processed/sequences_H.fasta \
        --metadata data/processed/metadata_clean.tsv

    # relabel + align only, skip the tree (fast, for eyeballing the alignment)
    python scripts/03_align_and_tree.py --fasta ... --metadata ... --stop-after align

Requires mafft and iqtree on PATH (see environment.yml). trimAl optional.

Tip label format
----------------
    ACCESSION|host_group|decimal_year          e.g.  KC1|wild_felid|2010.615

BEAST and TempEst both parse the trailing field as the sampling date, so the
year must be LAST and the separator consistent. Do not reorder these fields.

IMPORTANT — before you trust any tree from this script
------------------------------------------------------
Run it once on a published CDV dataset and confirm you recover the known lineage
structure (America-1/2, Asia-1/2/3/4, Europe wildlife, Arctic-like, etc.). If you
cannot reproduce a known result, you cannot trust a new one. This is the step
that replaces a supervisor.
"""

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError as exc:                       # noqa: F401
    # Raise rather than sys.exit: this module gets imported by the notebook,
    # and SystemExit inside a notebook produces a confusing traceback.
    raise ImportError(
        "pandas is required.\n"
        "  terminal: pip install pandas\n"
        "  notebook: %pip install pandas   (then restart the kernel)"
    ) from exc


def have(tool: str) -> bool:
    return shutil.which(tool) is not None


def run(cmd: list[str], log_path: Path | None = None) -> None:
    logging.info("$ %s", " ".join(cmd))
    if log_path:
        with log_path.open("w") as fh:
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    else:
        proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise SystemExit(f"command failed (exit {proc.returncode}): {' '.join(cmd)}")


def read_fasta(path: Path) -> dict[str, str]:
    seqs, acc, chunks = {}, None, []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(chunks)
                acc, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if acc is not None:
        seqs[acc] = "".join(chunks)
    return seqs


def write_fasta(path: Path, seqs: dict[str, str], width: int = 60) -> None:
    with path.open("w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def relabel(seqs: dict[str, str], meta: pd.DataFrame) -> tuple[dict[str, str], list[str]]:
    """Rename ACCESSION -> ACCESSION|host_group|decimal_year. Returns (seqs, dropped)."""
    lookup = meta.set_index("accession")[["host_group", "decimal_year"]].to_dict("index")
    out, dropped, seen = {}, [], set()
    for acc, seq in seqs.items():
        info = lookup.get(acc)
        if info is None or pd.isna(info["decimal_year"]):
            dropped.append(acc)
            continue
        base = acc.split(".")[0]                     # strip GenBank version suffix
        label = f"{base}|{info['host_group']}|{info['decimal_year']:.3f}"
        if label in seen:                            # defensive; should not happen
            label = f"{base}_dup|{info['host_group']}|{info['decimal_year']:.3f}"
        seen.add(label)
        out[label] = seq
    return out, dropped


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--outdir", default=Path("data/processed"), type=Path)
    ap.add_argument("--prefix", default="H", help="Output file prefix (default: %(default)s)")
    ap.add_argument("--threads", default="AUTO", help="Threads for MAFFT/IQ-TREE")
    ap.add_argument("--mafft-mode", default="auto", choices=["auto", "linsi", "ginsi"],
                    help="linsi is more accurate but slow above ~500 sequences")
    ap.add_argument("--trim", action="store_true", help="Run trimAl -automated1 after alignment")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--stop-after", choices=["relabel", "align", "tree"], default="tree")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s  %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    args.outdir.mkdir(parents=True, exist_ok=True)

    # --- 1. relabel -------------------------------------------------------
    seqs = read_fasta(args.fasta)
    meta = pd.read_csv(args.metadata, sep="\t")
    logging.info("%d sequences, %d metadata rows", len(seqs), len(meta))

    labelled, dropped = relabel(seqs, meta)
    if dropped:
        logging.warning("%d sequences had no usable metadata and were dropped: %s",
                        len(dropped), ", ".join(dropped[:8]) + (" ..." if len(dropped) > 8 else ""))
    if not labelled:
        logging.error("nothing left after relabelling — check that accessions match "
                      "between the FASTA and metadata_clean.tsv")
        return 1

    labelled_path = args.outdir / f"{args.prefix}_labelled.fasta"
    write_fasta(labelled_path, labelled)
    logging.info("wrote %s (%d sequences)", labelled_path, len(labelled))

    # TempEst / BEAST date file
    dates_path = args.outdir / f"{args.prefix}_dates.tsv"
    with dates_path.open("w") as fh:
        fh.write("taxon\tdate\n")
        for label in labelled:
            fh.write(f"{label}\t{label.rsplit('|', 1)[1]}\n")
    logging.info("wrote %s", dates_path)

    # composition summary — sanity check before spending compute
    counts = {}
    for label in labelled:
        counts[label.split("|")[1]] = counts.get(label.split("|")[1], 0) + 1
    logging.info("host composition of the alignment set:")
    for grp, cnt in sorted(counts.items(), key=lambda kv: -kv[1]):
        logging.info("    %-16s %5d", grp, cnt)
    thin = [g for g, c in counts.items() if c < 5]
    if thin:
        logging.warning("host groups with <5 sequences: %s", ", ".join(thin))
        logging.warning("  discrete trait analysis will be unreliable for these. "
                        "Consider merging them into a broader category or excluding them.")

    if args.stop_after == "relabel":
        return 0

    # --- 2. align ---------------------------------------------------------
    if not have("mafft"):
        logging.error("mafft not found on PATH. conda install -c bioconda mafft")
        return 1

    aln_path = args.outdir / f"{args.prefix}_aligned.fasta"
    mafft_cmd = ["mafft"]
    if args.mafft_mode == "linsi":
        mafft_cmd += ["--localpair", "--maxiterate", "1000"]
    elif args.mafft_mode == "ginsi":
        mafft_cmd += ["--globalpair", "--maxiterate", "1000"]
    else:
        mafft_cmd += ["--auto"]
    mafft_cmd += ["--adjustdirection", "--reorder",
                  "--thread", "-1" if args.threads == "AUTO" else str(args.threads),
                  str(labelled_path)]

    logging.info("aligning with MAFFT (%s)...", args.mafft_mode)
    with aln_path.open("w") as fh:
        proc = subprocess.run(mafft_cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        logging.error(proc.stderr[-2000:])
        return 1
    logging.info("wrote %s", aln_path)

    # --adjustdirection prefixes reverse-complemented seqs with _R_; strip it so
    # tip labels still match the date file.
    aln = read_fasta(aln_path)
    n_rev = sum(1 for k in aln if k.startswith("_R_"))
    if n_rev:
        logging.warning("%d sequences were reverse-complemented by MAFFT — "
                        "worth checking those records", n_rev)
        aln = {(k[3:] if k.startswith("_R_") else k): v for k, v in aln.items()}
        write_fasta(aln_path, aln)

    if args.trim:
        if not have("trimal"):
            logging.warning("trimal not found — skipping trimming")
        else:
            trimmed = args.outdir / f"{args.prefix}_aligned_trimmed.fasta"
            run(["trimal", "-in", str(aln_path), "-out", str(trimmed), "-automated1"])
            aln_path = trimmed

    logging.info("")
    logging.info(">>> STOP AND LOOK AT THE ALIGNMENT before building a tree.")
    logging.info(">>> Open %s in AliView or Jalview. Check the ends, check for", aln_path)
    logging.info(">>> frame shifts, check nothing is obviously misaligned.")
    logging.info("")

    if args.stop_after == "align":
        return 0

    # --- 3. ML tree -------------------------------------------------------
    if not have("iqtree") and not have("iqtree2"):
        logging.error("iqtree not found on PATH. conda install -c bioconda iqtree")
        return 1
    iqtree = "iqtree2" if have("iqtree2") else "iqtree"

    tree_prefix = args.outdir / f"{args.prefix}_ml"
    run([iqtree, "-s", str(aln_path), "-m", "MFP",
         "-B", str(args.bootstrap), "--alrt", "1000",
         "-T", str(args.threads), "--prefix", str(tree_prefix), "-redo"])

    logging.info("")
    logging.info("Tree: %s.treefile", tree_prefix)
    logging.info("Model selected: see %s.iqtree", tree_prefix)
    logging.info("")
    logging.info("NEXT:")
    logging.info("  1. Open the treefile in FigTree. Does the lineage structure look sane?")
    logging.info("  2. Check for a vaccine-strain clade — any field isolate sitting inside")
    logging.info("     it is a mislabelled vaccine sequence. Add it to config/vaccine_strains.txt")
    logging.info("     and re-run from step 02.")
    logging.info("  3. Load %s.treefile into TempEst with %s", tree_prefix, dates_path)
    logging.info("     and check the root-to-tip regression. R-squared and a positive slope")
    logging.info("     tell you whether tip-dated BEAST analysis is viable.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

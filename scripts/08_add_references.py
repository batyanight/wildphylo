#!/usr/bin/env python3
"""
08_add_references.py — fetch lineage reference sequences and build a classification set.

Downloads the reference accessions in config/lineage_references.tsv, extracts H gene
sequences, and writes a combined FASTA of your sequences plus references, ready for
re-alignment and tree building. The resulting tree tells you which named lineage your
clade falls in.

Protein accessions are resolved to their coding nucleotide record automatically via the
GenBank DBSOURCE / coded_by annotation; any that cannot be resolved are reported so you
can find the nucleotide equivalent by hand.

Reference tips are labelled:
    REF_<lineage>_<strain>|reference|<year or 0000>

The `reference` host group keeps them out of your host-state analysis — this set is for
lineage identification only, not for the phylodynamic run.

Usage
-----
    python scripts/08_add_references.py --email you@example.com \\
        --query-fasta data/processed/H_clade_3.fasta

    # include the global subsample too, for a fuller picture
    python scripts/08_add_references.py --email you@example.com \\
        --query-fasta data/processed/H_global400.fasta --prefix global400_with_refs

Then realign and rebuild:
    mafft --auto data/processed/clade3_with_refs.fasta > data/processed/clade3_with_refs_aln.fasta
    iqtree2 -s data/processed/clade3_with_refs_aln.fasta -m MFP -B 1000 --alrt 1000 \\
        --prefix data/processed/clade3_with_refs -redo

Open the tree and see which REF_ tips fall inside your clade.
"""

import argparse
import re
import sys
import time
from pathlib import Path

try:
    from Bio import Entrez, SeqIO
except ImportError as exc:
    raise ImportError(
        "Biopython is required.\n"
        "  terminal: pip install biopython\n"
        "  notebook: %pip install biopython   (then restart the kernel)"
    ) from exc

H_NAMES = {"h", "ha", "hemagglutinin", "haemagglutinin", "attachment protein"}
H_MIN, H_MAX = 1500, 1900


def load_references(path: Path):
    refs = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            print(f"  skipping malformed line: {line[:60]}", file=sys.stderr)
            continue
        refs.append({"accession": parts[0].strip(), "lineage": parts[1].strip(),
                     "strain": parts[2].strip(), "type": parts[3].strip()})
    return refs


def resolve_protein(acc: str) -> str | None:
    """Find the nucleotide accession coding for a protein accession."""
    try:
        handle = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode="text")
        text = handle.read()
        handle.close()
    except Exception as exc:                          # noqa: BLE001
        print(f"    could not fetch protein {acc}: {exc}")
        return None
    m = re.search(r'coded_by="(?:complement\()?([A-Z]+\d+(?:\.\d+)?)', text)
    if m:
        return m.group(1)
    m = re.search(r"DBSOURCE\s+.*?accession ([A-Z]+\d+(?:\.\d+)?)", text)
    return m.group(1) if m else None


def extract_H(record):
    """H gene from annotation, falling back to whole record if length is right."""
    for feat in record.features:
        if feat.type not in ("CDS", "gene", "mat_peptide"):
            continue
        label = ""
        for key in ("gene", "product", "note"):
            if key in feat.qualifiers:
                label = feat.qualifiers[key][0].strip().lower()
                break
        if label in H_NAMES or any(label.startswith(n + " ") for n in H_NAMES):
            try:
                seq = str(feat.extract(record.seq))
            except Exception:                          # noqa: BLE001
                continue
            if H_MIN <= len(seq) <= H_MAX:
                return seq, "annotation"
    if H_MIN <= len(record.seq) <= H_MAX:
        return str(record.seq), "whole_record"
    return "", ""


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, chunks = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def write_fasta(path: Path, seqs: dict[str, str], width: int = 60) -> None:
    with path.open("w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "-", text)[:40]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", required=True)
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--refs", default=Path("config/lineage_references.tsv"), type=Path)
    ap.add_argument("--query-fasta", required=True, type=Path,
                    help="Your sequences — UNALIGNED or aligned; gaps are stripped")
    ap.add_argument("--outdir", default=Path("data/processed"), type=Path)
    ap.add_argument("--prefix", default="clade3_with_refs")
    ap.add_argument("--cache", default=Path("data/raw/references.gb"), type=Path)
    args = ap.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    if not args.refs.is_file():
        print(f"Reference table not found: {args.refs}", file=sys.stderr)
        return 1
    if not args.query_fasta.is_file():
        print(f"Query FASTA not found: {args.query_fasta}", file=sys.stderr)
        return 1

    refs = load_references(args.refs)
    print(f"{len(refs)} reference accessions listed")
    by_lineage = {}
    for r in refs:
        by_lineage.setdefault(r["lineage"], 0)
        by_lineage[r["lineage"]] += 1
    for lin, n in sorted(by_lineage.items()):
        print(f"  {lin:<22} {n}")
    print()

    # --- resolve protein accessions to nucleotide ---
    nuc_ids, meta, unresolved = [], {}, []
    for r in refs:
        acc = r["accession"]
        if r["type"] == "prot":
            print(f"  resolving protein {acc} ...", end=" ")
            nuc = resolve_protein(acc)
            time.sleep(0.4)
            if not nuc:
                print("FAILED")
                unresolved.append(r); continue
            print(f"-> {nuc}")
            acc = nuc
        nuc_ids.append(acc)
        meta[acc.split(".")[0]] = r

    if unresolved:
        print(f"\n{len(unresolved)} protein accessions could not be resolved:")
        for r in unresolved:
            print(f"    {r['accession']:<12} {r['lineage']:<20} {r['strain']}")
        print("  Find the nucleotide equivalent on NCBI and add it to the table as type=nuc.\n")

    # --- fetch nucleotide records ---
    print(f"fetching {len(nuc_ids)} nucleotide records ...")
    args.cache.parent.mkdir(parents=True, exist_ok=True)
    try:
        handle = Entrez.efetch(db="nuccore", id=",".join(nuc_ids),
                               rettype="gb", retmode="text")
        args.cache.write_text(handle.read())
        handle.close()
    except Exception as exc:                           # noqa: BLE001
        print(f"fetch failed: {exc}", file=sys.stderr)
        return 1

    ref_seqs, failed = {}, []
    for record in SeqIO.parse(str(args.cache), "genbank"):
        base = record.id.split(".")[0]
        info = meta.get(base)
        if not info:
            continue
        seq, how = extract_H(record)
        if not seq:
            failed.append((record.id, info["lineage"], info["strain"], len(record.seq)))
            continue
        year = "0000"
        for feat in record.features:
            if feat.type == "source" and "collection_date" in feat.qualifiers:
                m = re.search(r"(\d{4})", feat.qualifiers["collection_date"][0])
                if m:
                    year = m.group(1)
                break
        label = (f"REF_{safe(info['lineage'])}_{safe(info['strain'])}_{base}"
                 f"|reference|{year}")
        ref_seqs[label] = seq
        print(f"  {record.id:<12} {info['lineage']:<22} {len(seq):>5} nt  ({how})")

    if failed:
        print(f"\n{len(failed)} records had no extractable H sequence:")
        for acc, lin, strain, ln in failed:
            print(f"    {acc:<12} {lin:<20} {strain:<20} record length {ln}")
        print("  Complete genomes without an H annotation need manual extraction.\n")

    if not ref_seqs:
        print("No reference sequences recovered — nothing to merge.", file=sys.stderr)
        return 1

    # --- merge with query sequences ---
    query = read_fasta(args.query_fasta)
    query = {k: v.replace("-", "").replace(".", "") for k, v in query.items()}
    combined = {**query, **ref_seqs}

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"{args.prefix}.fasta"
    write_fasta(out, combined)

    print(f"\n{len(query)} query + {len(ref_seqs)} reference = {len(combined)} sequences")
    print(f"wrote {out}")
    print(f"""
NEXT

  mafft --auto {out} > {args.outdir / (args.prefix + '_aln.fasta')}
  iqtree2 -s {args.outdir / (args.prefix + '_aln.fasta')} -m MFP -B 1000 --alrt 1000 \\
      -T AUTO --prefix {args.outdir / args.prefix} -redo

Then open {args.prefix}.treefile and answer two questions:

  1. Which REF_ tips fall INSIDE your clade? That is your lineage assignment.
     Expected for a North American raccoon-associated clade: America-2, possibly America-3.

  2. Do any of YOUR sequences fall inside the America-1/vaccine clade
     (REF_America-1-vaccine_Onderstepoort, REF_America-1-vaccine_Snyder-Hill)?
     Any that do are vaccine-derived sequences the name filter missed. Add them to
     config/vaccine_strains.txt and re-run the pipeline from step 02.

  3. Do the 4 Denmark sequences and the Meles meles record group with
     REF_Europe-wildlife_Danish-mink? That would explain them as genuine European
     wildlife lineage sequences rather than a clade boundary artefact.
""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

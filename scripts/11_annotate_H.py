#!/usr/bin/env python3
"""
11_annotate_H.py — add CDV H gene and domain annotations to an Auspice JSON.

The entropy panel is far more informative when the x-axis carries the protein's
domain structure: you can see immediately whether variable sites fall in the
receptor-binding head or the conserved transmembrane region.

RUN THIS BEFORE 10_add_entropy.py. This script establishes the coordinate system;
10 then reconstructs mutations into it. Running them the other way round leaves
the mutations in alignment coordinates and the panel will be wrong.

WHAT CHANGED, AND WHY
---------------------
Earlier versions wrote CDS features as *alignment column* spans. That is only
valid if the reference has no gaps inside the feature. It does: this alignment
carries insertions relative to the reference inside the head domain, so head came
out four columns too long, stopped being a multiple of 3, and Auspice discarded
it. The same thing happened to any feature spanning those columns, which is what
the old --slam and --max-feature-aa workarounds were really fighting. Splitting a
feature cannot fix a frame problem, and the amino acids downstream of an
insertion were being read out of frame.

Features are now written in reference nucleotide coordinates, where every CDS is
an exact multiple of 3 by construction. The old --max-feature-aa subdivision and
the SLAM-splitting logic are gone; Auspice 2.46+ stacks overlapping CDSs on
separate rows, so the SLAM sites can simply sit on top of the head domain.

Two things still have to be right:

  1. The alignment includes gap padding, so alignment column != H nucleotide
     position. This script maps between them using a reference sequence.
  2. The reading frame must be known. The script finds the longest ATG-initiated
     ORF and refuses to write an annotation if it does not look like H, rather
     than producing a plausible wrong one.

Domain boundaries (Beineke et al. / Sattler et al., aa positions in H protein)
live in cdv_h_coords.py.

Usage
-----
    python scripts/11_annotate_H.py \\
        --json auspice/cdv-phylodynamics.json \\
        --alignment data/processed/H_clade_3.fasta \\
        --reference "MK423844|domestic_dog|2018.500" \\
        --domains --slam \\
        --output auspice/cdv-phylodynamics.json

    # let the script pick the least-gapped sequence as reference
    python scripts/11_annotate_H.py --json ... --alignment ... --auto-reference \\
        --domains --output ...

Requires Auspice 2.46 or newer for overlapping features to render.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cdv_h_coords import (  # noqa: E402
    HFrame, read_fasta, pick_reference, resolve_features, build_annotations,
    backup_if_overwriting,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", required=True, type=Path)
    ap.add_argument("--alignment", required=True, type=Path)
    ap.add_argument("--reference", help="Taxon label to use as the coordinate reference")
    ap.add_argument("--auto-reference", action="store_true",
                    help="Use the least-gapped sequence")
    ap.add_argument("--domains", action="store_true",
                    help="Emit the domain breakdown as well as the full-length gene. "
                         "The gene feature is always written, so the amino-acid entropy "
                         "view works either way once 10_add_entropy.py has run.")
    ap.add_argument("--slam", action="store_true",
                    help="Include the two SLAM-binding sites as separate features. "
                         "They overlap the head domain, which Auspice 2.46+ renders as "
                         "a stacked row. Implies --domains.")
    ap.add_argument("--gene-name", default="H",
                    help="Name of the full-length CDS feature.")
    ap.add_argument("--offset", type=int, default=0,
                    help=argparse.SUPPRESS)
    ap.add_argument("--max-feature-aa", type=int, default=None,
                    help=argparse.SUPPRESS)
    ap.add_argument("--no-backup", action="store_true",
                    help="Do not keep a timestamped copy of the output file before "
                         "overwriting it.")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    for dead, msg in (("offset", "coordinates are now anchored to the reference ORF"),
                      ("max_feature_aa", "feature size was never the problem; see the "
                                         "module docstring")):
        if getattr(args, dead) not in (0, None):
            print(f"note: --{dead.replace('_', '-')} is ignored — {msg}", file=sys.stderr)

    if args.slam:
        args.domains = True

    if not args.json.is_file():
        print(f"not found: {args.json}", file=sys.stderr)
        return 1
    seqs = read_fasta(args.alignment)
    if not seqs:
        print("no sequences read", file=sys.stderr)
        return 1

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        print(f"alignment sequences differ in length {sorted(lengths)[:4]} — "
              "this file is not aligned", file=sys.stderr)
        return 1

    try:
        ref_name, ref_seq = pick_reference(seqs, args.reference, args.auto_reference)
    except KeyError as e:
        print(e, file=sys.stderr)
        return 1
    if args.auto_reference or not args.reference:
        print(f"reference (least gapped): {ref_name}")
    else:
        print(f"reference: {ref_name}")

    frame = HFrame(ref_name, ref_seq)
    print(f"reference length: {frame.ref_length} nt ungapped, "
          f"{frame.n_columns} alignment columns")

    # --- locate the coding sequence ---
    orf_start, orf_nt = frame.find_orf()
    stops_naive = {f: frame.internal_stops(f) for f in (0, 1, 2)}
    print(f"internal stops if the sequence began at base 1: {stops_naive}")

    if orf_start is None or orf_nt < 1500:
        print(f"\nERROR: the longest ATG-initiated ORF is only {orf_nt} nt "
              f"({orf_nt // 3} aa). CDV H is ~1824 nt (604-607 aa). This reference is "
              "either not a clean H coding sequence or is too fragmentary. Choose a "
              "different reference, e.g. add A75/17 (AF164967) to the alignment and "
              "pass --reference for it. No annotation written.", file=sys.stderr)
        return 1

    n_aa = frame.n_codons
    print(f"longest ORF: starts at reference base {orf_start + 1}, {orf_nt} nt "
          f"({n_aa} aa), frame {orf_start % 3}")
    if abs(n_aa - 605) > 25:
        print(f"\nWARNING: the ORF encodes {n_aa} aa, which is well away from the "
              "expected 604-607. Domain positions may be shifted. Verify before using "
              "this figure in a manuscript.", file=sys.stderr)
    if n_aa < 550:
        print(f"\nWARNING: the reference covers only ~{n_aa} aa of H. Domain "
              "annotations beyond that will be truncated or wrong. Consider adding a "
              "complete H reference sequence to the alignment.", file=sys.stderr)

    # How much of the CDS is affected by insertions in other sequences. This is
    # the quantity that used to silently corrupt the annotation.
    cds_cols = frame.col_for_refpos[orf_start + 1], frame.col_for_refpos[orf_start + orf_nt]
    insertions = (cds_cols[1] - cds_cols[0] + 1) - orf_nt
    if insertions:
        print(f"note: {insertions} alignment column(s) inside the CDS are insertions "
              f"relative to the reference; codons are assembled around them")

    # --- build features ---
    feats = resolve_features(frame, gene_name=args.gene_name,
                             domains=args.domains, slam=args.slam)
    ann = build_annotations(frame, feats, gene_name=args.gene_name)

    print("\nFEATURE -> REFERENCE NUCLEOTIDE COORDINATES")
    for name, a1, a2 in feats:
        f = ann[name]
        span = f["end"] - f["start"] + 1
        print(f"  {name:<16} aa {a1:>3}-{a2:<3}  ->  nt {f['start']}-{f['end']} "
              f"({span} nt, {span // 3} codons)")

    doc = json.loads(args.json.read_text())
    doc["meta"]["genome_annotations"] = ann
    if "entropy" not in doc["meta"].get("panels", []):
        doc["meta"].setdefault("panels", []).append("entropy")

    # Any mutation keys left over from an earlier annotation scheme refer to
    # features that no longer exist. Auspice ignores them, but they inflate the
    # file and make the JSON confusing to debug, so report them here.
    known = set(ann) | {"nuc"}
    stale = set()

    def scan(node):
        for k in node.get("branch_attrs", {}).get("mutations", {}):
            if k not in known:
                stale.add(k)
        for c in node.get("children", []):
            scan(c)

    scan(doc["tree"])
    if stale:
        print(f"\nnote: {len(stale)} stale mutation key(s) in the tree refer to features "
              f"that no longer exist: {', '.join(sorted(stale)[:6])}"
              f"{' ...' if len(stale) > 6 else ''}")
        print("      10_add_entropy.py will clear them when it rebuilds mutations.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    backup_if_overwriting(args.json, args.output, not args.no_backup)
    args.output.write_text(json.dumps(doc, indent=1))
    print(f"\nwrote {args.output} ({args.output.stat().st_size / 1024:.0f} KB)")
    print(f"{len(feats)} feature(s) added, all in reference nucleotide coordinates.")
    print("\nNow run 10_add_entropy.py with the same --alignment and --reference to "
          "rebuild mutations in these coordinates.")
    print("\nCaveat for the methods: domain boundaries are aa positions from the "
          "literature, mapped onto the reference ORF. State which reference you used.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

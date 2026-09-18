#!/usr/bin/env python3
"""
10_add_entropy.py — add the entropy panel to an Auspice JSON.

The entropy panel shows per-site nucleotide (or amino acid) diversity across the
gene, and lets you colour the tree by mutations at any site. It needs two things
the basic export does not have:

  1. a genome annotation, so Auspice knows the coordinate system
  2. per-node mutations, from ancestral sequence reconstruction

Ancestral states are reconstructed by Fitch parsimony. This is not maximum
likelihood (augur uses TreeTime) but for a single gene alignment with modest
divergence it gives essentially the same mutation assignments, and it has no
dependencies.

RUN 11_annotate_H.py FIRST. It establishes the coordinate system; this script
reconstructs mutations into it and will refuse to run without it.

WHAT CHANGED, AND WHY
---------------------
Mutations used to be emitted in *alignment column* coordinates while amino acids
were translated from raw column slices. Where the alignment carries an insertion
relative to the reference — and this one does, inside the head domain — the slice
picks up extra columns, the reading frame shifts, and every amino acid downstream
is wrong. Positions are now reference nucleotide coordinates and codons are
assembled from the columns holding the reference's own bases, so insertions are
stepped over rather than read through. Mutations falling in insertion columns
have no reference coordinate and are dropped, with a count reported; this is what
augur does too.

Amino acid mutations are emitted for *every* CDS in the annotation, keyed by
feature name and numbered from 1 within each feature, which is how Auspice looks
them up. With the full-length H feature present you get a whole-gene amino acid
view; with --domains on script 11 you also get one view per domain.

Mutation keys are now rebuilt from scratch on each run. Previously the script
wrote into whatever was already in the JSON, so keys from earlier annotation
schemes accumulated and never went away. Empty mutation lists are also omitted
rather than written for every feature on every node.

Usage
-----
    python scripts/11_annotate_H.py --json ... --alignment ... --domains --output ...
    python scripts/10_add_entropy.py \\
        --json auspice/cdv-phylodynamics.json \\
        --alignment data/processed/H_clade_3.fasta \\
        --reference "MK423844|domestic_dog|2018.500" \\
        --output auspice/cdv-phylodynamics.json

Use the same --alignment and --reference for both scripts. Tip names in the
alignment must match the tree exactly.
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cdv_h_coords import (  # noqa: E402
    AMBIG, BASES, CODON, DOMAINS, HFrame, read_fasta, pick_reference,
    annotation_feature_spans, backup_if_overwriting,
)


def collect(node, out):
    out.append(node)
    for c in node.get("children", []):
        collect(c, out)
    return out


def fitch_up(node, seqs, L):
    """Upward pass: build the candidate state set at every site."""
    if "children" not in node:
        s = seqs.get(node["name"])
        if s is None:
            return [set()] * L
        return [({c} if c not in AMBIG else set()) for c in s[:L]]

    child_sets = [fitch_up(c, seqs, L) for c in node["children"]]
    mine = []
    for i in range(L):
        sets = [cs[i] for cs in child_sets if cs[i]]
        if not sets:
            mine.append(set())
            continue
        inter = set.intersection(*sets)
        mine.append(inter if inter else set.union(*sets))
    node["_fitch"] = mine
    return mine


def assign_down(node, parent_state, L, states, nuc_muts, frame, counters):
    """Downward pass: resolve states, record nucleotide mutations in reference
    coordinates, and cache the resolved state for later translation."""
    if "children" not in node:
        s = node.get("_seq")
        # Ambiguities and gaps become None. Without this, a gap in a partial
        # sequence compares unequal to the parent base and is recorded as a
        # mutation, which manufactures huge numbers of spurious changes.
        state = [(c if c not in AMBIG else None) for c in s[:L]] if s else [None] * L
    else:
        f = node["_fitch"]
        state = []
        for i in range(L):
            if parent_state and parent_state[i] and parent_state[i] in f[i]:
                state.append(parent_state[i])
            elif f[i]:
                state.append(sorted(f[i])[0])
            else:
                state.append(parent_state[i] if parent_state else None)

    muts = []
    if parent_state:
        for i in range(L):
            a, b = parent_state[i], state[i]
            if a in BASES and b in BASES and a != b:
                refpos = frame.refpos_for_col.get(i + 1)
                if refpos is None:
                    # Column is an insertion relative to the reference, so the
                    # site has no coordinate in this dataset.
                    counters["dropped"] += 1
                    continue
                muts.append(f"{a}{refpos}{b}")

    states[node["name"]] = state
    nuc_muts[node["name"]] = muts

    for c in node.get("children", []):
        assign_down(c, state, L, states, nuc_muts, frame, counters)
    return state


def parent_map(node, parent=None, out=None):
    out = {} if out is None else out
    out[node["name"]] = parent
    for c in node.get("children", []):
        parent_map(c, node["name"], out)
    return out


def translate(state, triples):
    """Resolved state -> amino acids, one per codon-column triple."""
    aas = []
    for c1, c2, c3 in triples:
        codon = "".join(
            state[c - 1] if state[c - 1] in BASES else "N" for c in (c1, c2, c3)
        )
        aas.append(CODON.get(codon, "X"))
    return aas


def final_checks(doc, gene_name: str) -> bool:
    """Verify the two invariants whose violation caused the original bug.

    A CDS that is not a multiple of 3 gets silently discarded by Auspice. And if
    the per-domain amino acid counts do not sum to the whole-gene count, some
    feature is reading a different frame than the gene — which is exactly the
    failure that alignment-column coordinates produced, and it is invisible in
    the rendered panel.
    """
    ok = True
    ann = doc["meta"]["genome_annotations"]

    bad = [k for k, v in ann.items()
           if k != "nuc" and (v["end"] - v["start"] + 1) % 3]
    if bad:
        print(f"CHECK FAILED: CDS not a multiple of 3: {bad}", file=sys.stderr)
        ok = False
    else:
        print(f"  all {len(ann) - 1} CDS lengths divisible by 3")

    counts = Counter()

    def walk(n):
        for k, v in n.get("branch_attrs", {}).get("mutations", {}).items():
            counts[k] += len(v)
        for c in n.get("children", []):
            walk(c)

    walk(doc["tree"])

    present = [d for d, _, _ in DOMAINS if d in ann]
    if gene_name in ann and present:
        whole, parts = counts.get(gene_name, 0), sum(counts[d] for d in present)
        if whole == parts:
            print(f"  {gene_name} amino acid mutations ({whole}) == sum over "
                  f"{len(present)} domains")
        else:
            print(f"CHECK FAILED: {gene_name} has {whole} amino acid mutations but the "
                  f"domains sum to {parts}. A feature is out of frame with the gene.",
                  file=sys.stderr)
            ok = False

    stale = [k for k in counts if k != "nuc" and k not in ann]
    if stale:
        print(f"CHECK FAILED: mutation keys with no matching feature: {stale}",
              file=sys.stderr)
        ok = False
    else:
        print("  no orphaned mutation keys")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", required=True, type=Path)
    ap.add_argument("--alignment", required=True, type=Path)
    ap.add_argument("--reference", help="Taxon label defining the coordinate system. "
                                        "Must match the one given to 11_annotate_H.py.")
    ap.add_argument("--auto-reference", action="store_true",
                    help="Use the least-gapped sequence")
    ap.add_argument("--gene", default="H", help=argparse.SUPPRESS)
    ap.add_argument("--gene-start", type=int, default=None, help=argparse.SUPPRESS)
    ap.add_argument("--cds", nargs=2, type=int, metavar=("START", "END"),
                    help=argparse.SUPPRESS)
    ap.add_argument("--cds-from-json", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--no-backup", action="store_true",
                    help="Do not keep a timestamped copy of the output file before "
                         "overwriting it.")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    if args.gene_start is not None or args.cds:
        print("note: --gene-start and --cds are ignored — CDS extents now come from "
              "the annotation written by 11_annotate_H.py", file=sys.stderr)
    if args.cds_from_json:
        print("note: --cds-from-json is now the only behaviour and can be dropped",
              file=sys.stderr)

    for f in (args.json, args.alignment):
        if not f.is_file():
            print(f"not found: {f}", file=sys.stderr)
            return 1

    doc = json.loads(args.json.read_text())
    seqs = read_fasta(args.alignment)
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        print(f"alignment sequences differ in length {sorted(lengths)[:4]} — "
              "this file is not aligned", file=sys.stderr)
        return 1
    L = lengths.pop()

    ann = doc["meta"].get("genome_annotations")
    if not ann or set(ann) <= {"nuc"}:
        print("ERROR: the JSON has no CDS features. Run 11_annotate_H.py first — it "
              "establishes the coordinate system this script reconstructs into.",
              file=sys.stderr)
        return 1

    # --- rebuild the same coordinate frame script 11 used ---
    try:
        ref_name, ref_seq = pick_reference(seqs, args.reference, args.auto_reference)
    except KeyError as e:
        print(e, file=sys.stderr)
        return 1
    frame = HFrame(ref_name, ref_seq)
    frame.find_orf()
    print(f"reference: {ref_name} ({frame.ref_length} nt ungapped, {L} columns)")

    if frame.orf_start is None:
        print("ERROR: no ORF found in the reference.", file=sys.stderr)
        return 1
    if ann["nuc"]["end"] != frame.ref_length:
        print(f"ERROR: the annotation says nuc is {ann['nuc']['end']} nt but this "
              f"reference is {frame.ref_length} nt. The two scripts were given "
              "different --reference or --alignment values.", file=sys.stderr)
        return 1

    # --- map each annotated CDS back to codon-column triples ---
    regions = []
    for name, start, end in annotation_feature_spans(ann):
        span = end - start + 1
        offset = start - frame.orf_start - 1
        if offset < 0 or offset % 3 or span % 3:
            print(f"ERROR: feature {name} (nt {start}-{end}) is not in frame with the "
                  f"reference ORF starting at {frame.orf_start + 1}. Re-run "
                  "11_annotate_H.py.", file=sys.stderr)
            return 1
        a1 = offset // 3 + 1
        a2 = a1 + span // 3 - 1
        if a2 > frame.n_codons:
            print(f"ERROR: feature {name} extends to aa {a2}, past the {frame.n_codons} "
                  "aa ORF. Re-run 11_annotate_H.py.", file=sys.stderr)
            return 1
        regions.append((name, a1, a2, frame.feature_codon_columns(a1, a2)))
    print(f"translating {len(regions)} CDS feature(s) from the JSON annotation")

    nodes = collect(doc["tree"], [])
    tips = [n for n in nodes if "children" not in n]
    matched = sum(1 for t in tips if t["name"] in seqs)
    print(f"{len(tips)} tips in tree, {len(seqs)} in alignment, {matched} matched")
    if matched < len(tips) * 0.9:
        print("ERROR: fewer than 90% of tips matched. Tip names must be identical "
              "between the tree and the alignment.", file=sys.stderr)
        missing = [t["name"] for t in tips if t["name"] not in seqs][:3]
        print(f"  unmatched examples: {missing}", file=sys.stderr)
        return 1

    for t in tips:
        t["_seq"] = seqs.get(t["name"])

    # --- ancestral reconstruction, once ---
    print(f"reconstructing ancestral states over {L} sites...")
    sys.setrecursionlimit(10000)
    fitch_up(doc["tree"], seqs, L)
    states, nuc_muts, counters = {}, {}, {"dropped": 0}
    assign_down(doc["tree"], None, L, states, nuc_muts, frame, counters)

    total = sum(len(v) for v in nuc_muts.values())
    on_branches = sum(1 for v in nuc_muts.values() if v)
    per_branch = total / on_branches if on_branches else 0
    print(f"{total} nucleotide mutations assigned across {on_branches} branches "
          f"({per_branch:.1f} per branch)")
    if counters["dropped"]:
        print(f"{counters['dropped']} mutation(s) dropped: they fall in alignment "
              "columns that are insertions relative to the reference and so have no "
              "coordinate in this dataset")

    # A branch cannot plausibly change more sites than exist, and in practice
    # should change far fewer. This catches gap/ambiguity handling errors.
    n_var = sum(1 for i in range(L)
                if len({s[i] for s in seqs.values() if s[i] in BASES}) > 1)
    if per_branch > n_var * 0.25:
        print(f"\nWARNING: {per_branch:.0f} mutations per branch against only "
              f"{n_var} variable sites. That is implausible and usually means "
              f"gaps or ambiguity codes are being treated as real bases.",
              file=sys.stderr)

    # --- amino acid mutations, per feature ---
    pmap = parent_map(doc["tree"])
    translations = {}
    aa_muts = {}
    for name, a1, a2, triples in regions:
        cache = {}
        n_aa = 0
        for node_name, st in states.items():
            par = pmap.get(node_name)
            if par is None:
                continue
            if par not in cache:
                cache[par] = translate(states[par], triples)
            if node_name not in cache:
                cache[node_name] = translate(st, triples)
            a, b = cache[par], cache[node_name]
            ms = [f"{x}{i + 1}{y}" for i, (x, y) in enumerate(zip(a, b))
                  if x != y and x not in "X*" and y not in "X*"]
            if ms:
                aa_muts.setdefault(node_name, {})[name] = ms
                n_aa += len(ms)
        translations[name] = n_aa
        print(f"  {name:<16} aa {a1:>3}-{a2:<3} ({a2 - a1 + 1} codons): "
              f"{n_aa} amino-acid mutations")

    # --- write mutations, rebuilding the block from scratch ---
    known = {name for name, _, _, _ in regions}
    stale = set()
    for n in nodes:
        old = n.get("branch_attrs", {}).get("mutations", {})
        stale |= {k for k in old if k != "nuc" and k not in known}
        fresh = {}
        m = nuc_muts.get(n["name"], [])
        if m:
            fresh["nuc"] = m
        for gene, ms in aa_muts.get(n["name"], {}).items():
            fresh[gene] = ms
        n.setdefault("branch_attrs", {})["mutations"] = fresh
        for k in ("_fitch", "_seq"):
            n.pop(k, None)
    if stale:
        print(f"\ncleared {len(stale)} stale mutation key(s) left by earlier runs: "
              f"{', '.join(sorted(stale)[:6])}{' ...' if len(stale) > 6 else ''}")

    if "entropy" not in doc["meta"].get("panels", []):
        doc["meta"].setdefault("panels", []).append("entropy")

    # variable-site summary, for your own information
    var = 0
    for i in range(L):
        col = Counter(s[i] for s in seqs.values() if s[i] not in AMBIG)
        if len(col) > 1:
            var += 1
    print(f"{var} of {L} columns are variable ({var / L:.1%})")

    print("\nchecks:")
    if not final_checks(doc, args.gene):
        print("\nNot writing — see the failed check above.", file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    backup_if_overwriting(args.json, args.output, not args.no_backup)
    args.output.write_text(json.dumps(doc, indent=1))
    print(f"\nwrote {args.output} ({args.output.stat().st_size / 1024:.0f} KB)")
    print("The entropy panel will now appear below the tree in Auspice, with an "
          "amino-acid view for every annotated feature.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

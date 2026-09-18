#!/usr/bin/env python3
"""
09_make_auspice.py — convert a BEAST MCC tree directly to an Auspice v2 JSON.

Does what `augur import beast` + `augur export v2` would do, without the augur
dependency chain (which pulls in cvxopt and a compiler toolchain for features this
project does not use).

Reads the annotated MCC tree from TreeAnnotator, preserves node dates and the
reconstructed discrete trait, joins tip metadata, and writes a JSON that renders in
Auspice — locally or by dragging onto https://auspice.us

Usage
-----
    python scripts/09_make_auspice.py \\
        --mcc beast/clade3/run100M_fixed/clade3_mcc_fixed.tree \\
        --metadata data/processed/H_clade_3_metadata.tsv \\
        --full-metadata data/processed/metadata_clean.tsv \\
        --most-recent 2023.622 \\
        --output auspice/cdv_clade3.json

Then either drag the JSON onto auspice.us, or run `auspice view --datasetDir auspice`.
"""

import argparse
import datetime
import json
import re
import sys
from pathlib import Path


def parse_mcc(path: Path):
    raw = path.read_text()
    tr = {}
    m = re.search(r"Translate(.*?);", raw, re.S)
    if m:
        for line in m.group(1).replace("\n", " ").split(","):
            p = line.split()
            if len(p) >= 2:
                tr[p[0].strip()] = p[1].strip().strip("'\"")
    tl = [l for l in raw.splitlines() if l.strip().startswith("tree ")]
    if not tl:
        raise ValueError("no tree line found in the NEXUS file")
    nwk = tl[0][tl[0].index("("):].strip().rstrip(";")
    return tr, nwk


def parse_annotation(s):
    d = {}
    for a in re.finditer(r'([\w.%]+)=(?:"([^"]*)"|\{([^}]*)\}|([^,\]]+))', s):
        d[a.group(1)] = a.group(2) or a.group(3) or a.group(4)
    return d


class Node:
    __slots__ = ("children", "ann", "label", "length")
    def __init__(self):
        self.children = []; self.ann = {}; self.label = None; self.length = 0.0


def build(tr, nwk):
    i = [0]
    def parse():
        n = Node()
        if nwk[i[0]] == "(":
            i[0] += 1
            while True:
                n.children.append(parse())
                if nwk[i[0]] == ",": i[0] += 1; continue
                if nwk[i[0]] == ")": i[0] += 1; break
        else:
            j = i[0]
            while nwk[i[0]] not in "[:,)": i[0] += 1
            key = nwk[j:i[0]].strip()
            n.label = tr.get(key, key)
        if i[0] < len(nwk) and nwk[i[0]] == "[":
            j = nwk.index("]", i[0]); n.ann = parse_annotation(nwk[i[0]+1:j]); i[0] = j + 1
        if i[0] < len(nwk) and nwk[i[0]] == ":":
            i[0] += 1; j = i[0]
            while i[0] < len(nwk) and nwk[i[0]] not in ",)[": i[0] += 1
            try: n.length = float(nwk[j:i[0]])
            except ValueError: n.length = 0.0
            if i[0] < len(nwk) and nwk[i[0]] == "[":
                j = nwk.index("]", i[0]); i[0] = j + 1
        return n
    return parse()


def load_metadata(sub_path, full_path):
    """taxon label -> dict of attributes."""
    meta = {}
    if not sub_path or not Path(sub_path).is_file():
        return meta
    import csv
    with open(sub_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    extra = {}
    if full_path and Path(full_path).is_file():
        with open(full_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                extra[str(r.get("accession", "")).split(".")[0]] = r
    for r in rows:
        lab = r.get("label") or ""
        if not lab:
            continue
        acc = r.get("accession", "")
        e = extra.get(str(acc).split(".")[0], {})
        meta[lab] = {
            "accession": acc,
            "host": r.get("host_group", ""),
            "species": e.get("host_canonical", "") or r.get("host_canonical", ""),
            "country": e.get("country", ""),
            "date_precision": r.get("date_precision", ""),
        }
    return meta


def convert(node, most_recent, meta, trait_key, trait_name, div=0.0):
    """Recursively build the Auspice node dict."""
    out = {}
    height = float(node.ann.get("height", 0.0) or 0.0)
    num_date = most_recent - height
    div = div + (node.length or 0.0)

    attrs = {"div": round(div, 6), "num_date": {"value": round(num_date, 4)}}

    lo = node.ann.get("height_95%_HPD")
    if lo:
        try:
            a, b = [float(x) for x in re.split(r"[,\s]+", lo.strip()) if x]
            attrs["num_date"]["confidence"] = [round(most_recent - max(a, b), 4),
                                               round(most_recent - min(a, b), 4)]
        except (ValueError, TypeError):
            pass

    if trait_key in node.ann:
        entry = {"value": node.ann[trait_key]}
        # BEAST_CLASSIC writes the full posterior over states as
        # location.set={a,b,c} and location.set.prob={0.5,0.3,0.2}.
        # Use it: it shows real uncertainty rather than only the modal state.
        states = node.ann.get(f"{trait_key}.set")
        probs = node.ann.get(f"{trait_key}.set.prob")
        if states and probs:
            try:
                sl = [x.strip().strip('"') for x in states.split(",") if x.strip()]
                pl = [float(x) for x in probs.split(",") if x.strip()]
                if len(sl) == len(pl):
                    entry["confidence"] = {s: round(v, 4) for s, v in zip(sl, pl)}
            except ValueError:
                pass
        if "confidence" not in entry:
            pr = node.ann.get(f"{trait_key}.prob")
            if pr:
                try: entry["confidence"] = {node.ann[trait_key]: round(float(pr), 4)}
                except ValueError: pass
        attrs[trait_name] = entry

    post = node.ann.get("posterior")
    if post:
        try: attrs["posterior"] = {"value": round(float(post), 4)}
        except ValueError: pass

    if node.label:                                   # tip
        out["name"] = node.label
        m = meta.get(node.label, {})
        for k in ("country", "species", "accession"):
            if m.get(k):
                attrs[k] = {"value": m[k]}
        if m.get("host") and trait_name not in attrs:
            attrs[trait_name] = {"value": m["host"]}
        elif m.get("host"):
            attrs[trait_name] = {"value": m["host"]}   # tip state is observed, not inferred
    else:
        convert.counter += 1
        out["name"] = f"NODE_{convert.counter:07d}"

    out["node_attrs"] = attrs
    if node.children:
        out["children"] = [convert(c, most_recent, meta, trait_key, trait_name, div)
                           for c in node.children]
    return out


convert.counter = 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mcc", required=True, type=Path)
    ap.add_argument("--metadata", type=Path, help="subset metadata with a 'label' column")
    ap.add_argument("--full-metadata", type=Path, help="metadata_clean.tsv, for country/species")
    ap.add_argument("--most-recent", type=float, required=True,
                    help="Decimal date of the most recent sample, e.g. 2023.622")
    ap.add_argument("--trait-key", default="location",
                    help="Annotation key in the MCC tree (BEAST_CLASSIC writes 'location')")
    ap.add_argument("--trait-name", default="host",
                    help="What to call it in Auspice")
    ap.add_argument("--title", default="America-2 canine distemper virus in North American carnivores")
    ap.add_argument("--maintainer", default="")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    if not args.mcc.is_file():
        print(f"MCC tree not found: {args.mcc}", file=sys.stderr); return 1

    tr, nwk = parse_mcc(args.mcc)
    root = build(tr, nwk)
    meta = load_metadata(args.metadata, args.full_metadata)

    # verify the time axis before writing anything
    def tips(n):
        if not n.children: yield n
        for c in n.children: yield from tips(c)
    bad = 0
    checked = 0
    for t in tips(root):
        try: date = float(t.label.split("|")[-1])
        except (ValueError, AttributeError, IndexError): continue
        h = float(t.ann.get("height", 0) or 0)
        if abs((args.most_recent - date) - h) > 0.05:
            bad += 1
        checked += 1
    print(f"time axis check: {checked - bad}/{checked} tip heights consistent")
    if bad:
        print(f"ERROR: {bad} tips have heights inconsistent with their dates.\n"
              "The time axis is wrong — do not build a figure from this tree.",
              file=sys.stderr)
        return 1

    tree = convert(root, args.most_recent, meta, args.trait_key, args.trait_name)

    colorings = [{"key": args.trait_name, "title": "Host group", "type": "categorical"},
                 {"key": "country", "title": "Country", "type": "categorical"},
                 {"key": "species", "title": "Host species", "type": "categorical"},
                 {"key": "posterior", "title": "Clade posterior", "type": "continuous"},
                 {"key": "num_date", "title": "Date", "type": "continuous"}]

    doc = {
        "version": "v2",
        "meta": {
            "title": args.title,
            "updated": datetime.date.today().isoformat(),
            "build_url": "https://github.com/batyanight/cdv-phylodynamics",
            "panels": ["tree"],
            "colorings": colorings,
            "filters": [args.trait_name, "country", "species"],
            "display_defaults": {"color_by": args.trait_name,
                                 "branch_label": "none",
                                 "distance_measure": "num_date"},
            "description": (
                "Phylodynamics of an America-2 canine distemper virus lineage in North "
                "American carnivores. Branches coloured by reconstructed host state; "
                "deep nodes are poorly resolved and should not be over-interpreted. "
                "Amino acid positions follow the H CDS of A75/17 (AF164967)."),
        },
        "tree": tree,
    }
    if args.maintainer:
        doc["meta"]["maintainers"] = [{"name": args.maintainer}]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=1))

    n_tips = sum(1 for _ in tips(root))
    print(f"{n_tips} tips, {convert.counter} internal nodes")
    print(f"root date {args.most_recent - float(root.ann.get('height', 0)):.1f}")
    print(f"metadata joined for {sum(1 for t in tips(root) if t.label in meta)} tips")
    print(f"\nwrote {args.output}  ({args.output.stat().st_size/1024:.0f} KB)")
    print("\nDrag it onto https://auspice.us to view — no install needed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

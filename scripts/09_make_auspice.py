#!/usr/bin/env python3
"""
09_make_auspice.py — convert a BEAST MCC tree to an Auspice v2 JSON.

Does what `augur import beast` + `augur export v2` would do, without the augur
dependency chain (which pulls cvxopt and a compiler toolchain for features this
project does not use).

    python scripts/09_make_auspice.py --config config/pathogen/cdv-1200.yaml \\
        --mcc builds/.../mcc.tree --metadata builds/.../sub_metadata.tsv \\
        --full-metadata builds/.../metadata_clean.tsv \\
        --gates builds/.../convergence.json builds/.../trait_signal.json \\
        --output auspice/cdv.json

Everything that names the pathogen comes from the config
-------------------------------------------------------
Title, maintainer, build URL, which attributes are colourable and the trait key
are all read from `nextstrain:` and `beast.discrete_trait:`. This used to carry
CDV defaults -- a hardcoded title naming the America-2 lineage, and a
description citing the H CDS of A75/17. In a one-pathogen repository that was
merely untidy. In a pathogen-agnostic one it is a bug that publishes a BTV
build described as canine distemper, and nothing downstream would catch it,
because the JSON is valid and renders.

The time axis is checked before anything is written
---------------------------------------------------
Every tip's reconstructed height must agree with the date in its own label. If
they disagree the time axis is wrong, and a figure built from it is wrong in a
way that looks entirely normal -- dates on the x-axis, nodes in plausible
places, all of it shifted. That check refuses to write rather than warning.

`--most-recent` is derived from the tip labels rather than being required.
Passing it by hand is how the axis gets silently offset: supply last year's
value against this year's tree and every node moves together, so nothing looks
wrong.

Gate verdicts travel with the build
-----------------------------------
`--gates` takes the JSON written by any of the gates (temporal, convergence,
trait signal, date randomisation) and folds their verdicts into the Auspice
description. A published tree is read by people who did not run the pipeline,
and a trait colouring that failed its signal test is indistinguishable from one
that passed unless the build says so.
"""

import argparse
import datetime
import json
import re
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.traittrees import (                                    # noqa: E402
    NewickParser, strip_leading_comment, detect_tag, label_date,
    TREE_LINE_RE, TRANSLATE_ROW_RE,
)


def parse_mcc(path: Path, tag=None):
    """
    Read a TreeAnnotator MCC tree into the flat node list lib/traittrees
    produces, with every annotation kept.

    The parser is shared with 08f and 08g rather than reimplemented here. An
    MCC tree and a posterior tree file are the same grammar, and a fourth
    parser in this repo would be a fourth place for the same quoting and
    nesting edge cases to be got wrong differently.
    """
    raw = path.read_text()
    translate = {}
    in_translate = False
    newick = None
    for line in raw.splitlines():
        s = line.strip()
        low = s.lower()
        if low.startswith("translate"):
            in_translate = True
            continue
        if in_translate:
            if s.startswith(";") or low.startswith("tree "):
                in_translate = False
            else:
                m = TRANSLATE_ROW_RE.match(s)
                if m:
                    translate[m.group(1)] = m.group(2).strip().strip("'\"")
                if s.endswith(";"):
                    in_translate = False
                continue
        m = TREE_LINE_RE.match(s)
        if m and newick is None:
            newick = strip_leading_comment(m.group(1))
    if newick is None:
        raise ValueError(f"no tree line found in {path}")

    resolved = detect_tag(newick, tag)
    nodes = NewickParser(newick, resolved, keep_annotations=True).parse()
    for nd in nodes:
        if nd["is_tip"] and nd["label"] in translate:
            nd["label"] = translate[nd["label"]]
    return nodes, resolved


def load_metadata(sub_path, full_path, dta_states=None):
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
        # The tip's host is collapsed through hosts.dta_states, exactly as 07
        # does when it builds the trait alignment. Without this the tips carry
        # the raw host group (wild_felid) while every internal node carries the
        # model's state (felid), so Auspice renders them as different
        # categories and the colouring implies a transition at every tip.
        raw_host = r.get("host_group", "")
        host = (dta_states or {}).get(raw_host, raw_host) if dta_states else raw_host
        meta[lab] = {
            "accession": acc,
            "host": host,
            "host_group_raw": raw_host,
            "species": e.get("host_canonical", "") or r.get("host_canonical", ""),
            "country": e.get("country", ""),
            "date_precision": r.get("date_precision", ""),
        }
    return meta


def convert(nodes, idx, most_recent, meta, trait_key, trait_name, div=0.0):
    """Recursively build the Auspice node dict from the flat node list."""
    node = nodes[idx]
    out = {}
    height = float(node["ann"].get("height", 0.0) or 0.0)
    num_date = most_recent - height
    div = div + (node["length"] or 0.0)

    attrs = {"div": round(div, 6), "num_date": {"value": round(num_date, 4)}}

    lo = node["ann"].get("height_95%_HPD")
    if lo:
        try:
            a, b = [float(x) for x in re.split(r"[,\s]+", lo.strip()) if x]
            attrs["num_date"]["confidence"] = [round(most_recent - max(a, b), 4),
                                               round(most_recent - min(a, b), 4)]
        except (ValueError, TypeError):
            pass

    if trait_key in node["ann"]:
        entry = {"value": node["ann"][trait_key]}
        # BEAST_CLASSIC writes the full posterior over states as
        # location.set={a,b,c} and location.set.prob={0.5,0.3,0.2}.
        # Use it: it shows real uncertainty rather than only the modal state.
        states = node["ann"].get(f"{trait_key}.set")
        probs = node["ann"].get(f"{trait_key}.set.prob")
        if states and probs:
            try:
                sl = [x.strip().strip('"') for x in states.split(",") if x.strip()]
                pl = [float(x) for x in probs.split(",") if x.strip()]
                if len(sl) == len(pl):
                    entry["confidence"] = {s: round(v, 4) for s, v in zip(sl, pl)}
            except ValueError:
                pass
        if "confidence" not in entry:
            pr = node["ann"].get(f"{trait_key}.prob")
            if pr:
                try: entry["confidence"] = {node["ann"][trait_key]: round(float(pr), 4)}
                except ValueError: pass
        attrs[trait_name] = entry

    post = node["ann"].get("posterior")
    if post:
        try: attrs["posterior"] = {"value": round(float(post), 4)}
        except ValueError: pass

    if node["label"]:                                # tip
        out["name"] = node["label"]
        m = meta.get(node["label"], {})
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
    if node["children"]:
        out["children"] = [
            convert(nodes, c, most_recent, meta, trait_key, trait_name, div)
            for c in node["children"]]
    return out


convert.counter = 0


def gate_summary(paths):
    """
    Verdicts from whichever gate JSONs were given, as (name, verdict, note).

    A published tree is read by people who did not run the pipeline. A trait
    colouring that failed its signal test renders identically to one that
    passed, so the verdicts have to travel with the build or they are invisible
    exactly where they matter most.
    """
    out = []
    for p in paths or []:
        p = Path(p)
        if not p.is_file():
            continue
        try:
            doc = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        if not isinstance(doc, dict) or "verdict" not in doc:
            continue
        name = p.stem.replace("_", " ")
        note = ""
        if isinstance(doc.get("reasons"), list) and doc["reasons"]:
            note = str(doc["reasons"][0])
        elif isinstance(doc.get("checks"), dict):
            errs = doc["checks"].get("errors") or []
            if errs:
                note = str(errs[0])
        out.append((name, str(doc["verdict"]), note))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mcc", required=True, type=Path)
    ap.add_argument("--config", type=Path,
                    help="Pathogen config. Title, maintainer, build URL, "
                         "colourings and the trait key come from here; "
                         "nothing naming the pathogen is hardcoded.")
    ap.add_argument("--metadata", type=Path,
                    help="Subsampled metadata with a 'label' column.")
    ap.add_argument("--full-metadata", type=Path,
                    help="metadata_clean.tsv, for country and species.")
    ap.add_argument("--alignment", type=Path,
                    help="Subsampled alignment. Used only to report how many "
                         "sequences the tree was built from.")
    ap.add_argument("--gates", nargs="*", type=Path, default=[],
                    help="Gate JSONs whose verdicts belong in the description.")
    ap.add_argument("--most-recent", type=float, default=None,
                    help="Decimal date of the most recent sample. Derived from "
                         "the tip labels when omitted, which is safer: a value "
                         "passed by hand that is a year stale shifts every "
                         "node together, so nothing looks wrong.")
    ap.add_argument("--trait-key", default=None,
                    help="Annotation key in the MCC tree. Auto-detected, or "
                         "taken from beast.discrete_trait.trait in the config.")
    ap.add_argument("--trait-name", default=None,
                    help="What to call the trait in Auspice.")
    ap.add_argument("--title", default=None)
    ap.add_argument("--maintainer", default=None)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    if not args.mcc.is_file():
        print(f"MCC tree not found: {args.mcc}", file=sys.stderr)
        return 1

    cfg = {}
    if args.config:
        import yaml
        cfg = yaml.safe_load(args.config.read_text()) or {}
    ns = cfg.get("nextstrain") or {}
    dt = (cfg.get("beast") or {}).get("discrete_trait") or {}

    trait_key = args.trait_key or dt.get("trait")
    trait_name = args.trait_name or dt.get("trait") or "host"

    try:
        nodes, resolved_tag = parse_mcc(args.mcc, trait_key)
    except (ValueError, SystemExit) as e:
        print(f"ERROR reading {args.mcc}: {e}", file=sys.stderr)
        return 1
    trait_key = resolved_tag
    root = 0
    dta_states = (cfg.get("hosts") or {}).get("dta_states") or {}
    meta = load_metadata(args.metadata, args.full_metadata, dta_states)

    tips = [n for n in nodes if n["is_tip"]]

    # The time axis. Derived from the labels rather than asserted, then checked
    # against every tip before anything is written.
    dated = [(t, label_date(t["label"])) for t in tips]
    dated = [(t, d) for t, d in dated if d is not None]
    if args.most_recent is not None:
        most_recent = args.most_recent
    elif dated:
        most_recent = max(d for _, d in dated)
        print(f"most recent tip date {most_recent:.4f} (from tip labels)")
    else:
        print("ERROR: no tip label carries a date and --most-recent was not "
              "given, so the time axis cannot be established.", file=sys.stderr)
        return 1

    bad = 0
    for t, date in dated:
        h = float(t["ann"].get("height", 0) or 0)
        if abs((most_recent - date) - h) > 0.05:
            bad += 1
    print(f"time axis check: {len(dated) - bad}/{len(dated)} tip heights consistent")
    if bad:
        print(f"ERROR: {bad} tips have heights inconsistent with their dates.\n"
              "The time axis is wrong -- every node would be plotted in the "
              "wrong place, and the figure would look entirely normal. Refusing "
              "to write. Check --most-recent against the tree's own tip dates.",
              file=sys.stderr)
        return 1

    convert.counter = 0
    tree = convert(nodes, root, most_recent, meta, trait_key, trait_name)

    # Colourings: whatever the config asks to colour by, plus what is always
    # available from the tree itself.
    wanted = list(ns.get("colour_by") or ns.get("color_by") or [])
    titles = {trait_name: "Host group", "country": "Country",
              "species": "Host species", "host_group": "Host group",
              "lineage": "Lineage", "region": "Region"}
    colorings = []
    seen = set()
    for key in [trait_name, *wanted]:
        k = trait_name if key == "host_group" else key
        if k in seen:
            continue
        seen.add(k)
        colorings.append({"key": k, "title": titles.get(k, k.replace("_", " ").title()),
                          "type": "categorical"})
    for key, title in (("posterior", "Clade posterior"), ("num_date", "Date")):
        colorings.append({"key": key, "title": title, "type": "continuous"})

    gates = gate_summary(args.gates)
    desc = [ns.get("description") or
            f"Time-scaled phylogeny with ancestral {trait_name} states "
            f"reconstructed under a discrete-trait model. Branches are coloured "
            f"by the reconstructed state; deep nodes carry wide uncertainty and "
            f"should not be read as confident assignments."]
    if gates:
        desc.append("")
        desc.append("**Gate verdicts for this build**")
        for name, verdict, note in gates:
            line = f"- {name}: **{verdict}**"
            if note:
                line += f" — {note}"
            desc.append(line)
        failed = [n for n, v, _ in gates if v in ("fail", "regression")]
        if failed:
            desc.append("")
            desc.append(
                f"One or more gates did not pass ({', '.join(failed)}). "
                "Quantities downstream of a failed gate are not supported by "
                "this analysis, whatever the figure shows.")

    doc = {
        "version": "v2",
        "meta": {
            "title": args.title or ns.get("title") or cfg.get("id") or "phylogeny",
            "updated": datetime.date.today().isoformat(),
            "panels": ["tree"],
            "colorings": colorings,
            "filters": [c["key"] for c in colorings if c["type"] == "categorical"],
            "display_defaults": {"color_by": trait_name,
                                 "branch_label": "none",
                                 "distance_measure": "num_date"},
            "description": "\n".join(desc),
        },
        "tree": tree,
    }
    if ns.get("build_url"):
        doc["meta"]["build_url"] = ns["build_url"]
    if ns.get("geo_resolutions"):
        doc["meta"]["geo_resolutions"] = [{"key": g} for g in ns["geo_resolutions"]]
    maint = args.maintainer or ns.get("maintainer")
    if maint:
        doc["meta"]["maintainers"] = [{"name": maint}]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=1))

    if args.alignment and args.alignment.is_file():
        n_seq = args.alignment.read_text().count(">")
        print(f"alignment: {n_seq} sequences")
    print(f"{len(tips)} tips, {convert.counter} internal nodes")
    print(f"root date {most_recent - float(nodes[root]['ann'].get('height', 0) or 0):.1f}")
    print(f"metadata joined for {sum(1 for t in tips if t['label'] in meta)} tips")
    for name, verdict, _ in gates:
        print(f"gate {name}: {verdict}")
    print(f"\nwrote {args.output}  ({args.output.stat().st_size / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

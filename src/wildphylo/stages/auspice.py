"""auspice stage: BEAST MCC tree to an Auspice dataset, then the entropy panel.

``wildphylo auspice dataset`` converts a BEAST MCC tree to an Auspice v2 JSON.
It does what ``augur import beast`` + ``augur export v2`` would do, without the
augur dependency chain (which pulls cvxopt and a compiler toolchain for features
this project does not use).

``wildphylo auspice entropy`` adds the entropy panel to an existing Auspice
dataset: per-node mutations from Fitch-parsimony ancestral reconstruction,
expressed in coordinate-reference nucleotide positions, plus amino-acid
mutations for every CDS in the dataset's genome annotation.

Everything that names the pathogen comes from the pathogen config
------------------------------------------------------------------
Title, maintainer, build URL, which attributes are colourable and the trait key
are read from ``nextstrain:`` and ``beast.discrete_trait:``. Nothing naming a
pathogen is hardcoded: in a pathogen-agnostic pipeline that would publish a BTV
build described as canine distemper, and the JSON is valid and renders, so
nothing downstream would catch it.

The time axis is checked before anything is written
---------------------------------------------------
Every tip's reconstructed height must agree with the date in its own label. If
they disagree the time axis is wrong, and a figure built from it is wrong in a
way that looks entirely normal. The step refuses to write rather than warning.
``--most-recent`` is derived from the tip labels rather than being required:
supply last year's value against this year's tree and every node moves
together, so nothing looks wrong.

Gate verdicts travel with the build
-----------------------------------
``--gates`` takes the JSON written by any gate and folds its verdict into the
Auspice description. A trait colouring that failed its signal test is
indistinguishable from one that passed unless the build says so.

Coordinates for the entropy step
--------------------------------
Auspice reads a CDS in *coordinate reference* nucleotide positions, not
alignment columns. An alignment span is only a valid CDS if the reference has no
gaps inside it; an insertion in any other sequence would make the span a
non-multiple of 3 and Auspice silently drops the feature. So ``nuc`` runs
1..len(ungapped reference), codons are assembled from the columns holding the
reference's own bases, and mutations at insertion columns (which have no
coordinate) are dropped and counted. This is also what augur does.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import datetime
import json
import logging
import re
import sys
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml

from wildphylo.cli import step
from wildphylo.traittrees import (
    TRANSLATE_ROW_RE,
    TREE_LINE_RE,
    NewickParser,
    Node,
    detect_tag,
    label_date,
    strip_leading_comment,
)

if TYPE_CHECKING:
    from collections.abc import Iterable

log = logging.getLogger(__name__)

# An Auspice tree node: heterogeneous JSON (name, node_attrs, children, ...).
AuspiceNode = dict[str, Any]

# --- auspice dataset -------------------------------------------------------


def parse_mcc(path: Path, tag: str | None = None) -> tuple[list[Node], str]:
    """Read a TreeAnnotator MCC tree into the flat node list, keeping every annotation.

    The parser is shared with the jump-history and state-through-time steps: an
    MCC tree and a posterior tree file are the same grammar, and a further
    parser would be a further place for the same quoting and nesting edge cases
    to be got wrong differently.
    """
    raw = path.read_text()
    translate: dict[str, str] = {}
    in_translate = False
    newick: str | None = None
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


def load_metadata(
    sub_path: Path | None,
    full_path: Path | None,
    dta_states: dict[str, str] | None = None,
) -> dict[str, dict[str, str]]:
    """Map each taxon label to its attributes, joining subsampled and full metadata."""
    meta: dict[str, dict[str, str]] = {}
    if not sub_path or not Path(sub_path).is_file():
        return meta
    with Path(sub_path).open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    extra: dict[str, dict[str, str]] = {}
    if full_path and Path(full_path).is_file():
        with Path(full_path).open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                extra[str(r.get("accession", "")).split(".")[0]] = r
    for r in rows:
        lab = r.get("label") or ""
        if not lab:
            continue
        acc = r.get("accession", "")
        e = extra.get(str(acc).split(".")[0], {})
        # The tip's host is collapsed through hosts.dta_states, exactly as the
        # BEAST XML step does when it builds the trait alignment. Without this
        # the tips carry the raw host group (wild_felid) while every internal
        # node carries the trait state (felid), so Auspice renders them as
        # different categories and the colouring implies a transition at every tip.
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


def convert(
    nodes: list[Node],
    idx: int,
    most_recent: float,
    meta: dict[str, dict[str, str]],
    trait_key: str,
    trait_name: str,
    counter: list[int],
    div: float = 0.0,
) -> AuspiceNode:
    """Recursively build the Auspice node dict from the flat node list.

    ``counter`` is a one-element list holding the number of internal nodes named
    so far; it is how internal nodes get unique ``NODE_nnnnnnn`` names.
    """
    node = nodes[idx]
    out: AuspiceNode = {}
    height = float(node["ann"].get("height", 0.0) or 0.0)
    num_date = most_recent - height
    div = div + (node["length"] or 0.0)

    attrs: dict[str, Any] = {"div": round(div, 6), "num_date": {"value": round(num_date, 4)}}

    lo = node["ann"].get("height_95%_HPD")
    if lo:
        try:
            a, b = (float(x) for x in re.split(r"[,\s]+", lo.strip()) if x)
            attrs["num_date"]["confidence"] = [
                round(most_recent - max(a, b), 4),
                round(most_recent - min(a, b), 4),
            ]
        except (ValueError, TypeError):
            pass

    if trait_key in node["ann"]:
        entry: dict[str, Any] = {"value": node["ann"][trait_key]}
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
                    entry["confidence"] = {s: round(v, 4) for s, v in zip(sl, pl, strict=True)}
            except ValueError:
                pass
        if "confidence" not in entry:
            pr = node["ann"].get(f"{trait_key}.prob")
            if pr:
                with contextlib.suppress(ValueError):
                    entry["confidence"] = {node["ann"][trait_key]: round(float(pr), 4)}
        attrs[trait_name] = entry

    post = node["ann"].get("posterior")
    if post:
        with contextlib.suppress(ValueError):
            attrs["posterior"] = {"value": round(float(post), 4)}

    if node["label"]:  # tip
        out["name"] = node["label"]
        m = meta.get(node["label"], {})
        for k in ("country", "species", "accession"):
            if m.get(k):
                attrs[k] = {"value": m[k]}
        if m.get("host"):
            attrs[trait_name] = {"value": m["host"]}  # tip state is observed, not inferred
    else:
        counter[0] += 1
        out["name"] = f"NODE_{counter[0]:07d}"

    out["node_attrs"] = attrs
    if node["children"]:
        out["children"] = [
            convert(nodes, c, most_recent, meta, trait_key, trait_name, counter, div)
            for c in node["children"]
        ]
    return out


def gate_summary(paths: Iterable[Path | str] | None) -> list[tuple[str, str, str]]:
    """Verdicts from whichever gate JSONs were given, as (name, verdict, note).

    A published tree is read by people who did not run the pipeline. A trait
    colouring that failed its signal test renders identically to one that
    passed, so the verdicts have to travel with the build or they are invisible
    exactly where they matter most. Missing and malformed files are skipped.
    """
    out: list[tuple[str, str, str]] = []
    for p_raw in paths or []:
        p = Path(p_raw)
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


def add_arguments_dataset(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo auspice dataset``."""
    parser.add_argument("--mcc", required=True, type=Path)
    parser.add_argument(
        "--config",
        type=Path,
        help="Pathogen config. Title, maintainer, build URL, "
        "colourings and the trait key come from here; "
        "nothing naming the pathogen is hardcoded.",
    )
    parser.add_argument("--metadata", type=Path, help="Subsampled metadata with a 'label' column.")
    parser.add_argument(
        "--full-metadata", type=Path, help="metadata_clean.tsv, for country and species."
    )
    parser.add_argument(
        "--alignment",
        type=Path,
        help="Subsampled alignment. Used only to report how many "
        "sequences the tree was built from.",
    )
    parser.add_argument(
        "--gates",
        nargs="*",
        type=Path,
        default=[],
        help="Gate JSONs whose verdicts belong in the description.",
    )
    parser.add_argument(
        "--most-recent",
        type=float,
        default=None,
        help="Decimal date of the most recent sample. Derived from "
        "the tip labels when omitted, which is safer: a value "
        "passed by hand that is a year stale shifts every "
        "node together, so nothing looks wrong.",
    )
    parser.add_argument(
        "--trait-key",
        default=None,
        help="Annotation key in the MCC tree. Auto-detected, or "
        "taken from beast.discrete_trait.trait in the config.",
    )
    parser.add_argument("--trait-name", default=None, help="What to call the trait in Auspice.")
    parser.add_argument("--title", default=None)
    parser.add_argument("--maintainer", default=None)
    parser.add_argument("--output", required=True, type=Path)


@step("auspice", "dataset", "Convert a BEAST MCC tree to an Auspice v2 dataset")
def dataset(args: argparse.Namespace) -> int:
    """Write the Auspice dataset; non-zero if the MCC tree or its time axis is unusable."""
    if not args.mcc.is_file():
        log.error("MCC tree not found: %s", args.mcc)
        return 1

    cfg: dict[str, Any] = {}
    if args.config:
        cfg = yaml.safe_load(args.config.read_text()) or {}
    ns = cfg.get("nextstrain") or {}
    dt = (cfg.get("beast") or {}).get("discrete_trait") or {}

    trait_hint = args.trait_key or dt.get("trait")
    trait_name = args.trait_name or dt.get("trait") or "host"

    read_error: str | None = None
    try:
        nodes, trait_key = parse_mcc(args.mcc, trait_hint)
    except (ValueError, SystemExit) as e:
        read_error = str(e)
    if read_error is not None:
        log.error("ERROR reading %s: %s", args.mcc, read_error)
        return 1
    root = 0
    dta_states = (cfg.get("hosts") or {}).get("dta_states") or {}
    meta = load_metadata(args.metadata, args.full_metadata, dta_states)

    tips = [n for n in nodes if n["is_tip"]]

    # The time axis. Derived from the labels rather than asserted, then checked
    # against every tip before anything is written.
    dated: list[tuple[Node, float]] = []
    for t in tips:
        d = label_date(t["label"])
        if d is not None:
            dated.append((t, d))
    if args.most_recent is not None:
        most_recent = args.most_recent
    elif dated:
        most_recent = max(d for _, d in dated)
        log.info("most recent tip date %.4f (from tip labels)", most_recent)
    else:
        log.error(
            "ERROR: no tip label carries a date and --most-recent was not "
            "given, so the time axis cannot be established."
        )
        return 1

    bad = 0
    for t, date in dated:
        h = float(t["ann"].get("height", 0) or 0)
        if abs((most_recent - date) - h) > 0.05:
            bad += 1
    log.info("time axis check: %d/%d tip heights consistent", len(dated) - bad, len(dated))
    if bad:
        log.error(
            "ERROR: %d tips have heights inconsistent with their dates.\n"
            "The time axis is wrong -- every node would be plotted in the "
            "wrong place, and the figure would look entirely normal. Refusing "
            "to write. Check --most-recent against the tree's own tip dates.",
            bad,
        )
        return 1

    counter = [0]
    tree = convert(nodes, root, most_recent, meta, trait_key, trait_name, counter)

    # Colourings: whatever the config asks to colour by, plus what is always
    # available from the tree itself.
    wanted = list(ns.get("colour_by") or ns.get("color_by") or [])
    titles = {
        trait_name: "Host group",
        "country": "Country",
        "species": "Host species",
        "host_group": "Host group",
        "lineage": "Lineage",
        "region": "Region",
    }
    colorings: list[dict[str, str]] = []
    seen: set[str] = set()
    for key in [trait_name, *wanted]:
        k = trait_name if key == "host_group" else key
        if k in seen:
            continue
        seen.add(k)
        colorings.append(
            {"key": k, "title": titles.get(k, k.replace("_", " ").title()), "type": "categorical"}
        )
    for key, title in (("posterior", "Clade posterior"), ("num_date", "Date")):
        colorings.append({"key": key, "title": title, "type": "continuous"})

    gates = gate_summary(args.gates)
    desc = [
        ns.get("description")
        or f"Time-scaled phylogeny with ancestral {trait_name} states "
        f"reconstructed under a discrete-trait model. Branches are coloured "
        f"by the reconstructed state; deep nodes carry wide uncertainty and "
        f"should not be read as confident assignments."
    ]
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
                "this analysis, whatever the figure shows."
            )

    doc: dict[str, Any] = {
        "version": "v2",
        "meta": {
            "title": args.title or ns.get("title") or cfg.get("id") or "phylogeny",
            "updated": datetime.date.today().isoformat(),
            "panels": ["tree"],
            "colorings": colorings,
            "filters": [c["key"] for c in colorings if c["type"] == "categorical"],
            "display_defaults": {
                "color_by": trait_name,
                "branch_label": "none",
                "distance_measure": "num_date",
            },
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
        log.info("alignment: %d sequences", n_seq)
    log.info("%d tips, %d internal nodes", len(tips), counter[0])
    log.info("root date %.1f", most_recent - float(nodes[root]["ann"].get("height", 0) or 0))
    log.info("metadata joined for %d tips", sum(1 for t in tips if t["label"] in meta))
    for name, verdict, _ in gates:
        log.info("gate %s: %s", name, verdict)
    log.info("wrote %s  (%.0f KB)", args.output, args.output.stat().st_size / 1024)
    return 0


# --- coordinate frame ------------------------------------------------------

GAPS = set("-.~")
STOPS = {"TAA", "TAG", "TGA"}
BASES = set("ACGT")
AMBIG = set("NnXx?-.~RYKMSWBDHVryknmswbdhv")

def _codon_table() -> dict[str, str]:
    bases = "TCAG"
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    table: dict[str, str] = {}
    k = 0
    for a in bases:
        for b in bases:
            for c in bases:
                table[a + b + c] = aas[k]
                k += 1
    return table


CODON = _codon_table()


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA file into {first word of header: upper-cased sequence}."""
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def pick_coordinate_reference(
    seqs: dict[str, str], name: str | None, auto: bool
) -> tuple[str, str]:
    """Return (name, sequence) of the coordinate reference; falls back to the least-gapped one."""
    if name and not auto:
        if name not in seqs:
            raise KeyError(
                f"reference {name!r} is not in the alignment. "
                f"Available (first 5): {list(seqs)[:5]}"
            )
        return name, seqs[name]
    chosen = min(seqs, key=lambda k: sum(1 for c in seqs[k] if c in GAPS))
    return chosen, seqs[chosen]


class CoordinateFrame:
    """Coordinate-reference coordinate system for one gapped reference sequence.

    Attributes:
        ref_length: Length of the reference with gaps removed.
        n_columns: Number of alignment columns.
        orf_start: 0-based offset of the ORF's first base within the ungapped
            reference (None until :meth:`find_orf` has run).
        n_codons: Codons in the ORF, excluding the terminal stop.
        col_for_refpos: ``col_for_refpos[p]`` is the 1-based alignment column
            holding reference base ``p`` (index 0 is unused).
        refpos_for_col: The reverse map; insertion columns are absent.
    """

    def __init__(self, ref_name: str, ref_seq: str) -> None:
        """Index the gapped reference sequence."""
        self.ref_name = ref_name
        self.ref = ref_seq.upper()
        self.n_columns = len(self.ref)

        self.col_for_refpos: list[int] = [0]
        self.refpos_for_col: dict[int, int] = {}
        n = 0
        for col, ch in enumerate(self.ref, start=1):
            if ch not in GAPS:
                n += 1
                self.col_for_refpos.append(col)
                self.refpos_for_col[col] = n
        self.ref_length = n
        self.ungapped = "".join(c for c in self.ref if c not in GAPS)

        self.orf_start: int | None = None
        self.orf_nt = 0
        self.n_codons = 0

    def find_orf(self) -> tuple[int | None, int]:
        """Find the longest ATG-initiated ORF in the ungapped reference.

        Many GenBank records are mRNAs or genome fragments with untranslated
        flanks, so the CDS rarely starts at base 1 and the frame cannot be
        assumed.
        """
        s = self.ungapped
        best_len = 0
        best_start: int | None = None
        for f in (0, 1, 2):
            i = f
            while i < len(s) - 2:
                if s[i : i + 3] == "ATG":
                    j = i
                    while j < len(s) - 2:
                        if s[j : j + 3] in STOPS:
                            break
                        j += 3
                    if j - i > best_len:
                        best_len, best_start = j - i, i
                    i = j + 3
                else:
                    i += 3
        self.orf_start = best_start
        self.orf_nt = best_len
        self.n_codons = best_len // 3
        return best_start, best_len

    def codon_columns(self, codon_index: int) -> tuple[int, int, int]:
        """Map a 0-based codon index within the ORF to its three alignment columns.

        Insertion columns are skipped, which is what keeps the reading frame
        correct across indels.
        """
        if self.orf_start is None:
            raise RuntimeError("call find_orf() first")
        p = self.orf_start + 3 * codon_index  # 0-based ungapped offset
        cols = self.col_for_refpos
        return cols[p + 1], cols[p + 2], cols[p + 3]

    def feature_codon_columns(self, a1: int, a2: int) -> list[tuple[int, int, int]]:
        """Codon-column triples for amino acids ``a1``..``a2`` (1-based, inclusive)."""
        return [self.codon_columns(i - 1) for i in range(a1, a2 + 1)]


def annotation_feature_spans(ann: dict[str, Any]) -> list[tuple[str, int, int]]:
    """CDS features from a genome_annotations object, as (name, start_nt, end_nt)."""
    out: list[tuple[str, int, int]] = []
    for name, f in ann.items():
        if name == "nuc":
            continue
        if "segments" in f:
            segs = f["segments"]
            out.append((name, int(segs[0]["start"]), int(segs[-1]["end"])))
        else:
            out.append((name, int(f["start"]), int(f["end"])))
    out.sort(key=lambda r: (r[1], r[2]))
    return out


def backup_if_overwriting(output: Path, enabled: bool = True) -> Path | None:
    """Copy ``output`` to a timestamped ``.bak.json`` before it is overwritten.

    The dataset is edited in place, so a failure part way through leaves a JSON
    with new annotations and stale mutations, which looks plausible in Auspice
    and is not. ``*.bak.json`` belongs in ``.gitignore``.
    """
    output = Path(output)
    if not enabled or not output.is_file():
        return None
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    dest = output.with_suffix(f".{stamp}.bak.json")
    dest.write_bytes(output.read_bytes())
    log.info("backed up existing %s -> %s", output, dest)
    return dest


# --- auspice entropy -------------------------------------------------------


def collect(node: AuspiceNode, out: list[AuspiceNode]) -> list[AuspiceNode]:
    """Append ``node`` and all its descendants to ``out`` (pre-order) and return it."""
    out.append(node)
    for c in node.get("children", []):
        collect(c, out)
    return out


def fitch_up(node: AuspiceNode, seqs: dict[str, str], n_sites: int) -> list[set[str]]:
    """Upward Fitch pass: build the candidate state set at every site."""
    if "children" not in node:
        s = seqs.get(node["name"])
        if s is None:
            return [set() for _ in range(n_sites)]
        return [({c} if c not in AMBIG else set()) for c in s[:n_sites]]

    child_sets = [fitch_up(c, seqs, n_sites) for c in node["children"]]
    mine: list[set[str]] = []
    for i in range(n_sites):
        sets = [cs[i] for cs in child_sets if cs[i]]
        if not sets:
            mine.append(set())
            continue
        inter = set.intersection(*sets)
        mine.append(inter or set.union(*sets))
    node["_fitch"] = mine
    return mine


def assign_down(
    node: AuspiceNode,
    parent_state: list[str | None] | None,
    n_sites: int,
    states: dict[str, list[str | None]],
    nuc_muts: dict[str, list[str]],
    frame: CoordinateFrame,
    counters: dict[str, int],
) -> list[str | None]:
    """Downward Fitch pass.

    Resolves states, records nucleotide mutations in coordinate-reference
    positions, and caches the resolved state for later translation.
    """
    state: list[str | None]
    if "children" not in node:
        s = node.get("_seq")
        # Ambiguities and gaps become None. Without this, a gap in a partial
        # sequence compares unequal to the parent base and is recorded as a
        # mutation, which manufactures huge numbers of spurious changes.
        state = [(c if c not in AMBIG else None) for c in s[:n_sites]] if s else [None] * n_sites
    else:
        f = node["_fitch"]
        state = []
        for i in range(n_sites):
            ps = parent_state[i] if parent_state else None
            if ps and ps in f[i]:
                state.append(ps)
            elif f[i]:
                state.append(min(f[i]))
            else:
                state.append(ps)

    muts: list[str] = []
    if parent_state:
        for i in range(n_sites):
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
        assign_down(c, state, n_sites, states, nuc_muts, frame, counters)
    return state


def parent_map(
    node: AuspiceNode, parent: str | None = None, out: dict[str, str | None] | None = None
) -> dict[str, str | None]:
    """Map each node name to its parent's name (None for the root)."""
    out = {} if out is None else out
    out[node["name"]] = parent
    for c in node.get("children", []):
        parent_map(c, node["name"], out)
    return out


def translate(state: list[str | None], triples: list[tuple[int, int, int]]) -> list[str]:
    """Translate a resolved state into amino acids, one per codon-column triple."""
    aas: list[str] = []
    for c1, c2, c3 in triples:
        codon = "".join(
            b if (b := state[c - 1]) is not None and b in BASES else "N" for c in (c1, c2, c3)
        )
        aas.append(CODON.get(codon, "X"))
    return aas


def final_checks(doc: AuspiceNode, locus_name: str) -> bool:
    """Verify the two invariants whose violation caused the original coordinate bug.

    A CDS that is not a multiple of 3 gets silently discarded by Auspice. And if
    the per-domain amino-acid counts do not sum to the whole-locus count, some
    feature is reading a different frame than the locus, which is exactly the
    failure alignment-column coordinates produced and is invisible in the
    rendered panel.
    """
    ok = True
    ann = doc["meta"]["genome_annotations"]

    bad = [k for k, v in ann.items() if k != "nuc" and (v["end"] - v["start"] + 1) % 3]
    if bad:
        log.error("CHECK FAILED: CDS not a multiple of 3: %s", bad)
        ok = False
    else:
        log.info("  all %d CDS lengths divisible by 3", len(ann) - 1)

    counts: Counter[str] = Counter()

    def walk(n: AuspiceNode) -> None:
        for k, v in n.get("branch_attrs", {}).get("mutations", {}).items():
            counts[k] += len(v)
        for c in n.get("children", []):
            walk(c)

    walk(doc["tree"])

    # Every annotation other than the locus itself (and the nucleotide track) is
    # a domain of the locus; the domain table lives in the pathogen repo.
    present = [d for d in ann if d not in (locus_name, "nuc")]
    if locus_name in ann and present:
        whole, parts = counts.get(locus_name, 0), sum(counts[d] for d in present)
        if whole == parts:
            log.info(
                "  %s amino acid mutations (%d) == sum over %d domains",
                locus_name,
                whole,
                len(present),
            )
        else:
            log.error(
                "CHECK FAILED: %s has %d amino acid mutations but the domains sum to %d. "
                "A feature is out of frame with the locus.",
                locus_name,
                whole,
                parts,
            )
            ok = False

    stale = [k for k in counts if k != "nuc" and k not in ann]
    if stale:
        log.error("CHECK FAILED: mutation keys with no matching feature: %s", stale)
        ok = False
    else:
        log.info("  no orphaned mutation keys")
    return ok


def add_arguments_entropy(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo auspice entropy``."""
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument(
        "--reference",
        help="Taxon label defining the coordinate system. "
        "Must match the one given to 11_annotate_H.py.",
    )
    parser.add_argument(
        "--auto-reference", action="store_true", help="Use the least-gapped sequence"
    )
    parser.add_argument("--gene", default="H", help=argparse.SUPPRESS)
    parser.add_argument("--gene-start", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        "--cds", nargs=2, type=int, metavar=("START", "END"), help=argparse.SUPPRESS
    )
    parser.add_argument("--cds-from-json", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not keep a timestamped copy of the output file before overwriting it.",
    )
    parser.add_argument("--output", required=True, type=Path)


@step("auspice", "entropy", "Add the entropy panel and per-node mutations to an Auspice dataset")
def entropy(args: argparse.Namespace) -> int:
    """Reconstruct ancestral states and write mutations; non-zero on a failed check."""
    if args.gene_start is not None or args.cds:
        log.warning(
            "note: --gene-start and --cds are ignored — CDS extents now come from "
            "the annotation written by 11_annotate_H.py"
        )
    if args.cds_from_json:
        log.warning("note: --cds-from-json is now the only behaviour and can be dropped")

    for f in (args.json, args.alignment):
        if not f.is_file():
            log.error("not found: %s", f)
            return 1

    doc = json.loads(args.json.read_text())
    seqs = read_fasta(args.alignment)
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        log.error(
            "alignment sequences differ in length %s — this file is not aligned",
            sorted(lengths)[:4],
        )
        return 1
    n_columns = lengths.pop()

    ann = doc["meta"].get("genome_annotations")
    if not ann or set(ann) <= {"nuc"}:
        log.error(
            "ERROR: the JSON has no CDS features. Run 11_annotate_H.py first — it "
            "establishes the coordinate system this step reconstructs into."
        )
        return 1

    # --- rebuild the same coordinate frame the annotate step used ---
    reference_error: str | None = None
    try:
        ref_name, ref_seq = pick_coordinate_reference(seqs, args.reference, args.auto_reference)
    except KeyError as e:
        reference_error = str(e)
    if reference_error is not None:
        log.error("%s", reference_error)
        return 1
    frame = CoordinateFrame(ref_name, ref_seq)
    frame.find_orf()
    log.info("reference: %s (%d nt ungapped, %d columns)", ref_name, frame.ref_length, n_columns)

    if frame.orf_start is None:
        log.error("ERROR: no ORF found in the reference.")
        return 1
    if ann["nuc"]["end"] != frame.ref_length:
        log.error(
            "ERROR: the annotation says nuc is %s nt but this reference is %d nt. "
            "The two steps were given different --reference or --alignment values.",
            ann["nuc"]["end"],
            frame.ref_length,
        )
        return 1

    # --- map each annotated CDS back to codon-column triples ---
    regions: list[tuple[str, int, int, list[tuple[int, int, int]]]] = []
    for name, start, end in annotation_feature_spans(ann):
        span = end - start + 1
        offset = start - frame.orf_start - 1
        if offset < 0 or offset % 3 or span % 3:
            log.error(
                "ERROR: feature %s (nt %d-%d) is not in frame with the reference ORF "
                "starting at %d. Re-run 11_annotate_H.py.",
                name,
                start,
                end,
                frame.orf_start + 1,
            )
            return 1
        a1 = offset // 3 + 1
        a2 = a1 + span // 3 - 1
        if a2 > frame.n_codons:
            log.error(
                "ERROR: feature %s extends to aa %d, past the %d aa ORF. Re-run 11_annotate_H.py.",
                name,
                a2,
                frame.n_codons,
            )
            return 1
        regions.append((name, a1, a2, frame.feature_codon_columns(a1, a2)))
    log.info("translating %d CDS feature(s) from the JSON annotation", len(regions))

    nodes = collect(doc["tree"], [])
    tips = [n for n in nodes if "children" not in n]
    matched = sum(1 for t in tips if t["name"] in seqs)
    log.info("%d tips in tree, %d in alignment, %d matched", len(tips), len(seqs), matched)
    if matched < len(tips) * 0.9:
        missing = [t["name"] for t in tips if t["name"] not in seqs][:3]
        log.error(
            "ERROR: fewer than 90%% of tips matched. Tip names must be identical "
            "between the tree and the alignment.\n  unmatched examples: %s",
            missing,
        )
        return 1

    for t in tips:
        t["_seq"] = seqs.get(t["name"])

    # --- ancestral reconstruction, once ---
    log.info("reconstructing ancestral states over %d sites...", n_columns)
    sys.setrecursionlimit(10000)
    fitch_up(doc["tree"], seqs, n_columns)
    states: dict[str, list[str | None]] = {}
    nuc_muts: dict[str, list[str]] = {}
    counters = {"dropped": 0}
    assign_down(doc["tree"], None, n_columns, states, nuc_muts, frame, counters)

    total = sum(len(v) for v in nuc_muts.values())
    on_branches = sum(1 for v in nuc_muts.values() if v)
    per_branch = total / on_branches if on_branches else 0
    log.info(
        "%d nucleotide mutations assigned across %d branches (%.1f per branch)",
        total,
        on_branches,
        per_branch,
    )
    if counters["dropped"]:
        log.info(
            "%d mutation(s) dropped: they fall in alignment columns that are insertions "
            "relative to the reference and so have no coordinate in this dataset",
            counters["dropped"],
        )

    # A branch cannot plausibly change more sites than exist, and in practice
    # should change far fewer. This catches gap/ambiguity handling errors.
    n_var = sum(
        1 for i in range(n_columns) if len({s[i] for s in seqs.values() if s[i] in BASES}) > 1
    )
    if per_branch > n_var * 0.25:
        log.warning(
            "WARNING: %.0f mutations per branch against only %d variable sites. That is "
            "implausible and usually means gaps or ambiguity codes are being treated as "
            "real bases.",
            per_branch,
            n_var,
        )

    # --- amino acid mutations, per feature ---
    pmap = parent_map(doc["tree"])
    aa_muts: dict[str, dict[str, list[str]]] = {}
    for name, a1, a2, triples in regions:
        cache: dict[str, list[str]] = {}
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
            ms = [
                f"{x}{i + 1}{y}"
                for i, (x, y) in enumerate(zip(a, b, strict=True))
                if x != y and x not in "X*" and y not in "X*"
            ]
            if ms:
                aa_muts.setdefault(node_name, {})[name] = ms
                n_aa += len(ms)
        log.info(
            "  %-16s aa %3d-%-3d (%d codons): %d amino-acid mutations",
            name,
            a1,
            a2,
            a2 - a1 + 1,
            n_aa,
        )

    # --- write mutations, rebuilding the block from scratch ---
    known = {name for name, _, _, _ in regions}
    stale: set[str] = set()
    for n in nodes:
        old = n.get("branch_attrs", {}).get("mutations", {})
        stale |= {k for k in old if k != "nuc" and k not in known}
        fresh: dict[str, list[str]] = {}
        m = nuc_muts.get(n["name"], [])
        if m:
            fresh["nuc"] = m
        fresh.update(aa_muts.get(n["name"], {}))
        n.setdefault("branch_attrs", {})["mutations"] = fresh
        for k in ("_fitch", "_seq"):
            n.pop(k, None)
    if stale:
        log.info(
            "cleared %d stale mutation key(s) left by earlier runs: %s%s",
            len(stale),
            ", ".join(sorted(stale)[:6]),
            " ..." if len(stale) > 6 else "",
        )

    if "entropy" not in doc["meta"].get("panels", []):
        doc["meta"].setdefault("panels", []).append("entropy")

    # variable-site summary, for your own information
    var = 0
    for i in range(n_columns):
        col = Counter(s[i] for s in seqs.values() if s[i] not in AMBIG)
        if len(col) > 1:
            var += 1
    log.info("%d of %d columns are variable (%.1f%%)", var, n_columns, 100 * var / n_columns)

    log.info("checks:")
    if not final_checks(doc, args.gene):
        log.error("Not writing — see the failed check above.")
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    backup_if_overwriting(args.output, not args.no_backup)
    args.output.write_text(json.dumps(doc, indent=1))
    log.info("wrote %s (%.0f KB)", args.output, args.output.stat().st_size / 1024)
    log.info(
        "The entropy panel will now appear below the tree in Auspice, with an "
        "amino-acid view for every annotated feature."
    )
    return 0

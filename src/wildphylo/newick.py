"""newick.py -- the one parser for annotated Newick and NEXUS tree files.

Two analyses used to carry a parser each. The jump-history and
state-through-time steps (and the Auspice export) read ancestral-state trees
and needed the node annotations; the trait signal gate read thousands of
posterior trees and needed only the topology, as arrays. Two grammars for one
file format means two places for the same quoting and annotation-placement edge
cases to be got wrong differently, so there is now one tokenizer and one
recursive-descent parser here, and everything else is a view of its output.

Representation
--------------
The parser produces a flat **node list**: one dict per node, parents always
before their children, so a single forward pass accumulates depths. That is
the form `beast jumps` and `beast state-through-time` walk, and the only one
that can carry the per-node state, `.prob` and (optionally) every annotation.

The trait signal statistics need something else: integer arrays (`left`,
`right`, `postorder`) they can index with numpy across thousands of trees, with
multifurcations resolved to a binary tree. That is `Topology`, built from a
node list by `to_topology`. It is derived, not parsed a second time: the node
list is the representation, `Topology` is the array view the statistics
consume. Building it is a single O(nodes) pass per tree, which is small beside
the statistics' own per-node numpy work (R+1 labellings per node, per tree).

What makes this parser rather than Bio.Phylo
--------------------------------------------
BEAST writes node annotations as [&host_group="felid",host_group.prob=0.93],
and Bio.Phylo discards them. The annotation IS the data for the discrete-trait
steps; the topology alone says nothing about host. Bio.Phylo is also the wrong
tool for ten thousand trees: it builds a Python object per clade. Parsing
directly gets three things right that a regex over the whole string does not:

  * A ']' inside a quoted annotation value does not end the comment. BEAST
    writes set-valued annotations and quoted labels that contain brackets.
  * `annotation_value` anchors on start-of-string or a comma, so asking for
    `host_group` does not match `host_group.prob` or `host_group.set` and come
    back with a probability where a state was wanted.
  * Annotations may precede the colon, follow it, or both, depending on which
    logger wrote the file. All positions are absorbed.

On the .prob field
------------------
Each node carries a `prob` slot, filled from `<tag>.prob` when present. It is
usually absent, and that is a property of the XML rather than of the parser:
beastclassic's TreeWithTraitLogger with a single `<metadata>` reference logs
the modal reconstructed state only. Without it, every transition is counted at
equal weight regardless of how confident the reconstruction was at either end,
and a node split 35/33 between two states contributes exactly as much as one at
95%. See DR-010.

Reading tree files
------------------
`iter_trees` streams node lists from a NEXUS file one tree at a time (posterior
files get large; there is no reason to hold 2000 trees in memory). The trait
signal gate wants a burn-in fraction, thinning and a cap applied per file
before parsing, so `iter_topologies` collects the tree strings first, slices
them, and parses only the kept ones. Both sit on the same NEXUS scanner.
"""

from __future__ import annotations

import contextlib
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, TextIO

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

# One node of the flat list NewickParser builds. Values are heterogeneous (parent
# index, state string, float length, child list, optional annotation dict), so
# Any is the truthful value type.
Node = dict[str, Any]

__all__ = [
    "NUMBER_RE",
    "TRANSLATE_ROW_RE",
    "TREE_LINE_RE",
    "NewickParser",
    "Node",
    "Topology",
    "annotation_value",
    "detect_tag",
    "iter_topologies",
    "iter_trees",
    "label_date",
    "node_depths",
    "parse_all_annotations",
    "parse_newick",
    "parse_topology",
    "read_bracket",
    "read_posterior",
    "smart_open",
    "strip_comments",
    "strip_leading_comment",
    "to_topology",
]

# --------------------------------------------------------------------------- #
# Annotation scanning
# --------------------------------------------------------------------------- #

NUMBER_RE = re.compile(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?")
TREE_LINE_RE = re.compile(r"^\s*tree\s+\S+\s*(?:\[[^\]]*\]\s*)?=\s*(.*;)\s*$", re.IGNORECASE)
TRANSLATE_ROW_RE = re.compile(r"^\s*(\S+)\s+(.+?)[,;]?\s*$")


def smart_open(path: Path | str) -> TextIO:
    """Open a text file, transparently gunzipping a ``.gz`` one."""
    p = Path(path)
    if p.name.endswith(".gz"):
        return gzip.open(p, "rt")
    return p.open()


def read_bracket(s: str, i: int) -> tuple[str, int]:
    """Read a [...] comment starting at s[i] == '['. Returns (content, next_i).

    Quote-aware so that a ']' inside a quoted annotation value does not end it.
    """
    if s[i] != "[":
        raise ValueError(f"expected '[' at offset {i}")
    i += 1
    start = i
    in_quote = False
    while i < len(s):
        c = s[i]
        if c == '"':
            in_quote = not in_quote
        elif c == "]" and not in_quote:
            return s[start:i], i + 1
        i += 1
    raise ValueError("unterminated [ annotation ]")


def annotation_value(content: str, key: str) -> str | None:
    """Pull key=value out of a BEAST annotation body.

    Anchored on start-of-string or a comma so that `location=` does not also
    match `location.prob=` or `location.set=`.
    """
    pat = re.compile(
        r"(?:^|,)\s*&?\s*" + re.escape(key) + r'=\s*(?:"([^"]*)"|\{([^}]*)\}|([^,]+))'
    )
    m = pat.search(content)
    if not m:
        return None
    return (m.group(1) or m.group(2) or m.group(3) or "").strip()


ANN_PAIR = re.compile(
    r"(?:^|,)\s*&?\s*([A-Za-z_][\w.%+-]*)\s*=\s*"
    r'(?:"([^"]*)"|\{([^}]*)\}|([^,]+))'
)


def parse_all_annotations(content: str) -> dict[str, str]:
    """Every key=value in a BEAST annotation body, as strings.

    Set-valued entries ({1950,1972}) keep their braces stripped and are
    returned as the raw comma-joined text, because the caller decides whether
    a given key is a pair of HPD bounds or a list of states -- guessing here
    would turn a two-element HPD into a two-state set.
    """
    out: dict[str, str] = {}
    for key, q, brace, bare in ANN_PAIR.findall(content):
        out[key] = (q or brace or bare or "").strip()
    return out


def strip_leading_comment(newick: str) -> str:
    """Drop the [&R] / [&U] rooting comment that precedes the Newick string."""
    s = newick.lstrip()
    while s.startswith("["):
        _, j = read_bracket(s, 0)
        s = s[j:].lstrip()
    return s


def strip_comments(s: str) -> str:
    """Remove every [...] annotation block, keeping quoted labels outside them intact."""
    out: list[str] = []
    i, n = 0, len(s)
    quote: str | None = None
    while i < n:
        c = s[i]
        if quote:
            out.append(c)
            if c == quote:
                quote = None
        elif c in "'\"":
            quote = c
            out.append(c)
        elif c == "[":
            _, i = read_bracket(s, i)
            continue
        elif c != "]":
            out.append(c)
        i += 1
    return "".join(out)


def detect_tag(newick: str, preferred: str | None = None) -> str:
    """Find the annotation key that carries the discrete state.

    beastclassic writes whatever `tag=` was set on AncestralStateTreeLikelihood
    (default 'location'); BEAST 1 DTA writes the trait name. Pick the first key
    whose values look like labels rather than numbers-with-decimals.
    """
    if preferred:
        return preferred
    i = newick.find("[")
    seen: list[str] = []
    while i != -1 and len(seen) < 5:
        content, j = read_bracket(newick, i)
        seen.append(content)
        i = newick.find("[", j)
    for content in seen:
        for key, raw in re.findall(r'(?:^|,)\s*&?\s*([A-Za-z_][\w.]*)=\s*("?[^,"]*"?)', content):
            if key.endswith((".prob", ".set")) or "." in key:
                continue
            if key.lower() in {"rate", "length", "height", "posterior", "r"}:
                continue
            val = raw.strip('"')
            if NUMBER_RE.fullmatch(val) and "." in val:
                continue
            return str(key)
    raise SystemExit(
        "Could not auto-detect the state annotation key. Pass --tag explicitly "
        "(it is the `tag=` attribute on AncestralStateTreeLikelihood in the XML; "
        'this repo\'s DTA XML uses tag="location").'
    )


# --------------------------------------------------------------------------- #
# The parser
# --------------------------------------------------------------------------- #


class NewickParser:
    """Recursive-descent parser for BEAST-annotated Newick.

    Produces a flat node list. Parents are always created before their children,
    so a single forward pass over the list is enough to accumulate depths.
    """

    def __init__(self, s: str, tag: str | None, keep_annotations: bool = False) -> None:
        self.s = s
        self.i = 0
        # The annotation key carrying the discrete state. None means "topology
        # only": comments are skipped without being read, which is what the
        # trait signal gate wants over thousands of trees.
        self.tag = tag
        self.nodes: list[Node] = []
        # An MCC tree carries far more per node than the state: height, its
        # 95% HPD, clade posterior, rate. The posterior-tree analyses want one
        # key across thousands of trees and would waste memory keeping the
        # rest; the Auspice export wants all of it from a single tree. Same
        # grammar, so one parser with a flag rather than a second parser.
        self.keep_annotations = keep_annotations

    # -- low-level ---------------------------------------------------------- #
    def peek(self) -> str:
        while self.i < len(self.s) and self.s[self.i] in " \t\n\r":
            self.i += 1
        return self.s[self.i] if self.i < len(self.s) else ""

    def read_label(self) -> str | None:
        out = []
        if self.peek() in ("'", '"'):
            q = self.s[self.i]
            self.i += 1
            while self.i < len(self.s) and self.s[self.i] != q:
                out.append(self.s[self.i])
                self.i += 1
            self.i += 1
            return "".join(out)
        while self.i < len(self.s) and self.s[self.i] not in "[]():,;":
            out.append(self.s[self.i])
            self.i += 1
        return "".join(out).strip() or None

    def read_number(self) -> float:
        m = NUMBER_RE.match(self.s, self.i)
        if not m:
            return 0.0
        self.i = m.end()
        return float(m.group(0))

    def read_comments(self, idx: int) -> None:
        """Absorb every [...] comment at the cursor into node ``idx``."""
        while self.peek() == "[":
            content, self.i = read_bracket(self.s, self.i)
            self._absorb(idx, content)

    # -- grammar ------------------------------------------------------------ #
    def node(self, parent: int | None) -> int:
        idx = len(self.nodes)
        self.nodes.append(
            {
                "parent": parent,
                "state": None,
                "prob": None,
                "length": 0.0,
                "label": None,
                "is_tip": True,
                "children": [],
                **({"ann": {}} if self.keep_annotations else {}),
            }
        )
        if parent is not None:
            self.nodes[parent]["children"].append(idx)
        if self.peek() == "(":
            self.nodes[idx]["is_tip"] = False
            self.i += 1  # consume '('
            while True:
                self.node(idx)
                c = self.peek()
                if c == ",":
                    self.i += 1
                    continue
                if c == ")":
                    self.i += 1
                    break
                raise ValueError(f"expected ',' or ')' near offset {self.i}")

        label = self.read_label()
        if label:
            self.nodes[idx]["label"] = label

        # annotation may sit before the colon, after it, or both
        self.read_comments(idx)
        if self.peek() == ":":
            self.i += 1
            self.read_comments(idx)
            self.nodes[idx]["length"] = self.read_number()
            self.read_comments(idx)
        return idx

    def _absorb(self, idx: int, content: str) -> None:
        nd = self.nodes[idx]
        if self.keep_annotations:
            nd["ann"].update(parse_all_annotations(content))
        if self.tag is None:
            return
        if nd["state"] is None:
            val = annotation_value(content, self.tag)
            if val is not None:
                nd["state"] = val
        if nd["prob"] is None:
            p = annotation_value(content, self.tag + ".prob")
            if p is not None:
                with contextlib.suppress(ValueError):
                    nd["prob"] = float(p)

    def parse(self) -> list[Node]:
        self.node(None)
        return self.nodes


def parse_newick(
    newick: str, tag: str | None = None, keep_annotations: bool = False
) -> list[Node]:
    """Newick string -> flat node list. The one entry point to the parser.

    A leading [&R] / [&U] comment is dropped. ``tag`` names the annotation key
    that carries the discrete state; leave it None to read the topology only.
    """
    return NewickParser(strip_leading_comment(newick), tag, keep_annotations).parse()


# --------------------------------------------------------------------------- #
# The array view the trait signal statistics consume
# --------------------------------------------------------------------------- #


@dataclass
class Topology:
    """A rooted tree as flat arrays.

    postorder  node indices, children always before their parent
    left,right  child indices; -1 for a tip. Multifurcations are resolved into
                a left-deep chain of internal nodes with zero-length branches,
                which changes no statistic here: AI, Fitch PS and MC are all
                invariant to arbitrary resolution of a polytomy by zero-length
                edges only when the resolution introduces no state conflict --
                Fitch PS can differ on a polytomy, so `n_polytomies` is
                recorded and reported rather than silently absorbed.
    tip_index   label -> row in the state vector
    """

    postorder: np.ndarray
    left: np.ndarray
    right: np.ndarray
    tip_index: dict[str, int]
    n_tips: int
    n_polytomies: int = 0

    @property
    def n_nodes(self) -> int:
        return len(self.postorder)


def to_topology(nodes: list[Node]) -> Topology:
    """Node list -> Topology.

    Tips are numbered 0..n_tips-1 in order of first appearance (left to right);
    internal nodes follow. The postorder array is built as the nodes are
    created, since a depth-first walk emits children before parents by
    construction.
    """
    left: list[int] = []
    right: list[int] = []
    labels: list[str] = []
    postorder: list[int] = []
    n_poly = 0

    def new_node(lc: int, r: int, lab: str) -> int:
        left.append(lc)
        right.append(r)
        labels.append(lab)
        idx = len(left) - 1
        postorder.append(idx)
        return idx

    made = [-1] * len(nodes)  # node-list index -> index of the Topology node standing for it
    stack: list[tuple[int, bool]] = [(0, False)]
    while stack:
        i, expanded = stack.pop()
        nd = nodes[i]
        if nd["is_tip"]:
            if nd["label"] is None:
                raise ValueError("a tip has no label")
            made[i] = new_node(-1, -1, nd["label"])
        elif not expanded:
            stack.append((i, True))
            stack.extend((c, False) for c in reversed(nd["children"]))
        else:
            kids = [made[c] for c in nd["children"]]
            if len(kids) > 2:
                n_poly += 1
            # left-deep resolution of any multifurcation
            node = kids[0]
            for k in kids[1:-1]:
                node = new_node(node, k, "")
            made[i] = new_node(node, kids[-1], "")

    # Renumber so tips occupy 0..n_tips-1: the state vector is indexed by tip.
    is_tip = [lc == -1 for lc in left]
    n_tips = sum(is_tip)
    remap = np.empty(len(left), dtype=np.int64)
    t = 0
    i = n_tips
    for idx in range(len(left)):
        if is_tip[idx]:
            remap[idx] = t
            t += 1
        else:
            remap[idx] = i
            i += 1

    n = len(left)
    lft = np.full(n, -1, dtype=np.int64)
    rgt = np.full(n, -1, dtype=np.int64)
    for idx in range(n):
        if not is_tip[idx]:
            lft[remap[idx]] = remap[left[idx]]
            rgt[remap[idx]] = remap[right[idx]]

    tip_index = {labels[idx]: int(remap[idx]) for idx in range(n) if is_tip[idx]}
    if len(tip_index) != n_tips:
        raise ValueError("duplicate tip labels in tree")

    return Topology(
        postorder=remap[np.array(postorder, dtype=np.int64)],
        left=lft,
        right=rgt,
        tip_index=tip_index,
        n_tips=n_tips,
        n_polytomies=n_poly,
    )


def parse_topology(newick: str) -> Topology:
    """Newick string -> Topology; annotations are skipped, not read."""
    return to_topology(parse_newick(newick))


# --------------------------------------------------------------------------- #
# NEXUS tree files
# --------------------------------------------------------------------------- #


def _scan_trees(
    path: Path | str, plain_newick: bool = False
) -> Iterator[tuple[dict[str, str], str]]:
    """Yield (translate, newick) for each tree in a NEXUS file, streaming.

    The translate dict is shared and complete by the time the first tree is
    yielded (BEAST writes the Translate block ahead of the trees). With
    ``plain_newick`` a bare ``(...);`` line also counts as a tree.
    """
    translate: dict[str, str] = {}
    in_translate = False
    with smart_open(path) as fh:
        for line in fh:
            stripped = line.strip()
            low = stripped.lower()
            if low.startswith("translate"):
                in_translate = True
                continue
            if in_translate:
                if stripped.startswith(";") or low in {"tree", ""} or low.startswith("tree "):
                    in_translate = False
                else:
                    m = TRANSLATE_ROW_RE.match(stripped)
                    if m:
                        translate[m.group(1)] = m.group(2).strip().strip("'\"")
                    if stripped.endswith(";"):
                        in_translate = False
                    continue
            m = TREE_LINE_RE.match(stripped)
            if m:
                yield translate, m.group(1)
            elif plain_newick and stripped.startswith("(") and stripped.endswith(";"):
                yield translate, stripped


def _apply_translate(nodes: list[Node], translate: dict[str, str]) -> None:
    for nd in nodes:
        if nd["is_tip"] and nd["label"] in translate:
            nd["label"] = translate[nd["label"]]


def iter_trees(
    path: Path | str, tag: str | None = None, keep_annotations: bool = False
) -> Iterator[tuple[int, list[Node], str]]:
    """Yield (tree_index, node_list, tag) for each tree in a NEXUS tree file.

    Streams line by line; posterior tree files get large and there is no reason
    to hold 2000 trees in memory at once. The Translate block is applied to the
    tip labels, and every node gets an ``nchild`` count. If ``tag`` is None it
    is detected from the first tree.
    """
    tree_i = -1
    resolved_tag = tag
    for translate, raw in _scan_trees(path):
        newick = strip_leading_comment(raw)
        if resolved_tag is None:
            resolved_tag = detect_tag(newick, tag)
        nodes = NewickParser(newick, resolved_tag, keep_annotations).parse()
        for nd in nodes:
            nd["nchild"] = 0
        for nd in nodes:
            if nd["parent"] is not None:
                nodes[nd["parent"]]["nchild"] += 1
        _apply_translate(nodes, translate)
        tree_i += 1
        yield tree_i, nodes, resolved_tag


def iter_topologies(
    path: Path | str, burnin: float = 0.10, thin: int = 1, max_trees: int | None = None
) -> Iterator[Topology]:
    """Yield a Topology for each retained tree of a NEXUS .trees file or plain newick.

    Resolves the Translate block if present.

    burnin is a FRACTION of the trees in the file, applied per file. Callers
    passing several files must apply it per file, which is what read_posterior
    below does -- pooling first and then discarding a fraction would drop the
    whole of one chain rather than the head of each.
    """
    scanned = list(_scan_trees(path, plain_newick=True))
    newicks = [raw for _, raw in scanned]
    translate = scanned[-1][0] if scanned else {}

    n = len(newicks)
    start = round(burnin * n)
    kept = newicks[start :: max(1, thin)]
    if max_trees and len(kept) > max_trees:
        step = int(np.ceil(len(kept) / max_trees))
        kept = kept[::step][:max_trees]

    for nk in kept:
        nodes = parse_newick(nk)
        _apply_translate(nodes, translate)
        yield to_topology(nodes)


def read_posterior(
    paths: Sequence[Path | str], burnin: float = 0.10, thin: int = 1, max_trees: int | None = None
) -> list[Topology]:
    """Read several chains, applying burnin and thinning within each."""
    per_file = None if max_trees is None else max(1, max_trees // max(1, len(paths)))
    trees: list[Topology] = []
    for p in paths:
        trees.extend(iter_topologies(p, burnin=burnin, thin=thin, max_trees=per_file))
    if not trees:
        raise ValueError("no trees read; check burnin and the file format")
    return trees


# --------------------------------------------------------------------------- #
# tree arithmetic
# --------------------------------------------------------------------------- #


def node_depths(nodes: list[Node]) -> list[float]:
    """Root-to-node distance for every node (parents precede children)."""
    depth = [0.0] * len(nodes)
    for i, nd in enumerate(nodes):
        p = nd["parent"]
        depth[i] = 0.0 if p is None else depth[p] + nd["length"]
    return depth


def label_date(label: str | None) -> float | None:
    """Decimal year from ACCESSION|host|year tip labels."""
    if not label:
        return None
    parts = label.split("|")
    for token in reversed(parts):
        try:
            return float(token)
        except ValueError:
            continue
    return None

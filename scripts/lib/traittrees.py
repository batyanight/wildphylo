"""
traittrees.py — parsing BEAST ancestral-state tree files.

Why this is a library and not a script
--------------------------------------
Two analyses read these files: 08f_jump_history (what changed, and when) and
08g_state_through_time (how many lineages carried each state, year by year).
They need the same parser, and a parser that disagrees with itself between the
two would produce a transition count and a lineage share that cannot both be
true of the same posterior.

It previously lived inside the jump-history script, and the state-through-time
script reached in and borrowed it with importlib, because a module whose name
starts with a digit cannot be imported normally. That broke the moment the
scripts were renamed. The pipeline numbering is a human convention for reading
order; shared code should not depend on it.

What makes this parser rather than Bio.Phylo
--------------------------------------------
BEAST writes node annotations as [&host_group="felid",host_group.prob=0.93],
and Bio.Phylo discards them. The annotation IS the data here — the tree
topology on its own says nothing about host. So the parsing has to be
annotation-aware, which means doing it directly.

Three things it gets right that a regex over the whole string does not:

  * A ']' inside a quoted annotation value does not end the comment. BEAST
    writes set-valued annotations and quoted labels that contain brackets.
  * `annotation_value` anchors on start-of-string or a comma, so asking for
    `host_group` does not match `host_group.prob` or `host_group.set` and come
    back with a probability where a state was wanted.
  * Annotations may precede the colon, follow it, or both, depending on which
    logger wrote the file. All three positions are absorbed.

On the .prob field
------------------
Each node carries a `prob` slot, filled from `<tag>.prob` when present. It is
usually absent, and that is a property of the XML rather than of the parser:
beastclassic's TreeWithTraitLogger with a single `<metadata>` reference logs
the modal reconstructed state only. Without it, every transition is counted at
equal weight regardless of how confident the reconstruction was at either end,
and a node split 35/33 between two states contributes exactly as much as one at
95%. See DR-010.
"""

from __future__ import annotations

import gzip
import re

__all__ = [
    "smart_open", "read_bracket", "annotation_value", "NewickParser",
    "strip_leading_comment", "detect_tag", "iter_trees", "node_depths",
    "label_date", "NUMBER_RE", "TREE_LINE_RE", "TRANSLATE_ROW_RE",
]

# --------------------------------------------------------------------------- #
# NEXUS / annotated-Newick parsing
# --------------------------------------------------------------------------- #

NUMBER_RE = re.compile(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?")
TREE_LINE_RE = re.compile(
    r"^\s*tree\s+\S+\s*(?:\[[^\]]*\]\s*)?=\s*(.*;)\s*$", re.IGNORECASE
)
TRANSLATE_ROW_RE = re.compile(r"^\s*(\S+)\s+(.+?)[,;]?\s*$")


def smart_open(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_bracket(s, i):
    """Read a [...] comment starting at s[i] == '['. Returns (content, next_i).

    Quote-aware so that a ']' inside a quoted annotation value does not end it.
    """
    assert s[i] == "["
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


def annotation_value(content, key):
    """Pull key=value out of a BEAST annotation body.

    Anchored on start-of-string or a comma so that `location=` does not also
    match `location.prob=` or `location.set=`.
    """
    pat = re.compile(
        r'(?:^|,)\s*&?\s*' + re.escape(key) + r'=\s*(?:"([^"]*)"|\{([^}]*)\}|([^,]+))'
    )
    m = pat.search(content)
    if not m:
        return None
    return (m.group(1) or m.group(2) or m.group(3) or "").strip()


class NewickParser:
    """Recursive-descent parser for BEAST-annotated Newick.

    Produces a flat node list. Parents are always created before their children,
    so a single forward pass over the list is enough to accumulate depths.
    """

    def __init__(self, s, tag):
        self.s = s
        self.i = 0
        self.tag = tag
        self.nodes = []

    # -- low-level ---------------------------------------------------------- #
    def peek(self):
        while self.i < len(self.s) and self.s[self.i] in " \t\n\r":
            self.i += 1
        return self.s[self.i] if self.i < len(self.s) else ""

    def read_label(self):
        start = self.i
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
        return "".join(out).strip() or (None if self.i == start else None)

    def read_number(self):
        m = NUMBER_RE.match(self.s, self.i)
        if not m:
            return 0.0
        self.i = m.end()
        return float(m.group(0))

    # -- grammar ------------------------------------------------------------ #
    def node(self, parent):
        idx = len(self.nodes)
        self.nodes.append(
            {
                "parent": parent,
                "state": None,
                "prob": None,
                "length": 0.0,
                "label": None,
                "is_tip": True,
            }
        )
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
        if self.peek() == "[":
            content, self.i = read_bracket(self.s, self.i)
            self._absorb(idx, content)

        if self.peek() == ":":
            self.i += 1
            if self.peek() == "[":
                content, self.i = read_bracket(self.s, self.i)
                self._absorb(idx, content)
            self.nodes[idx]["length"] = self.read_number()
            if self.peek() == "[":
                content, self.i = read_bracket(self.s, self.i)
                self._absorb(idx, content)
        return idx

    def _absorb(self, idx, content):
        if self.nodes[idx]["state"] is None:
            val = annotation_value(content, self.tag)
            if val is not None:
                self.nodes[idx]["state"] = val
        if self.nodes[idx]["prob"] is None:
            p = annotation_value(content, self.tag + ".prob")
            if p is not None:
                try:
                    self.nodes[idx]["prob"] = float(p)
                except ValueError:
                    pass

    def parse(self):
        self.node(None)
        return self.nodes


def strip_leading_comment(newick):
    """Drop the [&R] / [&U] rooting comment that precedes the Newick string."""
    s = newick.lstrip()
    while s.startswith("["):
        _, j = read_bracket(s, 0)
        s = s[j:].lstrip()
    return s


def detect_tag(newick, preferred=None):
    """Find the annotation key that carries the discrete state.

    beastclassic writes whatever `tag=` was set on AncestralStateTreeLikelihood
    (default 'location'); BEAST 1 DTA writes the trait name. Pick the first key
    whose values look like labels rather than numbers-with-decimals.
    """
    if preferred:
        return preferred
    i = newick.find("[")
    seen = []
    while i != -1 and len(seen) < 5:
        content, j = read_bracket(newick, i)
        seen.append(content)
        i = newick.find("[", j)
    for content in seen:
        for key, raw in re.findall(r'(?:^|,)\s*&?\s*([A-Za-z_][\w.]*)=\s*("?[^,"]*"?)', content):
            if key.endswith(".prob") or key.endswith(".set") or "." in key:
                continue
            if key.lower() in {"rate", "length", "height", "posterior", "r"}:
                continue
            val = raw.strip('"')
            if NUMBER_RE.fullmatch(val) and "." in val:
                continue
            return key
    raise SystemExit(
        "Could not auto-detect the state annotation key. Pass --tag explicitly "
        "(it is the `tag=` attribute on AncestralStateTreeLikelihood in the XML; "
        "this repo's DTA XML uses tag=\"location\")."
    )


def iter_trees(path, tag=None):
    """Yield (tree_index, node_list) for each tree in a NEXUS tree file.

    Streams line by line; posterior tree files get large and there is no reason
    to hold 2000 trees in memory at once.
    """
    translate = {}
    in_translate = False
    tree_i = -1
    resolved_tag = tag
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
            if not m:
                continue
            newick = strip_leading_comment(m.group(1))
            if resolved_tag is None:
                resolved_tag = detect_tag(newick, tag)
            nodes = NewickParser(newick, resolved_tag).parse()
            for nd in nodes:
                nd["nchild"] = 0
            for nd in nodes:
                if nd["parent"] is not None:
                    nodes[nd["parent"]]["nchild"] += 1
            for nd in nodes:
                if nd["is_tip"] and nd["label"] in translate:
                    nd["label"] = translate[nd["label"]]
            tree_i += 1
            yield tree_i, nodes, resolved_tag


# --------------------------------------------------------------------------- #
# tree arithmetic
# --------------------------------------------------------------------------- #


def node_depths(nodes):
    """Root-to-node distance for every node (parents precede children)."""
    depth = [0.0] * len(nodes)
    for i, nd in enumerate(nodes):
        p = nd["parent"]
        depth[i] = 0.0 if p is None else depth[p] + nd["length"]
    return depth


def label_date(label):
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

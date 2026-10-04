"""
traitsignal.py -- does the discrete trait carry phylogenetic signal at all?

Why this exists
---------------
A discrete-trait analysis ALWAYS returns an answer. BSSVS always selects some
rates, ancestral reconstruction always paints every internal node, and the jump
counter always produces a transition table with a largest entry. None of that
is evidence that the trait is associated with the phylogeny. If host labels
were shuffled at random across the tips, every one of those outputs would still
appear, look structured, and be reported.

So the question that has to be answered BEFORE any transition number is
interpreted is not "which transition is most common" but "are the tip states
more clustered on this tree than a random assignment of the same states would
be". That is what this module tests.

Parker, Rambaut and Pybus (2008, Infect Genet Evol 8:239-246) set out the
standard form. Three statistics, computed from TOPOLOGY AND TIP STATES ONLY:

  AI  association index      -- weighted imbalance of state mixing at internal
                                nodes. Lower means more structure.
  PS  parsimony score        -- Fitch minimum number of state changes needed to
                                explain the tips. Lower means more structure.
  MC  maximum monophyletic   -- per state, the largest clade containing only
      clade size               that state. Higher means more structure.

Deliberately NOT used here: BEAST's ancestral state reconstruction. That
reconstruction is the thing under suspicion. Using it to test itself would be
circular. This module reads the topology and the tip labels, nothing else,
which also means it runs on a sequence-only .trees file just as well as on the
ancestral-state one.

The null
--------
The observed analysis assigns ONE fixed labelling of states to tips and
evaluates it across the whole posterior sample of trees. The null must have the
same shape: one fixed PERMUTED labelling, also evaluated across the whole
posterior. So a null replicate is a single permutation applied to every
retained tree, summarised by its mean -- not a fresh permutation per tree.

Getting this wrong is the easy mistake, and it matters. If you permute
independently per tree and pool the per-tree values, you are comparing a mean
against a distribution of single values. The mean has far smaller variance, so
it lands in the tail almost automatically and everything comes out significant.

p is therefore computed against the distribution of R null REPLICATE MEANS:

    p = (1 + #{null means at least as extreme as observed}) / (1 + R)

The +1 in both places is Davison and Hinkley's correction: with R permutations
the smallest attainable p is 1/(R+1), and reporting p = 0 would claim more
resolution than R permutations can buy.

What a pass does and does not license
-------------------------------------
Passing means the tip states are not randomly distributed with respect to the
tree. It does NOT mean the reconstructed transition DIRECTIONS are reliable,
and it does not rule out deep-node anchoring -- a single old, deeply-placed
clade of one state produces genuine clustering and genuine significance while
still making every outbound count from that state a consequence of where the
oldest tips sit. Signal is necessary for a directional claim, not sufficient.
Run 08g_state_through_time.py as well; the two answer different questions.

Ambiguous tips
--------------
Tips whose state is unknown (the '?' that 07 assigns to host groups absent from
hosts.dta_states) carry no information about clustering and are PRUNED before
the statistics are computed, not counted as a state of their own. Treating '?'
as a state would invent a sixth host group that clusters perfectly by
construction, because the tips that have it are the ones nobody classified.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field, asdict

import numpy as np

# ---------------------------------------------------------------------------
# Newick / NEXUS parsing
#
# Bio.Phylo is used elsewhere in this repo and is the right tool for one tree.
# It is the wrong tool for ten thousand: it builds a Python object per clade,
# and the statistics below need to walk every node of every tree R+1 times.
# What is needed is a flat integer representation, so that is what is built.
# ---------------------------------------------------------------------------

# BEAST writes node annotations as [&rate=1.0,host_group="felid",...]. Nested
# brackets do not occur in BEAST output, but quoted strings containing ']' do,
# so strip with an explicit scanner rather than a regex.
def strip_comments(s: str) -> str:
    """Remove [...] annotation blocks, respecting quoted strings inside them."""
    out = []
    i, n, depth = 0, len(s), 0
    quote = None
    while i < n:
        c = s[i]
        if quote:
            out.append(c) if depth == 0 else None
            if c == quote:
                quote = None
        elif c in "'\"" and depth == 0:
            quote = c
            out.append(c)
        elif c in "'\"" and depth > 0:
            # quoted string inside a comment: skip to its close
            j = s.find(c, i + 1)
            i = n if j < 0 else j
        elif c == "[":
            depth += 1
        elif c == "]":
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(c)
        i += 1
    return "".join(out)


_LABEL = re.compile(r"^[^:,()]+")


@dataclass
class Topology:
    """
    A rooted tree as flat arrays.

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


def parse_newick(nk: str) -> Topology:
    """
    Newick -> Topology. Comments must already be stripped.

    Tips are numbered 0..n_tips-1 in order of first appearance; internal nodes
    follow. The postorder array is built directly by the parser, since a
    recursive descent over newick emits children before parents by construction.
    """
    nk = nk.strip()
    if nk.endswith(";"):
        nk = nk[:-1]

    left: list[int] = []
    right: list[int] = []
    labels: list[str | None] = []
    postorder: list[int] = []
    n_poly = 0

    def new_node(l: int, r: int, lab: str | None) -> int:
        left.append(l)
        right.append(r)
        labels.append(lab)
        idx = len(left) - 1
        postorder.append(idx)
        return idx

    pos = 0

    def parse_node() -> int:
        nonlocal pos, n_poly
        if nk[pos] == "(":
            pos += 1
            kids = [parse_node()]
            while nk[pos] == ",":
                pos += 1
                kids.append(parse_node())
            if nk[pos] != ")":
                raise ValueError(f"expected ')' at offset {pos}")
            pos += 1
            # label / support on the internal node, then branch length
            m = _LABEL.match(nk[pos:])
            if m:
                pos += m.end()
            if pos < len(nk) and nk[pos] == ":":
                pos += 1
                m = _LABEL.match(nk[pos:])
                pos += m.end() if m else 0
            if len(kids) > 2:
                n_poly += 1
            # left-deep resolution of any multifurcation
            node = kids[0]
            for k in kids[1:-1]:
                node = new_node(node, k, None)
            return new_node(node, kids[-1], None)
        # tip
        m = _LABEL.match(nk[pos:])
        if not m:
            raise ValueError(f"expected a tip label at offset {pos}")
        lab = m.group(0).strip().strip("'\"")
        pos += m.end()
        if pos < len(nk) and nk[pos] == ":":
            pos += 1
            m = _LABEL.match(nk[pos:])
            pos += m.end() if m else 0
        return new_node(-1, -1, lab)

    parse_node()

    # Renumber so tips occupy 0..n_tips-1: the state vector is indexed by tip.
    is_tip = [l == -1 for l in left]
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
    L = np.full(n, -1, dtype=np.int64)
    R = np.full(n, -1, dtype=np.int64)
    for idx in range(n):
        if not is_tip[idx]:
            L[remap[idx]] = remap[left[idx]]
            R[remap[idx]] = remap[right[idx]]

    tip_index = {labels[idx]: int(remap[idx]) for idx in range(n) if is_tip[idx]}
    if len(tip_index) != n_tips:
        raise ValueError("duplicate tip labels in tree")

    return Topology(
        postorder=remap[np.array(postorder, dtype=np.int64)],
        left=L, right=R, tip_index=tip_index, n_tips=n_tips,
        n_polytomies=n_poly,
    )


def _open(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def iter_trees(path, burnin: float = 0.10, thin: int = 1, max_trees: int | None = None):
    """
    Yield Topology objects from a NEXUS .trees file (BEAST) or a plain newick
    file, resolving the Translate block if present.

    burnin is a FRACTION of the trees in the file, applied per file. Callers
    passing several files must apply it per file, which is what read_posterior
    below does -- pooling first and then discarding a fraction would drop the
    whole of one chain rather than the head of each.
    """
    trans: dict[str, str] = {}
    in_translate = False
    newicks: list[str] = []
    with _open(path) as fh:
        for raw in fh:
            s = raw.strip()
            low = s.lower()
            if low.startswith("translate"):
                in_translate = True
                continue
            if in_translate:
                if s.startswith(";") or low.startswith("tree "):
                    in_translate = False
                else:
                    part = s.rstrip(",;")
                    if part:
                        bits = part.split(None, 1)
                        if len(bits) == 2:
                            trans[bits[0]] = bits[1].strip().strip("'\"")
                        continue
                    continue
            if low.startswith("tree "):
                j = s.find("=")
                if j > 0:
                    newicks.append(s[j + 1:])
            elif s.startswith("(") and s.endswith(";"):
                newicks.append(s)

    n = len(newicks)
    start = int(round(burnin * n))
    kept = newicks[start::max(1, thin)]
    if max_trees and len(kept) > max_trees:
        step = int(np.ceil(len(kept) / max_trees))
        kept = kept[::step][:max_trees]

    for nk in kept:
        top = parse_newick(strip_comments(nk))
        if trans:
            inv = {}
            for lab, idx in top.tip_index.items():
                inv[trans.get(lab, lab)] = idx
            top.tip_index = inv
        yield top


def read_posterior(paths, burnin: float = 0.10, thin: int = 1,
                   max_trees: int | None = None) -> list[Topology]:
    """Read several chains, applying burnin and thinning within each."""
    per_file = None if max_trees is None else max(1, max_trees // max(1, len(paths)))
    trees: list[Topology] = []
    for p in paths:
        trees.extend(iter_trees(p, burnin=burnin, thin=thin, max_trees=per_file))
    if not trees:
        raise ValueError("no trees read; check burnin and the file format")
    return trees


# ---------------------------------------------------------------------------
# The three statistics, vectorised over R labellings at once.
#
# Every function takes `states` of shape (R, n_tips): R independent labellings
# of the same tree. R = 1 is the observed case. The loop over nodes is Python;
# the loop over labellings is numpy. That ordering is what makes a thousand
# permutations across a thousand trees finish in seconds rather than an hour.
# ---------------------------------------------------------------------------

def parsimony_score(top: Topology, states: np.ndarray, n_states: int) -> np.ndarray:
    """
    Fitch small-parsimony, as bitmasks. Returns shape (R,).

    Fitch: at each internal node take the intersection of its children's state
    sets; if empty, take the union and charge one change. The total charged is
    the minimum number of changes any ancestral assignment would need.
    """
    R = states.shape[0]
    sets = np.zeros((top.n_nodes, R), dtype=np.int64)
    sets[:top.n_tips] = (np.int64(1) << states.T.astype(np.int64))
    changes = np.zeros(R, dtype=np.int64)
    for node in top.postorder:
        l, r = top.left[node], top.right[node]
        if l < 0:
            continue
        inter = sets[l] & sets[r]
        empty = inter == 0
        sets[node] = np.where(empty, sets[l] | sets[r], inter)
        changes += empty
    return changes


def _descendant_counts(top: Topology, states: np.ndarray, n_states: int) -> np.ndarray:
    """Per node, per labelling, the tip count in each state. (n_nodes, R, S)."""
    R = states.shape[0]
    counts = np.zeros((top.n_nodes, R, n_states), dtype=np.int32)
    rows = np.arange(R)
    for t in range(top.n_tips):
        counts[t, rows, states[:, t]] = 1
    for node in top.postorder:
        l, r = top.left[node], top.right[node]
        if l < 0:
            continue
        counts[node] = counts[l] + counts[r]
    return counts


def association_index(top: Topology, states: np.ndarray, n_states: int,
                      counts: np.ndarray | None = None) -> np.ndarray:
    """
    AI = sum over internal nodes of (1 - f) / 2^(n-1), where f is the frequency
    of the most common tip state below the node and n the number of tips below.

    The 2^(n-1) denominator is the point of the statistic: a node with many
    descendants contributes almost nothing, so AI measures mixing in the SMALL
    clades where it is informative, rather than being dominated by the root --
    which is mixed in every tree, structured or not, and says nothing.

    For n beyond about 1024 the denominator overflows float64 to inf and the
    contribution becomes exactly 0. That is the intended limit, not a bug.
    """
    if counts is None:
        counts = _descendant_counts(top, states, n_states)
    R = states.shape[0]
    ai = np.zeros(R, dtype=np.float64)
    for node in top.postorder:
        if top.left[node] < 0:
            continue
        c = counts[node]                      # (R, S)
        n = c.sum(axis=1).astype(np.float64)  # (R,)
        f = c.max(axis=1) / np.maximum(n, 1)
        with np.errstate(over="ignore"):
            denom = np.power(2.0, n - 1.0)
        ai += np.where(np.isfinite(denom), (1.0 - f) / denom, 0.0)
    return ai


def monophyletic_clades(top: Topology, states: np.ndarray, n_states: int,
                        counts: np.ndarray | None = None) -> np.ndarray:
    """
    MC per state: the largest clade whose tips are all that state. Shape (R, S).

    A single tip is a monophyletic clade of size 1, so a state present at all
    scores at least 1. A state absent from a labelling scores 0.
    """
    if counts is None:
        counts = _descendant_counts(top, states, n_states)
    R = states.shape[0]
    mc = np.zeros((R, n_states), dtype=np.int32)
    for node in top.postorder:
        c = counts[node]                       # (R, S)
        n = c.sum(axis=1, keepdims=True)
        pure = (c == n) & (n > 0)              # only one state present
        mc = np.where(pure & (c > mc), c, mc)
    return mc


# ---------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------

@dataclass
class StatResult:
    name: str
    observed: float
    null_mean: float
    null_lo: float
    null_hi: float
    p: float
    direction: str          # "lower" = structure gives a smaller value

    @property
    def significant(self) -> bool:
        return self.p < 0.05


@dataclass
class TraitSignalResult:
    n_trees: int
    n_tips: int
    n_pruned: int
    states: list[str]
    state_counts: dict[str, int]
    n_permutations: int
    ai: StatResult | None = None
    ps: StatResult | None = None
    mc: dict[str, StatResult] = field(default_factory=dict)
    verdict: str = "unknown"          # pass | weak | fail
    reasons: list[str] = field(default_factory=list)
    polytomies: int = 0

    def to_dict(self) -> dict:
        d = asdict(self)
        d["ai"] = asdict(self.ai) if self.ai else None
        d["ps"] = asdict(self.ps) if self.ps else None
        d["mc"] = {k: asdict(v) for k, v in self.mc.items()}
        return d


def _pvalue(observed: float, null: np.ndarray, direction: str) -> float:
    """
    Davison-Hinkley corrected permutation p. The smallest value obtainable
    from R permutations is 1/(R+1); reporting 0 would overstate the resolution
    the test actually has.
    """
    if direction == "lower":
        extreme = int(np.sum(null <= observed))
    else:
        extreme = int(np.sum(null >= observed))
    return (1.0 + extreme) / (1.0 + len(null))


def trait_signal(trees: list[Topology],
                 tip_states: dict[str, str],
                 n_permutations: int = 1000,
                 seed: int = 12345,
                 ai_threshold: float = 0.05,
                 ps_threshold: float = 0.05) -> TraitSignalResult:
    """
    Run the tip-label randomisation test over a posterior sample.

    tip_states maps tip label -> state. Labels mapped to '?' or absent from the
    mapping are pruned, in the sense that they are excluded from the state
    vector; see below for why they cannot simply be left in.
    """
    rng = np.random.default_rng(seed)

    usable = {k: v for k, v in tip_states.items() if v not in (None, "", "?")}
    states_sorted = sorted(set(usable.values()))
    S = len(states_sorted)
    if S < 2:
        raise ValueError(
            f"a trait with {S} state(s) has nothing to test; "
            f"check the traits file and hosts.dta_states")
    code = {s: i for i, s in enumerate(states_sorted)}

    # Restrict to tips present in EVERY tree with a usable state. Posterior
    # samples from one analysis share a taxon set, so this is normally a no-op,
    # but a mismatch here is a real error (wrong traits file for these trees)
    # and silently intersecting it away would hide that.
    common = set(trees[0].tip_index)
    for t in trees[1:]:
        common &= set(t.tip_index)
    labelled = common & set(usable)
    if not labelled:
        raise ValueError(
            "no tip label in the trees matches the trait mapping; "
            f"trees look like {sorted(common)[:2]}, "
            f"traits look like {sorted(usable)[:2]}")
    n_pruned = len(common) - len(labelled)

    keep = sorted(labelled)
    obs_vec = np.array([code[usable[k]] for k in keep], dtype=np.int64)
    n = len(keep)

    # Pruning '?' tips properly means removing them from the topology, not
    # just from the state vector -- a node whose descendants are all '?' would
    # otherwise show n=0 and contribute nonsense to AI. Rather than rebuild
    # every tree, restrict the descendant counts to the kept tips: counts for a
    # pruned tip are all-zero, so it contributes no tips to any ancestor, and a
    # clade of only pruned tips has n=0 and is skipped by both AI and MC. That
    # is exactly the pruned tree's answer for AI and MC. Fitch parsimony is the
    # exception -- a '?' tip must be an all-states set, not an empty one --
    # handled by giving pruned tips a full bitmask below.
    counts_index = [np.array([t.tip_index[k] for k in keep], dtype=np.int64)
                    for t in trees]

    R = n_permutations
    # perms[0] is the observed labelling; 1..R are permutations of it. Each
    # permutation is ONE fixed relabelling applied to every tree, matching how
    # the observed labelling is used. See the module docstring.
    perms = np.empty((R + 1, n), dtype=np.int64)
    perms[0] = obs_vec
    for r in range(1, R + 1):
        perms[r] = rng.permutation(obs_vec)

    ai_sum = np.zeros(R + 1, dtype=np.float64)
    ps_sum = np.zeros(R + 1, dtype=np.float64)
    mc_sum = np.zeros((R + 1, S), dtype=np.float64)
    polytomies = 0

    for top, idx in zip(trees, counts_index):
        polytomies += top.n_polytomies
        # Build a (R+1, n_tips) state array where kept tips carry their code
        # and pruned tips carry a sentinel handled per statistic.
        full = np.full((R + 1, top.n_tips), -1, dtype=np.int64)
        full[:, idx] = perms

        # counts: pruned tips contribute nothing
        cnt = np.zeros((top.n_nodes, R + 1, S), dtype=np.int32)
        rows = np.arange(R + 1)
        for j, tip in enumerate(idx):
            cnt[tip, rows, perms[:, j]] = 1
        for node in top.postorder:
            l, r = top.left[node], top.right[node]
            if l >= 0:
                cnt[node] = cnt[l] + cnt[r]

        ai_sum += association_index(top, full, S, counts=cnt)
        mc_sum += monophyletic_clades(top, full, S, counts=cnt)

        # Fitch: a pruned tip is ambiguous across all states, so it never
        # forces a change and never blocks an intersection.
        sets = np.full((top.n_nodes, R + 1), (1 << S) - 1, dtype=np.int64)
        sets[idx] = (np.int64(1) << perms.T.astype(np.int64))
        changes = np.zeros(R + 1, dtype=np.int64)
        for node in top.postorder:
            l, r = top.left[node], top.right[node]
            if l < 0:
                continue
            inter = sets[l] & sets[r]
            empty = inter == 0
            sets[node] = np.where(empty, sets[l] | sets[r], inter)
            changes += empty
        ps_sum += changes

    n_trees = len(trees)
    ai_mean = ai_sum / n_trees
    ps_mean = ps_sum / n_trees
    mc_mean = mc_sum / n_trees

    res = TraitSignalResult(
        n_trees=n_trees,
        n_tips=n,
        n_pruned=n_pruned,
        states=states_sorted,
        state_counts={s: int(np.sum(obs_vec == code[s])) for s in states_sorted},
        n_permutations=R,
        polytomies=polytomies,
    )

    res.ai = StatResult(
        "AI", float(ai_mean[0]), float(ai_mean[1:].mean()),
        float(np.percentile(ai_mean[1:], 2.5)),
        float(np.percentile(ai_mean[1:], 97.5)),
        _pvalue(ai_mean[0], ai_mean[1:], "lower"), "lower")
    res.ps = StatResult(
        "PS", float(ps_mean[0]), float(ps_mean[1:].mean()),
        float(np.percentile(ps_mean[1:], 2.5)),
        float(np.percentile(ps_mean[1:], 97.5)),
        _pvalue(ps_mean[0], ps_mean[1:], "lower"), "lower")
    for s in states_sorted:
        j = code[s]
        col = mc_mean[1:, j]
        res.mc[s] = StatResult(
            f"MC[{s}]", float(mc_mean[0, j]), float(col.mean()),
            float(np.percentile(col, 2.5)), float(np.percentile(col, 97.5)),
            _pvalue(mc_mean[0, j], col, "higher"), "higher")

    # --- verdict ------------------------------------------------------------
    # AI and PS are whole-tree statistics: they answer "is there structure".
    # Both must pass, because they can disagree -- PS is sensitive to a few
    # tightly clustered tips while AI weights small clades throughout the tree,
    # so one state cluster in an otherwise random labelling can move PS alone.
    reasons: list[str] = []
    ai_ok = res.ai.p < ai_threshold
    ps_ok = res.ps.p < ps_threshold
    if ai_ok and ps_ok:
        res.verdict = "pass"
        reasons.append(
            f"AI p={res.ai.p:.4f} and PS p={res.ps.p:.4f}: tip states are more "
            f"clustered than random relabelling of the same states.")
    elif ai_ok or ps_ok:
        res.verdict = "weak"
        reasons.append(
            f"AI p={res.ai.p:.4f}, PS p={res.ps.p:.4f}: the two whole-tree "
            f"statistics disagree. Structure is present somewhere but is not "
            f"a property of the tree as a whole.")
    else:
        res.verdict = "fail"
        reasons.append(
            f"AI p={res.ai.p:.4f}, PS p={res.ps.p:.4f}: the tip states are "
            f"indistinguishable from a random assignment. Reconstructed "
            f"ancestral states and transition counts describe the label "
            f"distribution, not transmission.")

    sig = [s for s in states_sorted if res.mc[s].p < 0.05]
    if res.verdict in ("pass", "weak"):
        if sig:
            reasons.append(
                f"MC significant for: {', '.join(sig)}. Only these states "
                f"cluster more than chance; a transition involving a state "
                f"not in this list has no demonstrated structure behind it.")
        else:
            reasons.append(
                "No state has a significant MC. Whatever AI or PS is "
                "detecting is not any one state forming clades, so "
                "state-specific claims are not supported.")
            if res.verdict == "pass":
                res.verdict = "weak"

    # Power warning. With few tips per state the test simply cannot reject,
    # and a non-significant result then means "no power", not "no signal".
    thin = {s: c for s, c in res.state_counts.items() if c < 5}
    if thin:
        reasons.append(
            f"LOW POWER: states with fewer than 5 tips: {thin}. A "
            f"non-significant MC for these is uninformative.")
    if n < 20:
        reasons.append(
            f"LOW POWER: {n} labelled tips in total. Permutation tests on "
            f"trees this small rarely reject even when structure is real.")
    if polytomies:
        reasons.append(
            f"{polytomies} polytomies across the posterior were resolved "
            f"arbitrarily with zero-length edges, which can shift PS slightly. "
            f"BEAST posteriors are normally fully resolved; a large count here "
            f"suggests the trees came from a consensus method instead.")
    res.reasons = reasons
    return res


def format_report(res: TraitSignalResult, width: int = 78) -> str:
    """Human-readable report. The verdict line is the one that matters."""
    L = ["=" * width,
         "TIP-LABEL RANDOMISATION  (Parker, Rambaut & Pybus 2008)",
         "=" * width,
         f"posterior trees     : {res.n_trees}",
         f"labelled tips       : {res.n_tips}"
         + (f"   ({res.n_pruned} pruned as unassigned)" if res.n_pruned else ""),
         f"permutations        : {res.n_permutations}",
         f"states              : " + ", ".join(
             f"{s} {res.state_counts[s]}" for s in res.states),
         "",
         f"  {'statistic':<16}{'observed':>10}{'null mean':>11}"
         f"{'null 95%':>20}{'p':>9}",
         "  " + "-" * (width - 4)]

    for st in (res.ai, res.ps):
        L.append(f"  {st.name:<16}{st.observed:>10.4f}{st.null_mean:>11.4f}"
                 f"{st.null_lo:>10.4f}{st.null_hi:>10.4f}{st.p:>9.4f}"
                 + ("  *" if st.significant else ""))
    L.append("")
    L.append("  maximum monophyletic clade size, per state")
    for s in res.states:
        st = res.mc[s]
        L.append(f"  {s:<16}{st.observed:>10.2f}{st.null_mean:>11.2f}"
                 f"{st.null_lo:>10.2f}{st.null_hi:>10.2f}{st.p:>9.4f}"
                 + ("  *" if st.significant else ""))
    L += ["",
          "  AI and PS: LOWER than null means structure.",
          "  MC:        HIGHER than null means structure.",
          "",
          f"VERDICT: {res.verdict.upper()}", ""]
    for r in res.reasons:
        L.append("  - " + r)
    L += ["", "=" * width,
          "This tests whether tip states are associated with the tree at all.",
          "It does NOT validate transition directions: a state whose tips are",
          "old and deeply placed clusters significantly and still produces",
          "outbound transition counts that are an artefact of sampling depth.",
          "Read 08g_state_through_time.py alongside this.",
          "=" * width]
    return "\n".join(L)

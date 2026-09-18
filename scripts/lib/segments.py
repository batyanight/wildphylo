"""
segments.py — segment congruence and reassortment screening.

Why this exists
---------------
The CDV pipeline assumed one locus, one tree, one clock. For a segmented virus
that assumption is not merely imprecise, it is false in a way that manufactures
results:

  * Reassortment means segments have genuinely DIFFERENT topologies. Concatenating
    them averages over conflicting histories and produces a tree that matches no
    segment's actual ancestry.
  * Segments also have different TMRCAs. For BTV, Seg-10's diversity is far
    younger than Seg-6's -- fitting one clock to a concatenate splits the
    difference and misdates both.
  * A reassortant taxon sits in different clades depending on the segment. Under
    a discrete-trait host analysis, that single taxon generates a spurious
    "host transition" in every segment tree where its placement moved, and the
    transition is an artefact of reassortment, not of transmission.

So for a segmented virus, congruence is not a diagnostic to run afterwards. It
is a precondition: it decides whether concatenation is legitimate at all, and it
flags the taxa whose placement cannot be trusted in any downstream host
reconstruction.

What is computed
----------------
  1. Pairwise topological distance between segment trees (normalised
     Robinson-Foulds on the shared taxon set). High distance = decoupled
     histories = do not concatenate.
  2. Per-taxon incongruence: how much each taxon's set of neighbours changes
     between segment trees. High scorers are reassortant candidates.
  3. A concatenation verdict, with the reason stated.

Caveats worth carrying into methods
-----------------------------------
Robinson-Foulds is blunt. It counts split differences without weighting by
support, so a poorly-supported rearrangement counts the same as a
well-supported one. Segment trees built from short segments (BTV Seg-10 is
822 nt) will look incongruent partly from lack of resolution rather than from
real reassortment. `collapse_support` addresses this by ignoring
weakly-supported splits before comparing; use it, and report the threshold.

This screens for reassortment. It does not test it formally -- that needs a
coalescent-with-reassortment model. What it does is stop a naive concatenated
analysis from being run unnoticed, which is the common failure.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from itertools import combinations


@dataclass
class CongruenceResult:
    segments: list = field(default_factory=list)
    n_shared_taxa: int = 0
    pairwise_rf: dict = field(default_factory=dict)      # "seg2|seg6" -> normalised RF
    mean_rf: float = 0.0
    max_rf: float = 0.0
    reassortant_candidates: list = field(default_factory=list)
    concatenate: bool = False
    verdict: str = "unknown"          # congruent | mixed | decoupled
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def splits(tree, taxa: set[str], collapse_support: float | None = None) -> set[frozenset]:
    """
    Every non-trivial bipartition of an unrooted tree, as a frozenset of the
    smaller side. Restricted to `taxa` so trees with different sampling can be
    compared.

    collapse_support ignores clades whose support falls below the threshold,
    so that unresolved nodes are not counted as topological disagreement.
    """
    out: set[frozenset] = set()
    for clade in tree.get_nonterminals():
        if collapse_support is not None:
            conf = getattr(clade, "confidence", None)
            if conf is not None and float(conf) < collapse_support:
                continue
        side = {t.name for t in clade.get_terminals()} & taxa
        other = taxa - side
        if len(side) < 2 or len(other) < 2:
            continue
        out.add(frozenset(min(side, other, key=lambda s: (len(s), sorted(s)))))
    return out


def normalised_rf(tree_a, tree_b, taxa: set[str],
                  collapse_support: float | None = None) -> float:
    """
    Robinson-Foulds distance scaled to [0, 1]. 0 = identical topology on the
    shared taxa; 1 = no split in common.
    """
    sa = splits(tree_a, taxa, collapse_support)
    sb = splits(tree_b, taxa, collapse_support)
    total = len(sa) + len(sb)
    if total == 0:
        return 0.0
    return len(sa ^ sb) / total


def neighbour_sets(tree, taxa: set[str], k: int = 5) -> dict[str, set[str]]:
    """
    For each taxon, its k nearest other taxa by patristic distance. Comparing
    these across segments localises incongruence to individual taxa, which
    pairwise RF cannot do -- RF says the trees differ, not which tip moved.
    """
    terms = [t for t in tree.get_terminals() if t.name in taxa]
    out: dict[str, set[str]] = {}
    for t in terms:
        d = []
        for u in terms:
            if u is t:
                continue
            try:
                d.append((tree.distance(t, u), u.name))
            except Exception:
                continue
        d.sort()
        out[t.name] = {name for _, name in d[:k]}
    return out


def find_reassortant_candidates(trees: dict, taxa: set[str], k: int = 5,
                                threshold: float = 0.6) -> list[dict]:
    """
    Score each taxon by how unstable its neighbourhood is across segments.

    score = 1 - (mean Jaccard similarity of its neighbour set between all
                 pairs of segment trees)

    A taxon that keeps the same neighbours in every segment scores near 0. A
    reassortant, whose segments came from different parents, scores high.

    This is a screen, not a test. A short segment with poor resolution also
    scores high. Treat high scorers as "exclude from the host-transition
    analysis or justify keeping them", not as confirmed reassortants.
    """
    if len(trees) < 2:
        return []
    nbrs = {seg: neighbour_sets(t, taxa, k) for seg, t in trees.items()}
    scored = []
    for taxon in sorted(taxa):
        sims = []
        for a, b in combinations(sorted(trees), 2):
            na, nb = nbrs[a].get(taxon, set()), nbrs[b].get(taxon, set())
            union = na | nb
            if not union:
                continue
            sims.append(len(na & nb) / len(union))
        if not sims:
            continue
        score = 1 - sum(sims) / len(sims)
        if score >= threshold:
            scored.append({"taxon": taxon, "incongruence": round(score, 3)})
    return sorted(scored, key=lambda d: -d["incongruence"])


def assess_congruence(trees: dict, *, collapse_support: float | None = 70.0,
                      congruent_max_rf: float = 0.15,
                      decoupled_min_rf: float = 0.40,
                      neighbour_k: int = 5,
                      reassortant_threshold: float = 0.6) -> CongruenceResult:
    """
    Decide whether segment trees may be concatenated, and flag reassortant
    candidates.

    trees: {segment_name: Bio.Phylo tree}
    """
    res = CongruenceResult(segments=sorted(trees))

    if len(trees) < 2:
        res.verdict = "congruent"
        res.concatenate = True
        res.notes.append(
            "only one segment supplied; congruence cannot be assessed and "
            "concatenation is vacuous")
        return res

    taxon_sets = [{t.name for t in tr.get_terminals()} for tr in trees.values()]
    shared = set.intersection(*taxon_sets)
    res.n_shared_taxa = len(shared)

    if len(shared) < 4:
        res.verdict = "unknown"
        res.concatenate = False
        res.notes.append(
            f"only {len(shared)} taxa are present in every segment. Congruence "
            "cannot be assessed, so concatenation cannot be justified. For "
            "segmented viruses this usually means segments were sequenced from "
            "different isolates and cannot be joined at all")
        return res

    union = set().union(*taxon_sets)
    if len(shared) < 0.5 * len(union):
        res.notes.append(
            f"only {len(shared)}/{len(union)} taxa have all segments. "
            "Concatenating would either discard most of the data or introduce "
            "large blocks of missing sequence, which biases branch lengths")

    for a, b in combinations(sorted(trees), 2):
        rf = normalised_rf(trees[a], trees[b], shared, collapse_support)
        res.pairwise_rf[f"{a}|{b}"] = round(rf, 4)

    vals = list(res.pairwise_rf.values())
    res.mean_rf = round(sum(vals) / len(vals), 4)
    res.max_rf = round(max(vals), 4)

    res.reassortant_candidates = find_reassortant_candidates(
        trees, shared, neighbour_k, reassortant_threshold)

    if res.max_rf <= congruent_max_rf:
        res.verdict = "congruent"
        res.concatenate = True
        res.notes.append(
            f"segment topologies agree (max normalised RF {res.max_rf:.3f}); "
            "concatenation is defensible, though segments may still differ in "
            "rate and should keep separate clock partitions")
    elif res.max_rf >= decoupled_min_rf:
        res.verdict = "decoupled"
        res.concatenate = False
        res.notes.append(
            f"segment topologies are decoupled (max normalised RF {res.max_rf:.3f}). "
            "Concatenation would average over conflicting histories and produce "
            "a tree matching no segment's actual ancestry. Analyse segments "
            "separately, and treat any host transition that is not recovered in "
            "several segments as unsupported")
    else:
        res.verdict = "mixed"
        res.concatenate = False
        res.notes.append(
            f"segment topologies partly disagree (max normalised RF {res.max_rf:.3f}). "
            "Default is to analyse separately. Concatenating a congruent subset "
            "is defensible if the subset is chosen on the RF matrix and stated")

    if res.reassortant_candidates:
        names = ", ".join(c["taxon"] for c in res.reassortant_candidates[:5])
        res.notes.append(
            f"{len(res.reassortant_candidates)} reassortant candidates "
            f"(showing up to 5): {names}. Each sits among different relatives "
            "depending on the segment. In a discrete-trait host analysis such a "
            "taxon generates a host transition that reflects reassortment "
            "rather than transmission")

    if collapse_support is not None:
        res.notes.append(
            f"splits with support below {collapse_support} were ignored before "
            "comparison; report this threshold, since RF is sensitive to it and "
            "short segments lose resolution rather than gaining real conflict")
    return res

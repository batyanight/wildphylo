"""
Tests for segment congruence screening.

Built on synthetic trees with known structure, so the expected answer is known
independently of the implementation. The key behaviours: identical topologies
must pass, a swapped clade must fail, and a single reassortant taxon must be
localised rather than just registering as "the trees differ".
"""

import io
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lib.segments import (  # noqa: E402
    assess_congruence, normalised_rf, splits, find_reassortant_candidates,
)

from Bio import Phylo  # noqa: E402


def t(newick):
    return Phylo.read(io.StringIO(newick), "newick")


# Trees must be big enough for the thresholds to mean anything. With six taxa
# there are only three informative splits, so even a complete clade swap tops
# out at RF 0.33 -- below the 0.40 decoupling threshold. That is a property of
# the tree size, not of the method, and using tiny fixtures would have forced a
# threshold that is wrong for real data. Twelve taxa is the smallest size where
# these behave as they will in practice.

# Six pairs in two clades of three pairs each.
CONGRUENT = ("((((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1):2,"
             "(((G:1,H:1):1,(I:1,J:1):1):1,(K:1,L:1):1):2);")
# Fully decoupled: the same tree shape with taxa permuted so that no cherry
# survives. Exchanging whole subclades is NOT enough -- it preserves every tip
# pair and only moves the deep splits, giving RF 0.33. Real reassortment
# between distant lineages relabels the tree, and that is what this represents.
SWAPPED = ("((((A:1,G:1):1,(C:1,I:1):1):1,(E:1,K:1):1):2,"
           "(((B:1,H:1):1,(D:1,J:1):1):1,(F:1,L:1):1):2);")
# One tip relocated: partial disagreement.
ONE_MOVED = ("((((A:1,B:1):1,(C:1,L:1):1):1,(E:1,F:1):1):2,"
             "(((G:1,H:1):1,(I:1,J:1):1):1,(K:1,D:1):1):2);")
TAXA = set("ABCDEFGHIJKL")


# --- splits and RF ----------------------------------------------------------

def test_identical_trees_have_zero_rf():
    assert normalised_rf(t(CONGRUENT), t(CONGRUENT), TAXA, None) == 0.0


def test_swapped_clades_have_high_rf():
    assert normalised_rf(t(CONGRUENT), t(SWAPPED), TAXA, None) > 0.4


def test_splits_ignore_trivial_bipartitions():
    """A single tip versus the rest carries no topological information."""
    s = splits(t(CONGRUENT), TAXA, None)
    assert all(2 <= len(x) <= len(TAXA) - 2 for x in s)


def test_rf_is_symmetric():
    a, b = t(CONGRUENT), t(SWAPPED)
    assert normalised_rf(a, b, TAXA, None) == normalised_rf(b, a, TAXA, None)


def test_low_support_splits_can_be_collapsed():
    """
    A short segment gives poorly-supported splits. Counting those as real
    conflict would make every short segment look reassortant rather than
    merely unresolved.
    """
    weak = t("((((A:1,B:1)20:1,(C:1,D:1)15:1)30:1,(E:1,F:1)22:1)25:2,"
             "(((G:1,H:1)18:1,(I:1,J:1)12:1)28:1,(K:1,L:1)19:2)21:2);")
    assert len(splits(weak, TAXA, collapse_support=70.0)) == 0
    assert len(splits(weak, TAXA, collapse_support=None)) > 0


def test_unlabelled_nodes_are_kept_not_collapsed():
    """
    A node with no support annotation is not evidence of low support. Treating
    missing as zero would silently erase real structure from any tree whose
    builder does not write confidences.
    """
    bare = t(CONGRUENT)
    assert len(splits(bare, TAXA, collapse_support=70.0)) > 0


# --- the concatenation decision --------------------------------------------

def test_congruent_segments_may_be_concatenated():
    r = assess_congruence({"seg2": t(CONGRUENT), "seg6": t(CONGRUENT)},
                          collapse_support=None)
    assert r.verdict == "congruent" and r.concatenate is True


def test_decoupled_segments_must_not_be_concatenated():
    """
    This is the failure the module exists to prevent: concatenating segments
    with conflicting histories yields a tree matching none of them.
    """
    r = assess_congruence({"seg2": t(CONGRUENT), "seg6": t(SWAPPED)},
                          collapse_support=None)
    assert r.verdict == "decoupled" and r.concatenate is False
    assert any("conflicting histories" in n for n in r.notes)


def test_mixed_defaults_to_not_concatenating():
    r = assess_congruence({"seg2": t(CONGRUENT), "seg6": t(ONE_MOVED)},
                          collapse_support=None)
    assert r.concatenate is False


def test_single_segment_cannot_be_assessed():
    r = assess_congruence({"seg2": t(CONGRUENT)})
    assert r.concatenate is True
    assert any("vacuous" in n for n in r.notes)


def test_too_few_shared_taxa_refuses_to_concatenate():
    """
    Segments sequenced from different isolates share almost no taxa. Joining
    them is meaningless, and must not silently succeed.
    """
    a = t("((A:1,B:1):1,(C:1,D:1):1);")
    b = t("((W:1,X:1):1,(Y:1,Z:1):1);")
    r = assess_congruence({"seg2": a, "seg6": b}, collapse_support=None)
    assert r.concatenate is False
    assert r.n_shared_taxa == 0
    assert any("different isolates" in n for n in r.notes)


def test_partial_segment_coverage_warns():
    a = t(CONGRUENT)
    b = t("((((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1):2,(W:1,X:1):2);")
    r = assess_congruence({"seg2": a, "seg6": b}, collapse_support=None)
    assert any("missing sequence" in n for n in r.notes)


# --- localising reassortment ------------------------------------------------

def test_reassortant_taxon_is_localised_not_just_detected():
    """
    RF says the trees differ. It cannot say which tip moved. The neighbour-set
    score must name the taxon, because that is what decides which tips to drop
    from the host analysis.
    """
    stable = t("((((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1):1,(G:1,R:1):1);")
    moved = t("((((A:1,B:1):1,(C:1,D:1):1):1,(E:1,R:1):1):1,(G:1,F:1):1);")
    cands = find_reassortant_candidates(
        {"seg2": stable, "seg6": moved}, set("ABCDEFGR"), k=2, threshold=0.5)
    names = [c["taxon"] for c in cands]
    assert "R" in names
    assert "A" not in names and "B" not in names


def test_stable_taxa_score_zero():
    cands = find_reassortant_candidates(
        {"seg2": t(CONGRUENT), "seg6": t(CONGRUENT)}, TAXA, k=3, threshold=0.5)
    assert cands == []


def test_reassortant_candidates_reported_in_notes():
    r = assess_congruence({"seg2": t(CONGRUENT), "seg6": t(SWAPPED)},
                          collapse_support=None, neighbour_k=3,
                          reassortant_threshold=0.5)
    assert r.reassortant_candidates
    assert any("reflects reassortment" in n for n in r.notes)


def test_single_segment_yields_no_candidates():
    assert find_reassortant_candidates({"seg2": t(CONGRUENT)}, TAXA) == []


# --- reporting --------------------------------------------------------------

def test_support_threshold_is_always_reported():
    """RF is sensitive to this threshold; it must never be silent."""
    r = assess_congruence({"seg2": t(CONGRUENT), "seg6": t(CONGRUENT)},
                          collapse_support=70.0)
    assert any("support below 70" in n for n in r.notes)


def test_all_pairs_are_scored():
    trees = {"seg2": t(CONGRUENT), "seg6": t(SWAPPED), "seg10": t(ONE_MOVED)}
    r = assess_congruence(trees, collapse_support=None)
    assert len(r.pairwise_rf) == 3
    assert r.max_rf >= r.mean_rf


def test_result_is_json_serialisable():
    import json
    r = assess_congruence({"seg2": t(CONGRUENT), "seg6": t(SWAPPED)},
                          collapse_support=None)
    json.dumps(r.to_dict())

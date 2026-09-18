"""
Tests for the per-step verification layer.

A check that never fires provides no assurance. Each test here confirms the
check fires on the failure it exists to catch, and stays quiet otherwise.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lib.preflight import (  # noqa: E402
    check_fetch, check_curation, check_trait_states, check_clock_prior,
    check_convergence,
)

CFG = {
    "fetch": {"taxid": 11232},
    "virus": {"expected_rate_range": [1.0e-4, 1.0e-2]},
    "hosts": {
        "min_group_size": 5,
        "dta_states": {
            "procyonid": "procyonid", "wild_canid": "wild_canid",
            "mustelid": "mustelid", "wild_felid": "felid",
            "domestic_cat": "felid", "domestic_dog": "domestic_dog",
        },
        "dta_other": None,
    },
    "beast": {"min_ess": 200, "clock_rate_prior": 7.46e-4,
              "clock_rate_source": "Panzera et al.",
              "discrete_trait": {"symmetric": True}},
    "dates": {"imprecise_date_policy": "midpoint"},
    "nextstrain": {"update": {"max_new_accessions": 500}},
}

# The real clade-3 host composition.
CLADE3 = {"procyonid": 91, "wild_canid": 32, "mustelid": 23,
          "domestic_dog": 8, "wild_felid": 7, "ailurid": 1}


# --- fetch ------------------------------------------------------------------

def test_empty_fetch_is_fatal():
    c = check_fetch(CFG, 0)
    assert not c.ok and "no records" in c.errors[0]


def test_losing_records_between_builds_is_fatal():
    """GenBank rarely withdraws records; a shrinking dataset is a query bug."""
    c = check_fetch(CFG, 150, n_previous=162)
    assert not c.ok
    assert "LOST" in c.errors[0]


def test_bulk_submission_warns():
    c = check_fetch(CFG, 1200, n_previous=162)
    assert c.ok and any("max_new_accessions" in w for w in c.warnings)


def test_normal_growth_is_clean():
    c = check_fetch(CFG, 180, n_previous=162)
    assert c.ok and not c.warnings


# --- curation ---------------------------------------------------------------

def test_zero_retained_points_at_locus_aliases():
    c = check_curation(CFG, 900, 0, 0, {})
    assert not c.ok and "locus.aliases" in c.errors[0]


def test_very_low_retention_warns():
    c = check_curation(CFG, 900, 40, 0, {"no_H_gene": 700})
    assert c.ok and any("survived curation" in w for w in c.warnings)


def test_dominant_exclusion_reason_warns():
    c = check_curation(CFG, 900, 400, 0, {"no_date": 480})
    assert any("single exclusion reason" in w for w in c.warnings)


def test_pending_review_warns():
    c = check_curation(CFG, 900, 400, 12, {})
    assert any("manual review" in w for w in c.warnings)


# --- trait states -----------------------------------------------------------

def test_undeclared_state_map_warns():
    cfg = {**CFG, "hosts": {**CFG["hosts"], "dta_states": None}}
    c = check_trait_states(cfg, CLADE3)
    assert any("dta_states is unset" in w for w in c.warnings)


def test_declared_map_collapses_and_reports_rate_count():
    c = check_trait_states(CFG, CLADE3)
    # procyonid, wild_canid, mustelid, felid, domestic_dog = 5 states -> 10 rates
    assert any("5 trait states -> 10" in n for n in c.notes)


def test_unmapped_group_is_reported_not_silently_dropped():
    """ailurid is not in the map and dta_other is null: it must be named."""
    c = check_trait_states(CFG, CLADE3)
    assert any("ailurid" in n for n in c.notes)


def test_twelve_states_warns_about_rate_count():
    """
    The unreduced host table: 12 groups -> 66 rates against 162 tips. This is
    the design that left 8 of 10 rates indistinguishable and the root at p=0.25.
    """
    cfg = {**CFG, "hosts": {**CFG["hosts"], "dta_states": None}}
    groups = {f"g{i}": 13 for i in range(12)}
    c = check_trait_states(cfg, groups)
    assert any("66" in n for n in c.notes)
    assert any("order of magnitude" in w for w in c.warnings)


def test_thinness_is_judged_after_collapsing_not_before():
    """
    domestic_cat (2) collapses into 'felid' alongside wild_felid (7), giving 9
    -- above min_group_size. The check must not warn, because the state map has
    already fixed the thinness. Judging group sizes before the collapse would
    produce a warning about a problem that no longer exists.
    """
    c = check_trait_states(CFG, {**CLADE3, "domestic_cat": 2})
    assert not any("min_group_size" in w for w in c.warnings)


def test_thin_state_warns_when_collapse_does_not_rescue_it():
    """domestic_dog maps to its own state, so 2 tips stays 2 tips."""
    c = check_trait_states(CFG, {**CLADE3, "domestic_dog": 2})
    assert any("min_group_size" in w for w in c.warnings)


def test_single_state_is_fatal():
    cfg = {**CFG, "hosts": {**CFG["hosts"],
                            "dta_states": {"procyonid": "procyonid"}}}
    c = check_trait_states(cfg, {"procyonid": 91})
    assert not c.ok


def test_asymmetric_model_doubles_rates_and_warns():
    cfg = {**CFG, "beast": {**CFG["beast"],
                            "discrete_trait": {"symmetric": False}}}
    c = check_trait_states(cfg, CLADE3)
    assert any("5 trait states -> 20" in n for n in c.notes)
    assert any("asymmetric" in w for w in c.warnings)


# --- clock prior ------------------------------------------------------------

def test_prior_matching_the_data_is_clean():
    c = check_clock_prior(CFG, rtt_slope=8.54e-4)
    assert c.ok and not c.warnings


def test_prior_from_the_wrong_pathogen_is_fatal():
    """A 50x disagreement between prior and data must stop the run."""
    cfg = {**CFG, "beast": {**CFG["beast"], "clock_rate_prior": 4.0e-2}}
    c = check_clock_prior(cfg, rtt_slope=8.54e-4)
    assert not c.ok
    assert "different pathogen" in c.errors[0]


def test_moderate_prior_data_disagreement_warns_only():
    cfg = {**CFG, "beast": {**CFG["beast"], "clock_rate_prior": 4.0e-3}}
    c = check_clock_prior(cfg, rtt_slope=8.54e-4)
    assert c.ok and any("prior-influenced" in w for w in c.warnings)


def test_rate_outside_genome_class_range_warns():
    cfg = {**CFG, "beast": {**CFG["beast"], "clock_rate_prior": 1.0e-6,
                            "clock_rate_source": "x"}}
    c = check_clock_prior(cfg, rtt_slope=1.1e-6)
    assert any("genome class" in w for w in c.warnings)


def test_uncited_prior_warns():
    cfg = {**CFG, "beast": {**CFG["beast"], "clock_rate_source": None}}
    c = check_clock_prior(cfg, rtt_slope=8.54e-4)
    assert any("clock_rate_source" in w for w in c.warnings)


def test_null_prior_with_no_slope_is_fatal():
    cfg = {**CFG, "beast": {**CFG["beast"], "clock_rate_prior": None}}
    c = check_clock_prior(cfg, rtt_slope=None)
    assert not c.ok


def test_null_prior_with_slope_uses_the_data():
    cfg = {**CFG, "beast": {**CFG["beast"], "clock_rate_prior": None}}
    c = check_clock_prior(cfg, rtt_slope=8.54e-4)
    assert c.ok and any("root-to-tip slope" in n for n in c.notes)


# --- convergence ------------------------------------------------------------

def test_low_ess_is_fatal():
    c = check_convergence(CFG, {"posterior": 1200, "clockRate": 88}, 2)
    assert not c.ok and "clockRate" in c.errors[0]


def test_good_ess_passes():
    c = check_convergence(CFG, {"posterior": 1200, "clockRate": 640}, 2)
    assert c.ok


def test_disagreeing_chains_are_fatal():
    c = check_convergence(CFG, {"posterior": 900}, 2,
                          between_chain_overlap={"treeHeight": 0.2})
    assert not c.ok and "disagree" in c.errors[0]


def test_single_chain_warns():
    c = check_convergence(CFG, {"posterior": 900}, 1)
    assert any("single chain" in w for w in c.warnings)


def test_interval_policy_requires_tip_height_ess():
    """Sampled tip dates mix poorly; their ESS must be checked, not assumed."""
    cfg = {**CFG, "dates": {"imprecise_date_policy": "interval"}}
    c = check_convergence(cfg, {"posterior": 900, "clockRate": 600}, 2)
    assert any("tip-height ESS" in w for w in c.warnings)


def test_render_is_readable():
    c = check_fetch(CFG, 0)
    assert "FAIL" in c.render()


# --- segmented-pathogen config validation ----------------------------------

import yaml  # noqa: E402
from lib.config import validate  # noqa: E402


def seg_cfg(**over):
    c = {
        "id": "btv", "name": "BTV",
        "virus": {"segmented": True, "reassortment_expected": True,
                  "expected_rate_range": [1e-5, 1e-2]},
        "fetch": {"taxid": 40051},
        "loci": [
            {"name": "seg2", "aliases": ["VP2"], "analyse": True},
            {"name": "seg6", "aliases": ["VP5"], "analyse": True},
            {"name": "seg10", "aliases": ["NS3"], "analyse": True},
        ],
        "segments": {"concatenate": "never", "congruence_gate": True},
        "dates": {"min_year": 1940, "imprecise_date_policy": "interval"},
        "hosts": {"table": "config/host_groups_btv.tsv", "vector_policy": "annotate"},
        "tree": {"model": "MFP"},
        "beast": {"clock_rate_prior": None, "n_chains": 2, "seeds": [1, 2],
                  "discrete_trait": {"enabled": False}},
        "references": {"assign_from": "seg2"},
    }
    c.update(over)
    return c


def test_segmented_config_validates():
    e, _w = validate(seg_cfg())
    assert e == []


def test_locus_and_loci_together_is_fatal():
    c = seg_cfg(locus={"name": "H", "aliases": ["H"]})
    e, _w = validate(c)
    assert any("both 'locus' and 'loci'" in x for x in e)


def test_missing_both_locus_forms_is_fatal():
    c = seg_cfg()
    del c["loci"]
    e, _w = validate(c)
    assert any("'locus' (unsegmented) or 'loci'" in x for x in e)


def test_concatenating_reassorting_segments_is_fatal():
    """
    The specific failure this pathogen was chosen to exercise: joining segments
    with conflicting histories produces a tree matching none of them.
    """
    c = seg_cfg(segments={"concatenate": "always", "congruence_gate": True})
    e, _w = validate(c)
    assert any("averages over conflicting histories" in x for x in e)


def test_single_analysed_segment_warns_reassortment_invisible():
    c = seg_cfg()
    for l in c["loci"][1:]:
        l["analyse"] = False
    _e, w = validate(c)
    assert any("reassortment would go undetected" in x for x in w)


def test_two_segments_warns_about_pairwise_only():
    c = seg_cfg()
    c["loci"][2]["analyse"] = False
    _e, w = validate(c)
    assert any("odd one out" in x for x in w)


def test_serotype_assignment_must_name_a_real_locus():
    c = seg_cfg(references={"assign_from": "seg99"})
    e, _w = validate(c)
    assert any("assign_from" in x for x in e)


def test_duplicate_locus_names_fatal():
    c = seg_cfg()
    c["loci"][1]["name"] = "seg2"
    e, _w = validate(c)
    assert any("duplicate locus names" in x for x in e)


def test_vector_as_state_warns():
    c = seg_cfg(hosts={"table": "config/host_groups_btv.tsv",
                       "vector_policy": "state"})
    _e, w = validate(c)
    assert any("where sequences were sampled" in x for x in w)


def test_bad_vector_policy_fatal():
    c = seg_cfg(hosts={"table": "config/host_groups_btv.tsv",
                       "vector_policy": "ignore"})
    e, _w = validate(c)
    assert any("vector_policy" in x for x in e)


def test_per_locus_clock_prior_is_bounds_checked():
    c = seg_cfg()
    c["loci"][0]["clock_rate_prior"] = 5.0     # subs/site/yr, absurd
    e, _w = validate(c)
    assert any("seg2" in x and "outside" in x for x in e)


def test_per_locus_date_policy_validated():
    c = seg_cfg()
    c["loci"][0]["imprecise_date_policy"] = "average"
    e, _w = validate(c)
    assert any("imprecise_date_policy" in x for x in e)


def test_congruence_gate_off_warns():
    c = seg_cfg(segments={"concatenate": "never", "congruence_gate": False})
    _e, w = validate(c)
    assert any("reassortants will not be flagged" in x for x in w)

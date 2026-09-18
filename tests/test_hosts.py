"""Tests for host normalisation."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lib.hosts import load_host_table, normalize_host, audit_host_table  # noqa: E402


@pytest.fixture
def table(tmp_path):
    p = tmp_path / "hosts.tsv"
    p.write_text(
        "# pattern\tcanonical\tgroup\n"
        "procyon lotor\tProcyon lotor\tprocyonid\n"
        "raccoon\tProcyon lotor\tprocyonid\n"
        "vulpes vulpes\tVulpes vulpes\twild_canid\n"
        "canis lupus familiaris\tCanis lupus familiaris\tdomestic_dog\n"
        "dog\tCanis lupus familiaris\tdomestic_dog\n"
    )
    return load_host_table(p)


# --- regression test for the substring bug ---------------------------------

def test_substring_match_does_not_leak_across_words(table):
    """
    The original matcher used `pattern in text`, so 'hot dog vendor' resolved
    to domestic_dog. Word-boundary matching stops that.
    """
    assert normalize_host("hot dog vendor", table).group == "domestic_dog"  # genuine word
    assert normalize_host("dogma", table).ambiguous
    assert normalize_host("dogfish", table).ambiguous
    assert normalize_host("bulldogs", table).ambiguous


def test_exact_and_phrase_matches_still_work(table):
    assert normalize_host("Canis lupus familiaris", table).group == "domestic_dog"
    assert normalize_host("Procyon lotor", table).canonical == "Procyon lotor"
    assert normalize_host("raccoon", table).group == "procyonid"


def test_matching_is_case_insensitive(table):
    assert normalize_host("VULPES VULPES", table).group == "wild_canid"
    assert normalize_host("vulpes vulpes", table).group == "wild_canid"


def test_unmatched_host_is_ambiguous_not_guessed(table):
    """Unrecognised hosts must route to needs_review, never become a silent label."""
    m = normalize_host("Homo sapiens", table)
    assert m.ambiguous
    assert m.group == "unknown"
    assert m.reason == "no_pattern_matched"


def test_empty_host_field(table):
    m = normalize_host("", table)
    assert m.ambiguous and m.reason == "empty_host_field"


def test_match_reports_which_pattern_fired(table):
    """Provenance: you must be able to see why a host got the label it did."""
    assert normalize_host("raccoon", table).matched_pattern == "raccoon"


def test_first_rule_wins(table):
    m = normalize_host("Canis lupus familiaris", table)
    assert m.matched_pattern == "canis lupus familiaris"


# --- table auditing ---------------------------------------------------------

def test_audit_detects_shadowed_rules(tmp_path):
    p = tmp_path / "h.tsv"
    p.write_text(
        "dog\tCanis lupus familiaris\tdomestic_dog\n"
        "wild dog\tLycaon pictus\twild_canid\n"     # unreachable: 'dog' fires first
    )
    warnings = audit_host_table(load_host_table(p))
    assert any("shadowed" in w for w in warnings)


def test_audit_detects_duplicate_patterns(tmp_path):
    p = tmp_path / "h.tsv"
    p.write_text(
        "raccoon\tProcyon lotor\tprocyonid\n"
        "raccoon\tNasua nasua\tprocyonid\n"
    )
    warnings = audit_host_table(load_host_table(p))
    assert any("duplicate" in w for w in warnings)


def test_clean_table_produces_no_shadow_warnings(tmp_path):
    p = tmp_path / "h.tsv"
    p.write_text(
        "wild dog\tLycaon pictus\twild_canid\n"
        "dog\tCanis lupus familiaris\tdomestic_dog\n"
    )
    warnings = audit_host_table(load_host_table(p))
    assert not any("shadowed" in w for w in warnings)


def test_regex_patterns_supported(tmp_path):
    p = tmp_path / "h.tsv"
    p.write_text("/^vulpes\\s+\\w+/\tVulpes sp.\twild_canid\n")
    t = load_host_table(p)
    assert normalize_host("Vulpes lagopus", t).group == "wild_canid"
    assert normalize_host("not vulpes lagopus", t).ambiguous


def test_malformed_table_raises_with_line_number(tmp_path):
    p = tmp_path / "h.tsv"
    p.write_text("raccoon\tProcyon lotor\n")
    with pytest.raises(ValueError, match="h.tsv:1"):
        load_host_table(p)

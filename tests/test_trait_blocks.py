"""
Tests for lib.beastxml.trait_blocks.

The discrete trait was previously added by hand in BEAUti, which meant the
model actually run was recorded nowhere: not in the config, not in git. The
config declared symmetric while the script's printed instructions said
asymmetric, and nothing could tell you which had been clicked. These tests fix
the generated XML to the config.
"""

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from lib.beastxml import (  # noqa: E402
    ROBUST_EIGEN, n_trait_rates, trait_blocks,
)

STATES = ["domestic_dog", "mustelid", "procyonid", "wild_canid", "wild_felid"]
TIPS = {
    "ACC1|procyonid|2016.5": "procyonid",
    "ACC2|wild_felid|2015.8": "wild_felid",
    "ACC3|domestic_dog|2021.0": "domestic_dog",
}


def blocks(**kw):
    kw.setdefault("symmetric", True)
    kw.setdefault("bssvs", True)
    return trait_blocks("america2", "host", STATES, TIPS, **kw)


# --- rate counting ----------------------------------------------------------

@pytest.mark.parametrize("n,sym,expected", [
    (2, True, 1), (2, False, 2),
    (5, True, 10), (5, False, 20),
    (6, True, 15), (6, False, 30),
])
def test_rate_count(n, sym, expected):
    assert n_trait_rates(n, sym) == expected


def test_one_state_is_refused():
    with pytest.raises(ValueError, match="at least 2 states"):
        n_trait_rates(1, True)


def test_asymmetric_doubles_the_rates():
    assert blocks(symmetric=False)["n_rates"] == 2 * blocks(symmetric=True)["n_rates"]


# --- the symmetry flag reaches the XML --------------------------------------

def test_symmetric_flag_is_written_not_assumed():
    assert 'symmetric="true"' in blocks(symmetric=True)["likelihood"]
    assert 'symmetric="false"' in blocks(symmetric=False)["likelihood"]


def test_rate_dimension_matches_symmetry():
    for sym, dim in ((True, 10), (False, 20)):
        b = blocks(symmetric=sym)
        assert f'dimension="{dim}"' in b["state"]


def test_asymmetric_gets_robust_eigen():
    """Asymmetric matrices destabilise the default eigen-decomposition and the
    run dies at startup. The fix must be automatic, not remembered."""
    assert ROBUST_EIGEN in blocks(symmetric=False)["likelihood"]
    assert ROBUST_EIGEN not in blocks(symmetric=True)["likelihood"]


# --- BSSVS ------------------------------------------------------------------

def test_bssvs_adds_indicators_and_they_are_logged():
    b = blocks(bssvs=True)
    assert "rateIndicator.s:host" in b["state"]
    assert "nonZeroRates.s:host" in b["prior"]
    assert "BitFlipOperator" in b["operators"]
    assert 'spec="Poisson"' in b["prior"]
    assert "rateIndicator.s:host" in b["log"]


def test_bssvs_off_removes_every_trace():
    b = blocks(bssvs=False)
    for key in ("state", "prior", "operators", "log", "likelihood"):
        assert "rateIndicator" not in b[key], key
        assert "nonZeroRates" not in b[key], key


def test_poisson_offset_keeps_the_matrix_connected():
    """offset = n_states - 1 is the minimum number of rates that can connect
    every state. A smaller offset lets the prior return -Infinity at the start
    state and the chain never begins."""
    assert 'offset="4"' in blocks()["prior"]


def test_poisson_lambda_is_overridable():
    assert 'lambda="2"' in blocks(poisson_lambda=2.0)["prior"]


# --- trait values -----------------------------------------------------------

def test_every_tip_appears_in_the_trait_set():
    data = blocks()["data"]
    for lbl, state in TIPS.items():
        assert f"{lbl}={state}" in data


def test_unknown_state_is_ambiguous_not_state_zero():
    """A '?' mapping to every state is the difference between 'we don't know'
    and 'it's a dog'."""
    assert "?=0 1 2 3 4" in blocks()["data"]


def test_tip_state_outside_the_declared_set_is_refused():
    with pytest.raises(ValueError, match="absent from the declared states"):
        trait_blocks("america2", "host", STATES,
                     {"ACC9|pinniped|2011": "pinniped"})


def test_state_codes_are_deterministic():
    a = trait_blocks("america2", "host", STATES, TIPS)["data"]
    b = trait_blocks("america2", "host", STATES, TIPS)["data"]
    assert a == b


# --- structural validity ----------------------------------------------------

@pytest.mark.parametrize("sym", [True, False])
@pytest.mark.parametrize("bssvs", [True, False])
def test_each_fragment_is_well_formed_xml(sym, bssvs):
    b = blocks(symmetric=sym, bssvs=bssvs)
    for key in ("data", "likelihood", "treelog"):
        ET.fromstring(b[key])
    for key in ("state", "prior", "operators", "log"):
        ET.fromstring(f"<wrap>{b[key]}</wrap>")


def test_references_resolve_to_declared_ids():
    """Every @id referenced by the likelihood must be declared in the state or
    inside the likelihood itself — a dangling idref is a startup crash."""
    b = blocks(symmetric=False)
    declared = set()
    for key in ("state", "likelihood", "data"):
        for el in ET.fromstring(f"<wrap>{b[key]}</wrap>").iter():
            if "id" in el.attrib:
                declared.add(el.attrib["id"])
    declared.add("Tree.t:america2")
    declared.add("TaxonSet.america2")

    for key in ("likelihood", "operators", "prior"):
        for el in ET.fromstring(f"<wrap>{b[key]}</wrap>").iter():
            for v in el.attrib.values():
                if v.startswith("@"):
                    assert v[1:] in declared, f"{key}: dangling ref {v}"


# --- end to end through 07 --------------------------------------------------

import random  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402

import yaml  # noqa: E402

COMPOSITION = {"procyonid": 25, "mustelid": 14, "wild_felid": 10,
               "domestic_dog": 9, "wild_canid": 7}


@pytest.fixture(scope="module")
def mock_aln(tmp_path_factory):
    """Host composition of the real america2 clade, so rate counts match."""
    rng = random.Random(7)
    rows, i = [], 0
    for host, n in COMPOSITION.items():
        for _ in range(n):
            i += 1
            seq = "".join(rng.choice("ACGT") for _ in range(600))
            rows.append(f">ACC{i:05d}|{host}|{1992 + rng.random() * 31.6:.3f}\n{seq}\n")
    p = tmp_path_factory.mktemp("aln") / "america2.fasta"
    p.write_text("".join(rows))
    return p


def run07(tmp_path, cfg_path, aln, *extra):
    out = tmp_path / "out.xml"
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts/07_make_beast_xml.py"),
         "--config", str(cfg_path), "--aln", str(aln),
         "--out-xml", str(out), "--out-traits", str(tmp_path / "t.txt"),
         *extra],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return out.read_text(), r.stdout


def variant(tmp_path, **discrete_trait):
    c = yaml.safe_load((ROOT / "config/pathogen/cdv-1200.yaml").read_text())
    c["beast"]["discrete_trait"].update(discrete_trait)
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(c, sort_keys=False))
    return p


def test_generated_xml_is_wellformed_and_has_no_dangling_refs(tmp_path, mock_aln):
    xml, _ = run07(tmp_path, ROOT / "config/pathogen/cdv-1200.yaml", mock_aln)
    root = ET.fromstring(xml)
    ids = {e.attrib["id"] for e in root.iter() if "id" in e.attrib}
    refs = {m for e in root.iter() for v in e.attrib.values()
            for m in re.findall(r"@([^\s\"]+)", v)}
    assert not refs - ids


def test_config_symmetry_reaches_the_generated_xml(tmp_path, mock_aln):
    asym, _ = run07(tmp_path, variant(tmp_path, symmetric=False), mock_aln)
    assert 'symmetric="false"' in asym and 'dimension="20"' in asym
    sym, _ = run07(tmp_path, variant(tmp_path, symmetric=True), mock_aln)
    assert 'symmetric="true"' in sym and 'dimension="10"' in sym


def test_no_trait_flag_produces_a_sequence_only_xml(tmp_path, mock_aln):
    xml, _ = run07(tmp_path, ROOT / "config/pathogen/cdv-1200.yaml",
                   mock_aln, "--no-trait")
    assert "AncestralStateTreeLikelihood" not in xml
    assert "rateIndicator" not in xml


def test_drt_replicates_carry_no_trait(tmp_path, mock_aln):
    """The date-randomisation test asks about the clock rate. Paying for the
    trait across 20 replicates buys nothing."""
    xml, _ = run07(tmp_path, ROOT / "config/pathogen/cdv-1200.yaml", mock_aln,
                   "--replicate", "3", "--no-trait")
    assert "DATE-RANDOMISED REPLICATE 3" in xml
    assert "AncestralStateTreeLikelihood" not in xml


def test_thin_rate_support_is_warned_about(tmp_path, mock_aln):
    _, out = run07(tmp_path, variant(tmp_path, symmetric=False), mock_aln)
    assert "tips/rate" in out and "decisive BF support" in out


def test_unmapped_host_becomes_ambiguous_not_a_state(tmp_path, mock_aln):
    """A pinniped has no dta_states mapping. It must enter as '?', not as
    whichever state happens to sort first."""
    extra = mock_aln.read_text() + ">ACC99999|pinniped|2011.500\n" + "A" * 600 + "\n"
    p = tmp_path / "with_pinniped.fasta"
    p.write_text(extra)
    xml, out = run07(tmp_path, ROOT / "config/pathogen/cdv-1200.yaml", p)
    assert "pinniped" in out and "ambiguous" in out
    assert "ACC99999|pinniped|2011.500=?" in xml


def test_two_state_companion_config_generates(tmp_path, mock_aln):
    cfg = ROOT / "config/pathogen/cdv-1200-dogwild.yaml"
    if not cfg.exists():
        pytest.skip("two-state companion config not present")
    xml, out = run07(tmp_path, cfg, mock_aln)
    assert 'dimension="2"' in xml
    assert "rateIndicator" not in xml      # BSSVS off: nothing to select
    assert "2 states, 2 rates" in out

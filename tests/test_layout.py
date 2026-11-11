"""
The build layout is the single source of truth for where a build's files go.

The expected paths below are the ones existing builds (e.g.
``builds/cdv-smoke/20261003``) already use, so a failure here means a rename
that would orphan old builds, not just a test to update.
"""

from pathlib import Path

import pytest

from wildphylo.layout import BuildLayout

ROOT = Path("builds/cdv/20260101")


@pytest.fixture
def layout() -> BuildLayout:
    return BuildLayout(ROOT, "cdv", "host_group")


# (method name, call arguments, path relative to the build root)
CASES: list[tuple[str, tuple[str, ...], str]] = [
    ("provenance", (), "provenance.json"),
    ("log", ("01_fetch",), "logs/01_fetch.log"),
    ("raw_records", (), "raw/cdv.gb"),
    ("accession_snapshot", (), "raw/cdv.acc"),
    ("interim_dir", (), "interim"),
    ("metadata", (), "interim/metadata_clean.tsv"),
    ("needs_review", (), "interim/needs_review.tsv"),
    ("review_marker", (), "interim/.review_ok"),
    ("flagged", (), "interim/flagged.tsv"),
    ("exclusions", (), "interim/exclusions.tsv"),
    ("locus_sequences", ("H",), "interim/H.fasta"),
    ("anchors", ("H",), "interim/H_anchors.fasta"),
    ("alignment", ("H",), "aligned/H_aligned.fasta"),
    ("alignment_qc", ("H",), "qc/H_alignment_qc.tsv"),
    ("alignment_qc_plot", ("H",), "qc/H_coverage.png"),
    ("ml_tree", ("H",), "tree/H_ml.treefile"),
    ("ml_tree_prefix", ("H",), "tree/H_ml"),
    ("temporal_report", ("H",), "temporal/H_temporal_signal.json"),
    ("temporal_report_md", ("H",), "temporal/H_temporal_signal.md"),
    ("temporal_plot", ("H",), "temporal/H_root_to_tip.png"),
    ("congruence_report", (), "segments/congruence.json"),
    ("congruence_md", (), "segments/congruence.md"),
    ("reassortant_candidates", (), "segments/reassortant_candidates.tsv"),
    ("subsample", ("H",), "subsampled/H_sub.fasta"),
    ("subsample_metadata", ("H",), "subsampled/H_sub_metadata.tsv"),
    ("beast_locus_dir", ("H",), "beast/H"),
    ("beast_chain_dir", ("H", "12345"), "beast/H/seed12345"),
    ("beast_xml", ("H", "12345"), "beast/H/seed12345/cdv.xml"),
    ("beast_traits", ("H", "12345"), "beast/H/seed12345/traits.txt"),
    ("beast_log", ("H", "12345"), "beast/H/seed12345/cdv.log"),
    ("beast_trees", ("H", "12345"), "beast/H/seed12345/cdv.trees"),
    ("beast_trait_trees", ("H", "12345"), "beast/H/seed12345/cdv_host_group.trees"),
    ("convergence_report", ("H",), "beast/H/convergence.json"),
    ("combined_log", ("H",), "beast/H/combined.log"),
    ("combined_trees", ("H",), "beast/H/combined.trees"),
    ("mcc_tree", ("H",), "beast/H/mcc.tree"),
    ("trait_signal_report", ("H",), "beast/H/trait_signal.json"),
    ("trait_signal_md", ("H",), "beast/H/trait_signal.md"),
    ("drt_dir", ("H",), "drt/H"),
    ("drt_xml", ("H", "3"), "drt/H/rep3/cdv_drt.xml"),
    ("drt_log", ("H", "3"), "drt/H/rep3/cdv_drt.log"),
    ("drt_summary", ("H",), "drt/H/drt_summary.json"),
    ("drt_plot", ("H",), "drt/H/drt.png"),
    ("auspice_dataset", ("cdv-phylodynamics",), "auspice/cdv-phylodynamics.json"),
    ("auspice_dataset", ("btv", "seg2"), "auspice/btv_seg2.json"),
    ("build_comparison", (), "build_comparison.md"),
    ("build_comparison_json", (), "build_comparison.json"),
]


@pytest.mark.parametrize(("method", "args", "relative"), CASES, ids=str)
def test_method_returns_the_expected_relative_path(layout, method, args, relative):
    assert getattr(layout, method)(*args) == ROOT / relative


def test_every_public_method_is_covered_by_a_case():
    """A new path added to the layout without a pinned expectation is a failure."""
    public = {
        name
        for name in dir(BuildLayout)
        if not name.startswith("_")
        and callable(getattr(BuildLayout, name))
        and name not in {"from_config", "at"}
    }
    # pathogen_dir sits above the build root, so it has its own test below.
    assert public - {"pathogen_dir"} == {method for method, _, _ in CASES}


def test_snakemake_wildcards_pass_through(layout):
    """Wildcard patterns are plain strings, so rules can use them as-is."""
    assert str(layout.alignment("{locus}")) == str(ROOT / "aligned" / "{locus}_aligned.fasta")
    assert str(layout.beast_log("{locus}", "{seed}")) == str(
        ROOT / "beast" / "{locus}" / "seed{seed}" / "cdv.log"
    )
    assert str(layout.drt_xml("{locus}", "{rep}")) == str(
        ROOT / "drt" / "{locus}" / "rep{rep}" / "cdv_drt.xml"
    )


def test_pathogen_dir_is_the_parent_of_the_build(layout):
    assert layout.pathogen_dir() == Path("builds/cdv")


def test_trait_trees_need_a_trait():
    no_trait = BuildLayout(ROOT, "cdv", None)
    with pytest.raises(ValueError, match="no discrete trait"):
        no_trait.beast_trait_trees("H", "1")


def test_from_config_reads_id_and_trait():
    cfg = {"id": "btv", "beast": {"discrete_trait": {"enabled": True, "trait": "host_group"}}}
    built = BuildLayout.from_config(cfg, "20260102")
    assert built == BuildLayout(Path("builds/btv/20260102"), "btv", "host_group")


def test_from_config_trait_is_none_unless_enabled():
    off = {"id": "btv", "beast": {"discrete_trait": {"enabled": False, "trait": "host_group"}}}
    assert BuildLayout.from_config(off, "20260102").trait is None
    assert BuildLayout.from_config({"id": "btv"}, "20260102").trait is None


def test_from_config_takes_a_builds_dir(tmp_path):
    built = BuildLayout.from_config({"id": "btv"}, "20260102", builds_dir=tmp_path)
    assert built.root == tmp_path / "btv" / "20260102"


def test_shipped_configs_resolve(configs):
    """Every shipped config yields a layout whose root is builds/<id>/<date>."""
    from wildphylo import config as cfgmod

    assert configs
    for path in configs:
        cfg = cfgmod.load(path)
        built = BuildLayout.from_config(cfg, "20260102")
        assert built.root == Path("builds") / cfg["id"] / "20260102"


def test_at_takes_the_root_as_given(tmp_path):
    cfg = {"id": "demo", "beast": {"discrete_trait": {"enabled": True, "trait": "host"}}}
    layout = BuildLayout.at(tmp_path / "anywhere", cfg)
    assert layout.root == tmp_path / "anywhere"
    assert layout.pathogen == "demo"
    assert layout.trait == "host"
    assert layout.metadata() == tmp_path / "anywhere" / "interim" / "metadata_clean.tsv"

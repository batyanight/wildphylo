"""
Tests for `wildphylo fetch records`, `wildphylo curate metadata` and `wildphylo curate review`.

Curation is tested end to end against synthetic GenBank fixtures
(tests/fixtures/make_genbank.py) rather than a live download. Every fixture
record exists to exercise one specific path, so the expected answer is known
independently of what the code does — and the tests do not depend on NCBI
being reachable or on what GenBank happens to contain this week.
"""

import logging
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

pd = pytest.importorskip("pandas")
pytest.importorskip("Bio")

from wildphylo import cli  # noqa: E402
from wildphylo.config import load  # noqa: E402
from wildphylo.stages import fetch as fetch_mod  # noqa: E402


@pytest.fixture(autouse=True)
def _run_from_repo_root(monkeypatch):
    """Config paths (host table, references table) are relative to the repo root."""
    monkeypatch.chdir(ROOT)


@pytest.fixture(scope="session")
def fixtures(tmp_path_factory):
    """Generate the synthetic GenBank files once per session."""
    out = tmp_path_factory.mktemp("gb")
    script = ROOT / "tests" / "fixtures" / "make_genbank.py"
    r = subprocess.run(
        [sys.executable, str(script)], capture_output=True, text=True, cwd=ROOT, check=False
    )
    assert r.returncode == 0, f"fixture generation failed: {r.stderr}"
    src = ROOT / "tests" / "fixtures"
    return {"cdv": src / "cdv_test.gb", "btv": src / "btv_test.gb", "tmp": out}


# `python -m wildphylo.cli` runs cli.py as __main__, a second copy of the module
# from the one the stage modules register their steps with, so no steps appear.
# Entering through `main` is what the installed `wildphylo` script does.
_ENTRY = "import sys; from wildphylo.cli import main; sys.exit(main(sys.argv[1:]))"


def run_curate(cfg: str, gb: Path, out: Path, *extra) -> subprocess.CompletedProcess:
    """Run `wildphylo curate metadata` in a subprocess, so stdout and stderr are separate."""
    return subprocess.run(
        [
            sys.executable,
            "-c",
            _ENTRY,
            "curate",
            "metadata",
            "--config",
            cfg,
            "--gb",
            str(gb),
            "--out-dir",
            str(out),
            "--log",
            str(out / "curate.log"),
            *extra,
        ],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=False,
    )


# =============================================================================
# fetch records — query construction
# =============================================================================


def test_query_uses_organism_exp():
    """
    [Organism:exp] expands to every strain and subtaxon beneath the node.
    Without it the query silently misses most records.
    """
    q = fetch_mod.build_query(load(ROOT / "tests/fixtures/configs/cdv.yaml"))
    assert "txid11232[Organism:exp]" in q


def test_query_uses_config_length_bounds():
    """Length bounds are pathogen-specific; BTV segments are short."""
    q = fetch_mod.build_query(load(ROOT / "config/pathogen/btv.yaml"))
    assert "300:4500[Sequence Length]" in q
    assert "txid40051" in q


def test_query_filter_with_leading_boolean_is_not_double_joined():
    cfg = load(ROOT / "tests/fixtures/configs/cdv.yaml")
    cfg["fetch"]["query_filter"] = "NOT patent[Title]"
    q = fetch_mod.build_query(cfg)
    assert "AND NOT" not in q
    assert q.endswith("NOT patent[Title]")


def test_query_filter_without_boolean_gets_an_and():
    cfg = load(ROOT / "tests/fixtures/configs/cdv.yaml")
    cfg["fetch"]["query_filter"] = "complete genome[Title]"
    assert "AND complete genome[Title]" in fetch_mod.build_query(cfg)


def test_empty_query_filter_leaves_no_dangling_operator():
    cfg = load(ROOT / "tests/fixtures/configs/cdv.yaml")
    cfg["fetch"]["query_filter"] = ""
    q = fetch_mod.build_query(cfg)
    assert not q.rstrip().endswith(("AND", "NOT", "OR"))


def test_fetch_requires_email(monkeypatch, tmp_path):
    """NCBI blocks anonymous bulk download; failing late wastes the run."""
    monkeypatch.setenv("NCBI_EMAIL", "")
    argv = ["fetch", "records", "--config", "tests/fixtures/configs/cdv.yaml", "--dry-run"]
    assert cli.main([*argv, "--log", str(tmp_path / "fetch.log")]) == 2


# =============================================================================
# fetch records — end to end with Entrez stubbed
# =============================================================================


@pytest.fixture
def stub_entrez(monkeypatch, fixtures):
    """Replace the two network calls with the synthetic CDV GenBank file."""
    text = fixtures["cdv"].read_text()
    total = sum(1 for line in text.splitlines() if line.startswith("VERSION"))
    calls = {"batches": 0, "reported_total": total}

    def esearch(_query):
        return "webenv", "1", calls["reported_total"]

    def batch(_webenv, _key, start, _size):
        calls["batches"] += 1
        return text if start == 0 else ""

    monkeypatch.setattr(fetch_mod, "esearch_history", esearch)
    monkeypatch.setattr(fetch_mod, "fetch_batch", batch)
    monkeypatch.setattr(fetch_mod.time, "sleep", lambda _s: None)
    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")
    calls["total"] = total
    return calls


def _fetch_argv(tmp_path, *extra):
    return [
        "fetch",
        "records",
        "--config",
        "tests/fixtures/configs/cdv.yaml",
        "--out",
        str(tmp_path / "raw" / "cdv.gb"),
        "--acc",
        str(tmp_path / "raw" / "cdv.acc"),
        "--log",
        str(tmp_path / "fetch.log"),
        *extra,
    ]


def test_fetch_records_writes_flat_file_and_accession_snapshot(stub_entrez, tmp_path, caplog):
    caplog.set_level(logging.INFO)  # pytest's handlers stop cli.main configuring the root level
    assert cli.main(_fetch_argv(tmp_path)) == 0
    accs = (tmp_path / "raw" / "cdv.acc").read_text().split()
    assert len(accs) == stub_entrez["total"]
    assert "AF100001.1" in accs
    assert (tmp_path / "raw" / "cdv.gb").stat().st_size > 0
    assert "matched" in (tmp_path / "fetch.log").read_text()


def test_fetch_records_dry_run_downloads_nothing(stub_entrez, tmp_path):
    assert cli.main(_fetch_argv(tmp_path, "--dry-run")) == 0
    assert stub_entrez["batches"] == 0
    assert not (tmp_path / "raw" / "cdv.gb").exists()


def test_fetch_records_leaves_an_existing_download_alone_without_force(stub_entrez, tmp_path):
    assert cli.main(_fetch_argv(tmp_path)) == 0
    assert cli.main(_fetch_argv(tmp_path)) == 0
    assert stub_entrez["batches"] == 1
    assert cli.main(_fetch_argv(tmp_path, "--force")) == 0
    assert stub_entrez["batches"] == 2


def test_a_short_download_is_an_error(stub_entrez, tmp_path):
    """NCBI reports more records than were received: every downstream count would be wrong."""
    stub_entrez["reported_total"] += 3
    assert cli.main(_fetch_argv(tmp_path)) == 1


def _honour_batch_window(monkeypatch, text):
    """Make the fetch_batch stub return only records start..start+size, like efetch."""
    records = [r + "//\n" for r in text.split("//\n") if r.strip()]
    requested = []

    def batch(_webenv, _key, start, size):
        requested.append((start, size))
        return "".join(records[start : start + size])

    monkeypatch.setattr(fetch_mod, "fetch_batch", batch)
    return records, requested


def test_max_records_caps_the_fetch_and_the_accession_snapshot(
    stub_entrez, fixtures, monkeypatch, tmp_path, caplog
):
    caplog.set_level(logging.INFO)
    records, requested = _honour_batch_window(monkeypatch, fixtures["cdv"].read_text())
    cap = 3
    assert cap < len(records)
    cfg = load(Path("tests/fixtures/configs/cdv.yaml"))
    cfg["fetch"]["max_records"] = cap
    monkeypatch.setattr(fetch_mod, "load", lambda _p: cfg)

    assert cli.main(_fetch_argv(tmp_path)) == 0

    accs = (tmp_path / "raw" / "cdv.acc").read_text().split()
    first_n = [
        line.split()[1]
        for r in records[:cap]
        for line in r.splitlines()
        if line.startswith("VERSION")
    ]
    assert accs == first_n
    assert sum(size for _start, size in requested) == cap
    log_text = (tmp_path / "fetch.log").read_text()
    assert f"capped at fetch.max_records = {cap}" in log_text
    assert f"{stub_entrez['total']} available" in log_text


def test_max_records_above_the_match_count_changes_nothing(stub_entrez, monkeypatch, tmp_path):
    cfg = load(Path("tests/fixtures/configs/cdv.yaml"))
    cfg["fetch"]["max_records"] = stub_entrez["total"] + 100
    monkeypatch.setattr(fetch_mod, "load", lambda _p: cfg)
    assert cli.main(_fetch_argv(tmp_path)) == 0
    accs = (tmp_path / "raw" / "cdv.acc").read_text().split()
    assert len(accs) == stub_entrez["total"]


# =============================================================================
# curate metadata — single locus (CDV)
# =============================================================================


@pytest.fixture(scope="session")
def cdv_out(fixtures, tmp_path_factory):
    out = tmp_path_factory.mktemp("cdv_cur")
    r = run_curate("tests/fixtures/configs/cdv.yaml", fixtures["cdv"], out)
    assert r.returncode == 0, r.stderr[-2000:]
    return out


def _tsv(out: Path, name: str):
    return pd.read_csv(out / name, sep="\t")


def test_curation_writes_every_expected_output(cdv_out):
    for f in (
        "metadata_all.tsv",
        "metadata_clean.tsv",
        "exclusions.tsv",
        "needs_review.tsv",
        "H.fasta",
    ):
        assert (cdv_out / f).is_file(), f"{f} was not written"


def test_malformed_date_does_not_crash_the_run(cdv_out):
    """
    The fixture contains 2013-13-45. The original parser raised ValueError and
    killed the whole run; the record must instead be excluded with a reason.
    """
    ex = _tsv(cdv_out, "exclusions.tsv")
    row = ex[ex["accession"] == "AF100007.1"]
    assert len(row) == 1
    assert row.iloc[0]["exclude_reason"] == "no_parseable_collection_date"


def test_implausible_year_is_excluded(cdv_out):
    ex = _tsv(cdv_out, "exclusions.tsv")
    assert "AF100008.1" in set(ex["accession"])
    rv = _tsv(cdv_out, "needs_review.tsv")
    reason = rv[rv["accession"] == "AF100008.1"].iloc[0]["date_reason"]
    assert "year_out_of_range" in reason


def test_vaccine_strain_excluded(cdv_out):
    ex = _tsv(cdv_out, "exclusions.tsv")
    row = ex[ex["accession"] == "AF100006.1"]
    assert row.iloc[0]["exclude_reason"] == "vaccine_or_vaccine_derived"


def test_record_without_the_locus_excluded(cdv_out):
    """A nucleocapsid-only record carries no H and cannot be used."""
    ex = _tsv(cdv_out, "exclusions.tsv")
    row = ex[ex["accession"] == "AF100010.1"]
    assert row.iloc[0]["exclude_reason"] == "no_analysed_locus_found"


def test_unrecognised_host_is_excluded_and_flagged(cdv_out):
    """A grey squirrel is not in the host table. It must not be guessed at."""
    ex = _tsv(cdv_out, "exclusions.tsv")
    assert ex[ex["accession"] == "AF100009.1"].iloc[0]["exclude_reason"] == "host_unresolved"
    assert "AF100009.1" in set(_tsv(cdv_out, "needs_review.tsv")["accession"])


def test_missing_date_is_excluded_but_not_flagged_for_review(cdv_out):
    """
    No date at all is missing data — there is nothing for a human to
    adjudicate. An unparseable date that WAS present is a curation gap and
    does need review. The two must not be conflated.
    """
    ex = _tsv(cdv_out, "exclusions.tsv")
    assert "AF100011.1" in set(ex["accession"])
    assert "AF100011.1" not in set(_tsv(cdv_out, "needs_review.tsv")["accession"])


def test_locus_alias_spellings_are_matched(cdv_out):
    """GenBank spells the same gene several ways; aliases come from config."""
    clean = _tsv(cdv_out, "metadata_clean.tsv")
    assert "AF100003.1" in set(clean["accession"])  # 'haemagglutinin'


def test_length_heuristic_fallback_and_review_flag(cdv_out):
    clean = _tsv(cdv_out, "metadata_clean.tsv")
    full = clean[clean["accession"] == "AF100004.1"].iloc[0]
    assert full["H_source"] == "length_heuristic_full"
    part = clean[clean["accession"] == "AF100005.1"].iloc[0]
    assert part["H_source"] == "length_heuristic_partial"
    # A partial hit is a guess, so it must be visible — but it needs no
    # decision, so it belongs in flagged.tsv rather than the review pile.
    assert "AF100005.1" in set(_tsv(cdv_out, "flagged.tsv")["accession"])


def test_year_only_date_uses_midpoint(cdv_out):
    clean = _tsv(cdv_out, "metadata_clean.tsv")
    row = clean[clean["accession"] == "AF100002.1"].iloc[0]
    assert row["decimal_year"] == pytest.approx(2013.5)
    assert row["date_precision"] == "year"


def test_pinniped_resolves_to_its_own_species_and_group(cdv_out):
    """
    The shadowing regression. `sea lion` used to be unreachable under `lion`,
    and Caspian seal collapsed to Phoca sp. Both must now be right.
    """
    clean = _tsv(cdv_out, "metadata_clean.tsv")
    row = clean[clean["accession"] == "AF100012.1"].iloc[0]
    assert row["host_group"] == "pinniped"
    assert row["host_canonical"] == "Pusa caspica"


def test_tip_labels_follow_the_downstream_contract(cdv_out):
    """ACCESSION|host_group|decimal_year is what 04 and 06 parse."""
    labels = [
        line[1:].strip()
        for line in (cdv_out / "H.fasta").read_text().splitlines()
        if line.startswith(">")
    ]
    assert labels
    for lab in labels:
        parts = lab.split("|")
        assert len(parts) == 3, lab
        float(parts[2])


def test_nothing_is_dropped_silently(cdv_out):
    """Every parsed record is either clean or carries an exclusion reason."""
    allr = _tsv(cdv_out, "metadata_all.tsv")
    clean = _tsv(cdv_out, "metadata_clean.tsv")
    ex = _tsv(cdv_out, "exclusions.tsv")
    assert len(clean) + len(ex) == len(allr)
    assert (ex["exclude_reason"].astype(str).str.strip() != "").all()


def test_keep_vaccines_flag_retains_them(fixtures, tmp_path):
    r = run_curate("tests/fixtures/configs/cdv.yaml", fixtures["cdv"], tmp_path, "--keep-vaccines")
    assert r.returncode == 0
    clean = pd.read_csv(tmp_path / "metadata_clean.tsv", sep="\t")
    assert "AF100006.1" in set(clean["accession"])


# =============================================================================
# curate metadata — multi-locus (BTV)
# =============================================================================


@pytest.fixture(scope="session")
def btv_out(fixtures, tmp_path_factory):
    out = tmp_path_factory.mktemp("btv_cur")
    r = run_curate("config/pathogen/btv.yaml", fixtures["btv"], out)
    assert r.returncode == 0, r.stderr[-2000:]
    return out


def test_one_fasta_per_analysed_locus(btv_out):
    for seg in ("seg2", "seg6", "seg10"):
        assert (btv_out / f"{seg}.fasta").is_file()
    # seg7 and seg3 are analyse: false
    assert not (btv_out / "seg7.fasta").exists()


def test_segments_share_taxon_labels(btv_out):
    """
    THE structural requirement for segmented pathogens. GenBank gives each
    segment its own accession, so labelling tips by accession leaves the
    shared-taxon set across segment trees EMPTY — and congruence, the entire
    reason for analysing segments separately, cannot be assessed at all.
    Labels must key on the isolate.
    """

    def labels(seg):
        return {
            line[1:].strip()
            for line in (btv_out / f"{seg}.fasta").read_text().splitlines()
            if line.startswith(">")
        }

    a, b, c = labels("seg2"), labels("seg6"), labels("seg10")
    assert a and a == b == c, (
        "segment trees share no taxa; congruence screening would be impossible"
    )


def test_segment_qualifier_assigns_the_right_locus(btv_out):
    meta = pd.read_csv(btv_out / "metadata_clean.tsv", sep="\t")
    for _, r in meta.iterrows():
        seg = str(r["segment_raw"]).strip()
        if seg in ("2", "6", "10"):
            assert bool(r[f"has_seg{seg}"]), f"{r['accession']} lost segment {seg}"


def test_record_with_no_segment_and_no_annotation_is_excluded(btv_out):
    ex = pd.read_csv(btv_out / "exclusions.tsv", sep="\t")
    assert "MN299999.1" in set(ex["accession"])


def test_vector_host_is_labelled_not_dropped(btv_out):
    """
    vector_policy is 'annotate': Culicoides sequences inform the tree but are
    excluded from the trait analysis by the dta_states map, not here.
    """
    meta = pd.read_csv(btv_out / "metadata_clean.tsv", sep="\t")
    assert "vector" in set(meta["host_group"])


def test_wild_and_domestic_ruminants_separated(btv_out):
    """bighorn sheep must not collapse into domestic sheep."""
    meta = pd.read_csv(btv_out / "metadata_clean.tsv", sep="\t")
    groups = set(meta["host_group"])
    assert {"domestic_sheep", "wild_bovid", "wild_cervid"} <= groups


# =============================================================================
# Reference anchors
# =============================================================================
#
# Reference sequences NAME clades; they are not dated tips. All eight usable
# CDV nucleotide references lack /collection_date and /host, so requiring them
# to pass the analysis filters removed every one of them — and the ten clades
# recovered from the first real run were therefore unnameable.

REFS = "tests/fixtures/refs_test.tsv"


@pytest.fixture(scope="session")
def anchor_cfg(tmp_path_factory):
    import yaml

    cfg = yaml.safe_load((ROOT / "tests/fixtures/configs/cdv.yaml").read_text())
    cfg["references"]["table"] = REFS
    cfg["references"]["include_as_anchors"] = True
    p = tmp_path_factory.mktemp("cfg") / "anchor.yaml"
    p.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return p


@pytest.fixture(scope="session")
def anchored(fixtures, anchor_cfg, tmp_path_factory):
    out = tmp_path_factory.mktemp("anchored")
    r = run_curate(str(anchor_cfg), fixtures["cdv"], out)
    assert r.returncode == 0, r.stderr[-2000:]
    return {"dir": out, "stderr": r.stderr, "stdout": r.stdout}


def _labels(path: Path) -> set[str]:
    return {line[1:].strip() for line in path.read_text().splitlines() if line.startswith(">")}


def test_undated_references_survive_as_anchors(anchored):
    """
    The regression. AF164967 and Z47762 have no /host and no /collection_date,
    exactly like the real references, and were being dropped as
    no_parseable_collection_date.
    """
    labels = _labels(anchored["dir"] / "H_anchors.fasta")
    accs = {line.split("|")[0] for line in labels}
    assert {"AF164967", "Z47762"} <= accs


def test_anchor_labels_carry_a_non_numeric_date(anchored):
    """
    Every downstream parser reads an unparseable date as "no date" and drops
    the tip from the regression. That is how an anchor stays in the tree while
    staying out of the clock.
    """
    for lab in _labels(anchored["dir"] / "H_anchors.fasta"):
        parts = lab.split("|")
        assert len(parts) == 3
        with pytest.raises(ValueError, match="could not convert"):
            float(parts[2])


def test_anchor_labels_carry_the_lineage(anchored):
    for lab in _labels(anchored["dir"] / "H_anchors.fasta"):
        assert lab.split("|")[1].startswith("ref_")
    assert any("America-2" in line for line in _labels(anchored["dir"] / "H_anchors.fasta"))


def test_anchors_are_not_in_the_analysis_fasta(anchored):
    """An anchor must not become a dated tip in the clock analysis."""
    accs = {line.split("|")[0] for line in _labels(anchored["dir"] / "H.fasta")}
    assert not ({"AF164967", "Z47762"} & accs)


def test_protein_accessions_are_skipped_with_a_warning(anchored):
    out = anchored["stdout"] + anchored["stderr"]
    assert "protein accessions" in out


def test_absent_references_are_reported(anchored):
    """
    A reference in the table but not in the data means that lineage cannot be
    named. Silence there would be worse than the missing reference.
    """
    out = anchored["stdout"] + anchored["stderr"]
    assert "NOTFETCHED" in out


def test_anchor_file_is_written_even_when_disabled(fixtures, tmp_path):
    """The workflow declares it as a fixed input, so it must always exist."""
    r = run_curate("config/pathogen/btv.yaml", fixtures["btv"], tmp_path)
    assert r.returncode == 0
    for seg in ("seg2", "seg6", "seg10"):
        assert (tmp_path / f"{seg}_anchors.fasta").is_file()


def test_an_anchor_without_the_locus_is_still_dropped(anchored):
    """An anchor with no sequence anchors nothing; that exemption has a limit."""
    meta = pd.read_csv(anchored["dir"] / "metadata_clean.tsv", sep="\t")
    refs = meta[meta["is_reference"]]
    assert len(refs) == 2
    assert refs["H_length"].min() > 0


# =============================================================================
# Locus length validation
# =============================================================================
#
# Annotation-based extraction used to be trusted unconditionally: a feature
# labelled with the locus name was accepted at any length. On real CDV data
# that admitted sequences of 1946 and 1947 nt against an 1824 nt H CDS, and
# 29 of 162 published clade-3 sequences have lengths that are not multiples
# of 3 and therefore cannot be translated in frame.


def test_overlong_annotated_feature_is_flagged(cdv_out):
    """
    1946 nt against an 1824 nt CDS — a `gene` feature spanning UTR, or a
    mis-annotation. It is 6.7% over, so the 15% heuristic tolerance let it
    through; annotation needs a tighter one.
    """
    fl = _tsv(cdv_out, "flagged.tsv")
    row = fl[fl["accession"] == "AF100015.1"]
    assert len(row) == 1, "an over-length annotated feature was not flagged"
    assert row.iloc[0]["H_source"].endswith("_overlength")


def test_out_of_frame_length_is_flagged(cdv_out):
    """1823 nt is not a multiple of 3, so it cannot be translated in frame."""
    assert "AF100016.1" in set(_tsv(cdv_out, "flagged.tsv")["accession"])
    allr = _tsv(cdv_out, "metadata_all.tsv")
    row = allr[allr["accession"] == "AF100016.1"].iloc[0]
    assert row["H_in_frame"] == False  # noqa: E712 — pandas bool column


def test_partial_annotations_are_not_flagged(cdv_out):
    """
    A feature SHORTER than the CDS is a partial annotation — a real sequence
    covering part of the gene, which is most of what GenBank holds. Flagging
    those put 1209 of 1390 CDV records into needs_review, a pile nobody can
    adjudicate, and the commonest of them was the 1338 nt amplicon that
    dominates the published clade 3.
    """
    allr = _tsv(cdv_out, "metadata_all.tsv")
    partial = allr[allr["has_H"] & (allr["H_length"] < 1824) & (allr["H_length"] >= 400)]
    assert len(partial), "fixture has no partial-H record to test with"
    flagged = partial["H_source"].astype(str).str.endswith("_overlength")
    assert not flagged.any(), "a partial annotation was flagged as suspicious"


def test_correct_length_is_not_flagged(cdv_out):
    """The check must not fire on a clean full-length CDS."""
    allr = _tsv(cdv_out, "metadata_all.tsv")
    row = allr[allr["accession"] == "AF100001.1"].iloc[0]
    assert row["H_length"] == 1824
    assert row["H_in_frame"] == True  # noqa: E712
    assert row["H_source"] == "annotation"


def test_in_frame_column_exists_for_every_locus(btv_out):
    meta = pd.read_csv(btv_out / "metadata_all.tsv", sep="\t")
    for seg in ("seg2", "seg6", "seg10"):
        assert f"{seg}_in_frame" in meta.columns


def test_cds_span_must_be_a_multiple_of_three():
    """
    `cds: [21, 1836]` spans 1816 nt = 605.33 codons — arithmetically impossible
    for a coding sequence, and it shipped in the CDV config feeding the
    annotation scripts.
    """
    import yaml

    from wildphylo.config import validate

    cfg = yaml.safe_load((ROOT / "tests/fixtures/configs/cdv.yaml").read_text())
    cfg["_path"] = "x"
    errors, _ = validate(cfg)
    assert not any("multiple of 3" in e for e in errors), (
        "the shipped CDV cds span is not a multiple of 3"
    )

    cfg["locus"]["coordinate_reference"]["cds"] = [21, 1836]
    errors, _ = validate(cfg)
    assert any("multiple of 3" in e for e in errors)


def test_cds_span_disagreeing_with_expected_length_warns():
    import yaml

    from wildphylo.config import validate

    cfg = yaml.safe_load((ROOT / "tests/fixtures/configs/cdv.yaml").read_text())
    cfg["_path"] = "x"
    cfg["locus"]["coordinate_reference"]["cds"] = [1, 1200]
    _e, warnings = validate(cfg)
    assert any("expected_length" in w and "wrong" in w for w in warnings)


# =============================================================================
# needs_review vs flagged
# =============================================================================
#
# The two were conflated and the pile reached 1209 records on the real CDV
# dataset — a volume nobody adjudicates, which makes the gate worthless.


def test_review_holds_only_records_needing_a_decision(cdv_out):
    rv = _tsv(cdv_out, "needs_review.tsv")
    assert len(rv), "fixture produced no actionable records"
    assert set(rv["review_reason"]) <= {"host_unrecognised", "date_unparseable"}


def test_sequence_property_flags_go_to_flagged_not_review(cdv_out):
    """
    'This sequence is partial' needs no decision. Burying the actionable
    records among hundreds of these is how a gate stops being read.
    """
    rv = _tsv(cdv_out, "needs_review.tsv")
    fl = _tsv(cdv_out, "flagged.tsv")
    assert (cdv_out / "flagged.tsv").is_file()
    # the over-length record is informational, not actionable
    assert "AF100015.1" in set(fl["accession"])
    assert "AF100015.1" not in set(rv["accession"])


def test_review_and_flagged_do_not_overlap(cdv_out):
    rv = set(_tsv(cdv_out, "needs_review.tsv")["accession"])
    fl = set(_tsv(cdv_out, "flagged.tsv")["accession"])
    assert not (rv & fl)


def test_unrecognised_host_is_still_actionable(cdv_out):
    """The grey squirrel: a human adds a pattern or confirms the drop."""
    rv = _tsv(cdv_out, "needs_review.tsv")
    row = rv[rv["accession"] == "AF100009.1"]
    assert len(row) == 1
    assert row.iloc[0]["review_reason"] == "host_unrecognised"


# =============================================================================
# curate metadata — in-process through cli.main
# =============================================================================


def test_curate_metadata_runs_in_process_through_the_cli(fixtures, tmp_path, caplog):
    caplog.set_level(logging.INFO)  # pytest's handlers stop cli.main configuring the root level
    out = tmp_path / "interim"
    argv = [
        "curate",
        "metadata",
        "--config",
        "tests/fixtures/configs/cdv.yaml",
        "--gb",
        str(fixtures["cdv"]),
        "--out-dir",
        str(out),
        "--log",
        str(tmp_path / "curate.log"),
    ]
    assert cli.main(argv) == 0
    assert (out / "H.fasta").is_file()
    assert (out / "needs_review.tsv").is_file()
    assert "NEEDS HUMAN REVIEW" in (tmp_path / "curate.log").read_text()


def test_curate_metadata_rejects_a_missing_genbank_file(tmp_path):
    argv = [
        "curate",
        "metadata",
        "--config",
        "tests/fixtures/configs/cdv.yaml",
        "--gb",
        str(tmp_path / "absent.gb"),
        "--out-dir",
        str(tmp_path / "interim"),
        "--log",
        str(tmp_path / "curate.log"),
    ]
    assert cli.main(argv) == 2


# =============================================================================
# curate review — the review gate
# =============================================================================


def _review_argv(needs_review: Path, marker: Path, *extra: str) -> list[str]:
    return [
        "curate",
        "review",
        "--needs-review",
        str(needs_review),
        "--marker",
        str(marker),
        *extra,
    ]


def test_review_gate_passes_and_writes_the_marker_when_nothing_needs_review(tmp_path):
    needs_review = tmp_path / "needs_review.tsv"
    needs_review.write_text("accession\treview_reason\n")
    marker = tmp_path / "interim" / ".review_ok"
    assert cli.main(_review_argv(needs_review, marker)) == 0
    assert marker.is_file()


def test_review_gate_refuses_while_records_need_review(tmp_path, caplog):
    needs_review = tmp_path / "needs_review.tsv"
    needs_review.write_text("accession\treview_reason\nAF100009.1\thost_unrecognised\n")
    marker = tmp_path / ".review_ok"
    assert cli.main(_review_argv(needs_review, marker)) == 1
    assert not marker.exists()
    assert "1 records need manual review" in caplog.text
    assert "--accept-unreviewed" in caplog.text


def test_review_gate_refusal_removes_a_stale_marker(tmp_path):
    needs_review = tmp_path / "needs_review.tsv"
    needs_review.write_text("accession\treview_reason\nAF100009.1\thost_unrecognised\n")
    marker = tmp_path / ".review_ok"
    marker.touch()
    assert cli.main(_review_argv(needs_review, marker)) == 1
    assert not marker.exists()


def test_review_gate_accept_unreviewed_proceeds_with_a_warning(tmp_path, caplog):
    needs_review = tmp_path / "needs_review.tsv"
    needs_review.write_text("accession\treview_reason\nAF100009.1\thost_unrecognised\n")
    marker = tmp_path / ".review_ok"
    assert cli.main(_review_argv(needs_review, marker, "--accept-unreviewed")) == 0
    assert marker.is_file()
    assert "proceed unreviewed" in caplog.text


def test_review_gate_refuses_on_real_curation_output(cdv_out, tmp_path):
    """The CDV fixture produces actionable records (the grey squirrel), so the gate must stop."""
    marker = tmp_path / ".review_ok"
    assert cli.main(_review_argv(cdv_out / "needs_review.tsv", marker)) == 1
    assert not marker.exists()

"""
Tests for `wildphylo align qc` and `wildphylo beast subsample`.

Both previously wrote to a hardcoded data/processed/. A rule that
writes outside its declared outputs cannot be tracked, re-run or cleaned by
Snakemake, and the stray file silently goes stale — so "the file landed where
it was asked to, and nowhere else" is the property worth testing.
"""

import contextlib
import io
import logging
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from wildphylo import cli

ROOT = Path(__file__).resolve().parents[1]

pd = pytest.importorskip("pandas")
pytest.importorskip("Bio")


@pytest.fixture(scope="module")
def curated(tmp_path_factory):
    """Curate the CDV fixture, then pad to a fake alignment."""
    out = tmp_path_factory.mktemp("qc")
    subprocess.run(
        [sys.executable, str(ROOT / "tests/fixtures/make_genbank.py")],
        capture_output=True,
        cwd=ROOT,
        check=True,
    )
    r = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from wildphylo.cli import main; sys.exit(main(sys.argv[1:]))",
            "curate",
            "metadata",
            "--config",
            "tests/fixtures/configs/cdv.yaml",
            "--gb",
            "tests/fixtures/cdv_test.gb",
            "--out-dir",
            str(out),
        ],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=False,
    )
    assert r.returncode == 0, r.stderr[-1500:]

    from Bio import SeqIO

    recs = list(SeqIO.parse(str(out / "H.fasta"), "fasta"))
    width = max(len(x.seq) for x in recs)
    aln = out / "H_aligned.fasta"
    aln.write_text("".join(f">{x.id}\n{str(x.seq).ljust(width, '-')}\n" for x in recs))
    return {"dir": out, "aln": aln, "meta": out / "metadata_clean.tsv"}


# --- align qc ---------------------------------------------------------------


def run_qc(*args):
    """Run `wildphylo align qc` in-process; return returncode, stdout and stderr."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = cli.main(["align", "qc", *map(str, args)])
    return SimpleNamespace(returncode=code, stdout=out.getvalue(), stderr=err.getvalue())


def test_qc_writes_to_the_requested_paths(curated, tmp_path):
    tsv, png = tmp_path / "deep" / "qc.tsv", tmp_path / "deep" / "cov.png"
    r = run_qc(
        "--config",
        "tests/fixtures/configs/cdv.yaml",
        "--aln",
        curated["aln"],
        "--metadata",
        curated["meta"],
        "--out-tsv",
        tsv,
        "--out-png",
        png,
    )
    assert r.returncode == 0, r.stderr[-1500:]
    assert tsv.is_file() and png.is_file()


def test_qc_creates_missing_parent_directories(curated, tmp_path):
    """Snakemake does not always create a rule's output directory first."""
    tsv = tmp_path / "a" / "b" / "c" / "qc.tsv"
    r = run_qc(
        "--aln",
        curated["aln"],
        "--out-tsv",
        tsv,
        "--out-png",
        tmp_path / "a" / "b" / "c" / "p.png",
    )
    assert r.returncode == 0
    assert tsv.is_file()


def test_qc_reads_thresholds_from_config(curated, tmp_path):
    r = run_qc(
        "--config",
        "tests/fixtures/configs/cdv.yaml",
        "--aln",
        curated["aln"],
        "--out-tsv",
        tmp_path / "q.tsv",
        "--out-png",
        tmp_path / "q.png",
    )
    assert r.returncode == 0
    # cdv.yaml sets min_identity_to_consensus 0.80, max_gap_fraction 0.60
    assert "0.6" in r.stdout or "0.8" in r.stdout


def _threshold_fixture(tmp_path):
    """A tiny alignment with sequences between the default and a stricter threshold.

    ``g`` is 40% gaps: it passes the default (flag above 50%) but fails a config
    with ``alignment.max_gap_fraction`` 0.30. ``m`` is 90% identical to the
    consensus: it passes the default (flag below 80%) but fails a config with
    ``alignment.min_identity_to_consensus`` 0.95. ``n`` is 5% ambiguous: it
    passes the default (flag above 5%) but fails a config with
    ``alignment.max_ambiguity_fraction`` 0.02.
    """
    import yaml

    core = "ATGGCTGAAGCTGAAGCTGA"
    gappy = "ATGGCTGAAGCT" + "-" * 8
    mismatched = "ATGGCTGAAGCTGAAGCCCC"
    ambiguous = "ATGGCTGAAGNTGAAGCTGA"
    rows = {"a": core, "b": core, "c": core, "g": gappy, "m": mismatched, "n": ambiguous}
    aln = tmp_path / "tiny.fasta"
    aln.write_text("".join(f">{k}|dog|2000\n{v}\n" for k, v in rows.items()))

    def config(**alignment):
        cfg = yaml.safe_load(Path("tests/fixtures/configs/cdv.yaml").read_text())
        cfg["alignment"].update(alignment)
        path = tmp_path / "strict.yaml"
        path.write_text(yaml.safe_dump(cfg))
        return path

    return aln, config


def _flags(tmp_path, aln, *extra):
    """Run QC on ``aln`` and return ``{sequence name: flag string}``."""
    tsv = tmp_path / "qc.tsv"
    r = run_qc("--aln", aln, "--out-tsv", tsv, "--out-png", tmp_path / "qc.png", *extra)
    assert r.returncode == 0, r.stderr[-1500:]
    table = pd.read_csv(tsv, sep="\t", keep_default_na=False)
    return dict(zip(table["name"], table["flag"], strict=True))


def test_qc_flags_use_the_config_gap_threshold(tmp_path):
    """The config threshold must decide the flags, not only the lines on the plot."""
    aln, config = _threshold_fixture(tmp_path)
    assert "mostly_gaps" not in _flags(tmp_path, aln)["g|dog|2000"]  # default 0.5 passes 40%
    flags = _flags(tmp_path, aln, "--config", config(max_gap_fraction=0.30))
    assert "mostly_gaps" in flags["g|dog|2000"]
    assert flags["a|dog|2000"] == ""


def test_qc_flags_use_the_config_identity_threshold(tmp_path):
    aln, config = _threshold_fixture(tmp_path)
    assert "low_identity" not in _flags(tmp_path, aln)["m|dog|2000"]
    flags = _flags(tmp_path, aln, "--config", config(min_identity_to_consensus=0.95))
    assert "low_identity" in flags["m|dog|2000"]


def test_qc_flags_use_the_config_ambiguity_threshold(tmp_path):
    aln, config = _threshold_fixture(tmp_path)
    assert "ambiguous" not in _flags(tmp_path, aln)["n|dog|2000"]
    flags = _flags(tmp_path, aln, "--config", config(max_ambiguity_fraction=0.02))
    assert "ambiguous>2%" in flags["n|dog|2000"]
    assert "ambiguous" not in flags["a|dog|2000"]


def test_qc_command_line_ambiguity_threshold_beats_the_config(tmp_path):
    aln, config = _threshold_fixture(tmp_path)
    strict = config(max_ambiguity_fraction=0.02)
    flags = _flags(tmp_path, aln, "--config", strict, "--ambiguity-threshold", "0.05")
    assert "ambiguous" not in flags["n|dog|2000"]


def test_qc_command_line_threshold_beats_the_config(tmp_path):
    """An explicit flag wins even when it equals the built-in default."""
    aln, config = _threshold_fixture(tmp_path)
    strict = config(max_gap_fraction=0.30)
    for explicit in ("0.5", "0.45"):
        flags = _flags(tmp_path, aln, "--config", strict, "--gap-threshold", explicit)
        assert "mostly_gaps" not in flags["g|dog|2000"]


def test_qc_cross_checks_labels_against_metadata(curated, tmp_path):
    """
    The label/metadata cross-check is what catches a tip-label bug in 02. On
    correctly curated data it must report zero disagreements.
    """
    r = run_qc(
        "--aln",
        curated["aln"],
        "--metadata",
        curated["meta"],
        "--out-tsv",
        tmp_path / "q.tsv",
        "--out-png",
        tmp_path / "q.png",
    )
    assert r.returncode == 0
    for line in r.stdout.splitlines():
        if "disagrees with metadata" in line or "not in metadata" in line:
            assert line.rstrip().endswith("0"), line


def run_cli(*argv):
    """Run ``wildphylo <argv>`` in-process; return exit code, stdout and stderr."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = cli.main([str(a) for a in argv])
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
            if isinstance(exc.code, str):
                err.write(exc.code)
    return SimpleNamespace(returncode=code, stdout=out.getvalue(), stderr=err.getvalue())


# --- 06 ---------------------------------------------------------------------


def test_subsample_writes_only_where_asked(curated, tmp_path):
    """
    Regression. subsample wrote a *_dates.tsv to data/processed/ regardless of
    --out-meta, and created that directory on every run even when unused.
    """
    existed_before = (ROOT / "data").exists()
    aln, meta = tmp_path / "s" / "sub.fasta", tmp_path / "s" / "sub_metadata.tsv"
    r = run_cli(
        "beast",
        "subsample",
        "--config",
        str(ROOT / "tests/fixtures/configs/cdv.yaml"),
        "--mode",
        "global",
        "--aln",
        curated["aln"],
        "--metadata",
        curated["meta"],
        "--target",
        "4",
        "--out-aln",
        aln,
        "--out-meta",
        meta,
    )
    assert r.returncode == 0, r.stderr[-1500:]
    assert aln.is_file() and meta.is_file()
    assert (tmp_path / "s" / "sub_dates.tsv").is_file(), (
        "the dates file escaped the requested directory"
    )
    # Check what THIS run created, not what merely exists. A stale data/ left
    # by the un-ported scripts is not a failure of the port, and asserting on
    # absolute absence made this fail for anyone who had run the old code.
    assert existed_before or not (ROOT / "data").exists(), (
        "subsample created data/ even though explicit outputs were given"
    )


def test_wild_groups_come_from_config_not_the_module_default(curated, tmp_path):
    """
    THE reason this needed porting. The module default is CDV's carnivore
    range. Applied to a ruminant pathogen it classifies every wild host as not
    wild, silently breaking the wild/domestic balance in subsampling.
    """
    common = [
        "--mode",
        "global",
        "--aln",
        str(curated["aln"]),
        "--metadata",
        str(curated["meta"]),
        "--target",
        "4",
    ]
    cdv = run_cli(
        "beast",
        "subsample",
        "--config",
        str(ROOT / "tests/fixtures/configs/cdv.yaml"),
        *common,
        "--out-aln",
        tmp_path / "c.fasta",
        "--out-meta",
        tmp_path / "c_metadata.tsv",
    )
    btv = run_cli(
        "beast",
        "subsample",
        "--config",
        str(ROOT / "config/pathogen/btv.yaml"),
        *common,
        "--out-aln",
        tmp_path / "b.fasta",
        "--out-meta",
        tmp_path / "b_metadata.tsv",
    )

    def wild_line(out):
        return next(line for line in out.splitlines() if line.startswith("wild host groups:"))

    c, b = wild_line(cdv.stdout), wild_line(btv.stdout)
    assert c != b, "both pathogens got the same wild-group set"
    assert "procyonid" in c and "procyonid" not in b
    assert "wild_cervid" in b and "wild_cervid" not in c


def test_subsample_is_deterministic(curated, tmp_path):
    """A seeded subsample must be reproducible, or no build can be compared."""

    def labels(tag):
        a = tmp_path / f"{tag}.fasta"
        run_cli(
            "beast",
            "subsample",
            "--config",
            str(ROOT / "tests/fixtures/configs/cdv.yaml"),
            "--mode",
            "global",
            "--aln",
            curated["aln"],
            "--metadata",
            curated["meta"],
            "--target",
            "4",
            "--seed",
            "999",
            "--out-aln",
            a,
            "--out-meta",
            tmp_path / f"{tag}_metadata.tsv",
        )
        return [line for line in a.read_text().splitlines() if line.startswith(">")]

    assert labels("r1") == labels("r2")


def test_subsample_output_labels_keep_the_contract(curated, tmp_path):
    aln = tmp_path / "sub.fasta"
    run_cli(
        "beast",
        "subsample",
        "--config",
        str(ROOT / "tests/fixtures/configs/cdv.yaml"),
        "--mode",
        "global",
        "--aln",
        curated["aln"],
        "--metadata",
        curated["meta"],
        "--target",
        "4",
        "--out-aln",
        aln,
        "--out-meta",
        tmp_path / "sub_metadata.tsv",
    )
    for line in aln.read_text().splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            assert len(parts) == 3
            float(parts[2])


# --- reference anchors in QC ------------------------------------------------


@pytest.fixture(scope="module")
def curated_with_anchors(tmp_path_factory):
    """Curate with anchors on, then pad to a fake alignment including them."""
    import subprocess

    import yaml

    out = tmp_path_factory.mktemp("anc")
    subprocess.run(
        [sys.executable, str(ROOT / "tests/fixtures/make_genbank.py")],
        capture_output=True,
        cwd=ROOT,
        check=True,
    )
    cfg = yaml.safe_load((ROOT / "tests/fixtures/configs/cdv.yaml").read_text())
    cfg["references"]["table"] = "tests/fixtures/refs_test.tsv"
    cfg["references"]["include_as_anchors"] = True
    cfgp = out / "cfg.yaml"
    cfgp.write_text(yaml.safe_dump(cfg, sort_keys=False))
    r = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from wildphylo.cli import main; sys.exit(main(sys.argv[1:]))",
            "curate",
            "metadata",
            "--config",
            str(cfgp),
            "--gb",
            "tests/fixtures/cdv_test.gb",
            "--out-dir",
            str(out),
        ],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=False,
    )
    assert r.returncode == 0, r.stderr[-1500:]

    from Bio import SeqIO

    recs = list(SeqIO.parse(str(out / "H.fasta"), "fasta")) + list(
        SeqIO.parse(str(out / "H_anchors.fasta"), "fasta")
    )
    w = max(len(x.seq) for x in recs)
    aln = out / "aln.fasta"
    aln.write_text("".join(f">{x.id}\n{str(x.seq).ljust(w, '-')}\n" for x in recs))
    return {"aln": aln, "meta": out / "metadata_clean.tsv"}


def test_qc_survives_empty_metadata_cells(curated_with_anchors, tmp_path):
    """
    REGRESSION. `dtype=str` does not make missing cells strings — pandas leaves
    them as float NaN, so `.strip()` raised AttributeError. Anchors were the
    first records with an empty decimal_year, and the whole run died at QC
    after the alignment had already been built.
    """
    r = run_qc(
        "--aln",
        curated_with_anchors["aln"],
        "--metadata",
        curated_with_anchors["meta"],
        "--out-tsv",
        tmp_path / "q.tsv",
        "--out-png",
        tmp_path / "q.png",
    )
    assert r.returncode == 0, r.stderr[-1500:]


def test_anchors_are_not_reported_as_disagreements(curated_with_anchors, tmp_path):
    """
    An anchor carries its LINEAGE in the host field and a non-numeric date, by
    design. Flagging that as a mismatch would teach the reader to ignore the
    one section that catches real tip-label bugs.
    """
    r = run_qc(
        "--aln",
        curated_with_anchors["aln"],
        "--metadata",
        curated_with_anchors["meta"],
        "--out-tsv",
        tmp_path / "q.tsv",
        "--out-png",
        tmp_path / "q.png",
    )
    assert r.returncode == 0
    out = r.stdout
    assert "reference anchors (not cross-checked): 2" in out
    for line in out.splitlines():
        if "disagrees with metadata" in line:
            assert line.rstrip().endswith("0"), line


def test_qc_reports_a_missing_alignment_and_fails(tmp_path):
    r = run_qc("--aln", tmp_path / "nope.fasta", "--out-tsv", tmp_path / "q.tsv")
    assert r.returncode == 1
    assert "Alignment not found" in r.stderr


# --- subsample --build-dir ---------------------------------------------------


def _subsample_build(curated, tmp_path):
    import shutil

    import yaml

    from wildphylo.layout import BuildLayout

    cfg_path = ROOT / "tests/fixtures/configs/cdv.yaml"
    layout = BuildLayout.at(
        tmp_path / "builds" / "cdv" / "20260101", yaml.safe_load(cfg_path.read_text())
    )
    layout.alignment("H").parent.mkdir(parents=True)
    shutil.copy(curated["aln"], layout.alignment("H"))
    layout.metadata().parent.mkdir(parents=True)
    shutil.copy(curated["meta"], layout.metadata())
    return cfg_path, layout


def test_subsample_build_dir_reads_and_writes_the_layout(curated, tmp_path):
    cfg, layout = _subsample_build(curated, tmp_path)
    r = run_cli(
        "beast",
        *("subsample", "--config", cfg, "--build-dir", layout.root, "--locus", "H"),
        *("--mode", "global", "--target", "4"),
    )
    assert r.returncode == 0, r.stderr[-1500:]
    assert layout.subsample("H").is_file()
    assert layout.subsample_metadata("H").is_file()


def test_subsample_explicit_outputs_override_the_build_dir(curated, tmp_path):
    cfg, layout = _subsample_build(curated, tmp_path)
    aln = tmp_path / "elsewhere" / "sub.fasta"
    r = run_cli(
        "beast",
        *("subsample", "--config", cfg, "--build-dir", layout.root, "--locus", "H"),
        *("--mode", "global", "--target", "4", "--out-aln", aln),
    )
    assert r.returncode == 0, r.stderr[-1500:]
    assert aln.is_file()
    assert not layout.subsample("H").exists()
    assert layout.subsample_metadata("H").is_file()


def test_subsample_build_dir_needs_locus(curated, tmp_path, caplog):
    caplog.set_level(logging.ERROR, logger="wildphylo")
    cfg, layout = _subsample_build(curated, tmp_path)
    r = run_cli(
        "beast",
        *("subsample", "--config", cfg, "--build-dir", layout.root, "--mode", "global"),
    )
    assert r.returncode == 2
    assert "--locus" in caplog.text


def test_subsample_clade_mode_without_a_clade_is_a_usage_error(curated, caplog):
    caplog.set_level(logging.ERROR, logger="wildphylo")
    r = run_cli("beast", "subsample", "--aln", curated["aln"], "--mode", "clade")
    assert r.returncode == 2
    assert "--mode clade requires --clade" in caplog.text

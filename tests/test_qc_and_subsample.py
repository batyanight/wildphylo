"""
Tests for the ported 04_alignment_qc.py and 06_subsample.py.

Both scripts previously wrote to a hardcoded data/processed/. A rule that
writes outside its declared outputs cannot be tracked, re-run or cleaned by
Snakemake, and the stray file silently goes stale — so "the file landed where
it was asked to, and nowhere else" is the property worth testing.
"""

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

pd = pytest.importorskip("pandas")
pytest.importorskip("Bio")


@pytest.fixture(scope="module")
def curated(tmp_path_factory):
    """Curate the CDV fixture, then pad to a fake alignment."""
    out = tmp_path_factory.mktemp("qc")
    subprocess.run([sys.executable, str(ROOT / "tests/fixtures/make_genbank.py")],
                   capture_output=True, cwd=ROOT, check=True)
    r = subprocess.run(
        [sys.executable, "scripts/02_curate_metadata.py",
         "--config", "config/pathogen/cdv.yaml",
         "--gb", "tests/fixtures/cdv_test.gb", "--out-dir", str(out)],
        capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-1500:]

    from Bio import SeqIO
    recs = list(SeqIO.parse(str(out / "H.fasta"), "fasta"))
    width = max(len(x.seq) for x in recs)
    aln = out / "H_aligned.fasta"
    aln.write_text("".join(f">{x.id}\n{str(x.seq).ljust(width, '-')}\n" for x in recs))
    return {"dir": out, "aln": aln, "meta": out / "metadata_clean.tsv"}


def run(script, *args, cwd=None):
    return subprocess.run([sys.executable, f"scripts/{script}", *map(str, args)],
                          capture_output=True, text=True, cwd=cwd or ROOT)


# --- 04 ---------------------------------------------------------------------

def test_qc_writes_to_the_requested_paths(curated, tmp_path):
    tsv, png = tmp_path / "deep" / "qc.tsv", tmp_path / "deep" / "cov.png"
    r = run("04_alignment_qc.py", "--config", "config/pathogen/cdv.yaml",
            "--aln", curated["aln"], "--metadata", curated["meta"],
            "--out-tsv", tsv, "--out-png", png)
    assert r.returncode == 0, r.stderr[-1500:]
    assert tsv.is_file() and png.is_file()


def test_qc_creates_missing_parent_directories(curated, tmp_path):
    """Snakemake does not always create a rule's output directory first."""
    tsv = tmp_path / "a" / "b" / "c" / "qc.tsv"
    r = run("04_alignment_qc.py", "--aln", curated["aln"], "--out-tsv", tsv,
            "--out-png", tmp_path / "a" / "b" / "c" / "p.png")
    assert r.returncode == 0
    assert tsv.is_file()


def test_qc_reads_thresholds_from_config(curated, tmp_path):
    r = run("04_alignment_qc.py", "--config", "config/pathogen/cdv.yaml",
            "--aln", curated["aln"], "--out-tsv", tmp_path / "q.tsv",
            "--out-png", tmp_path / "q.png")
    assert r.returncode == 0
    # cdv.yaml sets min_identity_to_consensus 0.80, max_gap_fraction 0.60
    assert "0.6" in r.stdout or "0.8" in r.stdout


def test_qc_cross_checks_labels_against_metadata(curated, tmp_path):
    """
    The label/metadata cross-check is what catches a tip-label bug in 02. On
    correctly curated data it must report zero disagreements.
    """
    r = run("04_alignment_qc.py", "--aln", curated["aln"],
            "--metadata", curated["meta"], "--out-tsv", tmp_path / "q.tsv",
            "--out-png", tmp_path / "q.png")
    assert r.returncode == 0
    for line in r.stdout.splitlines():
        if "disagrees with metadata" in line or "not in metadata" in line:
            assert line.rstrip().endswith("0"), line


# --- 06 ---------------------------------------------------------------------

def test_subsample_writes_only_where_asked(curated, tmp_path):
    """
    Regression. 06 wrote a *_dates.tsv to data/processed/ regardless of
    --out-meta, and created that directory on every run even when unused.
    """
    existed_before = (ROOT / "data").exists()
    aln, meta = tmp_path / "s" / "sub.fasta", tmp_path / "s" / "sub_metadata.tsv"
    r = run("06_subsample.py", "--config", "config/pathogen/cdv.yaml",
            "--mode", "global", "--aln", curated["aln"],
            "--metadata", curated["meta"], "--target", "4",
            "--out-aln", aln, "--out-meta", meta)
    assert r.returncode == 0, r.stderr[-1500:]
    assert aln.is_file() and meta.is_file()
    assert (tmp_path / "s" / "sub_dates.tsv").is_file(), \
        "the dates file escaped the requested directory"
    # Check what THIS run created, not what merely exists. A stale data/ left
    # by the un-ported scripts is not a failure of the port, and asserting on
    # absolute absence made this fail for anyone who had run the old code.
    assert existed_before or not (ROOT / "data").exists(), \
        "06 created data/ even though explicit outputs were given"


def test_wild_groups_come_from_config_not_the_module_default(curated, tmp_path):
    """
    THE reason this needed porting. The module default is CDV's carnivore
    range. Applied to a ruminant pathogen it classifies every wild host as not
    wild, silently breaking the wild/domestic balance in subsampling.
    """
    common = ["--mode", "global", "--aln", str(curated["aln"]),
              "--metadata", str(curated["meta"]), "--target", "4"]
    cdv = run("06_subsample.py", "--config", "config/pathogen/cdv.yaml",
              *common, "--out-aln", tmp_path / "c.fasta",
              "--out-meta", tmp_path / "c_metadata.tsv")
    btv = run("06_subsample.py", "--config", "config/pathogen/btv.yaml",
              *common, "--out-aln", tmp_path / "b.fasta",
              "--out-meta", tmp_path / "b_metadata.tsv")

    def wild_line(out):
        return next(l for l in out.splitlines() if l.startswith("wild host groups:"))

    c, b = wild_line(cdv.stdout), wild_line(btv.stdout)
    assert c != b, "both pathogens got the same wild-group set"
    assert "procyonid" in c and "procyonid" not in b
    assert "wild_cervid" in b and "wild_cervid" not in c


def test_subsample_is_deterministic(curated, tmp_path):
    """A seeded subsample must be reproducible, or no build can be compared."""
    def labels(tag):
        a = tmp_path / f"{tag}.fasta"
        run("06_subsample.py", "--config", "config/pathogen/cdv.yaml",
            "--mode", "global", "--aln", curated["aln"],
            "--metadata", curated["meta"], "--target", "4", "--seed", "999",
            "--out-aln", a, "--out-meta", tmp_path / f"{tag}_metadata.tsv")
        return [l for l in a.read_text().splitlines() if l.startswith(">")]
    assert labels("r1") == labels("r2")


def test_subsample_output_labels_keep_the_contract(curated, tmp_path):
    aln = tmp_path / "sub.fasta"
    run("06_subsample.py", "--config", "config/pathogen/cdv.yaml",
        "--mode", "global", "--aln", curated["aln"],
        "--metadata", curated["meta"], "--target", "4",
        "--out-aln", aln, "--out-meta", tmp_path / "sub_metadata.tsv")
    for line in aln.read_text().splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            assert len(parts) == 3
            float(parts[2])


# --- reference anchors in QC ------------------------------------------------

@pytest.fixture(scope="module")
def curated_with_anchors(tmp_path_factory):
    """Curate with anchors on, then pad to a fake alignment including them."""
    import subprocess, yaml
    out = tmp_path_factory.mktemp("anc")
    subprocess.run([sys.executable, str(ROOT / "tests/fixtures/make_genbank.py")],
                   capture_output=True, cwd=ROOT, check=True)
    cfg = yaml.safe_load((ROOT / "config/pathogen/cdv.yaml").read_text())
    cfg["references"]["table"] = "tests/fixtures/refs_test.tsv"
    cfg["references"]["include_as_anchors"] = True
    cfgp = out / "cfg.yaml"; cfgp.write_text(yaml.safe_dump(cfg, sort_keys=False))
    r = subprocess.run(
        [sys.executable, "scripts/02_curate_metadata.py", "--config", str(cfgp),
         "--gb", "tests/fixtures/cdv_test.gb", "--out-dir", str(out)],
        capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-1500:]

    from Bio import SeqIO
    recs = list(SeqIO.parse(str(out / "H.fasta"), "fasta")) + \
           list(SeqIO.parse(str(out / "H_anchors.fasta"), "fasta"))
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
    r = run("04_alignment_qc.py", "--aln", curated_with_anchors["aln"],
            "--metadata", curated_with_anchors["meta"],
            "--out-tsv", tmp_path / "q.tsv", "--out-png", tmp_path / "q.png")
    assert r.returncode == 0, r.stderr[-1500:]


def test_anchors_are_not_reported_as_disagreements(curated_with_anchors, tmp_path):
    """
    An anchor carries its LINEAGE in the host field and a non-numeric date, by
    design. Flagging that as a mismatch would teach the reader to ignore the
    one section that catches real tip-label bugs.
    """
    r = run("04_alignment_qc.py", "--aln", curated_with_anchors["aln"],
            "--metadata", curated_with_anchors["meta"],
            "--out-tsv", tmp_path / "q.tsv", "--out-png", tmp_path / "q.png")
    assert r.returncode == 0
    out = r.stdout
    assert "reference anchors (not cross-checked): 2" in out
    for line in out.splitlines():
        if "disagrees with metadata" in line:
            assert line.rstrip().endswith("0"), line

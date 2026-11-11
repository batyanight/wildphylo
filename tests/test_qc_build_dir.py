"""
Tests for ``wildphylo align qc --build-dir``.

With ``--config --build-dir --locus`` the step takes every path it was not
given from ``BuildLayout``, so the Snakefile does not assemble paths of its
own. The fixtures here are laid out under ``tmp_path`` exactly where the layout
says, and the assertions ask the layout rather than repeating its strings.
"""

import contextlib
import io
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from wildphylo import cli
from wildphylo.config import load
from wildphylo.layout import BuildLayout

pytest.importorskip("pandas")

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "tests/fixtures/configs/cdv.yaml"
LABELS = [f"AB00000{i}|procyonid|2010.{i}00" for i in range(1, 7)]


@pytest.fixture
def build(tmp_path):
    """A build directory holding an alignment and metadata where the layout expects them."""
    layout = BuildLayout.at(tmp_path / "builds" / "cdv" / "20260101", load(CONFIG))
    aln = layout.alignment("H")
    aln.parent.mkdir(parents=True)
    # Same sequence everywhere bar one base, 60 columns, no stops in frame 0.
    aln.write_text("".join(f">{n}\n{'GCT' * 20}\n" for n in LABELS))
    meta = layout.metadata()
    meta.parent.mkdir(parents=True)
    meta.write_text(
        "accession\thost_group\tdecimal_year\n"
        + "".join(f"AB00000{i}\tprocyonid\t2010.{i}00\n" for i in range(1, 7))
    )
    return layout


def run_qc(*args):
    """Run ``wildphylo align qc`` in-process; return returncode, stdout and stderr."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = cli.main(["align", "qc", *map(str, args)])
    return SimpleNamespace(returncode=code, stdout=out.getvalue(), stderr=err.getvalue())


def test_build_dir_alone_lands_every_file_where_the_layout_says(build):
    r = run_qc("--config", CONFIG, "--build-dir", build.root, "--locus", "H")
    assert r.returncode == 0, r.stderr[-1500:]
    assert build.alignment_qc("H").is_file()
    assert build.alignment_qc_plot("H").is_file()
    # The metadata default was picked up too: the label check ran.
    assert "Labels vs metadata" in r.stdout
    written = {p.relative_to(build.root) for p in build.root.rglob("*") if p.is_file()}
    assert written == {
        build.alignment("H").relative_to(build.root),
        build.metadata().relative_to(build.root),
        build.alignment_qc("H").relative_to(build.root),
        build.alignment_qc_plot("H").relative_to(build.root),
    }


def test_explicit_flags_override_the_layout(build, tmp_path):
    tsv, png = tmp_path / "elsewhere" / "q.tsv", tmp_path / "elsewhere" / "q.png"
    r = run_qc(
        "--config", CONFIG, "--build-dir", build.root, "--locus", "H",
        "--out-tsv", tsv, "--out-png", png, "--metadata", tmp_path / "no_such.tsv",
    )  # fmt: skip
    assert r.returncode == 0, r.stderr[-1500:]
    assert tsv.is_file() and png.is_file()
    assert not build.alignment_qc("H").exists()
    assert not build.alignment_qc_plot("H").exists()
    assert "Labels vs metadata" not in r.stdout  # the explicit (missing) metadata won


def test_explicit_aln_needs_no_locus(build, tmp_path):
    tsv, png = tmp_path / "q.tsv", tmp_path / "q.png"
    r = run_qc(
        "--config", CONFIG, "--build-dir", build.root, "--aln", build.alignment("H"),
        "--out-tsv", tsv, "--out-png", png,
    )  # fmt: skip
    assert r.returncode == 0, r.stderr[-1500:]
    assert tsv.is_file()


def test_build_dir_without_locus_is_a_usage_error(build, caplog):
    caplog.set_level(logging.ERROR, logger="wildphylo")
    r = run_qc("--config", CONFIG, "--build-dir", build.root)
    assert r.returncode == 2
    assert "--locus" in caplog.text
    assert not build.alignment_qc("H").exists()


def test_build_dir_without_config_is_a_usage_error(build, caplog):
    caplog.set_level(logging.ERROR, logger="wildphylo")
    r = run_qc("--build-dir", build.root, "--locus", "H")
    assert r.returncode == 2
    assert "--config" in caplog.text


def test_no_aln_and_no_build_dir_is_a_usage_error(tmp_path, caplog):
    caplog.set_level(logging.ERROR, logger="wildphylo")
    r = run_qc("--out-tsv", tmp_path / "q.tsv")
    assert r.returncode == 2
    assert "--aln" in caplog.text


def test_locus_without_build_dir_is_a_usage_error(build, caplog):
    caplog.set_level(logging.ERROR, logger="wildphylo")
    r = run_qc("--aln", build.alignment("H"), "--locus", "H")
    assert r.returncode == 2
    assert "--build-dir" in caplog.text

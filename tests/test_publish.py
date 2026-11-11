"""``wildphylo publish dataset``: refuse a regressed build, otherwise copy, never run git."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest
import yaml

from wildphylo import cli
from wildphylo.layout import BuildLayout


@pytest.fixture
def staged(tmp_path: Path) -> dict[str, Path]:
    """A built dataset, its provenance, and an empty pathogen repo checkout."""
    dataset = tmp_path / "build" / "cdv.json"
    dataset.parent.mkdir()
    dataset.write_text(json.dumps({"version": "v2", "tree": {}}))
    provenance = tmp_path / "build" / "provenance.json"
    provenance.write_text(json.dumps({"pathogen": "cdv-1200", "build_date": "20260919"}))
    repo = tmp_path / "cdv-phylodynamics"
    repo.mkdir()
    return {
        "dataset": dataset,
        "provenance": provenance,
        "repo": repo,
        "comparison": tmp_path / "build" / "build_comparison.json",
    }


def _comparison(staged: dict[str, Path], verdict: str) -> None:
    staged["comparison"].write_text(json.dumps({"pathogen": "cdv-1200", "verdict": verdict}))


def _publish(staged: dict[str, Path], *extra: str) -> int:
    return cli.main(
        [
            "publish",
            "dataset",
            "--dataset",
            str(staged["dataset"]),
            "--comparison",
            str(staged["comparison"]),
            "--pathogen-repo",
            str(staged["repo"]),
            "--name",
            "cdv-phylodynamics",
            "--provenance",
            str(staged["provenance"]),
            *extra,
        ]
    )


def test_a_regression_is_refused_and_nothing_is_written(staged, caplog):
    _comparison(staged, "regression")
    with caplog.at_level(logging.ERROR):
        assert _publish(staged) == 1
    assert "regression" in caplog.text
    assert list(staged["repo"].iterdir()) == []


@pytest.mark.parametrize("verdict", ["ok", "review"])
def test_a_clean_or_reviewable_build_is_published(staged, verdict):
    _comparison(staged, verdict)
    assert _publish(staged) == 0
    published = staged["repo"] / "auspice" / "cdv-phylodynamics.json"
    assert published.read_bytes() == staged["dataset"].read_bytes()
    prov = staged["repo"] / "builds" / "20260919" / "provenance.json"
    assert prov.read_bytes() == staged["provenance"].read_bytes()


def test_force_publishes_despite_a_regression(staged, caplog):
    _comparison(staged, "regression")
    with caplog.at_level(logging.WARNING):
        assert _publish(staged, "--force") == 0
    assert "despite a regression" in caplog.text
    assert (staged["repo"] / "auspice" / "cdv-phylodynamics.json").is_file()


def test_a_comparison_without_a_verdict_is_refused(staged):
    staged["comparison"].write_text(json.dumps({"pathogen": "cdv-1200"}))
    assert _publish(staged) == 1
    assert list(staged["repo"].iterdir()) == []


def test_provenance_without_a_build_date_is_refused(staged):
    _comparison(staged, "ok")
    staged["provenance"].write_text(json.dumps({"pathogen": "cdv-1200"}))
    assert _publish(staged) == 1
    assert list(staged["repo"].iterdir()) == []


def test_a_missing_dataset_is_refused(staged):
    _comparison(staged, "ok")
    staged["dataset"].unlink()
    assert _publish(staged) == 1
    assert list(staged["repo"].iterdir()) == []


def test_publishing_does_not_touch_git(staged):
    _comparison(staged, "ok")
    assert _publish(staged) == 0
    assert not (staged["repo"] / ".git").exists()


# --- --build-dir defaults ---------------------------------------------------

CONFIG = Path(__file__).parent / "fixtures" / "configs" / "cdv-1200.yaml"


def _build(tmp_path: Path) -> BuildLayout:
    """A finished build laid out as BuildLayout expects, with an ``ok`` comparison."""
    layout = BuildLayout.at(tmp_path / "build", yaml.safe_load(CONFIG.read_text()))
    for path, body in (
        (layout.auspice_dataset("cdv-phylodynamics"), {"version": "v2", "tree": {}}),
        (layout.build_comparison_json(), {"verdict": "ok"}),
        (layout.provenance(), {"build_date": "20260919"}),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(body))
    return layout


def _publish_build_dir(tmp_path: Path, *extra: str) -> int:
    repo = tmp_path / "repo"
    repo.mkdir(exist_ok=True)
    return cli.main(
        [
            "publish",
            "dataset",
            "--config",
            str(CONFIG),
            "--build-dir",
            str(tmp_path / "build"),
            "--locus",
            "H",
            "--pathogen-repo",
            str(repo),
            "--name",
            "cdv-phylodynamics",
            *extra,
        ]
    )


def test_build_dir_defaults_dataset_comparison_and_provenance(tmp_path):
    layout = _build(tmp_path)
    assert _publish_build_dir(tmp_path) == 0
    repo = tmp_path / "repo"
    assert (repo / "auspice" / "cdv-phylodynamics.json").read_text() == layout.auspice_dataset(
        "cdv-phylodynamics"
    ).read_text()
    assert (repo / "builds" / "20260919" / "provenance.json").is_file()


def test_an_explicit_flag_overrides_the_build_dir_default(tmp_path):
    layout = _build(tmp_path)
    layout.build_comparison_json().write_text(json.dumps({"verdict": "regression"}))
    ok = tmp_path / "other_comparison.json"
    ok.write_text(json.dumps({"verdict": "ok"}))
    assert _publish_build_dir(tmp_path) == 1  # the default comparison refuses
    assert _publish_build_dir(tmp_path, "--comparison", str(ok)) == 0


def test_build_dir_without_a_locus_is_refused(tmp_path, caplog):
    _build(tmp_path)
    assert (
        cli.main(
            [
                "publish",
                "dataset",
                "--config",
                str(CONFIG),
                "--build-dir",
                str(tmp_path / "build"),
                "--pathogen-repo",
                str(tmp_path),
                "--name",
                "x",
            ]
        )
        == 2
    )
    assert "--locus" in caplog.text

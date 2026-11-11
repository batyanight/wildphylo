"""``wildphylo publish dataset``: refuse a regressed build, otherwise copy, never run git."""

from __future__ import annotations

import json
import logging
from typing import TYPE_CHECKING

import pytest

from wildphylo import cli

if TYPE_CHECKING:
    from pathlib import Path


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

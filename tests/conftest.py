"""Shared test fixtures and path setup."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return ROOT


@pytest.fixture(scope="session")
def configs(repo_root) -> list[Path]:
    """Every shipped pathogen config except the template."""
    return sorted(p for p in (repo_root / "config" / "pathogen").glob("*.yaml")
                  if not p.name.startswith("_"))

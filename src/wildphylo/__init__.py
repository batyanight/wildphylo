"""wildphylo: a config-driven phylodynamics pipeline for wildlife disease."""

from pathlib import Path

__version__ = "0.3.0"

__all__ = ["__version__", "workflow_path"]


def workflow_path() -> Path:
    """Absolute path to the Snakefile that defines the pipeline's DAG.

    A pathogen repo depends on this package by release tag and runs the
    pipeline through a wrapper Snakefile that ``include:``s this one, so the
    workflow has to travel with the wheel rather than sit beside the source
    checkout. It lives inside the package for that reason, and this function
    is how anything outside the package finds it: a wrapper that reaches for a
    sibling checkout works only on a developer's machine.
    """
    p = Path(__file__).resolve().parent / "workflow" / "Snakefile"
    if not p.is_file():
        raise FileNotFoundError(
            f"wildphylo's Snakefile is missing from the installed package at {p}; "
            "the install is broken or predates 0.3.0"
        )
    return p

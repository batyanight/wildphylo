"""config stage: validate pathogen configs without running anything."""

from __future__ import annotations

from typing import TYPE_CHECKING

from wildphylo import config
from wildphylo.cli import step

if TYPE_CHECKING:
    import argparse


def add_arguments_validate(parser: argparse.ArgumentParser) -> None:
    """``wildphylo config validate CONFIG [CONFIG ...] [--strict]``."""
    parser.add_argument("configs", nargs="+", help="pathogen config YAML files")
    parser.add_argument("--strict", action="store_true", help="treat warnings as failures")


@step(
    "config",
    "validate",
    "Validate pathogen configs; non-zero on errors (or warnings with --strict)",
)
def validate(args: argparse.Namespace) -> int:
    """Validate pathogen configs without running anything.

    Delegates to :func:`wildphylo.config.main` with the same arguments. Exit code is
    non-zero on errors, or on warnings too with ``--strict``.
    """
    argv = list(args.configs) + (["--strict"] if args.strict else [])
    return config.main(argv)

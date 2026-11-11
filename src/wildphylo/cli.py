"""The ``wildphylo`` command.

One subparser per stage (fetch, curate, align, tree, gates, beast, auspice,
publish, config); one subcommand per step within a stage. Stage modules under
``wildphylo.stages`` register their steps with :func:`step`; this module only
dispatches. See ``docs/PLAN-stages.md`` for the interface and ``GLOSSARY.md``
for the terms.

Every step is a function ``run(args: argparse.Namespace) -> int`` paired with
an ``add_arguments(parser)`` hook. The exit code is the gate verdict: 0 for
pass or warn, non-zero for fail, so Snakemake stops on a failing gate.
"""

from __future__ import annotations

import argparse
import importlib
import logging
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass

STAGES: tuple[str, ...] = (
    "fetch",
    "curate",
    "align",
    "tree",
    "gates",
    "beast",
    "auspice",
    "publish",
    "config",
)

AddArguments = Callable[[argparse.ArgumentParser], None]
Run = Callable[[argparse.Namespace], int]


@dataclass(frozen=True)
class Step:
    """One registered ``wildphylo <stage> <name>`` subcommand."""

    stage: str
    name: str
    help: str
    add_arguments: AddArguments
    run: Run


_REGISTRY: dict[str, dict[str, Step]] = {stage: {} for stage in STAGES}


def step(stage: str, name: str, help: str) -> Callable[[Run], Run]:  # noqa: A002 -- argparse's own name
    """Register ``run`` as the step ``<stage> <name>``.

    The decorated module must also define ``add_arguments_<name with _>`` or
    pass one explicitly through :func:`register`. Steps are looked up by the
    stage module's import, so a stage module is imported exactly once.
    """
    if stage not in _REGISTRY:
        msg = f"unknown stage {stage!r}; stages are {', '.join(STAGES)}"
        raise ValueError(msg)

    def decorate(run: Run) -> Run:
        module = sys.modules[run.__module__]
        hook = getattr(module, "add_arguments_" + name.replace("-", "_"), None)
        if hook is None:
            msg = f"{run.__module__} registers step {name!r} without add_arguments_{name}"
            raise AttributeError(msg)
        register(Step(stage, name, help, hook, run))
        return run

    return decorate


def register(item: Step) -> None:
    """Add a step to the registry; a duplicate (stage, name) is a programming error."""
    if item.name in _REGISTRY[item.stage]:
        msg = f"step {item.stage} {item.name} registered twice"
        raise ValueError(msg)
    _REGISTRY[item.stage][item.name] = item


def steps(stage: str) -> dict[str, Step]:
    """The steps registered for ``stage``, importing its module on first use."""
    importlib.import_module(f"wildphylo.stages.{stage}")
    return dict(_REGISTRY[stage])


def build_parser() -> argparse.ArgumentParser:
    """Build the full parser: every stage, every step."""
    parser = argparse.ArgumentParser(
        prog="wildphylo",
        description="Config-driven phylodynamics pipeline for wildlife disease.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    by_stage = parser.add_subparsers(dest="stage", metavar="STAGE", required=True)
    for stage in STAGES:
        stage_parser = by_stage.add_parser(stage, help=f"{stage} stage")
        by_step = stage_parser.add_subparsers(dest="step", metavar="STEP", required=True)
        for item in steps(stage).values():
            step_parser = by_step.add_parser(item.name, help=item.help, description=item.help)
            item.add_arguments(step_parser)
            step_parser.set_defaults(_run=item.run)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``wildphylo`` console script."""
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    )
    run: Run = args._run  # noqa: SLF001 -- set by build_parser, private to this module
    return run(args)


if __name__ == "__main__":
    sys.exit(main())

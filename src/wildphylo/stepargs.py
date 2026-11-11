"""``--build-dir`` for every step: the flags, and the one place that resolves them.

A step run inside a build says ``--config C --build-dir B --locus L`` and every
path it reads or writes, and every tuning value the config holds, follows from
the :class:`~wildphylo.layout.BuildLayout`. The same step run by hand gives the
paths explicitly. Both must work, and the three sources of a value must always
rank the same way:

    explicit flag  >  derived from the build layout / config  >  legacy default

This module owns that ranking and everything around it, so a step states only
*what* is derived (a table of one small function per flag) and never *how*:

* :func:`add_build_dir_arguments` registers ``--build-dir``, ``--config`` (when
  the step has none yet) and the wildcard flags ``--locus``, ``--seed``, ``--rep``.
* :func:`resolve` fills the flags the command line left out, in place on
  ``args``, validating on the way: ``--build-dir`` needs ``--config``, and each
  derived value needs the wildcard flags its path contains. Derived output
  paths get their parent directories created and every resolved value is
  logged at debug level.
* :class:`UsageError` is how any of that fails; :func:`fail` reports one the
  same way for every step (log, then a return code).

The tables themselves stay in the stage modules: which flag comes from which
layout method or config key is that step's decision::

    resolve(
        args,
        derived={"tree": lambda b: b.layout.ml_tree(b.need("locus"))},
        required=("tree",),
    )

Nothing here knows a stage by name.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml

from wildphylo.config import ConfigError
from wildphylo.config import load as load_config
from wildphylo.layout import BuildLayout

if TYPE_CHECKING:
    import argparse
    from collections.abc import Callable, Iterable, Mapping

log = logging.getLogger(__name__)

_WILDCARD_HELP = {
    "locus": "Locus whose files --build-dir selects (needed unless the paths are given).",
    "seed": "Chain seed: with --build-dir, picks the chain directory beast/<locus>/seed<SEED>/.",
    "rep": "Date-randomisation replicate: with --build-dir, picks drt/<locus>/rep<REP>/.",
}


class UsageError(Exception):
    """A flag combination the step cannot run with; the message names the fix."""


class MissingFlagError(UsageError):
    """A derived value needs a wildcard flag (``--locus``, ...) that was not given."""

    def __init__(self, flag: str) -> None:
        super().__init__(f"--{flag} is needed")
        self.flag = flag


@dataclass(frozen=True)
class Build:
    """What a derived value is computed from: the layout, the config and the flags given."""

    layout: BuildLayout
    cfg: dict[str, Any]
    args: argparse.Namespace

    def need(self, flag: str) -> str:
        """The wildcard flag ``flag`` (such as ``locus``); :class:`MissingFlagError` if absent."""
        value = getattr(self.args, flag, None)
        if value is None or value == "":
            raise MissingFlagError(flag)
        return str(value)

    def setting(self, *keys: str) -> Any:
        """``cfg[k1][k2]...``, or None when any level is absent."""
        node: Any = self.cfg
        for key in keys:
            if not isinstance(node, dict):
                return None
            node = node.get(key)
        return node

    def loci(self) -> list[str]:
        """Analysed loci as the Snakefile lists them: ``loci`` with ``analyse``, else ``locus``."""
        if "loci" in self.cfg:
            return [str(item["name"]) for item in self.cfg["loci"] if item.get("analyse")]
        return [str(self.cfg["locus"]["name"])]

    def seeds(self) -> list[str]:
        """``beast.seeds``, one BEAST chain each."""
        seeds = self.setting("beast", "seeds") or []
        if not seeds:
            msg = "beast.seeds is empty in the config"
            raise UsageError(msg)
        return [str(seed) for seed in seeds]


def add_build_dir_arguments(
    parser: argparse.ArgumentParser,
    *,
    locus: bool = True,
    seed: bool = False,
    rep: bool = False,
) -> None:
    """Add ``--build-dir``, ``--config`` (unless the step has defined it) and the wildcards.

    A step that words its own ``--config`` help must add it *before* calling
    this. ``locus``, ``seed`` and ``rep`` choose which wildcard flags the
    step's files depend on; ``--seed`` is an integer, ``--locus`` and ``--rep``
    are strings (they are path components).
    """
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=None,
        help="Build directory (builds/<id>/<YYYYMMDD>). Every path flag not given "
        "explicitly defaults to its place in the build layout; needs --config and "
        "the wildcard flags that file depends on.",
    )
    if "--config" not in parser._option_string_actions:  # noqa: SLF001 -- argparse has no public lookup
        parser.add_argument(
            "--config",
            type=Path,
            default=None,
            help="Pathogen config. Required with --build-dir, which takes tuning "
            "defaults from it too.",
        )
    if locus:
        parser.add_argument("--locus", default=None, help=_WILDCARD_HELP["locus"])
    if seed:
        parser.add_argument("--seed", type=int, default=None, help=_WILDCARD_HELP["seed"])
    if rep:
        parser.add_argument("--rep", default=None, help=_WILDCARD_HELP["rep"])


def _flag(name: str) -> str:
    return "--" + name.replace("_", "-")


def _read_config(path: Path, *, validate: bool) -> dict[str, Any]:
    if validate:
        return load_config(path)
    if not path.is_file():
        msg = f"config not found: {path}"
        raise ConfigError(msg)
    cfg = yaml.safe_load(path.read_text())
    if not isinstance(cfg, dict):
        msg = f"{path}: top level must be a mapping"
        raise ConfigError(msg)
    return cfg


def resolve(
    args: argparse.Namespace,
    *,
    derived: Mapping[str, Callable[[Build], Any]] | None = None,
    fallbacks: Mapping[str, Any] | None = None,
    required: Iterable[str] = (),
    outputs: Iterable[str] = (),
    build_dir_only: Iterable[str] = (),
    cfg: dict[str, Any] | None = None,
    validate_config: bool = True,
) -> Build | None:
    """Fill the flags the command line left out, in place; raise :class:`UsageError`.

    ``derived`` and ``fallbacks`` are keyed by ``args`` attribute (``out_json`` for
    ``--out-json``). For each, in order:

    1. a value given on the command line (anything but None) is never touched;
    2. with ``--build-dir``, ``derived[name](build)`` supplies it. A function may
       return None for "no default here", and may call ``build.need("locus")``;
       a missing wildcard flag is reported once for all the values that needed it;
    3. else ``fallbacks[name]``: the legacy default. A callable is called with
       ``args``, after every earlier name has been resolved.

    Afterwards every name in ``required`` must be set (this is what
    ``required=True`` on the argparse flag used to be). ``outputs`` names the
    derived values that are files to write: their parent directories are
    created, after all validation has passed. ``build_dir_only`` names wildcard
    flags that are meaningless without ``--build-dir``. ``cfg`` is the config
    when the step has loaded it already; otherwise ``--config`` is loaded when
    ``--build-dir`` needs it, validated unless ``validate_config`` is False (the
    YAML is then taken as written: for steps that read only a few keys and
    never needed a complete config).

    Returns the :class:`Build` (layout, config, args) used, or None without
    ``--build-dir``.
    """
    build: Build | None = None
    if args.build_dir is None:
        for flag in build_dir_only:
            if getattr(args, flag, None) is not None:
                msg = f"--{flag} only applies together with --build-dir"
                raise UsageError(msg)
    else:
        if cfg is None:
            if args.config is None:
                msg = "--build-dir needs --config (the layout depends on the pathogen id)"
                raise UsageError(msg)
            try:
                cfg = _read_config(args.config, validate=validate_config)
            except ConfigError as exc:
                raise UsageError(str(exc)) from exc
        build = Build(BuildLayout.at(args.build_dir, cfg), cfg, args)

    derived_outputs: list[Path] = []
    if build is not None:
        output_names = set(outputs)
        needing: dict[str, list[str]] = {}
        for name, make in (derived or {}).items():
            if getattr(args, name) is not None:
                continue
            try:
                value = make(build)
            except MissingFlagError as exc:
                needing.setdefault(exc.flag, []).append(_flag(name))
                continue
            except (UsageError, ValueError) as exc:
                msg = f"cannot default {_flag(name)} from --build-dir: {exc}"
                raise UsageError(msg) from exc
            if value is None:
                continue
            setattr(args, name, value)
            log.debug("%s = %s (from --build-dir)", _flag(name), value)
            if name in output_names:
                derived_outputs.append(Path(value))
        if needing:
            msg = "; ".join(
                f"--build-dir needs --{flag} to derive {', '.join(names)} (or give them)"
                for flag, names in needing.items()
            )
            raise UsageError(msg)

    for name, default in (fallbacks or {}).items():
        if getattr(args, name) is None:
            value = default(args) if callable(default) else default
            setattr(args, name, value)
            log.debug("%s = %s (default)", _flag(name), value)

    missing = [_flag(name) for name in required if getattr(args, name) is None]
    if missing:
        where = "without --build-dir" if args.build_dir is None else "(--build-dir has no default)"
        verb = "is required" if len(missing) == 1 else "required"
        msg = f"{', '.join(missing)} {verb} {where}"
        raise UsageError(msg)

    for path in derived_outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
    return build


def fail(error: UsageError, *, code: int = 2) -> int:
    """Report a usage error the way every step does: log the message, return the exit code."""
    log.error("%s", error)
    return code

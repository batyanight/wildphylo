"""gates stage: the checks that decide whether the data can support the next step.

Each gate inspects data, writes a report and returns a verdict (see GLOSSARY.md).
Steps:

* ``temporal``: root-to-tip regression and the clock gate (run before BEAST).
* ``congruence``: per-segment tree topologies, before concatenating segments.
* ``convergence``: ESS and between-chain agreement of the BEAST chains.
* ``trait-signal``: tip-label randomisation test for discrete-trait analyses.
* ``drt``: the date-randomisation test, summarised and judged.
* ``compare``: what changed since the previous build, and is it a finding.

Exit codes follow each gate's own rule (documented on the step function); the
temporal, trait-signal and drt gates return 0 on a ``fail`` verdict unless
``--exit-nonzero-on-fail`` is given, so the report always survives.

Every step takes ``--build-dir`` (with ``--config``): each input, output and
tuning flag not given explicitly then defaults to the build layout
(:class:`wildphylo.layout.BuildLayout`) or the config. See :func:`_resolve`.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from wildphylo import newick
from wildphylo import traitsignal as tsig
from wildphylo.beastlog import (
    NON_PARAMETERS,
    ess,
    find_clock_parameter,
    hpd,
    overlap,
    read_beast_log,
)
from wildphylo.cli import step
from wildphylo.config import ConfigError, load
from wildphylo.layout import BuildLayout
from wildphylo.preflight import Checks, check_convergence
from wildphylo.segments import assess_congruence
from wildphylo.temporal import assess as assess_temporal_signal

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

    from Bio.Phylo.BaseTree import Tree

log = logging.getLogger(__name__)


# ===========================================================================
# --build-dir: defaults from the build layout and the config
# ===========================================================================


class _UsageError(Exception):
    """A flag combination a step cannot act on; the message names the fix."""


@dataclass(frozen=True)
class _Build:
    """What ``--build-dir`` resolves against: the layout, the config and the flags given."""

    layout: BuildLayout
    cfg: dict[str, Any]
    args: argparse.Namespace

    def need(self, flag: str) -> str:
        """The value of a per-locus flag such as ``locus``; an error when it was not given."""
        value = getattr(self.args, flag, None)
        if not value:
            msg = f"--{flag} is needed"
            raise _UsageError(msg)
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
            raise _UsageError(msg)
        return [str(seed) for seed in seeds]


def _add_build_arguments(
    parser: argparse.ArgumentParser, *, config: bool = True, locus: bool = True
) -> None:
    """``--build-dir``, ``--config`` (unless the step already has one) and ``--locus``."""
    parser.add_argument(
        "--build-dir",
        type=Path,
        help="Build directory (builds/<id>/<YYYYMMDD>). Every path flag not given "
        "explicitly defaults to its place in the build layout; needs --config.",
    )
    if config:
        parser.add_argument(
            "--config",
            type=Path,
            help="Pathogen config. Required with --build-dir, which takes tuning "
            "defaults from it too.",
        )
    if locus:
        parser.add_argument(
            "--locus",
            help="Locus whose files --build-dir selects (needed unless paths are given).",
        )


def _given(value: object) -> bool:
    return value is not None and value != []


def _flag(name: str) -> str:
    return "--" + name.replace("_", "-")


def _resolve(
    args: argparse.Namespace,
    *,
    derived: dict[str, Callable[[_Build], Any]] | None = None,
    required: Iterable[str] = (),
    fallbacks: dict[str, Any] | None = None,
) -> bool:
    """Fill the flags a step was not given; False (error logged) when it cannot run.

    Precedence for each flag: the value on the command line; else, with
    ``--build-dir``, ``derived[name](build)`` (a layout path or a config
    value; None means "no default"); else ``fallbacks[name]``. A flag in
    ``required`` still unset at the end is an error, which is what
    ``required=True`` on the argparse flag used to be. Without ``--build-dir``
    nothing is derived, so behaviour is exactly the explicit flags.
    """
    try:
        if args.build_dir is not None:
            if args.config is None:
                msg = "--build-dir needs --config"
                raise _UsageError(msg)  # noqa: TRY301 -- one handler for every usage error
            try:
                cfg = load(args.config)
            except ConfigError as e:
                raise _UsageError(str(e)) from e
            build = _Build(BuildLayout.at(args.build_dir, cfg), cfg, args)
            for name, default in (derived or {}).items():
                if _given(getattr(args, name)):
                    continue
                try:
                    value = default(build)
                except (_UsageError, ValueError) as e:
                    msg = f"cannot default {_flag(name)} from --build-dir: {e}"
                    raise _UsageError(msg) from e
                if value is not None:
                    setattr(args, name, value)
                    log.debug("%s defaults to %s", _flag(name), value)
                    if name.startswith("out"):
                        # Snakemake made output directories for us; a bare --build-dir run must.
                        Path(value).parent.mkdir(parents=True, exist_ok=True)
        for name, value in (fallbacks or {}).items():
            if getattr(args, name) is None:
                setattr(args, name, value)
        missing = [_flag(name) for name in required if not _given(getattr(args, name))]
        if missing:
            msg = f"missing {', '.join(missing)} (or give --build-dir with --config)"
            raise _UsageError(msg)  # noqa: TRY301
    except _UsageError as e:
        log.error("%s", e)  # noqa: TRY400 -- message only, no traceback
        return False
    return True


# ===========================================================================
# gates temporal
# ===========================================================================


def root_to_tip_distances(tree: Tree) -> dict[str, float]:
    """Distance from root to every tip, in substitutions per site."""
    # Bio.Phylo's tree methods ship without annotations (see the mypy override in the report).
    return {
        t.name: tree.distance(tree.root, t)  # type: ignore[no-untyped-call]
        for t in tree.get_terminals()  # type: ignore[no-untyped-call]
    }


def date_from_label(label: str) -> float | None:
    """Decimal year from the last ``|``-separated field of a tip label, or None."""
    try:
        return float(label.rsplit("|", 1)[-1])
    except (ValueError, IndexError):
        return None


def midpoint_root(tree: Tree) -> Tree:
    """Midpoint-root ``tree``; keep it as given if rooting fails."""
    # Midpoint rooting can fail on degenerate trees; keep the tree as given.
    with contextlib.suppress(Exception):
        tree.root_at_midpoint()  # type: ignore[no-untyped-call]
    return tree


def add_arguments_temporal(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo gates temporal``."""
    _add_build_arguments(parser)
    parser.add_argument("--tree", help="Required unless --build-dir.")
    parser.add_argument("--tree-format", default="newick")
    parser.add_argument("--metadata", help="TSV with columns name/strain and decimal_year")
    parser.add_argument(
        "--permutations", type=int, help="Default 1000, or temporal_signal.permutations."
    )
    parser.add_argument(
        "--min-r2", type=float, help="Default 0.10, or temporal_signal.min_r_squared."
    )
    parser.add_argument(
        "--min-span", type=float, help="Default 5.0, or temporal_signal.min_span_years."
    )
    parser.add_argument("--no-midpoint-root", action="store_true")
    parser.add_argument("--out-json", help="Required unless --build-dir.")
    parser.add_argument("--out-md")
    # A FAIL verdict is a result, not a crash. Exiting non-zero made Snakemake
    # treat the rule as failed and DELETE the report -- which is the single
    # most valuable artefact precisely when the gate fails. The workflow reads
    # the verdict from the JSON and decides what to schedule; it does not need
    # the exit code. Standalone shell use can opt back in.
    parser.add_argument(
        "--exit-nonzero-on-fail",
        action="store_true",
        help="Return exit code 1 when the verdict is 'fail'. Off by "
        "default so the report always survives.",
    )
    parser.add_argument("--out-plot")


def _read_dates(metadata: str) -> tuple[dict[str, Any], dict[str, Any]]:
    """Decimal years and date precisions per tip name from a metadata TSV."""
    import pandas as pd

    dates: dict[str, Any] = {}
    precisions: dict[str, Any] = {}
    df = pd.read_csv(metadata, sep="\t")
    namecol = next(
        (c for c in ("name", "strain", "accession", "taxon") if c in df.columns), df.columns[0]
    )
    for _, r in df.iterrows():
        dates[str(r[namecol])] = r.get("decimal_year")
        precisions[str(r[namecol])] = r.get("date_precision", "unknown")
    return dates, precisions


def _plot_root_to_tip(
    path: str, xs: list[float], ys: list[float], r: Any, verdict: str, perm_p: float
) -> None:
    """Scatter of divergence against date with the regression line; skip on any error."""
    try:
        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.scatter(xs, ys, s=18, alpha=0.7, edgecolor="none")
        lo, hi = min(xs), max(xs)
        ax.plot([lo, hi], [r.intercept + r.slope * lo, r.intercept + r.slope * hi], lw=1.5)
        ax.set_xlabel("Sampling date (decimal year)")
        ax.set_ylabel("Root-to-tip divergence (subs/site)")
        ax.set_title(
            f"Root-to-tip regression — {verdict.upper()}  (R²={r.r_squared:.3f}, p={perm_p:.3f})"
        )
        fig.tight_layout()
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150)
        plt.close(fig)
    except Exception as e:  # noqa: BLE001 -- a missing plot must never fail the gate
        log.info("NOTE: plot skipped (%s)", e)


@step("gates", "temporal", "Root-to-tip regression and clock gate (run before BEAST)")
def temporal(args: argparse.Namespace) -> int:
    """Automated root-to-tip regression and clock gate.

    Replaces the manual TempEst step. Tip labels are expected as
    ACCESSION|group|decimal_year, or dates can be supplied via ``--metadata``.

    Exit code is 0 on pass/weak and also on fail unless ``--exit-nonzero-on-fail``
    is given; 1 when fewer than three tips carry a usable date; 2 for unusable flags.

    With ``--build-dir`` and ``--config`` (and ``--locus``) the tree, metadata,
    outputs and thresholds default to the build layout and ``temporal_signal``.
    """
    if not _resolve(
        args,
        derived={
            "tree": lambda b: b.layout.ml_tree(b.need("locus")),
            "metadata": lambda b: b.layout.metadata(),
            "out_json": lambda b: b.layout.temporal_report(b.need("locus")),
            "out_md": lambda b: b.layout.temporal_report_md(b.need("locus")),
            "out_plot": lambda b: b.layout.temporal_plot(b.need("locus")),
            "permutations": lambda b: b.setting("temporal_signal", "permutations"),
            "min_r2": lambda b: b.setting("temporal_signal", "min_r_squared"),
            "min_span": lambda b: b.setting("temporal_signal", "min_span_years"),
        },
        required=("tree", "out_json"),
        fallbacks={"permutations": 1000, "min_r2": 0.10, "min_span": 5.0},
    ):
        return 2

    from Bio import Phylo

    tree = Phylo.read(args.tree, args.tree_format)  # type: ignore[attr-defined,no-untyped-call]
    if not args.no_midpoint_root:
        tree = midpoint_root(tree)

    dists = root_to_tip_distances(tree)

    dates: dict[str, Any] = {}
    precisions: dict[str, Any] = {}
    if args.metadata:
        dates, precisions = _read_dates(args.metadata)

    xs: list[float] = []
    ys: list[float] = []
    ps: list[str] = []
    unresolved: list[str] = []
    for name, d in dists.items():
        dec = dates.get(name)
        if dec is None:
            dec = date_from_label(name)
        if dec is None:
            unresolved.append(name)
            continue
        xs.append(float(dec))
        ys.append(d)
        ps.append(precisions.get(name, "unknown"))

    if len(xs) < 3:
        log.error("only %d tips carry a usable date; cannot regress.", len(xs))
        return 1
    if unresolved:
        log.info(
            "NOTE: %d tips had no parseable date and were excluded from the regression (e.g. %s)",
            len(unresolved),
            unresolved[:3],
        )

    rep = assess_temporal_signal(
        xs,
        ys,
        precisions=ps,
        n_perm=args.permutations,
        min_r2=args.min_r2,
        min_span=args.min_span,
    )

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(rep.to_dict(), indent=2))

    r = rep.regressions[0]
    perm_p = rep.permutation_p if rep.permutation_p is not None else float("nan")
    lines = [
        "# Temporal signal assessment",
        "",
        f"**Verdict: {rep.verdict.upper()}**",
        "",
        "| Metric | Value |",
        "|---|---|",
        f"| Tips | {rep.n_tips} |",
        (
            f"| Sampling span | {rep.earliest:.2f} \u2013 {rep.latest:.2f} "
            f"({rep.span_years:.1f} yr) |"
        ),
        f"| Root-to-tip slope | {r.slope:.3e} subs/site/yr |",
        f"| R² | {r.r_squared:.4f} |",
        f"| Correlation | {r.correlation:.4f} |",
        f"| Implied root date | {r.x_intercept:.1f} |"
        if r.x_intercept
        else "| Implied root date | undefined |",
        f"| Date-shuffling p | {perm_p:.4f} ({rep.permutation_n} perms) |",
        f"| Year-bin occupancy | {rep.occupancy:.0%} |",
        f"| Densest 10-yr window | {rep.concentration:.0%} of tips |",
        f"| Year-only or range dates | {rep.frac_year_only:.0%} |",
        "",
    ]
    if rep.notes:
        lines += ["## Cautions", ""] + [f"- {n}" for n in rep.notes] + [""]
    if args.out_md:
        Path(args.out_md).write_text("\n".join(lines))
    log.info("\n%s", "\n".join(lines))

    if args.out_plot:
        _plot_root_to_tip(args.out_plot, xs, ys, r, rep.verdict, perm_p)

    if rep.verdict == "fail":
        log.warning(
            "Verdict is FAIL. The report above has still been written — read "
            "it rather than re-running. A failed gate is a finding about the "
            "data, not a broken step."
        )
        if args.exit_nonzero_on_fail:
            return 1
    return 0


# ===========================================================================
# gates congruence
# ===========================================================================


def add_arguments_congruence(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo gates congruence``."""
    _add_build_arguments(parser, locus=False)
    parser.add_argument(
        "--trees",
        nargs="+",
        help="Segment trees. Required unless --build-dir, which takes the ML tree of "
        "every analysed locus.",
    )
    parser.add_argument(
        "--names",
        nargs="+",
        help="Segment names, in the same order as --trees. Defaults to the tree "
        "filenames, or to the locus names when --trees came from --build-dir.",
    )
    parser.add_argument("--tree-format", default="newick")
    parser.add_argument(
        "--collapse-support",
        type=float,
        help="Ignore splits below this support before comparing. "
        "Set to a negative number to disable. Default 70, or segments.collapse_support.",
    )
    parser.add_argument("--congruent-max-rf", type=float, default=0.15)
    parser.add_argument("--decoupled-min-rf", type=float, default=0.40)
    parser.add_argument("--neighbour-k", type=int, default=5)
    parser.add_argument(
        "--reassortant-threshold",
        type=float,
        help="Default 0.6, or segments.reassortant_threshold.",
    )
    parser.add_argument(
        "--require-concatenation",
        action="store_true",
        help="Fail if the segments cannot legitimately be joined.",
    )
    parser.add_argument("--out-json", help="Required unless --build-dir.")
    parser.add_argument("--out-md")
    parser.add_argument("--out-tsv", help="Reassortant candidates, one per row.")


@step("gates", "congruence", "Compare per-segment tree topologies before concatenating")
def congruence(args: argparse.Namespace) -> int:
    """Compare segment topologies before treating them as one history.

    Run after the per-segment ML trees exist and before subsampling. For an
    unsegmented pathogen this step does not run at all.

    Exit code is 0 when concatenation is legitimate or the verdict is merely
    mixed, and 1 when segments are decoupled AND ``--require-concatenation`` is
    given. Being decoupled is not itself an error -- for BTV it is the expected
    result. Exit code 2 means unusable input.

    With ``--build-dir`` and ``--config`` the trees are the ML trees of the
    analysed loci (named by locus), the thresholds come from ``segments`` and the
    reports go under ``segments/``.
    """
    # Locus names label the trees only when the trees came from the layout too.
    trees_from_layout = args.trees is None
    if not _resolve(
        args,
        derived={
            "trees": lambda b: [b.layout.ml_tree(name) for name in b.loci()],
            "names": lambda b: b.loci() if trees_from_layout else None,
            "collapse_support": lambda b: b.setting("segments", "collapse_support"),
            "reassortant_threshold": lambda b: b.setting("segments", "reassortant_threshold"),
            "out_json": lambda b: b.layout.congruence_report(),
            "out_md": lambda b: b.layout.congruence_md(),
            "out_tsv": lambda b: b.layout.reassortant_candidates(),
        },
        required=("trees", "out_json"),
        fallbacks={"collapse_support": 70.0, "reassortant_threshold": 0.6},
    ):
        return 2
    names = args.names or [Path(t).stem for t in args.trees]
    if len(names) != len(args.trees):
        log.error("%d trees but %d names", len(args.trees), len(names))
        return 2
    if len(args.trees) < 2:
        log.error("at least two segment trees are needed; congruence cannot be assessed from one")
        return 2

    from Bio import Phylo

    trees = {}
    for name, path in zip(names, args.trees, strict=True):
        try:
            trees[name] = Phylo.read(path, args.tree_format)  # type: ignore[attr-defined,no-untyped-call]
        except Exception as e:  # noqa: BLE001 -- any parser failure is unusable input
            log.error("could not read %s: %s", path, e)  # noqa: TRY400 -- message only, no traceback
            return 2

    support = None if args.collapse_support < 0 else args.collapse_support
    res = assess_congruence(
        trees,
        collapse_support=support,
        congruent_max_rf=args.congruent_max_rf,
        decoupled_min_rf=args.decoupled_min_rf,
        neighbour_k=args.neighbour_k,
        reassortant_threshold=args.reassortant_threshold,
    )

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(res.to_dict(), indent=2))

    lines = [
        "# Segment congruence",
        "",
        (
            f"**Verdict: {res.verdict.upper()}** — "
            f"concatenation {'permitted' if res.concatenate else 'refused'}"
        ),
        "",
        "| Metric | Value |",
        "|---|---|",
        f"| Segments | {', '.join(res.segments)} |",
        f"| Taxa present in every segment | {res.n_shared_taxa} |",
        f"| Mean normalised RF | {res.mean_rf:.4f} |",
        f"| Max normalised RF | {res.max_rf:.4f} |",
        f"| Reassortant candidates | {len(res.reassortant_candidates)} |",
        "",
    ]
    if res.pairwise_rf:
        lines += ["## Pairwise topological distance", "", "| Pair | normalised RF |", "|---|---|"]
        # Pair keys contain "|", which would break the Markdown table.
        lines += [
            f"| {k.replace('|', ' vs ')} | {v:.4f} |" for k, v in sorted(res.pairwise_rf.items())
        ]
        lines.append("")
    if res.reassortant_candidates:
        lines += [
            "## Reassortant candidates",
            "",
            (
                "Each sits among different relatives depending on the "
                "segment. In a discrete-trait host analysis such a taxon "
                "produces a host transition that reflects reassortment "
                "rather than transmission."
            ),
            "",
            "| Taxon | Incongruence |",
            "|---|---|",
        ]
        lines += [
            f"| {c['taxon']} | {c['incongruence']:.3f} |" for c in res.reassortant_candidates[:25]
        ]
        if len(res.reassortant_candidates) > 25:
            lines.append(f"| … and {len(res.reassortant_candidates) - 25} more | |")
        lines.append("")
    if res.notes:
        lines += ["## Notes", ""] + [f"- {n}" for n in res.notes] + [""]

    # Write every output BEFORE reporting. If stdout is a pipe that closes
    # early (a `| head` in a terminal, a truncated CI log), a report that comes
    # first can take the process down and the files are never written --
    # silently, since the exit status looks fine.
    if args.out_md:
        Path(args.out_md).write_text("\n".join(lines))
    if args.out_tsv:
        Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
        with Path(args.out_tsv).open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["taxon", "incongruence"])
            for c in res.reassortant_candidates:
                w.writerow([c["taxon"], c["incongruence"]])

    log.info("\n%s", "\n".join(lines))

    if args.require_concatenation and not res.concatenate:
        log.error("concatenation was required but the segments do not support it.")
        return 1
    return 0


# ===========================================================================
# gates convergence
# ===========================================================================
#
# ESS per parameter follows Tracer's algorithm, so the numbers match what you
# would read off Tracer or loganalyser rather than being a second opinion: a
# gate that disagrees with the tool the user checks it against gets overridden.
#
# OVL is the overlapping coefficient between chains, per parameter: the area
# shared by their marginal posteriors, 1.0 for identical and 0.0 for disjoint.
# Chosen over HPD comparison because two bimodal posteriors can have
# near-identical 95% intervals while placing their mass in different places,
# and over Gelman-Rubin because R-hat assumes roughly normal marginals, which
# clock rates and tree heights are not.
#
# The combined log pools post-burn-in samples from every chain with the Sample
# column renumbered. It is written here rather than by logcombiner because the
# pooling must happen only AFTER the chains are shown to agree.


def analyse_convergence(paths: list[Path], burnin: float, min_ess: int) -> dict[str, Any]:
    """ESS, between-chain overlap and HPDs for the BEAST logs in ``paths``."""
    chains: list[dict[str, Any]] = []
    for p in paths:
        header, data = read_beast_log(p)
        cut = round(burnin * len(data))
        chains.append(
            {
                "path": str(p),
                "header": header,
                "data": data[cut:],
                "n_total": len(data),
                "n_kept": len(data) - cut,
            }
        )
        if len(data) - cut < 10:
            raise ValueError(
                f"{p}: only {len(data) - cut} samples survive a "
                f"{burnin:.0%} burn-in; the chain is too short to assess"
            )

    # Every chain must log the same parameters. A mismatch means the XMLs
    # differ, which is a different analysis, not a different seed -- exactly
    # the mistake of editing one seed's XML and forgetting the other.
    ref = [h for h in chains[0]["header"] if h not in NON_PARAMETERS]
    for c in chains[1:]:
        other = [h for h in c["header"] if h not in NON_PARAMETERS]
        if other != ref:
            only_a = sorted(set(ref) - set(other))
            only_b = sorted(set(other) - set(ref))
            raise ValueError(
                f"chains log different parameters, so they are not replicates "
                f"of one analysis. Only in {chains[0]['path']}: {only_a}; "
                f"only in {c['path']}: {only_b}"
            )

    cols = {name: i for i, name in enumerate(chains[0]["header"])}

    # ESS is reported on the POOLED post-burn-in samples, which is what
    # loganalyser on a combined log gives, plus per chain so a single bad
    # chain is visible rather than being rescued by the other.
    ess_pooled: dict[str, float] = {}
    ess_per_chain: dict[str, list[float]] = {}
    hpds: dict[str, list[float]] = {}
    for name in ref:
        j = cols[name]
        per = [ess(c["data"][:, j]) for c in chains]
        ess_per_chain[name] = [round(e, 1) for e in per]
        # The pooled ESS of independent chains is their sum. Concatenating and
        # running the autocorrelation estimator over the join would read the
        # discontinuity as signal and understate it.
        ess_pooled[name] = round(float(sum(per)), 1)
        pooled = np.concatenate([c["data"][:, j] for c in chains])
        lo, hi = hpd(pooled)
        hpds[name] = [round(float(pooled.mean()), 6), round(lo, 6), round(hi, 6)]

    ovl: dict[str, float] = {}
    if len(chains) > 1:
        for name in ref:
            j = cols[name]
            pairwise = [
                overlap(chains[a]["data"][:, j], chains[b]["data"][:, j])
                for a in range(len(chains))
                for b in range(a + 1, len(chains))
            ]
            ovl[name] = round(min(pairwise), 4)

    return {
        "chains": chains,
        "params": ref,
        "cols": cols,
        "ess": ess_pooled,
        "ess_per_chain": ess_per_chain,
        "overlap": ovl,
        "hpd": hpds,
        "min_ess": min_ess,
        "burnin": burnin,
    }


def write_combined_log(res: dict[str, Any], out: Path) -> int:
    """Pool post-burn-in samples with a renumbered Sample column; return the row count."""
    chains = res["chains"]
    header = chains[0]["header"]
    sample_col = 0 if header[0] in NON_PARAMETERS else None
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out.open("w") as fh:
        fh.write(
            "# pooled post-burn-in samples from: " + ", ".join(c["path"] for c in chains) + "\n"
        )
        fh.write(
            "# burn-in " + f"{res['burnin']:.0%}" + " applied within each chain before pooling\n"
        )
        fh.write("\t".join(header) + "\n")
        for c in chains:
            for row in c["data"]:
                vals = list(row)
                if sample_col is not None:
                    vals[sample_col] = n
                fh.write(
                    "\t".join(
                        (f"{int(v)}" if i == sample_col else f"{v:.6g}")
                        for i, v in enumerate(vals)
                    )
                    + "\n"
                )
                n += 1
    return n


def convergence_report(res: dict[str, Any], checks: Checks) -> str:
    """Fixed-width text report of the convergence analysis."""
    width = 78
    lines = [
        "=" * width,
        "CONVERGENCE",
        "=" * width,
        f"chains        : {len(res['chains'])}",
        f"burn-in       : {res['burnin']:.0%} within each chain",
    ]
    lines.extend(
        f"  {c['path']}: {c['n_kept']} of {c['n_total']} samples kept" for c in res["chains"]
    )
    lines += [
        "",
        f"  {'parameter':<28}{'ESS':>10}{'per chain':>20}{'OVL':>8}",
        "  " + "-" * (width - 4),
    ]
    worst = sorted(res["params"], key=lambda p: res["ess"][p])
    for name in worst:
        per = "/".join(f"{e:.0f}" for e in res["ess_per_chain"][name])
        o = res["overlap"].get(name)
        flag = ""
        if res["ess"][name] < res["min_ess"]:
            flag = "  <- ESS"
        elif o is not None and o < 0.5:
            flag = "  <- chains disagree"
        lines.append(
            f"  {name[:28]:<28}{res['ess'][name]:>10.0f}{per:>20}"
            f"{('' if o is None else f'{o:.2f}'):>8}{flag}"
        )
    lines += [
        "",
        "  ESS is the pooled figure (chains summed). OVL is the smallest",
        "  pairwise overlap of the chains' marginal posteriors: 1.0 means",
        "  they agree, below 0.5 means at least one chain is somewhere the",
        "  others never went.",
        "",
        checks.render(),
        "=" * width,
    ]
    return "\n".join(lines)


def add_arguments_convergence(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo gates convergence``."""
    _add_build_arguments(parser, config=False)
    parser.add_argument(
        "--logs",
        nargs="+",
        help="BEAST .log files, one per chain. Required unless --build-dir, which "
        "takes one per beast.seeds entry.",
    )
    parser.add_argument(
        "--burnin",
        type=float,
        help="Fraction discarded from the head of EACH chain. "
        "Applied per chain, not to the pooled set. Default 0.10, or "
        "beast.burnin_fraction with --build-dir.",
    )
    parser.add_argument(
        "--min-ess",
        type=int,
        help="Default 200. A config's beast.min_ess overrides this flag, except with "
        "--build-dir, where the flag wins.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="Pathogen config, for min_ess and the date-policy "
        "check in preflight.check_convergence.",
    )
    parser.add_argument("--out-json", type=Path)
    parser.add_argument("--out-combined", type=Path)
    parser.add_argument("--out-md", type=Path)
    parser.add_argument(
        "--allow-fail",
        action="store_true",
        help="Report and exit 0 even on failure. For looking at a chain that is still running.",
    )


@step("gates", "convergence", "ESS and between-chain agreement of the BEAST chains")
def convergence(args: argparse.Namespace) -> int:
    """ESS and between-chain agreement, checked automatically.

    Replaces ``loganalyser -b 10`` run by hand, and adds the check loganalyser
    does not do: whether the chains agree WITH EACH OTHER. A single chain can
    look immaculate while sitting in a region of parameter space the other
    chain never visits. ESS measures how well a chain sampled where it went,
    not whether it went to the right place.

    Exit code is non-zero when the gate fails, so a workflow stops rather than
    handing an unconverged chain to TreeAnnotator.

    With ``--build-dir`` and ``--config`` (and ``--locus``) the logs are the
    chains of ``beast.seeds``, the thresholds come from ``beast`` and the
    reports go to the locus's BEAST directory.
    """
    if not _resolve(
        args,
        derived={
            "logs": lambda b: [b.layout.beast_log(b.need("locus"), s) for s in b.seeds()],
            "out_json": lambda b: b.layout.convergence_report(b.need("locus")),
            "out_combined": lambda b: b.layout.combined_log(b.need("locus")),
            "burnin": lambda b: b.setting("beast", "burnin_fraction"),
        },
        required=("logs",),
        fallbacks={"burnin": 0.10},
    ):
        return 2
    paths = [Path(p) for p in args.logs]
    missing = [p for p in paths if not p.is_file()]
    if missing:
        log.error("no such log file: %s", missing[0])
        return 1

    cfg: dict[str, Any] = {}
    if args.config:
        try:
            cfg = load(args.config)
        except ConfigError as e:
            log.error("%s", e)  # noqa: TRY400 -- message only, no traceback
            return 1
    # Without --build-dir the config's value beats the flag, as it always did;
    # with it, an explicit flag wins like every other flag.
    config_min_ess = cfg.get("beast", {}).get("min_ess")
    order = (args.min_ess, config_min_ess) if args.build_dir else (config_min_ess, args.min_ess)
    cfg.setdefault("beast", {})["min_ess"] = next((v for v in order if v), 200)

    try:
        res = analyse_convergence(paths, args.burnin, cfg["beast"]["min_ess"])
    except ValueError as e:
        log.error("%s", e)  # noqa: TRY400 -- message only, no traceback
        return 1

    checks = check_convergence(cfg, res["ess"], len(res["chains"]), res["overlap"] or None)
    text = convergence_report(res, checks)
    log.info("\n%s", text)

    if args.out_combined:
        n = write_combined_log(res, args.out_combined)
        log.info("wrote %s  (%d pooled samples)", args.out_combined, n)

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(
            json.dumps(
                {
                    "verdict": "pass" if checks.ok else "fail",
                    "n_chains": len(res["chains"]),
                    "burnin": res["burnin"],
                    "min_ess": res["min_ess"],
                    "logs": [c["path"] for c in res["chains"]],
                    "samples_kept": [c["n_kept"] for c in res["chains"]],
                    "ess": res["ess"],
                    "ess_per_chain": res["ess_per_chain"],
                    "overlap": res["overlap"],
                    "hpd": res["hpd"],
                    "checks": checks.to_dict(),
                },
                indent=2,
            )
        )
        log.info("wrote %s", args.out_json)

    if args.out_md:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text("```\n" + text + "\n```\n")
        log.info("wrote %s", args.out_md)

    if not checks.ok and not args.allow_fail:
        log.error("exiting non-zero: the chains have not converged.")
        return 1
    return 0


# ===========================================================================
# gates trait-signal
# ===========================================================================


def states_from_traits(path: Path) -> dict[str, str]:
    """Read the traits file ``beast xml`` writes: a two-column TSV with a header."""
    out: dict[str, str] = {}
    with Path(path).open(newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    for row in rows:
        if len(row) < 2 or row[0].lower() in ("traits", "name", "strain", "taxon"):
            continue
        out[row[0].strip()] = row[1].strip()
    if not out:
        raise ValueError(f"{path} yielded no taxon/state pairs")
    return out


def states_from_labels(labels: list[str], field: int, sep: str) -> dict[str, str]:
    """Tip state from field ``field`` of each label split on ``sep``."""
    out: dict[str, str] = {}
    bad: list[str] = []
    for lab in labels:
        parts = lab.split(sep)
        if len(parts) > field:
            out[lab] = parts[field].strip()
        else:
            bad.append(lab)
    if bad:
        raise ValueError(
            f"{len(bad)} tip labels have no field {field} when split on "
            f"{sep!r}; first is {bad[0]!r}. Use --traits instead, or set "
            f"--label-field / --label-sep to match your labels."
        )
    return out


def add_arguments_trait_signal(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo gates trait-signal``."""
    _add_build_arguments(parser, config=False)
    parser.add_argument(
        "--trees",
        nargs="+",
        help="Posterior tree files, one per chain. Burn-in is "
        "applied within each file, not to the pooled set. Required unless "
        "--build-dir, which takes the trait trees of every beast.seeds chain.",
    )
    parser.add_argument(
        "--traits", type=Path, help="Two-column TSV of taxon and state (what 07 writes)."
    )
    parser.add_argument(
        "--label-field",
        type=int,
        default=1,
        help="Field of the tip label holding the state, 0-based. "
        "Used when --traits is not given. Default 1, which "
        "is the group in ACCESSION|group|year.",
    )
    parser.add_argument("--label-sep", default="|")
    parser.add_argument(
        "--config",
        type=Path,
        help="Pathogen config. hosts.dta_states is applied to "
        "collapse host groups onto the trait states BEAST "
        "actually estimated rates between.",
    )
    parser.add_argument(
        "--burnin", type=float, help="Default 0.10, or beast.burnin_fraction with --build-dir."
    )
    parser.add_argument("--thin", type=int, default=1)
    parser.add_argument(
        "--max-trees",
        type=int,
        default=1000,
        help="Cap on retained trees after burn-in and thinning. "
        "The statistics converge well before a full "
        "posterior; more trees buy precision on the mean, "
        "not resolution on p, which --permutations sets.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        help="Null replicates. The smallest p obtainable is 1/(permutations+1). "
        "Default 1000, or beast.discrete_trait.signal_permutations with --build-dir.",
    )
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--out-json", type=Path)
    parser.add_argument("--out-md", type=Path)
    parser.add_argument(
        "--exit-nonzero-on-fail",
        action=argparse.BooleanOptionalAction,
        help="Exit 1 when the verdict is fail, so a workflow "
        "stops rather than publishing an unsupported table. Off by default; with "
        "--build-dir, on unless beast.discrete_trait.on_no_signal is warn.",
    )


@step("gates", "trait-signal", "Tip-label randomisation gate for discrete-trait analyses")
def trait_signal(args: argparse.Namespace) -> int:
    """Tip-label randomisation gate for discrete-trait analyses.

    Run this before reporting ANY number out of a discrete-trait reconstruction.
    A DTA always returns a transition table; this decides whether that table
    describes the phylogeny or the label distribution.

    Tip states come, in order of preference, from ``--traits``, then from the
    tip label field given by ``--label-field`` (labels here are
    ACCESSION|group|year), and are then collapsed through ``hosts.dta_states``
    from ``--config`` so that the states tested are exactly the states BEAST
    estimated rates between. Testing the uncollapsed host groups would answer a
    different question from the one the analysis asked.

    The statistics use topology and tip labels only, never BEAST's
    reconstructed ancestral states; see :mod:`wildphylo.traitsignal`.

    With ``--build-dir`` and ``--config`` (and ``--locus``) the trees are the
    trait trees of every ``beast.seeds`` chain, the tip states come from the
    tip labels collapsed through ``hosts.dta_states``, and thresholds and
    outputs default from ``beast`` and the layout.
    """
    if not _resolve(
        args,
        derived={
            "trees": lambda b: [b.layout.beast_trait_trees(b.need("locus"), s) for s in b.seeds()],
            "out_json": lambda b: b.layout.trait_signal_report(b.need("locus")),
            "out_md": lambda b: b.layout.trait_signal_md(b.need("locus")),
            "burnin": lambda b: b.setting("beast", "burnin_fraction"),
            "permutations": lambda b: b.setting("beast", "discrete_trait", "signal_permutations"),
            "exit_nonzero_on_fail": lambda b: (
                (b.setting("beast", "discrete_trait", "on_no_signal") or "stop") == "stop"
            ),
        },
        required=("trees",),
        fallbacks={"burnin": 0.10, "permutations": 1000, "exit_nonzero_on_fail": False},
    ):
        return 2
    paths = [Path(p) for p in args.trees]
    missing = [p for p in paths if not p.is_file()]
    if missing:
        log.error("no such tree file: %s", missing[0])
        return 1

    trees = newick.read_posterior(
        paths, burnin=args.burnin, thin=args.thin, max_trees=args.max_trees
    )
    log.info(
        "read %d trees from %d file(s) (burnin %.0f%%, thin %d, cap %d)",
        len(trees),
        len(paths),
        args.burnin * 100,
        args.thin,
        args.max_trees,
    )

    labels = sorted(trees[0].tip_index)
    if args.traits:
        states = states_from_traits(args.traits)
        src = str(args.traits)
    else:
        states = states_from_labels(labels, args.label_field, args.label_sep)
        src = f"tip label field {args.label_field}"
    log.info("tip states from %s", src)

    # Collapse through hosts.dta_states so the tested states match the model's.
    # A host group absent from the map becomes '?', exactly as `beast xml` does
    # it, and those tips are pruned rather than made into a state of their own.
    if args.config:
        import yaml

        cfg = yaml.safe_load(args.config.read_text())
        dta = (cfg.get("hosts") or {}).get("dta_states") or {}
        if dta:
            before = sorted(set(states.values()))
            states = {k: dta.get(v, "?") for k, v in states.items()}
            after = sorted(set(states.values()) - {"?"})
            log.info(
                "collapsed via hosts.dta_states: %d group(s) -> %d state(s) %s",
                len(before),
                len(after),
                after,
            )
        else:
            log.info("config has no hosts.dta_states; using states as read")

    res = tsig.trait_signal(
        trees,
        states,
        n_permutations=args.permutations,
        seed=args.seed,
        ai_threshold=args.alpha,
        ps_threshold=args.alpha,
    )
    text = tsig.format_report(res)
    log.info("\n%s", text)

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(res.to_dict(), indent=2))
        log.info("wrote %s", args.out_json)
    if args.out_md:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text("```\n" + text + "\n```\n")
        log.info("wrote %s", args.out_md)

    if res.verdict == "fail" and args.exit_nonzero_on_fail:
        log.error("exiting non-zero: the trait carries no demonstrated signal.")
        return 1
    return 0


# ===========================================================================
# gates drt
# ===========================================================================
#
# Root-to-tip regression (gates temporal) is a screen and can be passed by a
# dataset with no clock signal at all. The date-randomisation test asks the
# question properly: shuffle the sampling dates among the tips, re-estimate the
# clock rate, and ask whether the real rate is distinguishable from rates
# estimated on randomised dates. If not, the rate is driven by tree shape, and
# every date in the analysis is unidentifiable.
#
# Two criteria are in use and they disagree, so both are reported.
#   STRICT (Duan et al. 2016; Murray et al. 2016): the real rate's 95% HPD must
#       not overlap the HPD of ANY replicate.
#   PERMISSIVE (Ramsden et al. 2008): the real point estimate must fall outside
#       the range of the replicate point estimates.
# The verdict uses the strict form.
#
# Pitfalls the summary exists to catch: replicates run with a shorter chain
# than the real analysis have wide posteriors that overlap everything, so poor
# mixing looks exactly like absent temporal signal (low-ESS replicates make the
# verdict `inconclusive`, never `fail`); and randomisation that did not
# randomise (a reused permutation or seed) makes the test pass trivially.


def rate_summary(path: Path, burnin: float, clock_param: str | None) -> dict[str, Any]:
    """Clock-rate posterior from one log: mean, HPD, ESS, retained samples."""
    header, data = read_beast_log(path)
    cut = round(burnin * len(data))
    kept = data[cut:]
    if len(kept) < 10:
        raise ValueError(
            f"{path}: {len(kept)} samples after a {burnin:.0%} burn-in; too short to summarise"
        )
    name = find_clock_parameter(header, clock_param)
    col = kept[:, header.index(name)]
    lo, hi = hpd(col)
    return {
        "path": str(path),
        "param": name,
        "mean": float(col.mean()),
        "median": float(np.median(col)),
        "hpd_low": lo,
        "hpd_high": hi,
        "ess": round(ess(col), 1),
        "n_kept": len(kept),
        "n_total": len(data),
    }


def intervals_overlap(a: dict[str, Any], b: dict[str, Any]) -> bool:
    """True when the HPD intervals of two rate summaries overlap."""
    return bool(a["hpd_low"] <= b["hpd_high"] and b["hpd_low"] <= a["hpd_high"])


def assess_drt(
    real: dict[str, Any], reps: list[dict[str, Any]], min_ess: int = 100
) -> dict[str, Any]:
    """Compare the real clock rate against the randomised replicates."""
    overlapping = [r for r in reps if intervals_overlap(real, r)]
    rep_means = [r["mean"] for r in reps]

    # The permissive criterion compares the real estimate against the ENVELOPE
    # of the replicate posteriors, not against the range of their point
    # estimates. The range of point estimates shrinks toward the standard
    # error as replicates are added -- with 20 replicates it is far narrower
    # than any single HPD -- so testing against it makes almost any real rate
    # look "outside", and a genuine failure gets rescued by the weaker test.
    lo_m = min(r["hpd_low"] for r in reps)
    hi_m = max(r["hpd_high"] for r in reps)
    outside_range = not (lo_m <= real["mean"] <= hi_m)

    low_ess = [r for r in reps if r["ess"] < min_ess]
    # A replicate posterior wider than the real one is the signature of poor
    # mixing rather than of a flat likelihood, and it is what makes an
    # unconverged replicate overlap everything.
    real_width = real["hpd_high"] - real["hpd_low"]
    wide = [r for r in reps if (r["hpd_high"] - r["hpd_low"]) > 3 * real_width]

    reasons: list[str] = []
    if not overlapping:
        verdict = "pass"
        reasons.append(
            f"the real rate's HPD overlaps none of the {len(reps)} replicates: "
            "the clock signal is not reproducible from randomised dates, which "
            "is what temporal signal means."
        )
    elif outside_range:
        verdict = "weak"
        reasons.append(
            f"{len(overlapping)} of {len(reps)} replicate HPDs overlap the "
            "real one, but the real point estimate falls outside the range of "
            "replicate estimates. This passes the permissive criterion and "
            "fails the strict one; dates estimated here carry wider real "
            "uncertainty than their HPDs show."
        )
    else:
        verdict = "fail"
        reasons.append(
            f"the real rate is inside the spread of randomised ones "
            f"({len(overlapping)} of {len(reps)} HPDs overlap). The rate is "
            "being driven by tree shape rather than by sampling dates, so the "
            "TMRCA and every node age derived from it are not identifiable."
        )

    # An unconverged replicate overlaps everything, which turns poor mixing
    # into an apparent absence of temporal signal. That must not read as fail.
    if low_ess and verdict in ("fail", "weak"):
        verdict = "inconclusive"
        reasons.append(
            f"{len(low_ess)} replicate(s) have clock-rate ESS below {min_ess}. "
            "A replicate that has not mixed has a wide posterior that overlaps "
            "everything, so this result cannot be distinguished from poor "
            "sampling. Lengthen the replicate chains "
            "(beast.date_randomisation.chain_length) and rerun before reading "
            "this as absent temporal signal."
        )
    elif low_ess:
        reasons.append(
            f"{len(low_ess)} replicate(s) have ESS below {min_ess}. The "
            "verdict stands — poor mixing widens replicate HPDs, which makes "
            "a pass harder rather than easier — but the replicates are "
            "under-sampled."
        )
    if wide:
        reasons.append(
            f"{len(wide)} replicate(s) have an HPD more than 3x the width of "
            "the real one, which usually means the chain never settled rather "
            "than that the randomised likelihood is flat."
        )

    # The reverse failure: randomisation that did not randomise -- the same
    # permutation reused, or the same seed, so the replicates are one replicate
    # counted many times and the test passes trivially.
    #
    # The test cannot be a fixed threshold on the spread of replicate means.
    # Independent replicates ALSO have tightly clustered means: the spread of a
    # mean goes as SD/sqrt(ESS), so a fixed cutoff measures sample size rather
    # than independence and fires on exactly the well-sampled runs it should
    # ignore.
    #
    # What separates the two is that genuinely different permutations give
    # different point estimates, so the between-replicate spread should EXCEED
    # the Monte Carlo error of a single replicate's mean. When it does not, the
    # replicates are statistically one draw.
    spread = float(np.std(rep_means)) if len(rep_means) > 1 else 0.0
    mc_errors = [((r["hpd_high"] - r["hpd_low"]) / 3.92) / np.sqrt(max(r["ess"], 1)) for r in reps]
    expected = float(np.median(mc_errors)) if mc_errors else 0.0
    independent = spread > expected if expected > 0 else True
    if len(reps) > 1 and not independent:
        reasons.append(
            f"the replicate estimates vary by less between replicates "
            f"({spread:.2e}) than the Monte Carlo error within one of them "
            f"({expected:.2e}). Different date permutations give different "
            "point estimates, so this is the signature of replicates sharing a "
            "permutation or a seed rather than of a tight null. A test whose "
            "replicates are one replicate passes trivially."
        )
        if verdict == "pass":
            verdict = "inconclusive"

    return {
        "verdict": verdict,
        "reasons": reasons,
        "n_replicates": len(reps),
        "n_overlapping": len(overlapping),
        "overlapping": [r["path"] for r in overlapping],
        "real_outside_replicate_range": bool(outside_range),
        "replicate_envelope": [lo_m, hi_m],
        "replicate_mean_spread": spread,
        "replicate_mc_error": expected,
        "low_ess_replicates": [r["path"] for r in low_ess],
        "min_ess": min_ess,
    }


def plot_drt(real: dict[str, Any], reps: list[dict[str, Any]], out: Path) -> bool:
    """Plot the real HPD against the replicate HPDs; False if matplotlib is absent."""
    try:
        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    fig, ax = plt.subplots(figsize=(7, 0.3 * len(reps) + 2.2))
    for i, r in enumerate(reps):
        ax.plot([r["hpd_low"], r["hpd_high"]], [i, i], lw=2, color="0.55", solid_capstyle="butt")
        ax.plot([r["mean"]], [i], "o", ms=3, color="0.3")
    y = len(reps) + 0.6
    ax.plot(
        [real["hpd_low"], real["hpd_high"]], [y, y], lw=3, color="#b2182b", solid_capstyle="butt"
    )
    ax.plot([real["mean"]], [y], "o", ms=5, color="#b2182b")
    ax.axvspan(real["hpd_low"], real["hpd_high"], color="#b2182b", alpha=0.08)
    ax.set_yticks([y, *list(range(len(reps)))])
    ax.set_yticklabels(["real"] + [f"rep {i + 1}" for i in range(len(reps))], fontsize=7)
    ax.set_xlabel(f"{real['param']}  (subs/site/year)")
    ax.set_title("Date-randomisation test: real clock rate vs randomised dates", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return True


def drt_report(real: dict[str, Any], reps: list[dict[str, Any]], res: dict[str, Any]) -> str:
    """Fixed-width text report of the date-randomisation test."""
    width = 78
    lines = [
        "=" * width,
        "DATE-RANDOMISATION TEST",
        "=" * width,
        f"clock parameter : {real['param']}",
        f"replicates      : {len(reps)}",
        "",
        f"  {'':<14}{'mean':>12}{'95% HPD':>24}{'ESS':>8}",
        "  " + "-" * (width - 4),
        (
            f"  {'REAL':<14}{real['mean']:>12.3e}"
            f"{real['hpd_low']:>12.3e}{real['hpd_high']:>12.3e}{real['ess']:>8.0f}"
        ),
        "",
    ]
    for i, r in enumerate(reps, 1):
        flag = ""
        if intervals_overlap(real, r):
            flag = "  <- overlaps real"
        if r["ess"] < res["min_ess"]:
            flag += "  <- low ESS"
        lines.append(
            f"  {'rep ' + str(i):<14}{r['mean']:>12.3e}"
            f"{r['hpd_low']:>12.3e}{r['hpd_high']:>12.3e}"
            f"{r['ess']:>8.0f}{flag}"
        )
    lines += [
        "",
        (
            f"  strict     : {res['n_overlapping']} of {res['n_replicates']} "
            f"replicate HPDs overlap the real HPD  "
            f"({'PASS' if res['n_overlapping'] == 0 else 'FAIL'})"
        ),
        (
            f"  permissive : the real estimate is "
            f"{'OUTSIDE' if res['real_outside_replicate_range'] else 'INSIDE'} "
            f"the replicate envelope "
            f"[{res['replicate_envelope'][0]:.3e}, "
            f"{res['replicate_envelope'][1]:.3e}]"
        ),
        "",
        f"VERDICT: {res['verdict'].upper()}",
        "",
    ]
    lines.extend("  - " + reason for reason in res["reasons"])
    lines += [
        "",
        "=" * width,
        "A pass means the timescale is identifiable from the data, not that",
        "it is correct. It says nothing about whether the clock model, the",
        "tree prior or the subsample are appropriate.",
        "=" * width,
    ]
    return "\n".join(lines)


def _drt_logs(b: _Build) -> list[Path]:
    """The date-randomised replicates' logs, replicates numbered from 1 as the Snakefile does."""
    n_reps = b.setting("beast", "date_randomisation", "replicates")
    if not n_reps:
        msg = "beast.date_randomisation.replicates is not set in the config"
        raise _UsageError(msg)
    return [b.layout.drt_log(b.need("locus"), str(rep)) for rep in range(1, int(n_reps) + 1)]


def add_arguments_drt(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo gates drt``."""
    _add_build_arguments(parser)
    parser.add_argument(
        "--real",
        type=Path,
        help="Combined post-burn-in log from the real analysis (what `gates convergence` "
        "writes). Required unless --build-dir.",
    )
    parser.add_argument(
        "--reps",
        nargs="+",
        type=Path,
        help="One log per date-randomised replicate. Required unless --build-dir, which "
        "takes replicates 1..beast.date_randomisation.replicates.",
    )
    parser.add_argument(
        "--burnin",
        type=float,
        help="Applied to the replicates. The combined real log is "
        "already post-burn-in, so 0 is used for it unless "
        "--real-burnin says otherwise. Default 0.10, or beast.burnin_fraction "
        "with --build-dir.",
    )
    parser.add_argument("--real-burnin", type=float, default=0.0)
    parser.add_argument(
        "--clock-param",
        default=None,
        help="Column holding the clock rate. Auto-detected from "
        "the usual BEAST names when not given.",
    )
    parser.add_argument(
        "--min-ess",
        type=int,
        default=100,
        help="Replicate ESS below this makes the verdict inconclusive rather than a failure.",
    )
    parser.add_argument("--out-json", type=Path)
    parser.add_argument("--out-plot", type=Path)
    parser.add_argument("--out-md", type=Path)
    parser.add_argument(
        "--exit-nonzero-on-fail",
        action=argparse.BooleanOptionalAction,
        help="Off by default; with --build-dir, on unless "
        "beast.date_randomisation.on_fail is warn.",
    )


@step("gates", "drt", "Date-randomisation test, summarised and judged")
def drt(args: argparse.Namespace) -> int:
    """Summarise and judge the date-randomisation test.

    Root-to-tip regression (``gates temporal``) asks whether genetic distance
    increases with sampling date. It is a screen, and it can be passed by a dataset
    with no clock signal at all: with enough tips, a trivial correlation is
    significant, and population structure alone produces one when older samples
    happen to belong to a diverged lineage. The CDV H tree passed by permutation at
    R^2 = 0.024.

    The date-randomisation test asks the question properly. Shuffle the sampling
    dates among the tips, so the tree and the sequences are untouched and only the
    date-to-tip assignment is destroyed, then re-estimate the clock rate. If the
    real rate is indistinguishable from rates estimated on randomised dates, the
    rate is being driven by tree shape rather than by time, and every date in the
    analysis (TMRCA, node ages, the whole timescale) is unidentifiable.

    Two criteria are in use and they disagree, so both are reported:

      STRICT (Duan et al. 2016; Murray et al. 2016)
          the real rate's 95% HPD must not overlap the HPD of ANY replicate.

      PERMISSIVE (Ramsden et al. 2008, as originally applied)
          the real point estimate must fall outside the range of the replicate
          point estimates.

    The verdict uses the strict form, because the permissive one is passed routinely
    by datasets that go on to produce unstable TMRCAs: comparing a point estimate
    against point estimates ignores that the replicate posteriors are often very
    wide, and a real rate sitting just outside a cloud of overlapping intervals is
    not distinguishable from them.

    The verdict is ``pass``, ``weak``, ``fail`` or ``inconclusive``; the exit code
    is 1 on ``fail`` only with ``--exit-nonzero-on-fail``, and 1 on unreadable
    input.

    With ``--build-dir`` and ``--config`` (and ``--locus``) the real log is the
    locus's combined log, the replicates are its DRT logs, and the burn-in,
    outputs and fail policy default from ``beast`` and the layout.
    """
    if not _resolve(
        args,
        derived={
            "real": lambda b: b.layout.combined_log(b.need("locus")),
            "reps": _drt_logs,
            "burnin": lambda b: b.setting("beast", "burnin_fraction"),
            "out_json": lambda b: b.layout.drt_summary(b.need("locus")),
            "out_plot": lambda b: b.layout.drt_plot(b.need("locus")),
            "exit_nonzero_on_fail": lambda b: (
                (b.setting("beast", "date_randomisation", "on_fail") or "stop") == "stop"
            ),
        },
        required=("real", "reps"),
        fallbacks={"burnin": 0.10, "exit_nonzero_on_fail": False},
    ):
        return 2
    missing = [p for p in [args.real, *args.reps] if not p.is_file()]
    if missing:
        log.error("no such log file: %s", missing[0])
        return 1
    if len(args.reps) < 5:
        log.warning(
            "only %d replicates. Ten is a usual minimum; with fewer, the "
            "replicate range is too poorly estimated for the permissive "
            "criterion to mean much.",
            len(args.reps),
        )

    try:
        real = rate_summary(args.real, args.real_burnin, args.clock_param)
        reps = [rate_summary(p, args.burnin, args.clock_param or real["param"]) for p in args.reps]
    except ValueError as e:
        log.error("%s", e)  # noqa: TRY400 -- message only, no traceback
        return 1

    res = assess_drt(real, reps, min_ess=args.min_ess)
    text = drt_report(real, reps, res)
    log.info("\n%s", text)

    if args.out_plot:
        if plot_drt(real, reps, args.out_plot):
            log.info("wrote %s", args.out_plot)
        else:
            log.warning("matplotlib not available; skipped the plot")
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps({**res, "real": real, "replicates": reps}, indent=2))
        log.info("wrote %s", args.out_json)
    if args.out_md:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text("```\n" + text + "\n```\n")
        log.info("wrote %s", args.out_md)

    if res["verdict"] == "fail" and args.exit_nonzero_on_fail:
        log.error("exiting non-zero: the timescale is not identifiable.")
        return 1
    return 0


# ===========================================================================
# gates compare
# ===========================================================================
#
# A scheduled rebuild produces a new tree every month, and a new tree always
# differs from the old one. A large change between builds is a CURATION
# REGRESSION until shown otherwise: new data arriving monthly moves a TMRCA by
# months, not decades, and a dataset that shrinks has almost always lost
# records to a changed query rather than to GenBank withdrawals.
#
# Compared, whatever both builds have: accession lists (count and WHICH
# accessions came and went), the temporal signal verdict, R^2 and slope, every
# gate verdict, and the TMRCA from the Auspice JSON. Tree topology is
# deliberately NOT compared: it is large, noisy and uninformative about whether
# anything is wrong.


@dataclass
class Finding:
    """One line of the comparison report."""

    level: str  # ok | note | warn | error
    area: str
    message: str


@dataclass
class Comparison:
    """All findings from comparing two builds, and the verdict they imply."""

    pathogen: str
    previous: str | None = None
    current: str | None = None
    findings: list[Finding] = field(default_factory=list)

    def add(self, level: str, area: str, message: str) -> None:
        """Append a finding."""
        self.findings.append(Finding(level, area, message))

    @property
    def verdict(self) -> str:
        """``regression`` if any error, ``review`` if any warning, else ``ok``."""
        if any(f.level == "error" for f in self.findings):
            return "regression"
        if any(f.level == "warn" for f in self.findings):
            return "review"
        return "ok"

    def to_dict(self) -> dict[str, Any]:
        """JSON-ready form including the verdict."""
        d = asdict(self)
        d["verdict"] = self.verdict
        return d


def dated_builds(root: Path) -> list[Path]:
    """Dated build directories, oldest first. Names are YYYYMMDD."""
    if not root.is_dir():
        return []
    return sorted(
        (d for d in root.iterdir() if d.is_dir() and d.name.isdigit()), key=lambda d: d.name
    )


def read_accessions(build: Path) -> set[str] | None:
    """The accession snapshot of a build, or None if it has none."""
    hits = sorted(build.glob("raw/*.acc"))
    if not hits:
        return None
    return {ln.strip() for ln in hits[0].read_text().splitlines() if ln.strip()}


def read_json(path: Path) -> Any:
    """Parsed JSON, or None if the file is missing or malformed."""
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def read_temporal(build: Path) -> dict[str, Any] | None:
    """The first temporal signal report found in a build, or None."""
    for p in sorted(build.glob("temporal/*_temporal_signal.json")):
        doc = read_json(p)
        if doc and isinstance(doc, dict):
            return doc
    return None


def gate_verdicts(build: Path) -> dict[str, str]:
    """Every gate verdict in a build, keyed by a readable name."""
    out: dict[str, str] = {}
    patterns = {
        "temporal": "temporal/*_temporal_signal.json",
        "convergence": "beast/*/convergence.json",
        "trait_signal": "clades/*/trait_signal.json",
        "drt": "drt/*/drt_summary.json",
    }
    for name, pat in patterns.items():
        for p in sorted(build.glob(pat)):
            doc = read_json(p)
            if isinstance(doc, dict) and "verdict" in doc:
                # keep the directory so per-clade gates stay distinguishable
                key = (
                    f"{name}:{p.parent.name}"
                    if "*" in pat.split("/")[0] or len(pat.split("/")) > 2
                    else name
                )
                out[key] = doc["verdict"]
    return out


def auspice_tmrca(doc: Any) -> float | None:
    """Root date from an Auspice v2 JSON.

    The root's num_date is the TMRCA. Augur writes it under node_attrs; a
    missing value means the tree was not time-scaled, which is a different
    thing from a TMRCA of zero and must not be compared as a number.
    """
    if not isinstance(doc, dict):
        return None
    tree = doc.get("tree")
    if isinstance(tree, list):
        tree = tree[0] if tree else None
    if not isinstance(tree, dict):
        return None
    nd = (tree.get("node_attrs") or {}).get("num_date")
    if isinstance(nd, dict) and "value" in nd:
        try:
            return float(nd["value"])
        except (TypeError, ValueError):
            return None
    return None


VERDICT_RANK = {
    "pass": 3,
    "ok": 3,
    "weak": 2,
    "review": 2,
    "inconclusive": 1,
    "fail": 0,
    "regression": 0,
    "unknown": 1,
}


def compare_builds(
    prev: Path,
    curr: Path,
    pathogen: str,
    max_tmrca_shift: float = 5.0,
    alert_on_loss: bool = True,
    current_auspice: Any = None,
    previous_auspice: Any = None,
) -> Comparison:
    """Compare two dated builds: accessions, temporal signal, gate verdicts, TMRCA."""
    c = Comparison(pathogen=pathogen, previous=prev.name, current=curr.name)

    # --- accessions ---------------------------------------------------------
    a_prev, a_curr = read_accessions(prev), read_accessions(curr)
    if a_prev is None or a_curr is None:
        c.add(
            "note",
            "accessions",
            "one of the builds has no accession list, so dataset "
            "composition could not be compared",
        )
    else:
        gained = a_curr - a_prev
        lost = a_prev - a_curr
        c.add(
            "note",
            "accessions",
            f"{len(a_prev)} -> {len(a_curr)} ({len(a_curr) - len(a_prev):+d}); "
            f"{len(gained)} added, {len(lost)} removed",
        )
        if lost:
            shown = ", ".join(sorted(lost)[:10])
            more = f" and {len(lost) - 10} more" if len(lost) > 10 else ""
            level = "error" if alert_on_loss else "warn"
            c.add(
                level,
                "accessions",
                f"{len(lost)} accession(s) present in {prev.name} are absent "
                f"from {curr.name}: {shown}{more}. GenBank withdraws records "
                "rarely; a changed query or a stricter curation rule is the "
                "usual cause, and either is a regression rather than a "
                "finding",
            )
        # A stable count hides a composition change completely, and a
        # composition change is exactly what a broken query produces.
        if not lost and not gained and a_prev == a_curr:
            c.add("note", "accessions", "the dataset is unchanged")
        elif lost and gained and len(a_prev) == len(a_curr):
            c.add(
                "warn",
                "accessions",
                f"the count is unchanged but {len(lost)} record(s) were "
                "swapped out. A count-only check would have reported nothing",
            )

    # --- temporal signal ----------------------------------------------------
    t_prev, t_curr = read_temporal(prev), read_temporal(curr)
    if t_prev and t_curr:
        for key, label, tol in (("r_squared", "R^2", 0.5), ("slope", "clock slope", 0.5)):
            pv, cv = t_prev.get(key), t_curr.get(key)
            if isinstance(pv, (int, float)) and isinstance(cv, (int, float)) and pv:
                # Signed, because the direction is the whole point: R^2
                # falling is a regression and R^2 rising is not, and an
                # absolute magnitude printed with a + reads as the opposite
                # of what happened.
                change = (cv - pv) / abs(pv)
                msg = f"{label} {pv:.4g} -> {cv:.4g} ({change:+.0%})"
                c.add("warn" if abs(change) > tol else "note", "temporal", msg)

    # --- every gate verdict -------------------------------------------------
    g_prev, g_curr = gate_verdicts(prev), gate_verdicts(curr)
    for name in sorted(set(g_prev) | set(g_curr)):
        pv, cv = g_prev.get(name), g_curr.get(name)
        if pv is None or cv is None:
            c.add(
                "note",
                "gates",
                f"{name}: {'absent' if pv is None else pv} -> {'absent' if cv is None else cv}",
            )
            continue
        if pv == cv:
            continue
        before, after = VERDICT_RANK.get(pv, 1), VERDICT_RANK.get(cv, 1)
        level = "error" if after < before else "note"
        direction = "regressed" if after < before else "improved"
        c.add(level, "gates", f"{name} {direction}: {pv} -> {cv}")

    # --- TMRCA --------------------------------------------------------------
    tm_prev = auspice_tmrca(previous_auspice)
    tm_curr = auspice_tmrca(current_auspice)
    if tm_prev is not None and tm_curr is not None:
        shift = abs(tm_curr - tm_prev)
        msg = f"TMRCA {tm_prev:.1f} -> {tm_curr:.1f} ({tm_curr - tm_prev:+.1f} yr)"
        if shift > max_tmrca_shift:
            c.add(
                "error",
                "tmrca",
                msg + f"; exceeds alert_on_tmrca_shift_years "
                f"({max_tmrca_shift}). A month of new sequences moves a "
                "TMRCA by months. A shift this size is a changed "
                "dataset, prior or clock model until shown otherwise",
            )
        else:
            c.add("note", "tmrca", msg)
    elif tm_curr is not None or tm_prev is not None:
        c.add(
            "note",
            "tmrca",
            "only one build has a time-scaled tree, so the TMRCA could not be compared",
        )

    return c


def render_comparison(c: Comparison) -> str:
    """Markdown report of a comparison."""
    icon = {"ok": "ok", "note": "  ", "warn": "! ", "error": "X "}
    lines = [
        f"# Build comparison — {c.pathogen}",
        "",
        f"`{c.previous}` → `{c.current}`",
        "",
        f"**Verdict: {c.verdict.upper()}**",
        "",
    ]
    if c.verdict == "regression":
        lines += [
            (
                "A change this size is a curation regression until shown "
                "otherwise. Check the query and the curation rules before "
                "treating it as a result."
            ),
            "",
        ]
    for area in ("accessions", "temporal", "gates", "tmrca"):
        rows = [f for f in c.findings if f.area == area]
        if not rows:
            continue
        lines.append(f"## {area}")
        lines.extend(f"- {icon.get(f.level, '  ')} {f.message}" for f in rows)
        lines.append("")
    lines += [
        "---",
        "",
        (
            "Tree topology is deliberately not compared: two trees built from "
            "overlapping data differ mostly by reordering equally-supported "
            "splits, so a topological distance is large, noisy and says nothing "
            "about whether anything is wrong."
        ),
    ]
    return "\n".join(lines)


def _current_dataset(b: _Build) -> Path | None:
    """The build's Auspice JSON, if it exists: one locus has the bare dataset name."""
    build_name = b.setting("nextstrain", "build_name")
    if not build_name:
        msg = "nextstrain.build_name is not set in the config"
        raise _UsageError(msg)
    loci = b.loci()
    path = b.layout.auspice_dataset(str(build_name), None if len(loci) == 1 else b.need("locus"))
    if not path.is_file():
        log.info("no Auspice dataset at %s; comparing without a TMRCA", path)
        return None
    return path


def add_arguments_compare(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo gates compare``."""
    _add_build_arguments(parser)
    parser.add_argument(
        "--builds-dir",
        type=Path,
        help="builds/<pathogen>, holding dated build directories. Required unless "
        "--build-dir, whose parent it is.",
    )
    parser.add_argument(
        "--pathogen",
        help="Pathogen id, for the report header. Defaults to the builds-dir name.",
    )
    parser.add_argument(
        "--current",
        type=Path,
        help="Auspice JSON from the current build, when there is "
        "one. CI stops at the temporal gate and has none. With --build-dir it "
        "defaults to the build's dataset when that exists (--locus picks one "
        "when several loci are analysed).",
    )
    parser.add_argument(
        "--previous",
        type=Path,
        help="Auspice JSON from the previous build, for the TMRCA "
        "comparison. Found in the previous build directory "
        "when not given.",
    )
    parser.add_argument(
        "--max-tmrca-shift",
        type=float,
        help="Default 5.0, or nextstrain.update.alert_on_tmrca_shift_years with --build-dir.",
    )
    parser.add_argument(
        "--alert-on-sequence-loss",
        action=argparse.BooleanOptionalAction,
        help="Default on, or nextstrain.update.alert_on_sequence_loss with --build-dir.",
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument("--out-json", type=Path)
    parser.add_argument(
        "--exit-zero",
        action="store_true",
        help="Always exit 0. For looking at a comparison without failing a workflow.",
    )


@step("gates", "compare", "What changed since the last build, and is it a finding")
def compare(args: argparse.Namespace) -> int:
    """Compare the newest dated build against the one before it: what changed, and is it a finding?

    A scheduled rebuild produces a new tree every month, and a new tree always
    differs from the old one. Without a comparison step the only signal available is
    whether the pipeline exited zero, and a pipeline exits zero just as happily when
    a curation rule change has silently dropped two hundred sequences or when the
    TMRCA has moved forty years.

    The framing that matters: a large change between builds is a CURATION REGRESSION
    until shown otherwise. New data arriving monthly moves a TMRCA by months, not
    decades. A dataset that shrinks has almost always lost records to a changed
    query rather than to GenBank withdrawals. Treating these as discoveries is how a
    pipeline publishes its own bugs.

    What is compared is whatever both builds have, because CI and a full local run
    leave different artefacts behind:

      accession lists   count, and WHICH accessions came and went. A count can be
                        unchanged while the composition turns over completely, and
                        that is invisible to a count-only check.
      temporal signal   the gate verdict, R^2 and the clock slope. A verdict that
                        drops from pass to fail between builds is the single most
                        informative regression available, because everything
                        downstream is conditioned on it.
      gate verdicts     convergence, trait signal, date randomisation, per clade.
      TMRCA             from the Auspice JSON when there is one, against
                        alert_on_tmrca_shift_years.

    Deliberately NOT compared: tree topology. Two trees built from overlapping data
    differ in ways that are mostly reordering of equally-supported splits, so a
    topological distance between builds is large, noisy, and uninformative about
    whether anything is wrong. The quantities above are the ones where a change has
    an interpretation.

    Exit code is non-zero when something regressed, so CI opens an issue rather than
    publishing quietly. A first build has nothing to compare against and exits 0;
    ``--exit-zero`` forces 0.

    With ``--build-dir`` and ``--config`` the builds directory is the build's
    parent, the pathogen its id, the current dataset the build's Auspice JSON,
    the thresholds come from ``nextstrain.update`` and the reports go in the build.
    """
    if not _resolve(
        args,
        derived={
            "builds_dir": lambda b: b.layout.pathogen_dir(),
            "pathogen": lambda b: b.layout.pathogen,
            "current": _current_dataset,
            "max_tmrca_shift": lambda b: b.setting(
                "nextstrain", "update", "alert_on_tmrca_shift_years"
            ),
            "alert_on_sequence_loss": lambda b: b.setting(
                "nextstrain", "update", "alert_on_sequence_loss"
            ),
            "out": lambda b: b.layout.build_comparison(),
            "out_json": lambda b: b.layout.build_comparison_json(),
        },
        required=("builds_dir",),
        fallbacks={"max_tmrca_shift": 5.0, "alert_on_sequence_loss": True},
    ):
        return 2
    pathogen = args.pathogen or args.builds_dir.name
    builds = dated_builds(args.builds_dir)

    # A first build has nothing to compare against. That is the normal state of
    # a new pathogen, not a failure, and must not block the first publication.
    if len(builds) < 2:
        text = (
            f"# Build comparison — {pathogen}\n\n"
            f"Only {len(builds)} build in `{args.builds_dir}`; nothing to "
            f"compare against yet.\n"
        )
        log.info("\n%s", text)
        if args.out:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(text)
        return 0

    prev, curr = builds[-2], builds[-1]

    cur_doc = read_json(args.current) if args.current else None
    if args.current and cur_doc is None:
        log.error("could not read %s as JSON", args.current)
        return 1
    prev_doc = read_json(args.previous) if args.previous else None
    if prev_doc is None:
        for p in sorted(prev.glob("auspice/*.json")):
            prev_doc = read_json(p)
            if prev_doc:
                break

    c = compare_builds(
        prev,
        curr,
        pathogen,
        max_tmrca_shift=args.max_tmrca_shift,
        alert_on_loss=args.alert_on_sequence_loss,
        current_auspice=cur_doc,
        previous_auspice=prev_doc,
    )

    text = render_comparison(c)
    log.info("\n%s", text)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n")
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(c.to_dict(), indent=2))

    if c.verdict == "regression" and not args.exit_zero:
        log.error("exiting non-zero: something regressed between builds.")
        return 1
    return 0

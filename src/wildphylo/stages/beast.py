"""beast stage: subsample, BEAST XML generation, and posterior trait-tree summaries.

Steps (see ``docs/PLAN-stages.md``):

* ``subsample``: stratified subsample of an alignment (or one clade) sized for BEAST.
* ``xml``: BEAST 2 XML from a labelled alignment and the pathogen config.
* ``jumps``: when each host-state transition happened in a posterior of BEAST trees.
* ``state-through-time``: lineages carrying each host state, year by year.

Tip labels are ``ACCESSION|host_group|decimal_year`` throughout.
"""

from __future__ import annotations

import glob as globmod
import json
import logging
import random
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any

from wildphylo.beastxml import (
    clock_blocks,
    resolve_clock_rate,
    tip_date_blocks,
    trait_blocks,
    tree_prior_blocks,
)
from wildphylo.cli import step
from wildphylo.config import ConfigError
from wildphylo.config import load as load_config
from wildphylo.newick import (
    Node,
    count_trees,
    iter_trees,
    label_date,
    node_depths,
)
from wildphylo.preflight import check_clock_prior
from wildphylo.stepargs import UsageError, add_build_dir_arguments, fail, resolve

if TYPE_CHECKING:
    import argparse
    from collections.abc import Callable, Iterable, Sequence

    import numpy as np
    import numpy.typing as npt
    import pandas as pd

    from wildphylo.stepargs import Build

log = logging.getLogger(__name__)

DEFAULT_SEED = 20260819
GAP = frozenset("-.")


class BeastStepError(Exception):
    """A beast step cannot continue with the data it was given; exit code 1.

    Usage errors (a missing or inconsistent flag) are not this: they raise
    :class:`~wildphylo.stepargs.UsageError`, which :func:`~wildphylo.stepargs.fail`
    logs with exit code 2.
    """


def _say(text: str = "") -> None:
    """Write one line of report output to stdout (the steps' report is their stdout)."""
    sys.stdout.write(text + "\n")


def _err(text: str) -> None:
    """Write one line to stderr."""
    sys.stderr.write(text + "\n")


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA into {first-word-of-header: sequence}."""
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def write_fasta(path: Path, seqs: dict[str, str], width: int = 60) -> None:
    """Write sequences as FASTA, wrapped at ``width`` columns."""
    with path.open("w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


# =========================================================================== #
# subsample (was 06_subsample.py)
#
# Two modes:
#   global  stratified subsample across the whole dataset, preserving clade,
#           host-group and temporal structure
#   clade   all tips of one clade, optionally subsampled
#
# Why stratify: a random subsample of a dataset that is 52% domestic dog will be
# 52% domestic dog, discarding most of the rare host groups a host-jump analysis
# depends on. Stratifying by clade x host group x year bin keeps the informative
# tips and drops the redundant ones.
#
# Sampling within a stratum is deliberate, not random:
#   1. sequences with day-precision dates before year-only dates
#   2. longer ungapped sequences before shorter ones
#   3. one representative per identical-sequence set
# =========================================================================== #

# Fallback wild groups when no config is given. The pathogen config's
# hosts.wild_groups overrides this: the set below is CDV's carnivore range and
# means nothing for a pathogen of ruminants or rodents.
DEFAULT_WILD_GROUPS = frozenset(
    {
        "wild_felid",
        "wild_canid",
        "mustelid",
        "procyonid",
        "pinniped",
        "ursid",
        "ailurid",
        "viverrid",
    }
)

PRECISION_RANK = {"day": 0, "month": 1, "range": 2, "year": 3, "unknown": 4, "none": 5}


def strip_all_gap_columns(seqs: dict[str, str]) -> tuple[dict[str, str], int]:
    """Drop columns that are gaps in every remaining sequence.

    A subset of an alignment keeps the full column count of the original, so any
    column that only ever held bases belonging to dropped sequences survives as
    an all-gap column. These are not harmless. Anything that treats an alignment
    column span as a coordinate (a CDS annotation, a codon partition, a site
    index quoted in a figure) is thrown off by them, and the offset is invisible
    because the columns look like ordinary gaps.
    """
    if not seqs:
        return seqs, 0
    n_cols = len(next(iter(seqs.values())))
    keep = [i for i in range(n_cols) if any(s[i] not in GAP for s in seqs.values())]
    if len(keep) == n_cols:
        return seqs, 0
    return {n: "".join(s[i] for i in keep) for n, s in seqs.items()}, n_cols - len(keep)


def parse_label(label: str) -> tuple[str, str, float | None]:
    """Split a tip label into (accession, host group, decimal year or None)."""
    parts = str(label).split("|")
    if len(parts) < 3:
        return (parts[0], "unparsed", None)
    try:
        year: float | None = float(parts[-1])
    except ValueError:
        year = None
    return (parts[0], parts[1], year)


def build_records(
    seqs: dict[str, str], clades_df: pd.DataFrame | None, meta_df: pd.DataFrame | None
) -> pd.DataFrame:
    """One row per aligned sequence, with clade, host group, year and quality fields."""
    import pandas as pd

    clade_of = (
        dict(zip(clades_df["accession"], clades_df["clade"], strict=True))
        if clades_df is not None
        else {}
    )
    prec_of = {}
    if meta_df is not None:
        base = meta_df["accession"].astype(str).str.split(".").str[0]
        prec_of = dict(
            zip(
                base,
                meta_df.get("date_precision", pd.Series(["year"] * len(meta_df))),
                strict=True,
            )
        )

    rows = []
    for label, seq in seqs.items():
        acc, host, year = parse_label(label)
        ungapped = "".join(c for c in seq if c not in GAP)
        rows.append(
            {
                "label": label,
                "accession": acc,
                "host_group": host,
                "decimal_year": year,
                "clade": clade_of.get(acc, "unassigned"),
                "date_precision": prec_of.get(acc, "unknown"),
                "ungapped_len": len(ungapped),
                "ambig": sum(1 for c in ungapped.upper() if c in "NX?"),
                "seq_key": ungapped.upper(),
            }
        )
    return pd.DataFrame(rows)


def rank_within_stratum(df: pd.DataFrame) -> pd.DataFrame:
    """Order a stratum best-first: precise date, few ambiguity codes, long sequence."""
    d = df.copy()
    d["prec_rank"] = d["date_precision"].map(lambda p: PRECISION_RANK.get(p, 4))
    d["ambig_frac"] = d["ambig"] / d["ungapped_len"].clip(lower=1)
    return d.sort_values(
        ["prec_rank", "ambig_frac", "ungapped_len"], ascending=[True, True, False]
    )


def drop_duplicates(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Keep one representative per identical sequence, preferring the best-dated."""
    ranked = rank_within_stratum(df)
    kept = ranked.drop_duplicates(subset="seq_key", keep="first")
    return kept, len(df) - len(kept)


def stratified_sample(
    df: pd.DataFrame,
    target: int,
    year_bin: int,
    wild_groups: frozenset[str],
    min_per_stratum: int = 1,
) -> pd.DataFrame:
    """Allocate ``target`` across clade x host group x year-bin strata.

    Every non-empty stratum gets at least ``min_per_stratum``, then the remaining
    budget is distributed proportional to the square root of stratum size, which
    keeps large strata dominant without letting them crowd out the rare host
    groups entirely. When there are more strata than budget, one tip is taken per
    stratum and wild-group strata come first.
    """
    import pandas as pd

    d = df.copy()
    d["year_bin"] = d["decimal_year"].apply(
        lambda y: int(y // year_bin) * year_bin if pd.notna(y) else -1
    )
    d["stratum"] = d["clade"] + "|" + d["host_group"] + "|" + d["year_bin"].astype(str)

    strata = dict(iter(d.groupby("stratum")))  # dict(groupby) would look up a "keys" column
    if len(strata) >= target:
        picks = []
        order = sorted(
            strata.items(),
            key=lambda kv: (kv[1]["host_group"].iloc[0] not in wild_groups, -len(kv[1])),
        )
        for _, g in order[:target]:
            picks.append(rank_within_stratum(g).head(1))
        return pd.concat(picks)

    allocation = dict.fromkeys(strata, min_per_stratum)
    remaining = target - sum(allocation.values())
    if remaining > 0:
        weights = {k: max(len(g) - min_per_stratum, 0) ** 0.5 for k, g in strata.items()}
        total_w = sum(weights.values())
        if total_w > 0:
            for k, w in weights.items():
                allocation[k] += int(remaining * w / total_w)

    # hand out any rounding leftovers to the largest under-sampled strata
    used = sum(min(allocation[k], len(g)) for k, g in strata.items())
    leftover = target - used
    if leftover > 0:
        for k, g in sorted(strata.items(), key=lambda kv: -len(kv[1])):
            room = len(g) - allocation[k]
            if room > 0:
                take = min(room, leftover)
                allocation[k] += take
                leftover -= take
            if leftover <= 0:
                break

    picks = [rank_within_stratum(g).head(min(allocation[k], len(g))) for k, g in strata.items()]
    return pd.concat(picks)


def summarise(df: pd.DataFrame, title: str, wild_groups: frozenset[str]) -> None:
    """Report host-group, year and clade composition of a set of sequences."""
    _say(f"\n{title}: {len(df)} sequences")
    if df.empty:
        return
    hosts = Counter(df["host_group"])
    n_wild = sum(n for h, n in hosts.items() if h in wild_groups)
    for host, n in hosts.most_common():
        _say(f"    {host:<16} {n:>5}{'  <- wild' if host in wild_groups else ''}")
    _say(f"    {'wild total':<16} {n_wild:>5}  ({n_wild / len(df):.1%})")
    yrs = df["decimal_year"].dropna()
    if len(yrs):
        _say(f"    years {yrs.min():.1f} - {yrs.max():.1f}")
    clades = Counter(df["clade"])
    if len(clades) > 1:
        _say("    clades: " + ", ".join(f"{c}={n}" for c, n in clades.most_common(8)))


def add_arguments_subsample(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo beast subsample``."""
    parser.add_argument("--mode", choices=["global", "clade"], required=True)
    parser.add_argument("--aln", type=Path, help="Alignment (default with --build-dir: aligned/)")
    parser.add_argument("--clades", type=Path, help="clade table from `wildphylo tree clades`")
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--clade", help="Clade name, required for --mode clade")
    parser.add_argument(
        "--target",
        type=int,
        default=None,
        help="Target sequence count (global default 400; clade default = keep all)",
    )
    parser.add_argument("--year-bin", type=int, default=5, help="Time bin width in years")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--keep-duplicates", action="store_true")
    parser.add_argument(
        "--keep-all-gap-columns",
        action="store_true",
        help="Keep columns that are all-gap in the subset. Off by default: "
        "such columns shift every downstream column coordinate and are "
        "invisible when you look at the alignment.",
    )
    parser.add_argument("--outdir", type=Path, default=Path("data/processed"))
    parser.add_argument("--prefix", default=None)
    # Explicit outputs, so Snakemake can track what this rule produces.
    parser.add_argument(
        "--out-aln",
        type=Path,
        default=None,
        help="Subsampled alignment. Overrides --outdir/--prefix.",
    )
    parser.add_argument(
        "--out-meta",
        type=Path,
        default=None,
        help="Subsampled metadata. Overrides --outdir/--prefix.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Pathogen config: supplies hosts.wild_groups, "
        "subsample.target_tips, year_bin and seed. Required with --build-dir.",
    )
    add_build_dir_arguments(parser)


@step(
    "beast",
    "subsample",
    "Stratified subsample of an alignment (or one clade) sized for BEAST",
)
def subsample(args: argparse.Namespace) -> int:
    """Reduce an alignment to a BEAST-tractable size, or pull one clade.

    Two modes:

      global   stratified subsample across the whole dataset, preserving clade,
               host-group and temporal structure
      clade    extract all tips belonging to one clade, optionally subsampled

    Why stratify: a random subsample of a dataset that is 52% domestic dog will be
    52% domestic dog, discarding most of the rare host groups your host-jump
    analysis depends on. Stratifying by clade x host x time bin keeps the
    informative tips and drops the redundant ones.

    Sampling within a stratum is deliberate, not random:

      1. sequences with day-precision dates before year-only dates
      2. longer ungapped sequences before shorter ones
      3. one representative per identical-sequence set

    Writes the subset FASTA, a metadata table for the subset, and a dates table.
    The selection is deterministic given the inputs (ranking, not random draws);
    ``--seed`` is accepted and read from the config but does not change it.
    """
    import pandas as pd

    try:
        resolve(
            args,
            derived={
                "aln": lambda b: b.layout.alignment(b.need("locus")),
                "metadata": lambda b: b.layout.metadata(),
                "out_aln": lambda b: b.layout.subsample(b.need("locus")),
                "out_meta": lambda b: b.layout.subsample_metadata(b.need("locus")),
            },
            required=("aln",),
            outputs=("out_aln", "out_meta"),
        )
    except UsageError as exc:
        return fail(exc)
    if args.mode == "clade" and not args.clade:
        return fail(UsageError("--mode clade requires --clade (e.g. --clade clade_3)"))
    if not args.aln.is_file():
        _err(f"Alignment not found: {args.aln}")
        return 1

    seqs = read_fasta(args.aln)
    # Config first, command line wins. The wild-group set in particular MUST
    # come from config: the module default is CDV's carnivore range and is
    # simply wrong for a pathogen of ruminants or rodents, where it would
    # silently classify every wild host as domestic.
    wild_groups = DEFAULT_WILD_GROUPS
    if args.config:
        try:
            cfg = load_config(args.config)
            configured = (cfg.get("hosts") or {}).get("wild_groups")
            if configured:
                wild_groups = frozenset(configured)
            sub_cfg = cfg.get("subsample", {})
            if args.target is None and "target_tips" in sub_cfg:
                args.target = int(sub_cfg["target_tips"])
            if args.year_bin == 5 and "year_bin" in sub_cfg:
                args.year_bin = int(sub_cfg["year_bin"])
            if args.seed == DEFAULT_SEED and "seed" in sub_cfg:
                args.seed = int(sub_cfg["seed"])
        except ConfigError as exc:
            return fail(UsageError(str(exc)))
        except Exception as exc:  # noqa: BLE001
            _err(f"NOTE: could not read settings from {args.config}: {exc}")
    _say(f"wild host groups: {', '.join(sorted(wild_groups)) or '(none)'}")

    clades_df = (
        pd.read_csv(args.clades, sep="\t") if args.clades and args.clades.is_file() else None
    )
    meta_df = (
        pd.read_csv(args.metadata, sep="\t") if args.metadata and args.metadata.is_file() else None
    )
    if clades_df is None and args.mode == "clade":
        return fail(UsageError("--mode clade needs --clades (run `wildphylo tree clades` first)"))

    df = build_records(seqs, clades_df, meta_df)
    _say(f"Loaded {len(df)} aligned sequences")
    if (df["clade"] == "unassigned").any():
        n = int((df["clade"] == "unassigned").sum())
        _say(
            f"  note: {n} sequences had no clade assignment "
            "(present in the alignment but not the tree)"
        )

    # ---- select pool ----------------------------------------------------
    if args.mode == "clade":
        pool = df[df["clade"] == args.clade].copy()
        if pool.empty:
            _err(f"No sequences in {args.clade}. Available: {sorted(df['clade'].unique())}")
            return 1
        title = args.clade
    else:
        pool = df.copy()
        title = "global"

    summarise(pool, f"Pool ({title})", wild_groups)

    # ---- deduplicate ------------------------------------------------------
    if not args.keep_duplicates:
        pool, n_dropped = drop_duplicates(pool)
        if n_dropped:
            _say(
                f"\nRemoved {n_dropped} exact duplicate sequences "
                "(kept the best-dated representative of each)"
            )

    # ---- subsample --------------------------------------------------------
    target = args.target
    if args.mode == "global" and target is None:
        target = 400
    if target is not None and len(pool) > target:
        pool = stratified_sample(pool, target, args.year_bin, wild_groups)
        method = f"stratified subsample to {target}"
    else:
        method = "all sequences retained"
    _say(f"\n{method}")

    summarise(pool, "FINAL SUBSET", wild_groups)

    # sanity checks the analysis actually depends on
    hosts = Counter(pool["host_group"])
    thin = {h: n for h, n in hosts.items() if h in wild_groups and n < 5}
    if thin:
        _say(f"\n  WARNING: wild host groups with <5 sequences: {thin}")
        _say("  Discrete trait analysis will be unreliable for these. Consider merging")
        _say("  them into a broader category, or raising --target.")
    yrs = pool["decimal_year"].dropna()
    if len(yrs) and (yrs.max() - yrs.min()) < 10:
        _say(
            f"\n  WARNING: temporal span is only {yrs.max() - yrs.min():.1f} years — "
            "weak clock signal likely."
        )

    # ---- write -------------------------------------------------------------
    # Only create the fallback directory when it will actually be used;
    # otherwise every run leaves an empty data/processed/ behind.
    if not (args.out_aln and args.out_meta):
        args.outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or (f"H_{args.clade}" if args.mode == "clade" else f"H_global{len(pool)}")

    keep_labels = list(pool["label"])
    subset = {lbl: seqs[lbl] for lbl in keep_labels}
    if not args.keep_all_gap_columns:
        subset, n_stripped = strip_all_gap_columns(subset)
        if n_stripped:
            n_before = len(next(iter(seqs.values())))
            _say(
                f"\nStripped {n_stripped} all-gap columns "
                f"({n_before} -> {n_before - n_stripped}); they held only bases from "
                "sequences this subset dropped"
            )
    out_fasta = args.out_aln or args.outdir / f"{prefix}.fasta"
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(out_fasta, subset)

    out_meta = args.out_meta or args.outdir / f"{prefix}_metadata.tsv"
    out_meta.parent.mkdir(parents=True, exist_ok=True)
    pool.drop(columns=["seq_key"]).to_csv(out_meta, sep="\t", index=False)

    # Follow --out-meta when it is given, so no file escapes the build
    # directory. A rule that writes outside its declared outputs cannot be
    # tracked or cleaned, and the stray file silently goes stale.
    out_dates = (
        out_meta.with_name(out_meta.stem.replace("_metadata", "") + "_dates.tsv")
        if args.out_meta
        else args.outdir / f"{prefix}_dates.tsv"
    )
    out_dates.parent.mkdir(parents=True, exist_ok=True)
    with out_dates.open("w") as fh:
        fh.write("taxon\tdate\n")
        for lbl in keep_labels:
            fh.write(f"{lbl}\t{lbl.rsplit('|', 1)[1]}\n")

    _say(f"\nwrote {out_fasta}")
    _say(f"wrote {out_meta}")
    _say(f"wrote {out_dates}")
    _say("\nNOTE: this FASTA is already aligned and already has BEAST-style tip labels.")
    _say("Do NOT realign it or feed it to a step that expects unaligned sequences with")
    _say("bare accessions. Run IQ-TREE directly:")
    _say(f"\n    iqtree2 -s {out_fasta} -m MFP -B 1000 --alrt 1000 -T AUTO \\")
    _say(f"        --prefix {args.outdir / prefix}_ml -redo")
    _say("\nAll-gap columns have already been stripped (--keep-all-gap-columns to keep them).")
    _say("Re-check TempEst — a subset with better temporal balance often has")
    _say("markedly stronger clock signal than the full tree did.")
    return 0


# =========================================================================== #
# xml (was 07_make_beast_xml.py)
#
# Standard tip-dated analysis: HKY+G4 substitution model, clock model and tree
# prior from the config, tip dates parsed from the labels. Also writes a traits
# file (taxon -> host group) and can generate date-randomised replicates for the
# temporal signal test (DRT).
#
# Open the generated XML in BEAUti before running it (File > Load). BEAUti
# reports any incompatibility with your BEAST version and lets you inspect every
# prior. Do not run a chain for days without that check. This targets BEAST
# 2.7.x namespaces.
# =========================================================================== #

# Namespaces BEAST searches when resolving an unqualified spec= name.
# beast.base.inference.operator is required for ScaleOperator / DeltaExchangeOperator;
# beast.base.evolution.operator for the tree operators. Omitting either produces
# "Class could not be found" errors at load time.
NS = (
    "beast.base.evolution.alignment:beast.base.evolution.tree.coalescent:"
    "beast.base.util:beast.base.evolution.nuc:"
    "beast.base.evolution.operator:beast.base.evolution.operator.kernel:"
    "beast.base.inference.operator:beast.base.inference.operator.kernel:"
    "beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:"
    "beast.base.evolution.likelihood:beast.base.inference:"
    "beast.base.inference.parameter:beast.base.inference.distribution:"
    "beast.base.evolution.branchratemodel:beast.base.evolution.speciation"
)

TipInfo = dict[str, tuple[str, str, float]]  # label -> (accession, host group, decimal year)


def measured_slope_from(report: Path | None) -> float | None:
    """Root-to-tip slope from the temporal gate's JSON report, or None if unusable.

    On a tree whose branch lengths are substitutions per site this IS a
    clock-rate estimate, measured from the dataset actually being analysed
    rather than inherited from a paper about another virus.

    Returns None for a time-scaled tree (slope ~1.0), because there the slope is
    a units check, not a rate: using it as a prior would centre the clock on 1.0
    subs/site/year.
    """
    if not report or not Path(report).is_file():
        return None
    try:
        data = json.loads(Path(report).read_text())
        regs = data.get("regressions") or []
        if not regs:
            return None
        slope = float(regs[0].get("slope", 0.0))
    except Exception:  # noqa: BLE001
        return None
    if slope <= 0:
        return None
    if abs(slope - 1.0) < 0.05:
        _say(
            "NOTE: the temporal report's slope is ~1.0, so it came from a "
            "time-scaled tree and is not a substitution rate. Ignoring it; "
            "supply a rate in the config or run the gate on the ML tree."
        )
        return None
    return slope


def uncertainty_map(metadata: Path | None, labels: Iterable[str]) -> dict[str, float]:
    """Per-tip date uncertainty in years, for tips whose date is imprecise."""
    if not metadata or not Path(metadata).is_file():
        return {}
    import pandas as pd

    df = pd.read_csv(metadata, sep="\t")
    if "date_uncertainty_years" not in df.columns:
        return {}
    by_acc: dict[str, float] = {}
    for _, r in df.iterrows():
        acc = str(r.get("accession", "")).split(".")[0]
        try:
            u = float(r.get("date_uncertainty_years") or 0.0)
        except (TypeError, ValueError):
            u = 0.0
        by_acc[acc] = u
    out: dict[str, float] = {}
    for lbl in labels:
        u = by_acc.get(lbl.split("|")[0], 0.0)
        if u > 0.01:  # day-precision tips need no sampling
            out[lbl] = u
    return out


def parse_labels(seqs: Iterable[str]) -> tuple[TipInfo, list[str]]:
    """Parse ``ACCESSION|host|year`` labels into ({label: (acc, host, year)}, bad labels)."""
    out: TipInfo = {}
    bad: list[str] = []
    for label in seqs:
        parts = label.split("|")
        if len(parts) < 3:
            bad.append(label)
            continue
        try:
            year = float(parts[-1])
        except ValueError:
            bad.append(label)
            continue
        out[label] = (parts[0], parts[1], year)
    return out, bad


def build_xml(
    seqs: dict[str, str],
    info: TipInfo,
    rate: float,
    chain: int,
    log_every: int,
    prefix: str,
    seed_note: str = "",
    *,
    tree_prior: str = "coalescent_constant",
    clock_model: str = "relaxed_lognormal",
    date_policy: str = "midpoint",
    uncertainties: dict[str, float] | None = None,
    rate_note: str = "",
    trait: dict[str, Any] | None = None,
) -> str:
    """Assemble the BEAST XML as text.

    The model is not fixed. ``tree_prior``, ``clock_model`` and the tip-date
    policy come from the pathogen config, because a constant-size coalescent and
    a relaxed lognormal clock are choices, not constants; applying them to, say,
    a clade containing a mass-mortality outbreak is simply the wrong model.
    """
    seq_blocks = "\n".join(
        f'        <sequence id="seq_{lbl}" spec="Sequence" taxon="{lbl}" '
        f'totalcount="4" value="{s}"/>'
        for lbl, s in seqs.items()
    )

    dates = ",".join(f"{lbl}={info[lbl][2]:.4f}" for lbl in seqs)
    n = len(seqs)

    tp = tree_prior_blocks(tree_prior, prefix, n)
    ck = clock_blocks(clock_model, prefix, rate, n)
    td = tip_date_blocks(date_policy, prefix, uncertainties or {}, {k: info[k][2] for k in seqs})

    tp_state, tp_prior, tp_operators, tp_log = (
        tp["state"],
        tp["prior"],
        tp["operators"],
        tp["log"],
    )
    ck_state, ck_prior, ck_branchrate, ck_operators, ck_log, ck_ref = (
        ck["state"],
        ck["prior"],
        ck["branchrate"],
        ck["operators"],
        ck["log"],
        ck["ref"],
    )
    td_operators = td.get("operators", "")
    td_priors = td.get("priors", "")
    td_log = td.get("log", "")
    sampled_note = (
        f" ({td['n_sampled']} tips sampled within their window)" if td.get("n_sampled") else ""
    )

    # Discrete trait. Absent -> every fragment is empty and the XML is exactly
    # what it was before, so a sequence-only analysis is unaffected.
    tr = trait or {}
    tr_data = tr.get("data", "")
    tr_state = tr.get("state", "")
    tr_prior = tr.get("prior", "")
    tr_likelihood = tr.get("likelihood", "")
    tr_operators = tr.get("operators", "")
    tr_log = tr.get("log", "")
    tr_treelog = tr.get("treelog", "").replace("__LOGEVERY__", str(log_every))
    trait_note = (
        f"discrete trait  : {tr['n_states']} states, {tr['n_rates']} rates"
        if tr
        else "discrete trait  : none (sequence only)"
    )

    return f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Generated by wildphylo beast xml
     {n} sequences
     clock rate prior : {rate:g} subs/site/year
     rate provenance  : {rate_note}
     clock model      : {clock_model}
     tree prior       : {tree_prior}
     tip dates        : {date_policy}{sampled_note}
     {trait_note}
     {seed_note}
     VERIFY IN BEAUti BEFORE RUNNING. -->
<beast beautitemplate="Standard" beautistatus="" namespace="{NS}" required="" version="2.7">

    <data id="{prefix}" spec="Alignment" name="alignment">
{seq_blocks}
    </data>

{tr_data}

    <map name="Uniform">beast.base.inference.distribution.Uniform</map>
    <map name="Exponential">beast.base.inference.distribution.Exponential</map>
    <map name="LogNormal">beast.base.inference.distribution.LogNormalDistributionModel</map>
    <map name="Normal">beast.base.inference.distribution.Normal</map>
    <map name="Gamma">beast.base.inference.distribution.Gamma</map>
    <map name="OneOnX">beast.base.inference.distribution.OneOnX</map>
    <map name="prior">beast.base.inference.distribution.Prior</map>
    <map name="Beta">beast.base.inference.distribution.Beta</map>

    <run id="mcmc" spec="beast.base.inference.MCMC" chainLength="{chain}">

        <state id="state" spec="State" storeEvery="5000">
            <tree id="Tree.t:{prefix}" spec="beast.base.evolution.tree.Tree" name="stateNode">
                <trait id="dateTrait.t:{prefix}" spec="beast.base.evolution.tree.TraitSet"
                       traitname="date-forward" value="{dates}">
                    <taxa id="TaxonSet.{prefix}" spec="TaxonSet">
                        <alignment idref="{prefix}"/>
                    </taxa>
                </trait>
                <taxonset idref="TaxonSet.{prefix}"/>
            </tree>
            <parameter id="kappa.s:{prefix}" spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>
            <parameter id="gammaShape.s:{prefix}" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            <parameter id="freqParameter.s:{prefix}" spec="parameter.RealParameter" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>
{tp_state}
{ck_state}
{tr_state}
        </state>

        <init id="RandomTree.t:{prefix}" spec="RandomTree" estimate="false" initial="@Tree.t:{prefix}" taxa="@{prefix}">
            <populationModel id="ConstantPopulation0.t:{prefix}" spec="ConstantPopulation">
                <parameter id="randomPopSize.t:{prefix}" spec="parameter.RealParameter" name="popSize">1.0</parameter>
            </populationModel>
        </init>

        <distribution id="posterior" spec="CompoundDistribution">

            <distribution id="prior" spec="CompoundDistribution">
{tp_prior}
                <prior id="KappaPrior.s:{prefix}" name="distribution" x="@kappa.s:{prefix}">
                    <LogNormal id="LogNormalDistributionModel.0" name="distr" M="1.0" S="1.25"/>
                </prior>
                <prior id="GammaShapePrior.s:{prefix}" name="distribution" x="@gammaShape.s:{prefix}">
                    <Exponential id="Exponential.0" name="distr" mean="1.0"/>
                </prior>
                <!-- Rate prior. S=1.0 is wide enough that the data can move it.
                     Provenance: {rate_note} -->
{ck_prior}
{td_priors}
{tr_prior}
            </distribution>

            <distribution id="likelihood" spec="CompoundDistribution" useThreads="true">
                <distribution id="treeLikelihood.{prefix}" spec="TreeLikelihood" data="@{prefix}" tree="@Tree.t:{prefix}">
                    <siteModel id="SiteModel.s:{prefix}" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape.s:{prefix}">
                        <parameter id="mutationRate.s:{prefix}" spec="parameter.RealParameter" estimate="false" lower="0.0" name="mutationRate">1.0</parameter>
                        <parameter id="proportionInvariant.s:{prefix}" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>
                        <substModel id="hky.s:{prefix}" spec="HKY" kappa="@kappa.s:{prefix}">
                            <frequencies id="estimatedFreqs.s:{prefix}" spec="Frequencies" frequencies="@freqParameter.s:{prefix}"/>
                        </substModel>
                    </siteModel>
{ck_branchrate}
                </distribution>
{tr_likelihood}
            </distribution>
        </distribution>

        <!-- Classic operators: stable across BEAST 2.6/2.7. The Bactrian kernel
             variants mix slightly better but live in packages whose paths have
             moved between versions. -->
        <operator id="KappaScaler.s:{prefix}" spec="ScaleOperator" parameter="@kappa.s:{prefix}" scaleFactor="0.5" weight="0.1"/>
        <operator id="gammaShapeScaler.s:{prefix}" spec="ScaleOperator" parameter="@gammaShape.s:{prefix}" scaleFactor="0.5" weight="0.1"/>
        <operator id="FrequenciesExchanger.s:{prefix}" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:{prefix}"/>
{tp_operators}
{ck_operators}
{td_operators}
{tr_operators}

        <operator id="CoalescentConstantTreeScaler.t:{prefix}" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:{prefix}" weight="3.0"/>
        <operator id="CoalescentConstantTreeRootScaler.t:{prefix}" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:{prefix}" weight="3.0"/>
        <operator id="CoalescentConstantUniformOperator.t:{prefix}" spec="Uniform" tree="@Tree.t:{prefix}" weight="30.0"/>
        <operator id="CoalescentConstantSubtreeSlide.t:{prefix}" spec="SubtreeSlide" tree="@Tree.t:{prefix}" weight="15.0"/>
        <operator id="CoalescentConstantNarrow.t:{prefix}" spec="Exchange" tree="@Tree.t:{prefix}" weight="15.0"/>
        <operator id="CoalescentConstantWide.t:{prefix}" spec="Exchange" isNarrow="false" tree="@Tree.t:{prefix}" weight="3.0"/>
        <operator id="CoalescentConstantWilsonBalding.t:{prefix}" spec="WilsonBalding" tree="@Tree.t:{prefix}" weight="3.0"/>

        <logger id="tracelog" spec="Logger" fileName="{prefix}.log" logEvery="{log_every}" model="@posterior" sanitiseHeaders="true" sort="smart">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood.{prefix}"/>
            <log id="TreeHeight.t:{prefix}" spec="beast.base.evolution.tree.TreeStatLogger" tree="@Tree.t:{prefix}"/>
            <log idref="kappa.s:{prefix}"/>
            <log idref="gammaShape.s:{prefix}"/>
            <log idref="freqParameter.s:{prefix}"/>
{tp_log}
{ck_log}
{td_log}
{tr_log}
            <log id="rate.c:{prefix}" spec="beast.base.evolution.RateStatistic" branchratemodel="@{ck_ref}" tree="@Tree.t:{prefix}"/>
        </logger>

        <logger id="screenlog" spec="Logger" logEvery="{log_every}">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
        </logger>

        <logger id="treelog.t:{prefix}" spec="Logger" fileName="{prefix}.trees" logEvery="{log_every}" mode="tree">
            <log id="TreeWithMetaDataLogger.t:{prefix}" spec="beast.base.evolution.TreeWithMetaDataLogger" branchratemodel="@{ck_ref}" tree="@Tree.t:{prefix}"/>
        </logger>

{tr_treelog}

        <operatorschedule id="OperatorSchedule" spec="OperatorSchedule"/>
    </run>
</beast>
"""  # noqa: E501  # BEAST XML template; wrapping would change the generated file


def add_arguments_xml(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo beast xml``."""
    parser.add_argument(
        "--aln", type=Path, help="Alignment (default with --build-dir: subsampled/)"
    )
    parser.add_argument("--out", type=Path, help="Output directory (legacy)")
    parser.add_argument("--out-xml", type=Path, help="Explicit XML path.")
    parser.add_argument("--out-traits", type=Path, help="Explicit traits path.")
    parser.add_argument(
        "--config",
        type=Path,
        help="Pathogen config: clock rate, clock model, tree prior, "
        "chain length and date policy all come from here. Required with --build-dir.",
    )
    add_build_dir_arguments(parser, seed=True, rep=True)
    parser.add_argument(
        "--metadata", type=Path, help="Curated metadata, for per-tip date uncertainty."
    )
    parser.add_argument(
        "--temporal-report",
        type=Path,
        help="Temporal gate JSON. When the config has no "
        "clock_rate_prior, the prior is centred on the "
        "root-to-tip slope measured from THIS dataset.",
    )
    parser.add_argument(
        "--replicate", type=int, default=None, help="Date-randomisation replicate number."
    )
    parser.add_argument("--prefix", default=None)
    parser.add_argument(
        "--rate", type=float, default=None, help="Override the clock rate prior, subs/site/year."
    )
    parser.add_argument("--chain", type=int, default=None)
    parser.add_argument(
        "--log-every", type=int, default=None, help="Default: chain/10000, giving ~10k samples"
    )
    parser.add_argument(
        "--randomize-dates",
        type=int,
        nargs="?",
        const=1,
        default=0,
        help="Generate N date-randomised replicates for the DRT (bare flag: 1; "
        "with --replicate/--rep it only marks the intent, the replicate is single)",
    )
    parser.add_argument(
        "--no-trait",
        action="store_true",
        help="Write a sequence-only XML even when the config "
        "enables the discrete trait. Used for the DRT, where "
        "only the clock rate matters.",
    )


def _resolve_xml(args: argparse.Namespace) -> bool:
    """Apply ``--build-dir`` defaults for ``xml``; True when it is a DRT replicate.

    Two layouts: a real chain (``beast/<locus>/seed<SEED>/``, selected by
    ``--seed``) and a date-randomised replicate (``drt/<locus>/rep<REP>/``,
    selected by ``--rep`` or ``--replicate``). A replicate writes no traits file:
    the layout declares none, and an undeclared output goes stale unnoticed.
    """
    if args.build_dir is None:
        return False
    if args.rep is None and args.replicate is not None:
        args.rep = str(args.replicate)
    if args.replicate is None and args.rep is not None:
        try:
            args.replicate = int(args.rep)
        except ValueError as exc:
            msg = f"--rep must be a replicate number, got {args.rep!r}"
            raise UsageError(msg) from exc
    is_drt = args.rep is not None
    if args.randomize_dates and not is_drt and args.out_xml is None:
        msg = "--build-dir with --randomize-dates needs --rep (one replicate per call)"
        raise UsageError(msg)

    derived: dict[str, Callable[[Build], Any]] = {
        "aln": lambda b: b.layout.subsample(b.need("locus")),
        "metadata": lambda b: b.layout.subsample_metadata(b.need("locus")),
        "temporal_report": lambda b: b.layout.temporal_report(b.need("locus")),
    }
    outputs: tuple[str, ...]
    if is_drt:
        derived["out_xml"] = lambda b: b.layout.drt_xml(b.need("locus"), b.need("rep"))
        outputs = ("out_xml",)
    else:
        derived["out_xml"] = lambda b: b.layout.beast_xml(b.need("locus"), b.need("seed"))
        derived["out_traits"] = lambda b: b.layout.beast_traits(b.need("locus"), b.need("seed"))
        outputs = ("out_xml", "out_traits")
    resolve(args, derived=derived, outputs=outputs)
    return is_drt


def _locus_config(cfg: dict[str, Any], wanted: str | None) -> dict[str, Any]:
    """The locus block to analyse: the named (or first analysed) locus of a multi-locus config."""
    if "loci" in cfg:
        for locus in cfg["loci"]:
            if not locus.get("analyse"):
                continue
            if wanted is None or locus["name"] == wanted:
                return dict(locus)
        return {}
    return dict(cfg.get("locus", {}))


def _discrete_trait(
    cfg: dict[str, Any],
    dt: dict[str, Any],
    seqs: dict[str, str],
    info: TipInfo,
    prefix: str,
) -> dict[str, Any]:
    """Build the discrete-trait XML blocks for the tips in ``seqs``; raises ValueError.

    Previously this was a manual BEAUti step, which meant the model actually run
    was recorded in neither the config nor git, and the config said symmetric
    while the printed instructions said asymmetric, with no way to tell
    afterwards which had been clicked.

    ``hosts.dta_states`` maps host group to trait state, so several groups can
    collapse onto one state. A host group absent from the map has no declared
    state: it becomes '?', which the codeMap makes ambiguous across all states.
    Assigning it to a state instead would invent an observation.
    """
    state_map = cfg.get("hosts", {}).get("dta_states") or {}
    tip_states: dict[str, str] = {}
    unmapped: dict[str, int] = {}
    for k in seqs:
        host = info[k][1]
        st = state_map.get(host) if state_map else host
        if st is None:
            unmapped[host] = unmapped.get(host, 0) + 1
            st = "?"
        tip_states[k] = st
    if unmapped:
        _say(f"NOTE: host groups with no dta_states mapping, entered as ambiguous: {unmapped}")

    observed = sorted({v for v in tip_states.values() if v != "?"})
    declared_all = sorted(set(state_map.values())) if state_map else observed
    dropped = [s for s in declared_all if s not in observed]
    if dropped:
        # A state with no tips contributes rates nothing can inform, and
        # asymmetric BSSVS is already rate-hungry.
        _say(f"NOTE: trait states with no tips in this alignment, excluded: {dropped}")
    trait = trait_blocks(
        prefix,
        dt.get("trait", "host"),
        observed,
        tip_states,
        symmetric=bool(dt.get("symmetric", True)),
        bssvs=bool(dt.get("bssvs", True)),
        poisson_lambda=dt.get("poisson_lambda"),
    )
    sym = "symmetric" if dt.get("symmetric", True) else "asymmetric"
    bs = "BSSVS on" if dt.get("bssvs", True) else "BSSVS off"
    _say(f"discrete trait: {trait['n_states']} states, {trait['n_rates']} rates ({sym}, {bs})")
    per_rate = len(seqs) / trait["n_rates"]
    if per_rate < 10:
        _say(
            f"  WARNING: {len(seqs)} tips for {trait['n_rates']} rates "
            f"({per_rate:.1f} tips/rate). Most rates will be "
            "unidentified; report only those with decisive BF support."
        )
    _say()
    return dict(trait)


def _randomised_xml(
    seqs: dict[str, str],
    info: TipInfo,
    years: list[float],
    rng: random.Random,
    number: int,
    settings: dict[str, Any],
    seed_note: str,
) -> str:
    """BEAST XML for date-randomised replicate ``number``: tip dates shuffled across tips.

    The replicate is sequence-only and fixes dates at their (shuffled) values:
    the DRT asks whether the REAL dates carry signal, and sampling shuffled dates
    within their windows would blur the very thing being tested.
    """
    shuffled = list(years)
    rng.shuffle(shuffled)
    info_r = {k: (info[k][0], info[k][1], shuffled[j]) for j, k in enumerate(seqs)}
    return build_xml(
        seqs,
        info_r,
        settings["rate"],
        settings["chain"],
        settings["log_every"],
        f"{settings['prefix']}_drt{number}",
        seed_note=seed_note,
        tree_prior=settings["tree_prior"],
        clock_model=settings["clock_model"],
        date_policy="midpoint",
        rate_note=settings["rate_note"],
    )


@step("beast", "xml", "Generate a BEAST 2 XML (and DRT replicates) from a labelled alignment")
def xml(args: argparse.Namespace) -> int:
    """Generate a BEAST 2 XML from a labelled alignment, plus a traits file.

    A tip-dated analysis with tip dates parsed from the labels
    (``ACCESSION|host_group|decimal_year``). The tree prior, clock model and
    substitution model come from the ``beast:`` block of ``--config`` (defaults:
    constant-size coalescent, uncorrelated relaxed lognormal clock), with the clock
    rate prior from the config or the measured root-to-tip slope. A traits file
    (taxon -> host group) is written for the discrete trait analysis, which
    ``--no-trait`` omits. ``--randomize-dates N`` writes N date-randomised
    replicates for the date-randomisation test.

    IMPORTANT: open the generated XML in BEAUti before running it (File > Load).
    BEAUti reports any incompatibility with your BEAST version and lets you inspect
    every prior. Do not run a chain for days without that check. This step targets
    BEAST 2.7.x namespaces.
    """
    try:
        return _xml(args)
    except UsageError as exc:
        return fail(exc)
    except BeastStepError as exc:
        _err(str(exc))
        return 1


def _xml(args: argparse.Namespace) -> int:
    drt_layout = _resolve_xml(args)
    if args.aln is None:
        msg = "--aln is required (or --build-dir with --config, --locus, --seed)"
        raise UsageError(msg)
    seed = DEFAULT_SEED if args.seed is None else args.seed
    # --- config -------------------------------------------------------------
    cfg: dict[str, Any] = {}
    locus_cfg: dict[str, Any] = {}
    if args.config:
        try:
            cfg = load_config(args.config)
        except ConfigError as e:
            raise UsageError(str(e)) from e
        locus_cfg = _locus_config(cfg, args.locus)

    beast_cfg = cfg.get("beast", {})
    tree_prior = beast_cfg.get("tree_prior", "coalescent_constant")
    clock_model = beast_cfg.get("clock_model", "relaxed_lognormal")
    date_policy = (
        locus_cfg.get("imprecise_date_policy")
        or cfg.get("dates", {}).get("imprecise_date_policy")
        or "midpoint"
    )
    if args.chain is None:
        args.chain = int(beast_cfg.get("chain_length", 100_000_000))

    if not args.aln.is_file():
        raise BeastStepError(f"Alignment not found: {args.aln}")

    seqs = read_fasta(args.aln)
    info, bad = parse_labels(seqs)
    if bad:
        raise BeastStepError(
            f"ERROR: {len(bad)} labels not in ACCESSION|host|year form, e.g. {bad[:3]}"
        )
    seqs = {k: v for k, v in seqs.items() if k in info}

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise BeastStepError(
            f"ERROR: sequences differ in length {sorted(lengths)[:4]} — not aligned."
        )

    prefix = args.prefix or args.aln.stem.replace("-", "_")
    log_every = args.log_every or max(args.chain // 10000, 1000)

    # --- clock rate ---------------------------------------------------------
    # An explicit override, else the config, else the slope measured from this
    # dataset. There is no literature default: a CDV H-gene rate applied to
    # whatever pathogen happened to be running is how an assumption becomes a result.
    slope = measured_slope_from(args.temporal_report)
    if args.rate is not None:
        rate, rate_note = args.rate, "command line --rate"
    else:
        try:
            rate, rate_note = resolve_clock_rate(cfg, slope, locus_cfg)
        except ValueError as e:
            raise BeastStepError(f"ERROR: {e}") from e

    if cfg:
        checks = check_clock_prior(cfg, slope)
        _say(checks.render())
        if not checks.ok:
            raise BeastStepError(
                "\nRefusing to write an XML whose clock prior disagrees with "
                "the data. Resolve this before committing days of compute."
            )

    xml_path = args.out_xml or (args.out / f"{prefix}.xml" if args.out else Path(f"{prefix}.xml"))
    traits_path: Path | None = args.out_traits or (
        None
        if drt_layout
        else (args.out / f"{prefix}_traits.txt" if args.out else Path(f"{prefix}_traits.txt"))
    )
    xml_path.parent.mkdir(parents=True, exist_ok=True)
    if traits_path is not None:
        traits_path.parent.mkdir(parents=True, exist_ok=True)
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)

    years = [info[k][2] for k in seqs]
    hosts: dict[str, int] = {}
    for k in seqs:
        hosts[info[k][1]] = hosts.get(info[k][1], 0) + 1

    _say(f"{len(seqs)} sequences, {lengths.pop()} sites")
    _say(f"date range: {min(years):.2f} - {max(years):.2f} ({max(years) - min(years):.1f} yrs)")
    _say(f"host groups: {dict(sorted(hosts.items(), key=lambda kv: -kv[1]))}")
    _say(f"rate prior mean: {rate:g}   ({rate_note})")
    _say(f"clock model: {clock_model}   tree prior: {tree_prior}   tip dates: {date_policy}")
    _say(f"chain: {args.chain:,}   log every: {log_every:,}\n")

    thin = {h: n for h, n in hosts.items() if n < 5}
    if thin:
        _say(f"WARNING: host groups with <5 sequences: {thin}")
        _say("  Discrete trait analysis will be unreliable for these — consider merging.\n")

    # --- traits file for BEAUti ------------------------------------------
    if traits_path is not None:
        with traits_path.open("w") as fh:
            fh.write("traits\thost\n")
            for k in seqs:
                fh.write(f"{k}\t{info[k][1]}\n")
        note = (
            ""
            if (beast_cfg.get("discrete_trait") or {}).get("enabled") and not args.no_trait
            else "   <- import in BEAUti if adding a trait by hand"
        )
        _say(f"wrote {traits_path}{note}")

    # --- discrete trait -----------------------------------------------------
    trait: dict[str, Any] | None = None
    dt = beast_cfg.get("discrete_trait") or {}
    if dt.get("enabled") and not args.no_trait:
        try:
            trait = _discrete_trait(cfg, dt, seqs, info, prefix)
        except ValueError as e:
            raise BeastStepError(f"ERROR: {e}") from e

    # --- main XML ---------------------------------------------------------
    uncertainties = uncertainty_map(args.metadata, seqs) if date_policy == "interval" else {}
    if date_policy == "interval" and not uncertainties:
        _say(
            "NOTE: date policy is 'interval' but no per-tip uncertainty was "
            "found (pass --metadata). Tip dates will be fixed at their "
            "midpoints, which reports an HPD narrower than the data support."
        )

    xml_text = build_xml(
        seqs,
        info,
        rate,
        args.chain,
        log_every,
        prefix,
        tree_prior=tree_prior,
        clock_model=clock_model,
        date_policy=date_policy,
        uncertainties=uncertainties,
        rate_note=rate_note,
        trait=trait,
    )
    main_path = xml_path
    main_path.write_text(xml_text)
    try:
        ET.fromstring(xml_text)
        _say(f"wrote {main_path}   (well-formed XML)")
    except ET.ParseError as exc:
        raise BeastStepError(f"wrote {main_path}   BUT XML IS MALFORMED: {exc}") from exc

    # --- date randomised replicates ---------------------------------------
    drt = {
        "rate": rate,
        "chain": args.chain,
        "log_every": log_every,
        "prefix": prefix,
        "tree_prior": tree_prior,
        "clock_model": clock_model,
        "rate_note": rate_note,
    }
    if args.replicate is not None:
        # One replicate, written to --out-xml. The workflow fans these out as
        # separate jobs rather than producing 20 files from one invocation.
        rng = random.Random(seed + args.replicate)
        main_path.write_text(
            _randomised_xml(
                seqs,
                info,
                years,
                rng,
                args.replicate,
                drt,
                f"DATE-RANDOMISED REPLICATE {args.replicate} — dates "
                "shuffled for the temporal signal test. NOT a real analysis.",
            )
        )
        _say(f"wrote date-randomised replicate {args.replicate} -> {main_path}")
        return 0

    if args.randomize_dates:
        rng = random.Random(seed)
        drt_dir = (args.out / "date_randomised") if args.out else xml_path.parent
        drt_dir.mkdir(parents=True, exist_ok=True)
        for i in range(1, args.randomize_dates + 1):
            (drt_dir / f"{prefix}_drt{i}.xml").write_text(
                _randomised_xml(
                    seqs,
                    info,
                    years,
                    rng,
                    i,
                    drt,
                    f"DATE-RANDOMISED REPLICATE {i} — dates shuffled, "
                    "for the temporal signal test. NOT a real analysis.",
                )
            )
        _say(f"wrote {args.randomize_dates} date-randomised replicates -> {drt_dir}")
        _say("  Run these, then compare the ucldMean 95% HPDs against the real analysis.")
        _say("  Real estimate outside the randomised distribution = genuine temporal signal.")

    trait_steps = (
        ""
        if trait is None
        else f"""
   The discrete trait is already in the XML: {trait["n_states"]} states,
   {trait["n_rates"]} rates. Nothing to add in BEAUti.
"""
    )
    _say(f"""
NEXT STEPS
{trait_steps}
1. Smoke test before committing days of compute:

     sed 's/chainLength="{args.chain}"/chainLength="1000000"/' {main_path.name} > smoke.xml
     beast -seed 1 -threads 2 -overwrite smoke.xml

   Confirm it starts, the likelihood climbs, and — if the trait is in —
   that {prefix}_{(trait or {}).get("tag", "host")}.trees appears with
   ancestral states on its nodes. Then delete smoke.xml and its outputs.

2. Full run, one directory per seed so nothing overwrites anything:

     mkdir -p seed12345 seed54321
     caffeinate -i beast -seed 12345 -threads 2 -overwrite -prefix seed12345/ {main_path.name}
     caffeinate -i beast -seed 54321 -threads 2 -overwrite -prefix seed54321/ {main_path.name}

3. Convergence: every parameter ESS > 200, and the two chains agreeing with
   each other, not merely each converged on its own.

     loganalyser -b 10 seed12345/{prefix}.log

4. Date randomisation test (--no-trait; only the clock rate matters there),
   then TreeAnnotator for the MCC tree.
""")
    return 0


# =========================================================================== #
# Shared by jumps and state-through-time. Both read the same posterior of BEAST
# discrete-trait trees through wildphylo.newick: a transition count and a
# lineage share computed by two different parsers could not both be true of the
# same posterior.
# =========================================================================== #


def quantile(sorted_vals: Sequence[float], q: float) -> float:
    """Linear-interpolated quantile of an already-sorted list."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def mean(vals: Sequence[float]) -> float:
    """Arithmetic mean, NaN for an empty list."""
    return sum(vals) / len(vals) if vals else float("nan")


def sd(vals: Sequence[float]) -> float:
    """Sample standard deviation, NaN for fewer than two values."""
    if len(vals) < 2:
        return float("nan")
    m = mean(vals)
    return float((sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5)


def make_decoder(code_map: str | None) -> Callable[[str | None], str | None]:
    """Decoder from numeric trait codes to state names, from a BEAST ``codeMap`` string.

    Only needed when the trees store numeric codes instead of names. A state
    that is not a known code passes through unchanged.
    """
    lookup: dict[str, str] = {}
    if code_map:
        for chunk in code_map.split(","):
            if "=" not in chunk:
                continue
            name, codes = chunk.split("=", 1)
            code_list = codes.strip().split()
            if len(code_list) == 1 and code_list[0].isdigit():
                lookup[code_list[0]] = name.strip()

    def decode(state: str | None) -> str | None:
        return None if state is None else lookup.get(state, state)

    return decode


def expand_tree_files(patterns: Iterable[str]) -> list[Path]:
    """Expand glob patterns and plain paths into tree files; a unmatched pattern is fatal."""
    files: list[Path] = []
    for pattern in patterns:
        if any(ch in pattern for ch in "*?["):
            # user patterns may be absolute; Path.glob rejects those
            hits = [Path(h) for h in sorted(globmod.glob(pattern))]  # noqa: PTH207
        else:
            hits = [Path(pattern)]
        if not hits:
            raise BeastStepError(f"no file matched: {pattern}")
        files.extend(hits)
    return files


def _year_at(i: int, most_recent: float, max_depth: float, depth: list[float]) -> float:
    """Calendar year of node ``i``: the most recent tip year minus its distance from the tips."""
    return most_recent - (max_depth - depth[i])


def _most_recent_year(nodes: list[Node], tip_idx: list[int], given: float | None) -> float:
    """Decimal year of the most recent tip: ``given``, else the latest date in the tip labels."""
    if given is not None:
        return given
    dates = [d for d in (label_date(nodes[i]["label"]) for i in tip_idx) if d is not None]
    if not dates:
        raise BeastStepError("Tip labels carry no parseable date and --most-recent was not given.")
    return max(dates)


def _add_tree_arguments(parser: argparse.ArgumentParser) -> None:
    """Flags shared by jumps and state-through-time (``--trees`` ... ``--most-recent``)."""
    parser.add_argument(
        "--trees",
        nargs="+",
        default=None,
        help="one or more .host.trees files (.gz ok); pass every seed. "
        "With --build-dir the default is every seed in beast.seeds for --locus.",
    )
    parser.add_argument(
        "--burnin",
        type=float,
        default=0.10,
        help="fraction of trees discarded from the START of each file (default 0.10)",
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="annotation key holding the state (default: auto-detect; "
        "this repo's XML uses 'location')",
    )
    parser.add_argument(
        "--code-map",
        default=None,
        help='BEAST codeMap string, e.g. "domestic_dog=0,mustelid=1,..." '
        "— only needed if the trees store numeric codes",
    )
    parser.add_argument(
        "--most-recent",
        type=float,
        default=None,
        help="decimal year of the most recent tip (default: read from tip labels)",
    )
    parser.add_argument(
        "--config", type=Path, default=None, help="Pathogen config; required with --build-dir."
    )
    add_build_dir_arguments(parser)


def _resolve_trees(args: argparse.Namespace) -> list[Path]:
    """Trait-tree files: ``--trees``, else every ``beast.seeds`` chain under --build-dir."""

    def default(b: Build) -> list[Path]:
        locus = b.need("locus")
        try:
            return [b.layout.beast_trait_trees(locus, seed) for seed in b.seeds()]
        except ValueError as exc:
            msg = f"{exc}; enable beast.discrete_trait or pass --trees explicitly"
            raise UsageError(msg) from exc

    resolve(args, derived={"trees": default})
    if args.trees is None:
        msg = "--trees is required (or --build-dir with --config and --locus)"
        raise UsageError(msg)
    return expand_tree_files([str(t) for t in args.trees])


# =========================================================================== #
# jumps (was 08f_jump_history.py)
#
# Counts host-state transitions in BEAST2 DTA trait trees and records WHEN each
# one happened. Counting transitions alone (X -> Y happened N times) cannot
# distinguish a real transmission pattern from a sampling artefact. If the oldest
# tips in a dataset all carry one host label, ancestral reconstruction pulls the
# deep nodes toward that label, and every lineage descending from those deep
# nodes is scored as an outbound transition from it. The count is identical
# either way; only the timing of the transitions separates the two cases:
#
#   - jumps concentrated near the root, around the oldest sample dates
#       -> consistent with deep-node anchoring by the oldest tips
#   - jumps spread through the tree, including recent internal nodes
#       -> consistent with ongoing transmission between those hosts
#
# It parses the .host.trees file written by beastclassic's TreeWithTraitLogger
# (AncestralStateTreeLikelihood), walks every branch of every retained posterior
# tree, and writes one row per state change with the calendar year of the parent
# node, the child node, and the branch midpoint.
#
# It is a branch-level counter, not a full Markov-jump history. A branch whose
# two endpoints differ is scored as ONE transition even if the underlying
# continuous-time process made three. Counts are therefore a lower bound and
# slightly biased against long branches.
# =========================================================================== #

# (file, tree index, from, to, parent year, child year, mid year, branch length,
#  parent prob, child prob, child is tip)
Jump = tuple[str, int, str, str, float, float, float, float, float | None, float | None, bool]


class _JumpScan:
    """Everything collected while walking the retained trees."""

    def __init__(self, tag: str | None) -> None:
        self.jumps: list[Jump] = []
        self.tip_census: list[tuple[str, float | None]] | None = None  # from the first tree
        self.per_tree_counts: list[int] = []
        self.root_states: Counter[str | None] = Counter()
        self.trees_seen: dict[str, int] = {}
        self.trees_used = 0
        self.date_errors: list[float] = []
        self.missing_states = 0
        self.tag = tag


def _scan_jumps(args: argparse.Namespace, files: list[Path]) -> _JumpScan:
    """Walk every retained tree of every file, collecting branches whose endpoints differ."""
    decode = make_decoder(args.code_map)
    scan = _JumpScan(args.tag)

    for path in files:
        n_trees = count_trees(path)
        if n_trees == 0:
            raise BeastStepError(f"{path}: no trees found — is this a tree file?")
        skip = round(args.burnin * n_trees)
        scan.trees_seen[str(path)] = n_trees

        for tree_i, nodes, tag in iter_trees(path, scan.tag):
            scan.tag = tag
            if tree_i < skip:
                continue
            scan.trees_used += 1

            depth = node_depths(nodes)
            tip_idx = [i for i, nd in enumerate(nodes) if nd["is_tip"]]
            max_depth = max(depth[i] for i in tip_idx)
            most_recent = _most_recent_year(nodes, tip_idx, args.most_recent)

            year = partial(_year_at, most_recent=most_recent, max_depth=max_depth, depth=depth)

            # parser sanity check: do reconstructed tip dates match the labels?
            if scan.trees_used <= 5:
                for i in tip_idx:
                    d = label_date(nodes[i]["label"])
                    if d is not None:
                        scan.date_errors.append(abs(year(i) - d))

            if scan.tip_census is None:
                scan.tip_census = []
                for i in tip_idx:
                    st = decode(nodes[i]["state"])
                    if st is None:  # fall back to the host field of the label
                        parts = (nodes[i]["label"] or "").split("|")
                        st = parts[1] if len(parts) > 2 else "?"
                    scan.tip_census.append((st, label_date(nodes[i]["label"])))

            # tolerate a unary stem above the true root (some loggers emit one)
            r, stem = 0, set()
            while nodes[r]["state"] is None and nodes[r].get("nchild") == 1:
                stem.add(r)
                r = next(i for i, nd in enumerate(nodes) if nd["parent"] == r)
            scan.root_states[decode(nodes[r]["state"])] += 1
            scan.missing_states += sum(
                1 for i, nd in enumerate(nodes) if nd["state"] is None and i != r and i not in stem
            )

            n_here = 0
            for i, nd in enumerate(nodes):
                p = nd["parent"]
                if p is None:
                    continue
                s_from, s_to = decode(nodes[p]["state"]), decode(nd["state"])
                if s_from is None or s_to is None or s_from == s_to:
                    continue
                if args.min_prob is not None:
                    pp, cp = nodes[p]["prob"], nd["prob"]
                    if pp is not None and pp < args.min_prob:
                        continue
                    if cp is not None and cp < args.min_prob:
                        continue
                y_p, y_c = year(p), year(i)
                scan.jumps.append(
                    (
                        path.name,
                        tree_i,
                        s_from,
                        s_to,
                        y_p,
                        y_c,
                        (y_p + y_c) / 2.0,
                        nd["length"],
                        nodes[p]["prob"],
                        nd["prob"],
                        nd["is_tip"],
                    )
                )
                n_here += 1
            scan.per_tree_counts.append(n_here)
    return scan


def _write_jump_table(out: Path, jumps: list[Jump]) -> None:
    """Write the per-jump TSV, one row per state change."""
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(
            "file\ttree\tfrom\tto\tparent_year\tchild_year\tmid_year\t"
            "branch_length\tparent_prob\tchild_prob\tchild_is_tip\n"
        )
        fh.writelines(
            "{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.6f}\t{}\t{}\t{}\n".format(
                row[0],
                row[1],
                row[2],
                row[3],
                row[4],
                row[5],
                row[6],
                row[7],
                "" if row[8] is None else f"{row[8]:.4f}",
                "" if row[9] is None else f"{row[9]:.4f}",
                int(row[10]),
            )
            for row in jumps
        )


REPORT_WIDTH = 78


def _report_focus_jump(
    focus: str, jumps: list[Jump], by_type: dict[tuple[str, str], list[float]]
) -> None:
    """Report one FROM:TO transition against all others, with a by-decade breakdown."""
    try:
        fa, fb = focus.split(":")
    except ValueError:
        msg = "--focus must look like FROM:TO, e.g. wild_felid:procyonid"
        raise UsageError(msg) from None
    ys = sorted(by_type.get((fa, fb), []))
    focus_rows = [j for j in jumps if (j[2], j[3]) == (fa, fb)]
    if not ys:
        _say(f"\n[focus] no {fa} -> {fb} transitions recorded.")
        return
    others = sorted(y for (a, b), v in by_type.items() if (a, b) != (fa, fb) for y in v)
    _say(f"\n--- FOCUS: {fa} -> {fb} " + "-" * (REPORT_WIDTH - 13 - len(fa) - len(fb)))
    _say(f"  n = {len(ys)}  ({len(ys) / len(jumps):.1%} of all transitions)")
    _say(
        f"  jump year   median {quantile(ys, 0.5):.1f}, "
        f"IQR {quantile(ys, 0.25):.1f}-{quantile(ys, 0.75):.1f}, "
        f"range {ys[0]:.1f}-{ys[-1]:.1f}"
    )
    _say(
        f"  all others  median {quantile(others, 0.5):.1f}, "
        f"IQR {quantile(others, 0.25):.1f}-{quantile(others, 0.75):.1f}"
    )
    oldest_decade = sum(1 for y in ys if y < ys[0] + 10) / len(ys)
    _say(f"  share within 10 yr of its own oldest jump: {oldest_decade:.0%}")
    n_term = sum(1 for j in focus_rows if j[10])
    _say(f"  on terminal branches: {n_term / len(focus_rows):.0%}")

    _say("\n  by decade (count, and share on terminal branches):")
    lo_d = int(min(ys) // 10 * 10)
    hi_d = int(max(ys) // 10 * 10)
    for d in range(lo_d, hi_d + 10, 10):
        rows = [j for j in focus_rows if d <= j[6] < d + 10]
        if not rows:
            continue
        bar = "#" * max(1, round(40 * len(rows) / len(focus_rows)))
        t = sum(1 for j in rows if j[10]) / len(rows)
        _say(
            f"    {d}s {len(rows):>7} ({len(rows) / len(focus_rows):>4.0%}) term {t:>4.0%}  {bar}"
        )
    _say("\n  Two separated humps mean two processes are being pooled: deep")
    _say("  jumps from root anchoring, and shallow ones from something else.")
    _say("  If the shallow hump is mostly terminal, it is tip labels, not spread.")


def _report_jumps(args: argparse.Namespace, scan: _JumpScan) -> None:
    """Print the jump-history report: tip census, root state, ranked transitions, focus."""
    jumps = scan.jumps
    w = REPORT_WIDTH
    _say("=" * w)
    _say("HOST-TRANSITION HISTORY  (branch-level, height-stratified)")
    _say("=" * w)
    _say(f"state annotation key : {scan.tag}")
    for f, n in scan.trees_seen.items():
        _say(f"  {f}: {n} trees, {round(args.burnin * n)} discarded as burn-in")
    _say(f"trees analysed       : {scan.trees_used}")
    _say(f"transitions recorded : {len(jumps)}")
    _say(
        f"per tree             : mean {mean(scan.per_tree_counts):.2f}, "
        f"sd {sd(scan.per_tree_counts):.2f}"
    )
    if scan.missing_states:
        _say(
            f"WARNING              : {scan.missing_states} nodes carried no "
            f"'{scan.tag}' state and were skipped — transition counts are "
            "incomplete. Check --tag."
        )
    if scan.date_errors:
        _say(
            f"tip-date check       : max |reconstructed - label| = {max(scan.date_errors):.4f} yr "
            f"({'OK' if max(scan.date_errors) < 0.01 else 'MISMATCH — heights are suspect'})"
        )

    census: dict[str, list[float | None]] = defaultdict(list)
    for st, yr in scan.tip_census or []:
        census[st].append(yr)
    n_tips = sum(len(v) for v in census.values())

    _say("\n--- TIP CENSUS " + "-" * (w - 15))
    _say(f"  {'state':<16}{'tips':>6}{'% tips':>8}{'oldest':>10}{'newest':>10}")
    oldest_overall: float | None = None
    oldest_states: list[str] = []
    for st, yrs in sorted(census.items(), key=lambda kv: -len(kv[1])):
        ys = [y for y in yrs if y is not None]
        lo = min(ys) if ys else float("nan")
        hi = max(ys) if ys else float("nan")
        _say(f"  {st!s:<16}{len(yrs):>6}{len(yrs) / n_tips:>8.1%}{lo:>10.1f}{hi:>10.1f}")
        if ys and (oldest_overall is None or lo < oldest_overall):
            oldest_overall = lo
    if oldest_overall is not None:
        for st, yrs in census.items():
            ys = [y for y in yrs if y is not None]
            if ys and min(ys) <= oldest_overall + 1.0:
                oldest_states.append(st)

    _say("\n--- ROOT STATE (posterior) " + "-" * (w - 27))
    total_roots = sum(scan.root_states.values())
    top_state, top_n = scan.root_states.most_common(1)[0]
    for state, n in scan.root_states.most_common():
        _say(f"  {state!s:<16} {n / total_roots:6.1%}  ({n})")
    tip_share = len(census.get(top_state or "", [])) / n_tips if n_tips else float("nan")
    _say(
        f"\n  {top_state} holds the root in {top_n / total_roots:.0%} of trees "
        f"while carrying {tip_share:.0%} of the tips."
    )
    if oldest_states:
        oldest_names = ", ".join(sorted(oldest_states))
        _say(f"  Oldest tips (within 1 yr of {oldest_overall:.1f}): {oldest_names}")
    if top_state in oldest_states and tip_share < 0.25:
        _say("  READ: a minority state that also anchors the oldest tips is winning")
        _say("  the root. Outbound transitions from it are inflated by that anchoring;")
        _say("  treat its rank in the table as an upper bound, not an estimate.")

    mids = sorted(j[6] for j in jumps)
    cutoff = args.cutoff if args.cutoff is not None else quantile(mids, 0.5)

    by_type: dict[tuple[str, str], list[float]] = defaultdict(list)
    term_type: dict[tuple[str, str], int] = defaultdict(int)
    for j in jumps:
        by_type[(j[2], j[3])].append(j[6])
        if j[10]:
            term_type[(j[2], j[3])] += 1

    _say(f"\n--- TRANSITIONS, ranked  (cutoff for 'deep' = {cutoff:.1f}) " + "-" * 12)
    _say(f"  {'transition':<30}{'n':>7}{'%':>7}{'median':>8}{'IQR':>16}{'%deep':>7}{'%term':>7}")
    for (a, b), yrs_by_type in sorted(by_type.items(), key=lambda kv: -len(kv[1])):
        ys = sorted(yrs_by_type)
        pct = len(ys) / len(jumps)
        deep = sum(1 for y in ys if y < cutoff) / len(ys)
        term = term_type[(a, b)] / len(ys)
        _say(
            f"  {a + ' -> ' + b:<30}{len(ys):>7}{pct:>7.1%}"
            f"{quantile(ys, 0.5):>8.1f}"
            f"{f'{quantile(ys, 0.25):.1f}-{quantile(ys, 0.75):.1f}':>16}"
            f"{deep:>7.0%}{term:>7.0%}"
        )
    _say("  %deep = share older than the cutoff. %term = share landing on a")
    _say("  terminal branch, i.e. explaining one tip's label rather than any")
    _say("  onward transmission.")

    if args.focus:
        _report_focus_jump(args.focus, jumps, by_type)

    _say("\n" + "=" * w)
    _say("Branch-level counts are a LOWER BOUND: a branch with two endpoints in")
    _say("different states is scored once regardless of how many jumps the")
    _say("continuous-time process actually made along it.")
    _say("=" * w)


def add_arguments_jumps(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo beast jumps``."""
    _add_tree_arguments(parser)
    parser.add_argument(
        "--min-prob",
        type=float,
        default=None,
        help="drop jumps where either endpoint's state probability is below this "
        "(only works if the logger wrote <tag>.prob)",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=None,
        help="year separating 'deep' from 'shallow' jumps (default: median of all jump midpoints)",
    )
    parser.add_argument(
        "--focus",
        default=None,
        help="transition to interrogate, as FROM:TO, e.g. wild_felid:procyonid",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="write the per-jump TSV here (default: no file, summary only)",
    )


def _require_jumps(scan: _JumpScan) -> None:
    """Fail when no transition was found: the tag or code map is probably wrong."""
    if not scan.jumps:
        raise BeastStepError(
            "No transitions found. Check --tag / --code-map against the tree file."
        )


@step("beast", "jumps", "Date every host-state transition in a posterior of BEAST trait trees")
def jumps(args: argparse.Namespace) -> int:
    """Count host-state transitions in BEAST trait trees and record WHEN each happened.

    Counting transitions alone (X -> Y happened N times) cannot distinguish a real
    transmission pattern from a sampling artefact. If the oldest tips in a dataset
    all carry one host label, ancestral reconstruction pulls the deep nodes toward
    that label, and every lineage descending from those deep nodes is scored as an
    outbound transition from it. The count is identical either way; only the
    *timing* of the transitions separates the two cases:

      - jumps concentrated near the root, around the oldest sample dates:
        consistent with deep-node anchoring by the oldest tips
      - jumps spread through the tree, including recent internal nodes:
        consistent with ongoing transmission between those hosts

    This parses the ``.host.trees`` file written by beastclassic's
    TreeWithTraitLogger (AncestralStateTreeLikelihood), walks every branch of every
    posterior tree, and writes one row per state change with the calendar year of
    the parent node, the child node, and the branch midpoint.

    It is a branch-level counter, not a full Markov-jump history. A branch whose two
    endpoints differ is scored as ONE transition even if the underlying
    continuous-time process made three. Counts are therefore a lower bound and
    slightly biased against long branches. For exact jump histories you need BEAST's
    complete history logger; this is the diagnostic you can run on output you
    already have.

    What it reports:

      * tip-date consistency check (parser sanity: reconstructed tip dates vs the
        dates embedded in the tip labels)
      * a TIP CENSUS: how many tips carry each state and over what date range. This
        is the comparison that makes the root state interpretable: a state held by
        few, old tips that nonetheless wins the root is anchoring.
      * posterior distribution of the ROOT state
      * transitions ranked by frequency, with median and IQR of jump year
      * for each transition, the share of its jumps older than a cutoff year and the
        share falling on TERMINAL branches. A jump on a terminal branch describes
        one tip's own label and carries no information about onward transmission; a
        transition that is mostly terminal is not a corridor.
      * optional focus block (``--focus FROM:TO``) comparing one transition against
        all others, with a by-decade breakdown that exposes bimodal (mixture) timing

    If the annotation stores numeric codes instead of names, pass the codeMap string
    straight out of the BEAST XML with ``--code-map``. The most recent tip date is
    taken from the tip labels unless ``--most-recent`` is given.
    """
    try:
        scan = _scan_jumps(args, _resolve_trees(args))
        _require_jumps(scan)
        if args.out:
            _write_jump_table(args.out, scan.jumps)
        _report_jumps(args, scan)
    except UsageError as exc:
        return fail(exc)
    except BeastStepError as exc:
        _err(str(exc))
        return 1
    return 0


# =========================================================================== #
# state-through-time (was 08g_state_through_time.py)
#
# How many lineages carry each host state, year by year. `jumps` counts
# transitions and dates them, but a count has no denominator: a state can emit
# many transitions early simply because that is where the tree is narrow and
# every emission is visible. This step supplies the denominator. It slices the
# posterior trees at fixed time points and counts, at each slice, how many
# lineages are painted with each state (the discrete-trait version of a
# lineages-through-time plot). Read it as:
#
#   * a state ENRICHED relative to its tip share only in the deep tree, decaying
#     to roughly its tip share once more lineages exist, is being driven by the
#     root, not by ongoing presence in the tree
#   * a state that holds enrichment above its tip share across decades is a
#     structural feature of the topology, which a randomization test then has to
#     rule in or out as signal
#
# Enrichment, not raw share, is the quantity to read: 20% of lineages against a
# 10% tip share is 2x enriched, against a 40% tip share it is a depletion.
#
# Slices where the tree holds only one or two lineages cannot support a
# proportion at all, so slices below --min-lineages are excluded from both the
# table and the verdict.
#
# Method: each branch spans [parent_year, child_year] and has a state at each
# end. For a slice in the older half of a branch the parent's state is used, for
# the younger half the child's (the midpoint convention; an approximation for
# branches whose endpoints differ). Counts are averaged over the retained
# posterior trees. Cost scales with trees x branches x slices; --thin 10 (the
# default) keeps every 10th tree, which loses little because adjacent logged
# trees are highly autocorrelated.
# =========================================================================== #

# (parent year, child year, parent state, child state)
Branch = tuple[float, float, str | None, str | None]


def _read_branches(
    args: argparse.Namespace, files: list[Path]
) -> tuple[list[list[Branch]], list[float], Counter[str], str | None]:
    """Parse retained trees into flat branch lists.

    Returns (branches per tree, root years, tip census, resolved annotation tag).
    """
    decode = make_decoder(args.code_map)
    per_tree: list[list[Branch]] = []
    root_years: list[float] = []
    tip_census: Counter[str] | None = None
    resolved_tag: str | None = args.tag

    for path in files:
        n_trees = count_trees(path)
        if n_trees == 0:
            raise BeastStepError(f"{path}: no trees found")
        skip = round(args.burnin * n_trees)

        for tree_i, nodes, tag in iter_trees(path, resolved_tag):
            resolved_tag = tag
            if tree_i < skip or (tree_i - skip) % args.thin:
                continue
            depth = node_depths(nodes)
            tip_idx = [i for i, nd in enumerate(nodes) if nd["is_tip"]]
            max_depth = max(depth[i] for i in tip_idx)
            most_recent = _most_recent_year(nodes, tip_idx, args.most_recent)

            if tip_census is None:
                tip_census = Counter()
                for i in tip_idx:
                    tip_census[decode(nodes[i]["state"]) or "?"] += 1

            branches: list[Branch] = []
            for i, nd in enumerate(nodes):
                p = nd["parent"]
                year_i = _year_at(i, most_recent, max_depth, depth)
                if p is None:
                    root_years.append(year_i)
                    continue
                branches.append(
                    (
                        _year_at(p, most_recent, max_depth, depth),
                        year_i,
                        decode(nodes[p]["state"]),
                        decode(nd["state"]),
                    )
                )
            if not branches:
                continue
            per_tree.append(branches)

    if not per_tree:
        raise BeastStepError("No trees retained — check --burnin and --thin.")
    return per_tree, root_years, tip_census or Counter(), resolved_tag


def _lineage_counts(
    per_tree: list[list[Branch]], step_years: float
) -> tuple[npt.NDArray[np.float64], list[str], dict[str, npt.NDArray[np.float64]]]:
    """Mean lineages in each state on a fixed year grid, averaged over the trees."""
    import numpy as np

    y_lo = min(min(b[0] for b in br) for br in per_tree)
    y_hi = max(max(b[1] for b in br) for br in per_tree)
    grid = np.arange(np.floor(y_lo), np.ceil(y_hi) + step_years, step_years)

    states = sorted({s for br in per_tree for b in br for s in (b[2], b[3]) if s})
    counts = {s: np.zeros(len(grid)) for s in states}

    for branches in per_tree:
        for y0, y1, s0, s1 in branches:
            if y1 <= y0:
                continue
            mid = (y0 + y1) / 2.0
            # older half carries the parent state, younger half the child state
            for lo, hi, st in ((y0, mid, s0), (mid, y1, s1)):
                if st is None or hi <= lo:
                    continue
                i0 = int(np.searchsorted(grid, lo, side="left"))
                i1 = int(np.searchsorted(grid, hi, side="right"))
                if i1 > i0:
                    counts[st][i0:i1] += 1.0

    for s in states:
        counts[s] /= len(per_tree)
    return grid, states, counts


def _report_focus_state(
    args: argparse.Namespace,
    grid: npt.NDArray[np.float64],
    states: list[str],
    counts: dict[str, npt.NDArray[np.float64]],
    total: npt.NDArray[np.float64],
    tip_census: Counter[str],
    tip_share: dict[str, float],
) -> None:
    """Report one state's enrichment in the earliest, middle and latest thirds of the tree."""
    import numpy as np

    f = args.focus
    if f not in states:
        msg = f"--focus {f} not among states: {', '.join(states)}"
        raise UsageError(msg)
    share = counts[f] / np.where(total > 0, total, np.nan)
    valid = [k for k in range(len(grid)) if total[k] >= args.min_lineages]
    if not valid:
        _say(f"\n[focus] no slice holds {args.min_lineages:.0f}+ lineages.")
        return
    tip_total = sum(tip_census.values())
    ts = tip_share[f]
    _say(f"\n--- FOCUS: {f} " + "-" * (REPORT_WIDTH - 12 - len(f)))
    _say(f"  tip share: {ts:.0%}  ({tip_census[f]} of {tip_total} tips)")
    _say(
        f"  interpretable window: {grid[valid[0]]:.0f}-{grid[valid[-1]]:.0f} "
        f"({len(valid)} slices with >= {args.min_lineages:.0f} lineages)"
    )

    # thirds of the interpretable window, weighted by lineages present
    third = max(1, len(valid) // 3)
    blocks = [
        ("earliest", valid[:third]),
        ("middle", valid[third : 2 * third]),
        ("latest", valid[2 * third :]),
    ]
    enr: dict[str, float] = {}
    _say(f"\n  {'window':<10}{'years':>14}{'mean lin':>10}{'share':>9}{'enrichment':>12}")
    for name, ks in blocks:
        if not ks:
            continue
        w = np.array([total[k] for k in ks])
        sh = float(np.average([share[k] for k in ks], weights=w))
        enr[name] = sh / ts if ts else float("nan")
        _say(
            f"  {name:<10}{f'{grid[ks[0]]:.0f}-{grid[ks[-1]]:.0f}':>14}"
            f"{float(np.mean(w)):>10.1f}{sh:>9.0%}{enr[name]:>12.1f}x"
        )

    e0, e1 = enr.get("earliest", float("nan")), enr.get("latest", float("nan"))
    _say()
    if e0 >= 1.5 and e1 < 1.25:
        _say(f"  READ: {e0:.1f}x enrichment in the deep tree decaying to {e1:.1f}x. This state's")
        _say("  prominence is carried by its oldest tips. Outbound transition")
        _say("  counts from it are a counting consequence of that anchoring,")
        _say("  not a rate. Do not report them as transition frequencies.")
    elif e1 >= 1.25:
        _say(f"  READ: enrichment persists ({e0:.1f}x early, {e1:.1f}x late).")
        _say("  Not explained by deep-node anchoring. Whether it is signal or")
        _say("  label noise is now a question for tip-label randomization.")
    else:
        _say(f"  READ: no meaningful enrichment at any depth ({e0:.1f}x early, {e1:.1f}x late).")
        _say("  The state tracks its tip frequency; there is nothing here to")
        _say("  explain as transmission structure.")


def _report_state_through_time(
    args: argparse.Namespace,
    n_kept: int,
    root_years: list[float],
    tip_census: Counter[str],
    tag: str | None,
    grid: npt.NDArray[np.float64],
    states: list[str],
    counts: dict[str, npt.NDArray[np.float64]],
    total: npt.NDArray[np.float64],
) -> None:
    """Print the state-through-time report: share and enrichment tables, then the focus block."""
    import numpy as np

    w = REPORT_WIDTH
    _say("=" * w)
    _say("STATE THROUGH TIME  (mean lineages per posterior tree)")
    _say("=" * w)
    _say(f"state annotation key : {tag}")
    _say(f"trees retained       : {n_kept}  (burnin {args.burnin:.0%}, thin {args.thin})")
    _say(
        f"root year            : mean {np.mean(root_years):.1f}, "
        f"95% range {np.percentile(root_years, 2.5):.1f}-"
        f"{np.percentile(root_years, 97.5):.1f}"
    )
    if tip_census:
        n = sum(tip_census.values())
        _say(
            "tips                 : "
            + ", ".join(f"{s} {c} ({c / n:.0%})" for s, c in tip_census.most_common())
        )

    tip_total = sum(tip_census.values())
    tip_share = {
        s: (tip_census.get(s, 0) / tip_total if tip_total else float("nan")) for s in states
    }

    header = "year".ljust(7) + "lin".rjust(6) + "   " + "".join(s[:8].rjust(10) for s in states)
    _say("\n  " + header)
    suppressed = 0
    for k, yr in enumerate(grid):
        if abs(yr - round(yr)) > 1e-9 or round(yr) % 5:
            continue
        if total[k] < args.min_lineages:
            if total[k] > 0:
                suppressed += 1
            continue
        row = "".join(f"{counts[s][k] / total[k]:>9.0%} " for s in states)
        _say(f"  {yr:<7.0f}{total[k]:>6.1f}   {row}")
    _say("\n  Columns are the share of lineages in each state at that year.")
    _say("  'lin' is the mean number of lineages alive in the tree then.")
    if suppressed:
        _say(
            f"  {suppressed} slice(s) hidden: fewer than {args.min_lineages:.0f} "
            "lineages, so a percentage there describes one or two branches."
        )

    _say("\n  ENRICHMENT (share of lineages / share of tips; 1.0 = as expected)")
    _say("  " + header)
    for k, yr in enumerate(grid):
        if abs(yr - round(yr)) > 1e-9 or round(yr) % 5:
            continue
        if total[k] < args.min_lineages:
            continue
        row = ""
        for s in states:
            e = (counts[s][k] / total[k]) / tip_share[s] if tip_share[s] else float("nan")
            row += f"{e:>9.1f} "
        _say(f"  {yr:<7.0f}{total[k]:>6.1f}   {row}")
    _say("  A state well above 1.0 only in the deep tree, falling to ~1.0 once")
    _say("  lineages accumulate, is anchored by its oldest tips rather than")
    _say("  persistent. Sustained values above 1.0 are the case for structure.")

    if args.focus:
        _report_focus_state(args, grid, states, counts, total, tip_census, tip_share)

    _say("\n" + "=" * w)
    _say("Branch states come from BEAST's sampled reconstruction; this step")
    _say("summarises them, it does not re-estimate them. A state can be")
    _say("persistent here and still be noise — only a tip-label randomization")
    _say("test can establish that the trait carries signal at all.")
    _say("=" * w)


def add_arguments_state_through_time(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo beast state-through-time``."""
    _add_tree_arguments(parser)
    parser.add_argument(
        "--thin", type=int, default=10, help="keep every Nth post-burnin tree (default 10)"
    )
    parser.add_argument(
        "--min-lineages",
        type=float,
        default=3.0,
        help="suppress slices holding fewer than this many lineages; "
        "proportions there are not interpretable (default 3)",
    )
    parser.add_argument(
        "--step", type=float, default=1.0, help="slice interval in years (default 1.0)"
    )
    parser.add_argument(
        "--focus", default=None, help="state to summarise in the verdict lines, e.g. felid"
    )
    parser.add_argument("--out", type=Path, default=None, help="write the year x state table here")


@step(
    "beast",
    "state-through-time",
    "Lineages carrying each host state, year by year, across a posterior of trait trees",
)
def state_through_time(args: argparse.Namespace) -> int:
    """Count how many lineages carry each host state, year by year.

    ``beast jumps`` counts transitions and dates them. It cannot tell you whether a
    state that dominates the deep tree is a root-era phenomenon or a persistent
    feature of the reconstruction, because a count has no denominator: a state can
    emit many transitions early simply because that is where the tree is narrow and
    every emission is visible.

    This supplies the denominator. It slices the posterior trees at fixed time
    points and counts, at each slice, how many lineages are painted with each state:
    the discrete-trait version of a lineages-through-time plot. Read it as:

      * a state ENRICHED relative to its tip share only in the deep tree, decaying to
        roughly its tip share once more lineages exist, is being driven by the root,
        not by ongoing presence in the tree
      * a state that holds enrichment above its tip share across decades is a
        structural feature of the topology, which a randomisation test then has to
        rule in or out as signal

    Enrichment, not raw share, is the quantity to read. A state holding 20% of
    lineages means nothing on its own; 20% against a 10% tip share is a 2x
    enrichment, and against a 40% tip share is a depletion.

    Slices where the tree holds only one or two lineages cannot support a proportion
    at all ("40% of 2 lineages" is one lineage), so slices below ``--min-lineages``
    are excluded from both the table and the verdict.

    Method: each branch spans [parent_year, child_year] and has a state at each end.
    For a slice falling in the older half of a branch the parent's state is used, for
    the younger half the child's. That is the standard midpoint convention and it is
    an approximation: a branch whose endpoints differ made its transition at some
    unknown point along its length, and the midpoint is the expectation under a
    uniform prior, not the truth. Branches with identical endpoints (the large
    majority) are unaffected. Counts are averaged over the retained posterior trees,
    so a "mean lineages" of 3.4 means that at that moment the average tree had 3.4
    lineages in that state.

    Performance: cost scales with trees x branches x slices. ``--thin 10`` (the
    default) keeps every 10th tree, which for a 10,000-tree run leaves about 900 per
    seed and is far more than enough: adjacent logged trees are highly
    autocorrelated and add almost nothing. Use ``--thin 1`` for the full file.
    ``--focus STATE`` focuses the verdict lines on one state.
    """
    import numpy as np

    try:
        per_tree, root_years, tip_census, tag = _read_branches(args, _resolve_trees(args))
        grid, states, counts = _lineage_counts(per_tree, args.step)
        total = np.sum([counts[s] for s in states], axis=0)

        if args.out:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            with args.out.open("w") as fh:
                fh.write("year\tstate\tmean_lineages\tproportion\n")
                for k, yr in enumerate(grid):
                    if total[k] <= 0:
                        continue
                    for s in states:
                        fh.write(
                            f"{yr:.1f}\t{s}\t{counts[s][k]:.4f}\t{counts[s][k] / total[k]:.4f}\n"
                        )

        _report_state_through_time(
            args, len(per_tree), root_years, tip_census, tag, grid, states, counts, total
        )
    except UsageError as exc:
        return fail(exc)
    except BeastStepError as exc:
        _err(str(exc))
        return 1
    return 0

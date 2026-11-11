"""align stage: look at an alignment before spending compute on a tree.

``wildphylo align qc`` is the only step; its docstring says what it checks.
"""

from __future__ import annotations

import logging
import sys
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from wildphylo.cli import step
from wildphylo.config import ConfigError
from wildphylo.config import load as load_config
from wildphylo.layout import BuildLayout

if TYPE_CHECKING:
    import argparse

log = logging.getLogger(__name__)

AMBIG = set("NnXx?")
GAP = set("-.")
STOPS = {"TAA", "TAG", "TGA"}

# Built-in thresholds, used when neither the command line nor the config sets one.
DEFAULT_GAP_THRESHOLD = 0.5
DEFAULT_IDENTITY_THRESHOLD = 0.80
DEFAULT_AMBIGUITY_THRESHOLD = 0.05
DEFAULT_COVERAGE_THRESHOLD = 0.5
DEFAULT_OUTDIR = Path("data/processed")
_CONFIG_KEYS = {"max_gap_fraction", "min_identity_to_consensus", "max_ambiguity_fraction"}


def _say(text: str = "", *, err: bool = False) -> None:
    """Write one line of the step's report to stdout (or stderr)."""
    print(text, file=sys.stderr if err else sys.stdout)


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA file into ``{first word of header: sequence}``."""
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip()
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
    """Write ``seqs`` as FASTA, wrapped at ``width`` columns."""
    with path.open("w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def consensus(seqs: dict[str, str]) -> str:
    """Majority base per column, ignoring gaps and ambiguity codes."""
    if not seqs:
        return ""
    length = len(next(iter(seqs.values())))
    out = []
    for i in range(length):
        col = Counter(s[i].upper() for s in seqs.values() if s[i] not in GAP and s[i] not in AMBIG)
        out.append(col.most_common(1)[0][0] if col else "-")
    return "".join(out)


def identity_to(seq: str, ref: str) -> float:
    """Fraction of comparable positions that match. Gaps in either are skipped."""
    same = comp = 0
    # strict=False: both are alignment rows, equal length in practice;
    # truncating to the shorter was the prior behaviour.
    for a, b in zip(seq.upper(), ref, strict=False):
        if a in GAP or b in GAP or a in AMBIG:
            continue
        comp += 1
        same += a == b
    return same / comp if comp else 0.0


def count_stops(seq: str, frame: int = 0) -> int:
    """Count internal stop codons of ``seq`` in ``frame``, ignoring a terminal stop."""
    ungapped = "".join(c for c in seq.upper() if c not in GAP)[frame:]
    codons = [ungapped[i : i + 3] for i in range(0, len(ungapped) - 2, 3)]
    return sum(1 for c in codons[:-1] if c in STOPS)  # ignore a terminal stop


def _report_missing_alignment(aln: Path) -> None:
    """Explain a missing alignment file and what to run to get one."""
    _say(f"Alignment not found: {aln}\n", err=True)
    parent = aln.parent
    if parent.is_dir():
        found = sorted(p.name for p in parent.glob("*.fasta"))
        _say(f"FASTA files in {parent}: {found or '(none)'}", err=True)
    else:
        _say(f"Directory {parent} doesn't exist either.", err=True)
    _say(
        "\nRun the alignment step first: the Snakefile's align rule "
        "(snakemake <build>/aligned/<locus>_aligned.fasta).",
        err=True,
    )


def _report_coverage(cov: list[float], threshold: float) -> list[int]:
    """Print the well-covered column span; return the indices of those columns."""
    n_cols = len(cov)
    good = [i for i, c in enumerate(cov) if c >= threshold]
    if good:
        first, last = good[0], good[-1]
        _say(f"Columns at >= {threshold:.0%} coverage: {first + 1} to {last + 1} (1-based)")
        if first > 0 or last < n_cols - 1:
            _say(f"  ragged: {first} leading and {n_cols - 1 - last} trailing columns are sparse")
            _say(f"  suggested:  --trim-to {first + 1} {last + 1} --write-trimmed")
    _say(f"Mean column coverage: {sum(cov) / n_cols:.1%}")
    return good


def _report_all_gap_columns(cov: list[float], good: list[int]) -> None:
    """Report columns that are a gap in every sequence.

    Such a column carries no data but still occupies a position. Anything that
    treats a column span as a coordinate (a CDS annotation, a codon partition,
    a site index in a figure caption) is silently shifted by it. These usually
    appear when an alignment is subset without re-stripping, so the columns
    held bases from sequences now gone.
    """
    n_cols = len(cov)
    empty = [i for i, c in enumerate(cov) if c == 0]
    if not empty:
        _say("No all-gap columns.")
        _say()
        return
    runs: list[tuple[int, int]] = []
    start, prev = empty[0], empty[0]
    for i in empty[1:]:
        if i != prev + 1:
            runs.append((start, prev))
            start = i
        prev = i
    runs.append((start, prev))
    interior = [
        (a, b)
        for a, b in runs
        if a > (good[0] if good else 0) and b < (good[-1] if good else n_cols - 1)
    ]
    _say(
        f"\n{len(empty)} all-gap columns in {len(runs)} run(s)"
        + (f", {sum(b - a + 1 for a, b in interior)} of them interior" if interior else "")
    )
    for a, b in runs[:8]:
        _say(f"    columns {a + 1}-{b + 1}  ({b - a + 1} nt)")
    if len(runs) > 8:
        _say(f"    ... and {len(runs) - 8} more")
    _say("  These shift every column coordinate after them. Remove with")
    _say("  --strip-all-gap --write-stripped, or re-run `wildphylo beast subsample`, which")
    _say("  now strips them when it writes a subset.")
    _say()


def _flag(df: Any, mask: Any, label: str) -> int:
    """Append ``label`` to the flag column of the rows in ``mask``; return how many."""
    df.loc[mask, "flag"] = (df.loc[mask, "flag"] + ";" + label).str.lstrip(";")
    return int(mask.sum())


def _sequence_table(seqs: dict[str, str], n_cols: int) -> Any:
    """Per-sequence ungapped length, gap and ambiguity fractions, identity to consensus."""
    cons = consensus(seqs)
    rows = []
    for name, s in seqs.items():
        ung = "".join(c for c in s if c not in GAP)
        rows.append(
            {
                "name": name,
                "ungapped_len": len(ung),
                "gap_frac": sum(1 for c in s if c in GAP) / n_cols,
                "ambig_frac": (sum(1 for c in ung if c in AMBIG) / len(ung)) if ung else 1.0,
                "identity_to_consensus": identity_to(s, cons),
            }
        )
    return pd.DataFrame(rows)


def _flag_sequences(
    df: Any, seqs: dict[str, str], gap_threshold: float, identity: float, ambiguity: float
) -> None:
    """Flag suspect sequences in ``df`` and print the flag counts.

    The alignment's reading frame is determined from the whole dataset: the
    frame with the fewest total internal stops. Sequences that still have stops
    in THAT frame are the genuinely suspect ones (frameshift, wrong strand,
    pseudogene, sequencing error). Testing each sequence against its own best
    frame is meaningless, because frames 1 and 2 of any real coding sequence
    are full of stops.
    """
    df["flag"] = ""
    n_gap = _flag(df, df["gap_frac"] > gap_threshold, "mostly_gaps")
    n_amb = _flag(df, df["ambig_frac"] > ambiguity, f"ambiguous>{ambiguity * 100:g}%")
    n_id = _flag(df, df["identity_to_consensus"] < identity, "low_identity")
    frame_totals = {f: sum(count_stops(s, f) for s in seqs.values()) for f in (0, 1, 2)}
    aln_frame = min(frame_totals, key=lambda f: frame_totals[f])
    df["stops_in_aln_frame"] = [count_stops(seqs[n], aln_frame) for n in df["name"]]
    n_st = _flag(df, df["stops_in_aln_frame"] > 0, "internal_stop")

    _say("Flagged sequences")
    _say(f"  mostly gaps (>{gap_threshold:.0%})        : {n_gap}")
    _say(f"  ambiguous bases (>{ambiguity:.0%})       : {n_amb}")
    _say(f"  identity to consensus <{identity:.0%}      : {n_id}")
    _say(f"  internal stops in frame {aln_frame}          : {n_st}")
    _say()
    _say(
        f"  Alignment reading frame inferred as {aln_frame} "
        f"(total internal stops per frame: {frame_totals})."
    )
    _say("  If frame 0 isn't the winner, the alignment starts mid-codon — normal for")
    _say("  GenBank records, but use this offset when setting up codon partitions.")
    _say("  Sequences with stops in the inferred frame are worth checking: frameshift,")
    _say("  wrong strand, or sequencing error.")
    _say()

    flagged = df[df["flag"] != ""].sort_values("identity_to_consensus")
    if not flagged.empty:
        _say("Worth looking at by eye (worst first):")
        _say(
            flagged[["name", "ungapped_len", "gap_frac", "identity_to_consensus", "flag"]]
            .head(20)
            .to_string(index=False)
        )
        _say()


def _report_duplicates(seqs: dict[str, str]) -> None:
    """Report sets of exactly identical (ungapped) sequences."""
    by_seq: dict[str, list[str]] = {}
    for name, s in seqs.items():
        by_seq.setdefault("".join(c for c in s.upper() if c not in GAP), []).append(name)
    dups = {k: v for k, v in by_seq.items() if len(v) > 1}
    if not dups:
        return
    _say(f"{len(dups)} sets of identical sequences ({sum(len(v) for v in dups.values())} seqs):")
    for names in list(dups.values())[:5]:
        _say("   " + " " + ", ".join(names[:4]) + (" ..." if len(names) > 4 else ""))
    _say("  Identical sequences are legitimate data, not errors — but they")
    _say("  contribute no phylogenetic signal and slow BEAST down. Consider")
    _say("  keeping one per host/year/country combination.\n")


def _report_label_check(seqs: dict[str, str], metadata: Path) -> None:
    """Cross-check tip labels (accession, host group, year) against the metadata table.

    Tip labels encode host group and date, so they are a copy of the metadata
    that can go stale independently of it. A label that disagrees with the
    table is a real error: it will silently propagate into the tree, the BEAST
    run and every figure, and nothing downstream re-derives it.
    """
    # dtype=str does NOT make missing cells strings: pandas leaves them as
    # float NaN, so any downstream .strip() raises AttributeError. Reference
    # anchors are the first records to carry an empty decimal_year, which is
    # how this surfaced. Fill so every cell really is a string.
    meta_df = pd.read_csv(metadata, sep="\t", dtype=str).fillna("")
    if "accession" not in meta_df.columns:
        return
    by_acc = {
        str(a).split(".")[0]: r
        for a, r in zip(meta_df["accession"], meta_df.to_dict("records"), strict=True)
    }
    unknown: list[str] = []
    host_bad: list[tuple[str, str, str]] = []
    date_bad: list[tuple[str, str, str]] = []
    anchors_seen: list[str] = []
    for name in seqs:
        parts = name.split("|")
        if len(parts) < 3:
            continue
        acc, host, year = parts[0].split(".")[0], parts[1], parts[-1]
        # Reference anchors carry their LINEAGE in the host field and a
        # non-numeric date, by design. Cross-checking them against the
        # host and date columns reports a disagreement that is not one,
        # which would train the reader to ignore this whole section.
        if host.startswith("ref_"):
            anchors_seen.append(name)
            continue
        row = by_acc.get(acc)
        if row is None:
            unknown.append(name)
            continue
        m_host = str(row.get("host_group") or "").strip()
        if m_host and m_host != host:
            host_bad.append((name, host, m_host))
        m_year = str(row.get("decimal_year") or "").strip()
        try:
            if m_year and abs(float(m_year) - float(year)) > 0.01:
                date_bad.append((name, year, m_year))
        except ValueError:
            pass

    _say("Labels vs metadata")
    if anchors_seen:
        _say(f"  reference anchors (not cross-checked): {len(anchors_seen)}")
    _say(f"  accession not in metadata          : {len(unknown)}")
    _say(f"  host group disagrees with metadata : {len(host_bad)}")
    _say(f"  date disagrees with metadata       : {len(date_bad)}")
    for name, in_label, in_meta in (host_bad + date_bad)[:10]:
        _say(f"    {name}")
        _say(f"      label says {in_label!r}, metadata says {in_meta!r}")
    for name in unknown[:5]:
        _say(f"    {name}  (no metadata row)")
    if host_bad or date_bad:
        _say("  Fix the label, not the metadata, unless the table is the stale one.")
        _say("  Relabelling changes tip names, so anything built from these labels")
        _say("  — trees, BEAST XML, exported datasets — has to be regenerated.")
    _say()


def _report_host_composition(seqs: dict[str, str]) -> None:
    """Print how many sequences each host group (second label field) contributes."""
    groups = Counter(n.split("|")[1] for n in seqs if "|" in n)
    if not groups:
        return
    _say("Host composition of the alignment:")
    for g, c in groups.most_common():
        _say(f"    {g:<16} {c:>5}")
    thin = [g for g, c in groups.items() if c < 5]
    if thin:
        _say(f"  WARNING: fewer than 5 sequences in: {', '.join(thin)}")
    _say()


def _resolve_thresholds(args: argparse.Namespace) -> None:
    """Settle ``args.*_threshold`` before anything uses them.

    Precedence per threshold: the command-line flag if it was given (even when
    it equals the built-in default), else ``alignment.*`` from ``--config``,
    else the built-in default. The flags and the plot read the settled values,
    so they cannot disagree.
    """
    from_config: dict[str, float] = {}
    if args.config:
        try:
            alignment_cfg = load_config(args.config).get("alignment", {})
            from_config = {k: float(v) for k, v in alignment_cfg.items() if k in _CONFIG_KEYS}
        except Exception as exc:  # noqa: BLE001 -- any unreadable config falls back to the defaults
            _say(f"NOTE: could not read thresholds from {args.config}: {exc}")
    for attr, key, default in (
        ("gap_threshold", "max_gap_fraction", DEFAULT_GAP_THRESHOLD),
        ("identity_threshold", "min_identity_to_consensus", DEFAULT_IDENTITY_THRESHOLD),
        ("ambiguity_threshold", "max_ambiguity_fraction", DEFAULT_AMBIGUITY_THRESHOLD),
        ("coverage_threshold", None, DEFAULT_COVERAGE_THRESHOLD),
    ):
        if getattr(args, attr) is None:
            setattr(args, attr, from_config.get(key, default) if key else default)


def _write_coverage_plot(cov: list[float], df: Any, args: argparse.Namespace, n_cols: int) -> None:
    """Write the per-column coverage and identity-to-consensus plot."""
    try:
        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pyplot as plt

        _fig, ax = plt.subplots(2, 1, figsize=(10, 5.5))
        ax[0].fill_between(range(1, n_cols + 1), cov, color="#4B3A79", lw=0)
        ax[0].axhline(args.coverage_threshold, color="#C0697F", ls="--", lw=1)
        ax[0].set_title("Per-column coverage (fraction of sequences with a base)")
        ax[0].set_xlabel("alignment column")
        ax[0].set_ylim(0, 1)
        ax[1].hist(df["identity_to_consensus"], bins=40, color="#6E5CA0")
        ax[1].axvline(args.identity_threshold, color="#C0697F", ls="--", lw=1)
        ax[1].set_title("Identity to consensus")
        plt.tight_layout()
        png = args.out_png or args.plot or args.outdir / (args.aln.stem + "_qc.png")
        png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(png, dpi=130)
        _say(f"coverage plot          -> {png}")
    except ImportError:
        _say("(matplotlib not installed — skipping the plot)")


def _strip_all_gap_columns(
    seqs: dict[str, str], cov: list[float], n_cols: int, args: argparse.Namespace
) -> None:
    """Optionally remove all-gap columns, writing the result with ``--write-stripped``."""
    keep = [i for i in range(n_cols) if cov[i] > 0]
    if len(keep) == n_cols:
        _say("\nNo all-gap columns to strip.")
        return
    stripped = {n: "".join(s[i] for i in keep) for n, s in seqs.items()}
    _say(f"\nStripped {n_cols - len(keep)} all-gap columns ({n_cols} -> {len(keep)})")
    if args.write_stripped:
        out = args.outdir / (args.aln.stem + "_stripped.fasta")
        write_fasta(out, stripped)
        _say(f"wrote {out}")
    else:
        _say("(add --write-stripped to save it)")


def _trim_columns(seqs: dict[str, str], args: argparse.Namespace) -> None:
    """Optionally trim to ``--trim-to`` columns, writing the result with ``--write-trimmed``."""
    s0, s1 = args.trim_to
    trimmed = {n: s[s0 - 1 : s1] for n, s in seqs.items()}
    empty = [n for n, s in trimmed.items() if set(s) <= GAP]
    if empty:
        _say(f"\n{len(empty)} sequences become all-gap after trimming and were dropped:")
        _say("    " + ", ".join(empty[:6]))
        trimmed = {n: s for n, s in trimmed.items() if n not in empty}
    _say(f"\nTrimmed to columns {s0}-{s1}: {len(trimmed)} sequences, {s1 - s0 + 1} columns")
    if args.write_trimmed:
        out = args.outdir / (args.aln.stem + "_trimmed.fasta")
        write_fasta(out, trimmed)
        _say(f"wrote {out}")
    else:
        _say("(add --write-trimmed to save it)")


class _UsageError(Exception):
    """A command line that cannot be resolved to paths; reported, exit code 2."""


def _resolve_paths(args: argparse.Namespace) -> None:
    """Fill in the paths the command line left out, in place.

    Without ``--build-dir`` nothing is derived: ``--aln`` is required and the
    other paths behave as before. With it, every path flag not given explicitly
    defaults to the matching :class:`BuildLayout` method, so
    ``--config C --build-dir B --locus L`` is a complete invocation. An explicit
    flag always wins. ``--plot`` counts as an explicit ``--out-png``.
    """
    if args.build_dir is None:
        if args.locus is not None:
            msg = "--locus only applies together with --build-dir"
            raise _UsageError(msg)
        if args.aln is None:
            msg = "--aln is required unless --build-dir is given"
            raise _UsageError(msg)
        if args.outdir is None:
            args.outdir = DEFAULT_OUTDIR
        return

    if args.config is None:
        msg = "--build-dir needs --config (the layout depends on the pathogen id)"
        raise _UsageError(msg)
    try:
        layout = BuildLayout.at(args.build_dir, load_config(args.config))
    except ConfigError as exc:
        raise _UsageError(str(exc)) from exc

    needs_locus = [
        flag
        for flag, value in (
            ("--aln", args.aln),
            ("--out-tsv", args.out_tsv),
            ("--out-png", args.out_png or args.plot),
        )
        if value is None
    ]
    if needs_locus and args.locus is None:
        msg = f"--build-dir needs --locus to derive {', '.join(needs_locus)} (or give them)"
        raise _UsageError(msg)

    locus = args.locus
    if args.aln is None:
        args.aln = layout.alignment(locus)
    if args.metadata is None:
        args.metadata = layout.metadata()
    if args.out_tsv is None:
        args.out_tsv = layout.alignment_qc(locus)
    if args.out_png is None and args.plot is None:
        args.out_png = layout.alignment_qc_plot(locus)
    if args.outdir is None:
        args.outdir = args.out_tsv.parent
    log.debug(
        "resolved: aln=%s metadata=%s out_tsv=%s out_png=%s outdir=%s",
        args.aln,
        args.metadata,
        args.out_tsv,
        args.out_png or args.plot,
        args.outdir,
    )


def add_arguments_qc(parser: argparse.ArgumentParser) -> None:
    """``wildphylo align qc (--aln ALN | --config YAML --build-dir DIR --locus L) [...]``."""
    parser.add_argument(
        "--aln", type=Path, default=None, help="Alignment FASTA (default with --build-dir: layout)"
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=None,
        help="Metadata TSV for the label check (default with --build-dir: layout)",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help=f"Where trimmed/stripped alignments go (default: {DEFAULT_OUTDIR}, or with "
        "--build-dir the directory of the QC table)",
    )
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=None,
        help="Build directory (builds/<id>/<date>). Every path not given explicitly "
        "defaults to its place in the build layout. Needs --config and --locus.",
    )
    parser.add_argument(
        "--locus", default=None, help="Locus whose files to use; required with --build-dir"
    )
    # Explicit output paths. Snakemake has to know exactly which files a rule
    # produces; a rule that writes to a hardcoded directory cannot be tracked,
    # re-run or cleaned, and its outputs quietly accumulate outside the build.
    parser.add_argument(
        "--out-tsv", type=Path, default=None, help="Per-sequence QC table. Overrides --outdir."
    )
    parser.add_argument(
        "--out-png", type=Path, default=None, help="Coverage plot. Overrides --outdir and --plot."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Pathogen config; supplies QC thresholds from "
        "alignment.* when they are not given on the command line.",
    )
    parser.add_argument(
        "--gap-threshold",
        type=float,
        default=None,
        help="Flag sequences more than this fraction gaps "
        f"(default: alignment.max_gap_fraction from --config, else {DEFAULT_GAP_THRESHOLD})",
    )
    parser.add_argument(
        "--identity-threshold",
        type=float,
        default=None,
        help="Flag sequences below this identity to consensus "
        f"(default: alignment.min_identity_to_consensus from --config, else "
        f"{DEFAULT_IDENTITY_THRESHOLD})",
    )
    parser.add_argument(
        "--ambiguity-threshold",
        type=float,
        default=None,
        help="Flag sequences with more than this fraction ambiguous bases "
        f"(default: alignment.max_ambiguity_fraction from --config, else "
        f"{DEFAULT_AMBIGUITY_THRESHOLD})",
    )
    parser.add_argument(
        "--coverage-threshold",
        type=float,
        default=None,
        help="Columns with less than this coverage count as ragged "
        f"(default {DEFAULT_COVERAGE_THRESHOLD})",
    )
    parser.add_argument(
        "--trim-to",
        nargs=2,
        type=int,
        metavar=("START", "END"),
        help="Trim to these 1-based column positions (inclusive)",
    )
    parser.add_argument("--write-trimmed", action="store_true")
    parser.add_argument(
        "--strip-all-gap",
        action="store_true",
        help="Remove columns that are gaps in every sequence",
    )
    parser.add_argument(
        "--write-stripped", action="store_true", help="Save the result of --strip-all-gap"
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=None,
        help="Write a coverage plot PNG (default: alongside the alignment)",
    )


@step(
    "align", "qc", "Inspect an alignment (gaps, ambiguity, identity, stops, labels) before treeing"
)
def qc(args: argparse.Namespace) -> int:
    """Inspect an alignment before spending compute on a tree.

    This is the "look at your alignment" step. On a desktop you would open the
    alignment in AliView or Jalview; on a hosted kernel you probably cannot, so
    this reports the same things numerically and as a coverage plot.

    What it checks:

      * per-sequence gap and ambiguity fractions, with outliers flagged
      * per-column coverage, and suggested trim points for ragged ends
      * identity of each sequence to the consensus (catches misaligned or
        wrong-gene records)
      * in-frame stop codons (catches frameshifts and wrong reading frames)
      * exact duplicate sequences
      * all-gap columns, which shift every downstream column coordinate
      * tip labels against the metadata table, so a stale host group cannot survive
      * host-group composition of what survives

    Thresholds: a flag given on the command line wins, else ``alignment.*`` from
    ``--config``, else the built-in default. The same values drive the flags and
    the lines on the plot.

    Nothing here decides anything for you. It tells you which sequences are worth
    looking at by eye and which columns are mostly gaps. Use ``--trim-to START END
    --write-trimmed`` once you have chosen trim points from the coverage plot.
    Paths: give ``--aln`` (and ``--out-tsv``, ``--out-png``) explicitly, or give
    ``--config`` with ``--build-dir`` and ``--locus`` and every path not given
    explicitly is taken from the build layout.

    Exit code is 1 when the alignment is missing, empty or not aligned, 2 when
    the command line cannot be resolved to paths.
    """
    try:
        _resolve_paths(args)
    except _UsageError as exc:
        _say(f"wildphylo align qc: {exc}", err=True)
        return 2

    if not args.aln.is_file():
        _report_missing_alignment(args.aln)
        return 1

    seqs = read_fasta(args.aln)
    if not seqs:
        _say(f"No sequences read from {args.aln}", err=True)
        return 1
    _resolve_thresholds(args)

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        _say(
            f"ERROR: sequences have different lengths {sorted(lengths)[:5]} — "
            "this file is not aligned.",
            err=True,
        )
        return 1
    n_cols = lengths.pop()
    n_seqs = len(seqs)
    _say(f"{n_seqs} sequences, {n_cols} alignment columns\n")

    cov = [sum(1 for s in seqs.values() if s[i] not in GAP) / n_seqs for i in range(n_cols)]
    good = _report_coverage(cov, args.coverage_threshold)
    _report_all_gap_columns(cov, good)

    df = _sequence_table(seqs, n_cols)
    _flag_sequences(
        df, seqs, args.gap_threshold, args.identity_threshold, args.ambiguity_threshold
    )
    _report_duplicates(seqs)
    if args.metadata and args.metadata.is_file():
        _report_label_check(seqs, args.metadata)
        _report_host_composition(seqs)

    qc_path = args.out_tsv or args.outdir / (args.aln.stem + "_qc.tsv")
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    df.sort_values("identity_to_consensus").to_csv(qc_path, sep="\t", index=False)
    _say(f"per-sequence QC table -> {qc_path}")

    _write_coverage_plot(cov, df, args, n_cols)
    if args.strip_all_gap:
        _strip_all_gap_columns(seqs, cov, n_cols, args)
    if args.trim_to:
        _trim_columns(seqs, args)

    _say("\nNext: if the alignment looks sound, build the tree (the Snakefile's ml_tree rule).")
    return 0

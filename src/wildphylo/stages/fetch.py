"""fetch stage: download records from GenBank for any taxon.

``wildphylo fetch records`` is the only step; its docstring says how it is driven.
"""

from __future__ import annotations

import datetime as dt
import logging
import os
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any

from Bio import Entrez

from wildphylo.cli import step
from wildphylo.config import ConfigError, load
from wildphylo.logfile import log_to_file
from wildphylo.preflight import check_fetch

if TYPE_CHECKING:
    import argparse

RETRY_LIMIT = 5
RETRY_BACKOFF = 3.0  # seconds, doubled each retry
log = logging.getLogger(__name__)


def build_query(cfg: dict[str, Any]) -> str:
    """Assemble the Entrez query from config.

    ``[Organism:exp]`` expands to every strain and subtaxon beneath the node,
    which is what you want and is easy to leave off by hand. Length bounds come
    from config because they are pathogen-specific: 200-20000 suits a
    morbillivirus and would silently truncate a large dsDNA genome.
    """
    f = cfg["fetch"]
    lo = f.get("min_length", 200)
    hi = f.get("max_length", 50000)
    query = f"txid{f['taxid']}[Organism:exp] AND {lo}:{hi}[Sequence Length]"
    extra = (f.get("query_filter") or "").strip()
    if extra:
        # A filter may carry its own boolean ("NOT patent[Title]"); if not,
        # join with AND.
        if extra.upper().startswith(("AND ", "NOT ", "OR ")):
            query = f"{query} {extra}"
        else:
            query = f"{query} AND {extra}"
    return query


def esearch_history(query: str) -> tuple[str, str, int]:
    """Run esearch with usehistory, so results are paged server-side."""
    handle = Entrez.esearch(  # type: ignore[no-untyped-call]  # Bio.Entrez ships no annotations
        db="nuccore", term=query, usehistory="y", retmax=0
    )
    record = Entrez.read(handle)  # type: ignore[no-untyped-call]  # Bio.Entrez ships no annotations
    handle.close()
    return record["WebEnv"], record["QueryKey"], int(record["Count"])


def _efetch_text(webenv: str, query_key: str, start: int, size: int) -> str:
    """Single efetch call; raises ValueError on an empty response."""
    handle = Entrez.efetch(  # type: ignore[no-untyped-call]  # Bio.Entrez ships no annotations
        db="nuccore",
        rettype="gb",
        retmode="text",
        retstart=start,
        retmax=size,
        webenv=webenv,
        query_key=query_key,
    )
    data: str = handle.read()
    handle.close()
    if not data.strip():
        msg = "empty response from NCBI"
        raise ValueError(msg)
    return data


def fetch_batch(webenv: str, query_key: str, start: int, size: int) -> str:
    """Fetch one batch of GenBank records, retrying on transient NCBI failures."""
    delay = RETRY_BACKOFF
    for attempt in range(1, RETRY_LIMIT + 1):
        try:
            data = _efetch_text(webenv, query_key, start, size)
        except Exception as exc:
            if attempt == RETRY_LIMIT:
                raise
            log.warning(
                "batch at %d failed (attempt %d/%d): %s -- retrying in %.0fs",
                start,
                attempt,
                RETRY_LIMIT,
                exc,
                delay,
            )
            time.sleep(delay)
            delay *= 2
        else:
            return data
    msg = "unreachable"
    raise RuntimeError(msg)


def previous_accession_count(out_path: Path, pathogen: str) -> int | None:
    """Count the accessions in the previous build of this pathogen.

    Lets ``check_fetch`` spot a dataset that shrank. Returns None on a first
    run or outside the ``builds/`` layout.
    """
    try:
        raw = out_path.resolve().parent  # .../builds/<id>/<date>/raw
        build_dir = raw.parent
        pathogen_dir = build_dir.parent
        if pathogen_dir.name != pathogen:
            return None
        earlier = sorted(
            (d for d in pathogen_dir.iterdir() if d.is_dir() and d.name < build_dir.name),
            reverse=True,
        )
        for d in earlier:
            for acc in sorted(d.glob("raw/*.acc")):
                with acc.open() as fh:
                    return sum(1 for line in fh if line.strip())
    except Exception:  # noqa: BLE001
        return None
    return None


def add_arguments_records(parser: argparse.ArgumentParser) -> None:
    """Register the flags of ``wildphylo fetch records``."""
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out", type=Path, help="GenBank flat file to write.")
    parser.add_argument("--acc", type=Path, help="Accession list to write.")
    parser.add_argument(
        "--outdir", type=Path, help="Instead of --out/--acc: write <id>_<date>.gb/.acc here."
    )
    parser.add_argument("--email", default=None, help="NCBI email (default: $NCBI_EMAIL).")
    parser.add_argument("--api-key", default=None, help="NCBI API key (default: $NCBI_API_KEY).")
    parser.add_argument("--query", help="Override the config-derived query entirely.")
    parser.add_argument("--query-filter", help="Override fetch.query_filter for this run only.")
    parser.add_argument("--taxid", type=int, help="Override fetch.taxid for this run only.")
    parser.add_argument("--log", type=Path)
    parser.add_argument(
        "--force", action="store_true", help="Re-download even if the output already exists."
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Report how many records match, then exit."
    )


@step("fetch", "records", "Download GenBank records for the configured taxon")
def records(args: argparse.Namespace) -> int:
    """Download the records matching the pathogen config's query and list their accessions.

    Config-driven. The Entrez query is built from ``fetch.taxid`` and
    ``fetch.query_filter``, so switching pathogens means switching config files
    rather than editing code.

        wildphylo fetch records --config config/pathogen/cdv.yaml --out cdv.gb --acc cdv.acc

        # how many records match, without downloading
        wildphylo fetch records --config config/pathogen/btv.yaml --dry-run

    Email is required by NCBI. Pass ``--email`` or set ``NCBI_EMAIL``. An API key
    (free, from https://www.ncbi.nlm.nih.gov/account/) raises the rate limit from
    3/s to 10/s; set ``NCBI_API_KEY`` or pass ``--api-key``. Neither belongs in a
    config file, which is why both are read from the environment.

    Exit codes: 0 success, 1 nothing matched or the download was short, 2 bad usage.
    """
    try:
        cfg = load(args.config)
    except ConfigError as e:
        log.error("%s", e)  # noqa: TRY400 -- a bad config is user error; no traceback
        return 2

    if args.taxid:
        cfg["fetch"]["taxid"] = args.taxid
    if args.query_filter is not None:
        cfg["fetch"]["query_filter"] = args.query_filter

    stamp = dt.date.today().strftime("%Y%m%d")
    pathogen = cfg["id"]
    if args.out:
        gb_path: Path = args.out
        acc_path: Path = args.acc or gb_path.with_suffix(".acc")
    else:
        outdir: Path = args.outdir or Path("data/raw")
        gb_path = outdir / f"{pathogen}_{stamp}.gb"
        acc_path = args.acc or outdir / f"{pathogen}_{stamp}.acc"

    log_path: Path = args.log or Path("logs") / f"01_fetch_{pathogen}_{stamp}.log"
    with log_to_file(log_path):
        return _download(args, cfg, gb_path, acc_path)


def _download(args: argparse.Namespace, cfg: dict[str, Any], gb_path: Path, acc_path: Path) -> int:
    """Query NCBI, check the match count, then write the flat file and accession list."""
    pathogen = cfg["id"]
    for w in cfg.get("_warnings", []):
        log.warning("config: %s", w)

    email = args.email or os.environ.get("NCBI_EMAIL")
    api_key = args.api_key or os.environ.get("NCBI_API_KEY")
    if not email:
        log.error(
            "NCBI requires an email address. Pass --email or set NCBI_EMAIL. "
            "It is read from the environment rather than the config because it "
            "is a personal credential, not an analysis parameter."
        )
        return 2

    Entrez.email = email  # type: ignore[assignment]  # Bio.Entrez types the module global as None
    if api_key:
        Entrez.api_key = api_key  # type: ignore[assignment]  # Bio.Entrez types the module global as None
        log.info("using NCBI API key (10 req/s)")
    else:
        log.info("no API key -- limited to 3 req/s. Registering is free.")

    query = args.query or build_query(cfg)
    log.info("pathogen: %s (%s)", cfg["name"], pathogen)
    log.info("query: %s", query)

    if gb_path.exists() and not args.force:
        log.info("%s already exists -- nothing to do (use --force)", gb_path)
        return 0

    webenv, query_key, total = esearch_history(query)
    log.info("matched %d records", total)

    checks = check_fetch(cfg, total, previous_accession_count(gb_path, pathogen))
    log.info("%s", checks.render())
    if not checks.ok:
        log.error("fetch checks failed; nothing downloaded")
        return 1

    if args.dry_run:
        log.info("dry run -- exiting without download")
        return 0

    batch = cfg["fetch"].get("batch_size", 200)
    gb_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with gb_path.open("w") as out:
        for start in range(0, total, batch):
            out.write(fetch_batch(webenv, query_key, start, batch))
            log.info("  %d / %d records", min(start + batch, total), total)
            time.sleep(0.12 if api_key else 0.40)

    with gb_path.open() as fh:
        accs = [line.split()[1] for line in fh if line.startswith("VERSION")]
    acc_path.parent.mkdir(parents=True, exist_ok=True)
    acc_path.write_text("\n".join(accs) + "\n")

    log.info("wrote %s (%.1f MB)", gb_path, gb_path.stat().st_size / 1e6)
    log.info("wrote %s (%d accessions)", acc_path, len(accs))
    log.info("done in %.1f s", time.time() - t0)

    if len(accs) != total:
        # A short download is an error, not a warning. Every downstream count,
        # and the build-to-build comparison, would otherwise be computed on a
        # truncated dataset without anyone noticing.
        log.error(
            "expected %d records but parsed %d accessions. The download is "
            "incomplete -- rerun with --force rather than proceeding.",
            total,
            len(accs),
        )
        return 1
    return 0

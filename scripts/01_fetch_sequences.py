#!/usr/bin/env python3
"""
01_fetch_sequences.py — download records from GenBank for any taxon.

Config-driven. The Entrez query is built from `fetch.taxid` and
`fetch.query_filter`, so switching pathogens means switching config files
rather than editing this script.

    python scripts/01_fetch_sequences.py --config config/pathogen/cdv.yaml \\
        --out data/raw/cdv.gb --acc data/raw/cdv.acc

    # how many records match, without downloading
    python scripts/01_fetch_sequences.py --config config/pathogen/btv.yaml --dry-run

Email is required by NCBI. Pass --email or set NCBI_EMAIL. An API key (free,
from https://www.ncbi.nlm.nih.gov/account/) raises the rate limit from 3/s to
10/s; set NCBI_API_KEY or pass --api-key. Neither belongs in a config file,
which is why both are read from the environment.

Exit codes: 0 success, 1 nothing matched or the download was short, 2 bad usage.
"""

import argparse
import datetime as dt
import logging
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config import load, ConfigError  # noqa: E402
from lib.preflight import check_fetch      # noqa: E402

try:
    from Bio import Entrez
except ImportError as exc:                  # pragma: no cover
    raise ImportError(
        "Biopython is required.\n"
        "  terminal: pip install biopython\n"
        "  notebook: %pip install biopython   (then restart the kernel)"
    ) from exc

RETRY_LIMIT = 5
RETRY_BACKOFF = 3.0       # seconds, doubled each retry


def setup_logging(log_path: Path | None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stdout)]
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        handlers=handlers, force=True)


def build_query(cfg: dict) -> str:
    """
    Assemble the Entrez query from config.

    [Organism:exp] expands to every strain and subtaxon beneath the node, which
    is what you want and is easy to leave off by hand. Length bounds come from
    config because they are pathogen-specific: 200-20000 suits a morbillivirus
    and would silently truncate a large dsDNA genome.
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
    """esearch with usehistory, so results are paged server-side."""
    handle = Entrez.esearch(db="nuccore", term=query, usehistory="y", retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return record["WebEnv"], record["QueryKey"], int(record["Count"])


def fetch_batch(webenv: str, query_key: str, start: int, size: int) -> str:
    """One batch of GenBank records, retrying on transient NCBI failures."""
    delay = RETRY_BACKOFF
    for attempt in range(1, RETRY_LIMIT + 1):
        try:
            handle = Entrez.efetch(
                db="nuccore", rettype="gb", retmode="text",
                retstart=start, retmax=size,
                webenv=webenv, query_key=query_key)
            data = handle.read()
            handle.close()
            if not data.strip():
                raise ValueError("empty response from NCBI")
            return data
        except Exception as exc:               # noqa: BLE001 — NCBI throws many types
            if attempt == RETRY_LIMIT:
                raise
            logging.warning(
                "batch at %d failed (attempt %d/%d): %s — retrying in %.0fs",
                start, attempt, RETRY_LIMIT, exc, delay)
            time.sleep(delay)
            delay *= 2
    raise RuntimeError("unreachable")


def previous_accession_count(out_path: Path, pathogen: str) -> int | None:
    """
    How many accessions the previous build of this pathogen had, so
    check_fetch can spot a dataset that shrank. Returns None on a first run or
    outside the builds/ layout.
    """
    try:
        raw = out_path.resolve().parent          # .../builds/<id>/<date>/raw
        build_dir = raw.parent
        pathogen_dir = build_dir.parent
        if pathogen_dir.name != pathogen:
            return None
        earlier = sorted((d for d in pathogen_dir.iterdir()
                          if d.is_dir() and d.name < build_dir.name), reverse=True)
        for d in earlier:
            for acc in sorted(d.glob("raw/*.acc")):
                return sum(1 for line in acc.open() if line.strip())
    except Exception:                            # noqa: BLE001
        return None
    return None


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--out", type=Path, help="GenBank flat file to write.")
    ap.add_argument("--acc", type=Path, help="Accession list to write.")
    ap.add_argument("--outdir", type=Path,
                    help="Instead of --out/--acc: write <id>_<date>.gb/.acc here.")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    ap.add_argument("--query", help="Override the config-derived query entirely.")
    ap.add_argument("--query-filter",
                    help="Override fetch.query_filter for this run only.")
    ap.add_argument("--taxid", type=int, help="Override fetch.taxid for this run only.")
    ap.add_argument("--log", type=Path)
    ap.add_argument("--force", action="store_true",
                    help="Re-download even if the output already exists.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report how many records match, then exit.")
    a = ap.parse_args()

    try:
        cfg = load(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2

    if a.taxid:
        cfg["fetch"]["taxid"] = a.taxid
    if a.query_filter is not None:
        cfg["fetch"]["query_filter"] = a.query_filter

    stamp = dt.date.today().strftime("%Y%m%d")
    pathogen = cfg["id"]
    if a.out:
        gb_path = a.out
        acc_path = a.acc or gb_path.with_suffix(".acc")
    else:
        outdir = a.outdir or Path("data/raw")
        gb_path = outdir / f"{pathogen}_{stamp}.gb"
        acc_path = a.acc or outdir / f"{pathogen}_{stamp}.acc"

    setup_logging(a.log or Path("logs") / f"01_fetch_{pathogen}_{stamp}.log")
    for w in cfg.get("_warnings", []):
        logging.warning("config: %s", w)

    if not a.email:
        logging.error(
            "NCBI requires an email address. Pass --email or set NCBI_EMAIL. "
            "It is read from the environment rather than the config because it "
            "is a personal credential, not an analysis parameter.")
        return 2

    Entrez.email = a.email
    if a.api_key:
        Entrez.api_key = a.api_key
        logging.info("using NCBI API key (10 req/s)")
    else:
        logging.info("no API key — limited to 3 req/s. Registering is free.")

    query = a.query or build_query(cfg)
    logging.info("pathogen: %s (%s)", cfg["name"], pathogen)
    logging.info("query: %s", query)

    if gb_path.exists() and not a.force:
        logging.info("%s already exists — nothing to do (use --force)", gb_path)
        return 0

    webenv, query_key, total = esearch_history(query)
    logging.info("matched %d records", total)

    checks = check_fetch(cfg, total, previous_accession_count(gb_path, pathogen))
    print(checks.render())
    if not checks.ok:
        logging.error("fetch checks failed; nothing downloaded")
        return 1

    if a.dry_run:
        logging.info("dry run — exiting without download")
        return 0

    batch = cfg["fetch"].get("batch_size", 200)
    gb_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with gb_path.open("w") as out:
        for start in range(0, total, batch):
            out.write(fetch_batch(webenv, query_key, start, batch))
            logging.info("  %d / %d records", min(start + batch, total), total)
            time.sleep(0.12 if a.api_key else 0.40)

    accs = [line.split()[1] for line in gb_path.open()
            if line.startswith("VERSION")]
    acc_path.parent.mkdir(parents=True, exist_ok=True)
    acc_path.write_text("\n".join(accs) + "\n")

    logging.info("wrote %s (%.1f MB)", gb_path, gb_path.stat().st_size / 1e6)
    logging.info("wrote %s (%d accessions)", acc_path, len(accs))
    logging.info("done in %.1f s", time.time() - t0)

    if len(accs) != total:
        # A short download is an error, not a warning. Every downstream count,
        # and the build-to-build comparison, would otherwise be computed on a
        # truncated dataset without anyone noticing.
        logging.error(
            "expected %d records but parsed %d accessions. The download is "
            "incomplete — rerun with --force rather than proceeding.",
            total, len(accs))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

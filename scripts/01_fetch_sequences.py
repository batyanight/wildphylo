#!/usr/bin/env python3
"""
01_fetch_sequences.py — pull all Canine distemper virus records from GenBank.

Downloads in batches, caches raw GenBank flatfiles to data/raw/, and is safe to
re-run: already-downloaded batches are skipped unless --force is given.

Usage
-----
    python scripts/01_fetch_sequences.py --email you@example.com
    python scripts/01_fetch_sequences.py --email you@example.com --api-key $NCBI_API_KEY

Register for an NCBI API key (free, 30 seconds) at
https://www.ncbi.nlm.nih.gov/account/ — it raises the rate limit from 3 to 10
requests/second and makes this roughly 3x faster.

Outputs
-------
    data/raw/cdv_YYYYMMDD.gb        concatenated GenBank flatfile
    data/raw/cdv_YYYYMMDD.acc       accession list (for provenance)
    logs/01_fetch_YYYYMMDD.log      query, counts, timings

Provenance note: GenBank changes. The dated filename and accession list are what
make your analysis reproducible later. Do not overwrite old pulls.
"""

import argparse
import datetime as dt
import logging
import os
import sys
import time
from pathlib import Path

try:
    from Bio import Entrez
except ImportError as exc:                       # noqa: F401
    # Raise rather than sys.exit: this module gets imported by the notebook,
    # and SystemExit inside a notebook produces a confusing traceback.
    raise ImportError(
        "Biopython is required.\n"
        "  terminal: pip install biopython\n"
        "  notebook: %pip install biopython   (then restart the kernel)"
    ) from exc

# --- CDV in NCBI Taxonomy ------------------------------------------------
# 11232 = Canine distemper virus, now formally Morbillivirus canis.
# [Organism:exp] expands to include all strains/subtaxa beneath the node.
# Verify at https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=11232
DEFAULT_QUERY = (
    'txid11232[Organism:exp] '
    'AND 200:20000[Sequence Length] '
    'NOT patent[Properties]'
)

BATCH_SIZE = 200          # records per efetch call
RETRY_LIMIT = 5
RETRY_BACKOFF = 3.0       # seconds, doubled each retry


def setup_logging(log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )


def esearch_history(query: str) -> tuple[str, str, int]:
    """Run esearch with usehistory so we can page through results server-side."""
    handle = Entrez.esearch(db="nuccore", term=query, usehistory="y", retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return record["WebEnv"], record["QueryKey"], int(record["Count"])


def fetch_batch(webenv: str, query_key: str, start: int, size: int) -> str:
    """Fetch one batch of GenBank records, retrying on transient NCBI failures."""
    delay = RETRY_BACKOFF
    for attempt in range(1, RETRY_LIMIT + 1):
        try:
            handle = Entrez.efetch(
                db="nuccore",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=size,
                webenv=webenv,
                query_key=query_key,
            )
            data = handle.read()
            handle.close()
            if not data.strip():
                raise ValueError("empty response from NCBI")
            return data
        except Exception as exc:                       # noqa: BLE001 - NCBI throws many types
            if attempt == RETRY_LIMIT:
                raise
            logging.warning(
                "batch at %d failed (attempt %d/%d): %s — retrying in %.0fs",
                start, attempt, RETRY_LIMIT, exc, delay,
            )
            time.sleep(delay)
            delay *= 2
    raise RuntimeError("unreachable")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", required=True,
                    help="Your email. NCBI requires this and will contact you before blocking.")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"),
                    help="NCBI API key (or set NCBI_API_KEY env var).")
    ap.add_argument("--query", default=DEFAULT_QUERY,
                    help="Entrez query. Default pulls all CDV nucleotide records.")
    ap.add_argument("--outdir", default="data/raw", type=Path)
    ap.add_argument("--logdir", default="logs", type=Path)
    ap.add_argument("--force", action="store_true",
                    help="Re-download even if today's file already exists.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report how many records match, then exit without downloading.")
    args = ap.parse_args()

    stamp = dt.date.today().strftime("%Y%m%d")
    setup_logging(args.logdir / f"01_fetch_{stamp}.log")

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
        logging.info("using NCBI API key (10 req/s)")
    else:
        logging.info("no API key — limited to 3 req/s. Consider registering for one.")

    args.outdir.mkdir(parents=True, exist_ok=True)
    gb_path = args.outdir / f"cdv_{stamp}.gb"
    acc_path = args.outdir / f"cdv_{stamp}.acc"

    if gb_path.exists() and not args.force:
        logging.info("%s already exists — nothing to do (use --force to re-download)", gb_path)
        return 0

    logging.info("query: %s", args.query)
    webenv, query_key, total = esearch_history(args.query)
    logging.info("matched %d records", total)

    if total == 0:
        logging.error("no records matched. Check the query and the taxid.")
        return 1

    if args.dry_run:
        logging.info("dry run — exiting without download")
        return 0

    t0 = time.time()
    written = 0
    with gb_path.open("w") as out:
        for start in range(0, total, BATCH_SIZE):
            data = fetch_batch(webenv, query_key, start, BATCH_SIZE)
            out.write(data)
            written += data.count("\nLOCUS       ") + data.startswith("LOCUS")
            logging.info("  %d / %d records", min(start + BATCH_SIZE, total), total)
            time.sleep(0.12 if args.api_key else 0.40)

    # accession list, for provenance
    accs = []
    with gb_path.open() as fh:
        for line in fh:
            if line.startswith("VERSION"):
                accs.append(line.split()[1])
    acc_path.write_text("\n".join(accs) + "\n")

    logging.info("wrote %s (%.1f MB)", gb_path, gb_path.stat().st_size / 1e6)
    logging.info("wrote %s (%d accessions)", acc_path, len(accs))
    logging.info("done in %.1f s", time.time() - t0)

    if len(accs) != total:
        logging.warning(
            "expected %d records but parsed %d accessions — investigate before proceeding",
            total, len(accs),
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

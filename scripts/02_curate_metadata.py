#!/usr/bin/env python3
"""
02_curate_metadata.py — parse GenBank records into a curated metadata table.

Does four things:
  1. Extracts accession, host, country, collection date, strain, gene content
  2. Normalizes host names to canonical taxa and functional groups
  3. Flags vaccine and vaccine-derived strains for exclusion
  4. Parses collection dates to decimal years for tip-dated phylogenetics

Everything excluded is logged with a reason. Everything ambiguous is routed to
needs_review.tsv for human adjudication — that file is your job, not the script's.

Usage
-----
    python scripts/02_curate_metadata.py --gb data/raw/cdv_20260819.gb

Outputs
-------
    data/interim/metadata_all.tsv       every record, with flags, nothing dropped
    data/interim/needs_review.tsv       ambiguous host or date — REVIEW BY HAND
    data/interim/exclusions.tsv         what was excluded and why
    data/processed/metadata_clean.tsv   analysis-ready records
    data/processed/sequences_H.fasta    H (hemagglutinin) gene sequences
    data/processed/sequences_all.fasta  all retained sequences
    logs/02_curate_YYYYMMDD.log

Note on the H gene: CDV hemagglutinin is annotated inconsistently across GenBank
as "H", "HA", "hemagglutinin" or "haemagglutinin", and sometimes not at all on
whole-genome records. The extractor checks feature qualifiers first, then falls
back to a length heuristic on records whose description mentions H. Records it
can't resolve go to needs_review rather than being silently dropped.
"""

import argparse
import datetime as dt
import logging
import re
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError as exc:                       # noqa: F401
    # Raise rather than sys.exit: this module gets imported by the notebook,
    # and SystemExit inside a notebook produces a confusing traceback.
    raise ImportError(
        "Biopython is required.\n"
        "  terminal: pip install biopython\n"
        "  notebook: %pip install biopython   (then restart the kernel)"
    ) from exc

try:
    import pandas as pd
except ImportError as exc:                       # noqa: F401
    # Raise rather than sys.exit: this module gets imported by the notebook,
    # and SystemExit inside a notebook produces a confusing traceback.
    raise ImportError(
        "pandas is required.\n"
        "  terminal: pip install pandas\n"
        "  notebook: %pip install pandas   (then restart the kernel)"
    ) from exc


# CDV H gene is ~1824 nt (607 aa). Allow generous slack for partials.
H_FULL_MIN, H_FULL_MAX = 1750, 1900
H_PARTIAL_MIN = 400          # shorter than this is too little signal to be useful

H_NAMES = {"h", "ha", "hemagglutinin", "haemagglutinin", "hemaglutinin",
           "h protein", "attachment protein", "hemagglutinin protein"}
F_NAMES = {"f", "fusion", "fusion protein", "f protein"}
N_NAMES = {"n", "np", "nucleocapsid", "nucleoprotein", "nucleocapsid protein"}

MONTHS = {m: i for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun",
     "jul", "aug", "sep", "oct", "nov", "dec"], start=1)}


# ----------------------------------------------------------------------------
# config loading
# ----------------------------------------------------------------------------

def load_vaccine_patterns(path: Path) -> list[str]:
    pats = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            pats.append(line.lower())
    return pats


def load_host_table(path: Path) -> list[tuple[str, str, str]]:
    """Ordered list of (pattern, canonical_host, host_group). First match wins."""
    rows = []
    for line in path.read_text().splitlines():
        line = line.rstrip()
        if not line or line.lstrip().startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 3:
            logging.warning("skipping malformed host_groups line: %r", line)
            continue
        rows.append((parts[0].strip().lower(), parts[1].strip(), parts[2].strip()))
    return rows


# ----------------------------------------------------------------------------
# field extraction
# ----------------------------------------------------------------------------

def first_qualifier(feature, keys) -> str:
    for key in keys:
        if key in feature.qualifiers:
            return feature.qualifiers[key][0]
    return ""


def extract_source_fields(record) -> dict:
    out = {"host": "", "country": "", "collection_date": "",
           "strain": "", "isolate": "", "lat_lon": ""}
    for feat in record.features:
        if feat.type == "source":
            out["host"] = first_qualifier(feat, ["host", "specific_host", "isolation_source"])
            out["country"] = first_qualifier(feat, ["geo_loc_name", "country"])
            out["collection_date"] = first_qualifier(feat, ["collection_date"])
            out["strain"] = first_qualifier(feat, ["strain"])
            out["isolate"] = first_qualifier(feat, ["isolate"])
            out["lat_lon"] = first_qualifier(feat, ["lat_lon"])
            break
    return out


def gene_label(feature) -> str:
    label = first_qualifier(feature, ["gene", "product", "note"]).strip().lower()
    return re.sub(r"\s+", " ", label)


def extract_gene_seq(record, name_set) -> tuple[str, str]:
    """Return (sequence, how_found) for the first CDS/gene matching name_set."""
    for feat in record.features:
        if feat.type not in ("CDS", "gene", "mat_peptide"):
            continue
        label = gene_label(feat)
        if not label:
            continue
        if label in name_set or any(label.startswith(n + " ") for n in name_set):
            try:
                seq = str(feat.extract(record.seq))
            except Exception:                          # noqa: BLE001
                continue
            if seq:
                return seq, "annotation"
    return "", ""


def extract_H(record) -> tuple[str, str]:
    """H gene by annotation, falling back to a whole-record length heuristic."""
    seq, how = extract_gene_seq(record, H_NAMES)
    if seq:
        return seq, how

    desc = record.description.lower()
    mentions_H = any(k in desc for k in
                     ("hemagglutinin", "haemagglutinin", " h gene", "h protein", "(h)"))
    if mentions_H and H_FULL_MIN <= len(record.seq) <= H_FULL_MAX:
        return str(record.seq), "length_heuristic_full"
    if mentions_H and H_PARTIAL_MIN <= len(record.seq) < H_FULL_MIN:
        return str(record.seq), "length_heuristic_partial"
    return "", ""


# ----------------------------------------------------------------------------
# normalization
# ----------------------------------------------------------------------------

def normalize_host(raw: str, host_table) -> tuple[str, str, bool]:
    """Return (canonical_host, host_group, is_ambiguous)."""
    if not raw or not raw.strip():
        return "", "unknown", True
    low = raw.strip().lower()
    for pattern, canonical, group in host_table:
        if pattern in low:
            return canonical, group, False
    return "", "unknown", True


def parse_collection_date(raw: str) -> tuple[str, float | None, str]:
    """
    Return (iso_date_or_partial, decimal_year, precision).
    precision is one of: day, month, year, range, none.

    Decimal year uses the midpoint of the known interval, which is the honest
    representation for month- and year-only dates in tip-dated phylogenetics.
    """
    if not raw or not raw.strip():
        return "", None, "none"
    s = raw.strip()

    # "2011/2013" or "2011-01-01/2011-12-31" — a range. Use the midpoint.
    if "/" in s:
        parts = [p.strip() for p in s.split("/") if p.strip()]
        vals = [parse_collection_date(p)[1] for p in parts]
        vals = [v for v in vals if v is not None]
        if vals:
            return s, sum(vals) / len(vals), "range"
        return s, None, "none"

    # ISO: 2013-08-14 or 2013-08
    m = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})", s)
    if m:
        y, mo, d = map(int, m.groups())
        return s, to_decimal_year(y, mo, d), "day"
    m = re.fullmatch(r"(\d{4})-(\d{2})", s)
    if m:
        y, mo = map(int, m.groups())
        return s, to_decimal_year(y, mo, None), "month"

    # GenBank style: 14-Aug-2013 / Aug-2013
    m = re.fullmatch(r"(\d{1,2})-([A-Za-z]{3})-(\d{4})", s)
    if m:
        d, mon, y = int(m.group(1)), m.group(2).lower(), int(m.group(3))
        if mon in MONTHS:
            return f"{y:04d}-{MONTHS[mon]:02d}-{d:02d}", to_decimal_year(y, MONTHS[mon], d), "day"
    m = re.fullmatch(r"([A-Za-z]{3})-(\d{4})", s)
    if m:
        mon, y = m.group(1).lower(), int(m.group(2))
        if mon in MONTHS:
            return f"{y:04d}-{MONTHS[mon]:02d}", to_decimal_year(y, MONTHS[mon], None), "month"

    # bare year
    m = re.fullmatch(r"(\d{4})", s)
    if m:
        y = int(m.group(1))
        return s, to_decimal_year(y, None, None), "year"

    return s, None, "none"


def to_decimal_year(year: int, month: int | None, day: int | None) -> float:
    """Decimal year; midpoint of the known interval when month/day are missing."""
    start = dt.date(year, 1, 1)
    days_in_year = (dt.date(year + 1, 1, 1) - start).days
    if month is None:
        return year + 0.5
    if day is None:
        first = dt.date(year, month, 1)
        nxt = dt.date(year + 1, 1, 1) if month == 12 else dt.date(year, month + 1, 1)
        mid = first + (nxt - first) / 2
        return year + (mid - start).days / days_in_year
    return year + (dt.date(year, month, day) - start).days / days_in_year


def is_vaccine(record, meta, patterns) -> tuple[bool, str]:
    haystack = " | ".join([
        record.description, record.annotations.get("organism", ""),
        meta.get("strain", ""), meta.get("isolate", ""), meta.get("host", ""),
    ]).lower()
    for pat in patterns:
        if pat in haystack:
            return True, pat
    return False, ""


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gb", required=True, type=Path, help="GenBank file from step 01")
    ap.add_argument("--config", default=Path("config"), type=Path)
    ap.add_argument("--interim", default=Path("data/interim"), type=Path)
    ap.add_argument("--processed", default=Path("data/processed"), type=Path)
    ap.add_argument("--logdir", default=Path("logs"), type=Path)
    ap.add_argument("--min-length", type=int, default=H_PARTIAL_MIN,
                    help="Drop H sequences shorter than this (default %(default)s nt)")
    ap.add_argument("--keep-vaccines", action="store_true",
                    help="Retain vaccine strains (only for testing the filter)")
    args = ap.parse_args()

    stamp = dt.date.today().strftime("%Y%m%d")
    args.logdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        handlers=[logging.FileHandler(args.logdir / f"02_curate_{stamp}.log"),
                  logging.StreamHandler(sys.stdout)],
    )

    vaccine_patterns = load_vaccine_patterns(args.config / "vaccine_strains.txt")
    host_table = load_host_table(args.config / "host_groups.tsv")
    logging.info("loaded %d vaccine patterns, %d host patterns",
                 len(vaccine_patterns), len(host_table))

    rows, h_seqs, all_seqs = [], {}, {}

    for record in SeqIO.parse(str(args.gb), "genbank"):
        meta = extract_source_fields(record)
        canonical, group, host_ambiguous = normalize_host(meta["host"], host_table)
        iso_date, dec_year, precision = parse_collection_date(meta["collection_date"])
        vaccine, vaccine_hit = is_vaccine(record, meta, vaccine_patterns)
        h_seq, h_how = extract_H(record)
        f_seq, _ = extract_gene_seq(record, F_NAMES)
        n_seq, _ = extract_gene_seq(record, N_NAMES)

        rows.append({
            "accession": record.id,
            "description": record.description,
            "length": len(record.seq),
            "host_raw": meta["host"],
            "host_canonical": canonical,
            "host_group": group,
            "host_ambiguous": host_ambiguous,
            "country_raw": meta["country"],
            "country": meta["country"].split(":")[0].strip() if meta["country"] else "",
            "collection_date_raw": meta["collection_date"],
            "collection_date": iso_date,
            "decimal_year": dec_year,
            "date_precision": precision,
            "strain": meta["strain"],
            "isolate": meta["isolate"],
            "lat_lon": meta["lat_lon"],
            "is_vaccine": vaccine,
            "vaccine_pattern": vaccine_hit,
            "has_H": bool(h_seq),
            "H_length": len(h_seq),
            "H_source": h_how,
            "has_F": bool(f_seq),
            "has_N": bool(n_seq),
            "is_complete_genome": "complete genome" in record.description.lower(),
        })
        if h_seq:
            h_seqs[record.id] = h_seq
        all_seqs[record.id] = str(record.seq)

    if not rows:
        logging.error("no records parsed from %s", args.gb)
        return 1

    df = pd.DataFrame(rows)
    args.interim.mkdir(parents=True, exist_ok=True)
    args.processed.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.interim / "metadata_all.tsv", sep="\t", index=False)
    logging.info("parsed %d records", len(df))

    # ---- exclusions -------------------------------------------------------
    df["exclude_reason"] = ""

    def mark(mask, reason):
        newly = mask & (df["exclude_reason"] == "")
        df.loc[newly, "exclude_reason"] = reason
        return int(newly.sum())

    if not args.keep_vaccines:
        n = mark(df["is_vaccine"], "vaccine_or_vaccine_derived")
        logging.info("excluded %d vaccine / vaccine-derived records", n)
    else:
        logging.warning("--keep-vaccines set: vaccine strains RETAINED (test mode only)")

    n = mark(~df["has_H"], "no_H_gene_sequence")
    logging.info("excluded %d records without an identifiable H sequence", n)

    n = mark(df["H_length"] < args.min_length, f"H_shorter_than_{args.min_length}nt")
    logging.info("excluded %d records with H shorter than %d nt", n, args.min_length)

    # Catches whole-genome records where the H annotation spans the entire record,
    # or mis-annotated features. These are recoverable by hand, not junk.
    n = mark(df["H_length"] > 2100, "H_length_implausible_check_annotation")
    logging.info("excluded %d records with implausibly long H (annotation problem)", n)

    n = mark(df["decimal_year"].isna(), "no_parseable_collection_date")
    logging.info("excluded %d records without a usable collection date", n)

    n = mark(df["host_group"] == "unknown", "host_unresolved")
    logging.info("excluded %d records with unresolved host", n)

    excluded = df[df["exclude_reason"] != ""]
    clean = df[df["exclude_reason"] == ""].copy()
    excluded.to_csv(args.interim / "exclusions.tsv", sep="\t", index=False)

    # ---- needs review -----------------------------------------------------
    review = df[
        (df["host_ambiguous"] & (df["host_raw"].str.strip() != ""))
        | ((df["decimal_year"].isna()) & (df["collection_date_raw"].str.strip() != ""))
        | (df["H_source"] == "length_heuristic_partial")
        | (df["H_length"] > 2100)
    ]
    review.to_csv(args.interim / "needs_review.tsv", sep="\t", index=False)

    # ---- outputs ----------------------------------------------------------
    clean.to_csv(args.processed / "metadata_clean.tsv", sep="\t", index=False)

    def write_fasta(path: Path, seqs: dict, keep: set) -> int:
        written = 0
        with path.open("w") as fh:
            for acc, seq in seqs.items():
                if acc in keep:
                    fh.write(f">{acc}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
                    written += 1
        return written

    keep = set(clean["accession"])
    n_h = write_fasta(args.processed / "sequences_H.fasta", h_seqs, keep)
    n_all = write_fasta(args.processed / "sequences_all.fasta", all_seqs, keep)

    # ---- summary ----------------------------------------------------------
    logging.info("-" * 62)
    logging.info("RETAINED: %d of %d records", len(clean), len(df))
    logging.info("H sequences written: %d   (all-locus: %d)", n_h, n_all)
    logging.info("")
    logging.info("By host group:")
    for grp, cnt in clean["host_group"].value_counts().items():
        logging.info("    %-16s %5d", grp, cnt)
    logging.info("")
    logging.info("By date precision:")
    for prec, cnt in clean["date_precision"].value_counts().items():
        logging.info("    %-16s %5d", prec, cnt)
    if not clean.empty:
        logging.info("")
        logging.info("Date range: %.2f – %.2f",
                     clean["decimal_year"].min(), clean["decimal_year"].max())
        logging.info("Countries represented: %d", clean["country"].nunique())
    logging.info("")
    logging.info("NEEDS HUMAN REVIEW: %d records -> %s",
                 len(review), args.interim / "needs_review.tsv")
    logging.info("-" * 62)
    logging.info("")
    logging.info("NEXT: open needs_review.tsv. For each ambiguous host, either add a")
    logging.info("pattern to config/host_groups.tsv or confirm the record should be")
    logging.info("dropped. Then re-run this script. Repeat until the file is empty or")
    logging.info("everything left is genuinely unusable.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

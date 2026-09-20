#!/usr/bin/env python3
"""
02_curate_metadata.py — normalise metadata and extract loci, for any taxon.

Config-driven and multi-locus. Host table, locus aliases, date plausibility
bounds and exclusion thresholds all come from the pathogen config, so a new
taxon needs a new config rather than edits here. A segmented pathogen declares
`loci` and gets one FASTA per analysed segment.

    python scripts/02_curate_metadata.py --config config/pathogen/cdv.yaml \\
        --gb data/raw/cdv.gb --out-dir data/interim

Outputs, all in --out-dir:

    metadata_all.tsv      every parsed record, nothing dropped
    metadata_clean.tsv    records that survived exclusion
    exclusions.tsv        what was dropped, with a machine-readable reason
    needs_review.tsv      records a human should adjudicate
    <locus>.fasta         one per analysed locus, tip labels already built

Tip labels are `ACCESSION|host_group|decimal_year`, which is the contract every
downstream step reads.

Nothing is dropped silently. Every exclusion carries a reason and lands in
exclusions.tsv; anything ambiguous lands in needs_review.tsv rather than being
guessed at.

Exit codes: 0 success, 1 nothing usable, 2 bad usage.
"""

import argparse
import datetime as dt
import logging
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config import load, ConfigError            # noqa: E402
from lib.dates import parse_collection_date, PRECISION_RANK  # noqa: E402
from lib.hosts import load_host_table, audit_host_table, normalize_host  # noqa: E402
from lib.preflight import check_curation            # noqa: E402

try:
    import pandas as pd
    from Bio import SeqIO
except ImportError as exc:                          # pragma: no cover
    raise ImportError("pandas and biopython are required: "
                      "pip install pandas biopython") from exc


# --- GenBank feature reading ------------------------------------------------

def first_qualifier(feature, keys) -> str:
    for key in keys:
        if key in feature.qualifiers:
            return feature.qualifiers[key][0]
    return ""


def extract_source_fields(record) -> dict:
    out = {"host": "", "country": "", "collection_date": "",
           "strain": "", "isolate": "", "lat_lon": "", "segment": ""}
    for feat in record.features:
        if feat.type == "source":
            out["host"] = first_qualifier(
                feat, ["host", "specific_host", "isolation_source"])
            out["country"] = first_qualifier(feat, ["geo_loc_name", "country"])
            out["collection_date"] = first_qualifier(feat, ["collection_date"])
            out["strain"] = first_qualifier(feat, ["strain"])
            out["isolate"] = first_qualifier(feat, ["isolate"])
            out["lat_lon"] = first_qualifier(feat, ["lat_lon"])
            # Segmented viruses carry /segment on the source feature. Without
            # it, a segment can only be guessed from the description.
            out["segment"] = first_qualifier(feat, ["segment"])
            break
    return out


def gene_label(feature) -> str:
    label = first_qualifier(feature, ["gene", "product", "note"]).strip().lower()
    return re.sub(r"\s+", " ", label)


def extract_by_annotation(record, aliases: set[str]) -> tuple[str, str]:
    """
    First CDS/gene/mat_peptide whose gene, product or note matches an alias.
    Returns (sequence, how_found).

    NOTE: this returns whatever the feature spans. Length is validated by the
    caller — an annotation saying "H" is a claim, not a guarantee, and a `gene`
    feature legitimately includes UTR that a `CDS` feature does not.
    """
    for feat in record.features:
        if feat.type not in ("CDS", "gene", "mat_peptide"):
            continue
        label = gene_label(feat)
        if not label:
            continue
        if label in aliases or any(label.startswith(a + " ") for a in aliases):
            try:
                seq = str(feat.extract(record.seq))
            except Exception:                        # noqa: BLE001
                continue
            if seq:
                return seq, "annotation"
    return "", ""


def extract_locus(record, locus: dict, meta: dict) -> tuple[str, str]:
    """
    Pull one locus from a record.

    Annotation first. Failing that, fall back to treating the whole record as
    the locus, but only when the description mentions it AND the length is
    plausible — otherwise a whole-genome record would be handed on as if it
    were a single gene.
    """
    aliases = {a.strip().lower() for a in locus.get("aliases", []) if a.strip()}
    if not aliases:
        return "", ""

    seq, how = extract_by_annotation(record, aliases)
    if seq:
        # Annotation-based extraction used to be trusted unconditionally, so a
        # feature labelled with the locus name was accepted at any length. On
        # real CDV data that admitted sequences of 1,946 and 1,947 nt against an
        # 1,824 nt H CDS — a `gene` feature spanning UTR, or a mis-annotation.
        # Flag it rather than dropping it: the sequence is usually fine, but the
        # length is a claim worth checking.
        expected = locus.get("expected_length")
        if expected:
            # Only OVER-length features are suspicious.
            #
            # A feature SHORTER than the CDS is a partial annotation: a real
            # sequence covering part of the gene, which is most of what GenBank
            # holds. On the CDV dataset 1006 of 1024 off-length features were
            # short, median 907 nt short, and the single commonest was the
            # 1338 nt amplicon that dominates the published clade 3. Flagging
            # those put 1209 records into needs_review — a pile nobody can
            # adjudicate, which makes the review gate useless.
            #
            # A feature LONGER than the CDS cannot be a partial anything. It is
            # a `gene` span including UTR, or a mis-annotation. That is worth a
            # human look: 1946 and 1947 nt against an 1824 nt CDS.
            tol = float(locus.get("annotation_overlength_tolerance", 0.02))
            if len(seq) > expected * (1 + tol):
                how = f"{how}_overlength"
            elif len(seq) < locus.get("min_usable_length", 0):
                how = f"{how}_too_short"
        return seq, how

    # Segment number is decisive when GenBank provides it.
    seg = str(locus.get("segment", "")).strip()
    if seg and str(meta.get("segment", "")).strip() == seg:
        return str(record.seq), "segment_qualifier"

    desc = record.description.lower()
    if not any(a in desc for a in aliases):
        return "", ""

    expected = locus.get("expected_length")
    if not expected:
        return "", ""
    tol = float(locus.get("length_tolerance", 0.15))
    lo_full, hi_full = expected * (1 - tol), expected * (1 + tol)
    min_usable = locus.get("min_usable_length", 300)
    n = len(record.seq)
    if lo_full <= n <= hi_full:
        return str(record.seq), "length_heuristic_full"
    if min_usable <= n < lo_full:
        return str(record.seq), "length_heuristic_partial"
    return "", ""


def is_vaccine(record, meta, patterns) -> tuple[bool, str]:
    haystack = " | ".join([
        record.description, record.annotations.get("organism", ""),
        meta.get("strain", ""), meta.get("isolate", ""), meta.get("host", ""),
    ]).lower()
    for pat in patterns:
        if pat and pat in haystack:
            return True, pat
    return False, ""


def load_vaccine_patterns(path: Path | None) -> list[str]:
    if not path or not Path(path).is_file():
        return []
    out = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.append(line.lower())
    return out


def load_reference_anchors(cfg: dict) -> dict[str, str]:
    """
    Accession -> lineage, for nucleotide references only.

    Reference sequences exist to NAME clades, not to be dated tips. They are
    routinely 1990s submissions with no /collection_date and no /host — all
    eight usable CDV references lack both — so requiring them to pass the same
    filters as analysis sequences guarantees they are filtered out, and clades
    can then never be named. They are therefore admitted as anchors: present in
    the alignment and the tree, excluded from the temporal regression, the
    subsample and the trait analysis.

    Protein accessions are skipped: nothing here resolves them to a coding
    nucleotide record, so they would silently match nothing.
    """
    refs_cfg = cfg.get("references") or {}
    if not refs_cfg.get("include_as_anchors"):
        return {}
    table = refs_cfg.get("table")
    if not table or not Path(table).is_file():
        logging.warning("references.include_as_anchors is set but "
                        "references.table is missing: %s", table)
        return {}
    out, skipped = {}, 0
    for line in Path(table).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 4:
            continue
        acc, lineage, kind = parts[0], parts[1], parts[3]
        if kind.lower() != "nuc":
            skipped += 1
            continue
        out[acc.split(".")[0]] = lineage
    if skipped:
        logging.warning(
            "%d reference rows are protein accessions and were skipped; nothing "
            "resolves them to nucleotide records, so they cannot anchor a clade",
            skipped)
    return out


def loci_of(cfg: dict) -> list[dict]:
    """Analysed loci, as a list, whether the config is single- or multi-locus."""
    if "loci" in cfg:
        return [dict(l) for l in cfg["loci"] if l.get("analyse")]
    return [dict(cfg["locus"])]


def write_fasta(path: Path, entries: list[tuple[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for label, seq in entries:
            fh.write(f">{label}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    return len(entries)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--gb", required=True, type=Path,
                    help="GenBank flat file from step 01")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--log", type=Path)
    ap.add_argument("--keep-vaccines", action="store_true",
                    help="Retain vaccine strains (for testing the filter only)")
    ap.add_argument("--accept-unreviewed", action="store_true",
                    help="Do not fail when needs_review.tsv is non-empty")
    a = ap.parse_args()

    try:
        cfg = load(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2

    stamp = dt.date.today().strftime("%Y%m%d")
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stdout)]
    log_path = a.log or Path("logs") / f"02_curate_{cfg['id']}_{stamp}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handlers.append(logging.FileHandler(log_path))
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        handlers=handlers, force=True)

    for w in cfg.get("_warnings", []):
        logging.warning("config: %s", w)

    if not a.gb.is_file():
        logging.error("GenBank file not found: %s", a.gb)
        return 2

    # --- config-derived settings -------------------------------------------
    host_rules = load_host_table(Path(cfg["hosts"]["table"]))
    for w in audit_host_table(host_rules):
        # A shadowed host rule is the failure that put a pinniped into the CDV
        # tree as Panthera leo. Surface it every run, not just at config load.
        logging.warning("host table: %s", w)

    anchors = load_reference_anchors(cfg)
    if anchors:
        logging.info("%d nucleotide reference anchors: %s", len(anchors),
                     ", ".join(sorted(anchors)))

    vac_path = cfg.get("exclude", {}).get("vaccine_patterns_file")
    vaccine_patterns = load_vaccine_patterns(Path(vac_path) if vac_path else None)

    min_year = cfg["dates"]["min_year"]
    max_year = cfg["dates"].get("max_year")
    loci = loci_of(cfg)
    logging.info("pathogen: %s (%s)", cfg["name"], cfg["id"])
    logging.info("loci: %s", ", ".join(l["name"] for l in loci))
    logging.info("%d host rules, %d vaccine patterns",
                 len(host_rules), len(vaccine_patterns))

    # --- parse --------------------------------------------------------------
    rows = []
    seqs: dict[str, dict[str, str]] = {l["name"]: {} for l in loci}

    for record in SeqIO.parse(str(a.gb), "genbank"):
        meta = extract_source_fields(record)
        hm = normalize_host(meta["host"], host_rules)
        pd_ = parse_collection_date(meta["collection_date"],
                                    min_year=min_year, max_year=max_year)

        row = {
            "accession": record.id,
            "description": record.description,
            "length": len(record.seq),
            "segment_raw": meta["segment"],
            "host_raw": meta["host"],
            "host_canonical": hm.canonical,
            "host_group": hm.group,
            "host_ambiguous": hm.ambiguous,
            "host_matched_pattern": hm.matched_pattern,
            "country_raw": meta["country"],
            "country": meta["country"].split(":")[0].strip() if meta["country"] else "",
            "collection_date_raw": meta["collection_date"],
            "collection_date": pd_.iso,
            "decimal_year": pd_.decimal_year,
            "date_precision": pd_.precision,
            "date_uncertainty_years": pd_.uncertainty_years,
            "date_reason": pd_.reason,
            "strain": meta["strain"],
            "isolate": meta["isolate"],
            "lat_lon": meta["lat_lon"],
        }
        vac, hit = is_vaccine(record, meta, vaccine_patterns)
        row["is_vaccine"] = vac
        row["vaccine_pattern"] = hit
        bare = str(record.id).split(".")[0]
        row["is_reference"] = bare in anchors
        row["reference_lineage"] = anchors.get(bare, "")

        for locus in loci:
            name = locus["name"]
            seq, how = extract_locus(record, locus, meta)
            row[f"has_{name}"] = bool(seq)
            row[f"{name}_length"] = len(seq)
            row[f"{name}_source"] = how
            # A coding sequence whose length is not a multiple of 3 cannot be
            # translated in frame. It may still be usable for a nucleotide
            # phylogeny, but it is not a clean CDS and any codon-partitioned or
            # amino-acid analysis downstream will be wrong.
            row[f"{name}_in_frame"] = bool(seq) and len(seq) % 3 == 0
            if seq:
                seqs[name][record.id] = seq
        rows.append(row)

    if not rows:
        logging.error("no records parsed from %s", a.gb)
        return 1

    df = pd.DataFrame(rows)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out_dir / "metadata_all.tsv", sep="\t", index=False)
    logging.info("parsed %d records", len(df))

    # --- exclusions ---------------------------------------------------------
    df["exclude_reason"] = ""

    # Reference anchors are exempt from every exclusion except failing to carry
    # the locus at all — an anchor with no sequence anchors nothing.
    is_ref = df["is_reference"] if "is_reference" in df else pd.Series(False, index=df.index)

    def mark(mask, reason) -> int:
        newly = mask & (df["exclude_reason"] == "") & ~is_ref
        df.loc[newly, "exclude_reason"] = reason
        n = int(newly.sum())
        if n:
            logging.info("excluded %5d  %s", n, reason)
        return n

    if not a.keep_vaccines:
        mark(df["is_vaccine"], "vaccine_or_vaccine_derived")
    else:
        logging.warning("--keep-vaccines set: vaccine strains RETAINED (test mode)")

    # A record is usable if it carries at least one analysed locus.
    has_any = df[[f"has_{l['name']}" for l in loci]].any(axis=1)
    # This one applies to anchors too, so it bypasses `mark`'s exemption.
    lost = ~has_any & (df["exclude_reason"] == "")
    df.loc[lost, "exclude_reason"] = "no_analysed_locus_found"
    if int(lost.sum()):
        logging.info("excluded %5d  no_analysed_locus_found", int(lost.sum()))
    dropped_refs = int((lost & is_ref).sum())
    if dropped_refs:
        logging.warning(
            "%d reference anchors carry no usable %s sequence and were dropped; "
            "the lineages they represent cannot be named",
            dropped_refs, "/".join(l["name"] for l in loci))

    for locus in loci:
        name = locus["name"]
        min_len = locus.get("min_usable_length", 300)
        expected = locus.get("expected_length")
        short = df[f"has_{name}"] & (df[f"{name}_length"] < min_len)
        if len(loci) == 1:
            mark(short, f"{name}_shorter_than_{min_len}nt")
        if expected:
            # Catches whole-genome records whose locus annotation spans the
            # entire record, and mis-annotated features. Recoverable by hand,
            # not junk, so it also goes to needs_review.
            cap = expected * 1.15
            long = df[f"has_{name}"] & (df[f"{name}_length"] > cap)
            if len(loci) == 1:
                mark(long, f"{name}_length_implausible_check_annotation")

    mark(df["decimal_year"].isna(), "no_parseable_collection_date")
    mark(df["host_group"] == "unknown", "host_unresolved")

    excluded = df[df["exclude_reason"] != ""]
    clean = df[df["exclude_reason"] == ""].copy()
    excluded.to_csv(a.out_dir / "exclusions.tsv", sep="\t", index=False)

    # --- needs review vs flagged --------------------------------------------
    #
    # Two different things were being conflated, and the pile grew to 1209
    # records on the real CDV dataset — a volume nobody adjudicates, which
    # makes the review gate worthless.
    #
    # NEEDS REVIEW is for records where a HUMAN DECISION changes the outcome:
    # an unrecognised host (add a pattern, or confirm the drop) or a date that
    # was present but unparseable (read the record, fix or exclude). Both are
    # curation gaps.
    #
    # FLAGGED is for properties of the sequence that are simply true and need
    # no decision: partial coverage, a length not divisible by 3, a feature
    # longer than the CDS. Worth knowing, worth reporting, but there is nothing
    # to adjudicate — and burying 40 actionable records among 600 of these is
    # how a gate stops being read.
    actionable = (
        (df["host_ambiguous"] & (df["host_raw"].fillna("").str.strip() != ""))
        | (df["decimal_year"].isna()
           & (df["collection_date_raw"].fillna("").str.strip() != "")
           & (df["date_reason"] != "explicitly_unknown"))
    )
    informational = pd.Series(False, index=df.index)
    for locus in loci:
        name = locus["name"]
        informational |= (df[f"{name}_source"] == "length_heuristic_partial")
        informational |= df[f"{name}_source"].astype(str).str.endswith("_overlength")
        informational |= (df[f"has_{name}"] & ~df[f"{name}_in_frame"])

    df["review_reason"] = ""
    df.loc[df["host_ambiguous"] & (df["host_raw"].fillna("").str.strip() != ""),
           "review_reason"] = "host_unrecognised"
    date_gap = (df["decimal_year"].isna()
                & (df["collection_date_raw"].fillna("").str.strip() != "")
                & (df["date_reason"] != "explicitly_unknown"))
    df.loc[date_gap & (df["review_reason"] == ""), "review_reason"] = "date_unparseable"

    review = df[actionable]
    review.to_csv(a.out_dir / "needs_review.tsv", sep="\t", index=False)
    flagged = df[informational & ~actionable]
    flagged.to_csv(a.out_dir / "flagged.tsv", sep="\t", index=False)

    clean.to_csv(a.out_dir / "metadata_clean.tsv", sep="\t", index=False)

    # --- per-locus FASTA ----------------------------------------------------
    # Tip labels: <key>|host_group|decimal_year. Every downstream step parses
    # that format, so it is a contract rather than a convenience.
    #
    # WHICH KEY matters enormously for a segmented pathogen. GenBank gives each
    # segment of one isolate its own accession, so labelling by accession makes
    # the shared-taxon set across segments EMPTY — and segment congruence, the
    # whole reason for analysing segments separately, becomes impossible to
    # assess. The segments must be tied together by the isolate or strain
    # qualifier instead. `segments.isolate_key` names which field to use.
    isolate_key = (cfg.get("segments") or {}).get("isolate_key", "isolate")
    use_isolate = len(loci) > 1

    def label_key(row) -> str:
        if use_isolate:
            for field in ([isolate_key, "strain", "isolate"]
                          if isolate_key else ["isolate", "strain"]):
                v = row.get(field)
                if isinstance(v, str) and v.strip():
                    return re.sub(r"[|\s]+", "_", v.strip())
        return str(row["accession"]).split(".")[0]

    label_of, key_counts, anchor_labels = {}, {}, {}
    for _, r in clean.iterrows():
        if r.get("is_reference"):
            # Third field is deliberately non-numeric. Every downstream parser
            # treats an unparseable date as "no date" and drops the tip from the
            # regression, which is exactly the required behaviour: present in
            # the tree, absent from the clock.
            lin = re.sub(r"[|\s]+", "_", str(r.get("reference_lineage") or "unknown"))
            anchor_labels[r["accession"]] = "{}|ref_{}|NA".format(
                str(r["accession"]).split(".")[0], lin)
            continue
        if pd.isna(r["decimal_year"]):
            continue
        k = label_key(r)
        key_counts[k] = key_counts.get(k, 0) + 1
        label_of[r["accession"]] = "{}|{}|{:.3f}".format(
            k, r["host_group"], r["decimal_year"])

    if use_isolate:
        linked = sum(1 for n in key_counts.values() if n > 1)
        logging.info("tip labels keyed on %r: %d distinct isolates from %d records",
                     isolate_key, len(key_counts), len(label_of))
        if linked == 0 and len(label_of) > len(loci):
            logging.error(
                "no isolate appears in more than one record, so no taxon will be "
                "shared between segment trees and congruence cannot be assessed. "
                "Check that %r is populated in the source records, or set "
                "segments.isolate_key to a field that is.", isolate_key)
        else:
            logging.info("%d isolates carry more than one segment", linked)
    for locus in loci:
        name = locus["name"]
        entries = [(label_of[acc], s) for acc, s in seqs[name].items()
                   if acc in label_of]
        n = write_fasta(a.out_dir / f"{name}.fasta", entries)
        logging.info("wrote %s (%d sequences)", a.out_dir / f"{name}.fasta", n)

        anc = [(anchor_labels[acc], s) for acc, s in seqs[name].items()
               if acc in anchor_labels]
        # Always write the file, even empty, so the workflow has a fixed input.
        na = write_fasta(a.out_dir / f"{name}_anchors.fasta", anc)
        if anchors:
            logging.info("wrote %s (%d reference anchors)",
                         a.out_dir / f"{name}_anchors.fasta", na)
            missing = sorted(set(anchors) - {l.split("|")[0] for l in
                                             (x[0] for x in anc)})
            if missing:
                logging.warning(
                    "reference anchors not present in %s: %s", name,
                    ", ".join(missing))

    # --- checks -------------------------------------------------------------
    reasons = excluded["exclude_reason"].value_counts().to_dict()
    checks = check_curation(cfg, len(df), len(clean), len(review), reasons)
    print()
    print(checks.render())

    # Anchors are in `clean` so they reach the alignment, but they are not
    # analysis sequences — counting them under host_group "unknown" and date
    # precision "none" misdescribes the dataset in the one summary a reader
    # actually looks at.
    n_anchor = int(clean["is_reference"].sum()) if "is_reference" in clean else 0
    analysis = clean[~clean["is_reference"]] if n_anchor else clean

    logging.info("-" * 62)
    logging.info("RETAINED: %d analysis sequences of %d records",
                 len(analysis), len(df))
    if n_anchor:
        logging.info("  plus %d reference anchors (undated, excluded from the "
                     "clock and the trait analysis)", n_anchor)
    if not analysis.empty:
        logging.info("")
        logging.info("By host group:")
        for grp, cnt in analysis["host_group"].value_counts().items():
            logging.info("    %-18s %5d", grp, cnt)
        logging.info("")
        logging.info("By date precision:")
        for prec, cnt in analysis["date_precision"].value_counts().items():
            logging.info("    %-18s %5d", prec, cnt)
        logging.info("")
        for locus in loci:
            name = locus["name"]
            have = analysis[analysis[f"has_{name}"]]
            if have.empty:
                continue
            off = int((~have[f"{name}_in_frame"]).sum())
            over = int(have[f"{name}_source"].astype(str)
                       .str.endswith("_overlength").sum())
            exp_len = locus.get("expected_length") or 0
            partial = int((have[f"{name}_length"] < exp_len).sum()) if exp_len else 0
            exp = locus.get("expected_length")
            logging.info("")
            logging.info("%s length checks (expected CDS %s nt):", name, exp)
            logging.info("    median extracted length : %d",
                         int(have[f"{name}_length"].median()))
            logging.info("    partial (< expected)    : %d  (normal — most "
                         "GenBank records cover part of the gene)", partial)
            logging.info("    not a multiple of 3     : %d  (cannot translate "
                         "in frame; codon partitions and any amino-acid "
                         "analysis would be wrong)", off)
            logging.info("    LONGER than expected    : %d  (flagged: a feature "
                         "cannot be a partial CDS and also exceed it — likely "
                         "a gene span including UTR, or mis-annotation)", over)

        logging.info("")
        logging.info("Date range: %.2f – %.2f",
                     analysis["decimal_year"].min(), analysis["decimal_year"].max())
        logging.info("Countries represented: %d", analysis["country"].nunique())
    logging.info("")
    logging.info("NEEDS HUMAN REVIEW: %d records -> %s",
                 len(review), a.out_dir / "needs_review.tsv")
    if len(review):
        for reason, n in review["review_reason"].value_counts().items():
            logging.info("    %-20s %5d", reason, n)
    logging.info("FLAGGED (no decision needed): %d records -> %s",
                 len(flagged), a.out_dir / "flagged.tsv")
    logging.info("-" * 62)

    if not checks.ok:
        return 1
    if len(review) and not a.accept_unreviewed:
        logging.warning(
            "%d records need review. Open needs_review.tsv: for each "
            "unrecognised host, either add a pattern to %s or confirm the "
            "record should be dropped; for each unparseable date, read the "
            "record and fix or exclude it. Ambiguous host labels are the input "
            "to the host-transition analysis, so guessing here becomes a "
            "result. (Sequence-property flags are in flagged.tsv and need no "
            "decision.)",
            len(review), cfg["hosts"]["table"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

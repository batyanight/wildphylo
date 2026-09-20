#!/usr/bin/env python3
"""
Build synthetic GenBank files for testing curation.

Real GenBank downloads are large, gitignored, and change under you. These
fixtures are small, deterministic, and every record exists to exercise one
specific path through 02_curate_metadata.py — so the expected answer is known
independently of what the code does.
"""

import random
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def record(acc, desc, length, seed, *, host="", date="", country="",
           strain="", isolate="", segment="", gene=None, gene_span=None,
           organism="Canine distemper virus"):
    rec = SeqRecord(Seq(dna(length, seed)), id=acc, name=acc.split(".")[0],
                    description=desc)
    rec.annotations["molecule_type"] = "RNA"
    rec.annotations["organism"] = organism
    rec.annotations["taxonomy"] = ["Viruses"]
    q = {"organism": [organism], "mol_type": ["genomic RNA"]}
    if host:    q["host"] = [host]
    if date:    q["collection_date"] = [date]
    if country: q["geo_loc_name"] = [country]
    if strain:  q["strain"] = [strain]
    if isolate: q["isolate"] = [isolate]
    if segment: q["segment"] = [segment]
    rec.features.append(
        SeqFeature(FeatureLocation(0, length), type="source", qualifiers=q))
    if gene:
        s, e = gene_span or (0, length)
        rec.features.append(SeqFeature(
            FeatureLocation(s, e), type="CDS",
            qualifiers={"gene": [gene], "product": [gene]}))
    return rec


def cdv_records():
    """CDV: unsegmented, H gene. Each record exercises one path."""
    r = []
    # 1. clean, annotated H, day-precision date
    r.append(record("AF100001.1", "Canine distemper virus strain A hemagglutinin (H) gene",
                    1824, 1, host="Procyon lotor (raccoon)",
                    date="2015-06-14", country="USA: New York",
                    gene="H", gene_span=(0, 1824)))
    # 2. clean, year-only date
    r.append(record("AF100002.1", "Canine distemper virus H gene", 1824, 2,
                    host="Canis lupus familiaris", date="2013",
                    country="USA", gene="H", gene_span=(0, 1824)))
    # 3. alias spelling 'haemagglutinin'
    r.append(record("AF100003.1", "CDV haemagglutinin", 1824, 3,
                    host="Mustela putorius furo", date="14-Aug-2011",
                    gene="haemagglutinin", gene_span=(0, 1824)))
    # 4. no annotation, description mentions H, full length -> heuristic
    r.append(record("AF100004.1", "Canine distemper virus hemagglutinin gene, complete cds",
                    1820, 4, host="Vulpes vulpes", date="2016-03"))
    # 5. no annotation, partial length -> heuristic_partial -> needs_review
    r.append(record("AF100005.1", "Canine distemper virus hemagglutinin partial cds",
                    600, 5, host="Procyon lotor", date="2018"))
    # 6. VACCINE -> excluded
    r.append(record("AF100006.1", "Canine distemper virus strain Onderstepoort H gene",
                    1824, 6, host="Canis lupus familiaris", date="2001",
                    strain="Onderstepoort", gene="H", gene_span=(0, 1824)))
    # 7. malformed date -> excluded, and NOT a crash
    r.append(record("AF100007.1", "Canine distemper virus H gene", 1824, 7,
                    host="Procyon lotor", date="2013-13-45",
                    gene="H", gene_span=(0, 1824)))
    # 8. implausible year -> excluded
    r.append(record("AF100008.1", "Canine distemper virus H gene", 1824, 8,
                    host="Procyon lotor", date="3013",
                    gene="H", gene_span=(0, 1824)))
    # 9. unrecognised host -> needs_review + excluded
    r.append(record("AF100009.1", "Canine distemper virus H gene", 1824, 9,
                    host="Sciurus carolinensis (grey squirrel)", date="2017",
                    gene="H", gene_span=(0, 1824)))
    # 10. no H at all (N gene only) -> excluded
    r.append(record("AF100010.1", "Canine distemper virus nucleocapsid gene", 1572, 10,
                    host="Procyon lotor", date="2014", gene="N",
                    gene_span=(0, 1572)))
    # 11. no date at all -> excluded, NOT needs_review (nothing to adjudicate)
    r.append(record("AF100011.1", "Canine distemper virus H gene", 1824, 11,
                    host="Procyon lotor", gene="H", gene_span=(0, 1824)))
    # 15. Annotated as H but 1946 nt against an 1824 nt CDS — the real case
    # (KU666057, MH810099). A `gene` feature spanning UTR, or a mis-annotation.
    # Trusted unconditionally before length validation existed.
    r.append(record("AF100015.1", "Canine distemper virus H gene", 1946, 66,
                    host="Procyon lotor", date="2012-11-23",
                    gene="H", gene_span=(0, 1946)))
    # 16. Annotated H, length not a multiple of 3: cannot translate in frame.
    r.append(record("AF100016.1", "Canine distemper virus H gene", 1823, 67,
                    host="Procyon lotor", date="2014-02-10",
                    gene="H", gene_span=(0, 1823)))
    # 13-14. REFERENCE-STYLE records: no /host, no /collection_date. This is
    # what real 1990s reference genomes look like, and why they were all being
    # filtered out before anchors existed.
    r.append(record("AF164967.1", "Canine distemper virus strain A75/17 complete genome",
                    1824, 64, gene="H", gene_span=(0, 1824)))
    r.append(record("Z47762.1", "Canine distemper virus H gene", 1824, 65,
                    gene="H", gene_span=(0, 1824)))
    # 12. the shadowing case: a pinniped
    r.append(record("AF100012.1", "Canine distemper virus H gene", 1824, 12,
                    host="Pusa caspica (Caspian seal)", date="2000-04-02",
                    gene="H", gene_span=(0, 1824)))
    return r


def btv_records():
    """BTV: segmented. Seg-2, Seg-6, Seg-10 across shared isolates."""
    r = []
    specs = [
        ("seg2", "2", "VP2", 2926), ("seg6", "6", "VP5", 1638),
        ("seg10", "10", "NS3", 822),
    ]
    hosts = [("Ovis aries", "2012-05-03"),
             ("Ovis canadensis (bighorn sheep)", "2015"),
             ("Odocoileus virginianus (white-tailed deer)", "2018-09-11"),
             ("Culicoides sonorensis", "2019-07")]
    n = 0
    for hi, (host, date) in enumerate(hosts):
        for name, segnum, prod, length in specs:
            n += 1
            r.append(record(
                f"MN{200000 + n}.1",
                f"Bluetongue virus isolate BTV-{hi+1} segment {segnum} {prod} gene",
                length, 100 + n, host=host, date=date, country="USA",
                isolate=f"BTV-{hi+1}", segment=segnum,
                gene=prod, gene_span=(0, length),
                organism="Bluetongue virus"))
    # one record with no segment qualifier and no annotation -> unusable
    r.append(record("MN299999.1", "Bluetongue virus partial sequence", 400, 999,
                    host="Ovis aries", date="2020", organism="Bluetongue virus"))
    return r


if __name__ == "__main__":
    out = Path(__file__).resolve().parent
    SeqIO.write(cdv_records(), out / "cdv_test.gb", "genbank")
    SeqIO.write(btv_records(), out / "btv_test.gb", "genbank")
    print(f"wrote {out/'cdv_test.gb'} and {out/'btv_test.gb'}")

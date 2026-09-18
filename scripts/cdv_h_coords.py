#!/usr/bin/env python3
"""
cdv_h_coords.py — shared H-gene coordinate handling for scripts 10 and 11.

NOT A PIPELINE STEP. This is a module imported by 10_add_entropy.py and
11_annotate_H.py. It has no main() and does nothing if run directly. It sits
unnumbered among the 00-11 scripts for that reason.

WHY THIS EXISTS
---------------
Scripts 10 and 11 must agree exactly on three things: which sequence defines the
coordinate system, where the reading frame starts, and how alignment columns map
to reference nucleotide positions. When they were computed separately in each
script they could and did disagree, so both now import this module.

THE COORDINATE SYSTEM
---------------------
Auspice datasets are written in *reference nucleotide coordinates*, not alignment
columns. This is not a stylistic choice — it is forced by how Auspice reads a CDS.

An alignment column span is only a valid CDS if the reference has no gaps inside
it. If any other sequence carries an insertion relative to the reference, that
insertion adds columns to the span, the span stops being a multiple of 3, and
Auspice silently drops the feature. Worse, translating the raw column slice
shifts the reading frame at every insertion, so the amino acids downstream are
wrong.

So: `nuc` runs 1..len(ungapped reference), every CDS is expressed in those
coordinates, and codons are assembled from the specific columns holding the
reference's bases (skipping insertion columns). Mutations at insertion columns
have no reference coordinate and are dropped, which is also what augur does.

DOMAIN BOUNDARIES
-----------------
Amino acid positions in the H protein, from Beineke et al. / Sattler et al.
These are literature positions mapped onto whichever reference you choose; say
which reference you used in the methods.
"""

from pathlib import Path

GAPS = set("-.~")
STOPS = {"TAA", "TAG", "TGA"}
BASES = set("ACGT")
AMBIG = set("NnXx?-.~RYKMSWBDHVryknmswbdhv")

# H protein domains: (name, first aa, last aa), 1-based inclusive.
DOMAINS = [
    ("cytoplasmic",   1,   35),
    ("transmembrane", 36,  58),
    ("stalk",         59,  154),
    ("connector",     155, 187),
    ("head",          188, 604),
]

# SLAM (CD150) contact residues. These sit inside the head domain. Modern Auspice
# stacks overlapping CDSs on separate rows, so they no longer need the head to be
# split around them.
SLAM_SITES = [
    ("SLAM_526_529", 526, 529),
    ("SLAM_547_548", 547, 548),
]

COLORS = {
    "H":              "#5097BA",
    "cytoplasmic":    "#4C90C0",
    "transmembrane":  "#8EBC66",
    "stalk":          "#E4B143",
    "connector":      "#DE8A5A",
    "head":           "#DB2823",
    "SLAM_526_529":   "#7F4CA5",
    "SLAM_547_548":   "#B14C8F",
}

DESCRIPTIONS = {
    "H":             "Hemagglutinin (attachment) glycoprotein, full CDS",
    "cytoplasmic":   "H protein cytoplasmic tail",
    "transmembrane": "H protein transmembrane anchor",
    "stalk":         "H protein stalk",
    "connector":     "H protein connecting region",
    "head":          "H protein receptor-binding head (beta-propeller)",
    "SLAM_526_529":  "SLAM/CD150 contact residues 526-529",
    "SLAM_547_548":  "SLAM/CD150 contact residues 547-548",
}


def _codon_table():
    bases = "TCAG"
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    table, k = {}, 0
    for a in bases:
        for b in bases:
            for c in bases:
                table[a + b + c] = aas[k]
                k += 1
    return table


CODON = _codon_table()


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, chunks = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def pick_reference(seqs: dict[str, str], name: str | None, auto: bool):
    """Return (name, sequence). Falls back to the least-gapped sequence."""
    if name and not auto:
        if name not in seqs:
            raise KeyError(
                f"reference {name!r} is not in the alignment. "
                f"Available (first 5): {list(seqs)[:5]}"
            )
        return name, seqs[name]
    chosen = min(seqs, key=lambda k: sum(1 for c in seqs[k] if c in GAPS))
    return chosen, seqs[chosen]


class HFrame:
    """Reference-based coordinate system for one gapped reference sequence.

    Attributes
    ----------
    ref_length   : length of the reference with gaps removed
    n_columns    : number of alignment columns
    orf_start    : 0-based offset of the ORF's first base within the ungapped reference
    n_codons     : codons in the ORF, excluding the terminal stop
    """

    def __init__(self, ref_name: str, ref_seq: str):
        self.ref_name = ref_name
        self.ref = ref_seq.upper()
        self.n_columns = len(self.ref)

        # col_for_refpos[p] is the 1-based alignment column holding reference base p.
        # refpos_for_col[c] is the reverse; insertion columns are simply absent.
        self.col_for_refpos = [None]
        self.refpos_for_col = {}
        n = 0
        for col, ch in enumerate(self.ref, start=1):
            if ch not in GAPS:
                n += 1
                self.col_for_refpos.append(col)
                self.refpos_for_col[col] = n
        self.ref_length = n
        self.ungapped = "".join(c for c in self.ref if c not in GAPS)

        self.orf_start = None
        self.orf_nt = 0
        self.n_codons = 0

    # -- reading frame ------------------------------------------------------

    def find_orf(self):
        """Longest ATG-initiated ORF in the ungapped reference.

        Many GenBank records are mRNAs or genome fragments with untranslated
        flanks, so the CDS rarely starts at base 1 and the frame cannot be
        assumed.
        """
        s = self.ungapped
        best_len, best_start = 0, None
        for f in (0, 1, 2):
            i = f
            while i < len(s) - 2:
                if s[i:i + 3] == "ATG":
                    j = i
                    while j < len(s) - 2:
                        if s[j:j + 3] in STOPS:
                            break
                        j += 3
                    if j - i > best_len:
                        best_len, best_start = j - i, i
                    i = j + 3
                else:
                    i += 3
        self.orf_start = best_start
        self.orf_nt = best_len
        self.n_codons = best_len // 3
        return best_start, best_len

    def internal_stops(self, frame: int) -> int:
        s = self.ungapped[frame:]
        codons = [s[i:i + 3] for i in range(0, len(s) - 2, 3)]
        return sum(1 for c in codons[:-1] if c in STOPS)

    # -- coordinate conversion ---------------------------------------------

    def aa_span_to_ref_nt(self, a1: int, a2: int) -> tuple[int, int]:
        """Amino acid span (1-based, inclusive) -> reference nucleotide span."""
        if self.orf_start is None:
            raise RuntimeError("call find_orf() first")
        return (self.orf_start + (a1 - 1) * 3 + 1, self.orf_start + a2 * 3)

    def codon_columns(self, codon_index: int) -> tuple[int, int, int]:
        """0-based codon index within the ORF -> its three alignment columns.

        Insertion columns are skipped, which is the whole point: this is what
        keeps the reading frame correct across indels.
        """
        p = self.orf_start + 3 * codon_index          # 0-based ungapped offset
        return tuple(self.col_for_refpos[p + k + 1] for k in range(3))

    def feature_codon_columns(self, a1: int, a2: int) -> list[tuple[int, int, int]]:
        return [self.codon_columns(i - 1) for i in range(a1, a2 + 1)]


def resolve_features(frame: HFrame, gene_name: str = "H",
                     domains: bool = False, slam: bool = False,
                     log=print) -> list[tuple[str, int, int]]:
    """Build the ordered feature list as (name, first_aa, last_aa).

    The full-length gene is always emitted first so the entropy panel has a
    feature spanning the whole CDS. Domains and SLAM sites overlap it; Auspice
    2.46+ stacks overlapping CDSs on separate rows rather than discarding them.
    """
    n = frame.n_codons
    feats = [(gene_name, 1, n)]
    if not domains:
        return feats

    for name, a1, a2 in DOMAINS:
        if a1 > n:
            log(f"  note: {name} (aa {a1}-{a2}) lies beyond the {n} aa ORF — skipped")
            continue
        end = a2
        if name == DOMAINS[-1][0]:
            # The head runs to the C-terminus; follow this reference's ORF rather
            # than the literature's exact length.
            if n != a2:
                log(f"  note: head extended to aa {n} to match this reference's ORF "
                    f"(literature boundary is {a2})")
            end = n
        elif a2 > n:
            log(f"  note: {name} truncated at aa {n} (ORF ends there)")
            end = n
        feats.append((name, a1, end))

    if slam:
        for name, a1, a2 in SLAM_SITES:
            if a2 > n:
                log(f"  note: {name} lies beyond the {n} aa ORF — skipped")
                continue
            feats.append((name, a1, a2))
    return feats


def build_annotations(frame: HFrame, feats, gene_name: str = "H") -> dict:
    """Feature list -> the `genome_annotations` object, in reference coordinates."""
    ann = {
        "nuc": {"start": 1, "end": frame.ref_length, "strand": "+", "type": "source"}
    }
    for name, a1, a2 in feats:
        start, end = frame.aa_span_to_ref_nt(a1, a2)
        span = end - start + 1
        # Guaranteed by construction, but this is the exact invariant whose
        # violation used to make features disappear, so assert it loudly.
        assert span % 3 == 0, f"{name}: {span} nt is not a multiple of 3"
        ann[name] = {
            "start": start,
            "end": end,
            "strand": "+",
            "type": "CDS",
            "gene": gene_name,
            "color": COLORS.get(name, "#666666"),
            "display_name": name,
            "description": DESCRIPTIONS.get(name, f"{gene_name} feature {name}"),
        }
    return ann


def backup_if_overwriting(source: Path, output: Path, enabled: bool = True,
                          log=print) -> Path | None:
    """Timestamped copy of `output` before it is overwritten.

    Both scripts edit the dataset in place, so a failure part way through leaves
    a JSON with new annotations and stale mutations — which looks plausible in
    Auspice and is not. Add *.bak.json to .gitignore.
    """
    output = Path(output)
    if not enabled or not output.is_file():
        return None
    from datetime import datetime
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    dest = output.with_suffix(f".{stamp}.bak.json")
    dest.write_bytes(output.read_bytes())
    log(f"backed up existing {output} -> {dest}")
    return dest


def annotation_feature_spans(ann: dict) -> list[tuple[str, int, int]]:
    """CDS features from a genome_annotations object, as (name, start_nt, end_nt)."""
    out = []
    for name, f in ann.items():
        if name == "nuc":
            continue
        if "segments" in f:
            segs = f["segments"]
            out.append((name, int(segs[0]["start"]), int(segs[-1]["end"])))
        else:
            out.append((name, int(f["start"]), int(f["end"])))
    out.sort(key=lambda r: (r[1], r[2]))
    return out


if __name__ == "__main__":
    import sys
    sys.exit("cdv_h_coords.py is a module imported by 10_add_entropy.py and "
             "11_annotate_H.py, not a pipeline step. Run those instead.")

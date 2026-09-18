# Porting status

The scripts carried over from `cdv-phylodynamics` work, but they take the
command-line arguments the old repo used, not the ones the new Snakefile
passes. **The DAG builds cleanly and the rules would then fail at runtime.**
This file records exactly where, so the gap is a task list rather than a
surprise.

Verified by comparing each script's `--help` against what `workflow/Snakefile`
invokes.

### Done

**`01_fetch_sequences.py`** — ported. Entrez query built from `fetch.taxid`,
`fetch.min_length`/`max_length` and `fetch.query_filter`; explicit `--out`/`--acc`;
runs `preflight.check_fetch` against the previous build's accession count before
downloading; a short download is now an error rather than a warning. Email and
API key are read from `NCBI_EMAIL`/`NCBI_API_KEY`, never from the config.

**`02_curate_metadata.py`** — ported, and now multi-locus. Host table, locus
aliases, date plausibility bounds and vaccine patterns all come from config.
Uses `lib.dates` (never raises on malformed input) and `lib.hosts` (word-boundary
matching, with the shadowing audit logged every run). Writes one FASTA per
analysed locus with tip labels already built.

One design change came out of testing it on BTV: **for a segmented pathogen the
tip label keys on the isolate, not the accession.** GenBank assigns each segment
its own accession, so accession-keyed labels leave the shared-taxon set across
segment trees empty — and congruence screening, the entire reason for analysing
segments separately, becomes impossible. `segments.isolate_key` names the field;
curation falls back to `strain`, and logs an error if no isolate appears in more
than one record.

### Still to do

| Script | Snakefile passes | Script accepts | Work |
|---|---|---|---|
| `04_alignment_qc.py` | `--aln --metadata --out-tsv --out-png` | `--aln --metadata --outdir --plot --trim-to` | Explicit output paths so Snakemake can track them. Verified to read `02`'s output correctly: label/metadata cross-check reports 0 disagreements |
| `06_subsample.py` | `--config --out-aln --out-meta` | `--outdir --prefix --clades` | Explicit outputs; read strata from `subsample.strata` rather than the hardcoded clade × host × time; drop the hardcoded `WILD` set in favour of `hosts.wild_groups` |
| `07_make_beast_xml.py` | `--config --metadata --temporal-report --out-xml --out-traits --replicate` | `--out --prefix --rate` | Centre the clock prior on the measured root-to-tip slope from `--temporal-report` when `clock_rate_prior` is null; honour `imprecise_date_policy: interval` by emitting sampled tip dates; read model choices from `beast.*` |
| `09_make_auspice.py` | `--config --alignment --output` | `--title --maintainer --most-recent --trait-key` | Read title/maintainer/colourings from `nextstrain.*`; generate the dataset description from the gate verdicts (see `docs/NEXTSTRAIN.md` §3) |

## Not yet written at all

Referenced by the workflow, no implementation:

- `08b_convergence.py` — ESS per parameter, between-chain overlap. Consumed by
  `lib.preflight.check_convergence`, which is written and tested.
- `08d_drt_summary.py` — date-randomisation summary: the real clock-rate HPD
  must not overlap the randomised replicates'.
- `12_compare_builds.py` — build-to-build comparison gate.
- `13_check_updates.py` — accession-delta check for the scheduled rebuild.

## Not used by the workflow

Carried over because they are useful and were part of the CDV analysis, but not
wired into the Snakefile and CDV-specific in places:

`03_align_and_tree.py`, `05_tree_summary.py`, `08_add_references.py`,
`10_add_entropy.py`, `11_annotate_H.py`, `cdv_h_coords.py`.

`11_annotate_H.py` and `cdv_h_coords.py` are hemagglutinin-specific and will
need generalising against `locus.coordinate_reference` before they apply to
another pathogen.

## Suggested order

1. ~~`01` and `02`~~ — done.
2. `04` and `06` — mechanical argument changes; `04` already consumes `02`'s
   output correctly.
3. `13_check_updates.py` — small, and it unblocks the scheduled rebuild.
4. `08b` — the convergence gate is already specified by its tests.
5. `07` — the largest piece, because of the clock-prior and tip-date work.
6. `09`, `12` — the publishing end.

## Testing without NCBI

`tests/fixtures/make_genbank.py` builds synthetic GenBank files for CDV
(unsegmented) and BTV (segmented). Each record exercises one specific path —
malformed date, implausible year, vaccine strain, unrecognised host, missing
locus, alias spelling, length-heuristic fallback, pinniped — so the expected
answer is known independently of what the code does, and the tests do not
depend on NCBI being reachable.

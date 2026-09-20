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

**`07_make_beast_xml.py`** — ported. The original baked one model into a single
f-string (HKY+G4, relaxed lognormal clock, constant-size coalescent) and
defaulted the clock rate to 7.46e-4 — a CDV H-gene value applied to whatever
pathogen was running.

Now: `beast.tree_prior` selects coalescent_constant, coalescent_skyline or
birth_death; `beast.clock_model` selects relaxed_lognormal or strict;
`imprecise_date_policy: interval` adds a tip-date operator so imprecise dates
are sampled within their window rather than fixed at the midpoint. Model blocks
live in `scripts/lib/beastxml.py`.

The clock prior resolves in order: `--rate`, then the per-locus config value,
then the pathogen config value, then the root-to-tip slope measured from THIS
dataset via `--temporal-report`. A slope of ~1.0 is rejected, because that comes
from a time-scaled tree and is a units check rather than a rate. The XML header
records which source was used.

`preflight.check_clock_prior` runs before writing: a prior disagreeing with the
measured slope by more than 10x is a hard refusal, because that is what a rate
copied from the wrong pathogen looks like and the XML would otherwise consume
days of compute producing an answer determined by the prior.

Why tree priors had to become configurable: a constant-size coalescent assumes a
stable, randomly sampled population. The CDV clade containing the 1994 Serengeti
epidemic has 29 of 82 tips from a single outbreak year. That is not a detail —
it is the wrong model, and it biases both the population size and the dates that
depend on it.

STILL MISSING: the discrete trait (DTA) block. `07` writes a traits file for
BEAUti but does not yet emit the trait partition, transition-rate matrix or
BSSVS in the XML, so the host-transition analysis is still a manual BEAUti step.
That is the next piece.

**`04_alignment_qc.py`** — ported. Takes `--out-tsv`/`--out-png` and creates
missing parent directories; reads QC thresholds from `alignment.*` when not
given on the command line. Verified to read `02`'s output correctly — the
label/metadata cross-check reports 0 disagreements.

**`06_subsample.py`** — ported. Takes `--out-aln`/`--out-meta`, and the
wild-group set now comes from `hosts.wild_groups` instead of a module-level
constant. That constant was CDV's carnivore range; applied to BTV it would have
classified every wild ruminant as not wild and silently broken the
wild/domestic balance in subsampling. CDV resolves to
`procyonid, mustelid, wild_canid, wild_felid, ailurid, ursid, viverrid`;
BTV to `wild_bovid, wild_camelid, wild_cervid, wild_other_ruminant`.

Both scripts previously wrote to a hardcoded `data/processed/`. `06` also
emitted a `*_dates.tsv` there regardless of `--out-meta`, and created the
directory on every run whether or not it was used. A file written outside a
rule's declared outputs cannot be tracked, re-run or cleaned by Snakemake, and
goes stale without anyone noticing.

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
| `09_make_auspice.py` | `--config --alignment --output` | `--title --maintainer --most-recent --trait-key` | Read title/maintainer/colourings from `nextstrain.*`; generate the dataset description from the gate verdicts (see `docs/NEXTSTRAIN.md` §3) |

## Reference anchors

Reference sequences NAME clades; they are not dated tips. All eight usable CDV
nucleotide references lack `/collection_date` and `/host` — five were dropped as
`no_parseable_collection_date`, three as `vaccine_or_vaccine_derived` — so the
ten clades recovered from the first real run could not be named at all.

`references.include_as_anchors: true` admits them as anchors. They are exempt
from every exclusion except failing to carry the locus, and are labelled
`ACC|ref_<lineage>|NA`. The non-numeric date field is load-bearing: every
downstream parser reads it as "no date" and drops the tip from the regression,
the subsample and the trait analysis, while it stays in the alignment and the
tree. `align` concatenates the anchors file, which is always written even when
empty so the rule input is fixed.

Protein accessions are skipped with a warning — 16 of the 24 CDV reference rows
are protein, and nothing resolves them to nucleotide records. References named
in the table but absent from the data are also reported, because a silently
missing reference means a silently unnameable lineage.

## Added, not yet in the workflow

**`05c_cut_clades.py`** — splits a tree into clades and tests each for temporal
signal. Written because the first real CDV run showed that a global tree across
all lineages is not a tip-dating dataset: root-to-tip on 2,008 tips gave
R² = 0.0003 with a negative slope, and after filtering to ≥1,700 nt it reached
only R² = 0.024 with an implied root of 1656. The tree is saturated (16.5
subs/site total length; IQ-TREE warns about long pairwise distances). Clock
signal lives within lineages.

Clades are cut from the tree rather than assigned from a reference panel,
because `config/lineage_references.tsv` has usable nucleotide references for
only three lineages (America-2, America-1/vaccine, Arctic-like) — the other 16
of 24 entries are protein accessions that nothing resolves, and several
lineages present in a global dataset have no reference at all. Name the clades
afterwards by seeing which references fall where.

Not yet wired into the Snakefile. Doing so needs a `clade` wildcard alongside
`locus`, resolved by a checkpoint, so each clade gets its own subsample, BEAST
run and temporal gate.

**`07_make_beast_xml.py`** — the discrete trait is now generated from
`beast.discrete_trait` rather than added by hand in BEAUti (DR-008).
`lib/beastxml.trait_blocks` emits the trait alignment, SVS substitution model,
BSSVS indicators, operators and the ancestral-state tree logger. `--no-trait`
forces a sequence-only XML, which is what the DRT replicates use.

The BEAST_CLASSIC class paths were verified against the installed package —
`examples/testDiscreteSmall.xml` and the source jar — not from documentation.
Both `SVSGeneralSubstitutionModel` and `RobustEigenSystem` live in
`beastclassic`, not `BEAST.base`; every tutorial written before the BEAST 2.7
package rename gives paths that fail to load.

Still unverified: whether the generated XML *runs*. It is well formed, every
idref resolves, and the class paths match the installed package, but no chain
has been started from it. The 1M smoke test in the printed NEXT STEPS is that
check, and it is the last step before the 100M pair.

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
2. ~~`04` and `06`~~ — done.
3. `13_check_updates.py` — small, and it unblocks the scheduled rebuild.
4. `08b` — the convergence gate is already specified by its tests.
5. ~~`07`~~ — done except the discrete-trait block.
6. `07` DTA block — trait partition, rate matrix, BSSVS.
7. `09`, `12` — the publishing end.

## Testing without NCBI

`tests/fixtures/make_genbank.py` builds synthetic GenBank files for CDV
(unsegmented) and BTV (segmented). Each record exercises one specific path —
malformed date, implausible year, vaccine strain, unrecognised host, missing
locus, alias spelling, length-heuristic fallback, pinniped — so the expected
answer is known independently of what the code does, and the tests do not
depend on NCBI being reachable.

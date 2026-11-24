# Build spec — EEEV, Northeast emergence and the enzootic cycle

**Status:** open questions resolved 2026-11-23 (DR-012, DR-013, DR-014). Pathogen
parameters are still unvalidated; the metadata probe results are in issue #2.
**Pipeline:** [wildphylo](https://github.com/batyanight/wildphylo) @ `8197e4b`
**Author:** Batya Nightingale
**Drafted:** 2026-11-23

---

## 1. Why this build

Eastern equine encephalitis virus (EEEV) is a mosquito-borne alphavirus
maintained in a passerine-bird / *Culiseta melanura* cycle in freshwater swamps,
with spillover into horses, humans and farmed ratites. Case fatality is high in
both horses and humans.

The reference phylodynamic analysis (Tan et al. 2018, *J Virol* 92:e00074-18)
sequenced 433 complete genomes spanning 1934–2014 and found a source–sink
structure: Florida enzootic with high diversity and long persistence, New York
and Massachusetts repeatedly re-seeded with low diversity and short persistence.

Two things have happened since that dataset closes:

1. The 2019 and 2024 Northeast seasons, both notable, are not in it.
2. Sequencing has continued, but no maintained public Nextstrain build for EEEV
   exists — unlike flu, RSV or SARS-CoV-2, there is nothing an outbreak
   responder can open and look at.

This build exists to produce that resource and to ask whether the 2018
source–sink conclusion still describes the virus.

## 2. Questions, and what each one costs

Two questions were considered. They are not equally supported by the pipeline as
it stands, and that asymmetry drives the whole plan.

### Q1 — Geography × time (primary scientific question)

*Does the Florida-source / Northeast-sink structure hold when post-2014
sequences are added, or has northern persistence increased?*

**Pipeline support: absent.** Three blockers, all real work, not config:

- `src/wildphylo/config.py:392` warns that `beast.discrete_trait.trait` other
  than `host_group` is "not populated by the curation step". The discrete-trait
  machinery is wired to host end-to-end: `curate` emits `host_group`,
  `traitsignal` prunes on `hosts.dta_states`, `beast xml` builds the trait block
  from it.
- `nextstrain.geo_resolutions` is `[country]` in both shipped configs. EEEV is
  essentially one country, so country resolution carries no information. The
  build needs `division` (US state), which nothing currently parses out of the
  GenBank geo field.
- Subsampling strata (`subsample.strata`) would need `division` as a stratum, or
  the sample is dominated by whichever state sequenced most.

**Estimated work:** a generalisation of the trait path from "host_group" to "a
named discrete trait column", plus a geo parser. This is the kind of change the
repo is structured for — one trait column, one curation step, one validator rule
— but it is a pipeline feature, not a pathogen config.

### Q2 — Host and the enzootic cycle (secondary, and feasible now)

*What do bird-derived sequences show about the enzootic cycle, relative to the
mosquito-pool and mammalian spillover sequences?*

**Pipeline support: native.** This is exactly what wildphylo's discrete-trait
path does, and the trait-signal gate (DR-009) is the right instrument: it asks
whether host associates with the phylogeny at all before any transition table is
believed.

**Feasibility risk: the metadata, not the code.** GenBank `host` is free text and
frequently blank. Most EEEV isolates are mosquito pools. The honest possibility
is that bird-derived sequences are too few to support any host state, in which
case the gate stops the analysis and that refusal is itself the finding.

### Decision

Build Q2 first, because it is runnable against the existing pipeline and its
feasibility resolves in hours rather than months. Treat Q1 as the capstone's
scientific core and the pipeline extension as part of the capstone's
contribution — a defensible framing: the build and the feature are one piece of
work.

**Do not start the Q1 extension before the Q2 metadata probe has run.** If host
coverage is unusable, the whole shape of the project changes.

## 3. Where this lives

DR-011: analysis configs live in the pathogen's own repository, and wildphylo
ships only the template, its smoke config, and configs for pathogens with no
repo yet.

- **Now:** `config/pathogen/eeev.yaml` in wildphylo, alongside `btv.yaml`.
- **On first real build:** move to `eeev-phylodynamics`, pinning a wildphylo
  release, with this spec and the builds beside it.
- `tests/fixtures/configs/eeev.yaml` gets a frozen copy if the test suite needs
  a single-locus vector-borne config to validate against.

## 4. Glossary additions

Terms this build introduces, for `GLOSSARY.md`. Stated here so code, config and
commits use one word each.

| Term | Meaning in this build |
|---|---|
| **enzootic cycle** | the passerine / *Cs. melanura* maintenance cycle; where the lineage persists between outbreak years |
| **bridge vector** | a mosquito feeding on both birds and mammals (*Coquillettidia perturbans*, *Aedes* spp.); the route out of the enzootic cycle |
| **spillover** | transmission out of the enzootic cycle into a dead-end host (equid, human, ratite) |
| **dead-end host** | a host that does not transmit onward; every EEEV mammalian host is one |
| **source–sink** | Florida as persistent source, Northeast as repeatedly re-seeded sink (Tan et al. 2018) |
| **division** | US state, in GenBank geo terms; the resolution Q1 needs |
| **ratite** | emu, ostrich, rhea — farmed birds with EEEV outbreaks; birds, but not enzootic hosts |

Note `dead-end host` against `reservoir_candidates` in the config: EEEV has no
mammalian reservoir, so that list contains only avian groups.

## 5. Pathogen parameters

To go in `config/pathogen/eeev.yaml`. Values marked **VERIFY** are from
literature or memory and must be confirmed against real records before a run.

```yaml
id: eeev
name: Eastern equine encephalitis virus
scientific_name: Alphavirus eastern  # VERIFY: ICTV renamed alphaviruses in 2023

virus:
  genome: ssRNA+
  segmented: false
  recombination_expected: false
  reassortment_expected: false
  expected_rate_range: [1.0e-5, 1.0e-3]  # VERIFY; EEEV is described as slow

fetch:
  taxid: 11021          # VERIFY: must be EEEV, not a parent taxon
  min_length: 300
  max_length: 12500     # genome ~11.7 kb
```

### The Madariaga problem — the main fetch hazard

South American EEEV was reclassified as **Madariaga virus**, a separate species,
on the grounds of distinct geography, diversity and human pathogenicity. Whether
taxid 11021 still returns Madariaga records is the single question that most
affects this dataset: if it does and they are not handled, the tree gains a deep
divergent clade, the root-to-tip slope is computed across two species, and the
clock calibration is wrong in a way that will not look wrong.

Three options:

- **exclude** — cleanest for a North American source–sink question.
- **retain as outgroup anchors** — `references.include_as_anchors: true`, in the
  tree for rooting, out of the clock and the trait analysis.
- **retain as ingroup** — only if the question becomes continental.

~~Default: **retain as anchors**~~ **Resolved: DR-014.** The probe showed taxid
1440170 is a sibling species already excluded by the EEEV query, and that only
one South American record under 11021 is a near-complete genome. Retaining what
arrives would give a single accidental anchor. Anchors are instead a pinned
accession list fetched from 1440170, and South American records under 11021 are
excluded by geo screen.

### Locus

EEEV is unsegmented, but wildphylo's `locus` block is modelled on a *named gene
feature* with GenBank aliases and a coordinate reference. A whole-genome build
does not fit that shape cleanly. **Open design question** — see §9.

**Resolved: DR-013.** The `locus` block gains a `whole_genome: true` form rather
than borrowing the polyprotein aliases, because the two CDSs omit the UTRs and
the junction region and `provenance.json` would then record a polyprotein build
over a genomic alignment. Filter on `nucleotide_completeness` against the
coordinate reference.

### Dates

```yaml
dates:
  min_year: 1930        # Ten Broeck, 1933 Virginia, is the oldest isolate
  max_year: null
  imprecise_date_policy: interval   # VERIFY against DR-001 decision rule
```

The decision rule in DR-001 is per-dataset and depends on substitutions per
genome per year. At a slow alphavirus rate over ~11.7 kb, λ should land well
above the midpoint floor, so `interval` is likely correct — but run the
recommender rather than assuming.

### Clock prior

`clock_rate_prior: null`, per the template and DR-005. The workflow centres the
prior on the measured root-to-tip slope. Do not copy a rate from the alphavirus
literature; the validator treats a >10× disagreement as a hard error precisely
to catch that.

## 6. The vector problem — worse here than in BTV

DR-006 set `annotate` as the default: vector-derived sequences stay in the
phylogeny, contribute no host transitions. For BTV that costs little, because
most BTV sequences come from ruminants.

**For EEEV the majority of sequences are mosquito pools.** Under `annotate`, the
discrete-trait analysis would run on the minority of tips, and the trait-signal
gate would likely fail for lack of states — not because host is unassociated
with the tree, but because the data were set aside.

*Measured 2026-11-23 (issue #2): mosquito-derived is 351 of 577 near-complete
genomes, so the majority claim holds. The predicted gate failure does not. The
remaining tips are equid 97, avian 58 and human 27, which is a workable DTA.
`annotate` would not have broken the analysis; it would have silently narrowed
it to spillover hosts, which is the weaker failure and the harder one to
notice.*

Three readings, none free:

| Policy | Consequence for EEEV |
|---|---|
| `annotate` (DR-006 default) | most tips excluded from the DTA; gate probably fails on state counts; the surviving analysis is about spillover hosts only |
| `state` | mosquito becomes a trait state; every rate into and out of it describes surveillance effort, as DR-006 argues |
| `exclude` | discards the majority of the dataset; not viable |

A fourth option DR-006 does not contemplate, and the one I think is right here:
**split the vector by feeding ecology.** *Cs. melanura* is ornithophilic and
sits inside the enzootic cycle; *Coquillettidia* and *Aedes* are bridge vectors
and sit on the spillover route. These are not one state, and the distinction is
the mechanism the question is about. A sequence from *Cs. melanura* is evidence
about the enzootic cycle in a way a sequence from a bridge vector is not.

This needs a decision record of its own — it amends DR-006 for ornithophilic
vectors rather than overriding it, and the argument has to be made explicitly
because it reintroduces exactly the sampling-effort confound DR-006 rejects.

**Resolved: DR-012.** The split is accepted for EEEV, with DR-006 left standing
as the default for every other vector-borne pathogen, and with the
sampling-effort confound accepted rather than answered. The reporting duties
DR-012 imposes are conditions of the decision, not suggestions.

## 7. Host table — `config/host_groups_eeev.tsv`

Groups: `passerine`, `other_wild_bird`, `farmed_bird`, `equid`, `human`,
`other_mammal`, `vector_ornithophilic`, `vector_bridge`, `vector_other`.

`reservoir_candidates: [passerine]`. `wild_groups` excludes `farmed_bird`.

Ordering hazards to audit (the shadowing audit will catch these; listed so a
reviewer knows they were considered):

- `wood duck`, `mallard` above `duck`
- `house sparrow`, `song sparrow` above `sparrow`
- `emu`, `ostrich`, `rhea` → `farmed_bird`, **not** wild bird. A naive
  wild/domestic split puts them with wild birds and they are the opposite of
  enzootic.
- `pheasant`, `chukar`, `turkey` → `farmed_bird`. Farmed gamebird EEEV outbreaks
  are well documented and would otherwise read as wild-bird signal.
- `culiseta melanura` and `culiseta` above any general mosquito pattern, and
  `coquillettidia` / `aedes` / `culex` mapped to bridge or other — this is the
  §6 decision expressed as table rows, so it cannot be written until §6 is
  settled.
- `horse`, `pony`, `donkey`, `mule`, `zebra` → `equid`
- GenBank host strings for pools are often free text like
  `Culiseta melanura pool` — confirm word-boundary matching handles the
  trailing token.

## 8. Gates — expected behaviour

| Gate | Expectation | If it fails |
|---|---|---|
| **review** | non-trivially dirty; free-text host and mosquito pool strings will fill `needs_review.tsv` | this is the work, not an obstacle; budget real time here |
| **temporal** | should pass comfortably — ~90-year span is unusually good | if it fails, suspect Madariaga contamination of the ingroup first |
| **congruence** | not applicable, unsegmented | — |
| **convergence** | standard 2-chain, ESS ≥ 200 | — |
| **trait signal** | **the real risk.** Depends entirely on §6 and host coverage | a `stop` here is a publishable negative result about surveillance sampling, if written up as one |
| **comparison** | first build has no predecessor | — |

Sampling bias is the honest limitation to state up front rather than discover:
EEEV sequence totals are in the hundreds, not thousands, and sampling is
concentrated where surveillance programmes are funded. Phylogeographic inference
on that is sensitive to sampling design. `temporal_signal.max_concentration`
exists to flag it; report it as a result rather than a caveat.

## 9. Open questions, resolved 2026-11-23

All four were settled after the metadata probe (issue #2). The probe also found
a fifth blocker that §2 did not list.

1. **§6 vector policy** → **DR-012.** Ornithophilic and bridge vectors are
   separate trait states for EEEV, amending DR-006 for this pathogen only. The
   sampling-effort confound DR-006 rejects is accepted here, because 351 of 577
   near-complete genomes are mosquito-derived and `annotate` would leave the
   enzootic cycle resting on 58 avian tips. DR-012 makes the per-state trapping
   counts and the contrast-only reading of the rates mandatory.
2. **Locus model** → **DR-013.** `locus` gains a `whole_genome: true` form. The
   polyprotein aliases would have shipped today and would have recorded a
   polyprotein build over a genomic alignment.
3. **Madariaga** → **DR-014.** Anchors are a pinned accession list from taxid
   1440170, and South American records under 11021 go out by geo screen. The
   premise in §5 was wrong: 1440170 is a sibling species, already excluded, and
   only one South American record under 11021 is a near-complete genome.
4. **Q1 extension scope** → **capstone deliverable.** The build and the
   pipeline extension are one piece of work: Q1 is the scientific core, the
   trait generalisation and geo parsing are the engineering contribution, Q2 is
   the first result. This needs the faculty host conversation and HPC access
   (§11) before the BEAST stage, not before the laptop steps.

### Trait states

`dta_states` is six: `vector_ornithophilic`, `vector_bridge`, `passerine`,
`other_wild_bird`, `farmed_bird`, `mammal_spillover`. Passerines stay separate
from other wild birds because only passerines are enzootic hosts, and lumping
them would let herons and ibises read as enzootic signal. `other_wild_bird` is
thin at 10 near-complete genomes and may not survive the rate-count check.
`equid`, `human` and `other_mammal` collapse into `mammal_spillover`: every one
is a dead-end host, so the distinctions among them carry no transmission
information. `vector_other` and `NO_HOST` stay annotated and out.

Six asymmetric states is 30 rates over roughly 542 host-bearing near-complete
tips. That is the count to put to the rate-count check before the first run.

### The fifth blocker, for Q1

§2 listed three blockers for the geography question. There is a fourth: of 764
US records with a geo string, only 301 name a state. 435 name a sub-state
locality (Oswego, Onondaga, Plymouth and Bristol are counties), and Florida and
New York appear at top level, so the field mixes resolutions within one column.
A `division` parser needs a county-to-state lookup for the majority of US
records, not a split on the colon.

## 10. First steps, in order

1. Metadata probe: pull every EEEV record's accession, date, geo and host.
   Resolves taxid, the Madariaga question and host coverage in one pass.
2. Host coverage report: what fraction of records have a usable host string, and
   what the mosquito/bird/mammal split is. This decides whether Q2 is viable.
3. Hand-check ~10 GenBank records for locus naming, per the new-pathogen
   checklist.
4. Draft `config/pathogen/eeev.yaml` and `config/host_groups_eeev.tsv`.
5. `snakemake --configfile config/pathogen/eeev.yaml -n` and work the warnings.

Steps 1–2 are an evening. Nothing downloads sequence data until step 5.

## 11. Compute

Steps 1–4 are a laptop. Alignment and the ML tree are a laptop. The BEAST stage
is not: `chain_length: 100000000` × 2 chains, plus 20 date-randomisation
replicates at 20M, is on the order of weeks of wall-clock on a MacBook.

On a cluster this is the natural shape of the job, and worth learning on this
build:

- Snakemake has a SLURM executor (`--executor slurm`), so the same Snakefile
  submits each rule as a job rather than running it locally. Per-rule
  `resources:` (threads, memory, walltime) become `sbatch` arguments.
- The 2 chains and 20 DRT replicates are independent — they are the textbook
  case for a **job array**, 22 jobs from one submission, rather than 22
  sequential runs.
- BEAST on a cluster usually comes from a `module load` rather than conda, which
  is the one place the environment stops being reproducible from
  `environment.yml`. Worth pinning the module version in `provenance.json`.
- The gates stay on the login node: they are minutes of Python on files that
  already exist, and gating *before* submitting a week of compute is the whole
  point of the architecture.

NYU has an HPC cluster; access for a capstone would run through the faculty
host, which is one more reason the host conversation comes before the build.

---

## Checklist (`.github/ISSUE_TEMPLATE/new-pathogen.md`)

- [ ] `config/pathogen/eeev.yaml` created from `_template.yaml`
- [ ] `config/host_groups_eeev.tsv` written, shadowing audit clean
- [ ] Locus aliases checked by hand against ~10 real GenBank records
- [ ] Coordinate reference accession chosen and pinned
- [ ] `dates.min_year` defensible (1930, Ten Broeck 1933)
- [ ] `clock_rate_prior` left null
- [ ] `dta_states` declared, rate count checked
- [ ] `snakemake --configfile config/pathogen/eeev.yaml -n` passes
- [ ] All validator warnings read
- [ ] §6 vector policy decided and recorded as a decision record
- [ ] Madariaga handling decided and recorded

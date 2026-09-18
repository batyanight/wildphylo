# wildphylo

**A config-driven phylodynamics pipeline for wildlife disease.**

Point it at a taxon. It curates public sequence data, tests whether that data
can support the analysis you are about to run, and refuses to proceed when it
cannot.

> **Status: early.** The library layer, the gates and the test suite are done
> and CI-verified. The acquisition and Nextstrain ends are being ported from
> the CDV implementation — see [Status](#status).

---

## What this is for

Most phylodynamics failures are not crashes. They are runs that complete and
return a plausible number from a wrong premise:

- a clock prior copied from a different virus
- a time axis silently inverted
- a host label that matched a substring of a longer word
- a trait model with more transition rates than observed transitions
- a TMRCA extrapolated from a decade of sampling and reported to the year

This pipeline is built around gates that catch those before compute is spent,
not after results are published.

## Pipeline

```
fetch → curate → [REVIEW GATE] → align → QC → ML tree → [TEMPORAL GATE]
                                                              ├─ fail → stop
                                                              └─ pass ↓
              [CONGRUENCE GATE, segmented pathogens only] → subsample
                            → BEAST ×N seeds → [CONVERGENCE GATE]
                            → MCC → Auspice → [COMPARISON GATE] → publish
```

**Review gate** — refuses while `needs_review.tsv` is non-empty. Ambiguous host
labels are the *input* to the host-transition analysis, not a detail.

**Temporal gate** — automated root-to-tip regression, date-shuffling
permutation test, and sampling-design diagnostics. On `fail` the Bayesian branch
is not scheduled at all. Replaces the manual TempEst step a user could skip.

**Congruence gate** — for segmented genomes, compares segment topologies and
refuses concatenation when histories are decoupled. Localises incongruence to
individual taxa, so reassortants can be excluded from the trait analysis rather
than generating host transitions that reflect reassortment instead of
transmission.

**Convergence gate** — ESS and between-chain agreement, checked automatically
before any number is reported.

**Comparison gate** — every rebuild is compared against the previous one. A lost
sequence or a TMRCA that moved is a curation regression until proven otherwise.

## Quick start

```bash
git clone https://github.com/batyanight/wildphylo.git
cd wildphylo
conda env create -f environment.yml
conda activate wildphylo
pip install -e ".[dev]"
pytest

# dry run validates the config before anything downloads
snakemake -s workflow/Snakefile --configfile config/pathogen/cdv.yaml --cores 8 -n
```

## Adding a pathogen

```bash
cp config/pathogen/_template.yaml config/pathogen/myvirus.yaml
```

Everything pathogen-specific lives in that one file: taxid, locus aliases, date
plausibility bounds, host table, trait-state map, clock prior, subsampling
strata, Nextstrain metadata. The template fails validation until you replace the
placeholder taxid, so it cannot be run unedited.

Two worked examples ship with the repo:

| | Genome | Exercises |
|---|---|---|
| `cdv.yaml` | ssRNA(−), unsegmented | the baseline: single locus, direct transmission, host-jump inference |
| `btv.yaml` | dsRNA, 10 segments | per-segment clocks and date policies, reassortment screening, vector-borne transmission, serotype structure |

## Some things the validator catches

Run against a real config, `lib/config` and `lib/preflight` report:

```
[trait_states] OK (0E 1W)
  note  5 trait states -> 10 symmetric rates, 162 tips
[clock_prior] OK (0E 0W)
  note  prior 7.460e-04 vs measured root-to-tip slope 8.540e-04 (ratio 0.87)

host table: line 109 'sea lion' -> pinniped is shadowed by line 65 'lion'
  -> wild_felid; move the specific rule above the general one
```

That last one is a real bug found in a real host table: `sea lion` was
unreachable, so a pinniped would have entered the tree as *Panthera leo*.

A clock prior disagreeing with the measured root-to-tip slope by more than 10×
is a hard error, not a warning — that is what a rate copied from the wrong
pathogen looks like.

## Repository layout

```
workflow/Snakefile           the DAG; gates are Snakemake checkpoints
config/pathogen/*.yaml       one file per pathogen
config/host_groups_*.tsv     host normalisation tables
scripts/
  04b_temporal_signal.py     the temporal gate, runnable standalone
  05b_segment_congruence.py  the congruence gate, runnable standalone
  lib/
  config.py                  config loading and validation
  dates.py                   date parsing; never raises on malformed input
  hosts.py                   word-boundary matching + shadowing audit
  temporal.py                root-to-tip, permutation test, verdict
  datepolicy.py              midpoint vs interval, decided per dataset
  segments.py                congruence screening and reassortment localisation
  preflight.py               per-step verification
tests/                       170 tests
docs/DECISIONS.md            why each default is what it is
docs/NEXTSTRAIN.md           organising many builds without chaos
docs/PORTING.md              what is not wired up yet, and exactly why
```

## Documentation

- **[docs/DECISIONS.md](docs/DECISIONS.md)** — decision records. Each one states
  the choice, the justification, what the alternative costs, and the pitfalls of
  what was chosen. Written so it can be cited in a methods section.
- **[docs/NEXTSTRAIN.md](docs/NEXTSTRAIN.md)** — one repo per pathogen, dataset
  naming, staging vs published, generated descriptions, archiving.
- **[docs/PORTING.md](docs/PORTING.md)** — the honest gap list: which scripts
  are not yet config-driven, and the exact argument mismatches.
- **[CONTRIBUTING.md](CONTRIBUTING.md)** — setup, testing conventions.

## Status

Done and CI-verified:

- [x] Config schema and validation, single-locus and segmented
- [x] Hardened date parsing with plausibility bounds
- [x] Host matching with shadowing audit
- [x] Temporal-signal gate with permutation test and inverted-axis detection
- [x] Per-dataset date-policy recommendation
- [x] Segment congruence and reassortment screening
- [x] Per-step preflight verification
- [x] Snakemake DAG with runtime gating, fanning out per locus
      (CDV schedules 9 jobs; BTV 18, across its 3 analysed segments)
- [x] CI: tests on two Python versions, config validation, DAG dry run

**Not yet working end to end.** Acquisition and curation are ported and tested;
four carried-over scripts still take the old repo's command-line arguments, and
four scripts are referenced but unwritten. Every gap is enumerated with its exact flag
mismatch in **[docs/PORTING.md](docs/PORTING.md)**.

- [x] Port `01_fetch` and `02_curate` onto the config and the new libraries
- [ ] Port `04`, `06`, `07`, `09`
- [ ] Write `08b_convergence.py`, `08d_drt_summary.py`, `12_compare_builds.py`,
      `13_check_updates.py`
- [ ] End-to-end run on CDV, then BTV
- [ ] Generated Auspice descriptions carrying gate verdicts

## Acknowledgements

Substantial parts of this repository — the gate architecture, the library
layer, the test suite, the decision records and the BTV config — were drafted
with Claude (Anthropic) in an extended back-and-forth with the author, who
specified the requirements, chose the design direction, and is responsible for
the scientific content.

Recorded here rather than in `CITATION.cff` because authorship carries
accountability for the work, which a model cannot hold; most publishers
(ICMJE, COPE, Nature, Science) say so explicitly. If you reuse this code and
that provenance matters to you, cite the repository and read this section.

## Citation

See [CITATION.cff](CITATION.cff). MIT licensed.

Built from [cdv-phylodynamics](https://github.com/batyanight/cdv-phylodynamics),
which remains the CDV pathogen repository.

## Author

Batya Nightingale — [ORCID 0000-0002-0706-8951](https://orcid.org/0000-0002-0706-8951)

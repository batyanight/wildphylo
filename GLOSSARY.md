# Glossary

One line per term. These are the words used in code, config keys, commit
messages and conversation. If a concept is missing, add it here before naming
it anywhere else.

## Pipeline

- **pathogen config** — the YAML under `config/pathogen/<id>.yaml` that holds everything pathogen-specific; code holds nothing pathogen-specific.
- **build** — one dated run of the workflow for one pathogen config; lives in `builds/<id>/<YYYYMMDD>/`. The newest date is the current build.
- **build date** — the `YYYYMMDD` stamp of a build; pinned with `--config build_date=` to resume or reproduce.
- **stage** — one of the eight top-level steps: fetch, curate, align, tree, gates, beast, auspice, publish. Each is a CLI subcommand group with defined inputs and outputs.
- **step** — a subcommand within a stage (e.g. `gates convergence`).
- **gate** — a step that inspects data and returns a verdict; a failing gate stops the run. Gates: review, temporal, congruence, convergence, trait signal, comparison.
- **verdict** — a gate's result: `pass`, `warn` or `fail`.
- **checks** — the structured list of named checks and their verdicts written by a step to `checks.json`.
- **provenance** — the `provenance.json` record of what a build asked for and got: query, counts, versions, config hash.
- **audit trail** — the subset of build outputs tracked in git: provenance, accession snapshots, QC and temporal reports, clade summaries, build comparison.
- **smoke build** — a tiny config (`cdv-smoke`) used to run the whole workflow quickly as a test.

## Data

- **accession** — a GenBank nucleotide accession; the unit of fetch.
- **accession snapshot** — the `raw/*.acc` list of accessions a build fetched.
- **record** — one GenBank entry with its sequence and metadata.
- **locus** — the gene or region analysed (CDV: H, hemagglutinin). Set in `locus.*`.
- **coordinate reference** — the accession whose CDS defines amino-acid numbering for the locus.
- **metadata** — the per-sequence table: accession, date, host, country, lineage, flags.
- **curation** — turning fetched records into clean metadata and sequences; writes `needs_review.tsv`.
- **needs_review** — rows curation could not resolve (ambiguous host, bad date). The review gate refuses to proceed while non-empty.
- **host** — the organism a sequence was sampled from, as reported.
- **host group** — the normalised host category from `config/host_groups.tsv` (e.g. procyonid, wild_canid).
- **wild group** — a host group counted as wildlife in `hosts.wild_groups`.
- **reservoir candidate** — a host group hypothesised to maintain the virus.
- **vaccine strain** — a sequence matching `config/vaccine_strains.txt`; excluded.
- **lab passaged** — a sequence derived from cell-culture passage; excluded when `exclude.lab_passaged` is true.
- **lineage** — the named CDV genotype (e.g. America-2) assigned by placement against lineage references.
- **lineage reference** — a published accession anchoring a lineage, from `config/lineage_references.tsv`.
- **anchor** — a lineage reference included in the alignment so the tree can be labelled.
- **clade** — a monophyletic group cut from the ML tree for analysis (e.g. clade 3).
- **clade cut** — the step that delineates clades from the ML tree.
- **segment** — one genome segment of a segmented pathogen; each is aligned and treed separately.

## Dates

- **collection date** — the sampling date parsed from the record.
- **imprecise date** — a collection date given to year or month only.
- **date policy** — how imprecise dates are handled: `midpoint`, `interval` or `exclude`. Decided per dataset by `datepolicy`.
- **date interval** — the (earliest, latest) bounds of an imprecise date, sampled in BEAST under the `interval` policy.
- **most recent sample** — the latest collection date; sets the time axis origin.

## Alignment and tree

- **alignment** — the MAFFT output for a locus, after QC.
- **alignment QC** — per-sequence gap, ambiguity and identity-to-consensus checks.
- **ML tree** — the IQ-TREE maximum-likelihood tree.
- **root-to-tip (RTT)** — regression of tip divergence on date; its slope is the measured clock rate.
- **temporal signal** — whether dates explain divergence: RTT R², date-shuffling permutation p, span and concentration.
- **date shuffling** — permuting dates across tips to get a null RTT distribution. Not the DRT.
- **congruence** — agreement between per-segment tree topologies; the congruence gate refuses concatenation when they disagree.
- **reassortant** — a taxon whose segments have different histories; excluded from trait analysis.
- **entropy** — per-site Shannon entropy of the alignment, exported as an Auspice colouring.

## Subsampling

- **subsample** — the stratified subset of sequences that goes to BEAST.
- **stratum** — one cell of the subsampling grid (host group × country × year bin).
- **cap** — the maximum sequences kept per stratum.

## BEAST

- **BEAST XML** — the generated BEAST 2 input; embeds the alignment.
- **clock rate prior** — the prior on substitution rate, from `beast.clock_rate_prior`, sanity-checked against the RTT slope.
- **chain** — one MCMC run; `beast.n_chains` chains with distinct seeds.
- **seed** — the random seed of a chain; part of its output path.
- **chain length** — number of MCMC states.
- **burn-in** — the leading fraction of states discarded.
- **ESS** — effective sample size of a logged parameter.
- **convergence** — all ESS above `min_ess` and between-chain agreement; checked by the convergence gate.
- **MCC tree** — the maximum-clade-credibility summary tree.
- **discrete trait analysis (DTA)** — reconstruction of host-state transitions on the tree.
- **trait** — the discrete character reconstructed (host).
- **trait state (DTA state)** — one of the host categories in `hosts.dta_states`, built from host groups.
- **BSSVS** — Bayesian stochastic search variable selection on transition rates.
- **trait signal** — whether trait is associated with the phylogeny (association index, parsimony score vs tip-shuffled null); checked by the trait signal gate.
- **jump** — one inferred trait transition on a branch.
- **jump history** — the table of jumps across the posterior, summarised per source→sink pair.
- **state through time** — the posterior proportion of lineages in each trait state per year.
- **date randomisation test (DRT)** — rerunning BEAST with dates permuted; the real-data rate must lie outside the randomised rates.

## Nextstrain

- **Auspice dataset** — the JSON under `auspice/` that nextstrain.org renders.
- **dataset name** — the JSON filename; determines the URL. Never contains a date.
- **description** — the Markdown under the tree, generated from gate verdicts and provenance.
- **colouring** — an Auspice colour-by field (host group, lineage, entropy).
- **build comparison** — the diff of the current build against the previous: sequence loss, TMRCA shift, verdict changes.
- **regression** — a build comparison alert; a regression until explained.
- **publish** — the step that runs the comparison gate and copies the Auspice dataset to the pathogen repo.
- **pathogen repo** — the per-pathogen GitHub repo (e.g. cdv-phylodynamics) that holds config, audit trail and the served Auspice dataset.

---
name: New pathogen
about: Checklist for adding a pathogen config
labels: new-pathogen
---

## Pathogen

- Name / scientific name:
- NCBI taxid:
- Genome class (ssRNA+/-, dsRNA, dsDNA, retro):
- Segmented? How many segments?
- Recombination or reassortment expected?

## Before opening a PR

- [ ] `config/pathogen/<id>.yaml` created from `_template.yaml`
- [ ] `config/host_groups_<id>.tsv` written, and the shadowing audit is clean
- [ ] Locus aliases checked by hand against ~10 real GenBank records
- [ ] Coordinate reference accession chosen and pinned
- [ ] `dates.min_year` set to something defensible for this pathogen
- [ ] `clock_rate_prior` left null, OR set with a `clock_rate_source` citation
- [ ] `dta_states` map declared, and the resulting rate count checked
- [ ] `snakemake --configfile config/pathogen/<id>.yaml -n` passes
- [ ] All validator warnings reviewed (not necessarily fixed, but read)

## Host table hazards

List any general/specific pattern pairs needing explicit ordering, e.g.
`bighorn sheep` above `sheep`. The audit catches these, but note them here so a
reviewer knows they were considered.

# CLAUDE.md — wildphylo

Context for Claude Code sessions. Keep this file short; detail lives in the docs
it points to. Personal working preferences live in `~/.claude/CLAUDE.md`, not here.

## What this is
A config-driven phylodynamics pipeline for wildlife disease, generalised from
the earlier `batyanight/cdv-phylodynamics` repo. Goal: works for any input
taxon, with gates that test whether the data can support the priors and
calibrations before compute is spent. Pipeline and gates: `README.md`.

Read before changing anything:
- `docs/PORTING.md` — which scripts from cdv-phylodynamics are ported to the new
  Snakefile CLI and which would still fail at runtime. Update it when you port one.
- `docs/DECISIONS.md` — decision records (DR-nnn). New pipeline-wide choices get
  a DR entry: what, why, the cost of the alternative, how to revisit per dataset.
- `docs/NEXTSTRAIN.md` — this repo is the pipeline; per-pathogen Nextstrain
  repos consume it.

## Layout
- `workflow/Snakefile` — the DAG
- `scripts/NN_*.py` — numbered steps; shared code in `scripts/lib/` (installed as `lib`)
- `config/pathogen/<name>.yaml` — everything pathogen-specific; start from
  `_template.yaml`. Shipped configs: cdv, cdv-1200, cdv-1200-dogwild,
  cdv-clades, cdv-smoke, btv.
- `builds/<config>/<YYYYMMDD>/` — build outputs; the newest date is current.
  Only the audit trail is tracked (provenance, checks, QC, temporal, clade summaries).

## Commands
    conda activate wildphylo
    pip install -e ".[dev]"
    pytest                                 # -m 'not slow' / 'not needs_tools' to skip heavy tests
    ruff check scripts/ tests/             # CI lint job; ruff comes from Homebrew, not the env
    snakemake -s workflow/Snakefile --configfile config/pathogen/cdv-smoke.yaml --cores 8 -n
NCBI CLIs are in the env: `esearch | efetch` (Entrez Direct) and `datasets`
(NCBI Datasets). Prefer them for quick GenBank/taxonomy lookups over writing
one-off Python, and over asking Batya to download files by hand.
BEAST2, Tracer, TempEst, FigTree are installed outside conda (Java).

## Reading this repo cheaply
Search and summarise; don't read large files whole. `rg`, `fd`, `jq`, `mlr` are
on PATH (Homebrew) and `seqkit` is in the env. If one is missing, say so rather
than falling back to `cat`.
- Snakefile (~700 lines): `rg -n '^rule ' workflow/Snakefile`, then read the rule.
- DECISIONS.md (~650 lines): `rg -n '^## DR-' docs/DECISIONS.md`, then read the one DR.
- Build audit files: `jq` on `provenance.json` / `checks.json`; `mlr --itsv` on QC TSVs.
- FASTA and alignments: `seqkit stats`, `seqkit fx2tab -n -l` — never print sequences.
- BEAST XMLs embed the alignment: `rg -n 'chainLength|<prior|<operator|clock' FILE.xml`.
- Posterior `.log` / `.trees`: never read directly. Use the gate scripts
  (`08b_convergence.py` etc.) or `head` for column names.

## Rules
- Nothing pathogen-specific in code: no hardcoded clock rates, taxids, host
  lists or date bounds. They go in the pathogen config.
- NCBI email/key come from `NCBI_EMAIL` / `NCBI_API_KEY`, never the config.
- A gate that fails stops the run. Don't add bypass flags without a DR entry.
- Every rebuild is compared to the previous one (`12_compare_builds.py`); a lost
  sequence or moved TMRCA is a regression until explained.

## Related
Powassan virus / Deer Tick Virus (lineage II) phylodynamics is a candidate
master's thesis that would be a new pathogen config here.

## Cluster
BEAST chains (`BEAST ×N seeds`) are the step that belongs on a cluster. When
discussing that, explain SLURM job arrays, modules and Snakemake cluster
profiles as new material; Batya is learning HPC.

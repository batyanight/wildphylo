# Plan: stage-level CLI and quality pass

Status: **steps 0 to 5 done on branch `stage-cli` (2026-11-10); the Newick
merge (step 4) done 2026-11-11.** Open items are listed at the end. Branch:
`stage-cli`. Decided 2026-11-10 by interview; the answers are recorded here so
the plan does not need re-litigating.

## Decisions

| Question | Decision |
|---|---|
| Scope | Quality pass plus verifying the Auspice end runs from the Snakefile for CDV |
| Structure | Stage-level entry points: one installed `wildphylo` command, stages as subpackages |
| Stages | fetch, curate, align, tree, gates, beast, auspice, publish |
| Lint | ruff `select = ["ALL"]` with a short, justified ignore list, pinned in `pyproject.toml` |
| Types | Hints on every public function; mypy in dev deps and CI |
| Pins | `environment.yml` keeps loose floors; `conda-lock.yml` carries exact versions |
| Glossary | `GLOSSARY.md` drafted by Claude, approved by Batya, then enforced in names |
| Unwired scripts | General ones become subcommands; CDV H-gene scripts leave for the pathogen repo; superseded ones are deleted |
| Phase 2 | cdv-phylodynamics becomes a pathogen repo on top of a tagged wildphylo release (own plan, after this one) |

## Target layout

```
src/wildphylo/
  __init__.py           version
  cli.py                `wildphylo` entry point: argparse with one subparser per stage
  config.py             (from lib/config)        load, validate, resolve_path
  preflight.py          (from lib/preflight)     Checks and check_* functions
  dates.py              (from lib/dates + datepolicy)
  hosts.py              (from lib/hosts)
  temporal.py           (from lib/temporal)
  segments.py           (from lib/segments)
  newick.py             one Newick/annotation parser (merges lib/traittrees + the parser half of lib/traitsignal)
  traitsignal.py        statistics only, on top of newick
  beastlog.py           (from lib/beastlog)
  beastxml.py           (from lib/beastxml)
  stages/
    fetch.py            01
    curate.py           02 (+ review gate)
    align.py            04 (alignment QC)
    tree.py             ml tree, 05c clade cut, 05b congruence helpers
    gates.py            04b temporal, 05b congruence, 08b convergence, 08e trait signal, 08d drt summary, 12 comparison
    beast.py            06 subsample, 07 xml, 08f jump history, 08g state through time
    auspice.py          09 make dataset, 10 entropy
    publish.py          12 compare + copy to pathogen repo, 13 check updates
src/wildphylo/workflow/Snakefile      calls `wildphylo <stage> <step> ...`
tests/                  unchanged files, imports updated; new tests per stage interface
```

`scripts/` is deleted when the last Snakefile rule is switched. `scripts/lib`
becomes `src/wildphylo` with a `src` layout so tests import the installed
package, not the working tree by accident.

Deleted outright (superseded by existing rules): `03_align_and_tree.py`,
`05_tree_summary.py`, `08_add_references.py`. Moved to cdv-phylodynamics in
phase 2: `11_annotate_H.py`, `cdv_h_coords.py`.

## CLI interface (for approval)

Every step takes `--config PATH` and `--build-dir PATH`; the Snakefile passes
both. Steps derive their default inputs and outputs from the build dir, so the
Snakefile only overrides when a rule needs a specific file. Each step writes
`checks.json` into its own subfolder and exits non-zero on a `fail` verdict
unless `--allow-fail` is given (gates only, and only via a DR entry). Usage
errors exit 2 and are logged; gate failures exit 1.

```
wildphylo fetch records        --config --build-dir [--force] [--dry-run]
wildphylo curate metadata      --config --build-dir [--accept-unreviewed]
wildphylo curate review        --config --build-dir            # review gate
wildphylo align qc             --config --build-dir --aln
wildphylo tree clades          --config --build-dir --tree
wildphylo gates temporal       --config --build-dir --tree --metadata
wildphylo gates congruence     --config --build-dir --trees ... --names ...
wildphylo gates convergence    --config --build-dir --logs ...
wildphylo gates trait-signal   --config --build-dir --trees ... --traits
wildphylo gates drt            --config --build-dir --real --reps ...
wildphylo gates compare        --config --build-dir --current --builds-dir
wildphylo beast subsample      --config --build-dir --aln --metadata [--clade]
wildphylo beast xml            --config --build-dir --aln [--randomize-dates --replicate N]
wildphylo beast jumps          --config --build-dir --trees ...
wildphylo beast state-through-time --config --build-dir --trees ...
wildphylo auspice dataset      --config --build-dir --mcc --metadata --gates ...
wildphylo auspice entropy      --config --build-dir --json --alignment
wildphylo publish dataset      --config --build-dir --pathogen-repo PATH
wildphylo publish check-updates --config-dir --builds-dir
wildphylo config validate      --config
```

Argument names stay as they are in the current scripts wherever the script
survives, so the Snakefile diff is mechanical. Script-specific tuning flags
(`--permutations`, `--burnin`, thresholds) stay but default from the config.

## Tooling

- `pyproject.toml`: `[tool.ruff] select = ["ALL"]`, `line-length = 99`,
  `target-version = "py311"`, ignore list with a one-line reason per rule.
  Expect D (docstrings), ANN (annotations) and PLR (complexity) to generate the
  most work; D1xx on tests is ignored.
- `[tool.mypy]` strict on `src/`, `ignore_missing_imports` for Bio and
  pandas until stubs are added. CI gets a `lint` job that runs both.
- `conda-lock -f environment.yml -p osx-arm64 -p linux-64` → `conda-lock.yml`.
  CI installs from the lock; `environment.yml` stays the human spec.
- `[project.scripts] wildphylo = "wildphylo.cli:main"`.

## Order of work

Each step: tests green, ruff clean on touched files, one commit. The Snakefile
and the smoke build must run after every step, so `main` is never broken.

0. Commit the pending `cdv-smoke.yaml` change on its own.
1. Glossary approved. Ruff config pinned; safe autofixes applied; the remaining
   findings fixed by hand in `lib` only. mypy added for `lib`.
2. `src/wildphylo` package with `cli.py` and an empty stage registry; `lib`
   modules move with import rewrites; tests pass against the installed package.
3. One stage at a time, in pipeline order: write the stage module by moving
   the script body behind the subcommand, add a test that invokes the
   subcommand on the fixture data, switch the Snakefile rule, run the smoke
   build, delete the script. fetch → curate → align → tree → gates → beast →
   auspice → publish.
4. Merge the two Newick parsers into `newick.py`; `traitsignal` and the
   jump/state steps consume it.
5. Lock file, CI lint job, docs: README layout section, CLAUDE.md commands,
   PORTING.md retired into this file's completion note.
6. Tag `v0.2.0`. Phase 2 starts from that tag.

## Known entropy this plan removes

- Two independent Newick parsers (`lib/traittrees.py`, `lib/traitsignal.py`).
- 96 root-logger calls across scripts instead of one named logger per module.
- 35 `subprocess.run` calls without `check=`.
- Ten scripts the Snakefile never calls, three of them superseded.
- `docs/PORTING.md` listing `09_make_auspice` as unported after it was ported.

## Phase 2: cdv-phylodynamics (summary; own plan later)

Tag the current code `v1.0.0` for the Zenodo record. Then: scripts removed,
the CDV analysis configs (`cdv`, `cdv-1200`, `cdv-1200-dogwild`, `cdv-clades`) move there from wildphylo (DR-011; frozen copies of two stay as test fixtures in `tests/fixtures/configs/`), environment
pins `wildphylo @ git+https://github.com/batyanight/wildphylo@v0.2.0`, a thin
Snakefile wrapper, tracked files follow the wildphylo audit-trail convention,
`METHODS.md` and `BEAST_MANIFEST.md` stay, notebook and setup notes leave.
Crashed and superseded BEAST runs go to a Zenodo deposit. The Auspice dataset
is regenerated through the pipeline and must pass the comparison gate against
the currently served JSON before it replaces it.

## Completion note (2026-11-10, updated 2026-11-11)

Done: glossary; ruff `ALL` pinned; mypy strict on `src/`; conda-lock;
blocking CI lint job; `src/wildphylo` package with the `wildphylo` command;
every numbered script ported into a stage module and deleted; one Newick
parser (`newick.py`) that `traitsignal` and the jump/state steps consume, and
the only module that knows the NEXUS tree-line grammar (tree counting
included); each step's `--help` shows its run docstring; `align qc` resolves
its thresholds (flag, else config, else default) before computing the flags, so
QC flags for configs that differ from the old CLI defaults differ from earlier
builds; the Snakefile
and both GitHub workflows call only `wildphylo <stage> <step>`; the review
gate is a step instead of inline Snakefile Python; `publish dataset` added;
`docs/PORTING.md` retired. `11_annotate_H.py` and `cdv_h_coords.py` were
removed here; their lint-cleaned versions are in git history at the commit
that deleted them and the originals are in the cdv-phylodynamics repo.

`--build-dir` (2026-11-11): `layout.BuildLayout` is the single source of build
paths; every step takes `--config --build-dir` plus wildcards and derives the
rest through `stepargs.resolve`; the Snakefile rules shrank by 69 lines. A
second offline smoke build (`builds/cdv-smoke/20261113`) through the
simplified rules produced the same QC table, temporal report, subsample and
BEAST XMLs as the explicit-path run, and the date-randomisation XML rule ran
for the first time (its `--randomize-dates` flag had never parsed).
The two builds' ML trees differed because IQ-TREE was unseeded; `tree.seed`
now fixes the search. Verified: two seeded 8-thread runs give the same
topology (140 of 140 splits) and log-likelihood, with branch lengths agreeing
to about 5e-5 (parallel floating-point); one thread is bit-identical.

Behaviour-change review (2026-11-11). The three changes flagged for review
before anything public is rebuilt:

- **QC thresholds now resolved before the flags are computed.** Safe. For CDV
  only `max_gap_fraction` differs from the built-in default (0.60 against
  0.50), so sequences with a gap fraction in that band are no longer flagged.
  Nothing downstream reads those flag columns: the only QC column that changes
  membership is `ungapped_len`, which `tree clades` uses to prune short tips
  and which no threshold touches. The change moves an audit table, not a
  result.
- **Several steps report through logging rather than stdout.** Safe. Every
  wildphylo call in the Snakefile is `2>&1 | tee {log}`, and the only rules
  that capture stdout to a file are the two mafft shell lines, which are tool
  output, not step reports. Log files are unchanged.
- **`beast subsample` refusing an unreadable config.** It was half-refusing.
  An invalid config already exited 2, but YAML that does not parse raises
  ParserError, which a broad `except Exception` swallowed, leaving the CDV
  wild-group set in place. Fixed here, along with the same hole in `align qc`.

Settled (2026-11-11), both found while rewriting the pathogen repo's METHODS:

- **`hosts.dta_states` mapped `domestic_cat`**, a group no host rule produces,
  so the merge never fired. `validate` now warns on such a key. The CDV
  configs keep the five published state names and leave `ailurid` unmapped,
  so the one red panda enters the trait as `?` rather than as a felid.
- **The H CDS bounds** `[21, 1844]` number the CDS from the first base of the
  H gene feature. AF164967.1 has gene H at 7059..9004 and CDS H at 7079..8902;
  the configs now carry the genome-relative pair, which is what the key says
  it holds.

Not done, in rough priority:

- **Smoke build passed against NCBI** (2026-11-11, `builds/cdv-smoke/20261111`,
  untracked, `NCBI_EMAIL` set): fetch capped at 150 of 3232 records, curate →
  review gate → align → qc → ml_tree (`tree.fast`) → temporal gate → subsample
  → BEAST XML ×2 seeds → DRT XML. The gate FAILS on a real capped fetch, and
  correctly: NCBI's first 150 CDV records are a 2011-2026 slice with R² =
  0.0005 and date-shuffling p = 0.39. `cdv-smoke` therefore declares
  `temporal_signal.on_fail: warn`, or nothing past the gate would ever be
  exercised against real data. Run it to the XML targets, not to rule `all`,
  which past the gate pulls the real 100M-state chains.
- **`drt_run` writes BEAST's posterior log into the working directory.** The
  rule's XML names its loggers relatively, so a chain started from the repo
  root drops `H_*_drt<N>.log` beside the Snakefile instead of under the
  replicate's build directory. Exposed when `drt_xml` first ran; the files are
  now gitignored, which is a plaster, not the fix. The rule should run with
  the replicate directory as its working directory, or the writer should
  absolutise the logger paths.
- **A smoke fixture that rule `all` can finish** would need chain lengths
  small enough to run and would then trip the convergence gate on ESS. Worth
  solving; not solved.
- **`publish dataset` is not wired into the Snakefile.** It needs a config key
  for the pathogen repo path, which phase 2 decides.
- **`beast.py`** is ~2000 lines and is the next consolidation target.
- **Report output channels changed** in several steps from print to logging
  (stderr). The Snakefile's `2>&1 | tee` keeps the log files the same.
- `v0.2.0` is tagged locally; push `main` and the tag so phase 2 can pin it.

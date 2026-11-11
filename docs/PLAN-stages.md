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
workflow/Snakefile      calls `wildphylo <stage> <step> ...`
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
unless `--allow-fail` is given (gates only, and only via a DR entry).

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

## Completion note (2026-11-10)

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

Not done, in rough priority:

- **Smoke build passed offline** (2026-11-11, `builds/cdv-smoke/20261112`,
  untracked): curate → review gate → align → qc → ml_tree (`tree.fast`) →
  temporal gate (verdict `weak`, 143 tips) → subsample → BEAST XML ×2 seeds,
  every step through the CLI. The raw input was the first 150 records of the
  2026-10-03 snapshot, standing in for `fetch.max_records: 150`; the fetch
  step itself was exercised only by its tests, because `NCBI_EMAIL` was unset.
  Run a fresh `cdv-smoke` build with the email set once before relying on
  the cap end to end.
- **`--build-dir` defaults:** the layout module is done
  (`src/wildphylo/layout.py`, `BuildLayout`; the Snakefile takes every path
  from it). Steps' `--build-dir` defaults, built on it, are the next task.
- **`publish dataset` is not wired into the Snakefile.** It needs a config key
  for the pathogen repo path, which phase 2 decides.
- **`beast.py`** is ~2000 lines and is the next consolidation target.
- **Report output channels changed** in several steps from print to logging
  (stderr). The Snakefile's `2>&1 | tee` keeps the log files the same.
- `v0.2.0` is tagged locally; push `main` and the tag so phase 2 can pin it.

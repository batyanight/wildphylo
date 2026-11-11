# Contributing

## Setup

```bash
git clone https://github.com/batyanight/wildphylo.git
cd wildphylo
conda env create -f environment.yml      # loose floors; or, for the exact pins:
                                         #   conda-lock install -n wildphylo conda-lock.yml
conda activate wildphylo
pip install -e ".[dev]"                  # installs the `wildphylo` command from src/
wildphylo --help                         # every stage and step
pytest -m 'not slow and not needs_tools' # the fast suite; plain `pytest` runs everything
```

`conda-lock.yml` carries the exact versions; `environment.yml` is the human
spec. After changing `environment.yml`, regenerate the lock and commit both:

```bash
conda-lock -f environment.yml -p osx-arm64 -p linux-64 --lockfile conda-lock.yml
```

## Layout

Code lives in `src/wildphylo/`. `cli.py` is the `wildphylo` command and
`stages/<stage>.py` holds that stage's steps; the other modules are shared
library code (`newick.py`, `temporal.py`, `beastlog.py`, ...). The full map is in
the README's repository layout. Use the terms in `GLOSSARY.md` for names, config
keys and commit messages; if a concept is missing, add it there first.

To add a step, write `run(args) -> int` in the stage module, decorate it with
`@step("<stage>", "<name>", "one-line help")` and give the module an
`add_arguments_<name>` hook. The run function's docstring is the step's
`wildphylo <stage> <name> --help` text, so put the rationale there: why the step
exists, what it checks, how to read its output. The exit code is the verdict.

## Adding a pathogen

```bash
cp config/pathogen/_template.yaml config/pathogen/myvirus.yaml
$EDITOR config/pathogen/myvirus.yaml
wildphylo config validate config/pathogen/myvirus.yaml
snakemake -s src/wildphylo/workflow/Snakefile --configfile config/pathogen/myvirus.yaml -n
```

The dry run validates the config before anything downloads. Work through its
warnings first. The template deliberately fails validation until you replace
the placeholder taxid.

Open a New Pathogen issue and work the checklist.

Once the pathogen has an analysis, its config belongs in that pathogen's own
repository, not here (DR-011). `config/pathogen/` keeps only `_template.yaml`,
`cdv-smoke.yaml` (the pipeline's own smoke fixture) and configs for pathogens
that have no repository yet. Configs the tests need as realistic fixtures go in
`tests/fixtures/configs/`.

## Testing and linting

Every behavioural change needs a test. Steps are tested in-process through
`cli.main([...])` on the small fixtures in `tests/fixtures/`, so the suite
finishes in seconds. Run lint and types after each change, not at the end; CI
runs the same three commands and a failure blocks the merge:

```bash
ruff check src/ tests/
ruff format --check src/ tests/
mypy
```

`ruff` and `mypy` come with the environment and the `dev` extra. The rule set is
pinned in `pyproject.toml` (`select = ["ALL"]` with a justified ignore list),
and mypy runs strict on `src/`.

Two rules specific to this project:

**A check that never fires provides no assurance.** If you add a validation
check, add a test that makes it fire, not only one that passes.

**A failing test on a synthetic fixture sometimes means the fixture is
unrealistic.** The segment congruence tests were first written on six-taxon
trees, which have only three informative splits — so a full clade swap topped
out at RF 0.33, below the 0.40 threshold. The tempting fix was lowering the
threshold, which would have been wrong for real data. Check the fixture before
you change a constant.

## Scientific changes

Anything that changes a number — a threshold, a default policy, a prior — needs
an entry in `docs/DECISIONS.md` recording the decision, the justification, the
alternative's cost, and the pitfalls of the choice made. "It's the standard
approach" is not a justification; say why it is right *here*.

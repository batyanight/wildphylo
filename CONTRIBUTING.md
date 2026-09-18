# Contributing

## Setup

```bash
git clone https://github.com/batyanight/wildphylo.git
cd wildphylo
conda env create -f environment.yml
conda activate wildphylo
pip install -e ".[dev]"
pytest
```

## Adding a pathogen

```bash
cp config/pathogen/_template.yaml config/pathogen/myvirus.yaml
$EDITOR config/pathogen/myvirus.yaml
snakemake -s workflow/Snakefile --configfile config/pathogen/myvirus.yaml -n
```

The dry run validates the config before anything downloads. Work through its
warnings first. The template deliberately fails validation until you replace
the placeholder taxid.

Open a New Pathogen issue and work the checklist.

## Testing

Every behavioural change needs a test. Two rules specific to this project:

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

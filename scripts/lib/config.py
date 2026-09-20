"""
config.py — load and validate a pathogen config.

Validation exists to catch the errors that are expensive later: a clock rate
copied from a different virus, a coordinate reference that is not in the
reference table, a host group named in `wild_groups` that no rule produces.
Each of these produces a plausible-looking result rather than a crash, which
is exactly why they need checking up front.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import yaml

REQUIRED = ["id", "name", "fetch", "dates", "hosts", "tree", "beast"]
# Exactly one of these must be present: `locus` for an unsegmented pathogen,
# `loci` for a segmented one.
REQUIRED_LOCUS = ("locus", "loci")

# Order-of-magnitude sanity bounds on a molecular clock, subs/site/year.
# Below: slower than most DNA viruses and host genomes. Above: faster than the
# fastest RNA viruses. A rate outside this is almost always a unit error.
RATE_FLOOR, RATE_CEILING = 1e-9, 1e-1


class ConfigError(Exception):
    pass


def load(path: str | Path) -> dict:
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"config not found: {path}")
    cfg = yaml.safe_load(path.read_text())
    if not isinstance(cfg, dict):
        raise ConfigError(f"{path}: top level must be a mapping")
    cfg["_path"] = str(path)
    errors, warnings = validate(cfg)
    if errors:
        raise ConfigError(
            f"{path} is not usable:\n" + "\n".join(f"  - {e}" for e in errors))
    cfg["_warnings"] = warnings
    return cfg


def validate(cfg: dict) -> tuple[list[str], list[str]]:
    """Return (errors, warnings). Errors block the run; warnings do not."""
    errors: list[str] = []
    warnings: list[str] = []

    for key in REQUIRED:
        if key not in cfg:
            errors.append(f"missing required section: {key}")
    has = [k for k in REQUIRED_LOCUS if k in cfg]
    if not has:
        errors.append("missing required section: 'locus' (unsegmented) or "
                      "'loci' (segmented)")
    elif len(has) == 2:
        errors.append("config defines both 'locus' and 'loci'; use 'loci' "
                      "alone for a segmented pathogen")
    if errors:
        return errors, warnings

    # --- fetch ---
    taxid = cfg["fetch"].get("taxid")
    if not isinstance(taxid, int) or taxid <= 0:
        errors.append(f"fetch.taxid must be a positive integer, got {taxid!r}")

    # --- dates ---
    d = cfg["dates"]
    min_year = d.get("min_year")
    max_year = d.get("max_year") or dt.date.today().year + 1
    if not isinstance(min_year, int):
        errors.append("dates.min_year must be an integer")
    elif min_year < 1800:
        warnings.append(
            f"dates.min_year = {min_year} is very permissive; it will admit "
            "typo'd years that can distort a calibration")
    if isinstance(min_year, int) and min_year >= max_year:
        errors.append(f"dates.min_year ({min_year}) >= max_year ({max_year})")
    policy = d.get("imprecise_date_policy", "midpoint")
    if policy not in {"midpoint", "interval", "exclude"}:
        errors.append(f"dates.imprecise_date_policy must be midpoint|interval|"
                      f"exclude, got {policy!r}")

    # --- segmented pathogens ---
    segmented = (cfg.get("virus") or {}).get("segmented", False)
    loci = cfg.get("loci")

    if loci is not None:
        if not isinstance(loci, list) or not loci:
            errors.append("loci must be a non-empty list")
            return errors, warnings
        names = [l.get("name") for l in loci]
        if len(set(names)) != len(names):
            errors.append(f"duplicate locus names in loci: {names}")
        analysed = [l for l in loci if l.get("analyse")]
        if not analysed:
            errors.append("no locus has analyse: true; nothing would be built")
        if not segmented:
            warnings.append(
                "config defines multiple loci but virus.segmented is false. If "
                "these really are segments, set it -- the congruence gate and "
                "the reassortment screen both key off that flag")

        seg_cfg = cfg.get("segments") or {}
        if segmented:
            if len(analysed) < 2:
                warnings.append(
                    f"only {len(analysed)} locus marked analyse: true on a "
                    "segmented pathogen. Congruence cannot be assessed from "
                    "one tree, so reassortment would go undetected")
            elif len(analysed) == 2:
                warnings.append(
                    "two analysed segments give a single pairwise RF value. "
                    "Three or more makes it possible to tell which segment is "
                    "the odd one out")
            if seg_cfg.get("concatenate") == "always":
                if (cfg.get("virus") or {}).get("reassortment_expected"):
                    errors.append(
                        "segments.concatenate is 'always' while "
                        "virus.reassortment_expected is true. Concatenating "
                        "reassorting segments averages over conflicting "
                        "histories and yields a tree matching no segment")
                else:
                    warnings.append(
                        "segments.concatenate is 'always'; congruence will not "
                        "be checked before joining segments")
            if not seg_cfg.get("congruence_gate", True):
                warnings.append(
                    "segments.congruence_gate is off; segment topologies will "
                    "not be compared and reassortants will not be flagged")

        # Per-locus clock priors get the same scrutiny as a global one.
        for l in analysed:
            lr = l.get("clock_rate_prior")
            if lr is not None:
                try:
                    lr = float(lr)
                except (TypeError, ValueError):
                    errors.append(f"loci[{l.get('name')}].clock_rate_prior is "
                                  f"not a number: {lr!r}")
                    continue
                if not (RATE_FLOOR <= lr <= RATE_CEILING):
                    errors.append(
                        f"loci[{l.get('name')}].clock_rate_prior {lr:g} is "
                        f"outside [{RATE_FLOOR:g}, {RATE_CEILING:g}]")
                elif not l.get("clock_rate_source"):
                    warnings.append(
                        f"loci[{l.get('name')}] sets a clock prior with no "
                        "clock_rate_source")
            pol = l.get("imprecise_date_policy")
            if pol and pol not in {"midpoint", "interval", "exclude"}:
                errors.append(
                    f"loci[{l.get('name')}].imprecise_date_policy must be "
                    f"midpoint|interval|exclude, got {pol!r}")

        ref = cfg.get("references") or {}
        af = ref.get("assign_from")
        if af and af not in names:
            errors.append(
                f"references.assign_from is {af!r}, which is not a defined "
                f"locus: {names}")

        vp = cfg.get("hosts", {}).get("vector_policy")
        if vp and vp not in {"exclude", "state", "annotate"}:
            errors.append(f"hosts.vector_policy must be exclude|state|annotate, "
                          f"got {vp!r}")
        if vp == "state":
            warnings.append(
                "hosts.vector_policy is 'state': every transition rate to or "
                "from the vector will reflect where sequences were sampled "
                "rather than where the virus is maintained")

        for l in loci:
            if not l.get("aliases"):
                errors.append(f"loci[{l.get('name')}].aliases must list at "
                              "least one name")
        return _finish(cfg, errors, warnings)

    if segmented:
        warnings.append(
            "virus.segmented is true but the config defines a single `locus`. "
            "Only that one segment will be analysed, and reassortment will be "
            "invisible")

    # --- locus (unsegmented) ---
    loc = cfg["locus"]
    if not loc.get("aliases"):
        errors.append("locus.aliases must list at least one name; GenBank "
                      "annotates the same gene inconsistently")
    ref = loc.get("coordinate_reference") or {}
    cds = ref.get("cds")
    if cds and len(cds) == 2:
        try:
            span = int(cds[1]) - int(cds[0]) + 1
            if span % 3:
                errors.append(
                    f"locus.coordinate_reference.cds {cds} spans {span} nt = "
                    f"{span/3:.2f} codons. A CDS must be a multiple of 3; this "
                    "one cannot be translated and would silently break any "
                    "codon-aware annotation downstream")
            exp = loc.get("expected_length")
            if exp and span != exp:
                warnings.append(
                    f"locus.coordinate_reference.cds spans {span} nt but "
                    f"locus.expected_length is {exp}; one of them is wrong")
        except (TypeError, ValueError):
            errors.append(f"locus.coordinate_reference.cds is not a pair of "
                          f"integers: {cds!r}")
    if not ref.get("accession"):
        warnings.append(
            "locus.coordinate_reference.accession is unset. Annotation "
            "positions will fall back to alignment columns, which shift every "
            "time the dataset changes and silently invalidate saved coordinates")

    return _finish(cfg, errors, warnings)


def _finish(cfg: dict, errors: list, warnings: list) -> tuple[list, list]:
    """Checks that apply to segmented and unsegmented configs alike."""
    # --- hosts ---
    h = cfg["hosts"]
    table_path = Path(h.get("table", ""))
    declared_groups: set[str] = set()
    if table_path.is_file():
        try:
            from lib.hosts import load_host_table, audit_host_table
            rules = load_host_table(table_path)
            declared_groups = {r.group for r in rules}
            warnings.extend(f"host table: {w}" for w in audit_host_table(rules))
        except Exception as e:                       # noqa: BLE001
            errors.append(f"hosts.table could not be parsed: {e}")
    else:
        errors.append(f"hosts.table not found: {table_path}")

    if declared_groups:
        for g in h.get("wild_groups", []):
            if g not in declared_groups:
                warnings.append(
                    f"hosts.wild_groups lists {g!r}, which no rule in the host "
                    "table produces; it will never match anything")
        for g in h.get("reservoir_candidates", []):
            if g not in declared_groups:
                warnings.append(
                    f"hosts.reservoir_candidates lists {g!r}, not produced by "
                    "any host rule")

    # --- beast ---
    b = cfg["beast"]
    rate = b.get("clock_rate_prior")
    if rate is not None:
        try:
            rate = float(rate)
        except (TypeError, ValueError):
            errors.append(f"beast.clock_rate_prior is not a number: {rate!r}")
            rate = None
    if rate is not None:
        if not (RATE_FLOOR <= rate <= RATE_CEILING):
            errors.append(
                f"beast.clock_rate_prior = {rate:g} subs/site/yr is outside "
                f"[{RATE_FLOOR:g}, {RATE_CEILING:g}]; check the units")
        if not b.get("clock_rate_source"):
            warnings.append(
                "beast.clock_rate_prior is set but clock_rate_source is empty. "
                "A clock prior with no citation is the single easiest way to "
                "carry an assumption from one pathogen into another")
    else:
        warnings.append(
            "beast.clock_rate_prior is null; the prior will be centred on the "
            "root-to-tip slope from your own data. This is the right default "
            "for a new pathogen but check the resulting prior before the run")

    if b.get("n_chains", 1) < 2:
        warnings.append(
            "beast.n_chains < 2: convergence cannot be assessed from a single "
            "chain, only non-convergence")
    seeds = b.get("seeds") or []
    if len(seeds) < b.get("n_chains", 1):
        errors.append(
            f"beast.seeds has {len(seeds)} entries but n_chains is "
            f"{b.get('n_chains')}; runs would not be reproducible")
    if len(set(seeds)) != len(seeds):
        errors.append("beast.seeds contains duplicates; chains would be identical")

    dt_cfg = b.get("discrete_trait") or {}
    if dt_cfg.get("enabled") and declared_groups:
        trait = dt_cfg.get("trait", "host_group")
        if trait != "host_group":
            warnings.append(
                f"beast.discrete_trait.trait is {trait!r}; only 'host_group' is "
                "populated by the curation step")
        # Count states AFTER the declared collapse, not raw host-table groups.
        # Counting before the map reports a rate count the analysis will never
        # use, and contradicts preflight.check_trait_states.
        mapping = h.get("dta_states") or {}
        if mapping:
            states = {mapping[g] for g in declared_groups if g in mapping}
            other = h.get("dta_other")
            if other:
                states.add(other)
            unmapped = sorted(g for g in declared_groups if g not in mapping)
            if unmapped and not other:
                warnings.append(
                    "host groups absent from hosts.dta_states and excluded "
                    f"from the trait analysis: {', '.join(unmapped)}")
            n_states = len(states)
        else:
            n_states = len(declared_groups)
            warnings.append(
                "hosts.dta_states is unset, so every host group becomes a "
                "trait state and the state count depends on which sequences "
                "happen to be subsampled")
        if dt_cfg.get("symmetric", True):
            n_rates = n_states * (n_states - 1) // 2
        else:
            n_rates = n_states * (n_states - 1)
        warnings.append(
            f"discrete trait has {n_states} states -> {n_rates} transition "
            "rates to estimate. Rates far outnumbering informative transitions "
            "is the usual cause of an unresolved root")

    # --- temporal gate ---
    ts = cfg.get("temporal_signal", {})
    if ts.get("on_fail") not in {"stop", "warn", None}:
        errors.append("temporal_signal.on_fail must be stop|warn")
    if ts.get("on_fail") == "warn":
        warnings.append(
            "temporal_signal.on_fail = warn: a dataset with no clock signal "
            "will still be handed to BEAST, which will return a confident "
            "posterior for a rate that is not identifiable")

    return errors, warnings


def main(argv: list[str] | None = None) -> int:
    """
    Validate configs from the command line.

        python -m lib.config config/pathogen/*.yaml

    Exists so CI does not have to embed Python inside a YAML block scalar --
    which is how ci.yml was silently made invalid, meaning no CI ran at all.
    Emits GitHub Actions annotations when running under Actions.
    """
    import argparse
    import os

    ap = argparse.ArgumentParser(description="Validate pathogen configs.")
    ap.add_argument("configs", nargs="+")
    ap.add_argument("--strict", action="store_true",
                    help="Treat warnings as failures.")
    a = ap.parse_args(argv)

    gha = bool(os.environ.get("GITHUB_ACTIONS"))
    failed = 0
    for path in a.configs:
        name = Path(path).name
        if name.startswith("_"):
            print(f"skip     {path} (template)")
            continue
        try:
            cfg = load(path)
        except ConfigError as e:
            failed += 1
            msg = str(e).replace("\n", " ")
            print(f"::error file={path}::{msg}" if gha else f"INVALID  {path}\n{e}")
            continue
        print(f"valid    {path}  ({cfg['name']})")
        for w in cfg.get("_warnings", []):
            print(f"::warning file={path}::{w}" if gha else f"  warn   {w}")
            if a.strict:
                failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main())

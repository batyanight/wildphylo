#!/usr/bin/env python3
"""
07_make_beast_xml.py — generate a BEAST2 XML from a labelled alignment.

Produces a standard tip-dated analysis:
    HKY+G4 substitution model
    uncorrelated relaxed lognormal clock, rate prior centred on --rate
    constant-population coalescent
    tip dates parsed from the labels

Also writes a traits file (taxon -> host group) for adding discrete trait analysis
in BEAUti, and can generate date-randomised replicates for the temporal signal test.

Usage
-----
    # main analysis
    python scripts/07_make_beast_xml.py \\
        --aln data/processed/H_clade_3_clean.fasta \\
        --out beast/clade3 --rate 7.46e-4 --chain 100000000

    # 20 date-randomised replicates for the DRT
    python scripts/07_make_beast_xml.py \\
        --aln data/processed/H_clade_3_clean.fasta \\
        --out beast/drt --rate 7.46e-4 --chain 20000000 --randomize-dates 20

Tip labels must be  ACCESSION|host_group|decimal_year

IMPORTANT
---------
Open the generated XML in BEAUti before running it (File > Load). BEAUti will
report any incompatibility with your BEAST version, and lets you inspect every
prior. Do not run a chain for days without that check. This script targets
BEAST 2.7.x namespaces.
"""

import argparse
import json
import random
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

# Namespaces BEAST searches when resolving an unqualified spec= name.
# beast.base.inference.operator is required for ScaleOperator / DeltaExchangeOperator;
# beast.base.evolution.operator for the tree operators. Omitting either produces
# "Class could not be found" errors at load time.
NS = ("beast.base.evolution.alignment:beast.base.evolution.tree.coalescent:"
      "beast.base.util:beast.base.evolution.nuc:"
      "beast.base.evolution.operator:beast.base.evolution.operator.kernel:"
      "beast.base.inference.operator:beast.base.inference.operator.kernel:"
      "beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:"
      "beast.base.evolution.likelihood:beast.base.inference:"
      "beast.base.inference.parameter:beast.base.inference.distribution:"
      "beast.base.evolution.branchratemodel:beast.base.evolution.speciation")


sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config import load as load_config, ConfigError  # noqa: E402
from lib.beastxml import (  # noqa: E402
    resolve_clock_rate, tree_prior_blocks, clock_blocks, tip_date_blocks,
    trait_blocks, n_trait_rates)
from lib.preflight import check_clock_prior  # noqa: E402


def measured_slope_from(report: Path | None) -> float | None:
    """
    Root-to-tip slope from the temporal gate's JSON report.

    On a tree whose branch lengths are substitutions per site this IS a
    clock-rate estimate, and it is measured from the dataset actually being
    analysed rather than inherited from a paper about another virus.

    Returns None for a time-scaled tree (slope ~1.0), because there the slope
    is a units check, not a rate — using it as a prior would centre the clock
    on 1.0 subs/site/year.
    """
    if not report or not Path(report).is_file():
        return None
    try:
        data = json.loads(Path(report).read_text())
        regs = data.get("regressions") or []
        if not regs:
            return None
        slope = float(regs[0].get("slope", 0.0))
    except Exception:                                    # noqa: BLE001
        return None
    if slope <= 0:
        return None
    if abs(slope - 1.0) < 0.05:
        print("NOTE: the temporal report's slope is ~1.0, so it came from a "
              "time-scaled tree and is not a substitution rate. Ignoring it; "
              "supply a rate in the config or run the gate on the ML tree.")
        return None
    return slope


def uncertainty_map(metadata: Path | None, labels) -> dict[str, float]:
    """Per-tip date uncertainty, for tips whose date is imprecise."""
    if not metadata or not Path(metadata).is_file():
        return {}
    try:
        import pandas as pd
    except ImportError:
        return {}
    df = pd.read_csv(metadata, sep="\t")
    if "date_uncertainty_years" not in df.columns:
        return {}
    by_acc = {}
    for _, r in df.iterrows():
        acc = str(r.get("accession", "")).split(".")[0]
        try:
            u = float(r.get("date_uncertainty_years") or 0.0)
        except (TypeError, ValueError):
            u = 0.0
        by_acc[acc] = u
    out = {}
    for lbl in labels:
        u = by_acc.get(lbl.split("|")[0], 0.0)
        if u > 0.01:                      # day-precision tips need no sampling
            out[lbl] = u
    return out


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, chunks = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def parse_labels(seqs):
    out = {}
    bad = []
    for label in seqs:
        parts = label.split("|")
        if len(parts) < 3:
            bad.append(label); continue
        try:
            year = float(parts[-1])
        except ValueError:
            bad.append(label); continue
        out[label] = (parts[0], parts[1], year)
    return out, bad


def build_xml(seqs, info, rate, chain, log_every, prefix, seed_note="",
              *, tree_prior="coalescent_constant", clock_model="relaxed_lognormal",
              date_policy="midpoint", uncertainties=None, rate_note="",
              trait=None):
    """
    Assemble the XML as text.

    The model is no longer fixed. tree_prior, clock_model and the tip-date
    policy come from the pathogen config, because a constant-size coalescent
    and a relaxed lognormal clock are choices, not constants — and applying
    them to, say, a clade containing a mass-mortality outbreak is simply the
    wrong model.
    """
    seq_blocks = "\n".join(
        f'        <sequence id="seq_{lbl}" spec="Sequence" taxon="{lbl}" '
        f'totalcount="4" value="{s}"/>'
        for lbl, s in seqs.items())

    dates = ",".join(f"{lbl}={info[lbl][2]:.4f}" for lbl in seqs)
    traits = ",".join(f"{lbl}={info[lbl][1]}" for lbl in seqs)
    n = len(seqs)

    tp = tree_prior_blocks(tree_prior, prefix, n)
    ck = clock_blocks(clock_model, prefix, rate, n)
    td = tip_date_blocks(date_policy, prefix, uncertainties or {},
                         {k: info[k][2] for k in seqs})

    tp_state, tp_prior, tp_operators, tp_log = (
        tp["state"], tp["prior"], tp["operators"], tp["log"])
    ck_state, ck_prior, ck_branchrate, ck_operators, ck_log, ck_ref = (
        ck["state"], ck["prior"], ck["branchrate"], ck["operators"],
        ck["log"], ck["ref"])
    td_operators = td.get("operators", "")
    td_priors = td.get("priors", "")
    td_log = td.get("log", "")
    sampled_note = (f" ({td['n_sampled']} tips sampled within their window)"
                    if td.get("n_sampled") else "")

    # Discrete trait. Absent -> every fragment is empty and the XML is exactly
    # what it was before, so a sequence-only analysis is unaffected.
    tr = trait or {}
    tr_data = tr.get("data", "")
    tr_state = tr.get("state", "")
    tr_prior = tr.get("prior", "")
    tr_likelihood = tr.get("likelihood", "")
    tr_operators = tr.get("operators", "")
    tr_log = tr.get("log", "")
    tr_treelog = tr.get("treelog", "").replace("__LOGEVERY__", str(log_every))
    trait_note = (f"discrete trait  : {tr['n_states']} states, {tr['n_rates']} rates"
                  if tr else "discrete trait  : none (sequence only)")

    return f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Generated by 07_make_beast_xml.py
     {n} sequences
     clock rate prior : {rate:g} subs/site/year
     rate provenance  : {rate_note}
     clock model      : {clock_model}
     tree prior       : {tree_prior}
     tip dates        : {date_policy}{sampled_note}
     {trait_note}
     {seed_note}
     VERIFY IN BEAUti BEFORE RUNNING. -->
<beast beautitemplate="Standard" beautistatus="" namespace="{NS}" required="" version="2.7">

    <data id="{prefix}" spec="Alignment" name="alignment">
{seq_blocks}
    </data>

{tr_data}

    <map name="Uniform">beast.base.inference.distribution.Uniform</map>
    <map name="Exponential">beast.base.inference.distribution.Exponential</map>
    <map name="LogNormal">beast.base.inference.distribution.LogNormalDistributionModel</map>
    <map name="Normal">beast.base.inference.distribution.Normal</map>
    <map name="Gamma">beast.base.inference.distribution.Gamma</map>
    <map name="OneOnX">beast.base.inference.distribution.OneOnX</map>
    <map name="prior">beast.base.inference.distribution.Prior</map>
    <map name="Beta">beast.base.inference.distribution.Beta</map>

    <run id="mcmc" spec="beast.base.inference.MCMC" chainLength="{chain}">

        <state id="state" spec="State" storeEvery="5000">
            <tree id="Tree.t:{prefix}" spec="beast.base.evolution.tree.Tree" name="stateNode">
                <trait id="dateTrait.t:{prefix}" spec="beast.base.evolution.tree.TraitSet"
                       traitname="date-forward" value="{dates}">
                    <taxa id="TaxonSet.{prefix}" spec="TaxonSet">
                        <alignment idref="{prefix}"/>
                    </taxa>
                </trait>
                <taxonset idref="TaxonSet.{prefix}"/>
            </tree>
            <parameter id="kappa.s:{prefix}" spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>
            <parameter id="gammaShape.s:{prefix}" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            <parameter id="freqParameter.s:{prefix}" spec="parameter.RealParameter" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>
{tp_state}
{ck_state}
{tr_state}
        </state>

        <init id="RandomTree.t:{prefix}" spec="RandomTree" estimate="false" initial="@Tree.t:{prefix}" taxa="@{prefix}">
            <populationModel id="ConstantPopulation0.t:{prefix}" spec="ConstantPopulation">
                <parameter id="randomPopSize.t:{prefix}" spec="parameter.RealParameter" name="popSize">1.0</parameter>
            </populationModel>
        </init>

        <distribution id="posterior" spec="CompoundDistribution">

            <distribution id="prior" spec="CompoundDistribution">
{tp_prior}
                <prior id="KappaPrior.s:{prefix}" name="distribution" x="@kappa.s:{prefix}">
                    <LogNormal id="LogNormalDistributionModel.0" name="distr" M="1.0" S="1.25"/>
                </prior>
                <prior id="GammaShapePrior.s:{prefix}" name="distribution" x="@gammaShape.s:{prefix}">
                    <Exponential id="Exponential.0" name="distr" mean="1.0"/>
                </prior>
                <!-- Rate prior. S=1.0 is wide enough that the data can move it.
                     Provenance: {rate_note} -->
{ck_prior}
{td_priors}
{tr_prior}
            </distribution>

            <distribution id="likelihood" spec="CompoundDistribution" useThreads="true">
                <distribution id="treeLikelihood.{prefix}" spec="TreeLikelihood" data="@{prefix}" tree="@Tree.t:{prefix}">
                    <siteModel id="SiteModel.s:{prefix}" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape.s:{prefix}">
                        <parameter id="mutationRate.s:{prefix}" spec="parameter.RealParameter" estimate="false" lower="0.0" name="mutationRate">1.0</parameter>
                        <parameter id="proportionInvariant.s:{prefix}" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>
                        <substModel id="hky.s:{prefix}" spec="HKY" kappa="@kappa.s:{prefix}">
                            <frequencies id="estimatedFreqs.s:{prefix}" spec="Frequencies" frequencies="@freqParameter.s:{prefix}"/>
                        </substModel>
                    </siteModel>
{ck_branchrate}
                </distribution>
{tr_likelihood}
            </distribution>
        </distribution>

        <!-- Classic operators: stable across BEAST 2.6/2.7. The Bactrian kernel
             variants mix slightly better but live in packages whose paths have
             moved between versions. -->
        <operator id="KappaScaler.s:{prefix}" spec="ScaleOperator" parameter="@kappa.s:{prefix}" scaleFactor="0.5" weight="0.1"/>
        <operator id="gammaShapeScaler.s:{prefix}" spec="ScaleOperator" parameter="@gammaShape.s:{prefix}" scaleFactor="0.5" weight="0.1"/>
        <operator id="FrequenciesExchanger.s:{prefix}" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:{prefix}"/>
{tp_operators}
{ck_operators}
{td_operators}
{tr_operators}

        <operator id="CoalescentConstantTreeScaler.t:{prefix}" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:{prefix}" weight="3.0"/>
        <operator id="CoalescentConstantTreeRootScaler.t:{prefix}" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:{prefix}" weight="3.0"/>
        <operator id="CoalescentConstantUniformOperator.t:{prefix}" spec="Uniform" tree="@Tree.t:{prefix}" weight="30.0"/>
        <operator id="CoalescentConstantSubtreeSlide.t:{prefix}" spec="SubtreeSlide" tree="@Tree.t:{prefix}" weight="15.0"/>
        <operator id="CoalescentConstantNarrow.t:{prefix}" spec="Exchange" tree="@Tree.t:{prefix}" weight="15.0"/>
        <operator id="CoalescentConstantWide.t:{prefix}" spec="Exchange" isNarrow="false" tree="@Tree.t:{prefix}" weight="3.0"/>
        <operator id="CoalescentConstantWilsonBalding.t:{prefix}" spec="WilsonBalding" tree="@Tree.t:{prefix}" weight="3.0"/>

        <logger id="tracelog" spec="Logger" fileName="{prefix}.log" logEvery="{log_every}" model="@posterior" sanitiseHeaders="true" sort="smart">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood.{prefix}"/>
            <log id="TreeHeight.t:{prefix}" spec="beast.base.evolution.tree.TreeStatLogger" tree="@Tree.t:{prefix}"/>
            <log idref="kappa.s:{prefix}"/>
            <log idref="gammaShape.s:{prefix}"/>
            <log idref="freqParameter.s:{prefix}"/>
{tp_log}
{ck_log}
{td_log}
{tr_log}
            <log id="rate.c:{prefix}" spec="beast.base.evolution.RateStatistic" branchratemodel="@{ck_ref}" tree="@Tree.t:{prefix}"/>
        </logger>

        <logger id="screenlog" spec="Logger" logEvery="{log_every}">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
        </logger>

        <logger id="treelog.t:{prefix}" spec="Logger" fileName="{prefix}.trees" logEvery="{log_every}" mode="tree">
            <log id="TreeWithMetaDataLogger.t:{prefix}" spec="beast.base.evolution.TreeWithMetaDataLogger" branchratemodel="@{ck_ref}" tree="@Tree.t:{prefix}"/>
        </logger>

{tr_treelog}

        <operatorschedule id="OperatorSchedule" spec="OperatorSchedule"/>
    </run>
</beast>
'''


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True, type=Path)
    ap.add_argument("--out", type=Path, help="Output directory (legacy)")
    ap.add_argument("--out-xml", type=Path, help="Explicit XML path.")
    ap.add_argument("--out-traits", type=Path, help="Explicit traits path.")
    ap.add_argument("--config", type=Path,
                    help="Pathogen config: clock rate, clock model, tree prior, "
                         "chain length and date policy all come from here.")
    ap.add_argument("--metadata", type=Path,
                    help="Curated metadata, for per-tip date uncertainty.")
    ap.add_argument("--temporal-report", type=Path,
                    help="Temporal gate JSON. When the config has no "
                         "clock_rate_prior, the prior is centred on the "
                         "root-to-tip slope measured from THIS dataset.")
    ap.add_argument("--locus", help="Which locus in a multi-locus config.")
    ap.add_argument("--replicate", type=int, default=None,
                    help="Date-randomisation replicate number.")
    ap.add_argument("--prefix", default=None)
    ap.add_argument("--rate", type=float, default=None,
                    help="Override the clock rate prior, subs/site/year.")
    ap.add_argument("--chain", type=int, default=None)
    ap.add_argument("--log-every", type=int, default=None,
                    help="Default: chain/10000, giving ~10k samples")
    ap.add_argument("--randomize-dates", type=int, default=0,
                    help="Generate N date-randomised replicates for the DRT")
    ap.add_argument("--no-trait", action="store_true",
                    help="Write a sequence-only XML even when the config "
                         "enables the discrete trait. Used for the DRT, where "
                         "only the clock rate matters.")
    ap.add_argument("--seed", type=int, default=20260819)
    args = ap.parse_args()

    # --- config -------------------------------------------------------------
    cfg, locus_cfg = {}, {}
    if args.config:
        try:
            cfg = load_config(args.config)
        except ConfigError as e:
            print(e, file=sys.stderr)
            return 2
        if "loci" in cfg:
            wanted = args.locus
            for l in cfg["loci"]:
                if not l.get("analyse"):
                    continue
                if wanted is None or l["name"] == wanted:
                    locus_cfg = l
                    break
        else:
            locus_cfg = cfg.get("locus", {})

    beast_cfg = cfg.get("beast", {})
    tree_prior = beast_cfg.get("tree_prior", "coalescent_constant")
    clock_model = beast_cfg.get("clock_model", "relaxed_lognormal")
    date_policy = (locus_cfg.get("imprecise_date_policy")
                   or cfg.get("dates", {}).get("imprecise_date_policy")
                   or "midpoint")
    if args.chain is None:
        args.chain = int(beast_cfg.get("chain_length", 100_000_000))

    if not args.aln.is_file():
        print(f"Alignment not found: {args.aln}", file=sys.stderr)
        return 1

    seqs = read_fasta(args.aln)
    info, bad = parse_labels(seqs)
    if bad:
        print(f"ERROR: {len(bad)} labels not in ACCESSION|host|year form, e.g. {bad[:3]}",
              file=sys.stderr)
        return 1
    seqs = {k: v for k, v in seqs.items() if k in info}

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        print(f"ERROR: sequences differ in length {sorted(lengths)[:4]} — not aligned.",
              file=sys.stderr)
        return 1

    prefix = args.prefix or args.aln.stem.replace("-", "_")
    log_every = args.log_every or max(args.chain // 10000, 1000)

    # --- clock rate ---------------------------------------------------------
    # The old default was 7.46e-4 — a CDV H-gene rate applied to whatever
    # pathogen happened to be running. Now: an explicit override, else the
    # config, else the slope measured from this dataset.
    slope = measured_slope_from(args.temporal_report)
    if args.rate is not None:
        rate, rate_note = args.rate, "command line --rate"
    else:
        try:
            rate, rate_note = resolve_clock_rate(cfg, slope, locus_cfg)
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1

    if cfg:
        checks = check_clock_prior(cfg, slope)
        print(checks.render())
        if not checks.ok:
            print("\nRefusing to write an XML whose clock prior disagrees with "
                  "the data. Resolve this before committing days of compute.",
                  file=sys.stderr)
            return 1

    xml_path = args.out_xml or (args.out / f"{prefix}.xml" if args.out
                                else Path(f"{prefix}.xml"))
    traits_path = args.out_traits or (args.out / f"{prefix}_traits.txt"
                                      if args.out else Path(f"{prefix}_traits.txt"))
    xml_path.parent.mkdir(parents=True, exist_ok=True)
    traits_path.parent.mkdir(parents=True, exist_ok=True)
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)

    years = [info[k][2] for k in seqs]
    hosts = {}
    for k in seqs:
        hosts[info[k][1]] = hosts.get(info[k][1], 0) + 1

    print(f"{len(seqs)} sequences, {lengths.pop()} sites")
    print(f"date range: {min(years):.2f} - {max(years):.2f} ({max(years)-min(years):.1f} yrs)")
    print(f"host groups: {dict(sorted(hosts.items(), key=lambda kv: -kv[1]))}")
    print(f"rate prior mean: {rate:g}   ({rate_note})")
    print(f"clock model: {clock_model}   tree prior: {tree_prior}   "
          f"tip dates: {date_policy}")
    print(f"chain: {args.chain:,}   log every: {log_every:,}\n")

    thin = {h: n for h, n in hosts.items() if n < 5}
    if thin:
        print(f"WARNING: host groups with <5 sequences: {thin}")
        print("  Discrete trait analysis will be unreliable for these — consider merging.\n")

    # --- traits file for BEAUti ------------------------------------------
    with traits_path.open("w") as fh:
        fh.write("traits\thost\n")
        for k in seqs:
            fh.write(f"{k}\t{info[k][1]}\n")
    note = ("" if (beast_cfg.get("discrete_trait") or {}).get("enabled")
            and not args.no_trait else "   <- import in BEAUti if adding a trait by hand")
    print(f"wrote {traits_path}{note}")

    # --- discrete trait -----------------------------------------------------
    # Previously this was a manual BEAUti step, which meant the model actually
    # run was recorded in neither the config nor git -- and the config said
    # symmetric while the instructions printed below said asymmetric, with no
    # way to tell afterwards which had been clicked.
    trait = None
    dt = beast_cfg.get("discrete_trait") or {}
    if dt.get("enabled") and not args.no_trait:
        # hosts.dta_states maps host group -> trait state, so several groups can
        # collapse onto one state (wild_felid and domestic_cat -> felid). A host
        # group absent from the map has no declared state: it becomes '?', which
        # the codeMap makes ambiguous across all states. Assigning it to a state
        # instead would invent an observation.
        state_map = cfg.get("hosts", {}).get("dta_states") or {}
        tip_states, unmapped = {}, {}
        for k in seqs:
            host = info[k][1]
            if state_map:
                st = state_map.get(host)
            else:
                st = host
            if st is None:
                unmapped[host] = unmapped.get(host, 0) + 1
                st = "?"
            tip_states[k] = st
        if unmapped:
            print(f"NOTE: host groups with no dta_states mapping, entered as "
                  f"ambiguous: {unmapped}")

        observed = sorted({v for v in tip_states.values() if v != "?"})
        declared_all = (sorted(set(state_map.values())) if state_map else observed)
        dropped = [s for s in declared_all if s not in observed]
        if dropped:
            # A state with no tips contributes rates nothing can inform, and
            # asymmetric BSSVS is already rate-hungry.
            print(f"NOTE: trait states with no tips in this alignment, "
                  f"excluded: {dropped}")
        declared = observed
        try:
            trait = trait_blocks(
                prefix, dt.get("trait", "host"), declared, tip_states,
                symmetric=bool(dt.get("symmetric", True)),
                bssvs=bool(dt.get("bssvs", True)),
                poisson_lambda=dt.get("poisson_lambda"))
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1
        sym = "symmetric" if dt.get("symmetric", True) else "asymmetric"
        bs = "BSSVS on" if dt.get("bssvs", True) else "BSSVS off"
        print(f"discrete trait: {trait['n_states']} states, "
              f"{trait['n_rates']} rates ({sym}, {bs})")
        per_rate = len(seqs) / trait["n_rates"]
        if per_rate < 10:
            print(f"  WARNING: {len(seqs)} tips for {trait['n_rates']} rates "
                  f"({per_rate:.1f} tips/rate). Most rates will be "
                  "unidentified; report only those with decisive BF support.")
        print()

    # --- main XML ---------------------------------------------------------
    uncertainties = uncertainty_map(args.metadata, seqs) if date_policy == "interval" else {}
    if date_policy == "interval" and not uncertainties:
        print("NOTE: date policy is 'interval' but no per-tip uncertainty was "
              "found (pass --metadata). Tip dates will be fixed at their "
              "midpoints, which reports an HPD narrower than the data support.")

    xml = build_xml(seqs, info, rate, args.chain, log_every, prefix,
                    tree_prior=tree_prior, clock_model=clock_model,
                    date_policy=date_policy, uncertainties=uncertainties,
                    rate_note=rate_note, trait=trait)
    main_path = xml_path
    main_path.write_text(xml)
    try:
        ET.fromstring(xml)
        print(f"wrote {main_path}   (well-formed XML)")
    except ET.ParseError as exc:
        print(f"wrote {main_path}   BUT XML IS MALFORMED: {exc}", file=sys.stderr)
        return 1

    # --- date randomised replicates ---------------------------------------
    if args.replicate is not None:
        # One replicate, written to --out-xml. The workflow fans these out as
        # separate jobs rather than producing 20 files from one invocation.
        rng = random.Random(args.seed + args.replicate)
        shuffled = list(years)
        rng.shuffle(shuffled)
        info_r = {k: (info[k][0], info[k][1], shuffled[j])
                  for j, k in enumerate(seqs)}
        xml_r = build_xml(
            seqs, info_r, rate, args.chain, log_every,
            f"{prefix}_drt{args.replicate}",
            seed_note=(f"DATE-RANDOMISED REPLICATE {args.replicate} — dates "
                       "shuffled for the temporal signal test. NOT a real analysis."),
            tree_prior=tree_prior, clock_model=clock_model,
            date_policy="midpoint", rate_note=rate_note)
        main_path.write_text(xml_r)
        print(f"wrote date-randomised replicate {args.replicate} -> {main_path}")
        return 0

    if args.randomize_dates:
        rng = random.Random(args.seed)
        drt_dir = (args.out / "date_randomised") if args.out else xml_path.parent
        drt_dir.mkdir(parents=True, exist_ok=True)
        for i in range(1, args.randomize_dates + 1):
            shuffled = list(years)
            rng.shuffle(shuffled)
            info_r = {k: (info[k][0], info[k][1], shuffled[j])
                      for j, k in enumerate(seqs)}
            xml_r = build_xml(seqs, info_r, rate, args.chain, log_every,
                              f"{prefix}_drt{i}",
                              seed_note=f"DATE-RANDOMISED REPLICATE {i} — dates shuffled, "
                                        "for the temporal signal test. NOT a real analysis.",
                              tree_prior=tree_prior, clock_model=clock_model,
                              date_policy="midpoint", rate_note=rate_note)
            (drt_dir / f"{prefix}_drt{i}.xml").write_text(xml_r)
        print(f"wrote {args.randomize_dates} date-randomised replicates -> {drt_dir}")
        print("  Run these, then compare the ucldMean 95% HPDs against the real analysis.")
        print("  Real estimate outside the randomised distribution = genuine temporal signal.")

    trait_steps = ("" if trait is None else f"""
   The discrete trait is already in the XML: {trait['n_states']} states,
   {trait['n_rates']} rates. Nothing to add in BEAUti.
""")
    print(f"""
NEXT STEPS
{trait_steps}
1. Smoke test before committing days of compute:

     sed 's/chainLength="{args.chain}"/chainLength="1000000"/' {main_path.name} > smoke.xml
     beast -seed 1 -threads 2 -overwrite smoke.xml

   Confirm it starts, the likelihood climbs, and — if the trait is in —
   that {prefix}_{(trait or {}).get('tag', 'host')}.trees appears with
   ancestral states on its nodes. Then delete smoke.xml and its outputs.

2. Full run, one directory per seed so nothing overwrites anything:

     mkdir -p seed12345 seed54321
     caffeinate -i beast -seed 12345 -threads 2 -overwrite -prefix seed12345/ {main_path.name}
     caffeinate -i beast -seed 54321 -threads 2 -overwrite -prefix seed54321/ {main_path.name}

3. Convergence: every parameter ESS > 200, and the two chains agreeing with
   each other, not merely each converged on its own.

     loganalyser -b 10 seed12345/{prefix}.log

4. Date randomisation test (--no-trait; only the clock rate matters there),
   then TreeAnnotator for the MCC tree.
""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

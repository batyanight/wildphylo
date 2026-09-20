"""
beastxml.py — model blocks for the BEAST XML, driven by config.

The original 07_make_beast_xml.py baked one model into a single f-string:
HKY+G4, uncorrelated relaxed lognormal clock, constant-size coalescent. Every
choice the pathogen config declares — tree_prior, clock_model,
substitution_model — was ignored, and the clock rate defaulted to a CDV
literature value regardless of which pathogen was being analysed.

This module supplies the pieces that differ, so the XML matches what the config
actually asks for. Each builder returns the XML fragments for one choice, and
the caller assembles them.

WHY THE TREE PRIOR MATTERS HERE. A constant-size coalescent assumes a stable
population sampled at random. The CDV clade containing the 1994 Serengeti
epidemic has 29 of 82 tips from a single outbreak year — a population that
crashed and a sample that is anything but random. Fitting a constant-size prior
to that is not a detail; it is the wrong model, and it biases both the
population size and the dates that depend on it.
"""

from __future__ import annotations


# --- clock rate prior -------------------------------------------------------

def resolve_clock_rate(cfg: dict, measured_slope: float | None,
                       locus: dict | None = None) -> tuple[float, str]:
    """
    (rate, provenance) for the clock prior.

    Config value wins when set, because it is a deliberate declaration with a
    citation attached. When it is null the prior is centred on the root-to-tip
    slope measured from THIS dataset — the right default for a pathogen with no
    published rate, and far better than inheriting one from another virus.
    """
    per_locus = (locus or {}).get("clock_rate_prior")
    if per_locus is not None:
        src = (locus or {}).get("clock_rate_source") or "config (per locus)"
        return float(per_locus), f"config: {src}"

    declared = (cfg.get("beast") or {}).get("clock_rate_prior")
    if declared is not None:
        src = (cfg.get("beast") or {}).get("clock_rate_source") or "config"
        return float(declared), f"config: {src}"

    if measured_slope and measured_slope > 0:
        return float(measured_slope), (
            "measured root-to-tip slope from this dataset "
            "(no clock_rate_prior in config)")

    raise ValueError(
        "no clock rate available: beast.clock_rate_prior is null and no usable "
        "root-to-tip slope was supplied. Run the temporal gate first and pass "
        "its report with --temporal-report, or set a rate in the config.")


# --- tree priors ------------------------------------------------------------

def tree_prior_blocks(kind: str, prefix: str, n_tips: int) -> dict[str, str]:
    """
    XML for one tree prior: state parameters, the prior distribution, its
    operators, and what to log.
    """
    if kind == "coalescent_constant":
        return {
            "state": f'''            <parameter id="popSize.t:{prefix}" spec="parameter.RealParameter" name="stateNode">0.3</parameter>''',
            "prior": f'''                <distribution id="CoalescentConstant.t:{prefix}" spec="Coalescent">
                    <populationModel id="ConstantPopulation.t:{prefix}" spec="ConstantPopulation" popSize="@popSize.t:{prefix}"/>
                    <treeIntervals id="TreeIntervals.t:{prefix}" spec="beast.base.evolution.tree.TreeIntervals" tree="@Tree.t:{prefix}"/>
                </distribution>
                <prior id="PopSizePrior.t:{prefix}" name="distribution" x="@popSize.t:{prefix}">
                    <OneOnX id="OneOnX.0" name="distr"/>
                </prior>''',
            "operators": f'''        <operator id="PopSizeScaler.t:{prefix}" spec="ScaleOperator" parameter="@popSize.t:{prefix}" scaleFactor="0.75" weight="3.0"/>''',
            "log": f'''            <log idref="popSize.t:{prefix}"/>
            <log idref="CoalescentConstant.t:{prefix}"/>''',
        }

    if kind == "coalescent_skyline":
        # Groups trade resolution against noise. 5 is a common default; with
        # very few tips even that over-parameterises, so scale it down.
        groups = max(2, min(5, n_tips // 15))
        return {
            "state": f'''            <parameter id="bPopSizes.t:{prefix}" spec="parameter.RealParameter" dimension="{groups}" lower="0.0" name="stateNode">380.0</parameter>
            <stateNode id="bGroupSizes.t:{prefix}" spec="parameter.IntegerParameter" dimension="{groups}">1</stateNode>''',
            "prior": f'''                <distribution id="BayesianSkyline.t:{prefix}" spec="BayesianSkyline" groupSizes="@bGroupSizes.t:{prefix}" popSizes="@bPopSizes.t:{prefix}">
                    <treeIntervals id="TreeIntervals.t:{prefix}" spec="beast.base.evolution.tree.TreeIntervals" tree="@Tree.t:{prefix}"/>
                </distribution>
                <distribution id="MarkovChainedPopSizes.t:{prefix}" spec="beast.base.inference.distribution.MarkovChainDistribution" jeffreys="true" parameter="@bPopSizes.t:{prefix}"/>''',
            "operators": f'''        <operator id="popSizesScaler.t:{prefix}" spec="ScaleOperator" parameter="@bPopSizes.t:{prefix}" scaleFactor="0.75" weight="15.0"/>
        <operator id="groupSizesDelta.t:{prefix}" spec="DeltaExchangeOperator" integer="true" weight="6.0" intparameter="@bGroupSizes.t:{prefix}"/>''',
            "log": f'''            <log idref="bPopSizes.t:{prefix}"/>
            <log idref="bGroupSizes.t:{prefix}"/>
            <log idref="BayesianSkyline.t:{prefix}"/>''',
        }

    if kind == "birth_death":
        return {
            "state": f'''            <parameter id="BDBirthRate.t:{prefix}" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="10000.0">1.0</parameter>
            <parameter id="BDDeathRate.t:{prefix}" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="1.0">0.5</parameter>''',
            "prior": f'''                <distribution id="BirthDeath.t:{prefix}" spec="beast.base.evolution.speciation.BirthDeathGernhard08Model" birthDiffRate="@BDBirthRate.t:{prefix}" relativeDeathRate="@BDDeathRate.t:{prefix}" tree="@Tree.t:{prefix}"/>
                <prior id="BirthRatePrior.t:{prefix}" name="distribution" x="@BDBirthRate.t:{prefix}">
                    <Uniform id="Uniform.bd" name="distr" upper="1000.0"/>
                </prior>
                <prior id="DeathRatePrior.t:{prefix}" name="distribution" x="@BDDeathRate.t:{prefix}">
                    <Uniform id="Uniform.dr" name="distr"/>
                </prior>''',
            "operators": f'''        <operator id="BirthRateScaler.t:{prefix}" spec="ScaleOperator" parameter="@BDBirthRate.t:{prefix}" scaleFactor="0.75" weight="3.0"/>
        <operator id="DeathRateScaler.t:{prefix}" spec="ScaleOperator" parameter="@BDDeathRate.t:{prefix}" scaleFactor="0.75" weight="3.0"/>''',
            "log": f'''            <log idref="BDBirthRate.t:{prefix}"/>
            <log idref="BDDeathRate.t:{prefix}"/>
            <log idref="BirthDeath.t:{prefix}"/>''',
        }

    raise ValueError(
        f"unknown tree prior {kind!r}; expected coalescent_constant, "
        "coalescent_skyline or birth_death")


# --- clock models -----------------------------------------------------------

def clock_blocks(kind: str, prefix: str, rate: float, n_tips: int) -> dict[str, str]:
    if kind == "strict":
        return {
            "state": f'''            <parameter id="clockRate.c:{prefix}" spec="parameter.RealParameter" name="stateNode">{rate:g}</parameter>''',
            "prior": f'''                <prior id="ClockPrior.c:{prefix}" name="distribution" x="@clockRate.c:{prefix}">
                    <LogNormal id="LogNormalDistributionModel.rate" name="distr" M="{rate:g}" S="1.0" meanInRealSpace="true"/>
                </prior>''',
            "branchrate": f'''                    <branchRateModel id="StrictClock.c:{prefix}" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:{prefix}"/>''',
            "operators": f'''        <operator id="StrictClockScaler.c:{prefix}" spec="ScaleOperator" parameter="@clockRate.c:{prefix}" scaleFactor="0.75" weight="3.0"/>''',
            "log": f'''            <log idref="clockRate.c:{prefix}"/>''',
            "ref": f"StrictClock.c:{prefix}",
        }

    # relaxed lognormal (default)
    return {
        "state": f'''            <parameter id="ucldMean.c:{prefix}" spec="parameter.RealParameter" name="stateNode">{rate:g}</parameter>
            <parameter id="ucldStdev.c:{prefix}" spec="parameter.RealParameter" lower="0.0" name="stateNode">0.1</parameter>
            <stateNode id="rateCategories.c:{prefix}" spec="parameter.IntegerParameter" dimension="{2 * n_tips - 2}">1</stateNode>''',
        "prior": f'''                <prior id="ucldMeanPrior.c:{prefix}" name="distribution" x="@ucldMean.c:{prefix}">
                    <LogNormal id="LogNormalDistributionModel.rate" name="distr" M="{rate:g}" S="1.0" meanInRealSpace="true"/>
                </prior>
                <prior id="ucldStdevPrior.c:{prefix}" name="distribution" x="@ucldStdev.c:{prefix}">
                    <Exponential id="Exponential.stdev" name="distr" mean="0.3337"/>
                </prior>''',
        "branchrate": f'''                    <branchRateModel id="RelaxedClock.c:{prefix}" spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel" clock.rate="@ucldMean.c:{prefix}" rateCategories="@rateCategories.c:{prefix}" tree="@Tree.t:{prefix}">
                        <LogNormal id="LogNormalDistributionModel.c:{prefix}" S="@ucldStdev.c:{prefix}" meanInRealSpace="true" name="distr" M="1.0"/>
                    </branchRateModel>''',
        "operators": f'''        <operator id="ucldMeanScaler.c:{prefix}" spec="ScaleOperator" parameter="@ucldMean.c:{prefix}" scaleFactor="0.5" weight="1.0"/>
        <operator id="ucldStdevScaler.c:{prefix}" spec="ScaleOperator" parameter="@ucldStdev.c:{prefix}" scaleFactor="0.5" weight="3.0"/>
        <operator id="CategoriesRandomWalk.c:{prefix}" spec="IntRandomWalkOperator" parameter="@rateCategories.c:{prefix}" weight="10.0" windowSize="1"/>
        <operator id="CategoriesSwapOperator.c:{prefix}" spec="SwapOperator" intparameter="@rateCategories.c:{prefix}" weight="10.0"/>
        <operator id="CategoriesUniform.c:{prefix}" spec="UniformOperator" parameter="@rateCategories.c:{prefix}" weight="10.0"/>''',
        "log": f'''            <log idref="ucldMean.c:{prefix}"/>
            <log idref="ucldStdev.c:{prefix}"/>''',
        "ref": f"RelaxedClock.c:{prefix}",
    }


# --- sampled tip dates ------------------------------------------------------

def tip_date_blocks(policy: str, prefix: str, uncertainties: dict[str, float],
                    dates: dict[str, float] | None = None) -> dict[str, str]:
    """
    Sample imprecise tip dates WITHIN their uncertainty window.

    Fixing a year-only date at the midpoint asserts precision the record does
    not have and reports a rate HPD narrower than the data support (DR-001).

    THE BOUNDING IS THE POINT. A bare TipDatesRandomWalker with a window size
    lets a tip's date wander without limit — over a long chain a year-only date
    can drift years away from the year actually recorded, which is worse than a
    midpoint because a midpoint is at least wrong in a known, bounded way. The
    first version of this function did exactly that, and BEAST's operator
    tuning gave it away by suggesting a 24-year proposal window for tips whose
    real uncertainty is half a year.

    So each sampled tip gets an MRCAPrior over its own single-taxon set, with a
    Uniform distribution bounded by that tip's window. Heights are measured back
    from the most recent sample, so a date d with uncertainty u becomes a height
    in [mrsd - (d+u), mrsd - (d-u)], clamped at zero because no tip can sit
    above the most recent one.

    The cost is still real: these are extra parameters that mix slowly, so
    convergence must be checked on the tip heights and not only on the rate.
    preflight.check_convergence warns when this policy is set and no tip-height
    ESS appears in the log.
    """
    if policy != "interval" or not uncertainties or not dates:
        return {"state": "", "operators": "", "priors": "", "log": "",
                "n_sampled": 0}

    mrsd = max(dates.values())
    taxonsets, priors, logs = [], [], []
    for lbl in sorted(uncertainties):
        if lbl not in dates:
            continue
        u = float(uncertainties[lbl])
        d = float(dates[lbl])
        lower = max(0.0, mrsd - (d + u))
        upper = max(lower + 1e-6, mrsd - (d - u))
        ident = "".join(ch if ch.isalnum() else "_" for ch in lbl)
        priors.append(
            f'''                <distribution id="tipdate.{ident}" spec="beast.base.evolution.tree.MRCAPrior"
                              tipsonly="true" tree="@Tree.t:{prefix}">
                    <taxonset id="tipset.{ident}" spec="TaxonSet">
                        <taxon id="{lbl}" spec="Taxon"/>
                    </taxonset>
                    <Uniform id="tipunif.{ident}" name="distr" lower="{lower:.6f}" upper="{upper:.6f}"/>
                </distribution>''')
        taxonsets.append(f'                    <taxon idref="{lbl}"/>')  # defined above
        logs.append(f'            <log idref="tipdate.{ident}"/>')

    if not priors:
        return {"state": "", "operators": "", "priors": "", "log": "",
                "n_sampled": 0}

    # Window sized to the typical uncertainty, not left to BEAST's tuner: the
    # bound is enforced by the priors above, and a wildly over-sized proposal
    # just wastes steps being rejected at the boundary.
    window = max(0.05, min(uncertainties.values()))
    rows = "\n".join(taxonsets)
    return {
        "state": "",
        "priors": "\n".join(priors),
        "operators": f'''        <operator id="TipDatesRandomWalker.t:{prefix}" spec="TipDatesRandomWalker" windowSize="{window:.4f}" tree="@Tree.t:{prefix}" weight="3.0">
            <taxonset id="SampledTips.{prefix}" spec="TaxonSet">
{rows}
            </taxonset>
        </operator>''',
        "log": "\n".join(logs),
        "n_sampled": len(priors),
    }

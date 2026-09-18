"""
preflight.py — per-step verification of inputs and settings.

You asked for checks and balances at each step. This is that layer: every rule
in the workflow calls the matching check before it does any work, and the check
either passes, warns, or refuses.

The organising idea is that the expensive errors in phylodynamics are not
crashes. They are runs that complete and produce a plausible number from a
wrong premise: a clock prior from the wrong virus, a trait analysis with more
rates than transitions, an alignment whose reference coordinates silently
shifted. Each check below targets one of those.

Every check returns Checks, which carries three tiers:

    errors    refuse to proceed; this cannot produce a valid result
    warnings  proceed, but the run carries a stated caveat
    notes     recorded in provenance so a future rebuild can be compared

Nothing here is pathogen-specific; thresholds come from the config.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict


@dataclass
class Checks:
    step: str
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def error(self, msg: str) -> "Checks":
        self.errors.append(msg); return self

    def warn(self, msg: str) -> "Checks":
        self.warnings.append(msg); return self

    def note(self, msg: str) -> "Checks":
        self.notes.append(msg); return self

    def to_dict(self) -> dict:
        return asdict(self)

    def render(self) -> str:
        lines = [f"[{self.step}] "
                 f"{'FAIL' if self.errors else 'OK'} "
                 f"({len(self.errors)}E {len(self.warnings)}W)"]
        lines += [f"  ERROR   {e}" for e in self.errors]
        lines += [f"  WARN    {w}" for w in self.warnings]
        lines += [f"  note    {n}" for n in self.notes]
        return "\n".join(lines)


# --- step 1: fetch ----------------------------------------------------------

def check_fetch(cfg: dict, n_records: int, n_previous: int | None = None) -> Checks:
    """
    Verify the download before curating it. An empty or wildly changed pull is
    almost always a query problem, not a discovery.
    """
    c = Checks("fetch")
    if n_records == 0:
        c.error(f"taxid {cfg['fetch']['taxid']} returned no records; check the "
                "taxid and query_filter before assuming the taxon is unsampled")
        return c
    if n_records < 20:
        c.warn(f"only {n_records} records; too few for clade structure or "
               "stratified subsampling to mean much")

    upd = cfg.get("nextstrain", {}).get("update", {})
    if n_previous is not None:
        delta = n_records - n_previous
        c.note(f"{delta:+d} records since the last build ({n_previous} -> {n_records})")
        if delta < 0:
            c.error(f"the dataset LOST {abs(delta)} records. GenBank records "
                    "are rarely withdrawn; this is almost always a changed "
                    "query or a truncated download, not real")
        elif delta > upd.get("max_new_accessions", 10**9):
            c.warn(f"{delta} new records exceeds max_new_accessions "
                   f"({upd['max_new_accessions']}); a bulk submission or a "
                   "changed query deserves eyes before days of compute")
    return c


# --- step 2: curation -------------------------------------------------------

def check_curation(cfg: dict, n_in: int, n_out: int, n_review: int,
                   exclusion_reasons: dict) -> Checks:
    """Verify that curation kept a defensible fraction and dropped for stated reasons."""
    c = Checks("curate")
    if n_out == 0:
        c.error("curation retained no records; check locus.aliases against a "
                "few GenBank records by hand -- an unrecognised gene name is "
                "the usual cause")
        return c
    kept = n_out / n_in if n_in else 0
    c.note(f"retained {n_out}/{n_in} ({kept:.0%})")
    if kept < 0.10:
        c.warn(f"only {kept:.0%} of records survived curation. Review "
               "exclusions.tsv: a dominant single reason usually means an "
               "over-strict filter rather than dirty data")
    top = sorted(exclusion_reasons.items(), key=lambda kv: -kv[1])[:1]
    if top and n_in and top[0][1] / n_in > 0.5:
        c.warn(f"a single exclusion reason ({top[0][0]}) accounts for "
               f"{top[0][1]/n_in:.0%} of all records")
    if n_review:
        c.warn(f"{n_review} records need manual review; their host or date "
               "labels are unverified and they feed the trait analysis directly")
    return c


# --- step 3: host groups and trait states ----------------------------------

def check_trait_states(cfg: dict, group_counts: dict) -> Checks:
    """
    Verify the discrete-trait design BEFORE a multi-day run, not after.

    The rate count is the thing that matters. n states gives n(n-1)/2 symmetric
    rates, and each needs observed transitions to be identifiable. A design
    with more rates than tips is not going to resolve, and no amount of chain
    length fixes it.
    """
    c = Checks("trait_states")
    hosts = cfg.get("hosts", {})
    mapping = hosts.get("dta_states") or {}
    other = hosts.get("dta_other")

    if not mapping:
        states = {g for g, n in group_counts.items() if n > 0}
        c.warn("hosts.dta_states is unset, so every host group becomes a trait "
               "state. Declare the collapse explicitly -- otherwise the number "
               "of states depends on which sequences happened to be subsampled")
    else:
        unmapped = [g for g in group_counts if g not in mapping and group_counts[g] > 0]
        if unmapped:
            if other is None:
                c.note(f"excluded from the trait analysis (unmapped, dta_other "
                       f"is null): {', '.join(sorted(unmapped))}")
            else:
                c.note(f"pooled into '{other}': {', '.join(sorted(unmapped))}")
        bad = [g for g in mapping if g not in group_counts]
        if bad:
            c.warn(f"hosts.dta_states maps groups absent from this dataset: "
                   f"{', '.join(sorted(bad))}")

        counts: dict[str, int] = {}
        for g, n in group_counts.items():
            st = mapping.get(g, other)
            if st:
                counts[st] = counts.get(st, 0) + n
        states = set(counts)

        min_n = hosts.get("min_group_size", 5)
        thin = {s: n for s, n in counts.items() if n < min_n}
        if thin:
            c.warn("trait states below min_group_size "
                   f"({min_n}): {thin}. A state with a handful of tips cannot "
                   "support an identifiable transition rate and will widen "
                   "every other rate's interval")
        if len(counts) < 2:
            c.error("fewer than 2 trait states after collapsing; there is no "
                    "transition to reconstruct")
            return c

    k = len(states)
    symmetric = cfg.get("beast", {}).get("discrete_trait", {}).get("symmetric", True)
    n_rates = k * (k - 1) // 2 if symmetric else k * (k - 1)
    n_tips = sum(group_counts.values())
    c.note(f"{k} trait states -> {n_rates} "
           f"{'symmetric' if symmetric else 'asymmetric'} rates, {n_tips} tips")

    if n_rates > n_tips / 10:
        c.warn(f"{n_rates} rates against {n_tips} tips. As a rule of thumb you "
               "want an order of magnitude more tips than rates; below that, "
               "most rates come back indistinguishable from the mean and the "
               "root state stays unresolved. Collapse further, or keep BSSVS on")
    if not symmetric and n_rates > n_tips / 20:
        c.warn("an asymmetric model doubles the rate count. Directional rates "
               "are worth having, but only where the transitions to estimate "
               "them from actually exist")
    return c


# --- step 4: clock prior ----------------------------------------------------

def check_clock_prior(cfg: dict, rtt_slope: float | None) -> Checks:
    """
    Verify the clock prior against the data and against the genome class.

    The failure this exists to prevent: a rate copied from a paper on a
    different virus, which produces a confident posterior centred on an
    assumption rather than on evidence.
    """
    c = Checks("clock_prior")
    beast = cfg.get("beast", {})
    prior = beast.get("clock_rate_prior")
    vr = (cfg.get("virus") or {}).get("expected_rate_range")

    if prior is None:
        if rtt_slope and rtt_slope > 0:
            c.note(f"clock prior will be centred on the measured root-to-tip "
                   f"slope, {rtt_slope:.3e} subs/site/yr")
        else:
            c.error("clock_rate_prior is null and no usable root-to-tip slope "
                    "is available; there is nothing to centre the prior on")
        return c

    prior = float(prior)
    if not beast.get("clock_rate_source"):
        c.warn("clock_rate_prior is set with no clock_rate_source. An "
               "uncited rate is how one pathogen's assumption ends up in "
               "another's result")

    if vr and not (vr[0] <= prior <= vr[1]):
        c.warn(f"clock_rate_prior {prior:.3e} sits outside the expected range "
               f"for this genome class ({vr[0]:.1e}-{vr[1]:.1e}). Defensible if "
               "you have a reason, but state it in clock_rate_source")

    if rtt_slope and rtt_slope > 0:
        ratio = prior / rtt_slope
        c.note(f"prior {prior:.3e} vs measured root-to-tip slope "
               f"{rtt_slope:.3e} (ratio {ratio:.2f})")
        if ratio > 10 or ratio < 0.1:
            c.error(f"the prior disagrees with the data by {ratio:.1f}x. Either "
                    "the prior belongs to a different pathogen or locus, or the "
                    "regression is being driven by something other than time. "
                    "Resolve this before committing compute")
        elif ratio > 3 or ratio < 0.33:
            c.warn(f"the prior and the measured slope differ by {ratio:.1f}x. A "
                   "relaxed clock can absorb this, but the posterior will be "
                   "prior-influenced -- check it against a prior-only run")
    return c


# --- step 5: convergence ----------------------------------------------------

def check_convergence(cfg: dict, ess: dict, n_chains: int,
                      between_chain_overlap: dict | None = None) -> Checks:
    """Verify the MCMC before any number from it is reported."""
    c = Checks("convergence")
    min_ess = cfg.get("beast", {}).get("min_ess", 200)
    low = {k: v for k, v in ess.items() if v < min_ess}
    if low:
        c.error(f"ESS below {min_ess} for: {low}. These parameters are not "
                "adequately sampled and their intervals are not trustworthy")
    else:
        c.note(f"all {len(ess)} parameters have ESS >= {min_ess}")

    if n_chains < 2:
        c.warn("a single chain can demonstrate non-convergence but never "
               "convergence")
    if between_chain_overlap:
        bad = [k for k, ov in between_chain_overlap.items() if ov < 0.5]
        if bad:
            c.error(f"chains disagree on: {bad}. Their marginal posteriors "
                    "overlap by less than half, which means at least one chain "
                    "has not found the same region of parameter space")

    policy = cfg.get("dates", {}).get("imprecise_date_policy")
    if policy == "interval" and not any("height" in k or "date" in k.lower()
                                        for k in ess):
        c.warn("imprecise_date_policy is 'interval', so tip dates are sampled "
               "parameters, but no tip-height ESS was found in the log. Sampled "
               "tip dates mix poorly and must be checked, not assumed")
    return c

# Decision records

Internal justification for choices the pipeline makes on your behalf. Each
records what was decided, why, what the alternative costs, and how the decision
should be revisited per dataset rather than inherited.

---

## DR-001 — Imprecise date handling: `interval` for CDV, decided per dataset thereafter

**Status:** accepted
**Applies to:** `dates.imprecise_date_policy`
**Decided by:** `scripts/lib/datepolicy.py`, run per dataset

### Decision

CDV clade 3 moves to `interval`. Other datasets are decided by the rule below,
not by inheriting this choice.

### Why this is not a global default

The usual advice is "use interval sampling when dates are imprecise". That is
correct in principle and unhelpful in practice, because interval sampling costs
chain length and mixing, and on many datasets it changes nothing. The honest
question is not *is midpoint imperfect* — it always is — but *can this
alignment tell the difference*.

Midpoint fails in two unrelated ways, which need separate tests because a
dataset can have one without the other.

#### Failure mode 1 — bias

Midpoint assumes true collection dates are uniform within the year. Wildlife
sampling rarely is: carcass recovery, outbreak response and field-season effort
all cluster seasonally. If real dates concentrate in autumn, every year-only tip
is placed months too early. That error is **directional**, so it does not
average out across tips — it accumulates into the rate and the TMRCA.

Testable, because the subset that *does* carry day-precision dates samples the
same seasonal process. Rayleigh test for circular uniformity on their
within-year positions.

**Measured on CDV clade 3** (n = 36 precise-dated tips):

| Quantity | Value |
|---|---|
| Rayleigh R | 0.274 |
| Rayleigh p | 0.067 |
| Mean within-year position | 0.541 |
| Implied systematic shift | **+0.041 yr (~15 days)** |

Non-significant, and the shift is small. The month histogram shows why: the
distribution is **bimodal** — a December–April peak and a June–July hole,
consistent with CDV outbreak seasonality in raccoons. Two modes roughly six
months apart largely cancel at the midpoint. So CDV sampling is genuinely
seasonal but *not* midpoint-biasing, and the bias test correctly returns null.

This is deliberate behaviour, not a miss: the test asks about bias, not about
seasonality. A unimodal autumn season would fire it.

#### Failure mode 2 — false precision

Even with zero bias, midpoint asserts a date the data did not supply. Whether
that matters depends on whether the alignment could have resolved it:

```
lambda      = clock_rate x alignment_length     substitutions per genome per year
resolvable  = lambda x window_width             substitutions across the window
```

- `resolvable << 1` — no amount of sequence data distinguishes January from
  December. The tip date is unidentifiable, and interval sampling adds a
  badly-mixing parameter for nothing. **Use midpoint.**
- `resolvable >~ 1` — the alignment carries roughly a substitution's worth of
  signal across the window. Midpoint discards it *and* reports a rate HPD
  narrower than the data support. **Use interval.**

**Measured on CDV clade 3:**

```
lambda     = 8.54e-4 x 1824  = 1.56 subs/genome/year
resolvable = 1.56 x 1.0 yr   = 1.56 substitutions
```

Comfortably above 1. This is what decides CDV, not the bias test.

### The decision rule, in order

1. No imprecise tips → `midpoint`, trivially.
2. Imprecise tips < 10% and no seasonal bias → `midpoint`; the choice has no
   leverage.
3. Significant seasonal bias with |shift| > 0.08 yr → `interval`, regardless of
   resolvability. A directional error propagates into the rate.
4. `resolvable` < 0.3 → `midpoint`. Tip dates are unidentifiable.
5. `resolvable` >= 1.0 and >= 10% imprecise → `interval`.
6. Anything between → `interval`, marginal confidence. Rationale: interval's
   failure mode is intervals that are too wide, which is visible. Midpoint's
   is intervals that are too narrow, which is not.

### Pitfalls, for citation in methods

**`midpoint`**
- Understates uncertainty in rate and TMRCA; HPDs are too narrow.
- Directionally biased under unimodal seasonal sampling.
- Cheap, stable, and comparable with most published analyses.

**`interval`**
- Wider, more honest HPDs.
- More parameters; needs longer chains. **ESS must be checked on the sampled
  tip heights, not just on the rate** — enforced by
  `preflight.check_convergence`, which warns if `interval` is set and no
  tip-height ESS appears in the log.
- Tip-date operators interact badly with strict clocks on shallow trees; mixing
  can collapse. CDV uses a relaxed lognormal clock, so this is not a concern
  here.
- **Not comparable with a previously published midpoint analysis.** This
  matters for you specifically: moving clade 3 to `interval` will likely widen
  the TMRCA HPD beyond the published 1962–1985. That is an improvement in
  honesty, not a new finding, and should be described as such.

**`exclude`**
- Only defensible when imprecise tips are a small minority *and* dropping them
  does not gut a host group or a time period. Almost never right for
  GenBank-derived wildlife datasets, where the imprecise tips are usually the
  **old** ones carrying the temporal span. For CDV, excluding year-only dates
  would remove 126 of 162 sequences including most of the pre-2010 record —
  collapsing the very span the TMRCA depends on.

### Known weakness of the bias test

It assumes date precision is **missing at random**. If the precise dates came
disproportionately from one submitter or one field season, the test describes
that subset rather than the sampling process. `n_precise` is always reported and
the caveat is always attached — `test_missing_at_random_caveat_always_present_when_bias_tested`
pins that it can never go silent.

For CDV, 36 of 162 precise dates is a thin basis. Worth checking whether they
cluster by submitting lab before leaning hard on the null result.

### How to revisit per run

The workflow runs `datepolicy.recommend()` after the ML tree exists (so a
measured root-to-tip slope is available) and writes its full reasoning to
`builds/<pathogen>/<date>/temporal/date_policy.json`. If the recommendation
disagrees with the config, the run warns rather than silently overriding —
the config is your declared intent, and the pipeline's job is to tell you when
the data disagrees with it, not to quietly win the argument.

---

## DR-002 — Discrete-trait states are declared, not emergent

**Status:** accepted
**Applies to:** `hosts.dta_states`, `hosts.dta_other`

### Decision

The collapse from host groups to trait states is an explicit map in config.
Groups absent from the map are excluded (`dta_other: null`) or pooled into a
named state.

### Why

The host table defines 12 functional groups. Under a symmetric model that is
`12 x 11 / 2 = 66` transition rates. The published clade-3 run used 5 states —
but that number was an emergent consequence of which sequences landed in the
subset, not a decision anyone made and could defend.

The arithmetic predicts the published outcome. 66 rates against 162 tips is
roughly 2.5 tips per rate; a workable rule of thumb wants an order of magnitude
more tips than rates. Below that, most rates come back indistinguishable from
the mean and the root state stays unresolved — which is precisely what the
README reports: eight of ten rates not individually distinguishable, root state
p = 0.25.

`preflight.check_trait_states` now computes the rate count before the run and
warns when it exceeds tips/10.

### Consequence worth noting

Collapsing can *fix* a thin state rather than merely hiding it. `domestic_cat`
(n=2) and `wild_felid` (n=7) both map to `felid`, giving 9 — above
`min_group_size`. Thinness is therefore evaluated **after** collapsing, not
before, or the pipeline would warn about a problem the map already solved.
Pinned by `test_thinness_is_judged_after_collapsing_not_before`.

### Pitfall

A collapse is a scientific claim: it asserts those hosts are epidemiologically
interchangeable for this pathogen. Pooling `domestic_cat` with `wild_felid` says
captive and free-ranging felids play the same role in transmission. Defensible
for the two felid outbreak clusters you recovered, both captive. Less defensible
if the question ever becomes about free-ranging felids specifically. The map
makes that claim visible and reviewable, which is the point.

---

## DR-003 — Snakemake over Nextflow

**Status:** accepted, reversible

Your scripts are Python and pass state through DataFrames. Snakemake keeps that
in one language and one process model. More decisively, its `checkpoint`
mechanism directly expresses the thing this pipeline most needs: a gate that
decides **at runtime** whether the expensive Bayesian branch is scheduled at
all. Verified — the dry run schedules 8 jobs, because nothing past the temporal
gate is planned until the gate resolves.

Nextflow would be the better call on a shared HPC scheduler or across multiple
languages. Revisit if either becomes true.

---

## DR-004 — Rebuilds are cron-triggered but accession-gated

**Status:** accepted
**Applies to:** `nextstrain.update.trigger`

Cron fires on schedule; the build proceeds only if GenBank has at least
`min_new_accessions` new records. Re-analysing an unchanged dataset burns days
of compute and produces a second, confusingly-dated copy of the same result.

`max_new_accessions` is the other half. A large jump is as suspicious as none —
it usually means the query changed or a bulk submission landed, both of which
warrant eyes before compute. A **negative** delta is a hard error:
GenBank rarely withdraws records, so a shrinking dataset is a query bug or a
truncated download, not a finding.

---

## DR-005 — Date policy is per-segment, not per-pathogen

**Status:** accepted
**Applies to:** `loci[].imprecise_date_policy`

### Decision

For segmented pathogens the date policy is set per locus. The pathogen-level
`dates.imprecise_date_policy` is a fallback only.

### Why

Running the DR-001 rule across all ten BTV segments, with published rates and
RefSeq segment lengths, gives different answers within the same virus:

| Seg | Protein | nt | λ (subs/genome/yr) | resolvable | policy |
|---|---|---|---|---|---|
| 1 | VP1 | 3944 | 1.97 | 1.97 | interval, strong |
| 2 | VP2 | 2926 | 1.46 | 1.46 | interval, strong |
| 3 | VP3 | 2772 | 1.39 | 1.39 | interval, strong |
| 4 | VP4 | 1981 | 0.99 | 0.99 | interval, marginal |
| 6 | VP5 | 1638 | 0.82 | 0.82 | interval, marginal |
| 7 | VP7 | 1156 | 0.58 | 0.58 | interval, marginal |
| 10 | NS3 | 822 | 0.35 | 0.35 | interval, marginal |
| — | CDV H | 1824 | 1.56 | 1.56 | interval, strong |

Segment length varies nearly fivefold, and rate varies too, so λ spans an order
of magnitude within one genome. Seg-10 at 0.35 sits just above the `midpoint`
floor of 0.3 — the shortest segment can barely resolve within-year timing at
all. A single pathogen-level policy would be wrong for most of the genome.

### Pitfall

`interval / marginal` on Seg-10 means the tip-date parameters may not mix.
`preflight.check_convergence` already warns when `interval` is set and no
tip-height ESS appears in the log. If Seg-10's tip heights come back with poor
ESS, drop that segment to `midpoint` and record it — a badly-mixing sampled
parameter is worse than an acknowledged approximation.

---

## DR-006 — The vector is annotated, not treated as a host state

**Status:** accepted
**Applies to:** `hosts.vector_policy`

### Decision

BTV sequences derived from *Culicoides* are retained in the phylogeny and the
metadata, but excluded from the discrete-trait host analysis. Default
`annotate`.

### Why

Every BTV host-to-host transition physically passes through a midge, so the
vector is not a host state in the sense that sheep or white-tailed deer are. A
midge sequence records **where the virus was sampled**, not where the lineage is
maintained. Treating `vector` as a discrete state would make every transition
rate into and out of it a description of entomological surveillance effort
rather than of transmission.

Alternatives and their costs:

- `exclude` — drops real data. Vector sequences are often the only evidence of
  local circulation in a region with no ruminant sampling, so discarding them
  loses geographic signal.
- `state` — simple, and the rates are uninterpretable.
- `annotate` — sequences inform tree topology and branch lengths, but
  contribute no host transitions they cannot speak to.

### Pitfall

`annotate` is not free either. Vector-derived tips still affect the coalescent
and the clock, and if vector sampling is concentrated in one country or one
season it will skew the tree even while contributing no trait states. Report the
count of vector-derived tips alongside the host composition.

### Note for CDV

This has no CDV analogue — direct transmission means the host label *is* the
maintenance host. The distinction only appears with a vector-borne pathogen,
which is part of why BTV was worth doing second.

---

## DR-007 — Segment congruence is a gate, not a diagnostic

**Status:** accepted
**Applies to:** `segments.concatenate`, `segments.congruence_gate`

### Decision

`concatenate: never` for BTV. The congruence screen runs regardless, and
`concatenate: always` combined with `virus.reassortment_expected: true` is a
**config error**, not a warning.

### Why

Three separate failures follow from concatenating reassorting segments:

1. The concatenated tree matches no segment's actual ancestry — it averages
   over conflicting histories.
2. Segments have genuinely different TMRCAs. BTV Seg-10's diversity is much
   younger than Seg-6's; one clock on a concatenate misdates both.
3. **The one that matters most here.** A reassortant taxon sits among different
   relatives depending on the segment. In a discrete-trait host analysis, that
   single taxon generates a host transition in every segment tree where its
   placement moved — and the transition reflects reassortment, not
   transmission. Given that host-jump inference is the whole point of this
   pipeline, an undetected reassortant is not a minor artefact; it is a false
   positive in the primary result.

`segments.py` therefore does two things pairwise RF alone cannot: it decides
whether concatenation is legitimate, and it **localises** incongruence to
individual taxa via neighbour-set instability, so those tips can be dropped from
the trait analysis. `exclude_reassortants_from_dta: true` is the default.

### Pitfalls

- Robinson-Foulds is blunt: it weights a poorly-supported rearrangement the
  same as a well-supported one. `collapse_support: 70` ignores weak splits
  before comparing, and the threshold is always reported, because RF is
  sensitive to it.
- Short segments lose resolution rather than gaining real conflict. Seg-10 at
  822 nt will look more incongruent than it is.
- This is a **screen, not a test**. A formal treatment needs a
  coalescent-with-reassortment model. What this prevents is a naive
  concatenated analysis running unnoticed, which is the common failure.

### Test-design note

The first version of the congruence tests used six-taxon trees, which have only
three informative splits — so even a complete clade swap topped out at RF 0.33,
below the 0.40 decoupling threshold. Tests were rewritten on twelve taxa. Had
the fixtures been left small, the fix would have been to lower the threshold,
which would have been wrong for real data. Worth remembering: a failing test on
a toy fixture is sometimes telling you the fixture is unrealistic, not that the
threshold is.

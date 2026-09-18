# Managing many Nextstrain builds

You have one build now. You will shortly have several per pathogen — BTV alone
wants one per analysed segment — and then several pathogens. The failure mode
is not technical. It is that in eighteen months you have fourteen datasets, four
of them stale, two contradicting each other, and no way to tell which is which
from the URL.

This is what to decide before that happens rather than after.

---

## 1. The single most important decision: one repo per pathogen

Nextstrain's own convention is **one GitHub repository per pathogen**, and you
should follow it even though it means this pipeline repo is not itself a
Nextstrain pathogen repo.

The split:

```
wildphylo/                   the pipeline. Config-driven, pathogen-agnostic.
                             Installed as a dependency. No Auspice JSONs.

cdv-phylodynamics/           one pathogen. Its config, its auspice/, its
                             provenance, its README, its DOI.
btv-phylodynamics/           another.
```

Why not one monorepo with everything in it:

- **Community URLs are per-repo.** `nextstrain.org/community/batyanight/<repo>`
  — every dataset in a monorepo shares one namespace, and the URL cannot tell
  a reader which pathogen they are looking at without reading the trailing
  segments.
- **Citation.** A DOI on a monorepo cites "my pipeline work". A DOI on
  `btv-phylodynamics` cites a specific analysis with a specific dataset, which
  is what a paper needs.
- **Blast radius.** A change to the pipeline should not force a version bump on
  every pathogen analysis you have ever published.
- **Repo size.** Auspice JSONs are hundreds of kB and are regenerated monthly.
  In a monorepo with twelve pathogens, git history becomes unusable within a
  year.

The cost is real: a pipeline fix has to be pulled into each pathogen repo. Pin
the pipeline version in each pathogen repo's config so the lag is explicit
rather than accidental:

```yaml
pipeline:
  version: "0.3.1"     # wildphylo release this build was produced with
```

---

## 2. Dataset naming determines your URLs

For Nextstrain Community, a repo `batyanight/btv-phylodynamics` serves:

| File in `auspice/` | URL |
|---|---|
| `btv-phylodynamics.json` | `/community/batyanight/btv-phylodynamics` |
| `btv-phylodynamics_seg2.json` | `…/btv-phylodynamics/seg2` |
| `btv-phylodynamics_seg2_wildlife.json` | `…/btv-phylodynamics/seg2/wildlife` |

The repo name prefix is mandatory, and underscores become path separators. So
**the filename is the URL**, and naming is not cosmetic.

Rules that will save you later:

- **Never put a date in the dataset name.** `btv_seg2_2026-03.json` creates a
  new URL every month, so every link you have ever shared points at a frozen
  build. Keep one stable URL per *analysis*, and archive dated copies under
  `builds/` where they are not served.
- **Order segments from general to specific.** `seg2/wildlife` not
  `wildlife/seg2` — Auspice's dataset selector nests on the path, and the first
  segment should be the thing a reader chooses first.
- **Keep the default dataset meaningful.** The bare repo URL should resolve to
  something, not 404. For BTV, make it the segment you would show someone who
  arrived with no context.
- **Build names must be `[A-Za-z0-9_-]`.** Periods break the path. There is a
  CI test enforcing this.

---

## 3. Every dataset needs a description, and it should be generated

Auspice renders `meta.description` as Markdown under the tree. This is the only
place a reader finds out what they are looking at, and it is the first thing
that goes stale.

Generate it from the build rather than writing it by hand:

```
This dataset was built on 2026-03-01 from 214 sequences (GenBank, accessions
listed in the repository).

Temporal signal: PASS (R² = 0.43, date-shuffling p = 0.001)
Clock prior: centred on the measured root-to-tip slope, 5.1 × 10⁻⁴
Date policy: interval (126/214 dates are year-only)
Trait states: 5, from 12 host groups (see config)
Convergence: all parameters ESS > 200, 2 chains

Known limitations: <from config>
Pipeline: wildphylo 0.3.1 · Config: btv.yaml (sha256 a3628b…)
```

A hand-written description drifts from the build within two rebuilds. A
generated one cannot. This matters more than it sounds: when someone emails you
in a year asking whether a conclusion still holds, the description is what tells
you which build they are reading.

---

## 4. Separate "published" from "exploratory" early

You will build things that are not ready to be seen. Three tiers:

| Tier | Where | Who sees it |
|---|---|---|
| Exploratory | `builds/<pathogen>/<date>/` only, never in `auspice/` | you |
| Staging | a `staging` branch, URL `/community/user/repo@staging` | collaborators you send the link to |
| Published | `main`, URL `/community/user/repo` | anyone |

Nextstrain Community serves any branch via `@branch`, which makes staging free.
Use it. The alternative — publishing to `main` and hoping nobody looks yet — is
how a half-finished analysis gets cited.

If you ever need genuinely private datasets, that is Nextstrain Groups, which is
a paid/managed offering. Community is public by construction: **do not put
anything embargoed in a Community repo**, including in a branch, since the repo
is public and branches are visible.

---

## 5. `auspice_config.json` per build, not one shared file

Each dataset gets its own Auspice config, because colourings differ. BTV Seg-2
should colour by serotype; Seg-10 should not, because serotype is a Seg-2
property and colouring a Seg-10 tree by it implies a correspondence that
reassortment breaks.

```
config/
  auspice/
    btv_seg2.json          colour_by: serotype, host_group, country
    btv_seg10.json         colour_by: host_group, country      (no serotype)
    cdv.json               colour_by: host_group, country, lineage
```

Things to set once and keep consistent across every build:

- **`colorings` with explicit `title`.** Auspice falls back to raw field names,
  so a reader sees `host_group` instead of "Host group".
- **A fixed colour map per host group**, shared across pathogens where the
  groups overlap. If `wild_cervid` is green in one dataset and orange in
  another, side-by-side comparison becomes actively misleading.
- **`filters`.** Without them a reader cannot subset by host or country.
- **`display_defaults`.** Decide what someone sees on arrival. Usually
  `colorBy: host_group`, `mapTriplicate: false`, `distanceMeasure: num_date`.
- **`build_url`.** Points the sidebar at the repo that produced the dataset.
  Cheap, and it is how someone finds your methods.

---

## 6. Archive every build; serve only the current one

```
auspice/                      served. current builds only. stable URLs.
builds/<pathogen>/<date>/     archived. never served. one per rebuild.
```

Keep `provenance.json`, the accession list, the temporal-signal report and the
QC table for every archived build — they are small. Drop the alignments, trees
and BEAST logs, which are large and regenerable.

`nextstrain.update.keep_builds` prunes the archive. Twelve is a reasonable
default: a year of monthly rebuilds.

The reason to archive at all: when a conclusion changes between rebuilds, the
only way to find out why is to diff the inputs. Without the accession lists you
cannot.

---

## 7. Compare every rebuild against the last one before publishing

A rebuild that silently changes a published TMRCA is the worst failure mode
here, because nothing errors — you just quietly start telling people something
different.

The `publish` rule refuses when:

- the dataset **lost** sequences (GenBank rarely withdraws records; this is a
  query bug)
- the TMRCA moved more than `alert_on_tmrca_shift_years`
- the temporal-signal gate changed verdict
- the number of trait states changed

None of these mean the new build is wrong. They mean a human should look before
it goes out. `min_new_accessions` / `max_new_accessions` gate the other end: too
few new records and the rebuild is pointless, too many and something changed
about the query or a bulk submission landed.

---

## 8. Narratives, when you have something to say

A narrative is a Markdown file that walks a reader through a dataset with the
view pinned at each step. `narratives/<repo>_<name>.md` in the same repo, served
at `/community/narratives/user/repo/name`.

This is the right home for the felid-outbreak story in your CDV work — the two
clusters 23 years apart, each matching a published investigation. That is
currently prose in a README describing a tree nobody is looking at while they
read it. As a narrative, the tree moves as the text does.

Worth knowing: narratives pin a dataset URL in their frontmatter. If you rename
a dataset, every narrative pointing at it breaks silently. One more reason
dataset names should be stable.

---

## 9. What to do right now, in order

1. Create `wildphylo` as the pipeline repo — this scaffold.
2. Keep `cdv-phylodynamics` as the CDV pathogen repo. Strip the pipeline code
   out of it once `wildphylo` is installable; leave the config, `auspice/`, the
   provenance, the DOI and the README. The existing Community URL keeps
   working, which matters because it is already published.
3. Create `btv-phylodynamics` when the BTV config is ready to run.
4. Set the shared host-group colour map before the second pathogen, not after.
   Retrofitting it means republishing everything.
5. Decide the CDV default dataset URL now, and do not change it again.

---

## Pitfalls specific to segmented pathogens

- **Do not publish a concatenated BTV tree**, even as a convenience view. A
  reader will treat it as the phylogeny, and it is not one. If congruence ever
  passes for a subset of segments, publish *that subset* with the RF matrix in
  the description.
- **One dataset per analysed segment**, each with its own temporal-signal
  verdict in its description. Segments can legitimately disagree, and a reader
  needs to see that rather than being shown whichever segment passed.
- **Do not colour non-Seg-2 trees by serotype.** It implies a correspondence
  that reassortment specifically breaks.
- **Flag reassortant tips in the metadata** and expose them as an Auspice
  filter, so a reader can remove them and see whether a host transition
  survives.

#!/usr/bin/env python3
"""
08e_trait_signal.py — tip-label randomisation gate for discrete-trait analyses.

Run this before reporting ANY number out of a discrete-trait reconstruction.
A DTA always returns a transition table; this decides whether that table
describes the phylogeny or the label distribution.

    python scripts/08e_trait_signal.py \\
        --trees builds/cdv-1200/20260919/clades/america2/seed*/america2_host_group.trees \\
        --config config/pathogen/cdv-1200.yaml \\
        --burnin 0.10 --max-trees 1000 --permutations 1000 \\
        --out-json builds/.../america2/trait_signal.json

Tip states come, in order of preference, from --traits, then from the tip
label field given by --label-field (labels here are ACCESSION|group|year), and
are then collapsed through hosts.dta_states from --config so that the states
tested are exactly the states BEAST estimated rates between. Testing the
uncollapsed host groups would answer a different question from the one the
analysis asked.

The statistics use topology and tip labels only, never BEAST's reconstructed
ancestral states — see scripts/lib/traitsignal.py for why that is deliberate.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib import traitsignal as tsig  # noqa: E402


def states_from_traits(path: Path) -> dict[str, str]:
    """Read the traits file 07 writes: a two-column TSV with a header."""
    out: dict[str, str] = {}
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    for row in rows:
        if len(row) < 2 or row[0].lower() in ("traits", "name", "strain", "taxon"):
            continue
        out[row[0].strip()] = row[1].strip()
    if not out:
        raise ValueError(f"{path} yielded no taxon/state pairs")
    return out


def states_from_labels(labels, field: int, sep: str) -> dict[str, str]:
    out, bad = {}, []
    for lab in labels:
        parts = lab.split(sep)
        if len(parts) > field:
            out[lab] = parts[field].strip()
        else:
            bad.append(lab)
    if bad:
        raise ValueError(
            f"{len(bad)} tip labels have no field {field} when split on "
            f"{sep!r}; first is {bad[0]!r}. Use --traits instead, or set "
            f"--label-field / --label-sep to match your labels.")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", nargs="+", required=True,
                    help="Posterior tree files, one per chain. Burn-in is "
                         "applied within each file, not to the pooled set.")
    ap.add_argument("--traits", type=Path,
                    help="Two-column TSV of taxon and state (what 07 writes).")
    ap.add_argument("--label-field", type=int, default=1,
                    help="Field of the tip label holding the state, 0-based. "
                         "Used when --traits is not given. Default 1, which "
                         "is the group in ACCESSION|group|year.")
    ap.add_argument("--label-sep", default="|")
    ap.add_argument("--config", type=Path,
                    help="Pathogen config. hosts.dta_states is applied to "
                         "collapse host groups onto the trait states BEAST "
                         "actually estimated rates between.")
    ap.add_argument("--burnin", type=float, default=0.10)
    ap.add_argument("--thin", type=int, default=1)
    ap.add_argument("--max-trees", type=int, default=1000,
                    help="Cap on retained trees after burn-in and thinning. "
                         "The statistics converge well before a full "
                         "posterior; more trees buy precision on the mean, "
                         "not resolution on p, which --permutations sets.")
    ap.add_argument("--permutations", type=int, default=1000,
                    help="Null replicates. The smallest p obtainable is "
                         "1/(permutations+1).")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--out-json", type=Path)
    ap.add_argument("--out-md", type=Path)
    ap.add_argument("--exit-nonzero-on-fail", action="store_true",
                    help="Exit 1 when the verdict is fail, so a workflow "
                         "stops rather than publishing an unsupported table.")
    a = ap.parse_args()

    paths = [Path(p) for p in a.trees]
    missing = [p for p in paths if not p.is_file()]
    if missing:
        print(f"ERROR: no such tree file: {missing[0]}", file=sys.stderr)
        return 1

    trees = tsig.read_posterior(paths, burnin=a.burnin, thin=a.thin,
                               max_trees=a.max_trees)
    print(f"read {len(trees)} trees from {len(paths)} file(s) "
          f"(burnin {a.burnin:.0%}, thin {a.thin}, cap {a.max_trees})")

    labels = sorted(trees[0].tip_index)
    if a.traits:
        states = states_from_traits(a.traits)
        src = str(a.traits)
    else:
        states = states_from_labels(labels, a.label_field, a.label_sep)
        src = f"tip label field {a.label_field}"
    print(f"tip states from {src}")

    # Collapse through hosts.dta_states so the tested states match the model's.
    # A host group absent from the map becomes '?', exactly as 07 does it, and
    # those tips are pruned rather than made into a state of their own.
    if a.config:
        import yaml
        cfg = yaml.safe_load(a.config.read_text())
        dta = (cfg.get("hosts") or {}).get("dta_states") or {}
        if dta:
            before = sorted(set(states.values()))
            states = {k: dta.get(v, "?") for k, v in states.items()}
            after = sorted(set(states.values()) - {"?"})
            print(f"collapsed via hosts.dta_states: "
                  f"{len(before)} group(s) -> {len(after)} state(s) {after}")
        else:
            print("config has no hosts.dta_states; using states as read")

    res = tsig.trait_signal(trees, states,
                            n_permutations=a.permutations, seed=a.seed,
                            ai_threshold=a.alpha, ps_threshold=a.alpha)
    report = tsig.format_report(res)
    print()
    print(report)

    if a.out_json:
        a.out_json.parent.mkdir(parents=True, exist_ok=True)
        a.out_json.write_text(json.dumps(res.to_dict(), indent=2))
        print(f"\nwrote {a.out_json}")
    if a.out_md:
        a.out_md.parent.mkdir(parents=True, exist_ok=True)
        a.out_md.write_text("```\n" + report + "\n```\n")
        print(f"wrote {a.out_md}")

    if res.verdict == "fail" and a.exit_nonzero_on_fail:
        print("\nexiting non-zero: the trait carries no demonstrated signal.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

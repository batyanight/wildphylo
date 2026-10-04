#!/usr/bin/env python3
"""
12_compare_builds.py — what changed since the last build, and is it a finding?

    # after a full build, with an Auspice JSON to compare
    python scripts/12_compare_builds.py --current auspice/cdv.json \\
        --builds-dir builds/cdv-1200 --max-tmrca-shift 5.0 --out comparison.md

    # in CI, which stops at the temporal gate and has no Auspice build
    python scripts/12_compare_builds.py --pathogen cdv-1200 \\
        --builds-dir builds/cdv-1200 --out comparison.md

Why this exists
---------------
A scheduled rebuild produces a new tree every month, and a new tree always
differs from the old one. Without a comparison step the only signal available
is whether the pipeline exited zero, and a pipeline exits zero just as happily
when a curation rule change has silently dropped two hundred sequences or when
the TMRCA has moved forty years.

The framing that matters: a large change between builds is a CURATION
REGRESSION until shown otherwise. New data arriving monthly moves a TMRCA by
months, not decades. A dataset that shrinks has almost always lost records to a
changed query rather than to GenBank withdrawals. Treating these as discoveries
is how a pipeline publishes its own bugs.

What is compared
----------------
Whatever both builds have, because CI and a full local run leave different
artefacts behind:

  accession lists   count, and WHICH accessions came and went. A count can be
                    unchanged while the composition turns over completely, and
                    that is invisible to a count-only check.
  temporal signal   the gate verdict, R^2 and the clock slope. A verdict that
                    drops from pass to fail between builds is the single most
                    informative regression available, because everything
                    downstream is conditioned on it.
  gate verdicts     convergence, trait signal, date randomisation, per clade.
  TMRCA             from the Auspice JSON when there is one, against
                    alert_on_tmrca_shift_years.

Deliberately NOT compared: tree topology. Two trees built from overlapping data
differ in ways that are mostly reordering of equally-supported splits, so a
topological distance between builds is large, noisy, and uninformative about
whether anything is wrong. The quantities above are the ones where a change has
an interpretation.

Exit code is non-zero when something regressed, so CI opens an issue rather
than publishing quietly.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path


@dataclass
class Finding:
    level: str            # ok | note | warn | error
    area: str
    message: str


@dataclass
class Comparison:
    pathogen: str
    previous: str | None = None
    current: str | None = None
    findings: list[Finding] = field(default_factory=list)

    def add(self, level, area, message):
        self.findings.append(Finding(level, area, message))

    @property
    def verdict(self) -> str:
        if any(f.level == "error" for f in self.findings):
            return "regression"
        if any(f.level == "warn" for f in self.findings):
            return "review"
        return "ok"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["verdict"] = self.verdict
        return d


def dated_builds(root: Path) -> list[Path]:
    """Dated build directories, oldest first. Names are YYYYMMDD."""
    if not root.is_dir():
        return []
    return sorted((d for d in root.iterdir() if d.is_dir() and d.name.isdigit()),
                  key=lambda d: d.name)


def read_accessions(build: Path) -> set[str] | None:
    hits = sorted(build.glob("raw/*.acc"))
    if not hits:
        return None
    return {l.strip() for l in hits[0].read_text().splitlines() if l.strip()}


def read_json(path: Path):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def read_temporal(build: Path) -> dict | None:
    for p in sorted(build.glob("temporal/*_temporal_signal.json")):
        doc = read_json(p)
        if doc:
            return doc
    return None


def gate_verdicts(build: Path) -> dict[str, str]:
    """Every gate verdict in a build, keyed by a readable name."""
    out: dict[str, str] = {}
    patterns = {
        "temporal": "temporal/*_temporal_signal.json",
        "convergence": "beast/*/convergence.json",
        "trait_signal": "clades/*/trait_signal.json",
        "drt": "drt/*/drt_summary.json",
    }
    for name, pat in patterns.items():
        for p in sorted(build.glob(pat)):
            doc = read_json(p)
            if isinstance(doc, dict) and "verdict" in doc:
                # keep the directory so per-clade gates stay distinguishable
                key = f"{name}:{p.parent.name}" if "*" in pat.split("/")[0] \
                    or len(pat.split("/")) > 2 else name
                out[key] = doc["verdict"]
    return out


def auspice_tmrca(doc) -> float | None:
    """
    Root date from an Auspice v2 JSON.

    The root's num_date is the TMRCA. Augur writes it under node_attrs; a
    missing value means the tree was not time-scaled, which is a different
    thing from a TMRCA of zero and must not be compared as a number.
    """
    if not isinstance(doc, dict):
        return None
    tree = doc.get("tree")
    if isinstance(tree, list):
        tree = tree[0] if tree else None
    if not isinstance(tree, dict):
        return None
    nd = (tree.get("node_attrs") or {}).get("num_date")
    if isinstance(nd, dict) and "value" in nd:
        try:
            return float(nd["value"])
        except (TypeError, ValueError):
            return None
    return None


# --- the comparisons --------------------------------------------------------

RANK = {"pass": 3, "ok": 3, "weak": 2, "review": 2, "inconclusive": 1,
        "fail": 0, "regression": 0, "unknown": 1}


def compare(prev: Path, curr: Path, pathogen: str,
            max_tmrca_shift: float = 5.0,
            alert_on_loss: bool = True,
            current_auspice=None,
            previous_auspice=None) -> Comparison:
    c = Comparison(pathogen=pathogen, previous=prev.name, current=curr.name)

    # --- accessions ---------------------------------------------------------
    a_prev, a_curr = read_accessions(prev), read_accessions(curr)
    if a_prev is None or a_curr is None:
        c.add("note", "accessions",
              "one of the builds has no accession list, so dataset "
              "composition could not be compared")
    else:
        gained = a_curr - a_prev
        lost = a_prev - a_curr
        c.add("note", "accessions",
              f"{len(a_prev)} -> {len(a_curr)} ({len(a_curr) - len(a_prev):+d}); "
              f"{len(gained)} added, {len(lost)} removed")
        if lost:
            shown = ", ".join(sorted(lost)[:10])
            more = f" and {len(lost) - 10} more" if len(lost) > 10 else ""
            level = "error" if alert_on_loss else "warn"
            c.add(level, "accessions",
                  f"{len(lost)} accession(s) present in {prev.name} are absent "
                  f"from {curr.name}: {shown}{more}. GenBank withdraws records "
                  "rarely; a changed query or a stricter curation rule is the "
                  "usual cause, and either is a regression rather than a "
                  "finding")
        # A stable count hides a composition change completely, and a
        # composition change is exactly what a broken query produces.
        if not lost and not gained and a_prev == a_curr:
            c.add("note", "accessions", "the dataset is unchanged")
        elif lost and gained and len(a_prev) == len(a_curr):
            c.add("warn", "accessions",
                  f"the count is unchanged but {len(lost)} record(s) were "
                  "swapped out. A count-only check would have reported nothing")

    # --- temporal signal ----------------------------------------------------
    t_prev, t_curr = read_temporal(prev), read_temporal(curr)
    if t_prev and t_curr:
        for key, label, tol in (("r_squared", "R^2", 0.5),
                                ("slope", "clock slope", 0.5)):
            pv, cv = t_prev.get(key), t_curr.get(key)
            if isinstance(pv, (int, float)) and isinstance(cv, (int, float)) and pv:
                # Signed, because the direction is the whole point: R^2
                # falling is a regression and R^2 rising is not, and an
                # absolute magnitude printed with a + reads as the opposite
                # of what happened.
                change = (cv - pv) / abs(pv)
                msg = f"{label} {pv:.4g} -> {cv:.4g} ({change:+.0%})"
                c.add("warn" if abs(change) > tol else "note", "temporal", msg)

    # --- every gate verdict -------------------------------------------------
    g_prev, g_curr = gate_verdicts(prev), gate_verdicts(curr)
    for name in sorted(set(g_prev) | set(g_curr)):
        pv, cv = g_prev.get(name), g_curr.get(name)
        if pv is None or cv is None:
            c.add("note", "gates",
                  f"{name}: {'absent' if pv is None else pv} -> "
                  f"{'absent' if cv is None else cv}")
            continue
        if pv == cv:
            continue
        before, after = RANK.get(pv, 1), RANK.get(cv, 1)
        level = "error" if after < before else "note"
        direction = "regressed" if after < before else "improved"
        c.add(level, "gates", f"{name} {direction}: {pv} -> {cv}")

    # --- TMRCA --------------------------------------------------------------
    tm_prev = auspice_tmrca(previous_auspice)
    tm_curr = auspice_tmrca(current_auspice)
    if tm_prev is not None and tm_curr is not None:
        shift = abs(tm_curr - tm_prev)
        msg = f"TMRCA {tm_prev:.1f} -> {tm_curr:.1f} ({tm_curr - tm_prev:+.1f} yr)"
        if shift > max_tmrca_shift:
            c.add("error", "tmrca",
                  msg + f"; exceeds alert_on_tmrca_shift_years "
                        f"({max_tmrca_shift}). A month of new sequences moves a "
                        "TMRCA by months. A shift this size is a changed "
                        "dataset, prior or clock model until shown otherwise")
        else:
            c.add("note", "tmrca", msg)
    elif tm_curr is not None or tm_prev is not None:
        c.add("note", "tmrca",
              "only one build has a time-scaled tree, so the TMRCA could not "
              "be compared")

    return c


def render(c: Comparison) -> str:
    icon = {"ok": "ok", "note": "  ", "warn": "! ", "error": "X "}
    L = [f"# Build comparison — {c.pathogen}", "",
         f"`{c.previous}` → `{c.current}`", "",
         f"**Verdict: {c.verdict.upper()}**", ""]
    if c.verdict == "regression":
        L += ["A change this size is a curation regression until shown "
              "otherwise. Check the query and the curation rules before "
              "treating it as a result.", ""]
    for area in ("accessions", "temporal", "gates", "tmrca"):
        rows = [f for f in c.findings if f.area == area]
        if not rows:
            continue
        L.append(f"## {area}")
        for f in rows:
            L.append(f"- {icon.get(f.level, '  ')} {f.message}")
        L.append("")
    L += ["---", "",
          "Tree topology is deliberately not compared: two trees built from "
          "overlapping data differ mostly by reordering equally-supported "
          "splits, so a topological distance is large, noisy and says nothing "
          "about whether anything is wrong."]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--builds-dir", type=Path, required=True,
                    help="builds/<pathogen>, holding dated build directories.")
    ap.add_argument("--pathogen", help="Pathogen id, for the report header. "
                                       "Defaults to the builds-dir name.")
    ap.add_argument("--current", type=Path,
                    help="Auspice JSON from the current build, when there is "
                         "one. CI stops at the temporal gate and has none.")
    ap.add_argument("--previous", type=Path,
                    help="Auspice JSON from the previous build, for the TMRCA "
                         "comparison. Found in the previous build directory "
                         "when not given.")
    ap.add_argument("--max-tmrca-shift", type=float, default=5.0)
    ap.add_argument("--no-alert-on-sequence-loss", action="store_true")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--out-json", type=Path)
    ap.add_argument("--exit-zero", action="store_true",
                    help="Always exit 0. For looking at a comparison without "
                         "failing a workflow.")
    a = ap.parse_args()

    pathogen = a.pathogen or a.builds_dir.name
    builds = dated_builds(a.builds_dir)

    # A first build has nothing to compare against. That is the normal state of
    # a new pathogen, not a failure, and must not block the first publication.
    if len(builds) < 2:
        text = (f"# Build comparison — {pathogen}\n\n"
                f"Only {len(builds)} build in `{a.builds_dir}`; nothing to "
                f"compare against yet.\n")
        print(text)
        if a.out:
            a.out.parent.mkdir(parents=True, exist_ok=True)
            a.out.write_text(text)
        return 0

    prev, curr = builds[-2], builds[-1]

    cur_doc = read_json(a.current) if a.current else None
    if a.current and cur_doc is None:
        print(f"ERROR: could not read {a.current} as JSON", file=sys.stderr)
        return 1
    prev_doc = read_json(a.previous) if a.previous else None
    if prev_doc is None:
        for p in sorted(prev.glob("auspice/*.json")):
            prev_doc = read_json(p)
            if prev_doc:
                break

    c = compare(prev, curr, pathogen,
                max_tmrca_shift=a.max_tmrca_shift,
                alert_on_loss=not a.no_alert_on_sequence_loss,
                current_auspice=cur_doc, previous_auspice=prev_doc)

    text = render(c)
    print(text)
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(text + "\n")
    if a.out_json:
        a.out_json.parent.mkdir(parents=True, exist_ok=True)
        a.out_json.write_text(json.dumps(c.to_dict(), indent=2))

    if c.verdict == "regression" and not a.exit_zero:
        print("\nexiting non-zero: something regressed between builds.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

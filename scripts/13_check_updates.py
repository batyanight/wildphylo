#!/usr/bin/env python3
"""
13_check_updates.py — should the scheduled rebuild actually run?

    python scripts/13_check_updates.py                 # human-readable
    python scripts/13_check_updates.py --github-output # for CI

Asks GenBank how many records the configured query matches NOW, compares that
against the accession list of the most recent build, and decides per pathogen
whether a rebuild is worth the compute.

Why a counting step exists at all
---------------------------------
A monthly rebuild that downloads, curates, aligns and rebuilds a tree for zero
new sequences produces a new dated build directory, a new provenance record and
a new Auspice JSON that are identical in content to last month's. That is not
harmless: it makes `builds/` a log of the cron schedule rather than a log of
when the dataset changed, and it buries the rebuilds that did change something.

The count is one esearch with retmax=0. It downloads nothing.

The query is imported from 01_fetch_sequences rather than reimplemented, so the
count can never be of a different query than the download. A second copy would
drift, and the symptom of drift here is a delta that looks like a real change in
GenBank.

Decisions
---------
delta < 0            LOST records. GenBank withdraws records very rarely, so
                     this is almost always a changed query or a truncated
                     previous download. Reported as an error and the rebuild
                     does NOT proceed: rebuilding would quietly replace a
                     larger dataset with a smaller one.
delta = 0            nothing to do.
0 < delta < min      below nextstrain.update.min_new_accessions: not worth a
                     rebuild, but recorded so the count is visible.
min <= delta <= max  proceed.
delta > max          proceed, flagged. A bulk submission is exactly when a
                     rebuild is most worth doing and most worth looking at.
no previous build    proceed. There is nothing to compare against, and the
                     first build has to happen somehow.

Environment
-----------
NCBI_EMAIL, NCBI_API_KEY   credentials, never read from a config file
FORCE=true                 proceed regardless of the delta
ONLY=<id>                  restrict the check to one pathogen
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lib.config import validate  # noqa: E402


def _fetch_module():
    """
    Import 01_fetch_sequences. Its name starts with a digit, so it is not a
    legal module name and cannot be imported normally. Loading it by path is
    deliberate: sharing build_query is the only way to guarantee the count and
    the download describe the same set of records.
    """
    spec = importlib.util.spec_from_file_location(
        "fetch01", HERE / "01_fetch_sequences.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def count_records(query: str) -> int:
    """esearch with retmax=0: returns the match count, downloads nothing."""
    from Bio import Entrez
    email = os.environ.get("NCBI_EMAIL")
    if not email:
        raise RuntimeError(
            "NCBI_EMAIL is not set. Entrez requires it, and it belongs in the "
            "environment rather than a config file, which gets committed.")
    Entrez.email = email
    if os.environ.get("NCBI_API_KEY"):
        Entrez.api_key = os.environ["NCBI_API_KEY"]
    handle = Entrez.esearch(db="nuccore", term=query, retmax=0)
    return int(Entrez.read(handle)["Count"])


def previous_build(builds_dir: Path, pathogen: str) -> tuple[str | None, int | None]:
    """
    The most recent dated build and its accession count.

    Reads the .acc file rather than the metadata, because .acc is what the
    query returned -- curation drops records afterwards, and comparing a
    post-curation count against a pre-curation one would report a change every
    time a curation rule was edited.
    """
    root = builds_dir / pathogen
    if not root.is_dir():
        return None, None
    dated = sorted((d for d in root.iterdir() if d.is_dir() and d.name.isdigit()),
                   key=lambda d: d.name)
    for d in reversed(dated):
        acc = d / "raw" / f"{pathogen}.acc"
        if acc.is_file():
            n = sum(1 for line in acc.read_text().splitlines() if line.strip())
            return d.name, n
    return (dated[-1].name if dated else None), None


@dataclass
class Decision:
    pathogen: str
    proceed: bool
    reason: str
    current: int | None = None
    previous: int | None = None
    previous_build: str | None = None
    delta: int | None = None
    level: str = "note"          # note | warn | error

    def line(self) -> str:
        mark = {"note": " ", "warn": "!", "error": "X"}[self.level]
        n = "-" if self.current is None else str(self.current)
        p = "-" if self.previous is None else str(self.previous)
        d = "" if self.delta is None else f"{self.delta:+d}"
        return (f" {mark} {self.pathogen:<22}{p:>8}{n:>8}{d:>8}   "
                f"{'REBUILD' if self.proceed else 'skip':<8} {self.reason}")


def decide(cfg: dict, pathogen: str, current: int,
           previous: int | None, prev_build: str | None,
           force: bool = False) -> Decision:
    upd = (cfg.get("nextstrain") or {}).get("update") or {}
    lo = upd.get("min_new_accessions", 1)
    hi = upd.get("max_new_accessions", 10 ** 9)

    if previous is None:
        return Decision(pathogen, True, "no previous build to compare against",
                        current, None, prev_build, None, "note")

    delta = current - previous
    if delta < 0:
        # Forcing past this is allowed, but it has to be a deliberate act: the
        # default must not be to overwrite a larger dataset with a smaller one.
        return Decision(
            pathogen, bool(force),
            f"dataset LOST {abs(delta)} records; check the query before "
            f"rebuilding" + (" (FORCE set)" if force else ""),
            current, previous, prev_build, delta, "error")
    if force:
        return Decision(pathogen, True, "FORCE is set",
                        current, previous, prev_build, delta, "note")
    if delta == 0:
        return Decision(pathogen, False, "no new records",
                        current, previous, prev_build, delta, "note")
    if delta < lo:
        return Decision(pathogen, False,
                        f"{delta} new, below min_new_accessions ({lo})",
                        current, previous, prev_build, delta, "note")
    if delta > hi:
        return Decision(pathogen, True,
                        f"{delta} new exceeds max_new_accessions ({hi}); "
                        f"a bulk submission or a changed query deserves eyes",
                        current, previous, prev_build, delta, "warn")
    return Decision(pathogen, True, f"{delta} new records",
                    current, previous, prev_build, delta, "note")


def load_configs(config_dir: Path, only: str | None) -> list[tuple[Path, dict]]:
    out = []
    for p in sorted(config_dir.glob("*.yaml")):
        if p.name.startswith("_"):
            continue
        cfg = yaml.safe_load(p.read_text()) or {}
        cfg["_path"] = str(p)
        if only and cfg.get("id") != only:
            continue
        errors, _ = validate(cfg)
        if errors:
            # Refuse rather than skip. A config that stopped validating is the
            # reason to look, and skipping it makes the rebuild silently narrow.
            raise SystemExit(
                f"{p} is not usable, so the scheduled rebuild cannot decide "
                f"whether to run it:\n" + "\n".join(f"  - {e}" for e in errors))
        out.append((p, cfg))
    if only and not out:
        raise SystemExit(f"no config has id {only!r}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config-dir", type=Path, default=Path("config/pathogen"))
    ap.add_argument("--builds-dir", type=Path, default=Path("builds"))
    ap.add_argument("--only", help="Restrict to one pathogen id. Overrides $ONLY.")
    ap.add_argument("--force", action="store_true",
                    help="Proceed regardless of the delta. Overrides $FORCE.")
    ap.add_argument("--github-output", action="store_true",
                    help="Append proceed= and pathogens= to $GITHUB_OUTPUT.")
    ap.add_argument("--out-json", type=Path)
    a = ap.parse_args()

    only = a.only or (os.environ.get("ONLY") or "").strip() or None
    force = a.force or os.environ.get("FORCE", "").lower() in ("1", "true", "yes")

    configs = load_configs(a.config_dir, only)
    fetch = _fetch_module()

    decisions: list[Decision] = []
    for path, cfg in configs:
        pid = cfg["id"]
        query = fetch.build_query(cfg)
        prev_build, previous = previous_build(a.builds_dir, pid)
        try:
            current = count_records(query)
        except Exception as e:                            # noqa: BLE001
            # A network or credential failure must not look like "no new
            # records". Skipping on an error would make an outage indefinitely
            # silent, and the rebuild would simply stop happening.
            decisions.append(Decision(
                pid, False, f"could not reach Entrez: {e}",
                None, previous, prev_build, None, "error"))
            continue
        decisions.append(decide(cfg, pid, current, previous, prev_build, force))

    W = 86
    print("=" * W)
    print("SCHEDULED REBUILD CHECK")
    print("=" * W)
    print(f"   {'pathogen':<22}{'previous':>8}{'current':>8}{'delta':>8}   "
          f"{'action':<8} reason")
    print("-" * W)
    for d in decisions:
        print(d.line())
    print("=" * W)

    proceed = [d.pathogen for d in decisions if d.proceed]
    errors = [d for d in decisions if d.level == "error"]
    if proceed:
        print(f"rebuilding: {', '.join(proceed)}")
    else:
        print("nothing to rebuild")

    if a.out_json:
        a.out_json.parent.mkdir(parents=True, exist_ok=True)
        a.out_json.write_text(json.dumps(
            {"proceed": bool(proceed),
             "pathogens": proceed,
             "decisions": [asdict(d) for d in decisions]}, indent=2))
        print(f"wrote {a.out_json}")

    if a.github_output:
        target = os.environ.get("GITHUB_OUTPUT")
        if not target:
            print("ERROR: --github-output given but $GITHUB_OUTPUT is unset",
                  file=sys.stderr)
            return 1
        with open(target, "a") as fh:
            fh.write(f"proceed={'true' if proceed else 'false'}\n")
            fh.write(f"pathogens={json.dumps(proceed)}\n")

    # An Entrez failure or a shrinking dataset is reported but does not fail
    # the step: the matrix job is gated on `proceed` instead. Failing here
    # would mean one unreachable pathogen cancels the rebuild of the others.
    for d in errors:
        print(f"note: {d.pathogen}: {d.reason}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

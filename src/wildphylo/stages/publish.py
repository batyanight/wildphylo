"""publish stage: decide whether to rebuild, and hand a finished dataset to its pathogen repo.

Steps: ``check-updates`` (should the scheduled rebuild run?) and ``dataset``
(copy a finished dataset into the pathogen repo).
"""

from __future__ import annotations

import json
import logging
import os
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml

from wildphylo.cli import step
from wildphylo.config import validate
from wildphylo.layout import dataset_path
from wildphylo.stages.fetch import build_query
from wildphylo.stepargs import Build, UsageError, add_build_dir_arguments, fail, resolve

if TYPE_CHECKING:
    import argparse
    from collections.abc import Callable


log = logging.getLogger(__name__)


# ===========================================================================
# publish check-updates
# ===========================================================================


def count_records(query: str) -> int:
    """Count matches with esearch retmax=0: returns the match count, downloads nothing."""
    from Bio import Entrez

    email = os.environ.get("NCBI_EMAIL")
    if not email:
        msg = (
            "NCBI_EMAIL is not set. Entrez requires it, and it belongs in the "
            "environment rather than a config file, which gets committed."
        )
        raise RuntimeError(msg)
    Entrez.email = email  # type: ignore[assignment]  # Bio.Entrez types the module global as None
    if os.environ.get("NCBI_API_KEY"):
        Entrez.api_key = os.environ["NCBI_API_KEY"]  # type: ignore[assignment]  # as above
    handle = Entrez.esearch(  # type: ignore[no-untyped-call]  # Bio.Entrez ships no annotations
        db="nuccore", term=query, retmax=0
    )
    record = Entrez.read(handle)  # type: ignore[no-untyped-call]  # Bio.Entrez ships no annotations
    return int(record["Count"])


def previous_build(builds_dir: Path, pathogen: str) -> tuple[str | None, int | None]:
    """The most recent dated build and its accession count.

    Reads the .acc file rather than the metadata, because .acc is what the
    query returned -- curation drops records afterwards, and comparing a
    post-curation count against a pre-curation one would report a change every
    time a curation rule was edited.
    """
    root = builds_dir / pathogen
    if not root.is_dir():
        return None, None
    dated = sorted(
        (d for d in root.iterdir() if d.is_dir() and d.name.isdigit()), key=lambda d: d.name
    )
    for d in reversed(dated):
        acc = d / "raw" / f"{pathogen}.acc"
        if acc.is_file():
            n = sum(1 for line in acc.read_text().splitlines() if line.strip())
            return d.name, n
    return (dated[-1].name if dated else None), None


@dataclass
class Decision:
    """Whether one pathogen's scheduled rebuild proceeds, and why."""

    pathogen: str
    proceed: bool
    reason: str
    current: int | None = None
    previous: int | None = None
    previous_build: str | None = None
    delta: int | None = None
    level: str = "note"  # note | warn | error

    def line(self) -> str:
        """One row of the report table."""
        mark = {"note": " ", "warn": "!", "error": "X"}[self.level]
        n = "-" if self.current is None else str(self.current)
        p = "-" if self.previous is None else str(self.previous)
        d = "" if self.delta is None else f"{self.delta:+d}"
        return (
            f" {mark} {self.pathogen:<22}{p:>8}{n:>8}{d:>8}   "
            f"{'REBUILD' if self.proceed else 'skip':<8} {self.reason}"
        )


def decide(  # ordered rule list, one exit per rule
    cfg: dict[str, Any],
    pathogen: str,
    current: int,
    previous: int | None,
    prev_build: str | None,
    force: bool = False,
) -> Decision:
    """Apply the decision table to one pathogen's record counts."""
    upd = (cfg.get("nextstrain") or {}).get("update") or {}
    lo = upd.get("min_new_accessions", 1)
    hi = upd.get("max_new_accessions", 10**9)

    if previous is None:
        return Decision(
            pathogen,
            True,
            "no previous build to compare against",
            current,
            None,
            prev_build,
            None,
            "note",
        )

    delta = current - previous
    if delta < 0:
        # Forcing past this is allowed, but it has to be a deliberate act: the
        # default must not be to overwrite a larger dataset with a smaller one.
        return Decision(
            pathogen,
            bool(force),
            f"dataset LOST {abs(delta)} records; check the query before "
            f"rebuilding" + (" (FORCE set)" if force else ""),
            current,
            previous,
            prev_build,
            delta,
            "error",
        )
    if force:
        return Decision(
            pathogen, True, "FORCE is set", current, previous, prev_build, delta, "note"
        )
    if delta == 0:
        return Decision(
            pathogen, False, "no new records", current, previous, prev_build, delta, "note"
        )
    if delta < lo:
        return Decision(
            pathogen,
            False,
            f"{delta} new, below min_new_accessions ({lo})",
            current,
            previous,
            prev_build,
            delta,
            "note",
        )
    if delta > hi:
        return Decision(
            pathogen,
            True,
            f"{delta} new exceeds max_new_accessions ({hi}); "
            f"a bulk submission or a changed query deserves eyes",
            current,
            previous,
            prev_build,
            delta,
            "warn",
        )
    return Decision(
        pathogen, True, f"{delta} new records", current, previous, prev_build, delta, "note"
    )


def load_configs(config_dir: Path, only: str | None) -> list[tuple[Path, dict[str, Any]]]:
    """Validated pathogen configs in ``config_dir``, optionally only the one with id ``only``."""
    out: list[tuple[Path, dict[str, Any]]] = []
    for p in sorted(config_dir.glob("*.yaml")):
        if p.name.startswith("_"):
            continue
        cfg: dict[str, Any] = yaml.safe_load(p.read_text()) or {}
        cfg["_path"] = str(p)
        if only and cfg.get("id") != only:
            continue
        if not only and (cfg.get("fetch") or {}).get("max_records"):
            # A capped fetch is a smoke dataset, not an analysis (see
            # fetch.max_records); the scheduled rebuild never runs it.
            continue
        errors, _ = validate(cfg)
        if errors:
            # Refuse rather than skip. A config that stopped validating is the
            # reason to look, and skipping it makes the rebuild silently narrow.
            raise SystemExit(
                f"{p} is not usable, so the scheduled rebuild cannot decide "
                f"whether to run it:\n" + "\n".join(f"  - {e}" for e in errors)
            )
        out.append((p, cfg))
    if only and not out:
        raise SystemExit(f"no config has id {only!r}")
    return out


def add_arguments_check_updates(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo publish check-updates``."""
    parser.add_argument("--config-dir", type=Path, default=Path("config/pathogen"))
    parser.add_argument("--builds-dir", type=Path, default=Path("builds"))
    parser.add_argument("--only", help="Restrict to one pathogen id. Overrides $ONLY.")
    parser.add_argument(
        "--force", action="store_true", help="Proceed regardless of the delta. Overrides $FORCE."
    )
    parser.add_argument(
        "--github-output",
        action="store_true",
        help="Append proceed= and pathogens= to $GITHUB_OUTPUT.",
    )
    parser.add_argument("--out-json", type=Path)


@step("publish", "check-updates", "Should the scheduled rebuild run? Counts new GenBank records")
def check_updates(args: argparse.Namespace) -> int:
    """Should the scheduled rebuild actually run? Count new GenBank records per pathogen.

    Asks GenBank how many records the configured query matches NOW, compares that
    against the accession list of the most recent build, and decides per pathogen
    whether a rebuild is worth the compute.

    Why a counting step exists at all: a monthly rebuild that downloads, curates,
    aligns and rebuilds a tree for zero new sequences produces a new dated build
    directory, a new provenance record and a new Auspice JSON that are identical in
    content to last month's. That makes ``builds/`` a log of the cron schedule
    rather than a log of when the dataset changed, and it buries the rebuilds that
    did change something. The count is one esearch with retmax=0. It downloads
    nothing.

    The query comes from :func:`wildphylo.stages.fetch.build_query` rather than
    being reimplemented, so the count can never be of a different query than the
    download. A second copy would drift, and the symptom of drift here is a delta
    that looks like a real change in GenBank.

    Decisions::

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

    Environment: ``NCBI_EMAIL`` and ``NCBI_API_KEY`` are credentials, never read
    from a config file. ``FORCE=true`` proceeds regardless of the delta.
    ``ONLY=<id>`` restricts the check to one pathogen.
    """
    only = args.only or (os.environ.get("ONLY") or "").strip() or None
    force = args.force or os.environ.get("FORCE", "").lower() in ("1", "true", "yes")

    configs = load_configs(args.config_dir, only)

    decisions: list[Decision] = []
    for _path, cfg in configs:
        pid = cfg["id"]
        query = build_query(cfg)
        prev_build, previous = previous_build(args.builds_dir, pid)
        try:
            current = count_records(query)
        except Exception as e:  # noqa: BLE001
            # A network or credential failure must not look like "no new
            # records". Skipping on an error would make an outage indefinitely
            # silent, and the rebuild would simply stop happening.
            decisions.append(
                Decision(
                    pid,
                    False,
                    f"could not reach Entrez: {e}",
                    None,
                    previous,
                    prev_build,
                    None,
                    "error",
                )
            )
            continue
        decisions.append(decide(cfg, pid, current, previous, prev_build, force))

    width = 86
    proceed = [d.pathogen for d in decisions if d.proceed]
    errors = [d for d in decisions if d.level == "error"]
    report = [
        "=" * width,
        "SCHEDULED REBUILD CHECK",
        "=" * width,
        f"   {'pathogen':<22}{'previous':>8}{'current':>8}{'delta':>8}   {'action':<8} reason",
        "-" * width,
        *(d.line() for d in decisions),
        "=" * width,
        f"rebuilding: {', '.join(proceed)}" if proceed else "nothing to rebuild",
    ]
    log.info("\n%s", "\n".join(report))

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(
            json.dumps(
                {
                    "proceed": bool(proceed),
                    "pathogens": proceed,
                    "decisions": [asdict(d) for d in decisions],
                },
                indent=2,
            )
        )
        log.info("wrote %s", args.out_json)

    if args.github_output:
        target = os.environ.get("GITHUB_OUTPUT")
        if not target:
            log.error("--github-output given but $GITHUB_OUTPUT is unset")
            return 1
        with Path(target).open("a") as fh:
            fh.write(f"proceed={'true' if proceed else 'false'}\n")
            fh.write(f"pathogens={json.dumps(proceed)}\n")

    # An Entrez failure or a shrinking dataset is reported but does not fail
    # the step: the matrix job is gated on `proceed` instead. Failing here
    # would mean one unreachable pathogen cancels the rebuild of the others.
    for d in errors:
        log.warning("%s: %s", d.pathogen, d.reason)
    return 0


# ===========================================================================
# publish dataset
# ===========================================================================


def _read_json_object(path: Path) -> dict[str, Any] | None:
    """The JSON object in ``path``, or None if it is missing, unparseable or not an object."""
    try:
        doc = json.loads(path.read_text())
    except (OSError, ValueError) as e:
        log.error("could not read %s as JSON: %s", path, e)  # noqa: TRY400 -- a message, not a traceback
        return None
    if not isinstance(doc, dict):
        log.error("%s is not a JSON object", path)
        return None
    return doc


_DATASET_DERIVED: dict[str, Callable[[Build], Any]] = {
    "dataset": lambda b: dataset_path(b.layout, b.cfg, b.need("locus")),
    "comparison": lambda b: b.layout.build_comparison_json(),
    "provenance": lambda b: b.layout.provenance(),
}


def add_arguments_dataset(parser: argparse.ArgumentParser) -> None:
    """Flags of ``wildphylo publish dataset``."""
    parser.add_argument(
        "--dataset",
        type=Path,
        help="Auspice JSON produced by `wildphylo auspice dataset`. Required without --build-dir.",
    )
    parser.add_argument(
        "--comparison",
        type=Path,
        help="build_comparison.json written by `wildphylo gates compare`. "
        "Required without --build-dir.",
    )
    parser.add_argument("--config", type=Path, help="Pathogen config. Required with --build-dir.")
    add_build_dir_arguments(parser)
    parser.add_argument(
        "--pathogen-repo",
        type=Path,
        required=True,
        help="Local checkout of the per-pathogen repo.",
    )
    parser.add_argument(
        "--name",
        required=True,
        help="Dataset name, e.g. cdv-phylodynamics. Written to "
        "<pathogen-repo>/auspice/<name>.json.",
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        help="provenance.json of the build; copied to "
        "<pathogen-repo>/builds/<build_date>/provenance.json. Required without --build-dir.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Publish even when the comparison verdict is 'regression'.",
    )


@step("publish", "dataset", "Copy a finished dataset and its provenance into the pathogen repo")
def dataset(args: argparse.Namespace) -> int:
    """Copy a finished Auspice dataset and its provenance into the pathogen repo.

    The copy goes into the local checkout of the per-pathogen repo, unless the build
    comparison says the build regressed; ``--force`` overrides a regression. It
    never runs git: committing and pushing the pathogen repo is a reviewed,
    deliberate act.

    With ``--build-dir`` (plus ``--config`` and ``--locus``) the dataset, comparison
    and provenance default to the build layout's files.
    """
    try:
        resolve(
            args,
            derived=_DATASET_DERIVED,
            required=tuple(_DATASET_DERIVED),
            validate_config=False,
        )
    except UsageError as e:
        return fail(e)
    for label, path in (
        ("dataset", args.dataset),
        ("comparison", args.comparison),
        ("provenance", args.provenance),
    ):
        if not path.is_file():
            log.error("%s file not found: %s", label, path)
            return 1
    if not args.pathogen_repo.is_dir():
        log.error("pathogen repo is not a directory: %s", args.pathogen_repo)
        return 1

    comparison = _read_json_object(args.comparison)
    provenance = _read_json_object(args.provenance)
    if comparison is None or provenance is None:
        return 1

    # A comparison with no verdict is not evidence of a clean build. Refuse
    # rather than treat "unreadable" as "fine".
    verdict = comparison.get("verdict")
    if not isinstance(verdict, str):
        log.error("%s has no 'verdict'; is it a `gates compare` output?", args.comparison)
        return 1
    if verdict == "regression":
        if not args.force:
            log.error(
                "refusing to publish %s: the build comparison reports a regression "
                "(%s). Review it, or pass --force to publish anyway.",
                args.name,
                args.comparison,
            )
            return 1
        log.warning("publishing %s despite a regression verdict (--force)", args.name)

    build_date = provenance.get("build_date")
    if not isinstance(build_date, str) or not build_date.isdigit():
        log.error("%s has no usable 'build_date' (expected YYYYMMDD)", args.provenance)
        return 1

    dataset_out = args.pathogen_repo / "auspice" / f"{args.name}.json"
    provenance_out = args.pathogen_repo / "builds" / build_date / "provenance.json"
    for src, dst in ((args.dataset, dataset_out), (args.provenance, provenance_out)):
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        log.info("wrote %s", dst)
    log.info(
        "%s is copied, not committed: review and commit it in %s yourself.",
        args.name,
        args.pathogen_repo,
    )
    return 0

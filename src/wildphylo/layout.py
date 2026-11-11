"""The build layout: every file a build writes, in one place.

A build lives in ``builds/<pathogen>/<YYYYMMDD>/`` (see GLOSSARY.md, "build").
This module is the only place that knows where inside it each file goes, so
the Snakefile, the steps and the tests all ask here rather than assembling
``f"{out}/tree/{locus}_ml.treefile"`` strings of their own. Moving a file means
changing one method.

Every method returns a :class:`~pathlib.Path` and does no I/O: nothing is
created, checked or read. ``locus``, ``seed`` and ``rep`` are plain strings, so
a Snakemake wildcard pattern such as ``"{locus}"`` passes straight through and
comes out as part of the path::

    >>> L = BuildLayout(Path("builds/cdv/20260101"), "cdv", None)
    >>> L.alignment("{locus}")
    PosixPath('builds/cdv/20260101/aligned/{locus}_aligned.fasta')

The relative paths here are the ones existing builds already use. Changing one
orphans every earlier build of that pathogen, so a rename is a decision, not a
tidy-up.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class BuildLayout:
    """Every file a build writes, by stage, under ``builds/<pathogen>/<build date>/``."""

    root: Path
    """The build directory: ``builds/<id>/<YYYYMMDD>``."""
    pathogen: str
    """The pathogen config's ``id``; part of several file names."""
    trait: str | None
    """The discrete trait reconstructed by BEAST (e.g. ``host_group``), or None
    when the config has no discrete trait analysis."""

    @classmethod
    def from_config(
        cls, cfg: dict[str, Any], build_date: str, builds_dir: Path = Path("builds")
    ) -> BuildLayout:
        """Layout of the build of ``cfg`` stamped ``build_date``.

        ``trait`` is set only when ``beast.discrete_trait.enabled`` is true;
        its name defaults to ``host_group`` as in the config template.
        """
        pathogen = str(cfg["id"])
        dta = (cfg.get("beast") or {}).get("discrete_trait") or {}
        trait = str(dta.get("trait", "host_group")) if dta.get("enabled") else None
        return cls(Path(builds_dir) / pathogen / str(build_date), pathogen, trait)

    @classmethod
    def at(cls, root: Path | str, cfg: dict[str, Any]) -> BuildLayout:
        """Layout of an existing build directory (``--build-dir``) for ``cfg``.

        The same rules as :meth:`from_config` for ``pathogen`` and ``trait``;
        the root is taken as given, so a build can live anywhere.
        """
        return cls(Path(root), str(cfg["id"]), cls.from_config(cfg, "0").trait)

    # --- the build as a whole -------------------------------------------------

    def pathogen_dir(self) -> Path:
        """``builds/<id>``: every dated build of this pathogen (the comparison gate scans it)."""
        return self.root.parent

    def provenance(self) -> Path:
        return self.root / "provenance.json"

    def log(self, name: str) -> Path:
        """``logs/<name>.log``; the numbered names (``01_fetch``) are just names."""
        return self.root / "logs" / f"{name}.log"

    # --- fetch ----------------------------------------------------------------

    def raw_records(self) -> Path:
        return self.root / "raw" / f"{self.pathogen}.gb"

    def accession_snapshot(self) -> Path:
        return self.root / "raw" / f"{self.pathogen}.acc"

    # --- curate ---------------------------------------------------------------

    def interim_dir(self) -> Path:
        """Where curation writes everything below."""
        return self.root / "interim"

    def metadata(self) -> Path:
        return self.interim_dir() / "metadata_clean.tsv"

    def needs_review(self) -> Path:
        return self.interim_dir() / "needs_review.tsv"

    def review_marker(self) -> Path:
        return self.interim_dir() / ".review_ok"

    def flagged(self) -> Path:
        return self.interim_dir() / "flagged.tsv"

    def exclusions(self) -> Path:
        return self.interim_dir() / "exclusions.tsv"

    def locus_sequences(self, locus: str) -> Path:
        return self.interim_dir() / f"{locus}.fasta"

    def anchors(self, locus: str) -> Path:
        return self.interim_dir() / f"{locus}_anchors.fasta"

    # --- align ----------------------------------------------------------------

    def alignment(self, locus: str) -> Path:
        return self.root / "aligned" / f"{locus}_aligned.fasta"

    def alignment_qc(self, locus: str) -> Path:
        return self.root / "qc" / f"{locus}_alignment_qc.tsv"

    def alignment_qc_plot(self, locus: str) -> Path:
        return self.root / "qc" / f"{locus}_coverage.png"

    # --- tree -----------------------------------------------------------------

    def ml_tree(self, locus: str) -> Path:
        return self.root / "tree" / f"{locus}_ml.treefile"

    def ml_tree_prefix(self, locus: str) -> Path:
        """IQ-TREE's ``--prefix``: every IQ-TREE output starts with this."""
        return self.root / "tree" / f"{locus}_ml"

    # --- gates: temporal signal, congruence -----------------------------------

    def temporal_report(self, locus: str) -> Path:
        return self.root / "temporal" / f"{locus}_temporal_signal.json"

    def temporal_report_md(self, locus: str) -> Path:
        return self.root / "temporal" / f"{locus}_temporal_signal.md"

    def temporal_plot(self, locus: str) -> Path:
        return self.root / "temporal" / f"{locus}_root_to_tip.png"

    def congruence_report(self) -> Path:
        return self.root / "segments" / "congruence.json"

    def congruence_md(self) -> Path:
        return self.root / "segments" / "congruence.md"

    def reassortant_candidates(self) -> Path:
        return self.root / "segments" / "reassortant_candidates.tsv"

    # --- beast: subsample, XML, chains ----------------------------------------

    def subsample(self, locus: str) -> Path:
        return self.root / "subsampled" / f"{locus}_sub.fasta"

    def subsample_metadata(self, locus: str) -> Path:
        return self.root / "subsampled" / f"{locus}_sub_metadata.tsv"

    def beast_locus_dir(self, locus: str) -> Path:
        return self.root / "beast" / locus

    def beast_chain_dir(self, locus: str, seed: str) -> Path:
        """One chain's directory; also BEAST's ``-prefix`` (with a trailing slash)."""
        return self.beast_locus_dir(locus) / f"seed{seed}"

    def beast_xml(self, locus: str, seed: str) -> Path:
        return self.beast_chain_dir(locus, seed) / f"{self.pathogen}.xml"

    def beast_traits(self, locus: str, seed: str) -> Path:
        return self.beast_chain_dir(locus, seed) / "traits.txt"

    def beast_log(self, locus: str, seed: str) -> Path:
        return self.beast_chain_dir(locus, seed) / f"{self.pathogen}.log"

    def beast_trees(self, locus: str, seed: str) -> Path:
        return self.beast_chain_dir(locus, seed) / f"{self.pathogen}.trees"

    def beast_trait_trees(self, locus: str, seed: str) -> Path:
        """The tree file carrying the reconstructed trait states; needs a trait."""
        if self.trait is None:
            raise ValueError("no discrete trait in this build: beast_trait_trees is undefined")
        return self.beast_chain_dir(locus, seed) / f"{self.pathogen}_{self.trait}.trees"

    # --- beast: summaries and gates -------------------------------------------

    def convergence_report(self, locus: str) -> Path:
        return self.beast_locus_dir(locus) / "convergence.json"

    def combined_log(self, locus: str) -> Path:
        return self.beast_locus_dir(locus) / "combined.log"

    def combined_trees(self, locus: str) -> Path:
        """Chains' trees after burn-in, merged; the input to the MCC tree."""
        return self.beast_locus_dir(locus) / "combined.trees"

    def mcc_tree(self, locus: str) -> Path:
        return self.beast_locus_dir(locus) / "mcc.tree"

    def trait_signal_report(self, locus: str) -> Path:
        return self.beast_locus_dir(locus) / "trait_signal.json"

    def trait_signal_md(self, locus: str) -> Path:
        return self.beast_locus_dir(locus) / "trait_signal.md"

    # --- date randomisation test ----------------------------------------------

    def drt_dir(self, locus: str) -> Path:
        return self.root / "drt" / locus

    def drt_xml(self, locus: str, rep: str) -> Path:
        return self.drt_dir(locus) / f"rep{rep}" / f"{self.pathogen}_drt.xml"

    def drt_log(self, locus: str, rep: str) -> Path:
        return self.drt_dir(locus) / f"rep{rep}" / f"{self.pathogen}_drt.log"

    def drt_summary(self, locus: str) -> Path:
        return self.drt_dir(locus) / "drt_summary.json"

    def drt_plot(self, locus: str) -> Path:
        return self.drt_dir(locus) / "drt.png"

    # --- auspice and comparison -----------------------------------------------

    def auspice_dataset(self, build_name: str, locus: str | None = None) -> Path:
        """The Auspice dataset JSON; ``<build_name>_<locus>.json`` when per-locus.

        The dataset name never contains a date (GLOSSARY.md, "dataset name").
        """
        stem = build_name if locus is None else f"{build_name}_{locus}"
        return self.root / "auspice" / f"{stem}.json"

    def build_comparison(self) -> Path:
        return self.root / "build_comparison.md"

    def build_comparison_json(self) -> Path:
        return self.root / "build_comparison.json"


def dataset_path(layout: BuildLayout, cfg: dict[str, Any], locus: str | None) -> Path:
    """Where ``auspice dataset`` writes for ``locus``, by the Snakefile's rule.

    A multi-locus config has one dataset per locus, named ``<build_name>_<locus>``;
    a single-locus config keeps the bare ``<build_name>`` so a published URL does
    not move (docs/NEXTSTRAIN.md). ``locus`` may be None only for a single-locus
    config. Raises ValueError when the config cannot name the dataset.
    """
    build_name = (cfg.get("nextstrain") or {}).get("build_name")
    if not build_name:
        raise ValueError("the config has no nextstrain.build_name")
    n_loci = sum(1 for item in cfg["loci"] if item.get("analyse")) if "loci" in cfg else 1
    if n_loci <= 1:
        return layout.auspice_dataset(str(build_name))
    if locus is None:
        raise ValueError("a multi-locus config has one dataset per locus; the locus is needed")
    return layout.auspice_dataset(str(build_name), locus)

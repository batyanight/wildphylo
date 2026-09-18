"""
hosts.py — host normalisation, taxon-agnostic.

The original matcher tested `pattern in raw.lower()`, an unbounded substring
test. That makes "hot dog vendor" a domestic dog and, more plausibly, makes
"prairie dog" one too. The CDV table only escapes this because someone
hand-ordered `Cynomys` above `dog`. A new pathogen's fresh host table has no
such protection, and the failure is silent: a mislabelled host becomes a
mislabelled tip becomes a spurious host-transition rate.

Two changes:

  1. Patterns match on word boundaries by default. `dog` matches "dog" and
     "wild dog" but not "hot dog vendor"... it *does* still match "prairie dog",
     because that is a genuine two-word phrase containing the word. Which is
     why:
  2. `audit_host_table` reports shadowing at load time — any pattern that can
     never win because an earlier pattern subsumes it, and any pattern that is
     a proper substring of another. You see the ambiguity before the run, not
     after the tree.

A pattern wrapped in slashes (/.../)  is treated as a regex, for the cases
where word boundaries are not enough.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


# Specimen descriptions that arrive via /isolation_source when /host is absent.
SAMPLE_TYPE_RE = re.compile(
    r"(?:[\w\s-]*\b(?:urine|blood|serum|swab|tissue|lung|brain|spleen|"
    r"faeces|feces|stool|saliva|csf|plasma|biopsy|necropsy|cell culture|"
    r"supernatant)\b[\w\s-]*)")


@dataclass(frozen=True)
class HostRule:
    pattern: str
    canonical: str
    group: str
    regex: re.Pattern
    is_regex: bool
    lineno: int


@dataclass(frozen=True)
class HostMatch:
    canonical: str
    group: str
    ambiguous: bool
    matched_pattern: str = ""
    reason: str = ""


def normalize_text(raw: str) -> str:
    """
    Tidy a /host string before matching.

    Submitters write `Canis_lupus_familiaris` with underscores, and underscore
    is a word character — so a word-boundary pattern can never match inside it.
    Also collapses whitespace and strips the surrounding punctuation GenBank
    accumulates.
    """
    s = str(raw).replace("_", " ")
    s = re.sub(r"[;,]", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def _compile(pattern: str) -> tuple[re.Pattern, bool]:
    if len(pattern) > 1 and pattern.startswith("/") and pattern.endswith("/"):
        return re.compile(pattern[1:-1], re.IGNORECASE), True
    # Word-boundary match. \b is wrong at non-word edges (e.g. a pattern ending
    # in "."), so guard with lookarounds that tolerate those.
    #
    # An optional plural suffix is allowed because submitters write "dogs" and
    # "Arctic foxes". Pure word-boundary matching rejected those, which was a
    # REGRESSION against the old substring matcher -- it cost ~100 domestic dog
    # records on the CDV dataset before this was added. The suffix is only
    # tried when the pattern does not already end in "s", so "vulpes" does not
    # acquire a spurious alternative.
    esc = re.escape(pattern)
    suffix = "" if pattern.endswith("s") else "(?:s|es)?"
    return re.compile(rf"(?<!\w){esc}{suffix}(?!\w)", re.IGNORECASE), False


def load_host_table(path: Path) -> list[HostRule]:
    """
    Load a TSV of  pattern <TAB> canonical_host <TAB> host_group.
    Blank lines and lines starting with # are ignored. Order matters: the
    first matching rule wins.
    """
    rules: list[HostRule] = []
    for i, line in enumerate(Path(path).read_text().splitlines(), start=1):
        line = line.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 3:
            raise ValueError(
                f"{path}:{i}: expected 3 tab-separated fields "
                f"(pattern, canonical, group), got {len(parts)}: {line!r}"
            )
        pat, canon, group = parts[0], parts[1], parts[2]
        rx, is_rx = _compile(pat.lower())
        rules.append(HostRule(pat.lower(), canon, group, rx, is_rx, i))
    if not rules:
        raise ValueError(f"{path}: no host rules found")
    return rules


def audit_host_table(rules: list[HostRule]) -> list[str]:
    """
    Return human-readable warnings about a host table. Called at load time by
    the curation step; failing to act on these is a choice, but an informed one.
    """
    warnings: list[str] = []

    seen: dict[str, HostRule] = {}
    for r in rules:
        if r.pattern in seen:
            warnings.append(
                f"duplicate pattern {r.pattern!r} at lines "
                f"{seen[r.pattern].lineno} and {r.lineno}; the later one is dead"
            )
        else:
            seen[r.pattern] = r

    # Shadowing: an earlier pattern that matches a later pattern's own text
    # means the later rule can never fire for that text.
    for i, later in enumerate(rules):
        if later.is_regex:
            continue
        for earlier in rules[:i]:
            if earlier.pattern == later.pattern:
                continue
            if earlier.regex.search(later.pattern):
                warnings.append(
                    f"line {later.lineno} {later.pattern!r} -> {later.group} is "
                    f"shadowed by line {earlier.lineno} {earlier.pattern!r} -> "
                    f"{earlier.group}; move the specific rule above the general one"
                )
                break

    # Groups with a single rule are often typos ("mustelid" vs "mustelidae").
    from collections import Counter
    gc = Counter(r.group for r in rules)
    singles = sorted(g for g, n in gc.items() if n == 1)
    if len(gc) > 3 and singles:
        warnings.append(
            "host groups defined by a single pattern (check for typos): "
            + ", ".join(singles)
        )
    return warnings


def normalize_host(raw: str, rules: list[HostRule]) -> HostMatch:
    """
    Map a free-text /host qualifier onto a canonical name and functional group.
    Unmatched input is returned as ambiguous rather than guessed at, so it
    lands in needs_review.tsv instead of silently becoming 'unknown' in a tree.
    """
    if not raw or not str(raw).strip():
        return HostMatch("", "unknown", True, reason="empty_host_field")
    text = normalize_text(raw)
    if not text:
        return HostMatch("", "unknown", True, reason="empty_host_field")

    # GenBank's /isolation_source is read as a fallback for /host, so clinical
    # sample types leak in. They describe the specimen, not the animal, and
    # must not be matched against the host table or sent to review as if a
    # pattern were missing.
    if SAMPLE_TYPE_RE.fullmatch(text.lower()):
        return HostMatch("", "unknown", True, reason="sample_type_not_a_host")
    for r in rules:
        if r.regex.search(text):
            return HostMatch(r.canonical, r.group, False, r.pattern)
    return HostMatch("", "unknown", True, reason="no_pattern_matched")

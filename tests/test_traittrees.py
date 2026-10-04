"""
Tests for the shared BEAST ancestral-state tree parser, and for the two
analyses that read it.

The parser is the piece both 08f and 08g depend on, so an error here does not
produce a crash — it produces a transition count and a lineage share that are
quietly wrong and agree with each other. The cases below are the ones where a
plausible-looking parser silently returns the wrong thing: a bracket inside a
quoted value, a `.prob` key shadowing the state key, and annotations appearing
on either side of the branch length depending on which logger wrote the file.

The end-to-end checks run both scripts on a posterior whose transitions were
constructed by hand, so the expected counts are known independently rather than
being whatever the code currently emits.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib import traittrees as tt  # noqa: E402


# --- annotation extraction --------------------------------------------------

def test_annotation_value_does_not_match_a_longer_key():
    """
    Asking for `host_group` must not return `host_group.prob`. Unanchored, the
    probability comes back where a state was wanted and every tip silently
    becomes the state "0.93".
    """
    content = '&host_group.prob=0.93,host_group="felid"'
    assert tt.annotation_value(content, "host_group") == "felid"
    assert tt.annotation_value(content, "host_group.prob") == "0.93"


def test_annotation_value_handles_set_valued_and_quoted_forms():
    content = '&host_group="felid",host_group.set={felid,procyonid},rate=1.0'
    assert tt.annotation_value(content, "host_group") == "felid"
    assert tt.annotation_value(content, "host_group.set") == "felid,procyonid"
    assert tt.annotation_value(content, "rate") == "1.0"


def test_annotation_value_returns_none_for_an_absent_key():
    assert tt.annotation_value('&host_group="felid"', "location") is None


def test_read_bracket_is_quote_aware():
    """A ']' inside a quoted value must not end the comment early."""
    s = '[&note="a]b",host_group="felid"]rest'
    content, nxt = tt.read_bracket(s, 0)
    assert tt.annotation_value(content, "host_group") == "felid"
    assert s[nxt:] == "rest"


def test_unterminated_bracket_is_an_error_not_a_truncation():
    with pytest.raises(ValueError, match="unterminated"):
        tt.read_bracket("[&host_group=felid", 0)


# --- the parser -------------------------------------------------------------

SIMPLE = ('((A[&host_group="felid"]:1.0,B[&host_group="felid"]:1.0)'
          '[&host_group="felid"]:2.0,C[&host_group="procyonid"]:3.0)'
          '[&host_group="procyonid"]:0.0;')


def parse(nwk, tag="host_group"):
    return tt.NewickParser(nwk, tag).parse()


def test_parser_builds_a_flat_list_with_parents_before_children():
    nodes = parse(SIMPLE)
    assert len(nodes) == 5
    for i, nd in enumerate(nodes):
        if nd["parent"] is not None:
            assert nd["parent"] < i, "a child precedes its parent"


def test_parser_reads_states_and_tip_flags():
    nodes = parse(SIMPLE)
    tips = [n for n in nodes if n["is_tip"]]
    assert len(tips) == 3
    assert sorted(n["state"] for n in tips) == ["felid", "felid", "procyonid"]
    assert nodes[0]["state"] == "procyonid"          # root


def test_annotation_is_absorbed_before_or_after_the_colon():
    """
    Different loggers put the annotation in different places. Both forms, and
    the one that splits across the colon, have to give the same answer.
    """
    before = '(A[&host_group="felid"]:1.0,B[&host_group="procyonid"]:1.0);'
    after = '(A:[&host_group="felid"]1.0,B:[&host_group="procyonid"]1.0);'
    trailing = '(A:1.0[&host_group="felid"],B:1.0[&host_group="procyonid"]);'
    want = ["felid", "procyonid"]
    for nwk in (before, after, trailing):
        got = [n["state"] for n in parse(nwk) if n["is_tip"]]
        assert got == want, nwk


def test_prob_is_read_when_present_and_none_when_absent():
    """
    The modal-state-only logger leaves .prob out entirely. That must read as
    None rather than as a default, so a downstream filter can tell "unknown"
    from "low".
    """
    with_prob = '(A[&host_group="felid",host_group.prob=0.93]:1.0,B[&host_group="felid"]:1.0);'
    nodes = [n for n in parse(with_prob) if n["is_tip"]]
    assert nodes[0]["prob"] == pytest.approx(0.93)
    assert nodes[1]["prob"] is None


def test_branch_lengths_and_depths():
    nodes = parse(SIMPLE)
    depths = tt.node_depths(nodes)
    assert depths[0] == 0.0
    assert max(depths) == pytest.approx(3.0)


def test_label_date_reads_the_year_from_a_tip_label():
    assert tt.label_date("MT932504|wild_felid|1992.7") == pytest.approx(1992.7)
    assert tt.label_date("no_date_here") is None
    assert tt.label_date(None) is None


def test_detect_tag_finds_the_state_key_and_skips_numeric_ones():
    nwk = '(A[&rate=1.0,posterior=0.99,host_group="felid"]:1.0,B:1.0);'
    assert tt.detect_tag(nwk) == "host_group"


def test_detect_tag_respects_an_explicit_preference():
    nwk = '(A[&location="x",host_group="felid"]:1.0,B:1.0);'
    assert tt.detect_tag(nwk, preferred="host_group") == "host_group"


# --- NEXUS reading ----------------------------------------------------------

def nexus(trees, translate=True):
    out = ["#NEXUS", "Begin trees;"]
    if translate:
        out += ["\tTranslate",
                "\t\t1 MT932504|wild_felid|1992,",
                "\t\t2 MW984525|procyonid|2015",
                ";"]
    out += [f"tree STATE_{i} = [&R] {t}" for i, t in enumerate(trees)]
    out.append("End;")
    return "\n".join(out) + "\n"


def test_translate_block_is_applied_to_tip_labels(tmp_path):
    p = tmp_path / "x.trees"
    p.write_text(nexus(['(1[&host_group="felid"]:1.0,2[&host_group="procyonid"]:1.0);']))
    _, nodes, tag = next(tt.iter_trees(p))
    labels = sorted(n["label"] for n in nodes if n["is_tip"])
    assert labels == ["MT932504|wild_felid|1992", "MW984525|procyonid|2015"]
    assert tag == "host_group"


def test_iter_trees_yields_every_tree(tmp_path):
    p = tmp_path / "x.trees"
    p.write_text(nexus([SIMPLE] * 5, translate=False))
    assert len([None for _ in tt.iter_trees(p)]) == 5


def test_gzipped_tree_files_are_read(tmp_path):
    import gzip
    p = tmp_path / "x.trees.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(nexus([SIMPLE] * 3, translate=False))
    assert len([None for _ in tt.iter_trees(p)]) == 3


# --- end to end, against transitions constructed by hand --------------------

KNOWN = ('((A|wild_felid|1990[&host_group="felid"]:5.0,'
         'B|wild_felid|1992[&host_group="felid"]:3.0)[&host_group="felid"]:10.0,'
         '(C|procyonid|2015[&host_group="procyonid"]:2.0,'
         'D|domestic_dog|2018[&host_group="domestic_dog"]:1.0)'
         '[&host_group="procyonid"]:12.0)[&host_group="felid"]:0.0;')
# Exactly two changes per tree: felid -> procyonid on an internal branch, and
# procyonid -> domestic_dog on a terminal one.


@pytest.fixture
def known_posterior(tmp_path):
    p = tmp_path / "known.trees"
    p.write_text(nexus([KNOWN] * 20, translate=False))
    return p


def test_jump_history_recovers_the_constructed_transitions(known_posterior, tmp_path):
    out = tmp_path / "jumps.tsv"
    r = subprocess.run(
        [sys.executable, "scripts/08f_jump_history.py", "--trees", str(known_posterior),
         "--burnin", "0.0", "--out", str(out)],
        cwd=ROOT, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "mean 2.00" in r.stdout, "two transitions per tree"
    assert "felid -> procyonid" in r.stdout
    assert "procyonid -> domestic_dog" in r.stdout
    rows = [l.split("\t") for l in out.read_text().splitlines()[1:] if l.strip()]
    assert len(rows) == 40, "20 trees x 2 transitions"


def test_a_jump_on_a_terminal_branch_is_flagged_as_terminal(known_posterior, tmp_path):
    """
    A transition on a terminal branch describes one tip's own label and says
    nothing about onward transmission. Reporting it alongside internal
    transitions without that distinction is how a tip census becomes a corridor.
    """
    r = subprocess.run(
        [sys.executable, "scripts/08f_jump_history.py", "--trees", str(known_posterior),
         "--burnin", "0.0", "--out", str(tmp_path / "j.tsv")],
        cwd=ROOT, capture_output=True, text=True, timeout=120)
    line = [l for l in r.stdout.splitlines() if "procyonid -> domestic_dog" in l][0]
    assert "100%" in line.split()[-1], "the dog jump sits on a terminal branch"
    internal = [l for l in r.stdout.splitlines() if "felid -> procyonid" in l][0]
    assert internal.split()[-1] == "0%"


def test_state_through_time_runs_on_the_same_file(known_posterior, tmp_path):
    out = tmp_path / "stt.tsv"
    r = subprocess.run(
        [sys.executable, "scripts/08g_state_through_time.py", "--trees",
         str(known_posterior), "--burnin", "0.0", "--thin", "1", "--out", str(out)],
        cwd=ROOT, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]
    assert out.is_file()
    header = out.read_text().splitlines()[0]
    assert "year" in header.lower()


# --- the regression that prompted the split ---------------------------------

def test_neither_script_imports_the_other_by_filename():
    """
    08g used to load the parser out of 08f with importlib, because a module
    whose name starts with a digit cannot be imported normally. That broke the
    moment the scripts were renamed out of the 12/13 slots. Shared code does
    not belong in a file whose name encodes its pipeline position.
    """
    for name in ("08f_jump_history.py", "08g_state_through_time.py"):
        text = (ROOT / "scripts" / name).read_text()
        assert "import_module" not in text, name
        assert "spec_from_file_location" not in text, name
        assert "from lib.traittrees import" in text, name


def test_both_scripts_and_the_library_agree_on_the_parser():
    """One parser, so a transition count and a lineage share describe one tree."""
    import importlib.util
    for name in ("08f_jump_history.py", "08g_state_through_time.py"):
        spec = importlib.util.spec_from_file_location(
            name[:3], ROOT / "scripts" / name)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        assert mod.iter_trees is tt.iter_trees, name
        assert mod.node_depths is tt.node_depths, name

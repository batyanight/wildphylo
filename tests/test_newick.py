"""
Tests for the shared BEAST ancestral-state tree parser, and for the two
analyses that read it.

The parser is the piece both jumps and state-through-time depend on, so an error here does not
produce a crash — it produces a transition count and a lineage share that are
quietly wrong and agree with each other. The cases below are the ones where a
plausible-looking parser silently returns the wrong thing: a bracket inside a
quoted value, a `.prob` key shadowing the state key, and annotations appearing
on either side of the branch length depending on which logger wrote the file.

The end-to-end checks run both steps (`wildphylo beast jumps` and
`wildphylo beast state-through-time`) on a posterior whose transitions were
constructed by hand, so the expected counts are known independently rather than
being whatever the code currently emits.
"""

from __future__ import annotations

import contextlib
import io
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
from wildphylo import cli  # noqa: E402
from wildphylo import newick as tt  # noqa: E402
from wildphylo.stages import beast  # noqa: E402

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

SIMPLE = (
    '((A[&host_group="felid"]:1.0,B[&host_group="felid"]:1.0)'
    '[&host_group="felid"]:2.0,C[&host_group="procyonid"]:3.0)'
    '[&host_group="procyonid"]:0.0;'
)


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
    assert nodes[0]["state"] == "procyonid"  # root


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
        out += [
            "\tTranslate",
            "\t\t1 MT932504|wild_felid|1992,",
            "\t\t2 MW984525|procyonid|2015",
            ";",
        ]
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


def test_count_trees_agrees_with_iter_trees(tmp_path):
    """One grammar: the count must be the number of trees ``iter_trees`` yields."""
    p = tmp_path / "x.trees"
    p.write_text(nexus([SIMPLE] * 5))
    assert tt.count_trees(p) == 5 == len([None for _ in tt.iter_trees(p)])


def test_count_trees_ignores_everything_that_is_not_a_tree_line(tmp_path):
    p = tmp_path / "x.trees"
    p.write_text(
        "#NEXUS\nBegin taxa;\n\tDimensions ntax=2;\nEnd;\nBegin trees;\n"
        "\tTranslate\n\t\t1 A|felid|1990,\n\t\t2 B|felid|1992\n;\n"
        "tree STATE_0 [&lnP=-1.5] = [&R] (1:1.0,2:1.0);\n"
        "[a comment (1,2);]\n"
        "tree STATE_1 = [&R] (1:1.0,2:1.0);\n"
        "End;\n"
    )
    assert tt.count_trees(p) == 2


def test_count_trees_reads_gzip_and_empty_files(tmp_path):
    import gzip

    gz = tmp_path / "x.trees.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(nexus([SIMPLE] * 3, translate=False))
    assert tt.count_trees(gz) == 3
    empty = tmp_path / "empty.trees"
    empty.write_text("#NEXUS\nBegin trees;\nEnd;\n")
    assert tt.count_trees(empty) == 0


def test_only_newick_knows_the_tree_line_grammar():
    """beast.py counts trees through newick, not with its own copy of the regex."""
    text = Path(beast.__file__).read_text()
    assert "TREE_LINE_RE" not in text
    assert beast.count_trees is tt.count_trees


# --- end to end, against transitions constructed by hand --------------------

KNOWN = (
    '((A|wild_felid|1990[&host_group="felid"]:5.0,'
    'B|wild_felid|1992[&host_group="felid"]:3.0)[&host_group="felid"]:10.0,'
    '(C|procyonid|2015[&host_group="procyonid"]:2.0,'
    'D|domestic_dog|2018[&host_group="domestic_dog"]:1.0)'
    '[&host_group="procyonid"]:12.0)[&host_group="felid"]:0.0;'
)
# Exactly two changes per tree: felid -> procyonid on an internal branch, and
# procyonid -> domestic_dog on a terminal one.


@pytest.fixture
def known_posterior(tmp_path):
    p = tmp_path / "known.trees"
    p.write_text(nexus([KNOWN] * 20, translate=False))
    return p


def run_cli(*argv):
    """Run ``wildphylo <argv>`` in-process; return exit code, stdout and stderr."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = cli.main([str(a) for a in argv])
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
            if isinstance(exc.code, str):
                err.write(exc.code)
    return SimpleNamespace(returncode=code, stdout=out.getvalue(), stderr=err.getvalue())


def test_jump_history_recovers_the_constructed_transitions(known_posterior, tmp_path):
    out = tmp_path / "jumps.tsv"
    r = run_cli(
        "beast",
        "jumps",
        "--trees",
        str(known_posterior),
        "--burnin",
        "0.0",
        "--out",
        str(out),
    )
    assert r.returncode == 0, r.stderr[-2000:]
    assert "mean 2.00" in r.stdout, "two transitions per tree"
    assert "felid -> procyonid" in r.stdout
    assert "procyonid -> domestic_dog" in r.stdout
    rows = [line.split("\t") for line in out.read_text().splitlines()[1:] if line.strip()]
    assert len(rows) == 40, "20 trees x 2 transitions"


def test_a_jump_on_a_terminal_branch_is_flagged_as_terminal(known_posterior, tmp_path):
    """
    A transition on a terminal branch describes one tip's own label and says
    nothing about onward transmission. Reporting it alongside internal
    transitions without that distinction is how a tip census becomes a corridor.
    """
    r = run_cli(
        "beast",
        "jumps",
        "--trees",
        str(known_posterior),
        "--burnin",
        "0.0",
        "--out",
        str(tmp_path / "j.tsv"),
    )
    line = next(line for line in r.stdout.splitlines() if "procyonid -> domestic_dog" in line)
    assert "100%" in line.split()[-1], "the dog jump sits on a terminal branch"
    internal = next(line for line in r.stdout.splitlines() if "felid -> procyonid" in line)
    assert internal.split()[-1] == "0%"


def test_state_through_time_runs_on_the_same_file(known_posterior, tmp_path):
    out = tmp_path / "stt.tsv"
    r = run_cli(
        "beast",
        "state-through-time",
        "--trees",
        str(known_posterior),
        "--burnin",
        "0.0",
        "--thin",
        "1",
        "--out",
        str(out),
    )
    assert r.returncode == 0, r.stderr[-2000:]
    assert out.is_file()
    header = out.read_text().splitlines()[0]
    assert "year" in header.lower()


# --- the regression that prompted the split ---------------------------------


def test_both_steps_use_the_one_shared_parser():
    """
    The 08g script used to load the parser out of 08f with importlib, because a
    module whose name starts with a digit cannot be imported normally. That
    broke the moment the scripts were renamed out of the 12/13 slots. Shared
    code does not belong in a file whose name encodes its pipeline position, and
    one parser means a transition count and a lineage share describe one tree.
    """
    text = Path(beast.__file__).read_text()
    assert "import_module" not in text
    assert "spec_from_file_location" not in text
    assert beast.iter_trees is tt.iter_trees
    assert beast.node_depths is tt.node_depths


# --- one parser, two views --------------------------------------------------
#
# The jump/state steps walk the node list; the trait signal gate consumes the
# Topology arrays. Both are built by the same parser, so the same tree must
# come out as the same shape, labels and states from either entry point.

TRAIT_NEXUS = """#NEXUS
Begin trees;
\tTranslate
\t\t1 MT932504|wild_felid|1992,
\t\t2 MT932511|wild_felid|1992,
\t\t3 MW984525|procyonid|2015,
\t\t4 KX024709|domestic_dog|2004
;
tree STATE_0 = [&R] ((1[&host_group="felid"]:0.5,2[&host_group="felid"]:0.5):0.3,(3:0.2,4:0.2):0.6):0.0;
tree STATE_1 = [&R] ((1:0.5,3:0.5):0.3,(2:0.2,4:0.2):0.6):0.0;
End;
"""  # noqa: E501 -- Nexus fixture; tree lines cannot wrap


def _node_clades(nodes):
    """The set of tip-label sets, one per internal node, from a node list."""
    below = [frozenset() for _ in nodes]
    for i in reversed(range(len(nodes))):  # children come after parents
        nd = nodes[i]
        if nd["is_tip"]:
            below[i] = frozenset([nd["label"]])
        if nd["parent"] is not None:
            below[nd["parent"]] = below[nd["parent"]] | below[i]
    return {below[i] for i, nd in enumerate(nodes) if not nd["is_tip"]}


def _topology_clades(top):
    """The same set of clades, read from the Topology arrays."""
    name = {v: k for k, v in top.tip_index.items()}
    below = {}
    clades = set()
    for raw in top.postorder:
        node = int(raw)
        left, right = int(top.left[node]), int(top.right[node])
        if left < 0:
            below[node] = frozenset([name[node]])
        else:
            below[node] = below[left] | below[right]
            clades.add(below[node])
    return clades


@pytest.mark.parametrize("which", [0, 1])
def test_both_views_of_a_nexus_tree_agree(tmp_path, which):
    p = tmp_path / "x.trees"
    p.write_text(TRAIT_NEXUS)
    _, nodes, _ = list(tt.iter_trees(p))[which]
    top = list(tt.iter_topologies(p, burnin=0.0))[which]
    assert sorted(top.tip_index) == sorted(n["label"] for n in nodes if n["is_tip"])
    assert _topology_clades(top) == _node_clades(nodes)
    assert top.n_tips == sum(n["is_tip"] for n in nodes)


def test_translated_labels_and_states_come_from_one_parse(tmp_path):
    p = tmp_path / "x.trees"
    p.write_text(TRAIT_NEXUS)
    _, nodes, tag = next(tt.iter_trees(p))
    top = next(tt.iter_topologies(p, burnin=0.0))
    assert tag == "host_group"
    states = {n["label"]: n["state"] for n in nodes if n["is_tip"]}
    assert states == {
        "MT932504|wild_felid|1992": "felid",
        "MT932511|wild_felid|1992": "felid",
        "MW984525|procyonid|2015": None,
        "KX024709|domestic_dog|2004": None,
    }
    assert set(top.tip_index) == set(states)


def test_reading_annotations_does_not_change_the_topology():
    plain = tt.parse_topology(tt.strip_comments(KNOWN))
    annotated = tt.parse_topology(KNOWN)
    assert annotated.tip_index == plain.tip_index
    assert annotated.postorder.tolist() == plain.postorder.tolist()
    assert annotated.left.tolist() == plain.left.tolist()
    assert annotated.right.tolist() == plain.right.tolist()
    # ... and the node view of the same string still has the states
    states = [n["state"] for n in tt.parse_newick(KNOWN, "host_group") if n["is_tip"]]
    assert states == ["felid", "felid", "procyonid", "domestic_dog"]


def test_a_polytomy_is_one_node_in_the_list_and_counted_in_the_topology():
    nk = "(A:1,B:1,C:1,D:1);"
    nodes = tt.parse_newick(nk)
    assert len(nodes[0]["children"]) == 4
    top = tt.parse_topology(nk)
    assert top.n_polytomies == 1
    assert top.n_tips == 4


# --- --build-dir defaults ---------------------------------------------------


def _trait_build(tmp_path):
    """A build dir whose two chains (beast.seeds of the config) each hold the known posterior."""
    import yaml

    from wildphylo.layout import BuildLayout

    cfg_path = ROOT / "tests/fixtures/configs/cdv-1200.yaml"
    cfg = yaml.safe_load(cfg_path.read_text())
    build = tmp_path / "builds" / cfg["id"] / "20260101"
    layout = BuildLayout.at(build, cfg)
    for seed in cfg["beast"]["seeds"]:
        f = layout.beast_trait_trees("H", str(seed))
        f.parent.mkdir(parents=True)
        f.write_text(nexus([KNOWN] * 20, translate=False))
    return cfg_path, build, layout


def test_jumps_default_to_every_seed_under_build_dir(tmp_path):
    cfg, build, _ = _trait_build(tmp_path)
    out = tmp_path / "jumps.tsv"
    r = run_cli(
        "beast",
        *("jumps", "--config", cfg, "--build-dir", build, "--locus", "H"),
        *("--burnin", "0.0", "--out", out),
    )
    assert r.returncode == 0, r.stderr[-2000:]
    rows = [line for line in out.read_text().splitlines()[1:] if line.strip()]
    assert len(rows) == 80, "2 seeds x 20 trees x 2 transitions"


def test_jumps_explicit_trees_override_build_dir(tmp_path):
    cfg, build, layout = _trait_build(tmp_path)
    out = tmp_path / "jumps.tsv"
    r = run_cli(
        "beast",
        *("jumps", "--config", cfg, "--build-dir", build, "--locus", "H"),
        *("--trees", layout.beast_trait_trees("H", "12345")),
        *("--burnin", "0.0", "--out", out),
    )
    assert r.returncode == 0, r.stderr[-2000:]
    rows = [line for line in out.read_text().splitlines()[1:] if line.strip()]
    assert len(rows) == 40, "one seed only"


def test_state_through_time_defaults_to_every_seed_under_build_dir(tmp_path):
    cfg, build, _ = _trait_build(tmp_path)
    out = tmp_path / "stt.tsv"
    r = run_cli(
        "beast",
        *("state-through-time", "--config", cfg, "--build-dir", build, "--locus", "H"),
        *("--burnin", "0.0", "--thin", "1", "--out", out),
    )
    assert r.returncode == 0, r.stderr[-2000:]
    assert "year" in out.read_text().splitlines()[0].lower()


def test_build_dir_without_locus_or_config_is_a_clear_error(tmp_path, caplog):
    caplog.set_level(logging.ERROR, logger="wildphylo")
    cfg, build, _ = _trait_build(tmp_path)
    no_locus = run_cli("beast", "jumps", "--config", cfg, "--build-dir", build)
    assert no_locus.returncode == 2
    assert "--locus" in caplog.text
    caplog.clear()
    no_cfg = run_cli("beast", "jumps", "--build-dir", build, "--locus", "H")
    assert no_cfg.returncode == 2
    assert "--config" in caplog.text


def test_jumps_without_trees_or_build_dir_is_a_usage_error(caplog):
    caplog.set_level(logging.ERROR, logger="wildphylo")
    r = run_cli("beast", "jumps")
    assert r.returncode == 2
    assert "--trees is required" in caplog.text


def test_a_malformed_focus_is_a_usage_error(tmp_path, caplog):
    caplog.set_level(logging.ERROR, logger="wildphylo")
    _, _, layout = _trait_build(tmp_path)
    trees = layout.beast_trait_trees("H", "12345")
    r = run_cli("beast", "jumps", "--trees", trees, "--burnin", "0.0", "--focus", "no_colon")
    assert r.returncode == 2
    assert "FROM:TO" in caplog.text

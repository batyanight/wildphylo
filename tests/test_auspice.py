"""
Tests for the Auspice export.

The failures that matter here are the ones that produce a VALID JSON which
renders normally and says something false: a build labelled with another
pathogen's title, a time axis offset so every node sits in the wrong place, or
a trait colouring where tips and internal nodes are in different categories.
None of those raise anything.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location(
        "aus09", ROOT / "scripts" / "09_make_auspice.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


MCC = """#NEXUS
Begin trees;
\tTranslate
\t\t1 MT932504|wild_felid|1992.0,
\t\t2 MW984525|procyonid|2015.0,
\t\t3 KX024709|domestic_dog|2020.0
;
tree TREE1 = [&R] ((1[&host_group="felid",host_group.prob=0.95,height=28.0]:12.0,\
2[&host_group="procyonid",height=5.0]:11.0)\
[&host_group="felid",height=16.0,height_95%_HPD={12.0,22.0},posterior=0.88]:14.0,\
3[&host_group="domestic_dog",height=0.0]:30.0)\
[&host_group="felid",height=30.0,posterior=1.0]:0.0;
End;
"""

META = ("label\taccession\thost_group\thost_canonical\n"
        "MT932504|wild_felid|1992.0\tMT932504\twild_felid\tLynx rufus\n"
        "MW984525|procyonid|2015.0\tMW984525\tprocyonid\tProcyon lotor\n"
        "KX024709|domestic_dog|2020.0\tKX024709\tdomestic_dog\tCanis familiaris\n")


@pytest.fixture
def build(tmp_path):
    (tmp_path / "mcc.tree").write_text(MCC)
    (tmp_path / "meta.tsv").write_text(META)
    return tmp_path


def run(build, out, *extra, config=None):
    cmd = [sys.executable, "scripts/09_make_auspice.py",
           "--mcc", str(build / "mcc.tree"),
           "--metadata", str(build / "meta.tsv"),
           "--output", str(out), *extra]
    if config:
        cmd += ["--config", str(config)]
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=120)


# --- the time axis ----------------------------------------------------------

def test_most_recent_is_derived_from_the_tip_labels(build, tmp_path):
    """
    Passing it by hand is how the axis gets silently offset: a value a year
    stale shifts every node together, so nothing looks wrong.
    """
    out = tmp_path / "a.json"
    r = run(build, out)
    assert r.returncode == 0, r.stderr[-1500:]
    assert "2020.0000" in r.stdout
    doc = json.loads(out.read_text())
    assert doc["tree"]["node_attrs"]["num_date"]["value"] == pytest.approx(1990.0)


def test_an_inconsistent_time_axis_refuses_to_write(build, tmp_path):
    """
    A wrong axis produces a figure that looks entirely normal, with every node
    in the wrong place. It must refuse rather than warn.
    """
    out = tmp_path / "a.json"
    r = run(build, out, "--most-recent", "2015.0")
    assert r.returncode == 1
    assert "time axis is wrong" in r.stderr
    assert not out.exists(), "nothing may be written when the axis is wrong"


def test_height_hpd_becomes_a_date_confidence_interval(build, tmp_path):
    out = tmp_path / "a.json"
    assert run(build, out).returncode == 0
    node = json.loads(out.read_text())["tree"]["children"][0]
    conf = node["node_attrs"]["num_date"]["confidence"]
    # heights 12 and 22 below a most-recent of 2020 -> 2008 and 1998
    assert conf == [pytest.approx(1998.0), pytest.approx(2008.0)]
    assert conf[0] < conf[1], "the interval must be ordered low to high"


# --- nothing names the pathogen unless the config does ----------------------

def test_the_title_comes_from_the_config(build, tmp_path):
    out = tmp_path / "a.json"
    r = run(build, out, config=ROOT / "config" / "pathogen" / "cdv-1200.yaml")
    assert r.returncode == 0, r.stderr[-1500:]
    meta = json.loads(out.read_text())["meta"]
    assert "CDV" in meta["title"]


def test_no_pathogen_name_is_hardcoded(build, tmp_path):
    """
    Regression. The script shipped a default title naming the America-2 CDV
    lineage and a description citing A75/17. In a pathogen-agnostic repo that
    publishes a BTV build described as canine distemper, and the JSON is valid
    and renders, so nothing downstream catches it.
    """
    out = tmp_path / "a.json"
    assert run(build, out).returncode == 0
    doc = out.read_text().lower()
    for word in ("america-2", "canine distemper", "a75/17", "af164967"):
        assert word not in doc, f"{word!r} leaked into a build with no config"
    # Check the executable code, not the module docstring: the docstring
    # quotes the removed strings in order to explain why they were removed,
    # and a test that forbids naming the bug makes the bug harder to fix.
    import ast
    src = (ROOT / "scripts" / "09_make_auspice.py").read_text()
    tree = ast.parse(src)
    body = tree.body[1:] if (tree.body and isinstance(tree.body[0], ast.Expr)
                             and isinstance(tree.body[0].value, ast.Constant)) \
        else tree.body
    code = "\n".join(ast.unparse(n) for n in body).lower()
    for word in ("a75/17", "america-2", "af164967"):
        assert word not in code, f"{word!r} is still hardcoded in the code"


def test_colourings_come_from_the_config(build, tmp_path):
    out = tmp_path / "a.json"
    run(build, out, config=ROOT / "config" / "pathogen" / "cdv-1200.yaml")
    keys = [c["key"] for c in json.loads(out.read_text())["meta"]["colorings"]]
    assert "country" in keys and "lineage" in keys
    assert keys[0] == "host_group", "the trait is the default colouring"


# --- the trait --------------------------------------------------------------

def test_tip_and_internal_states_share_a_category(build, tmp_path):
    """
    Regression. Tips took their host from metadata, which holds the RAW host
    group, while internal nodes carry the model's collapsed state. Auspice then
    renders wild_felid and felid as different colours and the tree implies a
    transition at every tip.
    """
    out = tmp_path / "a.json"
    run(build, out, config=ROOT / "config" / "pathogen" / "cdv-1200.yaml")
    tree = json.loads(out.read_text())["tree"]
    internal = tree["children"][0]["node_attrs"]["host_group"]["value"]
    tip = tree["children"][0]["children"][0]["node_attrs"]["host_group"]["value"]
    assert tip == internal == "felid"


def test_the_trait_key_is_auto_detected(mod, build):
    """The XML writes tag=host_group; the old default was 'location'."""
    nodes, tag = mod.parse_mcc(build / "mcc.tree")
    assert tag == "host_group"


def test_state_probability_becomes_a_confidence(build, tmp_path):
    out = tmp_path / "a.json"
    run(build, out, config=ROOT / "config" / "pathogen" / "cdv-1200.yaml")
    tree = json.loads(out.read_text())["tree"]
    tip = tree["children"][0]["children"][0]["node_attrs"]
    assert "species" in tip


# --- gate verdicts ----------------------------------------------------------

def test_gate_verdicts_appear_in_the_description(build, tmp_path):
    """
    A published tree is read by people who did not run the pipeline, and a
    trait colouring that failed its signal test renders identically to one
    that passed.
    """
    g1 = tmp_path / "trait_signal.json"
    g1.write_text(json.dumps({"verdict": "pass", "reasons": ["AI p=0.001"]}))
    g2 = tmp_path / "convergence.json"
    g2.write_text(json.dumps({"verdict": "fail",
                              "checks": {"errors": ["ESS below 200"]}}))
    out = tmp_path / "a.json"
    run(build, out, "--gates", str(g1), str(g2))
    desc = json.loads(out.read_text())["meta"]["description"]
    assert "trait signal" in desc and "**pass**" in desc
    assert "convergence" in desc and "**fail**" in desc
    assert "did not pass" in desc


def test_a_build_with_no_failed_gates_gets_no_warning(build, tmp_path):
    g = tmp_path / "trait_signal.json"
    g.write_text(json.dumps({"verdict": "pass", "reasons": ["AI p=0.001"]}))
    out = tmp_path / "a.json"
    run(build, out, "--gates", str(g))
    assert "did not pass" not in json.loads(out.read_text())["meta"]["description"]


def test_a_malformed_gate_file_is_skipped_not_fatal(mod, tmp_path):
    bad = tmp_path / "broken.json"
    bad.write_text("{not json")
    ok = tmp_path / "trait_signal.json"
    ok.write_text(json.dumps({"verdict": "pass"}))
    got = mod.gate_summary([bad, ok, tmp_path / "absent.json"])
    assert [g[1] for g in got] == ["pass"]


# --- the workflow's own invocation ------------------------------------------

def test_the_snakefile_flags_are_all_accepted(build, tmp_path):
    """
    Regression. The Snakefile passed --config/--alignment/--output while the
    script wanted --title/--maintainer/--most-recent/--trait-key, so the rule
    would have failed at runtime, potentially days into a build.
    """
    aln = tmp_path / "sub.fasta"
    aln.write_text(">a\nACGT\n>b\nACGT\n")
    out = tmp_path / "a.json"
    r = run(build, out, "--alignment", str(aln),
            "--full-metadata", str(build / "meta.tsv"),
            config=ROOT / "config" / "pathogen" / "cdv-1200.yaml")
    assert r.returncode == 0, r.stderr[-1500:]
    assert "2 sequences" in r.stdout


def test_the_snakefile_passes_only_flags_the_script_has():
    """Every --flag in the auspice rule must exist in the parser."""
    import re
    sf = (ROOT / "workflow" / "Snakefile").read_text()
    rule = sf.split("rule auspice:")[1].split("rule publish:")[0]
    called = set(re.findall(r"(--[a-z-]+)", rule))
    src = (ROOT / "scripts" / "09_make_auspice.py").read_text()
    declared = set(re.findall(r'ap\.add_argument\("(--[a-z-]+)"', src))
    assert called <= declared, f"Snakefile passes unknown flags: {called - declared}"


def test_cli_help_works():
    r = subprocess.run([sys.executable, "scripts/09_make_auspice.py", "--help"],
                       cwd=ROOT, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0

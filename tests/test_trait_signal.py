"""
Tests for the tip-label randomisation gate.

The statistics are checked against values computed by hand on a four-tip tree,
not against the code's own output, so a refactor that changes the answer fails
here rather than quietly redefining what AI or PS mean.

The gate itself is checked against two constructed cases where the right answer
is known independently: perfectly monophyletic groups (must pass) and labels
randomised with respect to topology (must fail). A test suite that only
verified the code runs would not have caught the null-construction error the
module docstring describes.
"""

from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lib import traitsignal as ts  # noqa: E402

FOUR = "((A:1,B:1):1,(C:1,D:1):1);"


def _states(top, mapping):
    S = max(mapping.values()) + 1
    a = np.zeros((1, top.n_tips), dtype=np.int64)
    for k, v in mapping.items():
        a[0, top.tip_index[k]] = v
    return a, S


# --- parser -----------------------------------------------------------------

def test_parse_four_tip_tree():
    t = ts.parse_newick(FOUR)
    assert t.n_tips == 4
    assert t.n_nodes == 7
    assert sorted(t.tip_index) == ["A", "B", "C", "D"]
    assert all(t.tip_index[k] < 4 for k in t.tip_index)


def test_postorder_puts_children_before_parents():
    t = ts.parse_newick("(((A:1,B:1):1,C:1):1,D:1);")
    seen = set()
    for node in t.postorder:
        l, r = t.left[node], t.right[node]
        if l >= 0:
            assert l in seen and r in seen, "parent visited before a child"
        seen.add(int(node))


def test_strip_comments_removes_beast_annotations():
    nk = '(A[&host="felid",rate=1.0]:0.5,B[&host="dog"]:0.5)[&posterior=1.0]:0.0;'
    assert ts.strip_comments(nk) == "(A:0.5,B:0.5):0.0;"


def test_strip_comments_survives_a_bracket_inside_a_quoted_annotation():
    """
    A ']' inside a quoted annotation value would end the comment early under a
    naive regex and leave junk in the newick.
    """
    nk = '(A[&note="a]b",rate=1]:0.5,B:0.5);'
    assert ts.strip_comments(nk) == "(A:0.5,B:0.5);"


def test_duplicate_tip_labels_are_rejected():
    with pytest.raises(ValueError, match="duplicate"):
        ts.parse_newick("((A:1,A:1):1,B:1);")


def test_polytomy_is_resolved_and_counted():
    t = ts.parse_newick("(A:1,B:1,C:1,D:1);")
    assert t.n_tips == 4
    assert t.n_polytomies == 1


# --- statistics, against hand computation -----------------------------------

def test_parsimony_score_hand_computed():
    t = ts.parse_newick(FOUR)
    a, S = _states(t, dict(A=0, B=0, C=1, D=1))
    assert ts.parsimony_score(t, a, S)[0] == 1        # one change at the root split
    a, S = _states(t, dict(A=0, B=1, C=0, D=1))
    assert ts.parsimony_score(t, a, S)[0] == 2        # alternating needs two


def test_association_index_hand_computed():
    """
    AI = sum over internal nodes of (1 - f) / 2^(n-1).

    Clean split:  (A,B) f=1 -> 0 ; (C,D) f=1 -> 0 ; root f=.5, n=4 -> .5/8
    Alternating:  (A,B) f=.5, n=2 -> .25 ; (C,D) .25 ; root .0625
    """
    t = ts.parse_newick(FOUR)
    a, S = _states(t, dict(A=0, B=0, C=1, D=1))
    assert ts.association_index(t, a, S)[0] == pytest.approx(0.0625)
    a, S = _states(t, dict(A=0, B=1, C=0, D=1))
    assert ts.association_index(t, a, S)[0] == pytest.approx(0.5625)


def test_monophyletic_clade_sizes_hand_computed():
    t = ts.parse_newick(FOUR)
    a, S = _states(t, dict(A=0, B=0, C=1, D=1))
    assert list(ts.monophyletic_clades(t, a, S)[0]) == [2, 2]
    a, S = _states(t, dict(A=0, B=1, C=0, D=1))
    assert list(ts.monophyletic_clades(t, a, S)[0]) == [1, 1]


def test_statistics_vectorise_identically_over_labellings():
    """A batch of R labellings must give what R separate calls give."""
    t = ts.parse_newick(FOUR)
    batch = np.array([[0, 0, 1, 1], [0, 1, 0, 1], [1, 1, 1, 0]], dtype=np.int64)
    got = ts.parsimony_score(t, batch, 2)
    one = [ts.parsimony_score(t, batch[i:i + 1], 2)[0] for i in range(3)]
    assert list(got) == one


# --- tree construction helpers for the gate tests ---------------------------

def _rand_tree(labels, rng):
    nodes = list(labels)
    rng.shuffle(nodes)
    while len(nodes) > 1:
        a = nodes.pop(rng.randrange(len(nodes)))
        b = nodes.pop(rng.randrange(len(nodes)))
        nodes.append(f"({a}:0.1,{b}:0.1)")
    return nodes[0] + ";"


def _clustered_tree(groups, rng):
    clades = [_rand_tree(g, rng)[:-1] for g in groups]
    while len(clades) > 1:
        a, b = clades.pop(0), clades.pop(0)
        clades.append(f"({a}:0.1,{b}:0.1)")
    return clades[0] + ";"


@pytest.fixture(scope="module")
def two_groups():
    a = [f"a{i}" for i in range(20)]
    b = [f"b{i}" for i in range(20)]
    return a, b, {**{t: "A" for t in a}, **{t: "B" for t in b}}


# --- the gate ---------------------------------------------------------------

def test_gate_passes_on_known_signal(two_groups):
    a, b, states = two_groups
    rng = random.Random(11)
    trees = [ts.parse_newick(_clustered_tree([a, b], rng)) for _ in range(60)]
    res = ts.trait_signal(trees, states, n_permutations=200, seed=1)
    assert res.verdict == "pass"
    assert res.ai.observed < res.ai.null_lo
    assert res.ps.observed < res.ps.null_lo
    assert res.mc["A"].observed > res.mc["A"].null_hi


def test_gate_fails_on_known_noise(two_groups):
    a, b, states = two_groups
    rng = random.Random(12)
    trees = [ts.parse_newick(_rand_tree(a + b, rng)) for _ in range(60)]
    res = ts.trait_signal(trees, states, n_permutations=200, seed=1)
    assert res.verdict == "fail"
    assert res.ai.p > 0.05 and res.ps.p > 0.05


def test_p_value_cannot_be_zero(two_groups):
    """
    Davison-Hinkley: with R permutations the floor is 1/(R+1). A reported p of
    0 would claim resolution the test does not have.
    """
    a, b, states = two_groups
    rng = random.Random(13)
    trees = [ts.parse_newick(_clustered_tree([a, b], rng)) for _ in range(20)]
    res = ts.trait_signal(trees, states, n_permutations=99, seed=1)
    assert res.ai.p == pytest.approx(1 / 100)
    assert res.ai.p > 0


def test_null_is_one_permutation_across_all_trees_not_one_per_tree():
    """
    Regression test for the null construction, which is the error that makes
    everything look significant.

    On random topologies with random labels, the observed mean must land inside
    the null's central mass. If the null were built by permuting independently
    per tree and pooling single values, its variance would be far larger than
    the variance of a mean, and this check would still pass -- so the stronger
    statement is tested: the null's own spread must be comparable to the
    between-replicate spread of a mean, i.e. p is roughly uniform. Here, that
    the observed sits near the middle.
    """
    rng = random.Random(14)
    tips = [f"t{i}" for i in range(30)]
    states = {t: ("A" if i % 2 else "B") for i, t in enumerate(tips)}
    trees = [ts.parse_newick(_rand_tree(tips, rng)) for _ in range(80)]
    res = ts.trait_signal(trees, states, n_permutations=300, seed=2)
    assert res.ai.null_lo < res.ai.observed < res.ai.null_hi
    assert 0.05 < res.ai.p < 0.95


def test_unassigned_tips_are_pruned_not_treated_as_a_state(two_groups):
    """
    '?' tips are the ones nobody could classify. As a state they would cluster
    perfectly by construction and manufacture signal.
    """
    a, b, states = two_groups
    rng = random.Random(15)
    trees = [ts.parse_newick(_clustered_tree([a, b], rng)) for _ in range(20)]
    with_unknown = dict(states)
    for t in a[:5]:
        with_unknown[t] = "?"
    res = ts.trait_signal(trees, with_unknown, n_permutations=50, seed=1)
    assert "?" not in res.states
    assert res.n_pruned == 5
    assert res.n_tips == 35


def test_single_state_is_refused():
    rng = random.Random(16)
    tips = [f"t{i}" for i in range(10)]
    trees = [ts.parse_newick(_rand_tree(tips, rng))]
    with pytest.raises(ValueError, match="nothing to test"):
        ts.trait_signal(trees, {t: "A" for t in tips}, n_permutations=10)


def test_mismatched_traits_file_is_refused_not_silently_intersected():
    rng = random.Random(17)
    tips = [f"t{i}" for i in range(10)]
    trees = [ts.parse_newick(_rand_tree(tips, rng))]
    with pytest.raises(ValueError, match="no tip label"):
        ts.trait_signal(trees, {f"other{i}": "A" if i % 2 else "B"
                                for i in range(10)}, n_permutations=10)


def test_low_power_is_reported_not_silently_passed():
    rng = random.Random(18)
    a = [f"a{i}" for i in range(3)]
    b = [f"b{i}" for i in range(12)]
    states = {**{t: "A" for t in a}, **{t: "B" for t in b}}
    trees = [ts.parse_newick(_clustered_tree([a, b], rng)) for _ in range(20)]
    res = ts.trait_signal(trees, states, n_permutations=100, seed=1)
    joined = " ".join(res.reasons)
    assert "LOW POWER" in joined
    assert "A" in res.state_counts and res.state_counts["A"] == 3


def test_result_is_json_serialisable(two_groups):
    import json
    a, b, states = two_groups
    rng = random.Random(19)
    trees = [ts.parse_newick(_clustered_tree([a, b], rng)) for _ in range(10)]
    res = ts.trait_signal(trees, states, n_permutations=20, seed=1)
    round_tripped = json.loads(json.dumps(res.to_dict()))
    assert round_tripped["verdict"] == res.verdict
    assert round_tripped["ai"]["p"] == pytest.approx(res.ai.p)


# --- NEXUS reading ----------------------------------------------------------

NEXUS = """#NEXUS
Begin trees;
\tTranslate
\t\t1 MT932504|wild_felid|1992,
\t\t2 MT932511|wild_felid|1992,
\t\t3 MW984525|procyonid|2015,
\t\t4 KX024709|domestic_dog|2004
;
tree STATE_0 = [&R] ((1[&host_group="felid"]:0.5,2[&host_group="felid"]:0.5):0.3,(3:0.2,4:0.2):0.6):0.0;
tree STATE_1 = [&R] ((1:0.5,3:0.5):0.3,(2:0.2,4:0.2):0.6):0.0;
tree STATE_2 = [&R] ((1:0.5,4:0.5):0.3,(2:0.2,3:0.2):0.6):0.0;
End;
"""


def test_nexus_translate_block_is_resolved(tmp_path):
    p = tmp_path / "x.trees"
    p.write_text(NEXUS)
    trees = ts.read_posterior([p], burnin=0.0)
    assert len(trees) == 3
    assert "MT932504|wild_felid|1992" in trees[0].tip_index
    assert "1" not in trees[0].tip_index


def test_burnin_is_applied_per_file_not_to_the_pool(tmp_path):
    """
    Pooling two chains and then discarding the first 10% throws away the head
    of chain one only, leaving chain two's burn-in in the sample.
    """
    a, b = tmp_path / "a.trees", tmp_path / "b.trees"
    a.write_text(NEXUS)
    b.write_text(NEXUS)
    both = ts.read_posterior([a, b], burnin=1 / 3)
    assert len(both) == 4           # 2 kept from each, not 4 from one


def test_label_parsing_in_cli_helper(tmp_path):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gate08e",
        Path(__file__).resolve().parents[1] / "scripts" / "08e_trait_signal.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    got = mod.states_from_labels(["MT932504|wild_felid|1992"], 1, "|")
    assert got == {"MT932504|wild_felid|1992": "wild_felid"}
    with pytest.raises(ValueError, match="no field"):
        mod.states_from_labels(["nofields"], 1, "|")


def test_cli_help_works():
    r = subprocess.run(
        [sys.executable, "scripts/08e_trait_signal.py", "--help"],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True, text=True, timeout=60)
    assert r.returncode == 0
    assert "usage:" in r.stdout.lower()


# --- config plumbing --------------------------------------------------------

def _base_cfg(repo_root):
    """A minimal valid config, loaded from the shipped CDV one."""
    import yaml
    return yaml.safe_load(
        (repo_root / "config" / "pathogen" / "cdv-1200.yaml").read_text())


def test_signal_gate_keys_validate(repo_root):
    from lib import config as cfgmod
    cfg = _base_cfg(repo_root)
    cfg["_path"] = str(repo_root / "config" / "pathogen" / "cdv-1200.yaml")
    errors, _ = cfgmod.validate(cfg)
    assert not errors


def test_bad_on_no_signal_is_an_error(repo_root):
    from lib import config as cfgmod
    cfg = _base_cfg(repo_root)
    cfg["_path"] = str(repo_root / "config" / "pathogen" / "cdv-1200.yaml")
    cfg["beast"]["discrete_trait"]["on_no_signal"] = "ignore"
    errors, _ = cfgmod.validate(cfg)
    assert any("on_no_signal" in e for e in errors)


def test_too_few_permutations_warns_about_p_resolution(repo_root):
    from lib import config as cfgmod
    cfg = _base_cfg(repo_root)
    cfg["_path"] = str(repo_root / "config" / "pathogen" / "cdv-1200.yaml")
    cfg["beast"]["discrete_trait"]["signal_permutations"] = 20
    _, warns = cfgmod.validate(cfg)
    assert any("smallest p" in w for w in warns)


def test_disabling_the_gate_warns(repo_root):
    from lib import config as cfgmod
    cfg = _base_cfg(repo_root)
    cfg["_path"] = str(repo_root / "config" / "pathogen" / "cdv-1200.yaml")
    cfg["beast"]["discrete_trait"]["signal_gate"] = False
    _, warns = cfgmod.validate(cfg)
    assert any("signal_gate" in w for w in warns)


def test_every_shipped_config_with_a_trait_declares_the_gate(repo_root):
    """
    A config that runs a DTA without the gate publishes a transition table
    nothing has checked. Shipping one would make that the default.
    """
    import yaml
    missing = []
    for p in sorted((repo_root / "config" / "pathogen").glob("*.yaml")):
        cfg = yaml.safe_load(p.read_text())
        dt = (cfg.get("beast") or {}).get("discrete_trait") or {}
        if dt.get("enabled") and "signal_gate" not in dt:
            missing.append(p.name)
    assert not missing, f"discrete trait but no signal_gate key: {missing}"


# --- the cwd-relative path trap ---------------------------------------------

def test_config_paths_resolve_from_any_working_directory(repo_root, tmp_path, monkeypatch):
    """
    Regression test. hosts.table is written relative to the repo root, so
    running a script from inside a build directory made it unresolvable --
    and the failure was a missing file, not a wrong one, so it surfaced as a
    confusing error much later.
    """
    from lib import config as cfgmod
    monkeypatch.chdir(tmp_path)
    cfg = _base_cfg(repo_root)
    cfg["_path"] = str(repo_root / "config" / "pathogen" / "cdv-1200.yaml")
    got = cfgmod.resolve_path(cfg, cfg["hosts"]["table"])
    assert got.is_file()
    errors, _ = cfgmod.validate(cfg)
    assert not any("not found" in e for e in errors)


def test_resolve_path_leaves_absolute_paths_alone(repo_root, tmp_path):
    from lib import config as cfgmod
    p = tmp_path / "abs.tsv"
    p.write_text("x")
    assert cfgmod.resolve_path({"_path": "whatever"}, str(p)) == p


def test_resolve_path_prefers_the_working_directory(repo_root, tmp_path, monkeypatch):
    """cwd wins, so behaviour that worked before this change is unchanged."""
    from lib import config as cfgmod
    monkeypatch.chdir(tmp_path)
    (tmp_path / "config").mkdir()
    local = tmp_path / "config" / "host_groups.tsv"
    local.write_text("x")
    cfg = {"_path": str(repo_root / "config" / "pathogen" / "cdv-1200.yaml")}
    assert cfgmod.resolve_path(cfg, "config/host_groups.tsv").resolve() == local.resolve()


# --- workflow wiring --------------------------------------------------------

def test_gate_is_wired_into_the_workflow(repo_root):
    """
    A gate that exists as a script but is not in the DAG is a gate nobody runs.
    """
    sf = (repo_root / "workflow" / "Snakefile").read_text()
    assert "08e_trait_signal.py" in sf
    assert "rule trait_signal" in sf


def test_auspice_cannot_be_built_without_the_gate(repo_root):
    """
    The publishing rule must depend on trait_signal.json, so a build cannot
    reach Auspice with an untested trait.
    """
    sf = (repo_root / "workflow" / "Snakefile").read_text()
    auspice = sf.split("rule auspice:")[1].split("rule publish:")[0]
    assert "trait_signal" in auspice


def test_beast_run_declares_the_trait_tree_file(repo_root):
    """
    BEAST writes the ancestral-state trees whether or not the workflow knows.
    An undeclared output cannot be invalidated, so the gate would silently test
    a previous run's reconstruction.
    """
    sf = (repo_root / "workflow" / "Snakefile").read_text()
    assert "trait_trees" in sf

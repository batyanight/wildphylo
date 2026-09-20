"""
Tests for 05c_cut_clades.py.

Built on synthetic trees with known structure and a known substitution rate, so
the expected answer exists independently of the implementation. The motivating
case is real: root-to-tip on the full CDV H tree gave R^2 = 0.024 and an implied
root of 1656, because the tree is saturated across lineages. Signal lives within
clades, and this script has to find them.
"""

import random
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

pytest.importorskip("Bio")
RATE = 8.5e-4


def clade(prefix, n, t0, span, host, stem, rng):
    """
    A balanced, clocklike clade.

    Both properties are needed. BALANCE, because a caterpillar subdivides into
    (n-1) tips plus a singleton, so a tighter threshold collapses it instead of
    splitting it — the subdivision test would then measure the fixture, not the
    algorithm. CLOCKLIKE, because the rate-recovery test needs root-to-tip
    distance to equal elapsed time times the rate exactly, so internal branch
    lengths must be carved out of each tip's path rather than added to it.
    """
    tips = []
    for i in range(n):
        d = t0 + rng.random() * span
        tips.append((f"{prefix}{i:03d}|{host}|{d:.3f}",
                     max((d - t0) * RATE, 1e-6)))     # target root-to-tip depth

    def build(items, parent_depth):
        if len(items) == 1:
            name, depth = items[0]
            return f"{name}:{max(depth - parent_depth, 1e-7):.8f}"
        # Put this node halfway between its parent and its shallowest tip, so
        # every descendant keeps its exact target depth.
        node_depth = parent_depth + 0.5 * (min(d for _, d in items) - parent_depth)
        mid = len(items) // 2
        return (f"({build(items[:mid], node_depth)},"
                f"{build(items[mid:], node_depth)})"
                f":{max(node_depth - parent_depth, 1e-7):.8f}")

    return f"{build(tips, 0.0)[:-len(build(tips, 0.0).split(':')[-1]) - 1]}:{stem}"


def structured_tree(path, stems=(2.5, 2.7, 2.6), seed=7):
    rng = random.Random(seed)
    parts = [
        clade("A", 60, 1990, 30, "procyonid", stems[0], rng),
        clade("B", 45, 1995, 28, "wild_canid", stems[1], rng),
        clade("C", 35, 2000, 24, "domestic_dog", stems[2], rng),
    ]
    path.write_text("(" + ",".join(parts) + ");")
    return path


def run(tree, out, *extra):
    return subprocess.run(
        [sys.executable, "scripts/05c_cut_clades.py", "--tree", str(tree),
         "--out-tsv", str(out), "--permutations", "200", *map(str, extra)],
        capture_output=True, text=True, cwd=ROOT)


def rows(path):
    import csv
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_finds_the_clades_that_are_there(tmp_path):
    t = structured_tree(tmp_path / "t.nwk")
    r = run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 20)
    assert r.returncode == 0, r.stderr[-1500:]
    assert len(rows(tmp_path / "c.tsv")) == 3


def test_recovers_the_simulated_rate_within_each_clade(tmp_path):
    """The point of splitting: within a clade the slope IS the clock rate."""
    t = structured_tree(tmp_path / "t.nwk")
    run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 20)
    for row in rows(tmp_path / "c.tsv"):
        assert float(row["slope"]) == pytest.approx(RATE, rel=0.15)
        assert float(row["r_squared"]) > 0.8
        assert row["verdict"] in ("pass", "weak")


def test_small_groups_are_excluded(tmp_path):
    t = structured_tree(tmp_path / "t.nwk")
    run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 50)
    got = rows(tmp_path / "c.tsv")
    assert len(got) == 1                      # only the 60-tip clade qualifies
    assert int(got[0]["n_tips"]) == 60


def test_a_tighter_threshold_splits_further(tmp_path):
    """
    Tightening the threshold subdivides — but only if --min-size allows it.
    Below some point the subgroups get too small to qualify and the result is
    ZERO clades rather than more of them, which is the honest answer: nothing
    here is both tight enough and large enough. Tested with a small min-size so
    the subdivision is visible.
    """
    t = structured_tree(tmp_path / "t.nwk")
    run(t, tmp_path / "loose.tsv", "--max-diversity", 0.10, "--min-size", 3)
    run(t, tmp_path / "tight.tsv", "--max-diversity", 0.02, "--min-size", 3)
    assert len(rows(tmp_path / "tight.tsv")) > len(rows(tmp_path / "loose.tsv"))


def test_tightening_past_what_min_size_allows_yields_nothing(tmp_path):
    """Fewer clades, not smaller ones, once min_size can no longer be met."""
    t = structured_tree(tmp_path / "t.nwk")
    run(t, tmp_path / "c.tsv", "--max-diversity", 0.02, "--min-size", 30)
    assert rows(tmp_path / "c.tsv") == []


def test_clade_diversity_stays_under_the_threshold(tmp_path):
    t = structured_tree(tmp_path / "t.nwk")
    run(t, tmp_path / "c.tsv", "--max-diversity", 0.06, "--min-size", 20)
    for row in rows(tmp_path / "c.tsv"):
        assert float(row["diversity"]) <= 0.06


def test_tip_lists_are_written(tmp_path):
    t = structured_tree(tmp_path / "t.nwk")
    run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 20,
        "--out-dir", tmp_path / "tips")
    files = sorted((tmp_path / "tips").glob("*_tips.txt"))
    assert len(files) == 3
    assert sum(1 for _ in files[0].open()) >= 20


def test_references_label_clades_when_they_fall_inside(tmp_path):
    t = structured_tree(tmp_path / "t.nwk")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "# comment header\n"
        "A000\tAmerica-2\tstrain\tnuc\tsource\n"
        "B000\tEurope\tstrain\tnuc\tsource\n"
        "ZZZZ\tAsia-4\tstrain\tprot\tsource\n")   # protein: must be ignored
    r = run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 20,
            "--references", refs)
    assert "2 usable nucleotide references loaded" in r.stdout
    labels = {row["clade"]: row["reference_lineages"]
              for row in rows(tmp_path / "c.tsv")}
    assert "America-2" in "".join(labels.values())
    assert "Europe" in "".join(labels.values())
    assert "Asia-4" not in "".join(labels.values())


def test_coverage_filter_prunes_tips(tmp_path):
    """
    Ragged coverage makes root-to-tip measure sequence length rather than time.
    Filtering has to happen before the cut, not after.
    """
    pd = pytest.importorskip("pandas")
    t = structured_tree(tmp_path / "t.nwk")
    from Bio import Phylo
    names = [x.name for x in Phylo.read(str(t), "newick").get_terminals()]
    qc = tmp_path / "qc.tsv"
    pd.DataFrame({"name": names,
                  "ungapped_len": [1800 if n.startswith("A") else 500
                                   for n in names]}).to_csv(qc, sep="\t", index=False)
    r = run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 20,
            "--qc", qc, "--min-ungapped", 1500)
    assert "pruned 80 tips" in r.stdout
    got = rows(tmp_path / "c.tsv")
    assert len(got) == 1 and int(got[0]["n_tips"]) == 60


def test_no_clades_when_everything_is_too_divergent(tmp_path):
    """
    A saturated tree must yield nothing rather than a group that looks clocklike
    by accident. Returning a confident clade from saturated data is the failure
    this whole step exists to avoid.
    """
    t = structured_tree(tmp_path / "t.nwk")
    r = run(t, tmp_path / "c.tsv", "--max-diversity", 0.0001, "--min-size", 20)
    assert r.returncode == 0
    assert rows(tmp_path / "c.tsv") == []
    assert "No clade passed" in r.stdout


# --- bootstrap support ------------------------------------------------------

def supported_tree(path, value=98, seed=7):
    """The same structured tree with a support value on every internal node."""
    import re
    src = structured_tree(path.with_suffix(".plain"), seed=seed).read_text()
    path.write_text(re.sub(r"\)(:)", lambda m: f"){value}" + m.group(1), src))
    return path


def test_support_is_reported_when_the_tree_has_it(tmp_path):
    t = supported_tree(tmp_path / "s.nwk")
    r = run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 20)
    assert r.returncode == 0, r.stderr[-1000:]
    for row in rows(tmp_path / "c.tsv"):
        assert float(row["stem_support"]) == 98.0


def test_absent_support_is_not_treated_as_zero(tmp_path):
    """
    A tree built without bootstraps is not a tree with zero support. Conflating
    them would flag every clade on such a tree as unreliable and train the
    reader to ignore the column.
    """
    t = structured_tree(tmp_path / "p.nwk")
    r = run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 20,
            "--min-support", 95)
    assert r.returncode == 0
    assert "carries no support values" in r.stdout
    for row in rows(tmp_path / "c.tsv"):
        assert row["stem_support"] == ""


def test_weak_stem_support_is_called_out(tmp_path):
    """
    A clade cut on diversity alone can be an artefact of the tree search.
    Fitting a clock to a group that may not exist is worth a warning before
    days of BEAST go into it.
    """
    t = supported_tree(tmp_path / "w.nwk", value=42)
    r = run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 20,
            "--min-support", 95)
    assert r.returncode == 0
    assert "stem support below 95" in r.stdout
    assert "may not exist" in r.stdout


def test_strong_support_produces_no_warning(tmp_path):
    t = supported_tree(tmp_path / "g.nwk", value=99)
    r = run(t, tmp_path / "c.tsv", "--max-diversity", 0.10, "--min-size", 20,
            "--min-support", 95)
    assert "stem support below" not in r.stdout

"""
Tests for the date-randomisation gate.

Each scenario is built so the right answer is known from the construction
rather than from what the code emits. Two of these caught real bugs while the
gate was being written, and both were cases where a wrong answer looked
entirely reasonable:

  * the permissive criterion originally compared the real point estimate
    against the RANGE OF REPLICATE POINT ESTIMATES. That range shrinks toward
    the standard error as replicates are added, so with 20 replicates almost
    any real rate falls outside it -- and a genuine failure got promoted to
    "weak" by the weaker of the two criteria.
  * the independence check originally used a fixed threshold on the spread of
    replicate means. Independent replicates also have tightly clustered means
    (spread goes as SD/sqrt(ESS)), so the check measured sample size and fired
    on well-sampled runs.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location(
        "drt08d", ROOT / "scripts" / "08d_drt_summary.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def write_log(path, mean, sd, n=2000, rho=0.0, seed=0, param="clockRate"):
    rng = np.random.default_rng(seed)
    x = rng.normal(mean, sd, n)
    if rho:                                   # autocorrelated = poorly mixed
        for i in range(1, n):
            x[i] = rho * x[i - 1] + (1 - rho) * x[i]
    with open(path, "w") as fh:
        fh.write("# BEAST v2.7.7\n")
        fh.write(f"Sample\tposterior\t{param}\n")
        for i, v in enumerate(x):
            fh.write(f"{i * 1000}\t{-900 + rng.normal():.4f}\t{v:.6e}\n")


def scenario(tmp_path, name, real_mean, rep_mean, rep_sd=2e-5,
             rep_jitter=2.5e-5, n_reps=20, rho=0.0, identical=False):
    """
    Independent replicates get their OWN true rate: different permutations
    produce different likelihood surfaces. identical=True reuses one rate,
    which is what a repeated permutation or seed looks like.
    """
    rng = np.random.default_rng(99)
    real = tmp_path / f"{name}_real.log"
    write_log(real, real_mean, 3e-5, seed=1)
    reps = []
    for i in range(n_reps):
        p = tmp_path / f"{name}_rep{i}.log"
        m = rep_mean if identical else rng.normal(rep_mean, rep_jitter)
        write_log(p, m, rep_sd, rho=rho, seed=100 + i)
        reps.append(p)
    return real, reps


def run(mod, real, reps, burnin=0.10, min_ess=100):
    r = mod.rate_summary(real, 0.0, None)
    rr = [mod.rate_summary(p, burnin, r["param"]) for p in reps]
    return r, rr, mod.assess(r, rr, min_ess=min_ess)


# --- the four scenarios -----------------------------------------------------

def test_real_temporal_signal_passes(mod, tmp_path):
    real, reps = scenario(tmp_path, "A", 7.0e-4, 2.0e-4)
    _, _, res = run(mod, real, reps)
    assert res["verdict"] == "pass"
    assert res["n_overlapping"] == 0


def test_no_temporal_signal_fails(mod, tmp_path):
    """
    Real rate inside the replicate cloud, every HPD overlapping. This is the
    case the original permissive criterion rescued to "weak".
    """
    real, reps = scenario(tmp_path, "B", 2.1e-4, 2.0e-4, rep_sd=3e-5)
    _, _, res = run(mod, real, reps)
    assert res["verdict"] == "fail"
    assert res["n_overlapping"] == len(reps)


def test_unmixed_replicates_do_not_turn_a_pass_into_a_failure(mod, tmp_path):
    """
    Replicates run on short chains have wide posteriors, and a wide posterior
    overlaps everything. Poor mixing must not read as absent temporal signal.
    """
    real, reps = scenario(tmp_path, "C", 7.0e-4, 3.0e-4, rep_sd=2.5e-4, rho=0.995)
    _, rr, res = run(mod, real, reps)
    assert res["verdict"] == "pass"
    assert res["low_ess_replicates"], "the low ESS should still be reported"
    assert any("ESS" in r for r in res["reasons"])


def test_reused_permutation_is_caught(mod, tmp_path):
    """
    Every replicate identical: one permutation counted twenty times. The test
    would pass trivially, so it must not be allowed to read as a pass.
    """
    real, reps = scenario(tmp_path, "D", 7.0e-4, 2.0e-4, identical=True)
    _, _, res = run(mod, real, reps)
    assert res["verdict"] == "inconclusive"
    assert any("permutation" in r for r in res["reasons"])


def test_overlap_in_the_bad_direction_is_still_caught(mod, tmp_path):
    """A real rate BELOW the randomised ones is as non-identifiable as above."""
    real, reps = scenario(tmp_path, "E", 5.0e-5, 2.0e-4)
    _, _, res = run(mod, real, reps)
    assert res["verdict"] == "pass"       # non-overlap is non-overlap
    assert res["real_outside_replicate_range"]


# --- the two regressions ----------------------------------------------------

def test_permissive_criterion_uses_the_envelope_not_the_range_of_means(mod, tmp_path):
    """
    Regression. The range of replicate MEANS shrinks toward the standard error
    as replicates are added; the envelope of their HPDs does not. Testing
    against the range made nearly every real rate look "outside" and rescued
    genuine failures.
    """
    real, reps = scenario(tmp_path, "F", 2.1e-4, 2.0e-4, rep_sd=3e-5)
    r, rr, res = run(mod, real, reps)
    rep_mean_range = (min(x["mean"] for x in rr), max(x["mean"] for x in rr))
    envelope = res["replicate_envelope"]
    assert envelope[0] < rep_mean_range[0] and envelope[1] > rep_mean_range[1]
    assert envelope[0] <= r["mean"] <= envelope[1]
    assert not res["real_outside_replicate_range"]


def test_independence_check_does_not_fire_on_well_sampled_replicates(mod, tmp_path):
    """
    Regression. Independent replicates have tightly clustered means too, so a
    fixed threshold on their spread measures sample size rather than
    independence. The comparison is against the Monte Carlo error of a single
    replicate's mean.
    """
    real, reps = scenario(tmp_path, "G", 7.0e-4, 2.0e-4, n_reps=20)
    _, _, res = run(mod, real, reps)
    assert res["replicate_mean_spread"] > res["replicate_mc_error"]
    assert not any("permutation" in r for r in res["reasons"])


# --- clock parameter detection ---------------------------------------------

def test_clock_parameter_is_detected_from_the_usual_names():
    from lib.beastlog import find_clock_parameter
    assert find_clock_parameter(["Sample", "posterior", "clockRate"]) == "clockRate"
    assert find_clock_parameter(["Sample", "ucldMean", "ucldStdev"]) == "ucldMean"
    assert find_clock_parameter(
        ["Sample", "clockRate.c:america2"]) == "clockRate.c:america2"


def test_ucld_stdev_is_never_mistaken_for_the_rate():
    """
    ucldStdev is the clock's variation, not its rate. Comparing it against
    randomised dates returns a confident verdict about the wrong parameter.
    """
    from lib.beastlog import find_clock_parameter
    assert find_clock_parameter(
        ["Sample", "ucldStdev", "ucldMean"]) == "ucldMean"


def test_a_log_with_no_clock_column_is_refused():
    from lib.beastlog import find_clock_parameter
    with pytest.raises(ValueError, match="no clock-rate column"):
        find_clock_parameter(["Sample", "posterior", "treeLikelihood"])


def test_explicit_clock_param_overrides_detection():
    from lib.beastlog import find_clock_parameter
    header = ["Sample", "clockRate", "myRate"]
    assert find_clock_parameter(header, preferred="myRate") == "myRate"


# --- CLI --------------------------------------------------------------------

def test_cli_writes_json_and_exits_nonzero_on_failure(tmp_path):
    rng = np.random.default_rng(7)
    real = tmp_path / "real.log"
    write_log(real, 2.1e-4, 3e-5, seed=1)
    reps = []
    for i in range(10):
        p = tmp_path / f"rep{i}.log"
        write_log(p, rng.normal(2.0e-4, 2.5e-5), 3e-5, seed=100 + i)
        reps.append(str(p))
    out = tmp_path / "drt.json"
    r = subprocess.run(
        [sys.executable, "scripts/08d_drt_summary.py", "--real", str(real),
         "--reps", *reps, "--out-json", str(out), "--exit-nonzero-on-fail"],
        cwd=ROOT, capture_output=True, text=True, timeout=180)
    assert r.returncode == 1, r.stdout[-1500:]
    doc = json.loads(out.read_text())
    assert doc["verdict"] == "fail"
    assert doc["real"]["param"] == "clockRate"
    assert len(doc["replicates"]) == 10


def test_cli_warns_when_there_are_too_few_replicates(tmp_path):
    real = tmp_path / "real.log"
    write_log(real, 7.0e-4, 3e-5, seed=1)
    reps = []
    for i in range(3):
        p = tmp_path / f"rep{i}.log"
        write_log(p, 2.0e-4, 2e-5, seed=200 + i)
        reps.append(str(p))
    r = subprocess.run(
        [sys.executable, "scripts/08d_drt_summary.py", "--real", str(real),
         "--reps", *reps],
        cwd=ROOT, capture_output=True, text=True, timeout=180)
    assert "only 3 replicates" in r.stderr


def test_cli_help_works():
    r = subprocess.run([sys.executable, "scripts/08d_drt_summary.py", "--help"],
                       cwd=ROOT, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0

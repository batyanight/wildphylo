"""
beastlog.py — reading BEAST .log files, and the statistics computed from them.

Shared by 08b_convergence (did the chains converge and agree?) and
08d_drt_summary (does the real clock rate sit outside the date-randomised
ones?). Both need the same parser, the same HPD and the same ESS, and a
convergence verdict computed one way against a DRT verdict computed another
would be two statements about one chain that cannot both be checked.

It lives here rather than inside 08b for the same reason traittrees.py does:
a module whose name starts with a digit cannot be imported normally, so shared
code in a numbered script has to be reached by importlib and breaks whenever
the pipeline numbering changes. The numbers are a reading order for humans.
"""

from __future__ import annotations

import numpy as np

__all__ = ["MAX_LAG", "NON_PARAMETERS", "read_beast_log", "ess", "hpd",
           "overlap", "find_clock_parameter", "CLOCK_CANDIDATES"]

# Tracer stops accumulating autocovariance at this lag. Kept identical to
# BEAST's ESS.java so the ESS reported here is the ESS the user sees in Tracer.
MAX_LAG = 2000

# Columns that are bookkeeping rather than parameters. An ESS for the sample
# index is meaningless, and posterior/likelihood/prior are diagnostics whose
# ESS is genuinely informative, so those stay.
NON_PARAMETERS = {"Sample", "state", "STATE"}


def read_beast_log(path: Path) -> tuple[list[str], np.ndarray]:
    """
    Parse a BEAST .log: '#' comment lines, then a tab-separated header whose
    first column is Sample, then numeric rows.

    A truncated final line is dropped rather than raising: a log from a chain
    still running, or killed mid-write, is a normal thing to look at, and
    failing on it would mean the only way to check progress is to wait.
    """
    header: list[str] | None = None
    rows: list[list[float]] = []
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n")
            if not s.strip() or s.lstrip().startswith("#"):
                continue
            parts = s.split("\t")
            if header is None:
                header = [p.strip() for p in parts]
                continue
            if len(parts) != len(header):
                continue                      # truncated or interleaved line
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue                      # a non-numeric column; skip row
    if header is None or not rows:
        raise ValueError(f"{path}: no usable samples found")
    return header, np.asarray(rows, dtype=np.float64)


def ess(x: np.ndarray) -> float:
    """
    Effective sample size, by Tracer's algorithm.

    Tracer accumulates the autocovariance gamma[lag], adding lag pairs to the
    variance estimate until a consecutive pair sums to a negative value, which
    is Geyer's initial-positive-sequence stopping rule. ESS is then
    N * gamma[0] / varStat.

    The autocovariance itself is computed by FFT rather than by Tracer's nested
    loop. The values are identical to floating point; the difference is that a
    10,000-sample log with thirty parameters takes under a second instead of
    several minutes, which decides whether this runs every build or gets
    skipped.
    """
    n = len(x)
    if n < 10:
        return float(n)
    v = x - x.mean()
    if not np.any(v):
        return float(n)                        # constant: every sample identical

    max_lag = min(n, MAX_LAG)
    # Linear (not circular) autocovariance via zero-padded FFT.
    size = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(v, size)
    acov = np.fft.irfft(f * np.conjugate(f), size)[:max_lag]
    acov /= (n - np.arange(max_lag))            # Tracer divides by n - lag

    var_stat = acov[0]
    if var_stat == 0:
        return float(n)
    for lag in range(2, max_lag, 2):
        pair = acov[lag - 1] + acov[lag]
        if pair <= 0:
            break
        var_stat += 2.0 * pair
    if var_stat <= 0:
        return float(n)
    return float(min(n, n * acov[0] / var_stat))


def hpd(x: np.ndarray, level: float = 0.95) -> tuple[float, float]:
    """Highest posterior density interval: the shortest interval holding `level`."""
    s = np.sort(x)
    n = len(s)
    k = max(1, int(np.floor(level * n)))
    if k >= n:
        return float(s[0]), float(s[-1])
    widths = s[k:] - s[:n - k]
    i = int(np.argmin(widths))
    return float(s[i]), float(s[i + k])


def overlap(a: np.ndarray, b: np.ndarray, bins: int = 100) -> float:
    """
    Overlapping coefficient of two samples: sum of min(p, q) over a shared
    histogram grid. 1.0 identical, 0.0 disjoint.

    The grid spans both samples, so two chains that explored different regions
    get a low value by construction rather than by a tuning choice.
    """
    lo = min(a.min(), b.min())
    hi = max(a.max(), b.max())
    if hi == lo:
        return 1.0
    edges = np.linspace(lo, hi, bins + 1)
    pa, _ = np.histogram(a, bins=edges, density=False)
    pb, _ = np.histogram(b, bins=edges, density=False)
    pa = pa / pa.sum()
    pb = pb / pb.sum()
    return float(np.minimum(pa, pb).sum())


# Names BEAST gives the clock rate, in the order worth trying. Which one
# appears depends on the clock model: a strict clock logs clockRate, an
# uncorrelated lognormal logs ucldMean (and ucldStdev, which is NOT the rate).
# Guessing wrong here means the date-randomisation test compares the wrong
# parameter and returns a confident verdict about something else.
CLOCK_CANDIDATES = (
    "clockRate", "clock.rate", "ucldMean", "ucedMean", "rate.mean",
    "meanRate", "rate", "clockRate.c", "ucldMean.c",
)


def find_clock_parameter(header, preferred=None):
    """
    Pick the column holding the clock rate.

    Matching is exact first, then by prefix, because BEAST appends the
    partition name (clockRate.c:america2). Explicitly refuses ucldStdev and
    the .prob/.indicator columns, which are the wrong parameter and would
    otherwise match a loose prefix rule.
    """
    if preferred:
        if preferred in header:
            return preferred
        hits = [h for h in header if h.startswith(preferred)]
        if len(hits) == 1:
            return hits[0]
        raise ValueError(
            f"--clock-param {preferred!r} matches {len(hits)} columns: {hits}")
    usable = [h for h in header
              if not h.endswith(("Stdev", "stdev", ".prob", ".indicator"))]
    for cand in CLOCK_CANDIDATES:
        if cand in usable:
            return cand
    for cand in CLOCK_CANDIDATES:
        hits = [h for h in usable if h.split(".")[0] == cand
                or h.split(":")[0].rstrip(".c") == cand]
        if hits:
            return sorted(hits, key=len)[0]
    raise ValueError(
        "no clock-rate column found. Tried " + ", ".join(CLOCK_CANDIDATES)
        + f"; the log has {sorted(usable)[:12]}. Pass --clock-param.")

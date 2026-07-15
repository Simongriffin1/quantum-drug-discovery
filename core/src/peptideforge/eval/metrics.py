"""Core correlation / error metrics for oracle and surrogate validity."""

from __future__ import annotations

import math
from collections.abc import Sequence


def _rankdata(values: Sequence[float]) -> list[float]:
    """Average ranks for ties (1-based ranks)."""
    indexed = sorted(enumerate(values), key=lambda iv: iv[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        # average rank of positions i..j (1-based: i+1 .. j+1)
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg
        i = j + 1
    return ranks


def pearson_r(x: Sequence[float], y: Sequence[float]) -> float:
    """Pearson correlation; raises if n < 2 or zero variance."""
    if len(x) != len(y):
        raise ValueError("x and y must have equal length")
    n = len(x)
    if n < 2:
        raise ValueError("need at least 2 pairs for correlation")
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    num = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y, strict=True))
    den_x = math.sqrt(sum((a - mean_x) ** 2 for a in x))
    den_y = math.sqrt(sum((b - mean_y) ** 2 for b in y))
    if den_x == 0.0 or den_y == 0.0:
        raise ValueError("zero variance — correlation undefined")
    return num / (den_x * den_y)


def spearman_rho(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rank correlation (tie-aware via average ranks)."""
    return pearson_r(_rankdata(x), _rankdata(y))


def rmse(predicted: Sequence[float], experimental: Sequence[float]) -> float:
    """Root mean squared error."""
    if len(predicted) != len(experimental):
        raise ValueError("predicted and experimental must have equal length")
    if not predicted:
        raise ValueError("need at least one pair for RMSE")
    mse = sum((p - e) ** 2 for p, e in zip(predicted, experimental, strict=True)) / len(predicted)
    return math.sqrt(mse)


def bootstrap_ci(
    x: Sequence[float],
    y: Sequence[float],
    *,
    statistic: str = "spearman",
    n_resamples: int = 1000,
    seed: int = 0,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    """Bootstrap percentile CI for Spearman or Pearson.

    Returns (point_estimate, ci_low, ci_high). A point estimate alone is
    insufficient to declare PASS/FAIL — always report the CI.
    """
    import random

    if len(x) != len(y):
        raise ValueError("x and y must have equal length")
    if len(x) < 3:
        raise ValueError("bootstrap_ci needs ≥ 3 pairs")
    if n_resamples < 100:
        raise ValueError("n_resamples must be ≥ 100 for a usable CI")
    if statistic not in {"spearman", "pearson"}:
        raise ValueError("statistic must be 'spearman' or 'pearson'")
    fn = spearman_rho if statistic == "spearman" else pearson_r
    point = fn(x, y)
    rng = random.Random(seed)
    n = len(x)
    samples: list[float] = []
    for _ in range(n_resamples):
        idx = [rng.randrange(n) for _ in range(n)]
        xs = [x[i] for i in idx]
        ys = [y[i] for i in idx]
        try:
            samples.append(fn(xs, ys))
        except ValueError:
            continue
    if len(samples) < max(100, n_resamples // 2):
        raise RuntimeError(
            f"bootstrap produced only {len(samples)} valid resamples "
            f"(need ≥{max(100, n_resamples // 2)}); check variance / ties"
        )
    samples.sort()
    lo_i = int(math.floor(alpha / 2 * (len(samples) - 1)))
    hi_i = int(math.ceil((1 - alpha / 2) * (len(samples) - 1)))
    return point, samples[lo_i], samples[hi_i]

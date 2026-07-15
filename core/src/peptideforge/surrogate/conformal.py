"""Split conformal prediction for calibrated prediction intervals.

Method: inductive (split) conformal (Vovk / Lei / Tibshirani). Fit on a training
partition, compute absolute residuals on a held-out calibration set, and take the
empirical (1−α) quantile as the interval half-width around the point prediction.
"""

from __future__ import annotations

import math
from collections.abc import Sequence


def conformal_quantile(scores: Sequence[float], *, coverage: float) -> float:
    """Finite-sample conformal quantile of absolute residuals.

    Uses the floor((n+1)·q)/n order-statistic convention (conservative).
    """
    if not scores:
        raise ValueError("conformal_quantile requires ≥ 1 score")
    if not 0.0 < coverage < 1.0:
        raise ValueError("coverage must be in (0, 1)")
    sorted_scores = sorted(scores)
    n = len(sorted_scores)
    # level q = ceil((n+1)·coverage) / n  → index = ceil((n+1)·coverage) - 1
    level = math.ceil((n + 1) * coverage) / n
    idx = min(n - 1, max(0, math.ceil(level * n) - 1))
    return sorted_scores[idx]


def absolute_residuals(predicted: Sequence[float], observed: Sequence[float]) -> list[float]:
    if len(predicted) != len(observed):
        raise ValueError("predicted and observed length mismatch")
    return [abs(p - o) for p, o in zip(predicted, observed, strict=True)]

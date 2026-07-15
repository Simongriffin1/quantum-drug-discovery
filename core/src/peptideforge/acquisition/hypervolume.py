"""Hypervolume utilities for multi-objective acquisition (P8).

2-objective exact HV is O(n log n); for M>2 we use a Monte Carlo estimate.
Objectives are treated as **maximization**; flip signs for minimize-oriented axes
before calling. Reference point must be dominated by attainable objectives.
"""

from __future__ import annotations

import random
from collections.abc import Sequence


def dominates(a: Sequence[float], b: Sequence[float]) -> bool:
    """True if a weakly dominates b and is strictly better on ≥1 objective (maximize)."""
    if len(a) != len(b):
        raise ValueError("objective vectors must have equal length")
    ge = all(x >= y for x, y in zip(a, b, strict=True))
    gt = any(x > y for x, y in zip(a, b, strict=True))
    return ge and gt


def nondominated_front(points: Sequence[Sequence[float]]) -> list[tuple[float, ...]]:
    """Pareto front under maximization."""
    front: list[tuple[float, ...]] = []
    for p in points:
        pt = tuple(float(x) for x in p)
        if any(dominates(f, pt) for f in front):
            continue
        front = [f for f in front if not dominates(pt, f)]
        front.append(pt)
    return front


def hypervolume_2d(front: Sequence[Sequence[float]], reference: Sequence[float]) -> float:
    """Exact 2-D hypervolume (maximization) vs a dominated reference point."""
    if len(reference) != 2:
        raise ValueError("hypervolume_2d requires a 2-D reference")
    pts = [
        (p[0], p[1])
        for p in nondominated_front(front)
        if p[0] > reference[0] and p[1] > reference[1]
    ]
    if not pts:
        return 0.0
    pts.sort(key=lambda p: p[0])
    # Compress to monotonically decreasing staircase (ascending x ⇒ descending y)
    filtered: list[tuple[float, float]] = []
    for x, y in pts:
        while filtered and filtered[-1][1] <= y:
            filtered.pop()
        filtered.append((x, y))
    hv = 0.0
    prev_x = reference[0]
    for x, y in filtered:
        hv += (x - prev_x) * (y - reference[1])
        prev_x = x
    return max(0.0, hv)


def hypervolume(
    front: Sequence[Sequence[float]],
    reference: Sequence[float],
    *,
    seed: int | None = None,
    n_samples: int = 8_000,
) -> float:
    """Hypervolume vs reference (maximize). Exact for 2-D; MC otherwise."""
    pts = nondominated_front(front)
    if not pts:
        return 0.0
    m = len(reference)
    if any(len(p) != m for p in pts):
        raise ValueError("front points must match reference dimension")
    if m == 2:
        return hypervolume_2d(pts, reference)
    rng = random.Random(seed)
    highs = [max(p[j] for p in pts) for j in range(m)]
    volume = 1.0
    for j in range(m):
        side = highs[j] - reference[j]
        if side <= 0.0:
            return 0.0
        volume *= side
    hits = 0
    for _ in range(n_samples):
        sample = [reference[j] + rng.random() * (highs[j] - reference[j]) for j in range(m)]
        if any(all(pj >= sj for pj, sj in zip(p, sample, strict=True)) for p in pts):
            hits += 1
    return volume * (hits / n_samples)


def hypervolume_improvement(
    front: Sequence[Sequence[float]],
    candidate: Sequence[float],
    reference: Sequence[float],
) -> float:
    """HV(front ∪ {candidate}) − HV(front)."""
    base = hypervolume(front, reference)
    new = hypervolume([*front, candidate], reference)
    return max(0.0, new - base)

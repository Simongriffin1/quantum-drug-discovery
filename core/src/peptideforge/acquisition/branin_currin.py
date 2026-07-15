"""Branin–Currin synthetic multi-objective benchmark (P8 acquisition validation).

Classic two-objective test functions used to validate BO acquisition independent
of biology. Objectives are **maximized** as (−Branin, −Currin). All fixtures
are ``synthetic_*`` — not physics labels.
"""

from __future__ import annotations

import math
import random
from uuid import UUID

from peptideforge.contracts.models import (
    CalibratedPrediction,
    Candidates,
    ObjectiveVector,
    PeptideCandidate,
)


def branin(x1: float, x2: float) -> float:
    """Branin function (minimized form); domain x1∈[-5,10], x2∈[0,15]."""
    a = 1.0
    b = 5.1 / (4.0 * math.pi**2)
    c = 5.0 / math.pi
    r = 6.0
    s = 10.0
    t = 1.0 / (8.0 * math.pi)
    return a * (x2 - b * x1**2 + c * x1 - r) ** 2 + s * (1.0 - t) * math.cos(x1) + s


def currin(x1: float, x2: float) -> float:
    """Currin exponential function (minimized form); domain x∈[0,1]²."""
    # Avoid division by zero at x2=0
    x2_safe = max(x2, 1e-6)
    term = 1.0 - math.exp(-1.0 / (2.0 * x2_safe))
    num = 2300 * x1**3 + 1900 * x1**2 + 2092 * x1 + 60
    den = 100 * x1**3 + 500 * x1**2 + 4 * x1 + 20
    return term * num / den


def branin_currin_objectives(x1: float, x2: float) -> tuple[float, float]:
    """Maximization objectives: (−Branin, −Currin) after mapping to Currin domain."""
    # Map Branin domain coords to [0,1] for Currin
    u = (x1 + 5.0) / 15.0
    v = x2 / 15.0
    u = min(1.0, max(0.0, u))
    v = min(1.0, max(0.0, v))
    return (-branin(x1, x2), -currin(u, v))


def make_branin_currin_pool(
    n: int,
    *,
    seed: int = 0,
) -> tuple[Candidates, dict[UUID, tuple[float, float]], dict[UUID, tuple[float, float]]]:
    """Discrete synthetic pool on the Branin domain.

    Returns (candidates, coords_by_id, objectives_by_id).
    Sequences are placeholder synthetic_* peptides (length ≥5).
    """
    rng = random.Random(seed)
    items: list[PeptideCandidate] = []
    coords: dict[UUID, tuple[float, float]] = {}
    objs: dict[UUID, tuple[float, float]] = {}
    aa = "ACDEFGHIKLMNPQRSTVWY"
    for i in range(n):
        x1 = rng.uniform(-5.0, 10.0)
        x2 = rng.uniform(0.0, 15.0)
        # Unique synthetic sequence tag
        tag = "".join(aa[(i * 7 + k) % 20] for k in range(8))
        cand = PeptideCandidate(
            sequence=tag,
            generation_method="synthetic_branin_currin",
            metadata={"x1": x1, "x2": x2, "idx": i},
        )
        items.append(cand)
        coords[cand.candidate_id] = (x1, x2)
        objs[cand.candidate_id] = branin_currin_objectives(x1, x2)
    return Candidates(items=tuple(items), seed=seed), coords, objs


def perfect_surrogate_predictions(
    pool: Candidates,
    objectives: dict[UUID, tuple[float, float]],
    *,
    epistemic_std: float = 0.5,
    observed: dict[UUID, tuple[float, float]] | None = None,
    observed_noise: float = 0.0,
) -> tuple[ObjectiveVector, ...]:
    """Near-oracle predictions for validating acquisition logic (synthetic_*).

    Unobserved points get true means + fixed epistemic_std.
    Observed points get noisy labels as means with reduced std.
    """
    out: list[ObjectiveVector] = []
    for cand in pool.items:
        true = objectives[cand.candidate_id]
        if observed is not None and cand.candidate_id in observed:
            mean = observed[cand.candidate_id]
            std = max(1e-3, observed_noise)
        else:
            mean = true
            std = epistemic_std
        preds = (
            CalibratedPrediction(
                candidate_id=cand.candidate_id,
                objective_name="neg_branin",
                mean=mean[0],
                lower=mean[0] - 2 * std,
                upper=mean[0] + 2 * std,
                epistemic_std=std,
                coverage_target=0.9,
            ),
            CalibratedPrediction(
                candidate_id=cand.candidate_id,
                objective_name="neg_currin",
                mean=mean[1],
                lower=mean[1] - 2 * std,
                upper=mean[1] + 2 * std,
                epistemic_std=std,
                coverage_target=0.9,
            ),
        )
        out.append(ObjectiveVector(candidate_id=cand.candidate_id, predictions=preds))
    return tuple(out)

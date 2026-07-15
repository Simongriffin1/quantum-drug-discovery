"""Random acquisition baseline — uniform over eligible candidates."""

from __future__ import annotations

import random

from peptideforge.contracts.models import AcquisitionBatch, Candidates, ObjectiveVector
from peptideforge.acquisition.utils import (
    align_predictions,
    apply_constraints,
    cost_per_candidate,
    pack_batch,
)


class RandomAcquisition:
    """Uniform random batch from the constraint-filtered pool (baseline)."""

    method = "random"

    def rank(
        self,
        pool: Candidates,
        predictions: tuple[ObjectiveVector, ...],
        *,
        batch_size: int,
        budget_remaining: float,
        constraints: dict[str, object] | None = None,
        seed: int | None = None,
    ) -> AcquisitionBatch:
        if batch_size < 1:
            raise ValueError("batch_size must be ≥ 1")
        if budget_remaining < 0.0:
            raise ValueError("budget_remaining must be ≥ 0")
        pred_map = align_predictions(pool, predictions)
        eligible, constrained = apply_constraints(pool, pred_map, constraints)
        cost = cost_per_candidate(constraints)
        max_affordable = int(budget_remaining // cost) if cost > 0 else 0
        n_select = min(batch_size, len(eligible), max_affordable)

        rng = random.Random(seed)
        chosen = rng.sample(eligible, n_select) if n_select else []
        selected = [(cid, 0.0) for cid in chosen]
        remaining = budget_remaining - n_select * cost
        ranked = pack_batch(
            selected=selected,
            constrained=constrained,
            batch_size=batch_size,
            budget_remaining=remaining,
            method=self.method,
            seed=seed,
        )
        return AcquisitionBatch(
            ranked=ranked,
            batch_size=batch_size,
            budget_remaining=remaining,
            method=self.method,
            seed=seed,
        )

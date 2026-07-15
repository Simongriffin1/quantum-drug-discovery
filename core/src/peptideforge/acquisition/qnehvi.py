"""qNEHVI-style batch acquisition for discrete candidate pools (P8).

Method: Monte Carlo estimate of Noisy Expected Hypervolume Improvement over a
discrete pool (Daulton et al. / BoTorch qNEHVI idea). Posterior samples are drawn
from N(mean, epistemic_std²) per objective (independent axes — a pragmatic
approximation when full covariances are unavailable).

Optional BoTorch path: set ``use_botorch=True``; raises if BoTorch is missing.
"""

from __future__ import annotations

import random
from uuid import UUID

from peptideforge.acquisition.hypervolume import hypervolume, nondominated_front
from peptideforge.acquisition.utils import (
    align_predictions,
    apply_constraints,
    cost_per_candidate,
    objective_means,
    objective_stds,
    pack_batch,
)
from peptideforge.contracts.models import AcquisitionBatch, Candidates, ObjectiveVector


class BoTorchUnavailableError(ImportError):
    """Raised when BoTorch qNEHVI is requested but not installed."""


class QNEHVIAcquisition:
    """Greedy Monte Carlo qNEHVI over a discrete candidate pool."""

    method = "qnehvi_mc"

    def __init__(
        self,
        *,
        reference_point: tuple[float, ...] | None = None,
        n_mc: int = 64,
        pareto_front: tuple[tuple[float, ...], ...] = (),
        use_botorch: bool = False,
    ) -> None:
        if n_mc < 8:
            raise ValueError("n_mc must be ≥ 8")
        self.reference_point = reference_point
        self.n_mc = n_mc
        self.pareto_front = pareto_front
        self.use_botorch = use_botorch

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
        if self.use_botorch:
            raise BoTorchUnavailableError(
                "BoTorch qNEHVI wiring is optional and not bundled in CI. "
                "Install botorch+torch and use a BoTorch-backed backend, or "
                "set use_botorch=False to use the stdlib MC qNEHVI."
            )

        pred_map = align_predictions(pool, predictions)
        eligible, constrained = apply_constraints(pool, pred_map, constraints)
        if not eligible:
            return AcquisitionBatch(
                ranked=pack_batch(
                    selected=[],
                    constrained=constrained,
                    batch_size=batch_size,
                    budget_remaining=budget_remaining,
                    method=self.method,
                    seed=seed,
                ),
                batch_size=batch_size,
                budget_remaining=budget_remaining,
                method=self.method,
                seed=seed,
            )

        # Infer objective dim + reference
        sample_vec = pred_map[eligible[0]]
        n_obj = len(sample_vec.predictions)
        ref = self.reference_point
        if ref is None:
            # Dominated corner: min mean − 1 on each axis
            means = [objective_means(pred_map[cid]) for cid in eligible]
            ref = tuple(min(m[j] for m in means) - 1.0 for j in range(n_obj))
        if len(ref) != n_obj:
            raise ValueError(
                f"reference_point dim {len(ref)} != n_objectives {n_obj}"
            )

        cost = cost_per_candidate(constraints)
        max_affordable = int(budget_remaining // cost)
        n_select = min(batch_size, len(eligible), max_affordable)

        rng = random.Random(seed)
        front: list[tuple[float, ...]] = list(self.pareto_front)
        remaining_ids = list(eligible)
        selected: list[tuple[UUID, float]] = []

        for _ in range(n_select):
            best_id: UUID | None = None
            best_score = float("-inf")
            for cid in remaining_ids:
                score = self._expected_hvi(pred_map[cid], front, ref, rng)
                if score > best_score:
                    best_score = score
                    best_id = cid
            if best_id is None:
                break
            selected.append((best_id, best_score))
            remaining_ids.remove(best_id)
            # Fantasize: add mean to front for sequential greedy q
            front = nondominated_front([*front, objective_means(pred_map[best_id])])

        remaining_budget = budget_remaining - len(selected) * cost
        ranked = pack_batch(
            selected=selected,
            constrained=constrained,
            batch_size=batch_size,
            budget_remaining=remaining_budget,
            method=self.method,
            seed=seed,
        )
        return AcquisitionBatch(
            ranked=ranked,
            batch_size=batch_size,
            budget_remaining=remaining_budget,
            method=self.method,
            seed=seed,
        )

    def _expected_hvi(
        self,
        vec: ObjectiveVector,
        front: list[tuple[float, ...]],
        reference: tuple[float, ...],
        rng: random.Random,
    ) -> float:
        means = objective_means(vec)
        stds = objective_stds(vec)
        base = hypervolume(front, reference)
        total = 0.0
        for _ in range(self.n_mc):
            sample = tuple(
                rng.gauss(m, s) for m, s in zip(means, stds, strict=True)
            )
            hv = hypervolume([*front, sample], reference)
            total += max(0.0, hv - base)
        return total / self.n_mc

"""Validate qNEHVI vs random on Branin–Currin by hypervolume (P8 acceptance)."""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict

from peptideforge.acquisition.branin_currin import (
    make_branin_currin_pool,
    perfect_surrogate_predictions,
)
from peptideforge.acquisition.hypervolume import hypervolume, nondominated_front
from peptideforge.acquisition.qnehvi import QNEHVIAcquisition
from peptideforge.acquisition.random_acq import RandomAcquisition
from peptideforge.contracts.models import Candidates
from uuid import UUID


class AcquisitionValidationReport(BaseModel):
    """Honest HV comparison — acquisition logic only, not biology."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    n_pool: int
    n_init: int
    n_rounds: int
    batch_size: int
    hv_qnehvi: float
    hv_random: float
    reference: tuple[float, float]
    passed: bool
    notes: str | None = None


def run_branin_currin_validation(
    *,
    n_pool: int = 80,
    n_init: int = 8,
    n_rounds: int = 6,
    batch_size: int = 2,
    seed: int = 0,
    epistemic_std: float = 1.0,
) -> AcquisitionValidationReport:
    """Sequential batch BO: qNEHVI must finish with higher HV than random."""
    pool, _coords, objectives = make_branin_currin_pool(n_pool, seed=seed)
    # Reference: slightly worse than empirical nadir of the pool
    all_obj = list(objectives.values())
    ref = (
        min(o[0] for o in all_obj) - 1.0,
        min(o[1] for o in all_obj) - 1.0,
    )

    hv_q = _campaign(
        pool,
        objectives,
        acquisition=QNEHVIAcquisition(
            reference_point=ref,
            n_mc=48,
            pareto_front=(),
        ),
        n_init=n_init,
        n_rounds=n_rounds,
        batch_size=batch_size,
        seed=seed,
        reference=ref,
        epistemic_std=epistemic_std,
    )
    hv_r = _campaign(
        pool,
        objectives,
        acquisition=RandomAcquisition(),
        n_init=n_init,
        n_rounds=n_rounds,
        batch_size=batch_size,
        seed=seed + 1,
        reference=ref,
        epistemic_std=epistemic_std,
    )
    return AcquisitionValidationReport(
        n_pool=n_pool,
        n_init=n_init,
        n_rounds=n_rounds,
        batch_size=batch_size,
        hv_qnehvi=hv_q,
        hv_random=hv_r,
        reference=ref,
        passed=hv_q > hv_r,
        notes="synthetic_branin_currin_perfect_surrogate",
    )


def _campaign(
    pool: Candidates,
    objectives: dict[UUID, tuple[float, float]],
    *,
    acquisition: QNEHVIAcquisition | RandomAcquisition,
    n_init: int,
    n_rounds: int,
    batch_size: int,
    seed: int,
    reference: tuple[float, float],
    epistemic_std: float,
) -> float:
    import random

    rng = random.Random(seed)
    ids = [c.candidate_id for c in pool.items]
    labeled: set[UUID] = set(rng.sample(ids, n_init))
    for _round in range(n_rounds):
        remaining = Candidates(
            items=tuple(c for c in pool.items if c.candidate_id not in labeled),
            seed=seed,
        )
        if not remaining.items:
            break
        observed = {cid: objectives[cid] for cid in labeled}
        front = nondominated_front(list(observed.values()))
        preds = perfect_surrogate_predictions(
            remaining,
            objectives,
            epistemic_std=epistemic_std,
            observed=None,
        )
        acq: QNEHVIAcquisition | RandomAcquisition
        if isinstance(acquisition, QNEHVIAcquisition):
            acq = QNEHVIAcquisition(
                reference_point=reference,
                n_mc=acquisition.n_mc,
                pareto_front=tuple(front),
            )
        else:
            acq = acquisition
        budget = float(batch_size)
        batch = acq.rank(
            remaining,
            preds,
            batch_size=batch_size,
            budget_remaining=budget,
            seed=seed + _round,
        )
        picked = [r.candidate_id for r in batch.ranked if not r.constrained_out]
        labeled.update(picked)

    front = nondominated_front([objectives[cid] for cid in labeled])
    return hypervolume(front, reference)

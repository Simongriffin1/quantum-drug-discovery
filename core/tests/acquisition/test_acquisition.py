"""Tests for P8 acquisition (qNEHVI + Branin–Currin validation)."""

from __future__ import annotations

import pytest

from peptideforge.acquisition import (
    BoTorchUnavailableError,
    QNEHVIAcquisition,
    RandomAcquisition,
    hypervolume,
    make_branin_currin_pool,
    nondominated_front,
    perfect_surrogate_predictions,
    run_branin_currin_validation,
)
from peptideforge.contracts.models import (
    CalibratedPrediction,
    Candidates,
    ObjectiveVector,
    PeptideCandidate,
)
from peptideforge.contracts.protocols import AcquisitionFunction


def _toy_pool() -> tuple[Candidates, tuple[ObjectiveVector, ...]]:
    cands = [
        PeptideCandidate(sequence="AAAAA", generation_method="synthetic_acq"),
        PeptideCandidate(sequence="CCCCC", generation_method="synthetic_acq"),
        PeptideCandidate(sequence="DDDDD", generation_method="synthetic_acq"),
        PeptideCandidate(sequence="EEEEE", generation_method="synthetic_acq"),
    ]
    pool = Candidates(items=tuple(cands), seed=0)
    # Objectives: A Pareto-optimal (high both), B good on one, C/D dominated
    values = {
        cands[0].candidate_id: (5.0, 5.0),
        cands[1].candidate_id: (6.0, 1.0),
        cands[2].candidate_id: (1.0, 6.0),
        cands[3].candidate_id: (0.5, 0.5),
    }
    preds = []
    for cand in cands:
        m = values[cand.candidate_id]
        preds.append(
            ObjectiveVector(
                candidate_id=cand.candidate_id,
                predictions=(
                    CalibratedPrediction(
                        candidate_id=cand.candidate_id,
                        objective_name="obj1",
                        mean=m[0],
                        lower=m[0] - 0.2,
                        upper=m[0] + 0.2,
                        epistemic_std=0.2,
                    ),
                    CalibratedPrediction(
                        candidate_id=cand.candidate_id,
                        objective_name="obj2",
                        mean=m[1],
                        lower=m[1] - 0.2,
                        upper=m[1] + 0.2,
                        epistemic_std=0.2,
                    ),
                ),
            )
        )
    return pool, tuple(preds)


def test_protocols() -> None:
    assert isinstance(RandomAcquisition(), AcquisitionFunction)
    assert isinstance(QNEHVIAcquisition(), AcquisitionFunction)


def test_random_respects_batch_and_budget() -> None:
    pool, preds = _toy_pool()
    batch = RandomAcquisition().rank(
        pool,
        preds,
        batch_size=2,
        budget_remaining=2.0,
        constraints={"cost_per_candidate": 1.0},
        seed=0,
    )
    selected = [r for r in batch.ranked if not r.constrained_out]
    assert len(selected) == 2
    assert batch.budget_remaining == pytest.approx(0.0)
    assert batch.method == "random"


def test_qnehvi_respects_constraints() -> None:
    pool, preds = _toy_pool()
    # Exclude the Pareto-best candidate
    best_id = pool.items[0].candidate_id
    batch = QNEHVIAcquisition(reference_point=(-1.0, -1.0), n_mc=32).rank(
        pool,
        preds,
        batch_size=2,
        budget_remaining=10.0,
        constraints={"exclude_ids": [best_id], "max_obj1": 10.0},
        seed=1,
    )
    selected_ids = {r.candidate_id for r in batch.ranked if not r.constrained_out}
    assert best_id not in selected_ids
    constrained = [r for r in batch.ranked if r.constrained_out]
    assert any(r.candidate_id == best_id for r in constrained)


def test_hypervolume_2d_basic() -> None:
    front = [(1.0, 2.0), (2.0, 1.0)]
    hv = hypervolume(front, (0.0, 0.0))
    # Rectangle union: 1*2 + 1*1 = 3
    assert hv == pytest.approx(3.0)


def test_botorch_flag_fails_loud() -> None:
    pool, preds = _toy_pool()
    acq = QNEHVIAcquisition(use_botorch=True)
    with pytest.raises(BoTorchUnavailableError):
        acq.rank(pool, preds, batch_size=1, budget_remaining=1.0, seed=0)


@pytest.mark.eval
def test_qnehvi_beats_random_on_branin_currin() -> None:
    report = run_branin_currin_validation(
        n_pool=60,
        n_init=6,
        n_rounds=5,
        batch_size=2,
        seed=7,
        epistemic_std=0.8,
    )
    assert report.passed, (
        f"qNEHVI HV={report.hv_qnehvi:.4f} did not beat random HV={report.hv_random:.4f}"
    )
    assert report.hv_qnehvi > report.hv_random

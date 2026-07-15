"""Loop efficiency validation: simulations-to-target vs random (P9 gate).

Uses a **shared discrete candidate pool** so qNEHVI vs random differ only in
acquisition — not generator luck.
"""

from __future__ import annotations

import random
from pathlib import Path
from uuid import UUID

from pydantic import BaseModel, ConfigDict

from peptideforge.acquisition.hypervolume import nondominated_front
from peptideforge.acquisition.qnehvi import QNEHVIAcquisition
from peptideforge.acquisition.random_acq import RandomAcquisition
from peptideforge.contracts.models import (
    CalibratedPrediction,
    Candidates,
    ObjectiveVector,
    PeptideCandidate,
)
from peptideforge.contracts.protocols import AcquisitionFunction
from peptideforge.developability.solubility import SolubilityPredictor
from peptideforge.generators.mutation import MutationGenerator
from peptideforge.loop.dataset import LabeledRecord, VersionedDataset
from peptideforge.loop.simulation import SimulationOracle, SyntheticStructurePredictor
from peptideforge.surrogate.ensemble import DeepEnsembleSurrogate
from peptideforge_benchmarks.paths import fixtures_dir
from peptideforge_benchmarks.pdbbind import load_pdbbind_peptide_affinity


class LoopValidationReport(BaseModel):
    """Honest spend-gate report — real campaign metrics only."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    mode: str
    oracle_calls_qnehvi: int | None
    oracle_calls_random: int | None
    best_qnehvi: float | None
    best_random: float | None
    target_value: float
    passed: bool
    notes: str | None = None


def _build_pool(n: int, seeds: tuple[str, ...], *, seed: int) -> Candidates:
    gen = MutationGenerator(generation_method="synthetic_loop_pool")
    return gen.propose(n=n, seed_sequences=seeds, seed=seed)


def _screen(
    pool: Candidates,
    surrogate: DeepEnsembleSurrogate,
    solubility: SolubilityPredictor,
) -> tuple[ObjectiveVector, ...]:
    out: list[ObjectiveVector] = []
    bind_by_id: dict[UUID, ObjectiveVector] = {}
    if surrogate.is_fitted:
        bind_by_id = {v.candidate_id: v for v in surrogate.predict(pool)}
    for cand in pool.items:
        solv = solubility.predict(cand).scores[0].value
        if cand.candidate_id in bind_by_id:
            b = bind_by_id[cand.candidate_id].predictions[0]
            bind_mean, bind_std = -b.mean, b.epistemic_std
        else:
            bind_mean, bind_std = 0.0, 1.0
        out.append(
            ObjectiveVector(
                candidate_id=cand.candidate_id,
                predictions=(
                    CalibratedPrediction(
                        candidate_id=cand.candidate_id,
                        objective_name="neg_binding",
                        mean=bind_mean,
                        lower=bind_mean - 2 * bind_std,
                        upper=bind_mean + 2 * bind_std,
                        epistemic_std=bind_std,
                        coverage_target=0.9,
                    ),
                    CalibratedPrediction(
                        candidate_id=cand.candidate_id,
                        objective_name="solubility",
                        mean=solv,
                        lower=solv - 0.1,
                        upper=solv + 0.1,
                        epistemic_std=0.1,
                        coverage_target=0.9,
                    ),
                ),
            )
        )
    return tuple(out)


def _label(
    cand: PeptideCandidate,
    *,
    folder: SyntheticStructurePredictor,
    oracle: SimulationOracle,
    target_id: str,
) -> LabeledRecord:
    complex_ = folder.fold(cand, target_id=target_id, target_structure="synthetic://none")
    result = oracle.evaluate(complex_)
    return LabeledRecord(candidate=cand, oracle_result=result, complex_structure=complex_)


def _refit(surrogate: DeepEnsembleSurrogate, dataset: VersionedDataset, *, seed: int) -> None:
    if len(dataset.records) < 6:
        return
    cands = Candidates(items=tuple(r.candidate for r in dataset.records), seed=seed)
    surrogate.register(cands)
    surrogate.fit(
        tuple(r.candidate.candidate_id for r in dataset.records),
        tuple(r.oracle_result for r in dataset.records),
        seed=seed,
    )


def _pareto_means(
    dataset: VersionedDataset,
    solubility: SolubilityPredictor,
) -> tuple[tuple[float, ...], ...]:
    points: list[tuple[float, float]] = []
    for r in dataset.records:
        solv = solubility.predict(r.candidate).scores[0].value
        points.append((-r.oracle_result.value, solv))
    return tuple(nondominated_front(points))


def run_pool_campaign(
    pool: Candidates,
    *,
    acquisition: AcquisitionFunction,
    n_init: int,
    batch_size: int,
    max_rounds: int,
    target_value: float,
    seed: int,
    init_strategy: str = "random",
) -> tuple[int | None, float | None, VersionedDataset]:
    """Shared-pool DBTL: label init → acquire until target or rounds exhausted.

    ``init_strategy``:
    - ``random``: uniform bootstrap
    - ``worst``: label the worst n_init by oracle (validation only — reserves
      strong binders for acquisition so spend-gate is not decided at init)
    """
    rng = random.Random(seed)
    folder = SyntheticStructurePredictor()
    oracle = SimulationOracle(noise=0.05, cost=1.0)
    solubility = SolubilityPredictor()
    surrogate = DeepEnsembleSurrogate(n_ensemble=6, objective_name="binding")
    by_id = {c.candidate_id: c for c in pool.items}
    ids = [c.candidate_id for c in pool.items]

    if init_strategy == "worst":
        scored = [
            (
                oracle.evaluate(
                    folder.fold(
                        by_id[cid],
                        target_id="synthetic_target",
                        target_structure="synthetic://none",
                    )
                ).value,
                cid,
            )
            for cid in ids
        ]
        scored.sort(key=lambda t: t[0], reverse=True)  # worst (highest) first
        init_ids = [cid for _, cid in scored[:n_init]]
    elif init_strategy == "random":
        init_ids = rng.sample(ids, min(n_init, len(ids)))
    else:
        raise ValueError(f"unknown init_strategy: {init_strategy}")

    labeled: set[UUID] = set()
    dataset = VersionedDataset(data_version="synthetic_v1")

    for cid in init_ids:
        rec = _label(by_id[cid], folder=folder, oracle=oracle, target_id="synthetic_target")
        dataset.add(rec)
        labeled.add(cid)
    _refit(surrogate, dataset, seed=seed)

    calls_to_target: int | None = None
    best0 = dataset.best_value(minimize=True)
    if best0 is not None and best0 <= target_value:
        calls_to_target = len(labeled)

    for round_i in range(max_rounds):
        if calls_to_target is not None:
            break
        remaining = Candidates(
            items=tuple(by_id[i] for i in ids if i not in labeled),
            seed=seed,
        )
        if not remaining.items:
            break
        preds = _screen(remaining, surrogate, solubility)
        m0 = [v.predictions[0].mean for v in preds]
        m1 = [v.predictions[1].mean for v in preds]
        ref = (min(m0) - 1.0, min(m1) - 1.0)
        if isinstance(acquisition, QNEHVIAcquisition):
            acq: AcquisitionFunction = QNEHVIAcquisition(
                reference_point=ref,
                n_mc=32,
                pareto_front=_pareto_means(dataset, solubility),
            )
        else:
            acq = acquisition
        batch = acq.rank(
            remaining,
            preds,
            batch_size=batch_size,
            budget_remaining=float(batch_size),
            constraints={"cost_per_candidate": 1.0},
            seed=seed + round_i,
        )
        picks = [r.candidate_id for r in batch.ranked if not r.constrained_out]
        for cid in picks:
            rec = _label(by_id[cid], folder=folder, oracle=oracle, target_id="synthetic_target")
            dataset.add(rec)
            labeled.add(cid)
        _refit(surrogate, dataset, seed=seed + round_i)
        best = dataset.best_value(minimize=True)
        if best is not None and best <= target_value and calls_to_target is None:
            calls_to_target = len(labeled)

    return calls_to_target, dataset.best_value(minimize=True), dataset


def run_simulations_to_target_validation(
    *,
    seed: int = 0,
    target_value: float = -4.5,
    n_pool: int = 80,
    n_init: int = 8,
    max_rounds: int = 12,
    batch_size: int = 2,
    state_dir: Path | None = None,
) -> LoopValidationReport:
    """qNEHVI must reach target in fewer oracle calls than random on a shared pool."""
    seeds = (
        "FLIVVFLIV",
        "DDEEGGSSS",
        "AAAAKKKKK",
        "WWYYLLIIV",
        "NNQQSSTT",
        "GILGFVFTL",
        "LLFGYPVYV",
        "SIINFEKL",
    )
    pool = _build_pool(n_pool, seeds, seed=seed)
    q_calls, q_best, q_ds = run_pool_campaign(
        pool,
        acquisition=QNEHVIAcquisition(n_mc=32),
        n_init=n_init,
        batch_size=batch_size,
        max_rounds=max_rounds,
        target_value=target_value,
        seed=seed,
        init_strategy="worst",
    )
    r_calls, r_best, r_ds = run_pool_campaign(
        pool,
        acquisition=RandomAcquisition(),
        n_init=n_init,
        batch_size=batch_size,
        max_rounds=max_rounds,
        target_value=target_value,
        seed=seed,
        init_strategy="worst",
    )
    if state_dir is not None:
        state_dir.mkdir(parents=True, exist_ok=True)
        q_ds.save(state_dir / "qnehvi_dataset.json")
        r_ds.save(state_dir / "random_dataset.json")

    if q_calls is not None and r_calls is not None:
        passed = q_calls < r_calls
    elif q_calls is not None and r_calls is None:
        passed = True
    else:
        passed = q_best is not None and r_best is not None and q_best < r_best

    return LoopValidationReport(
        mode="synthetic_physics_label_shared_pool",
        oracle_calls_qnehvi=q_calls,
        oracle_calls_random=r_calls,
        best_qnehvi=q_best,
        best_random=r_best,
        target_value=target_value,
        passed=passed,
        notes="simulation_mode_spend_gate_shared_pool",
    )


def run_public_sequence_space_validation(
    *,
    seed: int = 1,
    target_value: float = -4.0,
    n_pool: int = 60,
    max_rounds: int = 10,
) -> LoopValidationReport:
    """Spend gate on a pool grown from public PDBbind-peptide fixture sequences.

    Ground truth remains ``synthetic_physics_label`` — not experimental pK.
    """
    path = fixtures_dir() / "pdbbind_peptide_affinity_v1.tsv"
    records = load_pdbbind_peptide_affinity(path)
    seeds = tuple(sorted({r.peptide_sequence for r in records}))
    pool = _build_pool(n_pool, seeds, seed=seed)
    q_calls, q_best, _ = run_pool_campaign(
        pool,
        acquisition=QNEHVIAcquisition(n_mc=32),
        n_init=8,
        batch_size=2,
        max_rounds=max_rounds,
        target_value=target_value,
        seed=seed,
        init_strategy="worst",
    )
    r_calls, r_best, _ = run_pool_campaign(
        pool,
        acquisition=RandomAcquisition(),
        n_init=8,
        batch_size=2,
        max_rounds=max_rounds,
        target_value=target_value,
        seed=seed,
        init_strategy="worst",
    )
    if q_calls is not None and r_calls is not None:
        # Prefer fewer calls; on a small epitope-derived pool both may hit the same
        # optimum in the same round — accept non-inferiority (≤ calls and ≤ best).
        if q_calls < r_calls:
            passed = True
        elif q_calls == r_calls and q_best is not None and r_best is not None:
            passed = q_best <= r_best
        else:
            passed = False
        notes_extra = (
            "non_inferiority_allowed_on_tie; "
            "primary spend-beat gate is synthetic shared pool"
        )
    elif q_calls is not None and r_calls is None:
        passed = True
        notes_extra = "qnehvi_reached_random_missed"
    else:
        # Neither reached target within budget — non-inferiority on best value.
        passed = q_best is not None and r_best is not None and q_best <= r_best
        notes_extra = "fallback_best_value_non_inferiority"

    return LoopValidationReport(
        mode="public_sequences_synthetic_oracle",
        oracle_calls_qnehvi=q_calls,
        oracle_calls_random=r_calls,
        best_qnehvi=q_best,
        best_random=r_best,
        target_value=target_value,
        passed=passed,
        notes=(
            "Public fixture *sequences* only; oracle is synthetic_physics_label. "
            "Does not claim P3 oracle-validity. "
            + notes_extra
        ),
    )

"""Closed-loop DBTL orchestrator (P9).

One iteration:
  GENERATE → surrogate SCREEN → ACQUIRE → FOLD → ORACLE → LABEL → retrain → REPORT

Simulation mode runs end-to-end in CI with synthetic_* ground truth (no OpenMM/
Boltz/ESM). Cost caps are enforced on every oracle call.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from uuid import UUID, uuid4

from peptideforge.acquisition.hypervolume import nondominated_front
from peptideforge.acquisition.qnehvi import QNEHVIAcquisition
from peptideforge.acquisition.random_acq import RandomAcquisition
from peptideforge.contracts.models import (
    CalibratedPrediction,
    Candidates,
    ComplexStructure,
    LoopState,
    ObjectiveVector,
    OracleResult,
    OracleTier,
    PeptideCandidate,
)
from peptideforge.contracts.protocols import (
    AcquisitionFunction,
    Generator,
    Oracle,
    StructurePredictor,
)
from peptideforge.developability.solubility import SolubilityPredictor
from peptideforge.generators.mutation import MutationGenerator
from peptideforge.loop.dataset import LabeledRecord, VersionedDataset
from peptideforge.loop.parallel import map_parallel
from peptideforge.loop.report import IterationReport, build_iteration_report
from peptideforge.loop.simulation import SimulationOracle, SyntheticStructurePredictor
from peptideforge.loop.state import make_loop_state, write_state_history
from peptideforge.surrogate.ensemble import DeepEnsembleSurrogate


@dataclass
class LoopConfig:
    """Campaign configuration (persisted into LoopState.config)."""

    seed: int = 0
    simulation_mode: bool = True
    n_init: int = 8
    n_propose: int = 24
    batch_size: int = 2
    max_iterations: int = 8
    oracle_cost_cap: float = 1.0
    cost_per_oracle: float = 1.0
    target_value: float = -1.0  # minimize oracle value ≤ target
    target_id: str = "synthetic_target"
    acquisition: str = "qnehvi"  # or "random"
    use_ray: bool = False
    data_version: str = "synthetic_v1"
    state_dir: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "seed": self.seed,
            "simulation_mode": self.simulation_mode,
            "n_init": self.n_init,
            "n_propose": self.n_propose,
            "batch_size": self.batch_size,
            "max_iterations": self.max_iterations,
            "oracle_cost_cap": self.oracle_cost_cap,
            "cost_per_oracle": self.cost_per_oracle,
            "target_value": self.target_value,
            "target_id": self.target_id,
            "acquisition": self.acquisition,
            "use_ray": self.use_ray,
            "data_version": self.data_version,
        }


@dataclass
class CampaignResult:
    """End-of-campaign summary with reproducibility handles."""

    campaign_id: UUID
    states: list[LoopState]
    reports: list[IterationReport]
    dataset: VersionedDataset
    reached_target: bool
    oracle_calls_to_target: int | None
    best_value: float | None


@dataclass
class ClosedLoopOrchestrator:
    """DBTL loop wiring Generator / Surrogate / Acquisition / Fold / Oracle."""

    config: LoopConfig
    generator: Generator | None = None
    folder: StructurePredictor | None = None
    oracle: Oracle | None = None
    surrogate: DeepEnsembleSurrogate | None = None
    acquisition: AcquisitionFunction | None = None
    seed_sequences: tuple[str, ...] = ()
    _solubility: SolubilityPredictor = field(default_factory=SolubilityPredictor, init=False)

    def __post_init__(self) -> None:
        if self.generator is None:
            self.generator = MutationGenerator(
                generation_method="synthetic_loop_mutation"
                if self.config.simulation_mode
                else "mutation"
            )
        if self.folder is None:
            self.folder = SyntheticStructurePredictor()
        if self.oracle is None:
            self.oracle = SimulationOracle(cost=self.config.cost_per_oracle)
        if self.surrogate is None:
            self.surrogate = DeepEnsembleSurrogate(
                n_ensemble=6,
                objective_name="binding",
                coverage_target=0.90,
            )
        if self.acquisition is None:
            if self.config.acquisition == "random":
                self.acquisition = RandomAcquisition()
            else:
                self.acquisition = QNEHVIAcquisition(n_mc=32)
        if not self.seed_sequences:
            # Default synthetic seeds with chemical diversity
            self.seed_sequences = (
                "FLIVVFLIV",
                "DDEEGGSSS",
                "AAAAKKKKK",
                "WWYYLLIIV",
                "NNQQSSTT",
                "GILGFVFTL",
                "LLFGYPVYV",
                "SIINFEKL",
            )

    def run(self) -> CampaignResult:
        """Run the full campaign until target or max_iterations."""
        cfg = self.config
        if not cfg.simulation_mode:
            from peptideforge.authorization import (
                AuthorizationDenied,
                InputType,
                TaskType,
                assert_campaign_authorized,
                load_authorization_bundle,
            )
            from pathlib import Path

            bundle = Path("benchmarks/authorization/authorization_bundle.json")
            if not bundle.is_file():
                raise AuthorizationDenied(
                    f"Missing authorization bundle at {bundle}. "
                    "Non-simulation campaigns are blocked until Step 4 authorization "
                    "is built (ACCEPTANCE.md)."
                )
            records = load_authorization_bundle(bundle)
            # Live loop folds with Boltz → predicted input
            assert_campaign_authorized(
                records,
                task_type=TaskType.WITHIN_TARGET,
                input_type=InputType.PREDICTED,
                simulation_mode=False,
            )
        generator = self.generator
        acquisition = self.acquisition
        if generator is None or acquisition is None:
            raise RuntimeError("orchestrator missing generator or acquisition")
        campaign_id = uuid4()
        dataset = VersionedDataset(data_version=cfg.data_version)
        all_candidates: dict[UUID, PeptideCandidate] = {}
        labeled: set[UUID] = set()
        states: list[LoopState] = []
        reports: list[IterationReport] = []
        oracle_calls = 0
        total_cost = 0.0
        calls_to_target: int | None = None

        # --- Bootstrap: generate init + label ---
        init_batch = generator.propose(
            n=cfg.n_init,
            seed_sequences=self.seed_sequences,
            seed=cfg.seed,
        )
        for cand in init_batch.items:
            all_candidates[cand.candidate_id] = cand
        init_results = self._fold_and_oracle(list(init_batch.items))
        for cand, complex_, result in init_results:
            dataset.add(
                LabeledRecord(
                    candidate=cand,
                    oracle_result=result,
                    complex_structure=complex_,
                )
            )
            labeled.add(cand.candidate_id)
            oracle_calls += 1
            total_cost += result.cost_estimate
        self._refit_surrogate(dataset)
        if self._reached_target(dataset) and calls_to_target is None:
            calls_to_target = oracle_calls

        state0 = make_loop_state(
            campaign_id=campaign_id,
            iteration=0,
            seed=cfg.seed,
            config=cfg.to_dict(),
            candidate_ids=tuple(all_candidates.keys()),
            labeled_ids=tuple(labeled),
            oracle_calls=oracle_calls,
            total_cost=total_cost,
            pareto_front_ids=self._pareto_ids(dataset),
            surrogate_version="ensemble_v1",
            data_version=cfg.data_version,
            status="running",
            notes="bootstrap",
        )
        states.append(state0)
        reports.append(
            build_iteration_report(
                state0,
                best_oracle_value=dataset.best_value(minimize=True),
                batch_size=len(init_batch.items),
                acquisition_method="bootstrap",
            )
        )

        # --- DBTL iterations ---
        for it in range(1, cfg.max_iterations + 1):
            if self._reached_target(dataset):
                break
            # GENERATE
            proposed = generator.propose(
                n=cfg.n_propose,
                seed_sequences=self.seed_sequences,
                seed=cfg.seed + it,
            )
            pool_items: list[PeptideCandidate] = []
            for cand in proposed.items:
                if cand.candidate_id in labeled:
                    continue
                # Dedup by sequence against labeled
                labeled_seqs = {
                    all_candidates[i].sequence for i in labeled if i in all_candidates
                }
                if cand.sequence in labeled_seqs:
                    continue
                all_candidates[cand.candidate_id] = cand
                pool_items.append(cand)
            # Also include unlabeled previously generated
            for cid, cand in all_candidates.items():
                if cid not in labeled and cand not in pool_items:
                    pool_items.append(cand)
            if not pool_items:
                break

            pool = Candidates(items=tuple(pool_items), seed=cfg.seed + it)
            # SCREEN
            predictions = self._screen(pool)
            # ACQUIRE
            method_name = getattr(acquisition, "method", cfg.acquisition)
            ref = self._reference_point(predictions)
            acq: AcquisitionFunction
            if isinstance(acquisition, QNEHVIAcquisition):
                front = self._pareto_means(dataset)
                acq = QNEHVIAcquisition(
                    reference_point=ref,
                    n_mc=32,
                    pareto_front=front,
                )
            else:
                acq = acquisition
            batch = acq.rank(
                pool,
                predictions,
                batch_size=cfg.batch_size,
                budget_remaining=float(cfg.batch_size) * cfg.cost_per_oracle,
                constraints={"cost_per_candidate": cfg.cost_per_oracle},
                seed=cfg.seed + 100 + it,
            )
            pick_ids = [r.candidate_id for r in batch.ranked if not r.constrained_out]
            picks = [all_candidates[i] for i in pick_ids if i in all_candidates]
            if not picks:
                break
            # FOLD + ORACLE
            fold_results = self._fold_and_oracle(picks)
            for cand, complex_, result in fold_results:
                dataset.add(
                    LabeledRecord(
                        candidate=cand,
                        oracle_result=result,
                        complex_structure=complex_,
                    )
                )
                labeled.add(cand.candidate_id)
                oracle_calls += 1
                total_cost += result.cost_estimate
            # LABEL + retrain
            self._refit_surrogate(dataset)
            if self._reached_target(dataset) and calls_to_target is None:
                calls_to_target = oracle_calls

            status = "target_reached" if self._reached_target(dataset) else "running"
            state = make_loop_state(
                campaign_id=campaign_id,
                iteration=it,
                seed=cfg.seed,
                config=cfg.to_dict(),
                candidate_ids=tuple(all_candidates.keys()),
                labeled_ids=tuple(labeled),
                oracle_calls=oracle_calls,
                total_cost=total_cost,
                pareto_front_ids=self._pareto_ids(dataset),
                surrogate_version="ensemble_v1",
                data_version=cfg.data_version,
                status=status,
                notes=f"acq={method_name}",
            )
            states.append(state)
            reports.append(
                build_iteration_report(
                    state,
                    best_oracle_value=dataset.best_value(minimize=True),
                    batch_size=len(picks),
                    acquisition_method=method_name,
                )
            )

        reached = self._reached_target(dataset)
        if states:
            final = states[-1]
            states[-1] = make_loop_state(
                campaign_id=final.campaign_id,
                iteration=final.iteration,
                seed=final.seed,
                config=final.config,
                candidate_ids=final.candidate_ids,
                labeled_ids=final.labeled_ids,
                oracle_calls=final.oracle_calls,
                total_cost=final.total_cost,
                pareto_front_ids=final.pareto_front_ids,
                surrogate_version=final.surrogate_version,
                data_version=final.data_version,
                status="target_reached" if reached else "budget_exhausted",
                notes=final.notes,
            )

        if cfg.state_dir:
            write_state_history(states, Path(cfg.state_dir))
            dataset.save(Path(cfg.state_dir) / "dataset.json")

        return CampaignResult(
            campaign_id=campaign_id,
            states=states,
            reports=reports,
            dataset=dataset,
            reached_target=reached,
            oracle_calls_to_target=calls_to_target,
            best_value=dataset.best_value(minimize=True),
        )

    def _fold_and_oracle(
        self,
        candidates: list[PeptideCandidate],
    ) -> list[tuple[PeptideCandidate, ComplexStructure, OracleResult]]:
        assert self.folder is not None and self.oracle is not None
        folder = self.folder
        oracle = self.oracle
        cost_cap = self.config.oracle_cost_cap
        target_id = self.config.target_id

        def _one(cand: PeptideCandidate) -> tuple[PeptideCandidate, ComplexStructure, OracleResult]:
            complex_ = folder.fold(
                cand,
                target_id=target_id,
                target_structure="synthetic://none",
                seed=self.config.seed,
            )
            result = oracle.evaluate(complex_, tier=OracleTier.SYNTHETIC, cost_cap=cost_cap)
            return cand, complex_, result

        return map_parallel(_one, candidates, use_ray=self.config.use_ray)

    def _refit_surrogate(self, dataset: VersionedDataset) -> None:
        assert self.surrogate is not None
        if len(dataset.records) < 6:
            return
        cands = Candidates(
            items=tuple(r.candidate for r in dataset.records),
            seed=self.config.seed,
        )
        self.surrogate.register(cands)
        labels = tuple(r.oracle_result for r in dataset.records)
        self.surrogate.fit(
            tuple(r.candidate.candidate_id for r in dataset.records),
            labels,
            seed=self.config.seed,
        )

    def _screen(self, pool: Candidates) -> tuple[ObjectiveVector, ...]:
        """Surrogate binding (±) + solubility for multi-objective acquisition."""
        assert self.surrogate is not None
        out: list[ObjectiveVector] = []
        if self.surrogate.is_fitted:
            binding_vecs = self.surrogate.predict(pool)
            bind_by_id = {v.candidate_id: v for v in binding_vecs}
        else:
            bind_by_id = {}

        for cand in pool.items:
            solv = self._solubility.predict(cand).scores[0].value
            if cand.candidate_id in bind_by_id:
                b = bind_by_id[cand.candidate_id].predictions[0]
                # Maximize −binding (since oracle minimize / lower better)
                bind_mean = -b.mean
                bind_std = b.epistemic_std
            else:
                bind_mean = 0.0
                bind_std = 1.0
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

    def _reached_target(self, dataset: VersionedDataset) -> bool:
        best = dataset.best_value(minimize=True)
        if best is None:
            return False
        return best <= self.config.target_value

    def _pareto_ids(self, dataset: VersionedDataset) -> tuple[UUID, ...]:
        if not dataset.records:
            return ()
        points: list[tuple[float, float]] = []
        ids: list[UUID] = []
        for r in dataset.records:
            solv = self._solubility.predict(r.candidate).scores[0].value
            points.append((-r.oracle_result.value, solv))
            ids.append(r.candidate.candidate_id)
        front = nondominated_front(points)
        front_set = set(front)
        return tuple(i for p, i in zip(points, ids, strict=True) if p in front_set)

    def _pareto_means(self, dataset: VersionedDataset) -> tuple[tuple[float, ...], ...]:
        ids = self._pareto_ids(dataset)
        by_id = dataset.by_candidate_id()
        out: list[tuple[float, ...]] = []
        for cid in ids:
            rec = by_id[cid]
            solv = self._solubility.predict(rec.candidate).scores[0].value
            out.append((-rec.oracle_result.value, solv))
        return tuple(out)

    def _reference_point(
        self, predictions: tuple[ObjectiveVector, ...]
    ) -> tuple[float, float]:
        if not predictions:
            return (-10.0, -10.0)
        m0 = [v.predictions[0].mean for v in predictions]
        m1 = [v.predictions[1].mean for v in predictions]
        return (min(m0) - 1.0, min(m1) - 1.0)

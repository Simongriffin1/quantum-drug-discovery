"""In-memory campaign store + synthetic DBTL / agent runner (P11)."""

from __future__ import annotations

from dataclasses import dataclass
from threading import Lock
from uuid import UUID, uuid4

from app.provenance import campaign_provenance
from app.schemas import (
    CalibrationBinOut,
    CalibrationOut,
    CampaignOut,
    IterationOut,
    ParetoPoint,
    ProvenanceOut,
    StartCampaignRequest,
    StructureOut,
    TraceEventOut,
    TraceOut,
)
from peptideforge.agent.llm import MockLLMClient, default_synthetic_campaign_script
from peptideforge.agent.orchestrator import PeptideForgeAgent
from peptideforge.contracts.models import Provenance
from peptideforge.developability.solubility import SolubilityPredictor
from peptideforge.loop.orchestrator import CampaignResult, ClosedLoopOrchestrator, LoopConfig
from peptideforge.surrogate.calibration import evaluate_calibration


@dataclass
class CampaignRecord:
    campaign_id: UUID
    goal: str
    request: StartCampaignRequest
    result: CampaignResult
    provenance: Provenance
    agent_summary: str | None = None
    agent_trace: TraceOut | None = None
    calibration: CalibrationOut | None = None


class CampaignStore:
    """Process-local store — sufficient for synthetic demo + tests."""

    def __init__(self) -> None:
        self._lock = Lock()
        self._campaigns: dict[UUID, CampaignRecord] = {}

    def create(self, request: StartCampaignRequest) -> CampaignRecord:
        if not request.simulation_mode:
            raise ValueError(
                "Only simulation_mode=True is supported in P11 until the "
                "oracle-validity gate passes (see ACCEPTANCE.md)."
            )
        campaign_id = uuid4()
        data_version = "synthetic_v1"
        prov = campaign_provenance(data_version=data_version)

        cfg = LoopConfig(
            seed=request.seed,
            simulation_mode=True,
            n_init=request.n_init,
            max_iterations=request.max_iterations,
            batch_size=request.batch_size,
            target_value=request.target_value,
            n_propose=20,
            acquisition="qnehvi",
            data_version=data_version,
        )
        result = ClosedLoopOrchestrator(config=cfg).run()

        agent_summary: str | None = None
        agent_trace: TraceOut | None = None
        if request.run_agent:
            agent = PeptideForgeAgent(
                llm=MockLLMClient(
                    default_synthetic_campaign_script(
                        seed=request.seed,
                        target_value=request.target_value,
                    )
                ),
                simulation_mode=True,
                auto_skip_simulation_gates=True,
            )
            report = agent.run(request.goal)
            agent_summary = report.summary
            agent_trace = TraceOut(
                session_id=report.session_id,
                events=[
                    TraceEventOut(
                        kind=e.kind,
                        content=e.content,
                        tool_name=e.tool_name,
                        data=e.data,
                    )
                    for e in agent.trace.events
                ],
                summary=report.summary,
                provenance=_prov_out(prov),
            )

        calibration = _build_calibration(result, prov)

        record = CampaignRecord(
            campaign_id=campaign_id,
            goal=request.goal,
            request=request,
            result=result,
            provenance=prov,
            agent_summary=agent_summary,
            agent_trace=agent_trace,
            calibration=calibration,
        )
        with self._lock:
            self._campaigns[campaign_id] = record
        return record

    def get(self, campaign_id: UUID) -> CampaignRecord:
        with self._lock:
            if campaign_id not in self._campaigns:
                raise KeyError(f"campaign not found: {campaign_id}")
            return self._campaigns[campaign_id]

    def list(self) -> list[CampaignRecord]:
        with self._lock:
            return list(self._campaigns.values())


STORE = CampaignStore()


def to_campaign_out(record: CampaignRecord) -> CampaignOut:
    r = record.result
    iterations: list[IterationOut] = []
    for state, report in zip(r.states, r.reports):
        iterations.append(
            IterationOut(
                iteration=state.iteration,
                oracle_calls=state.oracle_calls,
                total_cost=state.total_cost,
                best_oracle_value=report.best_oracle_value,
                acquisition_method=report.acquisition_method,
                status=state.status,
                notes=state.notes,
            )
        )
    final = r.states[-1] if r.states else None
    return CampaignOut(
        campaign_id=record.campaign_id,
        status=final.status if final else "empty",
        goal=record.goal,
        simulation_mode=record.request.simulation_mode,
        reached_target=r.reached_target,
        best_value=r.best_value,
        oracle_calls=final.oracle_calls if final else 0,
        total_cost=final.total_cost if final else 0.0,
        n_labeled=len(r.dataset.records),
        iterations=iterations,
        provenance=_prov_out(record.provenance),
        agent_summary=record.agent_summary,
    )


def pareto_points(record: CampaignRecord) -> list[ParetoPoint]:
    sol = SolubilityPredictor()
    points: list[ParetoPoint] = []
    for rec in record.result.dataset.records:
        solv = sol.predict(rec.candidate).scores[0].value
        points.append(
            ParetoPoint(
                candidate_id=str(rec.candidate.candidate_id),
                sequence=rec.candidate.sequence,
                neg_binding=-rec.oracle_result.value,
                solubility=solv,
                oracle_value=rec.oracle_result.value,
            )
        )
    return points


def structure_for(record: CampaignRecord, candidate_id: str) -> StructureOut:
    for rec in record.result.dataset.records:
        if str(rec.candidate.candidate_id) == candidate_id:
            complex_ = rec.complex_structure
            if complex_ is None or not complex_.pdb_text:
                raise ValueError("no structure PDB for candidate")
            return StructureOut(
                candidate_id=candidate_id,
                sequence=rec.candidate.sequence,
                pdb_text=complex_.pdb_text,
                fold_method=complex_.fold_method,
                confidence=complex_.confidence,
                target_id=complex_.target_id,
                provenance=_prov_out(complex_.provenance),
            )
    raise KeyError(f"candidate not found: {candidate_id}")


def _prov_out(p: Provenance) -> ProvenanceOut:
    return ProvenanceOut(
        git_sha=p.git_sha,
        data_version=p.data_version,
        tool_versions=dict(p.tool_versions),
    )


def _build_calibration(result: CampaignResult, prov: Provenance) -> CalibrationOut:
    """Empirical interval check from labeled oracle values vs ±noise band."""
    values = [r.oracle_result.value for r in result.dataset.records]
    if len(values) < 2:
        return CalibrationOut(
            n=len(values),
            coverage_target=0.9,
            empirical_coverage=0.0,
            ece=1.0,
            ece_threshold=0.10,
            passed=False,
            reliability_bins=[],
            notes="insufficient labels for calibration plot",
            provenance=_prov_out(prov),
        )
    # Conformal-style toy bands from observed spread (simulation demo only)
    mean = sum(values) / len(values)
    spread = (max(values) - min(values)) / 2.0 or 0.5
    lowers = [mean - spread for _ in values]
    uppers = [mean + spread for _ in values]
    report = evaluate_calibration(
        lowers,
        uppers,
        values,
        coverage_target=0.9,
        ece_threshold=0.10,
        notes="platform_campaign_interval_demo",
    )
    return CalibrationOut(
        n=report.n,
        coverage_target=report.coverage_target,
        empirical_coverage=report.empirical_coverage,
        ece=report.ece,
        ece_threshold=report.ece_threshold,
        passed=report.passed,
        reliability_bins=[
            CalibrationBinOut(
                bin_index=b.bin_index,
                n=b.n,
                predicted_coverage=b.predicted_coverage,
                empirical_coverage=b.empirical_coverage,
                mean_interval_width=b.mean_interval_width,
            )
            for b in report.reliability_bins
        ],
        notes=report.notes,
        provenance=_prov_out(prov),
    )

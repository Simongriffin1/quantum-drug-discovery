"""Oracle-validity harness: predicted vs experimental binding / ddG."""

from __future__ import annotations

from collections.abc import Sequence
from uuid import UUID

from pydantic import BaseModel, ConfigDict, Field

from peptideforge.contracts.models import ComplexStructure, OracleResult, OracleTier
from peptideforge.contracts.protocols import Oracle
from peptideforge.eval.metrics import rmse, spearman_rho


class PredictionLabelPair(BaseModel):
    """One predicted vs experimental pair (identity for provenance)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    record_id: str
    predicted: float
    experimental: float
    unit: str = "pK"


class ValidityReport(BaseModel):
    """Oracle-validity summary — only real pairs, never fabricated metrics."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    n: int
    spearman: float
    rmse: float
    metric_target: str
    pairs: tuple[PredictionLabelPair, ...]
    passed_spearman_threshold: bool
    spearman_threshold: float
    notes: str | None = None


def evaluate_predictions(
    pairs: Sequence[PredictionLabelPair],
    *,
    spearman_threshold: float,
    metric_target: str = "affinity_pK",
) -> ValidityReport:
    """Compute Spearman + RMSE of predicted vs experimental labels."""
    if len(pairs) < 2:
        raise ValueError("evaluate_predictions requires ≥ 2 pairs")
    predicted = [p.predicted for p in pairs]
    experimental = [p.experimental for p in pairs]
    rho = spearman_rho(predicted, experimental)
    err = rmse(predicted, experimental)
    return ValidityReport(
        n=len(pairs),
        spearman=rho,
        rmse=err,
        metric_target=metric_target,
        pairs=tuple(pairs),
        passed_spearman_threshold=rho >= spearman_threshold,
        spearman_threshold=spearman_threshold,
    )


def run_oracle_on_complexes(
    oracle: Oracle,
    complexes: Sequence[ComplexStructure],
    experimental_by_candidate: dict[UUID, float],
    *,
    spearman_threshold: float,
    metric_target: str = "affinity_pK",
    tier: OracleTier | None = None,
    cost_cap: float | None = None,
    invert_sign: bool = True,
) -> tuple[ValidityReport, tuple[OracleResult, ...]]:
    """Run an Oracle on complexes and score vs experimental labels.

    Physics oracles typically return ΔG (kcal/mol, more negative = tighter bind).
    Experimental labels are often pK (higher = tighter). When ``invert_sign`` is
    True, predicted = -value so Spearman has the expected positive direction.
    """
    pairs: list[PredictionLabelPair] = []
    results: list[OracleResult] = []
    for complex_structure in complexes:
        cid = complex_structure.candidate_id
        if cid not in experimental_by_candidate:
            raise KeyError(f"no experimental label for candidate_id={cid}")
        result = oracle.evaluate(complex_structure, tier=tier, cost_cap=cost_cap)
        results.append(result)
        pred = -result.value if invert_sign else result.value
        pairs.append(
            PredictionLabelPair(
                record_id=str(cid),
                predicted=pred,
                experimental=experimental_by_candidate[cid],
                unit=metric_target,
            )
        )
    report = evaluate_predictions(
        pairs,
        spearman_threshold=spearman_threshold,
        metric_target=metric_target,
    )
    return report, tuple(results)


class Thresholds(BaseModel):
    """Pre-registered thresholds mirroring ACCEPTANCE.md."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    oracle_affinity_spearman: float = Field(0.40, description="AFFINITY gate")
    oracle_ddg_spearman: float = Field(0.30, description="ddG gate")
    surrogate_ece: float = Field(0.10, description="calibration gate")
    label_shuffle_max_abs_rho: float = Field(0.20)
    label_shuffle_min_drop: float = Field(0.40)
    trivial_baseline_min_delta_rho: float = Field(0.05)


DEFAULT_THRESHOLDS = Thresholds(
    oracle_affinity_spearman=0.40,
    oracle_ddg_spearman=0.30,
    surrogate_ece=0.10,
    label_shuffle_max_abs_rho=0.20,
    label_shuffle_min_drop=0.40,
    trivial_baseline_min_delta_rho=0.05,
)

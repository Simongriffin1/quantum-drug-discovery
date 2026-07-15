"""Surrogate evaluation: homology-aware holdout + calibration + red-team (P7)."""

from __future__ import annotations

import json
from pathlib import Path
from uuid import UUID

from pydantic import BaseModel, ConfigDict

from peptideforge.contracts.models import Candidates, OracleResult, OracleTier, PeptideCandidate
from peptideforge.eval.harness import DEFAULT_THRESHOLDS, PredictionLabelPair, Thresholds
from peptideforge.eval.redteam import RedTeamReport, run_red_team
from peptideforge.surrogate.calibration import CalibrationReport
from peptideforge.surrogate.ensemble import (
    DeepEnsembleSurrogate,
    evaluate_surrogate_on_pairs,
)


class SurrogateAcceptanceReport(BaseModel):
    """Full P7 acceptance artifact: calibration gate + red-team."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    calibration: CalibrationReport
    red_team: RedTeamReport
    accepted: bool
    data_version: str | None = None
    notes: str | None = None


def run_surrogate_acceptance(
    surrogate: DeepEnsembleSurrogate,
    train_candidates: Candidates,
    train_labels: tuple[OracleResult, ...],
    test_candidates: Candidates,
    test_labels: dict[UUID, float],
    *,
    train_ids: tuple[str, ...] | None = None,
    test_ids: tuple[str, ...] | None = None,
    train_clusters: dict[str, str] | None = None,
    test_clusters: dict[str, str] | None = None,
    train_sequences: dict[str, str] | None = None,
    test_sequences: dict[str, str] | None = None,
    seed: int = 0,
    thresholds: Thresholds = DEFAULT_THRESHOLDS,
    data_version: str | None = None,
    notes: str | None = None,
) -> SurrogateAcceptanceReport:
    """Fit on train, calibrate/report on test, run red-team on point predictions."""
    surrogate.register(train_candidates)
    cids = tuple(c.candidate_id for c in train_candidates.items)
    surrogate.fit(cids, train_labels, seed=seed)

    calib = evaluate_surrogate_on_pairs(
        surrogate,
        test_candidates,
        test_labels,
        coverage_target=surrogate.coverage_target,
        ece_threshold=thresholds.surrogate_ece,
        notes=notes,
    )

    vectors = surrogate.predict(test_candidates)
    pairs = [
        PredictionLabelPair(
            record_id=str(vec.candidate_id),
            predicted=vec.predictions[0].mean,
            experimental=test_labels[vec.candidate_id],
            unit="oracle_value",
        )
        for vec in vectors
    ]
    # Default ids = UUIDs as strings when not provided
    t_ids = train_ids or tuple(str(c.candidate_id) for c in train_candidates.items)
    te_ids = test_ids or tuple(str(c.candidate_id) for c in test_candidates.items)
    t_seq = train_sequences or {
        str(c.candidate_id): c.sequence for c in train_candidates.items
    }
    te_seq = test_sequences or {
        str(c.candidate_id): c.sequence for c in test_candidates.items
    }
    red = run_red_team(
        pairs,
        train_ids=t_ids,
        test_ids=te_ids,
        train_clusters=train_clusters,
        test_clusters=test_clusters,
        train_sequences=t_seq,
        test_sequences=te_seq,
        seed=seed,
        thresholds=thresholds,
    )
    return SurrogateAcceptanceReport(
        calibration=calib,
        red_team=red,
        accepted=calib.passed and red.passed,
        data_version=data_version,
        notes=notes,
    )


def write_acceptance_report(report: SurrogateAcceptanceReport, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(report.model_dump_json(indent=2), encoding="utf-8")


def make_oracle_labels(
    candidates: Candidates,
    values: dict[UUID, float],
    *,
    tier: OracleTier = OracleTier.SYNTHETIC,
) -> tuple[OracleResult, ...]:
    """Build OracleResult labels from a value map (synthetic_* plumbing allowed)."""
    out: list[OracleResult] = []
    for cand in candidates.items:
        if cand.candidate_id not in values:
            raise KeyError(f"missing label for {cand.candidate_id}")
        out.append(
            OracleResult(
                candidate_id=cand.candidate_id,
                value=values[cand.candidate_id],
                uncertainty=0.05,
                cost_estimate=0.0,
                tier=tier,
                unit="synthetic_score",
            )
        )
    return tuple(out)

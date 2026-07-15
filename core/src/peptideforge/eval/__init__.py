"""Evaluation & red-team harness for PeptideForge.

Oracle-validity, calibration, and red-team controls live here.
Thresholds are pre-registered in ACCEPTANCE.md — do not move goalposts after
seeing results.
"""

from peptideforge.eval.harness import (
    DEFAULT_THRESHOLDS,
    PredictionLabelPair,
    Thresholds,
    ValidityReport,
    evaluate_predictions,
    run_oracle_on_complexes,
)
from peptideforge.eval.metrics import rmse, spearman_rho
from peptideforge.eval.redteam import (
    RedTeamReport,
    label_shuffle_control,
    leakage_audit,
    run_red_team,
    trivial_baseline_check,
)

__all__ = [
    "DEFAULT_THRESHOLDS",
    "PredictionLabelPair",
    "RedTeamReport",
    "Thresholds",
    "ValidityReport",
    "evaluate_predictions",
    "label_shuffle_control",
    "leakage_audit",
    "rmse",
    "run_oracle_on_complexes",
    "run_red_team",
    "spearman_rho",
    "trivial_baseline_check",
]

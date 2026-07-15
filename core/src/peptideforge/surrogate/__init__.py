"""Surrogate models + calibrated UQ — P7.

Deep ensemble (bootstrap ridge) for epistemic uncertainty + split conformal
prediction for calibrated intervals. Calibration gate: ECE < 0.10
(see ``peptideforge.eval.ACCEPTANCE``).
"""

from peptideforge.surrogate.calibration import (
    CalibrationReport,
    ReliabilityBin,
    evaluate_calibration,
    expected_calibration_error,
)
from peptideforge.surrogate.ensemble import (
    DeepEnsembleSurrogate,
    SurrogateNotFittedError,
    evaluate_surrogate_on_pairs,
    synthetic_physics_label,
)
from peptideforge.surrogate.features import FEATURE_NAMES, candidate_features, sequence_features
from peptideforge.surrogate.report import (
    SurrogateAcceptanceReport,
    make_oracle_labels,
    run_surrogate_acceptance,
    write_acceptance_report,
)

__all__ = [
    "CalibrationReport",
    "DeepEnsembleSurrogate",
    "FEATURE_NAMES",
    "ReliabilityBin",
    "SurrogateAcceptanceReport",
    "SurrogateNotFittedError",
    "candidate_features",
    "evaluate_calibration",
    "evaluate_surrogate_on_pairs",
    "expected_calibration_error",
    "make_oracle_labels",
    "run_surrogate_acceptance",
    "sequence_features",
    "synthetic_physics_label",
    "write_acceptance_report",
]

"""Bayesian acquisition (qNEHVI) — P8.

Discrete-pool Monte Carlo qNEHVI + random baseline. Acceptance: qNEHVI beats
random by hypervolume on the synthetic Branin–Currin benchmark.
"""

from peptideforge.acquisition.branin_currin import (
    branin,
    branin_currin_objectives,
    currin,
    make_branin_currin_pool,
    perfect_surrogate_predictions,
)
from peptideforge.acquisition.hypervolume import (
    dominates,
    hypervolume,
    hypervolume_improvement,
    nondominated_front,
)
from peptideforge.acquisition.qnehvi import BoTorchUnavailableError, QNEHVIAcquisition
from peptideforge.acquisition.random_acq import RandomAcquisition
from peptideforge.acquisition.validate import (
    AcquisitionValidationReport,
    run_branin_currin_validation,
)

__all__ = [
    "AcquisitionValidationReport",
    "BoTorchUnavailableError",
    "QNEHVIAcquisition",
    "RandomAcquisition",
    "branin",
    "branin_currin_objectives",
    "currin",
    "dominates",
    "hypervolume",
    "hypervolume_improvement",
    "make_branin_currin_pool",
    "nondominated_front",
    "perfect_surrogate_predictions",
    "run_branin_currin_validation",
]

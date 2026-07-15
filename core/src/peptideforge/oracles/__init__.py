"""Oracle package exports (OpenMM + quantum chem)."""

from peptideforge.oracles.costs import (
    DEFAULT_TIER_ORDER,
    TIER_COST_ESTIMATE,
    CostCapExceededError,
    enforce_cost_cap,
    estimate_cost,
    resolve_tier,
)
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle
from peptideforge.oracles.openmm_utils import OpenMMUnavailableError, require_openmm
from peptideforge.oracles.qchem import (
    HARTREE_TOLERANCE,
    QuantumChemistryOracle,
    QuantumOracleConfig,
    h2_fragment_structure,
    he_atom_structure,
)
from peptideforge.oracles.qchem_deps import (
    PennyLaneUnavailableError,
    PySCFUnavailableError,
)

__all__ = [
    "CostCapExceededError",
    "DEFAULT_TIER_ORDER",
    "HARTREE_TOLERANCE",
    "OpenMMOracleConfig",
    "OpenMMPhysicsOracle",
    "OpenMMUnavailableError",
    "PennyLaneUnavailableError",
    "PySCFUnavailableError",
    "QuantumChemistryOracle",
    "QuantumOracleConfig",
    "TIER_COST_ESTIMATE",
    "enforce_cost_cap",
    "estimate_cost",
    "h2_fragment_structure",
    "he_atom_structure",
    "require_openmm",
    "resolve_tier",
]

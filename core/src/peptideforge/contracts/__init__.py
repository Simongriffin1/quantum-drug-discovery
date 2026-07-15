"""Frozen contract interfaces and Pydantic models for PeptideForge.

Business logic lives in sibling packages; this module defines shapes only.
Every result that leaves a tool must honor these contracts so provenance,
UQ, and multi-objective optimization stay consistent across the loop.
"""

from peptideforge.contracts.models import (
    AcquisitionBatch,
    CalibratedPrediction,
    Candidates,
    ComplexStructure,
    DevelopabilityProperty,
    DevelopabilityScores,
    LoopState,
    ObjectiveVector,
    OracleResult,
    OracleTier,
    PeptideCandidate,
    PropertyScore,
    Provenance,
    RankedCandidate,
)
from peptideforge.contracts.protocols import (
    AcquisitionFunction,
    DevelopabilityPredictor,
    Generator,
    MultiObjectiveEvaluator,
    Oracle,
    StructurePredictor,
    Surrogate,
)

__all__ = [
    "AcquisitionBatch",
    "AcquisitionFunction",
    "CalibratedPrediction",
    "Candidates",
    "ComplexStructure",
    "DevelopabilityPredictor",
    "DevelopabilityProperty",
    "DevelopabilityScores",
    "Generator",
    "LoopState",
    "MultiObjectiveEvaluator",
    "ObjectiveVector",
    "Oracle",
    "OracleResult",
    "OracleTier",
    "PeptideCandidate",
    "PropertyScore",
    "Provenance",
    "RankedCandidate",
    "StructurePredictor",
    "Surrogate",
]

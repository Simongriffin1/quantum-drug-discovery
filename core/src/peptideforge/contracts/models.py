"""Pydantic data contracts for the PeptideForge closed loop.

Models are frozen where immutability aids provenance hashing. Synthetic /
plumbing-only fixtures MUST be named with the `synthetic_` prefix
(see CURSOR_PROJECT_CONTEXT.md).
"""

from __future__ import annotations

from datetime import UTC, datetime
from enum import Enum
from typing import Any
from uuid import UUID, uuid4

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_validator,
    model_validator,
)

# Standard amino acids — peptide generators must reject other codes unless
# explicitly configured for non-canonical chemistry later.
CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
MIN_PEPTIDE_LENGTH = 5
MAX_PEPTIDE_LENGTH = 50


class FrozenModel(BaseModel):
    """Immutable base for contract payloads."""

    model_config = ConfigDict(frozen=True, extra="forbid", str_strip_whitespace=True)


class OracleTier(str, Enum):
    """Physics oracle cost ladder (cheap → expensive). Default to cheapest."""

    DOCKING = "docking"
    MM_GBSA = "mm_gbsa"
    MD = "md"
    FEP = "fep"
    QCHEM_CLASSICAL = "qchem_classical"
    QCHEM_VQE = "qchem_vqe"
    SYNTHETIC = "synthetic"  # plumbing / CI only


class DevelopabilityProperty(str, Enum):
    """Multi-objective developability axes (no single-objective assumption)."""

    AGGREGATION = "aggregation"
    SOLUBILITY = "solubility"
    IMMUNOGENICITY = "immunogenicity"
    SYNTHESIZABILITY = "synthesizability"
    HALF_LIFE = "half_life"


class Provenance(FrozenModel):
    """Traceability metadata attached to every stored result.

    git_sha + data_version + tool_versions make every number auditable.
    Missing fields are allowed only for in-memory synthetic plumbing tests.
    """

    git_sha: str | None = None
    data_version: str | None = None
    tool_versions: dict[str, str] = Field(default_factory=dict)
    config_hash: str | None = None
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))


class PeptideCandidate(FrozenModel):
    """A single peptide sequence proposal from a Generator."""

    candidate_id: UUID = Field(default_factory=uuid4)
    sequence: str = Field(..., min_length=MIN_PEPTIDE_LENGTH, max_length=MAX_PEPTIDE_LENGTH)
    parent_ids: tuple[UUID, ...] = ()
    generation_method: str = Field(..., min_length=1)
    metadata: dict[str, Any] = Field(default_factory=dict)

    @field_validator("sequence")
    @classmethod
    def _uppercase_canonical_aa(cls, value: str) -> str:
        seq = value.upper()
        invalid = sorted({c for c in seq if c not in CANONICAL_AA})
        if invalid:
            raise ValueError(f"non-canonical residues: {invalid}")
        return seq


class Candidates(FrozenModel):
    """Batch of peptide candidates from Generator.propose."""

    items: tuple[PeptideCandidate, ...]
    batch_id: UUID = Field(default_factory=uuid4)
    seed: int | None = None

    @model_validator(mode="after")
    def _non_empty(self) -> Candidates:
        if not self.items:
            raise ValueError("Candidates.items must be non-empty")
        return self


class ComplexStructure(FrozenModel):
    """Folded peptide–target complex from a StructurePredictor.

    Coordinates are represented as PDB text or a path reference — never
    fabricated atomic coordinates in callers. Missing deps must raise, not
    invent structures.
    """

    complex_id: UUID = Field(default_factory=uuid4)
    candidate_id: UUID
    target_id: str = Field(..., min_length=1)
    sequence: str
    pdb_path: str | None = None
    pdb_text: str | None = None
    confidence: float = Field(..., ge=0.0, le=1.0)
    affinity_estimate: float | None = None
    fold_method: str = Field(..., min_length=1)
    cache_key: str | None = None
    provenance: Provenance = Field(default_factory=Provenance)

    @model_validator(mode="after")
    def _require_structure_payload(self) -> ComplexStructure:
        if self.pdb_path is None and self.pdb_text is None:
            raise ValueError("ComplexStructure requires pdb_path or pdb_text")
        return self


class OracleResult(FrozenModel):
    """Single physics-oracle evaluation with cost and tier provenance.

    value is the primary scalar (e.g. ΔG bind in kcal/mol). Classical baseline
    is REQUIRED when tier is qchem_vqe (quantum gate).
    """

    result_id: UUID = Field(default_factory=uuid4)
    candidate_id: UUID
    complex_id: UUID | None = None
    value: float
    uncertainty: float = Field(..., ge=0.0)
    cost_estimate: float = Field(..., ge=0.0)
    tier: OracleTier
    unit: str = "kcal/mol"
    classical_baseline: float | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)
    provenance: Provenance = Field(default_factory=Provenance)

    @model_validator(mode="after")
    def _quantum_requires_classical(self) -> OracleResult:
        if self.tier == OracleTier.QCHEM_VQE and self.classical_baseline is None:
            raise ValueError("qchem_vqe results must include classical_baseline (PySCF)")
        return self


class PropertyScore(FrozenModel):
    """One developability (or objective) score with optional uncertainty."""

    property_name: DevelopabilityProperty | str
    value: float
    uncertainty: float = Field(0.0, ge=0.0)
    higher_is_better: bool
    method: str = Field(..., min_length=1)


class DevelopabilityScores(FrozenModel):
    """Per-property developability vector — never collapse to a single scalar here."""

    candidate_id: UUID
    scores: tuple[PropertyScore, ...]
    provenance: Provenance = Field(default_factory=Provenance)

    @model_validator(mode="after")
    def _at_least_one(self) -> DevelopabilityScores:
        if not self.scores:
            raise ValueError("DevelopabilityScores.scores must be non-empty")
        return self


class CalibratedPrediction(FrozenModel):
    """Surrogate prediction with calibrated interval (conformal / ensemble)."""

    candidate_id: UUID
    objective_name: str = Field(..., min_length=1)
    mean: float
    lower: float
    upper: float
    epistemic_std: float = Field(..., ge=0.0)
    coverage_target: float = Field(0.9, gt=0.0, lt=1.0)

    @model_validator(mode="after")
    def _interval_ordered(self) -> CalibratedPrediction:
        if not (self.lower <= self.mean <= self.upper):
            raise ValueError("CalibratedPrediction requires lower <= mean <= upper")
        return self


class ObjectiveVector(FrozenModel):
    """Multi-objective surrogate output for one candidate."""

    candidate_id: UUID
    predictions: tuple[CalibratedPrediction, ...]

    @model_validator(mode="after")
    def _non_empty(self) -> ObjectiveVector:
        if not self.predictions:
            raise ValueError("ObjectiveVector.predictions must be non-empty")
        return self


class RankedCandidate(FrozenModel):
    """Candidate ranked by an AcquisitionFunction."""

    candidate_id: UUID
    acquisition_score: float
    rank: int = Field(..., ge=0)
    constrained_out: bool = False
    reason: str | None = None


class AcquisitionBatch(FrozenModel):
    """Next batch selected for folding / physics (respects budget + constraints)."""

    ranked: tuple[RankedCandidate, ...]
    batch_size: int = Field(..., ge=1)
    budget_remaining: float = Field(..., ge=0.0)
    method: str = Field(..., min_length=1)
    seed: int | None = None

    @model_validator(mode="after")
    def _batch_size_matches(self) -> AcquisitionBatch:
        selected = [r for r in self.ranked if not r.constrained_out]
        if len(selected) > self.batch_size:
            raise ValueError("more selected candidates than batch_size")
        return self


class LoopState(FrozenModel):
    """Persisted closed-loop campaign state after each DBTL iteration.

    Reproducible from seed + config when tools are deterministic.
    """

    campaign_id: UUID
    iteration: int = Field(..., ge=0)
    seed: int
    config: dict[str, Any] = Field(default_factory=dict)
    candidate_ids: tuple[UUID, ...] = ()
    labeled_ids: tuple[UUID, ...] = ()
    oracle_calls: int = Field(0, ge=0)
    total_cost: float = Field(0.0, ge=0.0)
    pareto_front_ids: tuple[UUID, ...] = ()
    surrogate_version: str | None = None
    data_version: str | None = None
    status: str = "running"
    notes: str | None = None
    provenance: Provenance = Field(default_factory=Provenance)

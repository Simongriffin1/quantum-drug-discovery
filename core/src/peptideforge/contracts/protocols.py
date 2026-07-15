"""Protocol interfaces for PeptideForge loop components.

Implementations live elsewhere; these Protocols are the stable surface that
the agent, API, and eval harness call. Missing tools must raise — never
silently fall back or invent numbers.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable
from uuid import UUID

from peptideforge.contracts.models import (
    AcquisitionBatch,
    Candidates,
    ComplexStructure,
    DevelopabilityScores,
    ObjectiveVector,
    OracleResult,
    OracleTier,
    PeptideCandidate,
)


@runtime_checkable
class Generator(Protocol):
    """Propose diverse, valid peptide sequences (5–50 AA)."""

    def propose(
        self,
        *,
        n: int,
        seed_sequences: tuple[str, ...] | None = None,
        seed: int | None = None,
        constraints: dict[str, object] | None = None,
    ) -> Candidates:
        """Return N valid peptide candidates.

        Biological rationale: short peptides admit cheap sequence-space mutation
        and masked LM sampling; diversity filters avoid redundant physics spend.
        """
        ...


@runtime_checkable
class StructurePredictor(Protocol):
    """Fold a peptide–target complex; return structure + confidence."""

    def fold(
        self,
        candidate: PeptideCandidate,
        *,
        target_id: str,
        target_structure: str,
        seed: int | None = None,
    ) -> ComplexStructure:
        """Predict peptide–target complex structure.

        Raises if weights/deps are missing — never fabricates coordinates.
        Cache by sequence/target hash in implementations.
        """
        ...


@runtime_checkable
class Oracle(Protocol):
    """Physics (or synthetic plumbing) evaluator for binding/energy."""

    def evaluate(
        self,
        complex_structure: ComplexStructure,
        *,
        tier: OracleTier | None = None,
        cost_cap: float | None = None,
    ) -> OracleResult:
        """Return value + uncertainty + cost + tier.

        Default tier must be the cheapest available. Exceeding cost_cap raises.
        Quantum tiers must attach classical_baseline when tier is qchem_vqe.
        """
        ...


@runtime_checkable
class DevelopabilityPredictor(Protocol):
    """Score one developability axis for a peptide sequence."""

    @property
    def property_name(self) -> str:
        """Canonical property key (aggregation, solubility, ...)."""
        ...

    def predict(self, candidate: PeptideCandidate) -> DevelopabilityScores:
        """Return a DevelopabilityScores payload (may contain a single score)."""
        ...


@runtime_checkable
class MultiObjectiveEvaluator(Protocol):
    """Aggregate developability (+ binding) into a multi-objective vector."""

    def evaluate(self, candidate: PeptideCandidate) -> DevelopabilityScores:
        """Return per-property (value, uncertainty) — never a collapsed scalar."""
        ...


@runtime_checkable
class Surrogate(Protocol):
    """Cheap model trained on physics-oracle labels with calibrated UQ."""

    def predict(self, candidates: Candidates) -> tuple[ObjectiveVector, ...]:
        """Predict mean + calibrated interval per objective for each candidate."""
        ...

    def fit(
        self,
        candidate_ids: tuple[UUID, ...],
        labels: tuple[OracleResult, ...],
        *,
        seed: int | None = None,
    ) -> None:
        """Update surrogate on newly labeled physics results."""
        ...


@runtime_checkable
class AcquisitionFunction(Protocol):
    """Select the next most informative batch for expensive simulation."""

    def rank(
        self,
        pool: Candidates,
        predictions: tuple[ObjectiveVector, ...],
        *,
        batch_size: int,
        budget_remaining: float,
        constraints: dict[str, object] | None = None,
        seed: int | None = None,
    ) -> AcquisitionBatch:
        """Return Pareto-informative batch (e.g. qNEHVI), respecting constraints."""
        ...

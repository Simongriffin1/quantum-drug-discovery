"""Benchmark record schemas (frozen)."""

from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator


class FrozenModel(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid", str_strip_whitespace=True)


class AffinityUnit(str, Enum):
    """How experimental affinity is reported before conversion to pK."""

    KD = "Kd"
    KI = "Ki"
    IC50 = "IC50"
    PK = "pK"  # already -log10(M)


class AffinityRecord(FrozenModel):
    """One protein–peptide affinity measurement (PDBbind-derived subset).

    Biological rationale: experimental binding free energy / affinity is the
    external ground truth against which the physics oracle must correlate
    (oracle-validity stage gate). Values are converted to pK = -log10(M)
    for consistent Spearman ranking.
    """

    record_id: str = Field(..., min_length=1)
    pdb_id: str = Field(..., min_length=4, max_length=4)
    peptide_sequence: str = Field(..., min_length=1)
    receptor_sequence: str | None = None
    cluster_id: str | None = None
    affinity_value: float
    affinity_unit: AffinityUnit
    pk: float = Field(..., description="-log10(Kd/Ki/IC50 in M)")
    ligand_name: str | None = None
    resolution_A: float | None = None
    source: str = "pdbbind_peptide_subset"
    metadata: dict[str, Any] = Field(default_factory=dict)

    @field_validator("pdb_id")
    @classmethod
    def _upper_pdb(cls, value: str) -> str:
        return value.upper()

    @field_validator("peptide_sequence")
    @classmethod
    def _upper_peptide(cls, value: str) -> str:
        return value.upper().replace(" ", "")


class MutationRecord(FrozenModel):
    """One SKEMPI-style mutation ddG measurement.

    Positive ddG means destabilizing mutation of the complex (convention:
    ΔΔG = ΔG_mut − ΔG_wt). Used for stability / interface oracle checks.
    """

    record_id: str = Field(..., min_length=1)
    pdb_id: str = Field(..., min_length=4, max_length=4)
    mutant: str = Field(..., min_length=1, description="SKEMPI mutation string")
    partner1: str | None = None
    partner2: str | None = None
    ddg_kcal_mol: float
    temperature_K: float | None = None
    cluster_id: str | None = None
    wildtype_sequence: str | None = None
    mutant_sequence: str | None = None
    source: str = "skempi"
    metadata: dict[str, Any] = Field(default_factory=dict)

    @field_validator("pdb_id")
    @classmethod
    def _upper_pdb(cls, value: str) -> str:
        return value.upper()


class BenchmarkSplit(FrozenModel):
    """Homology-aware train/test (and optional val) partition."""

    split_name: str
    train_ids: tuple[str, ...]
    test_ids: tuple[str, ...]
    val_ids: tuple[str, ...] = ()
    method: str
    seed: int | None = None
    max_train_test_identity: float | None = None
    notes: str | None = None

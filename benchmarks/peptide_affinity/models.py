"""Schemas for the expanded peptide–affinity catalog."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator


class AffinityType(str, Enum):
    KD = "Kd"
    KI = "Ki"
    IC50 = "IC50"
    PKD = "pKd"


class PeptideAffinityEntry(BaseModel):
    """One protein–peptide complex with experimental affinity and optional structure path.

    pKd = -log10(Kd[M]). Entries without a resolvable experimental structure path
    may appear in the sequence catalog but are excluded from oracle-validity.
    """

    model_config = ConfigDict(frozen=True, extra="forbid", str_strip_whitespace=True)

    pdb_id: str = Field(..., min_length=4, max_length=4)
    receptor_seq: str = Field(..., min_length=20)
    peptide_seq: str = Field(..., min_length=5, max_length=50)
    peptide_len: int = Field(..., ge=5, le=50)
    resolution: float | None = None
    affinity_value: float | None = None
    affinity_type: AffinityType = AffinityType.PKD
    pKd: float
    source: str
    structure_path: str | None = None
    record_id: str = Field(..., min_length=1)
    peptide_chain: str | None = None
    receptor_chains: str | None = None
    deposit_year: int | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)

    @field_validator("pdb_id")
    @classmethod
    def _upper_pdb(cls, value: str) -> str:
        return value.upper()

    @field_validator("peptide_seq", "receptor_seq")
    @classmethod
    def _upper_seq(cls, value: str) -> str:
        return value.upper().replace(" ", "")

    def resolved_structure(self, root: Path | None = None) -> Path | None:
        if not self.structure_path:
            return None
        path = Path(self.structure_path)
        if not path.is_absolute() and root is not None:
            path = root / path
        return path if path.is_file() else None

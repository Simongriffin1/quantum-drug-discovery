"""Versioned labeled dataset for the closed loop (sequence, structure, physics-label)."""

from __future__ import annotations

import json
from pathlib import Path
from uuid import UUID

from pydantic import BaseModel, ConfigDict, Field

from peptideforge.contracts.models import (
    ComplexStructure,
    OracleResult,
    PeptideCandidate,
    Provenance,
)


class LabeledRecord(BaseModel):
    """One (candidate, optional fold, oracle label) triple."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    candidate: PeptideCandidate
    oracle_result: OracleResult
    complex_structure: ComplexStructure | None = None


class VersionedDataset(BaseModel):
    """In-memory + JSON-persistable labeled set (core campaign asset)."""

    model_config = ConfigDict(extra="forbid")

    data_version: str = "v0"
    records: list[LabeledRecord] = Field(default_factory=list)
    provenance: Provenance = Field(default_factory=Provenance)

    def add(self, record: LabeledRecord) -> None:
        self.records.append(record)

    def by_candidate_id(self) -> dict[UUID, LabeledRecord]:
        return {r.candidate.candidate_id: r for r in self.records}

    def candidate_ids(self) -> tuple[UUID, ...]:
        return tuple(r.candidate.candidate_id for r in self.records)

    def oracle_by_id(self) -> dict[UUID, OracleResult]:
        return {r.candidate.candidate_id: r.oracle_result for r in self.records}

    def values_by_id(self) -> dict[UUID, float]:
        return {r.candidate.candidate_id: r.oracle_result.value for r in self.records}

    def best_value(self, *, minimize: bool = True) -> float | None:
        if not self.records:
            return None
        vals = [r.oracle_result.value for r in self.records]
        return min(vals) if minimize else max(vals)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> VersionedDataset:
        return cls.model_validate_json(path.read_text(encoding="utf-8"))

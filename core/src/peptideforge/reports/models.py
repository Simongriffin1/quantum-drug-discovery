"""Traceable report models — every number cites a source (MLflow / JSON artifact)."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class TraceableNumber(BaseModel):
    """One reported figure with mandatory provenance."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    value: float | None
    unit: str | None = None
    source: str  # e.g. path, mlflow run_id, or "live:<function>"
    data_version: str | None = None
    notes: str | None = None


class ReportSection(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    title: str
    status: str  # PASS | FAIL | NOT_RUN | INFRA_ONLY
    summary: str
    numbers: tuple[TraceableNumber, ...] = ()
    details: dict[str, Any] = Field(default_factory=dict)


class BenchmarkReport(BaseModel):
    """Full PeptideForge credibility report (P12)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    title: str = "PeptideForge Benchmark Report"
    generated_at: str
    git_sha: str | None = None
    sections: tuple[ReportSection, ...]
    caveats: tuple[str, ...] = ()

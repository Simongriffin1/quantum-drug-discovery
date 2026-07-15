"""API schemas for campaign platform (P11)."""

from __future__ import annotations

from typing import Any
from uuid import UUID

from pydantic import BaseModel, Field


class StartCampaignRequest(BaseModel):
    goal: str = Field(..., min_length=1)
    seed: int = 0
    target_value: float = -4.0
    n_init: int = Field(8, ge=4, le=32)
    max_iterations: int = Field(3, ge=1, le=20)
    batch_size: int = Field(2, ge=1, le=8)
    simulation_mode: bool = True
    run_agent: bool = True


class ProvenanceOut(BaseModel):
    git_sha: str | None = None
    data_version: str | None = None
    tool_versions: dict[str, str] = Field(default_factory=dict)


class ParetoPoint(BaseModel):
    candidate_id: str
    sequence: str
    neg_binding: float
    solubility: float
    oracle_value: float


class CalibrationBinOut(BaseModel):
    bin_index: int
    n: int
    predicted_coverage: float
    empirical_coverage: float
    mean_interval_width: float


class CalibrationOut(BaseModel):
    n: int
    coverage_target: float
    empirical_coverage: float
    ece: float
    ece_threshold: float
    passed: bool
    reliability_bins: list[CalibrationBinOut]
    notes: str | None = None
    provenance: ProvenanceOut


class StructureOut(BaseModel):
    candidate_id: str
    sequence: str
    pdb_text: str
    fold_method: str
    confidence: float
    target_id: str
    provenance: ProvenanceOut


class TraceEventOut(BaseModel):
    kind: str
    content: str
    tool_name: str | None = None
    data: dict[str, Any] = Field(default_factory=dict)


class TraceOut(BaseModel):
    session_id: str
    events: list[TraceEventOut]
    summary: str | None = None
    provenance: ProvenanceOut


class IterationOut(BaseModel):
    iteration: int
    oracle_calls: int
    total_cost: float
    best_oracle_value: float | None
    acquisition_method: str
    status: str
    notes: str | None = None


class CampaignOut(BaseModel):
    campaign_id: UUID
    status: str
    goal: str
    simulation_mode: bool
    reached_target: bool
    best_value: float | None
    oracle_calls: int
    total_cost: float
    n_labeled: int
    iterations: list[IterationOut]
    provenance: ProvenanceOut
    agent_summary: str | None = None


class CampaignListOut(BaseModel):
    campaigns: list[CampaignOut]

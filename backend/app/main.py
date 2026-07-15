"""FastAPI application for PeptideForge campaigns (P11)."""

from __future__ import annotations

from uuid import UUID

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from app.campaigns import STORE, pareto_points, structure_for, to_campaign_out
from app.schemas import (
    CalibrationOut,
    CampaignListOut,
    CampaignOut,
    ParetoPoint,
    StartCampaignRequest,
    StructureOut,
    TraceOut,
)
from peptideforge import __version__ as core_version
from peptideforge.contracts.export_schemas import SCHEMA_MODELS
from pydantic import BaseModel, Field

app = FastAPI(
    title="PeptideForge API",
    description="Physics-grounded closed-loop peptide design platform",
    version=core_version,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://127.0.0.1:3000",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class HealthResponse(BaseModel):
    status: str = "ok"
    core_version: str
    contract_models: list[str] = Field(default_factory=list)


@app.get("/health", response_model=HealthResponse)
def health() -> HealthResponse:
    """Liveness probe — no fabricated campaign results."""
    return HealthResponse(
        status="ok",
        core_version=core_version,
        contract_models=sorted(SCHEMA_MODELS.keys()),
    )


@app.get("/contracts")
def list_contracts() -> dict[str, list[str]]:
    """Return registered contract model names (schemas only)."""
    return {"models": sorted(SCHEMA_MODELS.keys())}


@app.post("/campaigns", response_model=CampaignOut)
def start_campaign(body: StartCampaignRequest) -> CampaignOut:
    """Start a synthetic-mode DBTL campaign (+ optional agent trace)."""
    try:
        record = STORE.create(body)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return to_campaign_out(record)


@app.get("/campaigns", response_model=CampaignListOut)
def list_campaigns() -> CampaignListOut:
    return CampaignListOut(campaigns=[to_campaign_out(c) for c in STORE.list()])


@app.get("/campaigns/{campaign_id}", response_model=CampaignOut)
def get_campaign(campaign_id: UUID) -> CampaignOut:
    try:
        return to_campaign_out(STORE.get(campaign_id))
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/campaigns/{campaign_id}/pareto", response_model=list[ParetoPoint])
def get_pareto(campaign_id: UUID) -> list[ParetoPoint]:
    try:
        return pareto_points(STORE.get(campaign_id))
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/campaigns/{campaign_id}/structures/{candidate_id}", response_model=StructureOut)
def get_structure(campaign_id: UUID, candidate_id: str) -> StructureOut:
    try:
        return structure_for(STORE.get(campaign_id), candidate_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/campaigns/{campaign_id}/calibration", response_model=CalibrationOut)
def get_calibration(campaign_id: UUID) -> CalibrationOut:
    try:
        record = STORE.get(campaign_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    if record.calibration is None:
        raise HTTPException(status_code=404, detail="no calibration artifact")
    return record.calibration


@app.get("/campaigns/{campaign_id}/trace", response_model=TraceOut)
def get_trace(campaign_id: UUID) -> TraceOut:
    try:
        record = STORE.get(campaign_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    if record.agent_trace is None:
        raise HTTPException(status_code=404, detail="no agent trace (run_agent=false?)")
    return record.agent_trace

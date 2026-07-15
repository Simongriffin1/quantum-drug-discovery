"""Export Pydantic models to JSON Schema for API / cross-language consumers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from peptideforge.contracts.models import (
    AcquisitionBatch,
    CalibratedPrediction,
    Candidates,
    ComplexStructure,
    DevelopabilityScores,
    LoopState,
    ObjectiveVector,
    OracleResult,
    PeptideCandidate,
    Provenance,
)

SCHEMA_MODELS: dict[str, type[Any]] = {
    "Provenance": Provenance,
    "PeptideCandidate": PeptideCandidate,
    "Candidates": Candidates,
    "ComplexStructure": ComplexStructure,
    "OracleResult": OracleResult,
    "DevelopabilityScores": DevelopabilityScores,
    "CalibratedPrediction": CalibratedPrediction,
    "ObjectiveVector": ObjectiveVector,
    "AcquisitionBatch": AcquisitionBatch,
    "LoopState": LoopState,
}


def model_json_schema(name: str) -> dict[str, Any]:
    """Return JSON Schema for a named contract model."""
    if name not in SCHEMA_MODELS:
        raise KeyError(f"unknown schema model: {name}; known={sorted(SCHEMA_MODELS)}")
    schema: dict[str, Any] = SCHEMA_MODELS[name].model_json_schema()
    return schema


def export_all_schemas(output_dir: Path) -> list[Path]:
    """Write one JSON Schema file per contract model. Returns written paths."""
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name, model in SCHEMA_MODELS.items():
        path = output_dir / f"{name}.json"
        path.write_text(json.dumps(model.model_json_schema(), indent=2) + "\n", encoding="utf-8")
        written.append(path)
    return written


def main() -> None:
    """CLI: export schemas to peptideforge/contracts/schemas/."""
    here = Path(__file__).resolve().parent
    out = here / "schemas"
    paths = export_all_schemas(out)
    print(f"exported {len(paths)} schemas → {out}")


if __name__ == "__main__":
    main()

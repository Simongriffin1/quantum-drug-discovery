"""Tests for JSON Schema export."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from peptideforge.contracts.export_schemas import (
    SCHEMA_MODELS,
    export_all_schemas,
    model_json_schema,
)


@pytest.mark.contract
def test_all_schema_names_exportable() -> None:
    for name in SCHEMA_MODELS:
        schema = model_json_schema(name)
        assert "properties" in schema or "$defs" in schema or "title" in schema


@pytest.mark.contract
def test_unknown_schema_raises() -> None:
    with pytest.raises(KeyError):
        model_json_schema("NotARealModel")


@pytest.mark.contract
def test_export_writes_files(tmp_path: Path) -> None:
    written = export_all_schemas(tmp_path)
    assert len(written) == len(SCHEMA_MODELS)
    for path in written:
        data = json.loads(path.read_text(encoding="utf-8"))
        assert isinstance(data, dict)

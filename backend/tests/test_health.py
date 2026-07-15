"""Backend health contract tests (no campaign business logic yet)."""

from __future__ import annotations

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_health() -> None:
    response = client.get("/health")
    assert response.status_code == 200
    body = response.json()
    assert body["status"] == "ok"
    assert "LoopState" in body["contract_models"]


def test_contracts_list() -> None:
    response = client.get("/contracts")
    assert response.status_code == 200
    assert "OracleResult" in response.json()["models"]

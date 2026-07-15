"""Campaign lifecycle API tests (P11) — synthetic mode only."""

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


def test_campaign_lifecycle_synthetic() -> None:
    start = client.post(
        "/campaigns",
        json={
            "goal": "Find a strong binder in simulation mode",
            "seed": 0,
            "target_value": -4.0,
            "n_init": 8,
            "max_iterations": 2,
            "batch_size": 2,
            "simulation_mode": True,
            "run_agent": True,
        },
    )
    assert start.status_code == 200, start.text
    campaign = start.json()
    cid = campaign["campaign_id"]
    assert campaign["simulation_mode"] is True
    assert campaign["n_labeled"] >= 8
    assert campaign["provenance"]["data_version"] == "synthetic_v1"
    assert "peptideforge" in campaign["provenance"]["tool_versions"]
    assert campaign["iterations"]

    got = client.get(f"/campaigns/{cid}")
    assert got.status_code == 200
    assert got.json()["campaign_id"] == cid

    pareto = client.get(f"/campaigns/{cid}/pareto")
    assert pareto.status_code == 200
    points = pareto.json()
    assert len(points) >= 1
    assert "sequence" in points[0]
    assert "neg_binding" in points[0]

    structure = client.get(f"/campaigns/{cid}/structures/{points[0]['candidate_id']}")
    assert structure.status_code == 200
    body = structure.json()
    assert "ATOM" in body["pdb_text"] or "HEADER" in body["pdb_text"]
    assert body["provenance"]["data_version"] is not None or body["fold_method"]

    cal = client.get(f"/campaigns/{cid}/calibration")
    assert cal.status_code == 200
    assert "ece" in cal.json()
    assert cal.json()["provenance"]["tool_versions"]

    trace = client.get(f"/campaigns/{cid}/trace")
    assert trace.status_code == 200
    events = trace.json()["events"]
    assert any(e["kind"] == "gate_pause" for e in events)
    assert trace.json()["provenance"]["data_version"] == "synthetic_v1"


def test_reject_non_simulation() -> None:
    response = client.post(
        "/campaigns",
        json={
            "goal": "live",
            "simulation_mode": False,
            "run_agent": False,
        },
    )
    assert response.status_code == 400


def test_unknown_campaign_404() -> None:
    response = client.get("/campaigns/00000000-0000-0000-0000-000000000000")
    assert response.status_code == 404

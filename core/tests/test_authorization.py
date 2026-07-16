"""Tests for Step 4 authorization + claim ceiling."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from peptideforge.authorization import (
    AuthorizationDenied,
    AuthorizationRecord,
    InputType,
    TaskType,
    assert_campaign_authorized,
    build_authorization_bundle,
    filter_by_fold_confidence,
    write_authorization_bundle,
    load_authorization_bundle,
)
from peptideforge.claims import (
    PointAffinityClaimError,
    assert_no_point_affinity_claim,
    claim_ceiling_from_authorization,
    enrich_interval_width,
)


def test_authorization_roundtrip(tmp_path: Path) -> None:
    exp = {
        "gate_pass": True,
        "spearman": 0.381,
        "spearman_ci_low": 0.189,
        "spearman_ci_high": 0.556,
        "n": 100,
        "artifact": "skempi_powered.json",
    }
    predicted = {
        "status": "BLOCKED_BOLTZ_UNAVAILABLE",
        "mode_a": {"gate_pass": False, "n": 0},
        "block_reason": "Boltz unavailable",
        "artifact": "fold_degradation.json",
    }
    records = build_authorization_bundle(
        experimental_skempi=exp,
        predicted_degradation=predicted,
        split_id="skempi_powered_holdout_v1",
    )
    path = tmp_path / "auth.json"
    write_authorization_bundle(records, path)
    loaded = load_authorization_bundle(path)
    assert len(loaded) == len(records)
    exp_rec = next(
        r
        for r in loaded
        if r.task_type == TaskType.WITHIN_TARGET
        and r.input_type == InputType.EXPERIMENTAL
    )
    assert exp_rec.authorized is True
    pred_rec = next(
        r
        for r in loaded
        if r.task_type == TaskType.WITHIN_TARGET and r.input_type == InputType.PREDICTED
    )
    assert pred_rec.authorized is False


def test_reject_cross_target_and_predicted_when_blocked(tmp_path: Path) -> None:
    records = build_authorization_bundle(
        experimental_skempi={
            "gate_pass": True,
            "spearman": 0.381,
            "spearman_ci_low": 0.19,
            "spearman_ci_high": 0.55,
            "n": 100,
        },
        predicted_degradation={
            "mode_a": {"gate_pass": False, "n": 0},
            "block_reason": "blocked",
        },
        split_id="test",
    )
    with pytest.raises(AuthorizationDenied):
        assert_campaign_authorized(
            records,
            task_type=TaskType.CROSS_TARGET,
            input_type=InputType.EXPERIMENTAL,
        )
    with pytest.raises(AuthorizationDenied):
        assert_campaign_authorized(
            records,
            task_type=TaskType.WITHIN_TARGET,
            input_type=InputType.PREDICTED,
        )
    ok = assert_campaign_authorized(
        records,
        task_type=TaskType.WITHIN_TARGET,
        input_type=InputType.EXPERIMENTAL,
    )
    assert ok.authorized


def test_simulation_always_allowed() -> None:
    rec = assert_campaign_authorized(
        [],
        task_type=TaskType.WITHIN_TARGET,
        input_type=InputType.PREDICTED,
        simulation_mode=True,
    )
    assert rec.input_type == InputType.SIMULATION


def test_confidence_filter() -> None:
    kept, excluded = filter_by_fold_confidence(
        {"a": 0.9, "b": 0.5, "c": 0.8}, threshold=0.8
    )
    assert set(kept) == {"a", "c"}
    assert excluded == ["b"]


def test_claim_ceiling_no_point_affinity() -> None:
    rec = AuthorizationRecord(
        task_type=TaskType.WITHIN_TARGET,
        input_type=InputType.EXPERIMENTAL,
        structure_source="crystal",
        validated_rho=0.38,
        authorized=True,
        reason="test",
    )
    ceiling = claim_ceiling_from_authorization(rec)
    assert ceiling.allow_point_affinity is False
    assert "enrichment" in ceiling.max_claim.lower()
    with pytest.raises(PointAffinityClaimError):
        assert_no_point_affinity_claim(
            text="The predicted affinity is 5 nM",
            task_type=TaskType.WITHIN_TARGET,
            input_type=InputType.EXPERIMENTAL,
            authorization=rec,
        )
    assert_no_point_affinity_claim(
        text="Top candidates ranked for enrichment over random",
        task_type=TaskType.WITHIN_TARGET,
        input_type=InputType.EXPERIMENTAL,
        authorization=rec,
    )
    wide = enrich_interval_width(1.0, 0.38)
    assert wide > 1.5

"""Tests for bootstrap CIs and extended affinity validity gate."""

from __future__ import annotations

import pytest

from peptideforge.eval.affinity_validity import (
    evaluate_affinity_with_ci,
    run_trivial_baseline_battery,
)
from peptideforge.eval.harness import PredictionLabelPair
from peptideforge.eval.metrics import bootstrap_ci, spearman_rho
from peptideforge.eval.redteam import leakage_audit, run_red_team


@pytest.mark.eval
def test_bootstrap_ci_width_reported() -> None:
    x = [float(i) for i in range(30)]
    y = [0.8 * v + 0.1 for v in x]
    point, lo, hi = bootstrap_ci(x, y, statistic="spearman", n_resamples=200, seed=0)
    assert point == pytest.approx(spearman_rho(x, y), abs=1e-9)
    assert lo <= point <= hi
    assert hi - lo >= 0.0


@pytest.mark.eval
def test_gate_requires_min_n_and_ci() -> None:
    pairs = [
        PredictionLabelPair(record_id=f"r{i}", predicted=float(i), experimental=float(i))
        for i in range(10)
    ]
    report = evaluate_affinity_with_ci(pairs, min_n=30, n_bootstrap=200, seed=0)
    assert report.measurable is False
    assert report.passed is False
    assert report.notes is not None


@pytest.mark.eval
def test_gate_pass_on_strong_correlation_large_n() -> None:
    pairs = [
        PredictionLabelPair(
            record_id=f"r{i}", predicted=float(i) + 0.01, experimental=float(i)
        )
        for i in range(40)
    ]
    report = evaluate_affinity_with_ci(pairs, min_n=30, n_bootstrap=300, seed=0)
    assert report.measurable
    assert report.spearman >= 0.40
    assert report.spearman_ci_low > 0
    assert report.passed


@pytest.mark.eval
def test_trivial_baseline_battery_halts_when_baseline_wins() -> None:
    # Weak/anti-correlated oracle vs a strongly ranked length baseline
    pairs = [
        PredictionLabelPair(
            record_id=f"r{i}", predicted=float(12 - i), experimental=float(i)
        )
        for i in range(12)
    ]
    baselines = {"peptide_length": [float(i) for i in range(12)]}
    ok, rhos, msg = run_trivial_baseline_battery(pairs, baselines)
    assert ok is False
    assert msg is not None
    assert "peptide_length" in (msg or "")


@pytest.mark.eval
def test_leakage_audit_flags_high_identity() -> None:
    ok, findings = leakage_audit(
        ["a", "b"],
        ["c"],
        train_sequences={"a": "ACDEFGHIKL", "b": "XXXXXXXXXX"},
        test_sequences={"c": "ACDEFGHIKM"},
        max_identity=0.30,
    )
    assert ok is False
    assert any("identity" in f for f in findings)


@pytest.mark.eval
def test_deliberate_leak_flagged_by_red_team() -> None:
    pairs = [
        PredictionLabelPair(record_id="t1", predicted=1.0, experimental=1.0),
        PredictionLabelPair(record_id="t2", predicted=2.0, experimental=2.0),
        PredictionLabelPair(record_id="t3", predicted=3.0, experimental=3.0),
    ]
    report = run_red_team(
        pairs,
        train_ids=["tr1"],
        test_ids=["t1"],
        train_clusters={"tr1": "c0"},
        test_clusters={"t1": "c0"},
        train_sequences={"tr1": "ACDE"},
        test_sequences={"t1": "ACDE"},
    )
    assert report.leakage_passed is False
    assert report.passed is False

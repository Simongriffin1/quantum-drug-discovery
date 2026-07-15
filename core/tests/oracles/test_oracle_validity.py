"""Oracle-validity integration: real Spearman + MLflow log (slow)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from peptideforge.eval.harness import DEFAULT_THRESHOLDS
from peptideforge.oracles.openmm_utils import OpenMMUnavailableError, require_openmm
from peptideforge.oracles.validate_affinity import SUBSET_NAME, run_affinity_validity

STRUCTURES = Path(__file__).resolve().parents[3] / "benchmarks" / "fixtures" / "structures"


@pytest.mark.slow
@pytest.mark.golden
def test_oracle_validity_spearman_logged(tmp_path: Path) -> None:
    """Produce a real Spearman on the named interface subset and log to MLflow.

    Does NOT assert the stage-gate threshold — gate status is reported honestly.
    Threshold is pre-registered in ACCEPTANCE.md (ρ ≥ 0.40).
    """
    try:
        require_openmm()
    except OpenMMUnavailableError as exc:
        pytest.skip(str(exc))

    manifest = STRUCTURES / "structure_manifest_v1.tsv"
    if not manifest.is_file():
        pytest.skip("structure manifest missing — run prepare_interfaces_from_full")

    report, extras = run_affinity_validity(
        structures_dir=STRUCTURES,
        minimize_max_iterations=0,
        log_mlflow=True,
        mlflow_tracking_uri=f"sqlite:///{(tmp_path / 'oracle_validity.db').resolve()}",
        seed=0,
        platform="CPU",
    )
    assert extras["subset_name"] == SUBSET_NAME
    assert report.n >= 3
    assert report.spearman == report.spearman  # finite
    assert "mlflow" in extras and extras["mlflow"]["run_id"]
    assert report.spearman_threshold == DEFAULT_THRESHOLDS.oracle_affinity_spearman

    out = STRUCTURES / "oracle_validity_last_run.json"
    out.write_text(
        json.dumps(
            {
                "subset_name": extras["subset_name"],
                "spearman": report.spearman,
                "rmse": report.rmse,
                "n": report.n,
                "passed_threshold": report.passed_spearman_threshold,
                "spearman_threshold": report.spearman_threshold,
                "details": extras["details"],
                "mlflow": extras["mlflow"],
                "protocol": extras["protocol"],
                "failures": extras["failures"],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        f"\nORACLE-VALIDITY: Spearman={report.spearman:.4f} N={report.n} "
        f"passed_gate={report.passed_spearman_threshold} → {out}"
    )

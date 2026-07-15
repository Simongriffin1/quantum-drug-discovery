"""Oracle-validity: MM-GBSA vs experimental pK on a named structure subset.

Produces a **real** Spearman (never fabricated). Logs to MLflow when requested.
Threshold comparison uses pre-registered ACCEPTANCE.md values — gate may FAIL;
that result is reported honestly.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any
from uuid import uuid4

from peptideforge_benchmarks.pdbbind import load_pdbbind_peptide_affinity

from peptideforge.contracts.models import ComplexStructure, OracleTier, PeptideCandidate
from peptideforge.eval.harness import (
    DEFAULT_THRESHOLDS,
    PredictionLabelPair,
    ValidityReport,
    evaluate_predictions,
)
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle

# Named subset shown in reports — do not rename without updating ACCEPTANCE refs.
SUBSET_NAME = "pdbbind_peptide_affinity_v1_structures_openmm"


def _require_mlflow() -> Any:
    try:
        import mlflow
    except ImportError as exc:
        raise ImportError(
            "mlflow is required for oracle-validity logging but is not installed. "
            "Install with: pip install mlflow. Refusing to silently skip logging."
        ) from exc
    return mlflow


def load_structure_manifest(path: Path) -> list[dict[str, str]]:
    """Load structure manifest TSV (supports pdb_file or pdb_path columns)."""
    if not path.is_file():
        raise FileNotFoundError(f"structure manifest not found: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "record_id" not in reader.fieldnames:
            raise ValueError(f"invalid structure manifest: {path}")
        rows = [
            row for row in reader if row.get("record_id") and not row["record_id"].startswith("#")
        ]
    if not rows:
        raise ValueError(f"empty structure manifest: {path}")
    return rows


def _pdb_filename(row: dict[str, str]) -> str:
    name = row.get("pdb_file") or row.get("pdb_path")
    if not name:
        raise KeyError(f"manifest row {row.get('record_id')} missing pdb_file/pdb_path")
    return name


def _experimental_pk(row: dict[str, str], pk_by_id: dict[str, float]) -> float:
    if row.get("experimental_pk"):
        return float(row["experimental_pk"])
    rid = row["record_id"]
    if rid not in pk_by_id:
        raise KeyError(f"no experimental_pk for {rid}; add column or affinity fixture record")
    return pk_by_id[rid]


def _peptide_sequence(row: dict[str, str], seq_by_id: dict[str, str]) -> str:
    seq = row.get("peptide_sequence") or row.get("epitope") or seq_by_id.get(row["record_id"])
    if not seq:
        raise KeyError(f"no peptide sequence for {row['record_id']}")
    seq = seq.upper()
    if len(seq) < 5:
        seq = (seq + "AAAAA")[:5]
    return seq[:50]


def run_affinity_validity(
    *,
    structures_dir: Path,
    manifest_name: str = "structure_manifest_v1.tsv",
    tier: OracleTier = OracleTier.MM_GBSA,
    minimize_max_iterations: int = 100,
    md_steps: int = 0,
    seed: int = 0,
    platform: str | None = "CPU",
    log_mlflow: bool = False,
    mlflow_tracking_uri: str | None = None,
    mlflow_experiment: str = "peptideforge-oracle-validity",
) -> tuple[ValidityReport, dict[str, Any]]:
    """Evaluate OpenMM oracle on every row in the structure manifest.

    Returns (ValidityReport, extras). Raises if OpenMM missing or <2 successes.
    """
    manifest_path = structures_dir / manifest_name
    rows = load_structure_manifest(manifest_path)

    affinity = load_pdbbind_peptide_affinity()
    pk_by_id = {r.record_id: r.pk for r in affinity}
    seq_by_id = {r.record_id: r.peptide_sequence for r in affinity}

    pairs: list[PredictionLabelPair] = []
    details: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []

    for row in rows:
        pdb_file = structures_dir / _pdb_filename(row)
        if not pdb_file.is_file():
            failures.append({"record_id": row["record_id"], "error": f"missing {pdb_file}"})
            continue
        try:
            experimental_pk = _experimental_pk(row, pk_by_id)
            peptide_seq = _peptide_sequence(row, seq_by_id)
        except KeyError as exc:
            failures.append({"record_id": row["record_id"], "error": str(exc)})
            continue

        peptide_chain = row.get("peptide_chain") or "C"
        cand = PeptideCandidate(
            candidate_id=uuid4(),
            sequence=peptide_seq,
            generation_method="benchmark_fixture",
            metadata={"record_id": row["record_id"], "peptide_chain": peptide_chain},
        )
        complex_structure = ComplexStructure(
            candidate_id=cand.candidate_id,
            target_id=row.get("pdb_id") or row["record_id"],
            sequence=cand.sequence,
            pdb_path=str(pdb_file.resolve()),
            confidence=1.0,
            fold_method="experimental_pdb_fixture",
        )
        try:
            oracle = OpenMMPhysicsOracle(
                OpenMMOracleConfig(
                    minimize_max_iterations=minimize_max_iterations,
                    docking_minimize_iterations=min(40, minimize_max_iterations),
                    md_steps=md_steps if md_steps > 0 else 500,
                    seed=seed,
                    platform=platform,
                    peptide_chain_ids=(peptide_chain, "P", "L", "C"),
                )
            )
            result = oracle.evaluate(complex_structure, tier=tier)
        except Exception as exc:  # noqa: BLE001 — collect; fail if none succeed
            failures.append({"record_id": row["record_id"], "error": str(exc)})
            continue

        predicted = -float(result.value)
        pairs.append(
            PredictionLabelPair(
                record_id=row["record_id"],
                predicted=predicted,
                experimental=experimental_pk,
                unit="pK_proxy_neg_dG",
            )
        )
        details.append(
            {
                "record_id": row["record_id"],
                "pdb_file": _pdb_filename(row),
                "experimental_pk": experimental_pk,
                "mm_gbsa_kcal_mol": result.value,
                "predicted_neg_dG": predicted,
                "tier": result.tier.value,
                "cost_estimate": result.cost_estimate,
            }
        )

    if len(pairs) < 2:
        raise RuntimeError(
            f"oracle-validity needs ≥2 successful evaluations; got {len(pairs)}. "
            f"failures={failures}"
        )

    report = evaluate_predictions(
        pairs,
        spearman_threshold=DEFAULT_THRESHOLDS.oracle_affinity_spearman,
        metric_target="affinity_pK",
    )

    extras: dict[str, Any] = {
        "subset_name": SUBSET_NAME,
        "tier": tier.value,
        "n_manifest": len(rows),
        "n_success": len(pairs),
        "n_failed": len(failures),
        "failures": failures,
        "details": details,
        "protocol": {
            "method": "OpenMM MM-GBSA (amber14 + OBC2), ΔG ≈ E_c − E_r − E_l after minimize",
            "invert_sign": True,
            "minimize_max_iterations": minimize_max_iterations,
            "platform": platform,
            "seed": seed,
            "manifest": str(manifest_path),
            "structures": "trimmed public MHC–peptide interface PDBs (see REMARK headers)",
        },
    }

    if log_mlflow:
        mlflow = _require_mlflow()
        if not mlflow_tracking_uri:
            raise ValueError(
                "log_mlflow=True requires mlflow_tracking_uri "
                "(prefer sqlite:///…/oracle_validity.db — file store is maintenance-mode)"
            )
        mlflow.set_tracking_uri(mlflow_tracking_uri)
        mlflow.set_experiment(mlflow_experiment)
        with mlflow.start_run(run_name=f"{SUBSET_NAME}-{tier.value}") as run:
            mlflow.log_param("subset_name", SUBSET_NAME)
            mlflow.log_param("tier", tier.value)
            mlflow.log_param("minimize_max_iterations", minimize_max_iterations)
            mlflow.log_param("platform", platform or "default")
            mlflow.log_param("seed", seed)
            mlflow.log_metric("spearman", report.spearman)
            mlflow.log_metric("rmse", report.rmse)
            mlflow.log_metric("n", float(report.n))
            mlflow.log_metric("passed_spearman_threshold", float(report.passed_spearman_threshold))
            mlflow.log_metric("spearman_threshold", report.spearman_threshold)
            extras["mlflow"] = {
                "run_id": run.info.run_id,
                "experiment_id": run.info.experiment_id,
                "tracking_uri": mlflow.get_tracking_uri(),
            }

    return report, extras

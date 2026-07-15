"""CLI: MM-GBSA oracle-validity on prepared experimental structures with CIs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from peptideforge.eval.affinity_validity import (
    EXTENDED_THRESHOLDS,
    SUBSET_NAME_V2,
    charge_from_sequence,
    evaluate_affinity_with_ci,
    report_to_dict,
    run_trivial_baseline_battery,
    score_prepared_structures,
    trivial_baselines_for_structures,
    write_validity_artifact,
)
from peptideforge.eval.redteam import run_red_team
from peptideforge.oracles.validate_affinity import _require_mlflow


def _load_prep_manifest(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("ok") in {"1", "true", "True"} and row.get("complex_path"):
                rows.append(row)
    if not rows:
        raise ValueError(f"no successful prep rows in {path}")
    # Prefer scoreable (water-stripped) paths when present
    for row in rows:
        if row.get("scoreable_path"):
            row["complex_path"] = row["scoreable_path"]
    return rows


def _load_catalog(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {r["record_id"]: r for r in csv.DictReader(handle, delimiter="\t")}


def _load_splits(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--prep-manifest",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--splits",
        type=Path,
        required=True,
        help="splits_v2.json from peptide_affinity.splits (held-out test used ONCE)",
    )
    parser.add_argument(
        "--partition",
        choices=("test", "train", "val", "all"),
        default="test",
        help="Which split to score. Tuning must NEVER use test.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json"),
    )
    parser.add_argument("--minimize-max-iterations", type=int, default=0)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--mlflow", action="store_true")
    parser.add_argument("--mlflow-uri", type=str, default=None)
    args = parser.parse_args()

    catalog = _load_catalog(args.catalog)
    prep_rows = _load_prep_manifest(args.prep_manifest)
    splits = _load_splits(args.splits)
    cold = splits["cold_start"]
    cluster_of = splits["clusters"]

    if args.partition == "test":
        keep = set(cold["test_ids"])
    elif args.partition == "train":
        keep = set(cold["train_ids"])
    elif args.partition == "val":
        keep = set(cold["val_ids"])
    else:
        keep = set(catalog)

    manifest_rows: list[dict[str, str]] = []
    for row in prep_rows:
        rid = row["record_id"]
        if rid not in keep or rid not in catalog:
            continue
        entry = catalog[rid]
        manifest_rows.append(
            {
                "record_id": rid,
                "pdb_id": entry["pdb_id"],
                "pdb_path": row["complex_path"],
                "peptide_chain": entry.get("peptide_chain") or "C",
                "peptide_sequence": entry["peptide_seq"],
                "experimental_pk": entry["pKd"],
            }
        )

    if len(manifest_rows) < 3:
        raise SystemExit(
            f"only {len(manifest_rows)} prepared structures in partition={args.partition}"
        )

    pk_by_id = {r["record_id"]: float(r["experimental_pk"]) for r in manifest_rows}
    seq_by_id = {r["record_id"]: r["peptide_sequence"] for r in manifest_rows}
    pairs, details, failures = score_prepared_structures(
        manifest_rows,
        pk_by_id=pk_by_id,
        seq_by_id=seq_by_id,
        minimize_max_iterations=args.minimize_max_iterations,
        seed=args.seed,
        platform=args.platform,
    )
    if len(pairs) < 3:
        raise SystemExit(f"scoring failed for nearly all rows: {failures[:5]}")

    lengths = {rid: len(seq_by_id[rid]) for rid in seq_by_id}
    charges = {rid: charge_from_sequence(seq_by_id[rid]) for rid in seq_by_id}
    baselines = trivial_baselines_for_structures(
        pairs, peptide_lengths=lengths, net_charges=charges
    )
    base_ok, baseline_rhos, base_msg = run_trivial_baseline_battery(pairs, baselines)
    if not base_ok:
        # Halt and report — do not declare PASS
        print(f"HALT: trivial baseline battery failed: {base_msg}")

    # Leakage audit train vs this partition when scoring test
    train_ids = tuple(cold["train_ids"])
    test_ids = tuple(p.record_id for p in pairs)
    train_seqs = {
        rid: catalog[rid]["peptide_seq"] for rid in train_ids if rid in catalog
    }
    test_seqs = {p.record_id: seq_by_id[p.record_id] for p in pairs}
    train_clusters = {rid: cluster_of[rid] for rid in train_ids if rid in cluster_of}
    test_clusters = {rid: cluster_of[rid] for rid in test_ids if rid in cluster_of}
    red = run_red_team(
        pairs,
        train_ids=train_ids if args.partition == "test" else None,
        test_ids=test_ids if args.partition == "test" else None,
        train_clusters=train_clusters if args.partition == "test" else None,
        test_clusters=test_clusters if args.partition == "test" else None,
        train_sequences=train_seqs if args.partition == "test" else None,
        test_sequences=test_seqs if args.partition == "test" else None,
        baseline_predicted=baselines["peptide_length"],
        max_identity=0.30 if args.partition == "test" else None,
        seed=args.seed,
    )
    # Force red-team fail if trivial battery failed
    if not base_ok:
        from peptideforge.eval.redteam import RedTeamReport

        red = RedTeamReport(
            passed=False,
            label_shuffle_passed=red.label_shuffle_passed,
            label_shuffle_rho=red.label_shuffle_rho,
            trivial_baseline_passed=False,
            model_rho=red.model_rho,
            baseline_rho=red.baseline_rho,
            leakage_passed=red.leakage_passed,
            leakage_findings=red.leakage_findings + (base_msg or "baseline_fail",),
            notes=base_msg,
        )

    report = evaluate_affinity_with_ci(
        pairs,
        thresholds=EXTENDED_THRESHOLDS,
        min_n=EXTENDED_THRESHOLDS.min_n,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        red_team=red,
    )

    extras: dict[str, Any] = {
        "subset_name": SUBSET_NAME_V2,
        "data_version": SUBSET_NAME_V2,
        "partition": args.partition,
        "n_manifest": len(manifest_rows),
        "n_success": len(pairs),
        "n_failed": len(failures),
        "failures": failures,
        "details": details,
        "baseline_rhos": baseline_rhos,
        "protocol": {
            "method": "OpenMM MM-GBSA on prepared experimental structures",
            "minimize_max_iterations": args.minimize_max_iterations,
            "platform": args.platform,
            "seed": args.seed,
            "n_bootstrap": args.n_bootstrap,
            "catalog": str(args.catalog),
            "prep_manifest": str(args.prep_manifest),
            "splits": str(args.splits),
        },
    }

    if args.mlflow:
        mlflow = _require_mlflow()
        uri = args.mlflow_uri
        if not uri:
            raise SystemExit("--mlflow requires --mlflow-uri")
        mlflow.set_tracking_uri(uri)
        mlflow.set_experiment("peptideforge-oracle-validity-v2")
        with mlflow.start_run(run_name=f"{SUBSET_NAME_V2}-{args.partition}") as run:
            mlflow.log_param("subset_name", SUBSET_NAME_V2)
            mlflow.log_param("partition", args.partition)
            mlflow.log_param("data_version", SUBSET_NAME_V2)
            mlflow.log_metric("spearman", report.spearman)
            mlflow.log_metric("spearman_ci_low", report.spearman_ci_low)
            mlflow.log_metric("spearman_ci_high", report.spearman_ci_high)
            mlflow.log_metric("pearson", report.pearson)
            mlflow.log_metric("n", float(report.n))
            mlflow.log_metric("passed", float(report.passed))
            extras["mlflow"] = {
                "run_id": run.info.run_id,
                "experiment_id": run.info.experiment_id,
                "tracking_uri": mlflow.get_tracking_uri(),
            }

    payload = report_to_dict(report, extras)
    write_validity_artifact(args.out, payload)
    print(
        f"partition={args.partition} N={report.n} Spearman={report.spearman:.4f} "
        f"CI=[{report.spearman_ci_low:.4f},{report.spearman_ci_high:.4f}] "
        f"measurable={report.measurable} PASSED={report.passed}"
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

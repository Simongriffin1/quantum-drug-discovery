"""CLI: run oracle-validity and print / log Spearman (never fabricated)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from peptideforge.oracles.validate_affinity import SUBSET_NAME, run_affinity_validity


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structures-dir",
        type=Path,
        default=Path(__file__).resolve().parents[4] / "benchmarks" / "fixtures" / "structures",
    )
    parser.add_argument("--minimize", type=int, default=0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--platform", type=str, default="CPU")
    parser.add_argument("--mlflow", action="store_true")
    parser.add_argument(
        "--mlflow-uri",
        type=str,
        default=None,
        help="MLflow tracking URI (default: sqlite:///…/oracle_validity_mlflow.db)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Write JSON report path (default: structures/oracle_validity_last_run.json)",
    )
    args = parser.parse_args()

    out = args.out or (args.structures_dir / "oracle_validity_last_run.json")
    mlflow_uri = args.mlflow_uri
    if args.mlflow and mlflow_uri is None:
        mlflow_uri = f"sqlite:///{(args.structures_dir / 'oracle_validity_mlflow.db').resolve()}"

    report, extras = run_affinity_validity(
        structures_dir=args.structures_dir,
        minimize_max_iterations=args.minimize,
        seed=args.seed,
        platform=args.platform,
        log_mlflow=args.mlflow,
        mlflow_tracking_uri=mlflow_uri,
    )
    payload = {
        "subset_name": SUBSET_NAME,
        "spearman": report.spearman,
        "rmse": report.rmse,
        "n": report.n,
        "passed_threshold": report.passed_spearman_threshold,
        "spearman_threshold": report.spearman_threshold,
        "details": extras["details"],
        "failures": extras["failures"],
        "protocol": extras["protocol"],
        "mlflow": extras.get("mlflow"),
    }
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        f"ORACLE-VALIDITY subset={SUBSET_NAME} N={report.n} "
        f"Spearman={report.spearman:.4f} RMSE={report.rmse:.4f} "
        f"passed_gate={report.passed_spearman_threshold} → {out}"
    )


if __name__ == "__main__":
    main()

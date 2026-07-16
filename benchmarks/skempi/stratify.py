"""Step 4A.4 — confidence-stratified recovery of predicted-fold within-target ρ."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from peptideforge.eval.affinity_validity import evaluate_affinity_with_ci
from peptideforge.eval.harness import PredictionLabelPair


def stratified_rho_curve(
    pairs: list[PredictionLabelPair],
    quality_by_id: dict[str, float],
    *,
    thresholds: tuple[float, ...] = (0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    min_surviving_n: int = 30,
    seed: int = 0,
) -> dict[str, Any]:
    """ρ vs confidence threshold; identify min threshold that recovers the gate."""
    curve: list[dict[str, Any]] = []
    best_auth: dict[str, Any] | None = None
    for thr in thresholds:
        keep = [p for p in pairs if quality_by_id.get(p.record_id, -1.0) >= thr]
        if len(keep) < 3:
            curve.append(
                {
                    "threshold": thr,
                    "n": len(keep),
                    "status": "TOO_FEW",
                    "gate_pass": False,
                }
            )
            continue
        report = evaluate_affinity_with_ci(
            keep, min_n=min(min_surviving_n, 3), n_bootstrap=800, seed=seed
        )
        gate = bool(
            report.spearman >= 0.30
            and report.spearman_ci_low > 0.0
            and len(keep) >= min_surviving_n
        )
        row = {
            "threshold": thr,
            "n": len(keep),
            "spearman": report.spearman,
            "spearman_ci_low": report.spearman_ci_low,
            "spearman_ci_high": report.spearman_ci_high,
            "gate_pass": gate,
            "meets_n30": len(keep) >= min_surviving_n,
        }
        curve.append(row)
        if gate and best_auth is None:
            best_auth = {
                "metric": "fold_confidence",
                "threshold": thr,
                "gate_pass": True,
                **{k: row[k] for k in ("n", "spearman", "spearman_ci_low", "spearman_ci_high")},
            }
    return {
        "curve": curve,
        "stratified_authorization": best_auth
        or {
            "gate_pass": False,
            "reason": "no threshold recovers ρ≥0.30 with CI_low>0 and N≥30",
        },
        "min_surviving_n": min_surviving_n,
    }


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--degradation",
        type=Path,
        default=Path("benchmarks/skempi/data/fold_degradation_last_run.json"),
    )
    parser.add_argument(
        "--quality",
        type=Path,
        default=Path("benchmarks/skempi/data/predicted_folds/fold_quality.json"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/skempi/data/fold_stratify_last_run.json"),
    )
    args = parser.parse_args()
    if not args.degradation.is_file():
        raise SystemExit(f"missing {args.degradation}")
    deg = json.loads(args.degradation.read_text(encoding="utf-8"))
    if deg.get("status") in {"BLOCKED_BOLTZ_UNAVAILABLE", "BLOCKED"}:
        payload = {
            "status": deg.get("status"),
            "curve": [],
            "stratified_authorization": {
                "gate_pass": False,
                "reason": "no predicted scores — Boltz unavailable or blocked",
            },
        }
        args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(json.dumps(payload["stratified_authorization"], indent=2))
        return

    quality = {}
    if args.quality.is_file():
        q = json.loads(args.quality.read_text(encoding="utf-8"))
        for row in q.get("rows") or []:
            # Prefer dockq_proxy, else confidence
            val = row.get("dockq_proxy")
            if val is None:
                val = row.get("confidence")
            if val is not None:
                quality[row["record_id"]] = float(val)

    details = (deg.get("mode_a") or {}).get("details") or []
    pairs = [
        PredictionLabelPair(
            record_id=d["record_id"],
            predicted=float(d["predicted_ddg"]),
            experimental=float(d["experimental_ddg"]),
            unit="kcal/mol",
        )
        for d in details
        if "predicted_ddg" in d
    ]
    payload = stratified_rho_curve(pairs, quality)
    payload["status"] = "OK"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload["stratified_authorization"], indent=2))


if __name__ == "__main__":
    main()

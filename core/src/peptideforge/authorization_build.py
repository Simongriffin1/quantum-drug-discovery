"""Build and install the Step 4B authorization bundle from measured artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from peptideforge.authorization import (
    build_authorization_bundle,
    write_authorization_bundle,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--experimental",
        type=Path,
        default=Path("benchmarks/skempi/data/skempi_ddg_powered_last_run.json"),
    )
    parser.add_argument(
        "--degradation",
        type=Path,
        default=Path("benchmarks/skempi/data/fold_degradation_last_run.json"),
    )
    parser.add_argument(
        "--stratify",
        type=Path,
        default=Path("benchmarks/skempi/data/fold_stratify_last_run.json"),
    )
    parser.add_argument(
        "--holdout",
        type=Path,
        default=Path("benchmarks/skempi/data/skempi_powered_holdout_v1.json"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/authorization/authorization_bundle.json"),
    )
    args = parser.parse_args()

    exp = json.loads(args.experimental.read_text(encoding="utf-8"))
    exp["artifact"] = str(args.experimental)
    holdout = json.loads(args.holdout.read_text(encoding="utf-8"))
    split_id = holdout.get("split_id", "skempi_powered_holdout_v1")

    predicted = None
    if args.degradation.is_file():
        predicted = json.loads(args.degradation.read_text(encoding="utf-8"))
        if args.stratify.is_file():
            strat = json.loads(args.stratify.read_text(encoding="utf-8"))
            predicted["stratified_authorization"] = strat.get("stratified_authorization")
        if predicted.get("status") in {"BLOCKED_BOLTZ_UNAVAILABLE", "BLOCKED"}:
            # Keep structure for build_authorization_bundle — mode_a gate_pass false
            predicted.setdefault("mode_a", {"gate_pass": False, "n": 0})
            predicted["block_reason"] = predicted.get("block_reason") or predicted.get(
                "verdict"
            )

    records = build_authorization_bundle(
        experimental_skempi=exp,
        predicted_degradation=predicted,
        split_id=split_id,
    )
    write_authorization_bundle(records, args.out)
    summary = [
        {
            "task": r.task_type.value,
            "input": r.input_type.value,
            "authorized": r.authorized,
            "reason": r.reason[:120],
            "record_id": str(r.record_id),
        }
        for r in records
    ]
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

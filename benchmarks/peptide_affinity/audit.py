"""Attrition & selection-bias audit for the peptide-affinity funnel.

Quantifies drop reasons and whether survivors are charge-/length-skewed
relative to the PepBench source pool — resolving whether net_charge winning
is a selection artifact or a genuine oracle failure.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


def charge_from_sequence(seq: str) -> float:
    return float(
        sum(1 for c in seq.upper() if c in "KR") - sum(1 for c in seq.upper() if c in "DE")
    )

def _mean(xs: list[float]) -> float:
    return sum(xs) / len(xs) if xs else float("nan")


def _std(xs: list[float]) -> float:
    if len(xs) < 2:
        return float("nan")
    m = _mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def _summarize(xs: list[float]) -> dict[str, float]:
    if not xs:
        return {"n": 0.0}
    xs_sorted = sorted(xs)
    return {
        "n": float(len(xs)),
        "mean": _mean(xs),
        "std": _std(xs),
        "min": xs_sorted[0],
        "p50": xs_sorted[len(xs) // 2],
        "max": xs_sorted[-1],
    }


def run_audit(
    *,
    pepbench_csv: Path,
    catalog_tsv: Path,
    prep_manifest: Path,
    scoreable_ids_path: Path | None,
    oracle_artifact: Path | None,
) -> dict[str, Any]:
    pepbench = list(csv.DictReader(pepbench_csv.open(encoding="utf-8")))
    catalog = list(csv.DictReader(catalog_tsv.open(encoding="utf-8"), delimiter="\t"))
    prep = list(csv.DictReader(prep_manifest.open(encoding="utf-8"), delimiter="\t"))

    n_pepbench = len(pepbench)
    matched_ids = {r["record_id"] for r in catalog}
    prep_ok = {r["record_id"] for r in prep if r.get("ok") in {"1", "true", "True"}}
    prep_fail = {
        r["record_id"]: (r.get("error") or "unknown")
        for r in prep
        if r.get("ok") not in {"1", "true", "True"}
    }
    scoreable: set[str] = set()
    if scoreable_ids_path and scoreable_ids_path.is_file():
        scoreable = {x.strip() for x in scoreable_ids_path.read_text().split() if x.strip()}
    else:
        scoreable = {
            r["record_id"]
            for r in prep
            if r.get("ok") in {"1", "true", "True"} and r.get("scoreable_path")
        }

    drop_reasons = Counter()
    drop_reasons["pepbench_total"] = n_pepbench
    drop_reasons["matched_catalog"] = len(matched_ids)
    drop_reasons["no_rcsb_exact_chain_match"] = n_pepbench - len(matched_ids)
    drop_reasons["prep_failed"] = len(prep_fail)
    drop_reasons["prep_ok"] = len(prep_ok)
    drop_reasons["scoreable"] = len(scoreable)
    drop_reasons["prep_ok_not_scoreable"] = len(prep_ok - scoreable)

    fail_categories: Counter[str] = Counter()
    for err in prep_fail.values():
        key = err.split("—")[0].split(":")[0][:60]
        fail_categories[key] += 1

    # Survivor vs PepBench distributions
    def feats_from_pepbench(rows: list[dict[str, str]]) -> dict[str, list[float]]:
        charges, lengths, pkds = [], [], []
        for r in rows:
            pep = r["pep_seq"].strip().upper()
            if not (5 <= len(pep) <= 50):
                continue
            charges.append(charge_from_sequence(pep))
            lengths.append(float(len(pep)))
            pkds.append(float(r["label"]))
        return {"net_charge": charges, "peptide_len": lengths, "pKd": pkds}

    def feats_from_catalog(rows: list[dict[str, str]], keep: set[str]) -> dict[str, list[float]]:
        charges, lengths, pkds, res = [], [], [], []
        for r in rows:
            if r["record_id"] not in keep:
                continue
            pep = r["peptide_seq"].strip().upper()
            charges.append(charge_from_sequence(pep))
            lengths.append(float(len(pep)))
            pkds.append(float(r["pKd"]))
            if r.get("resolution"):
                try:
                    res.append(float(r["resolution"]))
                except ValueError:
                    pass
        return {
            "net_charge": charges,
            "peptide_len": lengths,
            "pKd": pkds,
            "resolution": res,
        }

    all_pb = feats_from_pepbench(pepbench)
    survivors = feats_from_catalog(catalog, scoreable)

    # Homogeneous-subset check: |net_charge| ≤ 1
    neutralish = [
        r
        for r in catalog
        if r["record_id"] in scoreable
        and abs(charge_from_sequence(r["peptide_seq"])) <= 1.0
    ]

    oracle_note = None
    if oracle_artifact and oracle_artifact.is_file():
        art = json.loads(oracle_artifact.read_text(encoding="utf-8"))
        oracle_note = {
            "spearman": art.get("spearman"),
            "baseline_rhos": art.get("baseline_rhos"),
            "n": art.get("n"),
        }

    # Call: flawed task vs bad oracle
    surv_charge_mean = _mean(survivors["net_charge"])
    pb_charge_mean = _mean(all_pb["net_charge"])
    charge_skew = abs(surv_charge_mean - pb_charge_mean) > 0.5
    call = (
        "selection_skew_contributes"
        if charge_skew
        else "oracle_genuinely_weak_vs_charge_baseline"
    )
    if charge_skew:
        narrative = (
            f"Survivors are charge-shifted vs PepBench "
            f"(mean net_charge {surv_charge_mean:.2f} vs {pb_charge_mean:.2f}). "
            "Part of the net_charge baseline win may be a selection artifact; "
            "still require the oracle to beat charge within homogeneous subsets."
        )
    else:
        narrative = (
            f"Survivor mean net_charge ({surv_charge_mean:.2f}) is close to PepBench "
            f"({pb_charge_mean:.2f}); net_charge winning is unlikely to be pure "
            "selection bias — treat as oracle electrostatic/solvation failure."
        )

    return {
        "funnel": {
            "pepbench_n": n_pepbench,
            "matched_n": len(matched_ids),
            "prep_ok_n": len(prep_ok),
            "scoreable_n": len(scoreable),
            "retention_vs_pepbench": len(scoreable) / n_pepbench if n_pepbench else 0.0,
        },
        "drop_reasons": dict(drop_reasons),
        "prep_fail_categories": dict(fail_categories.most_common()),
        "distributions": {
            "pepbench_all": {k: _summarize(v) for k, v in all_pb.items()},
            "scoreable_survivors": {k: _summarize(v) for k, v in survivors.items()},
        },
        "homogeneous_neutralish_n": len(neutralish),
        "charge_skew_detected": charge_skew,
        "call": call,
        "narrative": narrative,
        "oracle_artifact_snapshot": oracle_note,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pepbench",
        type=Path,
        default=Path("benchmarks/peptide_affinity/data/raw/PpI_ba_raw.csv"),
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path("benchmarks/peptide_affinity/data/peptide_affinity_catalog_v2.tsv"),
    )
    parser.add_argument(
        "--prep-manifest",
        type=Path,
        default=Path("benchmarks/peptide_affinity/data/prepared/prep_manifest.tsv"),
    )
    parser.add_argument(
        "--scoreable-ids",
        type=Path,
        default=Path("benchmarks/peptide_affinity/data/prepared/scoreable_success_ids.txt"),
    )
    parser.add_argument(
        "--oracle-artifact",
        type=Path,
        default=Path("benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/peptide_affinity/data/attrition_audit.json"),
    )
    args = parser.parse_args()
    if not args.pepbench.is_file():
        raise SystemExit(
            f"PepBench CSV missing: {args.pepbench}. See benchmarks/peptide_affinity/README.md"
        )
    report = run_audit(
        pepbench_csv=args.pepbench,
        catalog_tsv=args.catalog,
        prep_manifest=args.prep_manifest,
        scoreable_ids_path=args.scoreable_ids,
        oracle_artifact=args.oracle_artifact,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("funnel", "call", "narrative")}, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

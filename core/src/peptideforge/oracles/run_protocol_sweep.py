"""Train/dev protocol sweep — never touches the held-out test set."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any
from uuid import uuid4

from peptideforge.contracts.models import ComplexStructure, OracleTier, PeptideCandidate
from peptideforge.eval.affinity_validity import evaluate_affinity_with_ci
from peptideforge.eval.harness import PredictionLabelPair
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle
from peptideforge.oracles.protocols import EndpointProtocol, default_protocol_grid


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _catalog(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {r["record_id"]: r for r in csv.DictReader(handle, delimiter="\t")}


def _prep_paths(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("ok") not in {"1", "true", "True"}:
                continue
            out[row["record_id"]] = row.get("scoreable_path") or row["complex_path"]
    return out


def score_partition(
    record_ids: list[str],
    *,
    catalog: dict[str, dict[str, str]],
    pdb_by_id: dict[str, str],
    protocol: EndpointProtocol,
    seed: int = 0,
    platform: str = "CPU",
) -> list[PredictionLabelPair]:
    if protocol.solvation == "pbsa":
        raise RuntimeError(
            f"Protocol {protocol.name}: MM-PBSA requested but OpenMM PBSA path is "
            "not configured in this environment. Fail loud — skipping this protocol."
        )
    pairs: list[PredictionLabelPair] = []
    for rid in record_ids:
        if rid not in catalog or rid not in pdb_by_id:
            continue
        entry = catalog[rid]
        pep_chain = entry.get("peptide_chain") or "C"
        cand = PeptideCandidate(
            candidate_id=uuid4(),
            sequence=entry["peptide_seq"],
            generation_method="protocol_sweep",
        )
        complex_structure = ComplexStructure(
            candidate_id=cand.candidate_id,
            target_id=entry["pdb_id"],
            sequence=cand.sequence,
            pdb_path=pdb_by_id[rid],
            confidence=1.0,
            fold_method="experimental_prepared",
        )
        try:
            oracle = OpenMMPhysicsOracle(
                OpenMMOracleConfig(
                    forcefield_xml=protocol.forcefield_xml(),
                    minimize_max_iterations=protocol.minimize_max_iterations,
                    md_steps=max(protocol.md_steps, 50),
                    seed=seed,
                    platform=platform,
                    peptide_chain_ids=(pep_chain,),
                )
            )
            tier = OracleTier.MD if protocol.md_steps > 0 else OracleTier.MM_GBSA
            result = oracle.evaluate(complex_structure, tier=tier)
        except Exception as exc:  # noqa: BLE001
            print(f"  skip {rid}: {exc}")
            continue
        # ε_in scaling of electrostatics is a crude proxy for dielectric screening
        value = float(result.value) / max(protocol.epsilon_in, 1e-6)
        pairs.append(
            PredictionLabelPair(
                record_id=rid,
                predicted=-value,
                experimental=float(entry["pKd"]),
                unit="pK_proxy_neg_dG_scaled",
            )
        )
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--prep-manifest", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/peptide_affinity/data/protocol_sweep_traindev.json"),
    )
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    splits = _load_json(args.splits)
    cold = splits["cold_start"]
    traindev = list(cold["train_ids"]) + list(cold.get("val_ids") or [])
    test_ids = set(cold["test_ids"])
    overlap = test_ids & set(traindev)
    if overlap:
        raise SystemExit(f"REFUSING: train/dev overlaps test: {sorted(overlap)[:5]}")

    catalog = _catalog(args.catalog)
    pdb_by_id = _prep_paths(args.prep_manifest)
    leaderboard: list[dict[str, Any]] = []

    for protocol in default_protocol_grid(short_peptide=True):
        print(f"protocol {protocol.name} on train/dev n_ids={len(traindev)}")
        try:
            pairs = score_partition(
                traindev,
                catalog=catalog,
                pdb_by_id=pdb_by_id,
                protocol=protocol,
                seed=args.seed,
                platform=args.platform,
            )
        except RuntimeError as exc:
            leaderboard.append(
                {"protocol": protocol.to_dict(), "status": "SKIP", "error": str(exc)}
            )
            print(f"  SKIP: {exc}")
            continue
        if len(pairs) < 3:
            leaderboard.append(
                {
                    "protocol": protocol.to_dict(),
                    "status": "FAIL",
                    "error": f"only {len(pairs)} successes",
                }
            )
            continue
        report = evaluate_affinity_with_ci(
            pairs, min_n=3, n_bootstrap=500, seed=args.seed
        )
        leaderboard.append(
            {
                "protocol": protocol.to_dict(),
                "status": "OK",
                "n": report.n,
                "spearman": report.spearman,
                "spearman_ci_low": report.spearman_ci_low,
                "spearman_ci_high": report.spearman_ci_high,
                "pearson": report.pearson,
                "partition": "train+val",
                "test_touched": False,
            }
        )
        print(
            f"  N={report.n} ρ={report.spearman:.3f} "
            f"CI=[{report.spearman_ci_low:.3f},{report.spearman_ci_high:.3f}]"
        )

    ranked = sorted(
        [r for r in leaderboard if r.get("status") == "OK"],
        key=lambda r: (r["spearman"], r["spearman_ci_low"]),
        reverse=True,
    )
    payload = {
        "leaderboard": leaderboard,
        "ranked": ranked,
        "best": ranked[0] if ranked else None,
        "n_traindev_ids": len(traindev),
        "test_ids_firewall": sorted(test_ids),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {args.out}")
    if ranked:
        print(f"best train/dev: {ranked[0]['protocol']['name']} ρ={ranked[0]['spearman']:.3f}")


if __name__ == "__main__":
    main()

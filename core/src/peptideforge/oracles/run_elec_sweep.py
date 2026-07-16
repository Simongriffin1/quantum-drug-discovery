"""Electrostatics-led train/dev sweep (ε_in, salt, GB model) — test firewalled."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any
from uuid import uuid4

from peptideforge.contracts.models import ComplexStructure, OracleTier, PeptideCandidate
from peptideforge.eval.affinity_validity import (
    charge_from_sequence,
    evaluate_affinity_with_ci,
    run_trivial_baseline_battery,
    trivial_baselines_for_structures,
)
from peptideforge.eval.harness import PredictionLabelPair
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle


def _gb_xml(model: str) -> str:
    return {
        "obc1": "implicit/obc1.xml",
        "obc2": "implicit/obc2.xml",
        "gbn2": "implicit/gbn2.xml",
    }[model]


def elec_grid_dicts() -> list[dict[str, Any]]:
    """Small physically motivated grid — not a product over all axes."""
    out: list[dict[str, Any]] = []
    for eps in (1.0, 2.0, 4.0, 8.0):
        for salt in (0.0, 0.15):
            for gb in ("obc2", "gbn2"):
                out.append(
                    {
                        "name": f"gbsa_{gb}_eps{eps:g}_salt{salt:g}_min0",
                        "gb_model": gb,
                        "solute_dielectric": eps,
                        "salt_conc_M": salt,
                        "minimize_max_iterations": 0,
                    }
                )
    return out


def score_ids(
    record_ids: list[str],
    *,
    catalog: dict[str, dict[str, str]],
    pdb_by_id: dict[str, str],
    cfg: dict[str, Any],
    platform: str = "CPU",
    seed: int = 0,
) -> list[PredictionLabelPair]:
    pairs: list[PredictionLabelPair] = []
    gb = _gb_xml(str(cfg["gb_model"]))
    for rid in record_ids:
        if rid not in catalog or rid not in pdb_by_id:
            continue
        entry = catalog[rid]
        pep = entry.get("peptide_chain") or "C"
        cand = PeptideCandidate(
            candidate_id=uuid4(),
            sequence=entry["peptide_seq"],
            generation_method="elec_sweep",
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
                    forcefield_xml=("amber14-all.xml", gb),
                    minimize_max_iterations=int(cfg["minimize_max_iterations"]),
                    platform=platform,
                    seed=seed,
                    peptide_chain_ids=(pep,),
                    solute_dielectric=float(cfg["solute_dielectric"]),
                    salt_conc_M=float(cfg["salt_conc_M"]),
                )
            )
            result = oracle.evaluate(complex_structure, tier=OracleTier.MM_GBSA)
        except Exception as exc:  # noqa: BLE001
            print(f"  skip {rid}: {exc}", flush=True)
            continue
        pairs.append(
            PredictionLabelPair(
                record_id=rid,
                predicted=-float(result.value),
                experimental=float(entry["pKd"]),
                unit="pK_proxy_neg_dG",
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
        default=Path("benchmarks/peptide_affinity/data/elec_sweep_traindev.json"),
    )
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    splits = json.loads(args.splits.read_text(encoding="utf-8"))
    cold = splits["cold_start"]
    traindev = list(cold["train_ids"]) + list(cold.get("val_ids") or [])
    test_ids = set(cold["test_ids"])
    if set(traindev) & test_ids:
        raise SystemExit("REFUSING: train/dev overlaps test")

    catalog = {
        r["record_id"]: r
        for r in csv.DictReader(args.catalog.open(encoding="utf-8"), delimiter="\t")
    }
    pdb_by_id: dict[str, str] = {}
    for row in csv.DictReader(args.prep_manifest.open(encoding="utf-8"), delimiter="\t"):
        if row.get("ok") in {"1", "true", "True"}:
            pdb_by_id[row["record_id"]] = row.get("scoreable_path") or row["complex_path"]

    leaderboard: list[dict[str, Any]] = []
    for cfg in elec_grid_dicts():
        print(f"RUN {cfg['name']} n_ids={len(traindev)}", flush=True)
        pairs = score_ids(
            traindev,
            catalog=catalog,
            pdb_by_id=pdb_by_id,
            cfg=cfg,
            platform=args.platform,
            seed=args.seed,
        )
        if len(pairs) < 3:
            leaderboard.append({"config": cfg, "status": "FAIL", "n": len(pairs)})
            continue
        report = evaluate_affinity_with_ci(pairs, min_n=3, n_bootstrap=400, seed=args.seed)
        lengths = {p.record_id: len(catalog[p.record_id]["peptide_seq"]) for p in pairs}
        charges = {
            p.record_id: charge_from_sequence(catalog[p.record_id]["peptide_seq"]) for p in pairs
        }
        baselines = trivial_baselines_for_structures(
            pairs, peptide_lengths=lengths, net_charges=charges
        )
        beat_ok, base_rhos, base_msg = run_trivial_baseline_battery(pairs, baselines)
        leaderboard.append(
            {
                "config": cfg,
                "status": "OK",
                "n": report.n,
                "spearman": report.spearman,
                "spearman_ci_low": report.spearman_ci_low,
                "spearman_ci_high": report.spearman_ci_high,
                "pearson": report.pearson,
                "beats_trivial": beat_ok,
                "baseline_rhos": base_rhos,
                "baseline_msg": base_msg,
                "partition": "train+val",
                "test_touched": False,
            }
        )
        print(
            f"  ρ={report.spearman:.3f} CI=[{report.spearman_ci_low:.3f},"
            f"{report.spearman_ci_high:.3f}] beats_trivial={beat_ok}",
            flush=True,
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
        "n_traindev": len(traindev),
        "test_ids_firewall": sorted(test_ids),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {args.out}")
    if ranked:
        b = ranked[0]
        print(f"BEST {b['config']['name']} ρ={b['spearman']:.3f}")


if __name__ == "__main__":
    main()

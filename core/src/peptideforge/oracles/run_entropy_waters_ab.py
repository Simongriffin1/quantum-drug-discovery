"""Train/dev A/B for entropy + interface-water flags (test firewalled)."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from uuid import uuid4

from peptideforge.contracts.models import ComplexStructure, OracleTier, PeptideCandidate
from peptideforge.eval.affinity_validity import evaluate_affinity_with_ci
from peptideforge.eval.harness import PredictionLabelPair
from peptideforge.oracles.entropy import (
    apply_entropy_correction,
    truncated_nma_entropy_proxy,
)
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle
from peptideforge.oracles.waters import merge_waters_into_receptor_pdb


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--prep-manifest", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    catalog = {
        r["record_id"]: r
        for r in csv.DictReader(args.catalog.open(newline="", encoding="utf-8"), delimiter="\t")
    }
    prep = {
        r["record_id"]: r
        for r in csv.DictReader(
            args.prep_manifest.open(newline="", encoding="utf-8"), delimiter="\t"
        )
        if r.get("ok") in {"1", "true", "True"}
    }
    splits = json.loads(args.splits.read_text(encoding="utf-8"))
    traindev = list(splits["cold_start"]["train_ids"]) + list(
        splits["cold_start"]["val_ids"]
    )

    configs = [
        {"name": "baseline", "entropy": False, "waters": False},
        {"name": "entropy_on", "entropy": True, "waters": False},
        {"name": "waters_on", "entropy": False, "waters": True},
        {"name": "entropy_and_waters", "entropy": True, "waters": True},
    ]
    leaderboard = []
    work = args.out.parent / "ab_water_pdbs"
    work.mkdir(parents=True, exist_ok=True)

    for cfg in configs:
        pairs: list[PredictionLabelPair] = []
        for rid in traindev:
            if rid not in catalog or rid not in prep:
                continue
            entry = catalog[rid]
            row = prep[rid]
            pdb_path = row.get("scoreable_path") or row.get("complex_path")
            if not pdb_path or not Path(pdb_path).is_file():
                continue
            if cfg["waters"]:
                # Prep complex already retains flagged interface waters; scoreable strips them.
                base = row.get("complex_path") or ""
                water_src = row.get("interface_water_path") or None
                if base and Path(base).is_file():
                    pdb_path = base
                if water_src and Path(water_src).is_file() and Path(water_src).stat().st_size > 200:
                    merged = work / f"{rid}_with_waters.pdb"
                    try:
                        merge_waters_into_receptor_pdb(
                            Path(pdb_path),
                            Path(water_src),
                            merged,
                        )
                        pdb_path = str(merged)
                    except Exception as exc:  # noqa: BLE001
                        print(f"  water merge fail {rid}: {exc}", flush=True)
                        # Fall back to complex (may already include interface waters)
                        if not base or not Path(base).is_file():
                            continue
                        pdb_path = base
                elif not base or not Path(base).is_file():
                    print(f"  no water-bearing complex for {rid}", flush=True)
                    continue
            pep = entry.get("peptide_chain") or "C"
            cand = PeptideCandidate(
                candidate_id=uuid4(),
                sequence=entry["peptide_seq"],
                generation_method="ab",
            )
            cs = ComplexStructure(
                candidate_id=cand.candidate_id,
                target_id=entry["pdb_id"],
                sequence=cand.sequence,
                pdb_path=pdb_path,
                confidence=1.0,
                fold_method="experimental_prepared",
            )
            try:
                oracle = OpenMMPhysicsOracle(
                    OpenMMOracleConfig(
                        forcefield_xml=("amber14-all.xml", "implicit/gbn2.xml"),
                        minimize_max_iterations=0,
                        platform=args.platform,
                        seed=args.seed,
                        peptide_chain_ids=(pep,),
                        solute_dielectric=1.0,
                        salt_conc_M=0.0,
                    )
                )
                dg = float(oracle.evaluate(cs, tier=OracleTier.MM_GBSA).value)
                if cfg["entropy"]:
                    # Cheap proxy: ~peptide length as interface heavy-atom stand-in
                    n_heavy = max(10, 4 * len(entry["peptide_seq"]))
                    ent = truncated_nma_entropy_proxy(n_heavy_atoms_interface=n_heavy)
                    dg = apply_entropy_correction(dg, ent)
            except Exception as exc:  # noqa: BLE001
                print(f"  skip {rid}: {exc}", flush=True)
                continue
            pairs.append(
                PredictionLabelPair(
                    record_id=rid,
                    predicted=-dg,
                    experimental=float(entry["pKd"]),
                    unit="pK_proxy_neg_dG",
                )
            )
        if len(pairs) < 3:
            leaderboard.append({"config": cfg, "status": "TOO_FEW", "n": len(pairs)})
            continue
        report = evaluate_affinity_with_ci(
            pairs, min_n=3, n_bootstrap=400, seed=args.seed
        )
        leaderboard.append(
            {
                "config": cfg,
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
            f"{cfg['name']}: n={report.n} ρ={report.spearman:.3f} "
            f"CI=[{report.spearman_ci_low:.3f},{report.spearman_ci_high:.3f}]",
            flush=True,
        )

    leaderboard.sort(
        key=lambda r: (-(r.get("spearman") or -999), -(r.get("n") or 0))
    )
    payload = {
        "leaderboard": leaderboard,
        "chosen_hint": leaderboard[0] if leaderboard else None,
        "notes": (
            "Entropy uses truncated-NMA proxy (not MD Interaction Entropy). "
            "Waters use prep-flagged interface waters merged into complex."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"best": payload["chosen_hint"]}, indent=2))


if __name__ == "__main__":
    main()

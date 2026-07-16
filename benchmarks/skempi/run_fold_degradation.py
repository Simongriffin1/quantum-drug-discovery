"""Step 4A.3 — fold→score degradation: re-score predicted structures vs experimental ref.

Experimental ρ=0.381 is the REFERENCE and is never overwritten.
Reports MODE A / MODE B ρ with CIs and bootstrap paired Δρ.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any

from peptideforge.eval.affinity_validity import evaluate_affinity_with_ci
from peptideforge.eval.harness import PredictionLabelPair
from peptideforge.eval.metrics import spearman_rho
from peptideforge.eval.redteam import run_red_team
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle
from peptideforge.contracts.models import ComplexStructure, OracleTier, PeptideCandidate
from uuid import uuid4


EXPERIMENTAL_REFERENCE = {
    "spearman": 0.3806909215904269,
    "spearman_ci_low": 0.18889322486716895,
    "spearman_ci_high": 0.5556712698520563,
    "n": 100,
    "structures": "experimental_crystal_WT",
    "artifact": "benchmarks/skempi/data/skempi_ddg_powered_last_run.json",
}


def _score_ddg(
    wt_pdb: Path,
    mut_pdb: Path,
    *,
    peptide_chain: str,
    cfg: OpenMMOracleConfig,
) -> float:
    from skempi.run_skempi_ddg import score_ddg_pair

    pred, _ = score_ddg_pair(
        wt_pdb=wt_pdb, mut_pdb=mut_pdb, peptide_chain=peptide_chain, cfg=cfg
    )
    return pred


def bootstrap_paired_delta(
    exp_pairs: list[PredictionLabelPair],
    pred_pairs: list[PredictionLabelPair],
    *,
    n_bootstrap: int = 1000,
    seed: int = 0,
) -> tuple[float, float, float]:
    """Bootstrap CI for (ρ_experimental − ρ_predicted) on matched record_ids."""
    by_exp = {p.record_id: p for p in exp_pairs}
    by_pred = {p.record_id: p for p in pred_pairs}
    ids = sorted(set(by_exp) & set(by_pred))
    if len(ids) < 5:
        raise ValueError("need ≥5 matched ids for paired Δρ")
    rho_e = spearman_rho(
        [by_exp[i].predicted for i in ids], [by_exp[i].experimental for i in ids]
    )
    rho_p = spearman_rho(
        [by_pred[i].predicted for i in ids], [by_pred[i].experimental for i in ids]
    )
    point = rho_e - rho_p
    rng = random.Random(seed)
    samples: list[float] = []
    for _ in range(n_bootstrap):
        boot = [ids[rng.randrange(len(ids))] for _ in ids]
        re = spearman_rho(
            [by_exp[i].predicted for i in boot], [by_exp[i].experimental for i in boot]
        )
        rp = spearman_rho(
            [by_pred[i].predicted for i in boot], [by_pred[i].experimental for i in boot]
        )
        samples.append(re - rp)
    samples.sort()
    lo = samples[int(0.025 * len(samples))]
    hi = samples[int(0.975 * len(samples))]
    return point, lo, hi


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--holdout",
        type=Path,
        default=Path("benchmarks/skempi/data/skempi_powered_holdout_v1.json"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("benchmarks/skempi/data/predicted_folds/predicted_folds_manifest.json"),
    )
    parser.add_argument(
        "--experimental-artifact",
        type=Path,
        default=Path("benchmarks/skempi/data/skempi_ddg_powered_last_run.json"),
    )
    parser.add_argument(
        "--skempi-tsv",
        type=Path,
        default=Path("benchmarks/skempi/data/skempi_v2.tsv"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/skempi/data/fold_degradation_last_run.json"),
    )
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-pairs", type=int, default=100)
    args = parser.parse_args()

    holdout = json.loads(args.holdout.read_text(encoding="utf-8"))
    exp_art = json.loads(args.experimental_artifact.read_text(encoding="utf-8"))
    # Never overwrite experimental reference numbers
    experimental_ref = {
        **EXPERIMENTAL_REFERENCE,
        "spearman": exp_art.get("spearman", EXPERIMENTAL_REFERENCE["spearman"]),
        "spearman_ci_low": exp_art.get(
            "spearman_ci_low", EXPERIMENTAL_REFERENCE["spearman_ci_low"]
        ),
        "spearman_ci_high": exp_art.get(
            "spearman_ci_high", EXPERIMENTAL_REFERENCE["spearman_ci_high"]
        ),
        "n": exp_art.get("n", EXPERIMENTAL_REFERENCE["n"]),
    }

    if not args.manifest.is_file():
        payload = {
            "status": "BLOCKED",
            "block_reason": f"manifest missing: {args.manifest}",
            "experimental_reference": experimental_ref,
            "mode_a": None,
            "mode_b": None,
            "pre_registration": "ACCEPTANCE.md predicted-fold rule (committed before verdict)",
        }
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        raise SystemExit(payload["block_reason"])

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if manifest.get("status") == "BLOCKED_BOLTZ_UNAVAILABLE":
        payload = {
            "status": "BLOCKED_BOLTZ_UNAVAILABLE",
            "block_reason": manifest.get("error")
            or "Boltz-2 unavailable — predicted-fold ρ not measured",
            "experimental_reference": experimental_ref,
            "mode_a": {"gate_pass": False, "n": 0, "note": "not scored"},
            "mode_b": {"gate_pass": False, "n": 0, "note": "not scored"},
            "paired_delta_exp_minus_pred": None,
            "pre_registration": "ACCEPTANCE.md predicted-fold rule (committed before verdict)",
            "verdict": (
                "Predicted-fold campaigns remain BLOCKED. Experimental within-target "
                "authorization (ρ=0.381) is unchanged and is NOT re-cited as predicted ρ."
            ),
        }
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(json.dumps({"status": payload["status"], "verdict": payload["verdict"]}, indent=2))
        return

    # Score MODE A pairs when paths exist
    from peptideforge_benchmarks.skempi import load_skempi_ddg

    recs = {r.record_id: r for r in load_skempi_ddg(args.skempi_tsv)}
    cfg = OpenMMOracleConfig(
        forcefield_xml=("amber14-all.xml", "implicit/gbn2.xml"),
        minimize_max_iterations=0,
        platform=args.platform,
        seed=args.seed,
        solute_dielectric=1.0,
        salt_conc_M=0.0,
    )
    mode_a_pairs: list[PredictionLabelPair] = []
    details: list[dict[str, Any]] = []
    n_done = 0
    for prov in manifest.get("provenance") or []:
        if prov.get("mode") != "mutate_in_place":
            continue
        if prov.get("error") or not prov.get("wt_fold_path") or not prov.get("mutant_fold_path"):
            continue
        rid = prov["record_id"]
        rec = recs.get(rid)
        if rec is None:
            continue
        pep = (rec.partner2 or "B")[0]
        cfg_i = OpenMMOracleConfig(
            forcefield_xml=cfg.forcefield_xml,
            minimize_max_iterations=0,
            platform=args.platform,
            seed=args.seed,
            solute_dielectric=1.0,
            salt_conc_M=0.0,
            peptide_chain_ids=(pep, "P", "B"),
        )
        try:
            pred = _score_ddg(
                Path(prov["wt_fold_path"]),
                Path(prov["mutant_fold_path"]),
                peptide_chain=pep,
                cfg=cfg_i,
            )
        except Exception as exc:  # noqa: BLE001
            details.append({"record_id": rid, "error": str(exc)[:200]})
            continue
        # experimental ΔΔG from hold-out artifact details
        exp_ddg = None
        for d in exp_art.get("details") or []:
            if d.get("record_id") == rid and "experimental_ddg" in d:
                exp_ddg = float(d["experimental_ddg"])
                break
        if exp_ddg is None:
            exp_ddg = float(rec.ddg_kcal_mol)
        mode_a_pairs.append(
            PredictionLabelPair(
                record_id=rid, predicted=pred, experimental=exp_ddg, unit="kcal/mol"
            )
        )
        details.append(
            {"record_id": rid, "predicted_ddg": pred, "experimental_ddg": exp_ddg}
        )
        n_done += 1
        if n_done >= args.max_pairs:
            break

    if len(mode_a_pairs) < 3:
        payload = {
            "status": "BLOCKED",
            "block_reason": f"scored only {len(mode_a_pairs)} MODE A pairs",
            "experimental_reference": experimental_ref,
            "mode_a": {"n": len(mode_a_pairs), "gate_pass": False, "details": details},
            "mode_b": {"n": 0, "gate_pass": False, "note": "not primary"},
            "pre_registration": "ACCEPTANCE.md predicted-fold rule",
        }
        args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        raise SystemExit(payload["block_reason"])

    report = evaluate_affinity_with_ci(
        mode_a_pairs, min_n=30, n_bootstrap=1000, seed=args.seed
    )
    red = run_red_team(mode_a_pairs, seed=args.seed)
    gate = bool(
        report.spearman >= 0.30
        and report.spearman_ci_low > 0.0
        and red.passed
        and len(mode_a_pairs) >= 30
    )

    # Experimental pairs for paired Δ (use experimental *predictions* from crystal run)
    exp_pairs = [
        PredictionLabelPair(
            record_id=d["record_id"],
            predicted=float(d["predicted_ddg"]),
            experimental=float(d["experimental_ddg"]),
            unit="kcal/mol",
        )
        for d in exp_art.get("details") or []
        if "predicted_ddg" in d
    ]
    try:
        delta, dlo, dhi = bootstrap_paired_delta(exp_pairs, mode_a_pairs, seed=args.seed)
        paired = {
            "delta_exp_minus_pred": delta,
            "ci_low": dlo,
            "ci_high": dhi,
            "n_matched": len(set(p.record_id for p in mode_a_pairs) & set(p.record_id for p in exp_pairs)),
        }
    except ValueError as exc:
        paired = {"error": str(exc)}

    mode_a = {
        "spearman": report.spearman,
        "spearman_ci_low": report.spearman_ci_low,
        "spearman_ci_high": report.spearman_ci_high,
        "n": report.n,
        "gate_pass": gate,
        "red_team": {"passed": red.passed},
        "details": details,
        "mode": "mutate_in_place",
    }
    payload = {
        "status": "OK",
        "experimental_reference": experimental_ref,
        "mode_a": mode_a,
        "mode_b": {
            "note": "denovo comparison not scored in this run (MODE A primary)",
            "gate_pass": False,
            "n": 0,
        },
        "paired_delta_exp_minus_pred": paired,
        "table": {
            "experimental_rho": experimental_ref["spearman"],
            "predicted_mode_a_rho": report.spearman,
            "predicted_mode_a_ci": [report.spearman_ci_low, report.spearman_ci_high],
            "predicted_mode_b_rho": None,
            "paired_delta": paired,
        },
        "pre_registration": "ACCEPTANCE.md predicted-fold rule (committed before verdict)",
        "split_id": holdout.get("split_id"),
        "artifact": str(args.out),
        "verdict": (
            "PASS predicted-fold unconditional"
            if gate
            else "BLOCK predicted-fold unconditional — check stratification"
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "experimental_rho": experimental_ref["spearman"],
                "mode_a": {
                    k: mode_a[k]
                    for k in (
                        "n",
                        "spearman",
                        "spearman_ci_low",
                        "spearman_ci_high",
                        "gate_pass",
                    )
                },
                "paired": paired,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

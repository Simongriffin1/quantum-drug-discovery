"""SKEMPI v2 within-target ΔΔG validation harness.

Loads a local SKEMPI TSV (downloaded offline per README — never fabricates).
Scores WT vs mutant complexes with a fixed MM-GBSA protocol and correlates
predicted ΔΔG with experimental ΔΔG on a homology-aware held-out split.
"""

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
from peptideforge.eval.redteam import run_red_team
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle
from peptideforge_benchmarks.skempi import load_skempi_ddg
from peptideforge_benchmarks.splits import homology_aware_split


def _mutate_pdb_ca_proxy(pdb_path: Path, mutant_code: str, out_path: Path) -> None:
    """Best-effort single-residue rename for end-point scoring.

    SKEMPI mutant codes look like ``AI38A`` (chain A, Ile38→Ala). Without a full
    side-chain rebuild we only rename the residue in the PDB and rely on OpenMM
    hydrogens/templates — fail loud if residue not found.
    """
    # Parse: chain (1) + wt (1) + resnum + mut (1) — classic SKEMPI
    if len(mutant_code) < 4:
        raise ValueError(f"unparseable mutant code: {mutant_code}")
    chain = mutant_code[0]
    wt = mutant_code[1]
    mut = mutant_code[-1]
    num = mutant_code[2:-1]
    aa1to3 = {
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "D": "ASP",
        "C": "CYS",
        "Q": "GLN",
        "E": "GLU",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "V": "VAL",
    }
    mut3 = aa1to3.get(mut)
    if mut3 is None:
        raise ValueError(f"unknown mutant AA: {mut}")

    text = pdb_path.read_text(encoding="utf-8", errors="replace")
    found = False
    lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 26:
            cid = line[21]
            resseq = line[22:26].strip()
            if cid == chain and resseq == num:
                found = True
                line = line[:17] + f"{mut3:3s}" + line[20:]
        lines.append(line)
    if not found:
        raise ValueError(f"residue {chain}:{wt}{num} not found in {pdb_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def score_ddg_pair(
    *,
    wt_pdb: Path,
    mut_pdb: Path,
    peptide_chain: str,
    cfg: OpenMMOracleConfig,
) -> float:
    """Return predicted ΔΔG = ΔG_mut − ΔG_wt (kcal/mol)."""
    oracle = OpenMMPhysicsOracle(cfg)

    def score(path: Path) -> float:
        cand = PeptideCandidate(
            candidate_id=uuid4(), sequence="X", generation_method="skempi"
        )
        cs = ComplexStructure(
            candidate_id=cand.candidate_id,
            target_id=path.stem,
            sequence=cand.sequence,
            pdb_path=str(path),
            confidence=1.0,
            fold_method="experimental",
        )
        return float(oracle.evaluate(cs, tier=OracleTier.MM_GBSA).value)

    return score(mut_pdb) - score(wt_pdb)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skempi-tsv",
        type=Path,
        required=True,
        help="Local SKEMPI-like TSV (download offline; refuse missing).",
    )
    parser.add_argument(
        "--structure-dir",
        type=Path,
        required=True,
        help="Directory of WT PDBs named {pdb_id}.pdb",
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--test-fraction", type=float, default=0.3)
    parser.add_argument("--gb-model", default="gbn2", choices=("obc2", "gbn2", "obc1"))
    parser.add_argument("--solute-dielectric", type=float, default=1.0)
    parser.add_argument("--salt-conc", type=float, default=0.0)
    parser.add_argument("--max-pairs", type=int, default=80)
    args = parser.parse_args()

    if not args.skempi_tsv.is_file():
        raise SystemExit(
            f"Missing {args.skempi_tsv}. Download SKEMPI v2.0 CSV/TSV offline "
            "(https://life.bsc.es/pid/skempi2/) into benchmarks/skempi/data/ — "
            "refusing to fabricate mutation labels."
        )
    if not args.structure_dir.is_dir():
        raise SystemExit(
            f"Missing structure dir {args.structure_dir}. Place WT crystal PDBs "
            "named {{pdb_id}}.pdb — refusing to invent coordinates."
        )

    records = load_skempi_ddg(args.skempi_tsv)
    # Keep entries with usable WT PDB on disk
    usable = []
    for rec in records:
        wt = args.structure_dir / f"{rec.pdb_id}.pdb"
        if wt.is_file():
            usable.append(rec)
    if len(usable) < 6:
        raise SystemExit(
            f"only {len(usable)} SKEMPI rows have WT PDBs in {args.structure_dir} "
            f"(need ≥6). Download structures for the fixture PDB ids."
        )

    # Cluster by pdb_id if cluster_id missing
    clusters = {
        r.record_id: (r.cluster_id or r.pdb_id) for r in usable[: args.max_pairs]
    }
    split = homology_aware_split(
        tuple(clusters),
        clusters,
        test_fraction=args.test_fraction,
        seed=args.seed,
    )
    test_ids = set(split.test_ids)
    work_dir = args.out.parent / "skempi_mutant_pdbs"
    work_dir.mkdir(parents=True, exist_ok=True)

    gb = {
        "obc1": "implicit/obc1.xml",
        "obc2": "implicit/obc2.xml",
        "gbn2": "implicit/gbn2.xml",
    }[args.gb_model]
    cfg = OpenMMOracleConfig(
        forcefield_xml=("amber14-all.xml", gb),
        minimize_max_iterations=0,
        platform=args.platform,
        seed=args.seed,
        solute_dielectric=args.solute_dielectric,
        salt_conc_M=args.salt_conc,
        peptide_chain_ids=("B", "P", "L", "C", "D", "E"),
    )

    pairs: list[PredictionLabelPair] = []
    details: list[dict[str, Any]] = []
    for rec in usable[: args.max_pairs]:
        if rec.record_id not in test_ids:
            continue
        wt_pdb = args.structure_dir / f"{rec.pdb_id}.pdb"
        mut_pdb = work_dir / f"{rec.record_id}_mut.pdb"
        try:
            _mutate_pdb_ca_proxy(wt_pdb, rec.mutant, mut_pdb)
            # partner2 often peptide-like; use as peptide chain hint
            pep = (rec.partner2 or "B")[0]
            cfg_i = OpenMMOracleConfig(
                forcefield_xml=cfg.forcefield_xml,
                minimize_max_iterations=0,
                platform=args.platform,
                seed=args.seed,
                solute_dielectric=args.solute_dielectric,
                salt_conc_M=args.salt_conc,
                peptide_chain_ids=(pep,),
            )
            pred = score_ddg_pair(
                wt_pdb=wt_pdb, mut_pdb=mut_pdb, peptide_chain=pep, cfg=cfg_i
            )
        except Exception as exc:  # noqa: BLE001
            details.append({"record_id": rec.record_id, "error": str(exc)})
            continue
        pairs.append(
            PredictionLabelPair(
                record_id=rec.record_id,
                predicted=pred,
                experimental=float(rec.ddg_kcal_mol),
                unit="kcal/mol",
            )
        )
        details.append(
            {
                "record_id": rec.record_id,
                "predicted_ddg": pred,
                "experimental_ddg": float(rec.ddg_kcal_mol),
            }
        )

    if len(pairs) < 3:
        raise SystemExit(f"scored only {len(pairs)} held-out pairs — aborting")

    red = run_red_team(
        pairs,
        train_ids=list(split.train_ids),
        test_ids=list(split.test_ids),
        train_clusters={i: clusters[i] for i in split.train_ids if i in clusters},
        test_clusters={i: clusters[i] for i in split.test_ids if i in clusters},
        seed=args.seed,
    )
    report = evaluate_affinity_with_ci(
        pairs, min_n=3, n_bootstrap=1000, seed=args.seed, red_team=red
    )
    # Within-target gate (ACCEPTANCE.md): ρ≥0.30 and CI_low>0 — separate from
    # affinity min_n=40. Do not loosen after seeing results.
    gate_pass = bool(
        report.spearman >= 0.30
        and report.spearman_ci_low > 0.0
        and red.passed
    )

    payload = {
        "subset": "skempi_v2_within_target_heldout",
        "n": len(pairs),
        "spearman": report.spearman,
        "spearman_ci_low": report.spearman_ci_low,
        "spearman_ci_high": report.spearman_ci_high,
        "pearson": report.pearson,
        "gate_threshold": 0.30,
        "gate_pass": gate_pass,
        "red_team": {
            "passed": red.passed,
            "label_shuffle_passed": red.label_shuffle_passed,
            "trivial_baseline_passed": red.trivial_baseline_passed,
            "leakage_passed": red.leakage_passed,
            "model_rho": red.model_rho,
            "baseline_rho": red.baseline_rho,
        },
        "protocol": {
            "gb_model": args.gb_model,
            "solute_dielectric": args.solute_dielectric,
            "salt_conc_M": args.salt_conc,
        },
        "split": {
            "train": len(split.train_ids),
            "test": len(split.test_ids),
        },
        "details": details,
        "test_touched_affinity_set": False,
        "notes": (
            "CA-proxy mutation (residue rename only) — documents structural "
            "approximation; full side-chain rebuild is a follow-up."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({k: payload[k] for k in ("n", "spearman", "spearman_ci_low", "spearman_ci_high", "gate_pass", "red_team")}, indent=2))


if __name__ == "__main__":
    main()

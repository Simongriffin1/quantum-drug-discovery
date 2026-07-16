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
    """Mutate by stripping side chain beyond CB and renaming; hydrogens via PDBFixer.

    SKEMPI codes: ``LI18G`` = Leu→Gly on chain I residue 18 (wt+chain+num+mut).
    """
    try:
        from openmm.app import PDBFile, Modeller
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise ImportError(
            "SKEMPI mutate requires openmm+pdbfixer. Fail loud — no silent skip."
        ) from exc

    if len(mutant_code) < 4:
        raise ValueError(f"unparseable mutant code: {mutant_code}")
    wt = mutant_code[0]
    chain = mutant_code[1]
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

    keep_backbone = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3"}
    # ALA keeps CB; Gly keeps none beyond backbone
    keep_ala = keep_backbone | {"CB", "HB1", "HB2", "HB3"}

    pdb = PDBFile(str(pdb_path))
    modeller = Modeller(pdb.topology, pdb.positions)
    target = None
    for residue in modeller.topology.residues():
        if residue.chain.id == chain and residue.id.strip() == num:
            target = residue
            break
    if target is None:
        raise ValueError(f"residue {chain}:{wt}{num} not found in {pdb_path}")

    keep = keep_backbone if mut3 == "GLY" else keep_ala
    to_delete = [a for a in target.atoms() if a.name.strip() not in keep]
    if to_delete:
        modeller.delete(to_delete)
    # Rename remaining residue
    for residue in modeller.topology.residues():
        if residue.chain.id == chain and residue.id.strip() == num:
            residue.name = mut3
            break

    tmp = out_path.with_suffix(".tmp.pdb")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tmp.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle)

    fixer = PDBFixer(filename=str(tmp))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(f"REMARK mutate {mutant_code} -> {mut3}\n")
        PDBFile.writeFile(fixer.topology, fixer.positions, handle)
    tmp.unlink(missing_ok=True)


def _prepare_skempi_pdb(pdb_path: Path, out_path: Path, *, ph: float = 7.4) -> Path:
    """PDBFixer missing atoms/hydrogens for SKEMPI WT/mutant PDBs."""
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise ImportError(
            "SKEMPI prep requires openmm+pdbfixer. Fail loud — no silent skip."
        ) from exc
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    # Do not rebuild large missing loops — keep termini gaps as-is
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(f"REMARK skempi_prep source={pdb_path.name} ph={ph}\n")
        PDBFile.writeFile(fixer.topology, fixer.positions, handle)
    return out_path


def score_ddg_pair(
    *,
    wt_pdb: Path,
    mut_pdb: Path,
    peptide_chain: str,
    cfg: OpenMMOracleConfig,
    wt_energy: float | None = None,
    wt_prep_cache: dict[str, Path] | None = None,
    wt_energy_cache: dict[str, float] | None = None,
) -> tuple[float, float]:
    """Return (predicted ΔΔG, wt_energy) = (ΔG_mut − ΔG_wt, ΔG_wt)."""
    oracle = OpenMMPhysicsOracle(cfg)
    prep_dir = mut_pdb.parent / "prepared"
    cache_key = f"{wt_pdb.stem}::{peptide_chain}"

    if wt_prep_cache is not None and cache_key in wt_prep_cache:
        wt_prep = wt_prep_cache[cache_key]
    else:
        wt_prep = _prepare_skempi_pdb(wt_pdb, prep_dir / f"{wt_pdb.stem}_prep.pdb")
        if wt_prep_cache is not None:
            wt_prep_cache[cache_key] = wt_prep

    mut_prep = _prepare_skempi_pdb(mut_pdb, prep_dir / f"{mut_pdb.stem}_prep.pdb")

    def score(path: Path) -> float:
        cand = PeptideCandidate(
            candidate_id=uuid4(), sequence="AAAAA", generation_method="skempi"
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

    if wt_energy is not None:
        e_wt = wt_energy
    elif wt_energy_cache is not None and cache_key in wt_energy_cache:
        e_wt = wt_energy_cache[cache_key]
    else:
        e_wt = score(wt_prep)
        if wt_energy_cache is not None:
            wt_energy_cache[cache_key] = e_wt
    e_mut = score(mut_prep)
    return e_mut - e_wt, e_wt


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
    parser.add_argument("--test-fraction", type=float, default=0.35)
    parser.add_argument("--gb-model", default="gbn2", choices=("obc2", "gbn2", "obc1"))
    parser.add_argument("--solute-dielectric", type=float, default=1.0)
    parser.add_argument("--salt-conc", type=float, default=0.0)
    parser.add_argument("--max-pairs", type=int, default=800)
    parser.add_argument(
        "--min-test-n",
        type=int,
        default=100,
        help="Require ≥ this many successfully scored held-out pairs (powered gate).",
    )
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
    usable = [r for r in records if (args.structure_dir / f"{r.pdb_id}.pdb").is_file()]
    usable = usable[: args.max_pairs]
    if len(usable) < args.min_test_n:
        raise SystemExit(
            f"only {len(usable)} SKEMPI rows have WT PDBs (need ≥{args.min_test_n} "
            "pool for a powered held-out test). Download more structures."
        )

    clusters = {r.record_id: (r.cluster_id or r.pdb_id) for r in usable}
    split = homology_aware_split(
        tuple(clusters),
        clusters,
        test_fraction=args.test_fraction,
        seed=args.seed,
        split_name="skempi_homology_pdb_holdout",
    )
    test_ids = set(split.test_ids)
    n_test_assigned = sum(1 for r in usable if r.record_id in test_ids)
    print(
        f"pool={len(usable)} pdbs={len({r.pdb_id for r in usable})} "
        f"test_assigned={n_test_assigned} train={len(split.train_ids)}",
        flush=True,
    )
    if n_test_assigned < args.min_test_n:
        raise SystemExit(
            f"test assignment has only {n_test_assigned} rows (need ≥{args.min_test_n}). "
            "Increase --test-fraction or expand the mutation catalog."
        )

    work_dir = args.out.parent / "skempi_mutant_pdbs"
    work_dir.mkdir(parents=True, exist_ok=True)

    gb = {
        "obc1": "implicit/obc1.xml",
        "obc2": "implicit/obc2.xml",
        "gbn2": "implicit/gbn2.xml",
    }[args.gb_model]

    pairs: list[PredictionLabelPair] = []
    details: list[dict[str, Any]] = []
    wt_prep_cache: dict[str, Path] = {}
    wt_energy_cache: dict[str, float] = {}
    n_fail = 0

    for rec in usable:
        if rec.record_id not in test_ids:
            continue
        wt_pdb = args.structure_dir / f"{rec.pdb_id}.pdb"
        mut_pdb = work_dir / f"{rec.record_id}_mut.pdb"
        pep = (rec.partner2 or "B")[0]
        try:
            _mutate_pdb_ca_proxy(wt_pdb, rec.mutant, mut_pdb)
            cfg_i = OpenMMOracleConfig(
                forcefield_xml=("amber14-all.xml", gb),
                minimize_max_iterations=0,
                platform=args.platform,
                seed=args.seed,
                solute_dielectric=args.solute_dielectric,
                salt_conc_M=args.salt_conc,
                peptide_chain_ids=(pep,),
            )
            pred, _ = score_ddg_pair(
                wt_pdb=wt_pdb,
                mut_pdb=mut_pdb,
                peptide_chain=pep,
                cfg=cfg_i,
                wt_prep_cache=wt_prep_cache,
                wt_energy_cache=wt_energy_cache,
            )
        except Exception as exc:  # noqa: BLE001
            n_fail += 1
            details.append({"record_id": rec.record_id, "error": str(exc)[:200]})
            if n_fail <= 10 or n_fail % 25 == 0:
                print(f"  fail {rec.record_id}: {exc}", flush=True)
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
        if len(pairs) % 10 == 0:
            print(f"  scored {len(pairs)} / fails {n_fail}", flush=True)
        if len(pairs) >= args.min_test_n:
            print(
                f"  reached min_test_n={args.min_test_n}; stopping further scoring "
                f"(held-out membership fixed before scoring)",
                flush=True,
            )
            break

    if len(pairs) < args.min_test_n:
        raise SystemExit(
            f"scored only {len(pairs)} held-out pairs (need ≥{args.min_test_n}); "
            f"failures={n_fail}. Expand structures or fix mutation parse — "
            "refusing under-powered gate claim."
        )

    red = run_red_team(
        pairs,
        train_ids=list(split.train_ids),
        test_ids=[p.record_id for p in pairs],
        train_clusters={i: clusters[i] for i in split.train_ids if i in clusters},
        test_clusters={p.record_id: clusters[p.record_id] for p in pairs},
        seed=args.seed,
    )
    report = evaluate_affinity_with_ci(
        pairs, min_n=args.min_test_n, n_bootstrap=1000, seed=args.seed, red_team=red
    )
    gate_pass = bool(
        report.spearman >= 0.30
        and report.spearman_ci_low > 0.0
        and red.passed
        and len(pairs) >= args.min_test_n
    )

    payload = {
        "subset": "skempi_v2_within_target_heldout_powered",
        "n": len(pairs),
        "spearman": report.spearman,
        "spearman_ci_low": report.spearman_ci_low,
        "spearman_ci_high": report.spearman_ci_high,
        "pearson": report.pearson,
        "gate_threshold": 0.30,
        "min_test_n": args.min_test_n,
        "gate_pass": gate_pass,
        "n_failures": n_fail,
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
            "structures": "experimental_crystal_WT",
        },
        "split": {
            "train": len(split.train_ids),
            "test_assigned": n_test_assigned,
            "test_scored": len(pairs),
            "test_fraction": args.test_fraction,
            "clustering": "pdb_id_holdout",
        },
        "details": details,
        "test_touched_affinity_set": False,
        "prior_underpowered_run": "INVALIDATED_N16_CI_WIDTH",
        "notes": (
            "Side-chain strip to CB/backbone + PDBFixer rebuild; WT energies cached "
            "per complex. Experimental crystal WT structures only."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                k: payload[k]
                for k in (
                    "n",
                    "spearman",
                    "spearman_ci_low",
                    "spearman_ci_high",
                    "gate_pass",
                    "n_failures",
                    "red_team",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

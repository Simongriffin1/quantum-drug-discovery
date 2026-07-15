#!/usr/bin/env python3
"""Rebuild OpenMM-ready interface PDBs from shipped public ``*_full.pdb`` files.

Uses PDBFixer + OpenMM (no live network). Trims receptor residues far from the
peptide for CPU-friendly MM-GBSA while keeping a physically connected interface.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from openmm import Vec3, app, unit
from openmm.app import PDBFile
from pdbfixer import PDBFixer

# record_id → (full_pdb, receptor_chains, peptide_chain)
COMPLEXES: dict[str, tuple[str, tuple[str, ...], str]] = {
    "PP001": ("2VLL_full.pdb", ("A",), "C"),
    "PP003": ("1BD2_full.pdb", ("A",), "C"),
    "PP005": ("1DUZ_full.pdb", ("A",), "C"),
    "PP009": ("1QSE_full.pdb", ("A",), "C"),
    "PP011": ("2CLR_full.pdb", ("A",), "C"),
}


def _dist(a: Vec3, b: Vec3) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _residue_com(residue: app.topology.Residue, positions: list) -> Vec3 | None:
    coords = [positions[atom.index] for atom in residue.atoms()]
    if not coords:
        return None
    n = len(coords)
    return Vec3(
        sum(c.x for c in coords) / n,
        sum(c.y for c in coords) / n,
        sum(c.z for c in coords) / n,
    )


def trim_topology(
    topology: app.Topology,
    positions: list,
    *,
    keep_chains: set[str],
    peptide_chain: str,
    cutoff_nm: float,
) -> app.Modeller:
    """Keep peptide chain fully + receptor residues with COM within cutoff of peptide."""
    pos_nm = [Vec3(float(p.x), float(p.y), float(p.z)) for p in positions]
    pep_atoms = [
        pos_nm[a.index]
        for a in topology.atoms()
        if a.residue.chain.id == peptide_chain
    ]
    if not pep_atoms:
        raise ValueError(f"no atoms on peptide chain {peptide_chain}")

    # Contiguous receptor window covering all interface residues (avoids mid-chain cuts).
    chain_residues: dict[str, list] = {}
    for residue in topology.residues():
        chain_residues.setdefault(residue.chain.id, []).append(residue)

    keep_residues: set[int] = {r.index for r in chain_residues.get(peptide_chain, [])}

    for cid in keep_chains:
        residues = chain_residues.get(cid, [])
        if not residues:
            continue
        near_indices = []
        for i, residue in enumerate(residues):
            com = _residue_com(residue, pos_nm)
            if com is None:
                continue
            if any(_dist(com, pa) <= cutoff_nm for pa in pep_atoms):
                near_indices.append(i)
        if not near_indices:
            continue
        lo = max(0, min(near_indices) - 2)
        hi = min(len(residues) - 1, max(near_indices) + 2)
        for i in range(lo, hi + 1):
            keep_residues.add(residues[i].index)

    modeller = app.Modeller(topology, positions)
    to_delete = [r for r in topology.residues() if r.index not in keep_residues]
    modeller.delete(to_delete)
    return modeller


def prepare_one(
    full_path: Path,
    out_path: Path,
    *,
    receptor_chains: tuple[str, ...],
    peptide_chain: str,
    cutoff_A: float,
) -> dict[str, int]:
    fixer = PDBFixer(filename=str(full_path))
    keep = set(receptor_chains) | {peptide_chain}
    fixer.removeChains(chainIds=[c.id for c in fixer.topology.chains() if c.id not in keep])

    fixer.findMissingResidues()
    # Do not build missing terminal loops — clear missing residues far from ends
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    modeller = trim_topology(
        fixer.topology,
        fixer.positions,
        keep_chains=set(receptor_chains),
        peptide_chain=peptide_chain,
        cutoff_nm=cutoff_A / 10.0,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)

    topology = modeller.topology
    n_res = sum(1 for _ in topology.residues())
    n_atom = sum(1 for _ in topology.atoms())
    chains: dict[str, int] = {}
    for res in topology.residues():
        chains[res.chain.id] = chains.get(res.chain.id, 0) + 1
    return {"n_residues": n_res, "n_atoms": n_atom, **{f"chain_{k}": v for k, v in chains.items()}}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structures-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "fixtures" / "structures",
    )
    parser.add_argument("--cutoff", type=float, default=12.0)
    args = parser.parse_args()

    affinity_path = args.structures_dir.parent / "pdbbind_peptide_affinity_v1.tsv"
    pk: dict[str, str] = {}
    epitope: dict[str, str] = {}
    pdb_id: dict[str, str] = {}
    with affinity_path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pk[row["record_id"]] = row["pk"]
            epitope[row["record_id"]] = row["peptide_sequence"]
            pdb_id[row["record_id"]] = row["pdb_id"]

    manifest_rows: list[dict[str, str]] = []
    for record_id, (full_name, receptor_chains, peptide_chain) in COMPLEXES.items():
        full_path = args.structures_dir / full_name
        if not full_path.is_file():
            raise FileNotFoundError(full_path)
        out_name = f"{record_id}_{full_name.replace('_full.pdb', '')}_interface.pdb"
        out_path = args.structures_dir / out_name
        stats = prepare_one(
            full_path,
            out_path,
            receptor_chains=receptor_chains,
            peptide_chain=peptide_chain,
            cutoff_A=args.cutoff,
        )
        print(record_id, out_name, stats)
        manifest_rows.append(
            {
                "record_id": record_id,
                "pdb_id": full_name.replace("_full.pdb", ""),
                "pdb_path": out_name,
                "peptide_chain": peptide_chain,
                "receptor_chains": ",".join(receptor_chains),
                "epitope": epitope.get(record_id, ""),
                "experimental_pk": pk.get(record_id, ""),
            }
        )

    manifest = args.structures_dir / "structure_manifest_v1.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "record_id",
                "pdb_id",
                "pdb_path",
                "peptide_chain",
                "receptor_chains",
                "epitope",
                "experimental_pk",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    print("manifest →", manifest)


if __name__ == "__main__":
    main()

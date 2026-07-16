"""Interface crystallographic / bridging waters for MM-GBSA scoring.

Largest documented lever for PPI/peptide MM-GBSA when waters are retained
at the interface. Selection is deterministic and logged.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class WaterSelection:
    n_waters: int
    residue_indices: tuple[int, ...]
    cutoff_A: float
    method: str
    metadata: dict[str, Any]


def select_interface_waters(
    topology: Any,
    positions: Any,
    *,
    peptide_chain_ids: set[str],
    receptor_chain_ids: set[str],
    cutoff_A: float = 3.5,
) -> WaterSelection:
    """Keep HOH whose O is within ``cutoff_A`` of both peptide and receptor heavy atoms.

    Deterministic: residues ordered by topology index.
    """

    def heavy_coords(chain_ids: set[str]) -> list[tuple[float, float, float]]:
        out: list[tuple[float, float, float]] = []
        for atom in topology.atoms():
            if atom.residue.chain.id not in chain_ids:
                continue
            if atom.element is not None and atom.element.symbol == "H":
                continue
            if atom.residue.name in {"HOH", "WAT", "H2O"}:
                continue
            p = positions[atom.index]
            out.append((float(p.x), float(p.y), float(p.z)))
        return out

    pep = heavy_coords(peptide_chain_ids)
    rec = heavy_coords(receptor_chain_ids)
    if not pep or not rec:
        return WaterSelection(
            n_waters=0,
            residue_indices=(),
            cutoff_A=cutoff_A,
            method="interface_bridging",
            metadata={"reason": "empty peptide or receptor heavy atoms"},
        )

    def near(pt: tuple[float, float, float], cloud: list[tuple[float, float, float]]) -> bool:
        # positions are nm in OpenMM; convert cutoff Å → nm
        cut = cutoff_A / 10.0
        return any(
            math.sqrt((pt[0] - c[0]) ** 2 + (pt[1] - c[1]) ** 2 + (pt[2] - c[2]) ** 2) <= cut
            for c in cloud
        )

    kept: list[int] = []
    for residue in topology.residues():
        if residue.name not in {"HOH", "WAT", "H2O"}:
            continue
        o_atoms = [a for a in residue.atoms() if a.name.strip() in {"O", "OW"}]
        if not o_atoms:
            continue
        p = positions[o_atoms[0].index]
        pt = (float(p.x), float(p.y), float(p.z))
        if near(pt, pep) and near(pt, rec):
            kept.append(residue.index)

    kept.sort()
    return WaterSelection(
        n_waters=len(kept),
        residue_indices=tuple(kept),
        cutoff_A=cutoff_A,
        method="interface_bridging",
        metadata={},
    )


def merge_waters_into_receptor_pdb(
    complex_pdb: Path,
    interface_water_pdb: Path | None,
    out_pdb: Path,
    *,
    receptor_chain_ids: set[str] | None = None,
) -> WaterSelection:
    """Append interface waters from a prep artifact into the complex for scoring.

    Waters are left on their crystallographic chain; scoring treats non-peptide
    chains (including water) as receptor when ``include_interface_waters`` is set.
    Fail loud if complex missing; empty water file → zero waters (logged).
    """
    try:
        from openmm.app import PDBFile, Modeller
    except ImportError as exc:
        raise ImportError(
            "merge_waters_into_receptor_pdb requires OpenMM. Fail loud — no silent skip."
        ) from exc

    if not complex_pdb.is_file():
        raise FileNotFoundError(f"complex PDB missing: {complex_pdb}")

    pdb = PDBFile(str(complex_pdb))
    modeller = Modeller(pdb.topology, pdb.positions)
    n_added = 0
    if interface_water_pdb is not None and interface_water_pdb.is_file():
        w = PDBFile(str(interface_water_pdb))
        if w.topology.getNumAtoms() > 0:
            modeller.add(w.topology, w.positions)
            n_added = sum(1 for r in w.topology.residues() if r.name in {"HOH", "WAT", "H2O"})

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    with out_pdb.open("w", encoding="utf-8") as handle:
        handle.write(
            f"REMARK   1 interface_waters_merged n={n_added} "
            f"from={interface_water_pdb}\n"
        )
        PDBFile.writeFile(modeller.topology, modeller.positions, handle)

    return WaterSelection(
        n_waters=n_added,
        residue_indices=(),
        cutoff_A=3.5,
        method="prep_interface_waters_merge",
        metadata={
            "complex_pdb": str(complex_pdb),
            "interface_water_pdb": str(interface_water_pdb) if interface_water_pdb else None,
            "out_pdb": str(out_pdb),
            "receptor_chain_ids": sorted(receptor_chain_ids) if receptor_chain_ids else None,
        },
    )

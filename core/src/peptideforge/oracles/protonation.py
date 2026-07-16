"""Protonation at pH 7.4 for electrostatic fidelity.

Wrong His/Asp/Glu/Lys states corrupt Coulomb + GB terms — the components
diagnosed as non-informative in DIAGNOSIS.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ProtonationResult:
    ok: bool
    out_pdb: str | None
    method: str
    error: str | None = None


def protonate_pdbfixer(pdb_path: Path, out_pdb: Path, *, ph: float = 7.4) -> ProtonationResult:
    """PDBFixer addMissingHydrogens(ph) — always available with OpenMM stack."""
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise ImportError(
            "protonate_pdbfixer requires openmm + pdbfixer. Fail loud — no silent skip."
        ) from exc
    try:
        fixer = PDBFixer(filename=str(pdb_path))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        with out_pdb.open("w", encoding="utf-8") as handle:
            handle.write(f"REMARK   1 protonation method=pdbfixer ph={ph}\n")
            PDBFile.writeFile(fixer.topology, fixer.positions, handle)
        return ProtonationResult(ok=True, out_pdb=str(out_pdb), method="pdbfixer")
    except Exception as exc:  # noqa: BLE001
        return ProtonationResult(ok=False, out_pdb=None, method="pdbfixer", error=str(exc))


def protonate_pdb2pqr(pdb_path: Path, out_pqr: Path, *, ph: float = 7.4) -> ProtonationResult:
    """PDB2PQR + PROPKA — fail loud if not installed when requested."""
    try:
        from pdb2pqr.main import main as pdb2pqr_main  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "pdb2pqr/PROPKA not installed. pip install pdb2pqr  OR use method=pdbfixer. "
            "Refusing to silently fall back."
        ) from exc
    try:
        out_pqr.parent.mkdir(parents=True, exist_ok=True)
        pdb2pqr_main(
            [
                "--ff=AMBER",
                f"--with-ph={ph}",
                "--drop-water",
                str(pdb_path),
                str(out_pqr),
            ]
        )
        return ProtonationResult(ok=True, out_pdb=str(out_pqr), method="pdb2pqr_propka")
    except Exception as exc:  # noqa: BLE001
        return ProtonationResult(ok=False, out_pdb=None, method="pdb2pqr_propka", error=str(exc))

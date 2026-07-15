"""OpenMM import + PDB helpers. Missing deps fail loud — never fabricate energies."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any


class OpenMMUnavailableError(ImportError):
    """Raised when OpenMM (or openmm.app) cannot be imported."""


def require_openmm() -> tuple[Any, Any, Any]:
    """Import ``openmm``, ``openmm.app``, ``openmm.unit`` or raise loudly."""
    try:
        import openmm
        import openmm.app as app
        import openmm.unit as unit
    except ImportError as exc:
        raise OpenMMUnavailableError(
            "OpenMM is required for the physics oracle but is not installed. "
            "Install with: poetry install -E openmm   (or pip install openmm). "
            "Refusing to fabricate energies."
        ) from exc
    return openmm, app, unit


def resolve_pdb_path(pdb_path: str | None, pdb_text: str | None) -> Path:
    """Return a filesystem path to a PDB; write ``pdb_text`` to a temp file if needed."""
    if pdb_path is not None:
        path = Path(pdb_path)
        if not path.is_file():
            raise FileNotFoundError(f"complex PDB not found: {path}")
        return path
    if pdb_text is None:
        raise ValueError("ComplexStructure requires pdb_path or pdb_text")
    tmp = tempfile.NamedTemporaryFile(  # noqa: SIM115 — caller owns lifecycle via path
        mode="w",
        suffix=".pdb",
        delete=False,
        encoding="utf-8",
    )
    try:
        tmp.write(pdb_text)
        tmp.flush()
    finally:
        tmp.close()
    return Path(tmp.name)


def load_pdb(pdb_path: Path) -> tuple[Any, Any]:
    """Load topology + positions from a PDB file via OpenMM."""
    _, app, _ = require_openmm()
    pdb = app.PDBFile(str(pdb_path))
    return pdb.topology, pdb.positions


def all_chain_ids(topology: Any) -> list[str]:
    """Ordered unique chain IDs in a topology."""
    seen: list[str] = []
    for chain in topology.chains():
        cid = chain.id
        if cid not in seen:
            seen.append(cid)
    return seen


def chain_atom_indices(topology: Any, chain_ids: set[str]) -> list[int]:
    """Atom indices belonging to the given chain IDs."""
    indices: list[int] = []
    for atom in topology.atoms():
        if atom.residue.chain.id in chain_ids:
            indices.append(atom.index)
    return indices

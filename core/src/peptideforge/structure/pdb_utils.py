"""PDB helpers for structure prediction (sequence extraction, confidence)."""

from __future__ import annotations

from pathlib import Path

# Three-letter to one-letter map (standard amino acids).
_AA3_TO_1: dict[str, str] = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}


def read_pdb_text(path: str | Path) -> str:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"PDB not found: {p}")
    return p.read_text(encoding="utf-8", errors="replace")


def extract_chain_sequence(pdb_text: str, chain_id: str) -> str:
    """Extract one-letter sequence from ATOM records in order of appearance."""
    residues: list[tuple[int, str]] = []
    seen: set[tuple[str, int]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) < 22:
            continue
        if line[21].strip() != chain_id:
            continue
        resname = line[17:20].strip()
        if resname not in _AA3_TO_1:
            continue
        try:
            resseq = int(line[22:26].strip())
        except ValueError:
            continue
        key = (chain_id, resseq)
        if key in seen:
            continue
        seen.add(key)
        residues.append((resseq, _AA3_TO_1[resname]))
    residues.sort(key=lambda x: x[0])
    if not residues:
        raise ValueError(f"no standard residues found for chain {chain_id!r}")
    return "".join(aa for _, aa in residues)


def mean_plddt_from_pdb(pdb_text: str, *, chain_id: str | None = None) -> float:
    """Mean CA B-factor scaled to [0, 1] (Boltz/AF-style pLDDT storage)."""
    bfactors: list[float] = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) < 66:
            continue
        atom = line[12:16].strip()
        if atom != "CA":
            continue
        if chain_id is not None and line[21].strip() != chain_id:
            continue
        try:
            bfactors.append(float(line[60:66]))
        except ValueError:
            continue
    if not bfactors:
        raise ValueError("cannot derive confidence: no CA B-factors in PDB")
    mean_b = sum(bfactors) / len(bfactors)
    # pLDDT stored 0–100 in B-factor column
    return max(0.0, min(1.0, mean_b / 100.0))


def list_protein_chains(pdb_text: str) -> tuple[str, ...]:
    chains: set[str] = set()
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and len(line) >= 22:
            chains.add(line[21].strip())
    return tuple(sorted(chains))

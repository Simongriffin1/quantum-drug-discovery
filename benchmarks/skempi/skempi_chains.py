"""SKEMPI partner chain resolution for structure prediction."""

from __future__ import annotations

from peptideforge_benchmarks.models import MutationRecord
from peptideforge.structure.boltz2 import boltz_pdb_template_chain_id
from peptideforge.structure.pdb_utils import extract_chain_sequence, list_protein_chains

__all__ = [
    "boltz_pdb_template_chain_id",
    "skempi_chain_sequences",
    "skempi_partner_chain_ids",
]


def skempi_partner_chain_ids(rec: MutationRecord, pdb_text: str) -> list[str]:
    """Return ordered PDB chain ids for the SKEMPI complex (mutated partner first)."""
    if len(rec.mutant) < 2:
        raise ValueError(f"unparseable mutant code: {rec.mutant}")
    mut_chain = rec.mutant[1]
    p1 = rec.partner1 or ""
    p2 = rec.partner2 or ""
    present = set(list_protein_chains(pdb_text))

    ordered: list[str] = []
    if mut_chain in p2:
        ordered.extend(c for c in p2 if c in present)
        ordered.extend(c for c in p1 if c in present and c not in ordered)
    elif mut_chain in p1:
        ordered.extend(c for c in p1 if c in present)
        ordered.extend(c for c in p2 if c in present and c not in ordered)
    else:
        ordered.append(mut_chain)
        for c in p1 + p2:
            if c in present and c not in ordered:
                ordered.append(c)

    if mut_chain not in ordered:
        raise ValueError(
            f"mutation chain {mut_chain} not in SKEMPI partners for {rec.record_id}"
        )
    if len(ordered) < 2:
        raise ValueError(f"need ≥2 protein chains for {rec.record_id}, got {ordered}")
    return ordered


def skempi_chain_sequences(
    rec: MutationRecord, pdb_text: str
) -> list[tuple[str, str]]:
    """(pdb_chain_id, one_letter_sequence) tuples for Boltz input."""
    chains = skempi_partner_chain_ids(rec, pdb_text)
    return [(c, extract_chain_sequence(pdb_text, c)) for c in chains]


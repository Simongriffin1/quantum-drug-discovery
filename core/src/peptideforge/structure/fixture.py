"""Fixture structure predictor for CI / benchmark plumbing.

Returns precomputed interface PDBs from ``benchmarks/fixtures/structures/`` when
the candidate sequence matches a manifest epitope. Coordinates are experimental
fixtures — not fabricated. Unknown sequence/target pairs raise loudly.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from peptideforge.contracts.models import ComplexStructure, PeptideCandidate, Provenance
from peptideforge.structure.cache import FoldCache, FoldCacheKey
from peptideforge.structure.pdb_utils import read_pdb_text


@dataclass(frozen=True)
class ManifestRow:
    record_id: str
    pdb_id: str
    pdb_path: str
    peptide_chain: str
    epitope: str


def default_structures_dir() -> Path:
    """``benchmarks/fixtures/structures`` relative to repo layout."""
    here = Path(__file__).resolve()
    # core/src/peptideforge/structure -> repo root
    repo = here.parents[4]
    return repo / "benchmarks" / "fixtures" / "structures"


def load_manifest(manifest_path: Path | None = None) -> dict[str, ManifestRow]:
    """Load ``structure_manifest_v1.tsv`` keyed by pdb_id (uppercase)."""
    path = manifest_path or (default_structures_dir() / "structure_manifest_v1.tsv")
    if not path.is_file():
        raise FileNotFoundError(
            f"structure manifest missing: {path}. "
            "Benchmark fixtures must ship with CI (no live download)."
        )
    rows: dict[str, ManifestRow] = {}
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            entry = ManifestRow(
                record_id=row["record_id"],
                pdb_id=row["pdb_id"].upper(),
                pdb_path=row["pdb_path"],
                peptide_chain=row["peptide_chain"],
                epitope=row["epitope"].upper(),
            )
            rows[entry.pdb_id] = entry
    return rows


class FixtureStructurePredictor:
    """Return cached interface PDB when sequence matches benchmark epitope."""

    fold_method = "fixture_interface_pdb"

    def __init__(
        self,
        *,
        structures_dir: Path | None = None,
        manifest_path: Path | None = None,
        cache: FoldCache | None = None,
    ) -> None:
        self.structures_dir = structures_dir or default_structures_dir()
        self.manifest = load_manifest(manifest_path)
        self.cache = cache

    def fold(
        self,
        candidate: PeptideCandidate,
        *,
        target_id: str,
        target_structure: str,
        seed: int | None = None,
    ) -> ComplexStructure:
        key = FoldCacheKey(
            sequence=candidate.sequence,
            target_id=target_id,
            target_structure=target_structure,
        )
        if self.cache is not None:
            hit = self.cache.get(key)
            if hit is not None:
                return hit

        tid = target_id.upper()
        if tid not in self.manifest:
            raise ValueError(
                f"target_id {target_id!r} not in structure manifest; "
                "FixtureStructurePredictor only serves benchmark fixtures"
            )
        row = self.manifest[tid]
        if candidate.sequence.upper() != row.epitope:
            raise ValueError(
                f"sequence {candidate.sequence!r} does not match manifest epitope "
                f"{row.epitope!r} for target {target_id!r}. "
                "Refusing to return an unrelated experimental structure."
            )
        pdb_path = self.structures_dir / row.pdb_path
        if not pdb_path.is_file():
            raise FileNotFoundError(f"interface PDB missing: {pdb_path}")
        pdb_text = read_pdb_text(pdb_path)

        structure = ComplexStructure(
            candidate_id=candidate.candidate_id,
            target_id=tid,
            sequence=candidate.sequence,
            pdb_path=str(pdb_path.resolve()),
            pdb_text=pdb_text,
            confidence=1.0,
            fold_method=self.fold_method,
            cache_key=key.digest(),
            provenance=Provenance(
                data_version="structure_manifest_v1",
                tool_versions={"fixture_predictor": "0.1.0", "peptide_chain": row.peptide_chain},
            ),
        )
        if self.cache is not None:
            self.cache.put(key, structure)
        return structure

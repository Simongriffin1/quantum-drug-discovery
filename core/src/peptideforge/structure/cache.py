"""Disk cache for folded peptide–target complexes keyed by sequence hash."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from peptideforge.contracts.models import ComplexStructure, PeptideCandidate


@dataclass(frozen=True)
class FoldCacheKey:
    sequence: str
    target_id: str
    target_structure: str

    def digest(self) -> str:
        target_path = Path(self.target_structure)
        target_bytes = target_path.read_bytes() if target_path.is_file() else self.target_structure.encode()
        payload = f"{self.sequence}|{self.target_id}|".encode() + target_bytes
        return hashlib.sha256(payload).hexdigest()


class FoldCache:
    """Persist ``ComplexStructure`` payloads by fold hash (P6 acceptance: cache hits)."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)

    def path_for(self, key: FoldCacheKey) -> Path:
        return self.root / f"{key.digest()}.json"

    def get(self, key: FoldCacheKey) -> ComplexStructure | None:
        path = self.path_for(key)
        if not path.is_file():
            return None
        data = json.loads(path.read_text(encoding="utf-8"))
        return ComplexStructure.model_validate(data)

    def put(self, key: FoldCacheKey, structure: ComplexStructure) -> None:
        path = self.path_for(key)
        path.write_text(structure.model_dump_json(indent=2), encoding="utf-8")

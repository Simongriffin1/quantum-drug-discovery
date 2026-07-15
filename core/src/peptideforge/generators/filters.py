"""Validity and diversity filters for peptide sequence generation."""

from __future__ import annotations

import random
from collections.abc import Iterable, Sequence

from peptideforge.contracts.models import (
    CANONICAL_AA,
    MAX_PEPTIDE_LENGTH,
    MIN_PEPTIDE_LENGTH,
)


def normalize_allowed_residues(
    allowed: Iterable[str] | None,
) -> frozenset[str]:
    if allowed is None:
        return CANONICAL_AA
    upper = frozenset(aa.upper() for aa in allowed)
    invalid = upper - CANONICAL_AA
    if invalid:
        raise ValueError(f"non-canonical allowed residues: {sorted(invalid)}")
    return upper


def is_valid_sequence(
    seq: str,
    *,
    min_length: int = MIN_PEPTIDE_LENGTH,
    max_length: int = MAX_PEPTIDE_LENGTH,
    allowed_residues: frozenset[str] = CANONICAL_AA,
    fixed_positions: dict[int, str] | None = None,
) -> bool:
    """Return True if sequence passes length, alphabet, and fixed-position checks."""
    if not (min_length <= len(seq) <= max_length):
        return False
    if any(aa not in allowed_residues for aa in seq):
        return False
    if fixed_positions:
        for idx, aa in fixed_positions.items():
            if idx < 0 or idx >= len(seq) or seq[idx] != aa.upper():
                return False
    return True


def sequence_identity(a: str, b: str) -> float:
    """Fraction of identical positions (pad to max length with gaps counted as mismatch)."""
    if len(a) != len(b):
        n = max(len(a), len(b))
        a = a.ljust(n, "-")
        b = b.ljust(n, "-")
    if not a:
        return 1.0
    matches = sum(x == y for x, y in zip(a, b, strict=True))
    return matches / len(a)


def is_diverse_enough(
    candidate: str,
    pool: Sequence[str],
    *,
    max_identity: float = 0.95,
) -> bool:
    """Reject near-duplicates so physics budget is not spent redundantly."""
    return all(sequence_identity(candidate, existing) < max_identity for existing in pool)


def pick_mutation_positions(
    length: int,
    n_mutations: int,
    rng: random.Random,
    *,
    fixed_positions: frozenset[int] | None = None,
) -> list[int]:
    """Choose distinct mutable indices respecting fixed positions."""
    fixed = fixed_positions or frozenset()
    mutable = [i for i in range(length) if i not in fixed]
    if not mutable:
        return []
    k = min(n_mutations, len(mutable))
    return rng.sample(mutable, k)

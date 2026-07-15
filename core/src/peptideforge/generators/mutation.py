"""Constrained sequence-space mutation around seed peptides.

Biological rationale: local mutations explore binding-competent neighborhoods
without the cost of de novo generation; diversity filters avoid redundant
physics spend on near-identical sequences.
"""

from __future__ import annotations

import random
from typing import Any

from peptideforge.contracts.models import (
    CANONICAL_AA,
    MAX_PEPTIDE_LENGTH,
    MIN_PEPTIDE_LENGTH,
    Candidates,
    PeptideCandidate,
)
from peptideforge.generators.filters import (
    is_diverse_enough,
    is_valid_sequence,
    normalize_allowed_residues,
    pick_mutation_positions,
)


class MutationGenerator:
    """Valid, diverse peptide candidates via constrained point mutation."""

    def __init__(
        self,
        *,
        min_length: int = 5,
        max_length: int = 50,
        allowed_residues: frozenset[str] | None = None,
        mutations_per_child: tuple[int, int] = (1, 3),
        max_identity: float = 0.95,
        generation_method: str = "mutation",
        max_attempts_per_child: int = 200,
    ) -> None:
        self.min_length = min_length
        self.max_length = max_length
        self.allowed_residues = normalize_allowed_residues(
            allowed_residues if allowed_residues is not None else None
        )
        lo, hi = mutations_per_child
        if lo < 1 or hi < lo:
            raise ValueError("mutations_per_child must be (lo, hi) with 1 <= lo <= hi")
        self.mutations_per_child = mutations_per_child
        self.max_identity = max_identity
        self.generation_method = generation_method
        self.max_attempts_per_child = max_attempts_per_child

    def propose(
        self,
        *,
        n: int,
        seed_sequences: tuple[str, ...] | None = None,
        seed: int | None = None,
        constraints: dict[str, object] | None = None,
    ) -> Candidates:
        if n < 1:
            raise ValueError("n must be >= 1")
        rng = random.Random(seed)
        merged = dict(constraints or {})
        merged.setdefault("min_length", self.min_length)
        merged.setdefault("max_length", self.max_length)
        cfg = _parse_constraints(merged, default_allowed=self.allowed_residues)
        seeds = _normalize_seeds(seed_sequences, cfg, rng)
        selected: list[str] = []
        items: list[PeptideCandidate] = []

        while len(items) < n:
            parent = rng.choice(seeds)
            child = self._mutate_once(parent, rng, cfg)
            if child is None:
                continue
            if not is_diverse_enough(child, selected, max_identity=self.max_identity):
                continue
            selected.append(child)
            items.append(
                PeptideCandidate(
                    sequence=child,
                    generation_method=self.generation_method,
                    metadata={"parent_sequence": parent},
                )
            )

        return Candidates(items=tuple(items), seed=seed)

    def _mutate_once(
        self,
        parent: str,
        rng: random.Random,
        cfg: _MutationConfig,
    ) -> str | None:
        seq = list(parent)
        n_mut = rng.randint(*self.mutations_per_child)
        positions = pick_mutation_positions(
            len(seq),
            n_mut,
            rng,
            fixed_positions=cfg.fixed_positions,
        )
        for pos in positions:
            current = seq[pos]
            choices = [aa for aa in cfg.allowed_residues if aa != current]
            if not choices:
                return None
            seq[pos] = rng.choice(choices)
        candidate = "".join(seq)
        if not is_valid_sequence(
            candidate,
            min_length=cfg.min_length,
            max_length=cfg.max_length,
            allowed_residues=cfg.allowed_residues,
            fixed_positions=cfg.fixed_positions_dict,
        ):
            return None
        return candidate


class _MutationConfig:
    __slots__ = (
        "allowed_residues",
        "fixed_positions",
        "fixed_positions_dict",
        "max_length",
        "min_length",
    )

    def __init__(
        self,
        *,
        min_length: int,
        max_length: int,
        allowed_residues: frozenset[str],
        fixed_positions: dict[int, str] | None,
    ) -> None:
        self.min_length = min_length
        self.max_length = max_length
        self.allowed_residues = allowed_residues
        self.fixed_positions_dict = fixed_positions or {}
        self.fixed_positions = frozenset(self.fixed_positions_dict.keys())


def _parse_constraints(
    constraints: dict[str, object] | None,
    *,
    default_allowed: frozenset[str],
) -> _MutationConfig:
    if constraints is None:
        return _MutationConfig(
            min_length=MIN_PEPTIDE_LENGTH,
            max_length=MAX_PEPTIDE_LENGTH,
            allowed_residues=default_allowed,
            fixed_positions=None,
        )
    min_raw = constraints.get("min_length", MIN_PEPTIDE_LENGTH)
    max_raw = constraints.get("max_length", MAX_PEPTIDE_LENGTH)
    if not isinstance(min_raw, int) or not isinstance(max_raw, int):
        raise TypeError("min_length and max_length must be int")
    min_length = min_raw
    max_length = max_raw
    allowed_raw = constraints.get("allowed_residues")
    allowed: frozenset[str]
    if allowed_raw is None:
        allowed = default_allowed
    elif isinstance(allowed_raw, (list, tuple, set, frozenset)):
        allowed = normalize_allowed_residues(allowed_raw)
    else:
        raise TypeError("allowed_residues must be a sequence of one-letter codes")
    fixed_raw = constraints.get("fixed_positions")
    fixed: dict[int, str] | None = None
    if fixed_raw is not None:
        if not isinstance(fixed_raw, dict):
            raise TypeError("fixed_positions must be a dict[int, str]")
        fixed = {int(k): str(v).upper() for k, v in fixed_raw.items()}
    return _MutationConfig(
        min_length=min_length,
        max_length=max_length,
        allowed_residues=allowed,
        fixed_positions=fixed,
    )


def _normalize_seeds(
    seed_sequences: tuple[str, ...] | None,
    cfg: _MutationConfig,
    rng: random.Random,
) -> list[str]:
    if seed_sequences:
        seeds = [s.upper() for s in seed_sequences]
    else:
        length = rng.randint(cfg.min_length, min(cfg.max_length, 12))
        seeds = [
            "".join(rng.choice(sorted(cfg.allowed_residues)) for _ in range(length))
        ]
    for seq in seeds:
        if not is_valid_sequence(
            seq,
            min_length=cfg.min_length,
            max_length=cfg.max_length,
            allowed_residues=cfg.allowed_residues,
            fixed_positions=cfg.fixed_positions_dict,
        ):
            raise ValueError(f"invalid seed sequence: {seq!r}")
    return seeds

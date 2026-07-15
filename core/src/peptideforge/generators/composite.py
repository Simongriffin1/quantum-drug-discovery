"""Composite generator: ESM-2 masked sampling + mutation fallback."""

from __future__ import annotations

from peptideforge.contracts.models import Candidates, PeptideCandidate
from peptideforge.generators.deps import ESM2UnavailableError
from peptideforge.generators.esm2 import ESM2Generator
from peptideforge.generators.filters import is_diverse_enough
from peptideforge.generators.mutation import MutationGenerator


class PeptideGenerator:
    """Primary Generator: ESM-2 proposals plus mutation fill for diversity.

    When fair-esm is unavailable, falls back to mutation-only (still valid for CI
    if ``require_esm=False``). Never invents sequences without a seed or RNG path.
    """

    def __init__(
        self,
        *,
        require_esm: bool = False,
        esm_fraction: float = 0.5,
        mutation: MutationGenerator | None = None,
        esm: ESM2Generator | None = None,
    ) -> None:
        if not 0.0 <= esm_fraction <= 1.0:
            raise ValueError("esm_fraction must be in [0, 1]")
        self.require_esm = require_esm
        self.esm_fraction = esm_fraction
        self.mutation = mutation or MutationGenerator(generation_method="mutation")
        self.esm = esm or ESM2Generator()

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
        n_esm = int(round(n * self.esm_fraction)) if seed_sequences else 0
        n_mut = n - n_esm
        items: list[PeptideCandidate] = []

        if n_esm > 0:
            try:
                esm_batch = self.esm.propose(
                    n=n_esm,
                    seed_sequences=seed_sequences,
                    seed=seed,
                    constraints=constraints,
                )
                items.extend(esm_batch.items)
            except ESM2UnavailableError:
                if self.require_esm:
                    raise
                n_mut = n

        if n_mut > 0:
            mut_batch = self.mutation.propose(
                n=n_mut,
                seed_sequences=seed_sequences,
                seed=seed,
                constraints=constraints,
            )
            pool = [c.sequence for c in items]
            for cand in mut_batch.items:
                if is_diverse_enough(cand.sequence, pool, max_identity=0.95):
                    items.append(cand)
                    pool.append(cand.sequence)

        # Top up if ESM path failed partially
        attempts = 0
        while len(items) < n and attempts < n * 100:
            attempts += 1
            extra = self.mutation.propose(
                n=1,
                seed_sequences=seed_sequences,
                seed=None if seed is None else seed + attempts,
                constraints=constraints,
            )
            cand = extra.items[0]
            if is_diverse_enough(cand.sequence, [c.sequence for c in items], max_identity=0.95):
                items.append(cand)

        if len(items) < n:
            raise RuntimeError(f"PeptideGenerator produced {len(items)}/{n} candidates")

        # Deduplicate exact sequences while preserving order
        seen: set[str] = set()
        unique: list[PeptideCandidate] = []
        for cand in items:
            if cand.sequence in seen:
                continue
            seen.add(cand.sequence)
            unique.append(cand)
        if len(unique) < n:
            raise RuntimeError("insufficient unique candidates after deduplication")

        return Candidates(items=tuple(unique[:n]), seed=seed)

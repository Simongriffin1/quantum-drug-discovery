"""Tests for P6 sequence generators."""

from __future__ import annotations

import pytest

from peptideforge.contracts.models import PeptideCandidate
from peptideforge.contracts.protocols import Generator
from peptideforge.generators import (
    ESM2Generator,
    ESM2UnavailableError,
    MutationGenerator,
    PeptideGenerator,
)
from peptideforge.generators.filters import sequence_identity


@pytest.fixture
def seed_epitope() -> str:
    return "LLFGYPVYV"


def test_mutation_generator_honors_protocol(seed_epitope: str) -> None:
    gen = MutationGenerator(generation_method="synthetic_mutation")
    assert isinstance(gen, Generator)
    batch = gen.propose(n=8, seed_sequences=(seed_epitope,), seed=42)
    assert len(batch.items) == 8
    for cand in batch.items:
        assert len(cand.sequence) == len(seed_epitope)
        assert cand.generation_method == "synthetic_mutation"


def test_mutation_diversity(seed_epitope: str) -> None:
    gen = MutationGenerator(max_identity=0.9)
    batch = gen.propose(n=10, seed_sequences=(seed_epitope,), seed=0)
    seqs = [c.sequence for c in batch.items]
    assert len(set(seqs)) == len(seqs)
    for seq in seqs:
        assert 0.0 < sequence_identity(seq, seed_epitope) < 1.0


def test_mutation_respects_fixed_positions(seed_epitope: str) -> None:
    gen = MutationGenerator()
    batch = gen.propose(
        n=5,
        seed_sequences=(seed_epitope,),
        seed=1,
        constraints={"fixed_positions": {0: "L", 4: "Y"}},
    )
    for cand in batch.items:
        assert cand.sequence[0] == "L"
        assert cand.sequence[4] == "Y"


def test_peptide_generator_mutation_fallback(seed_epitope: str) -> None:
    gen = PeptideGenerator(require_esm=False, esm_fraction=0.5)
    batch = gen.propose(n=6, seed_sequences=(seed_epitope,), seed=7)
    assert len(batch.items) == 6


def test_esm2_raises_when_unavailable(seed_epitope: str) -> None:
    gen = ESM2Generator()
    try:
        import esm  # noqa: F401
    except ImportError:
        with pytest.raises(ESM2UnavailableError):
            gen.propose(n=2, seed_sequences=(seed_epitope,), seed=0)
    else:
        pytest.skip("fair-esm installed; unavailable-path test skipped")


def test_invalid_seed_rejected() -> None:
    gen = MutationGenerator()
    with pytest.raises(ValueError, match="invalid seed"):
        gen.propose(n=1, seed_sequences=("ACDE",), seed=0)

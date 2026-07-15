"""Tests for P6 structure predictors."""

from __future__ import annotations

from pathlib import Path

import pytest

from peptideforge.contracts.models import PeptideCandidate
from peptideforge.contracts.protocols import StructurePredictor
from peptideforge.structure import (
    Boltz2StructurePredictor,
    Boltz2UnavailableError,
    FixtureStructurePredictor,
    FoldCache,
    require_boltz_cli,
)
from peptideforge.structure.pdb_utils import extract_chain_sequence, mean_plddt_from_pdb

STRUCTURES = Path(__file__).resolve().parents[3] / "benchmarks" / "fixtures" / "structures"
INTERFACE = STRUCTURES / "PP003_1BD2_interface.pdb"


@pytest.fixture
def epitope_candidate() -> PeptideCandidate:
    return PeptideCandidate(sequence="LLFGYPVYV", generation_method="synthetic_fold")


def test_fixture_predictor_honors_protocol(epitope_candidate: PeptideCandidate) -> None:
    if not INTERFACE.is_file():
        pytest.skip(f"missing fixture {INTERFACE}")
    predictor = FixtureStructurePredictor(structures_dir=STRUCTURES)
    assert isinstance(predictor, StructurePredictor)
    folded = predictor.fold(
        epitope_candidate,
        target_id="1BD2",
        target_structure=str(INTERFACE),
    )
    assert folded.pdb_text is not None
    assert folded.confidence == 1.0
    assert folded.fold_method == "fixture_interface_pdb"
    assert "ATOM" in folded.pdb_text


def test_fixture_cache_hit(epitope_candidate: PeptideCandidate, tmp_path: Path) -> None:
    if not INTERFACE.is_file():
        pytest.skip(f"missing fixture {INTERFACE}")
    cache = FoldCache(tmp_path / "fold_cache")
    predictor = FixtureStructurePredictor(structures_dir=STRUCTURES, cache=cache)
    first = predictor.fold(
        epitope_candidate,
        target_id="1BD2",
        target_structure=str(INTERFACE),
    )
    second = predictor.fold(
        epitope_candidate,
        target_id="1BD2",
        target_structure=str(INTERFACE),
    )
    assert first.cache_key == second.cache_key
    assert first.pdb_path == second.pdb_path


def test_fixture_rejects_mismatched_sequence(epitope_candidate: PeptideCandidate) -> None:
    if not INTERFACE.is_file():
        pytest.skip(f"missing fixture {INTERFACE}")
    predictor = FixtureStructurePredictor(structures_dir=STRUCTURES)
    bad = PeptideCandidate(sequence="AAAAAAAAA", generation_method="synthetic_fold")
    with pytest.raises(ValueError, match="does not match manifest epitope"):
        predictor.fold(bad, target_id="1BD2", target_structure=str(INTERFACE))


def test_boltz2_raises_when_cli_missing(epitope_candidate: PeptideCandidate) -> None:
    try:
        require_boltz_cli()
    except Boltz2UnavailableError:
        pass
    else:
        pytest.skip("boltz CLI present; missing-deps test skipped")
    if not INTERFACE.is_file():
        pytest.skip(f"missing fixture {INTERFACE}")
    predictor = Boltz2StructurePredictor()
    with pytest.raises(Boltz2UnavailableError):
        predictor.fold(
            epitope_candidate,
            target_id="1BD2",
            target_structure=str(INTERFACE),
        )


def test_pdb_sequence_extraction() -> None:
    if not INTERFACE.is_file():
        pytest.skip(f"missing fixture {INTERFACE}")
    text = INTERFACE.read_text(encoding="utf-8")
    pep = extract_chain_sequence(text, "C")
    assert pep == "LLFGYPVYV"
    conf = mean_plddt_from_pdb(text, chain_id="C")
    assert 0.0 <= conf <= 1.0

"""Golden / contract tests for the OpenMM physics oracle."""

from __future__ import annotations

from pathlib import Path

import pytest

from peptideforge.contracts.models import ComplexStructure, OracleTier, PeptideCandidate
from peptideforge.oracles.costs import CostCapExceededError
from peptideforge.oracles.energy_conservation import run_energy_conservation_check
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle
from peptideforge.oracles.openmm_utils import OpenMMUnavailableError, require_openmm

STRUCTURES = Path(__file__).resolve().parents[3] / "benchmarks" / "fixtures" / "structures"
INTERFACE = STRUCTURES / "PP003_1BD2_interface.pdb"


@pytest.fixture(scope="module")
def openmm_ok() -> None:
    try:
        require_openmm()
    except OpenMMUnavailableError as exc:
        pytest.skip(str(exc))


@pytest.mark.golden
def test_energy_conservation_alanine_dipeptide(openmm_ok: None) -> None:
    result = run_energy_conservation_check(
        n_steps=500,
        timestep_fs=1.0,
        tolerance_kcal=5.0,
        platform="CPU",
        seed=0,
    )
    assert result.passed, (
        f"energy drift too large: max_abs_drift={result.max_abs_drift_kcal:.3f} "
        f"tol={result.tolerance_kcal}"
    )


@pytest.mark.golden
def test_mm_gbsa_oracle_contract_and_default_tier(openmm_ok: None) -> None:
    if not INTERFACE.is_file():
        pytest.skip(f"missing interface fixture {INTERFACE}")

    oracle = OpenMMPhysicsOracle(
        OpenMMOracleConfig(platform="CPU", peptide_chain_ids=("C",), seed=0)
    )
    cand = PeptideCandidate(sequence="LLFGYPVYV", generation_method="synthetic_test")
    complex_structure = ComplexStructure(
        candidate_id=cand.candidate_id,
        target_id="1BD2",
        sequence=cand.sequence,
        pdb_path=str(INTERFACE.resolve()),
        confidence=1.0,
        fold_method="experimental_pdb_fixture",
    )
    # Default tier = cheapest (docking)
    result = oracle.evaluate(complex_structure)
    assert result.tier == OracleTier.DOCKING
    assert result.unit == "kcal/mol"
    assert result.cost_estimate == pytest.approx(0.05)
    assert result.value == result.value  # finite

    mm = oracle.evaluate(complex_structure, tier=OracleTier.MM_GBSA)
    assert mm.tier == OracleTier.MM_GBSA
    assert mm.cost_estimate == pytest.approx(1.0)


@pytest.mark.golden
def test_cost_cap_enforced_on_evaluate(openmm_ok: None) -> None:
    if not INTERFACE.is_file():
        pytest.skip(f"missing interface fixture {INTERFACE}")
    oracle = OpenMMPhysicsOracle(
        OpenMMOracleConfig(platform="CPU", peptide_chain_ids=("C",), seed=0)
    )
    cand = PeptideCandidate(sequence="LLFGYPVYV", generation_method="synthetic_test")
    complex_structure = ComplexStructure(
        candidate_id=cand.candidate_id,
        target_id="1BD2",
        sequence=cand.sequence,
        pdb_path=str(INTERFACE.resolve()),
        confidence=1.0,
        fold_method="experimental_pdb_fixture",
    )
    with pytest.raises(CostCapExceededError):
        oracle.evaluate(complex_structure, tier=OracleTier.MM_GBSA, cost_cap=0.1)


@pytest.mark.golden
def test_missing_openmm_fails_loud(monkeypatch: pytest.MonkeyPatch) -> None:
    import peptideforge.oracles.mm_gbsa as mm_mod

    def _boom() -> tuple[object, object, object]:
        raise OpenMMUnavailableError("forced missing openmm")

    monkeypatch.setattr(mm_mod, "require_openmm", _boom)
    with pytest.raises(OpenMMUnavailableError):
        OpenMMPhysicsOracle()

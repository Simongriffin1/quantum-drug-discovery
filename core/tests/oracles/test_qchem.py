"""Golden / contract tests for the quantum chemistry oracle (P4)."""

from __future__ import annotations

import pytest

from peptideforge.contracts.models import OracleTier
from peptideforge.oracles.costs import CostCapExceededError
from peptideforge.oracles.qchem import (
    H2_BOND_ANGSTROM,
    HARTREE_TOLERANCE,
    HE_HF_STO3G_LITERATURE_HARTREE,
    PennyLaneUnavailableError,
    PySCFUnavailableError,
    QuantumChemistryOracle,
    QuantumOracleConfig,
    h2_fragment_structure,
    he_atom_structure,
)
from peptideforge.oracles.qchem_deps import require_pennylane, require_pyscf


@pytest.fixture(scope="module")
def pyscf_ok() -> None:
    try:
        require_pyscf()
    except PySCFUnavailableError as exc:
        pytest.skip(str(exc))


@pytest.fixture(scope="module")
def pennylane_ok(pyscf_ok: None) -> None:
    try:
        require_pennylane()
    except PennyLaneUnavailableError as exc:
        pytest.skip(str(exc))


@pytest.mark.golden
def test_he_hf_matches_literature(pyscf_ok: None) -> None:
    """PySCF HF/STO-3G He atom matches a published reference within tolerance."""
    oracle = QuantumChemistryOracle(QuantumOracleConfig(classical_method="hf", seed=0))
    result = oracle.evaluate(he_atom_structure(), tier=OracleTier.QCHEM_CLASSICAL)
    assert result.tier == OracleTier.QCHEM_CLASSICAL
    assert result.unit == "hartree"
    assert abs(result.value - HE_HF_STO3G_LITERATURE_HARTREE) <= HARTREE_TOLERANCE


@pytest.mark.golden
def test_h2_classical_fci_finite(pyscf_ok: None) -> None:
    oracle = QuantumChemistryOracle(QuantumOracleConfig(classical_method="fci", seed=0))
    result = oracle.evaluate(
        h2_fragment_structure(bond_angstrom=H2_BOND_ANGSTROM),
        tier=OracleTier.QCHEM_CLASSICAL,
    )
    # STO-3G FCI H2 ≈ −1.137 hartree at ~0.74 Å
    assert result.value < -1.0
    assert result.value > -1.3


@pytest.mark.golden
def test_h2_vqe_matches_classical_baseline(pennylane_ok: None) -> None:
    """Quantum gate: VQE within 0.01 hartree of PySCF classical; baseline always present."""
    oracle = QuantumChemistryOracle(
        QuantumOracleConfig(classical_method="fci", seed=0, vqe_steps=120, vqe_stepsize=0.2)
    )
    result = oracle.evaluate(
        h2_fragment_structure(bond_angstrom=H2_BOND_ANGSTROM),
        tier=OracleTier.QCHEM_VQE,
    )
    assert result.tier == OracleTier.QCHEM_VQE
    assert result.classical_baseline is not None
    assert result.unit == "hartree"
    delta = abs(result.value - result.classical_baseline)
    assert delta <= HARTREE_TOLERANCE, (
        f"VQE–classical |ΔE|={delta} > {HARTREE_TOLERANCE}; "
        f"VQE={result.value}, classical={result.classical_baseline}"
    )
    assert result.metadata.get("vqe_matched_classical") is True
    assert result.metadata.get("use_hardware") is False


@pytest.mark.golden
def test_vqe_never_without_classical_baseline_contract(pennylane_ok: None) -> None:
    oracle = QuantumChemistryOracle(QuantumOracleConfig(classical_method="fci", seed=0))
    result = oracle.evaluate(h2_fragment_structure(), tier=OracleTier.QCHEM_VQE)
    # Contract-level: QCHEM_VQE without classical_baseline would fail validation
    dumped = result.model_dump()
    assert dumped["classical_baseline"] is not None


@pytest.mark.golden
def test_qchem_cost_cap(pyscf_ok: None) -> None:
    oracle = QuantumChemistryOracle(QuantumOracleConfig(classical_method="hf", seed=0))
    with pytest.raises(CostCapExceededError):
        oracle.evaluate(he_atom_structure(), tier=OracleTier.QCHEM_CLASSICAL, cost_cap=1.0)


@pytest.mark.golden
def test_default_tier_is_classical(pyscf_ok: None) -> None:
    oracle = QuantumChemistryOracle(QuantumOracleConfig(classical_method="hf", seed=0))
    result = oracle.evaluate(he_atom_structure())
    assert result.tier == OracleTier.QCHEM_CLASSICAL


@pytest.mark.golden
def test_hardware_flag_fails_loud(pyscf_ok: None) -> None:
    with pytest.raises(NotImplementedError):
        QuantumChemistryOracle(QuantumOracleConfig(use_hardware=True))


@pytest.mark.golden
def test_missing_pyscf_fails_loud(monkeypatch: pytest.MonkeyPatch) -> None:
    import peptideforge.oracles.qchem as qmod

    def _boom() -> object:
        raise PySCFUnavailableError("forced missing pyscf")

    monkeypatch.setattr(qmod, "require_pyscf", _boom)
    with pytest.raises(PySCFUnavailableError):
        QuantumChemistryOracle()

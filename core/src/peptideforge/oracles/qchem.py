"""Quantum chemistry oracle: PySCF classical + optional PennyLane VQE.

Physical rationale
------------------
Quantum chem is a **precision scalpel** for small fragments / active-site models —
not a blanket advantage claim over classical MD. Every VQE result ships with its
PySCF classical baseline (quantum gate). Hardware backends are config-gated for
later without changing callers; default is a local simulator.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any
from uuid import UUID, uuid4

from peptideforge.contracts.models import (
    ComplexStructure,
    OracleResult,
    OracleTier,
    PeptideCandidate,
    Provenance,
)
from peptideforge.oracles.costs import enforce_cost_cap, resolve_tier
from peptideforge.oracles.qchem_deps import (
    PennyLaneUnavailableError,
    PySCFUnavailableError,
    require_pennylane,
    require_pyscf,
)

# H2 / He golden-test constants (documented for ACCEPTANCE.md).
H2_BOND_ANGSTROM = 0.7414
# STO-3G FCI H2 energy at ~0.74 Å is ~−1.137 hartree (basis-limited FCI, not
# the CBS exact −1.174). Tests compare VQE vs *this run's* PySCF baseline.
HE_HF_STO3G_LITERATURE_HARTREE = -2.80774745  # He atom, HF/STO-3G reference
HARTREE_TOLERANCE = 0.01


@dataclass(frozen=True)
class QuantumOracleConfig:
    """Quantum-chem hyperparameters. ``use_hardware`` reserved for later backends."""

    basis: str = "sto-3g"
    # Classical methods: hf | ccsd | fci | dft
    classical_method: str = "fci"
    dft_xc: str = "b3lypg"  # PySCF B3LYP
    vqe_layers: int = 1
    vqe_steps: int = 120
    vqe_stepsize: float = 0.2
    seed: int = 0
    # Simulator name for PennyLane; ignored when use_hardware=True (not wired yet).
    simulator_device: str = "default.qubit"
    use_hardware: bool = False
    max_atoms: int = 6  # hard cap — fail loud on oversized fragments


def geometry_from_complex(
    complex_structure: ComplexStructure,
) -> tuple[list[str], list[list[float]]]:
    """Extract (symbols, coords Å) from XYZ text/path on ComplexStructure.

    Quantum fragments are represented as XYZ in ``pdb_text`` / ``pdb_path``
    (not peptide PDB). Fail loud if geometry cannot be parsed.
    """
    text = complex_structure.pdb_text
    if text is None and complex_structure.pdb_path:
        from pathlib import Path

        path = Path(complex_structure.pdb_path)
        if not path.is_file():
            raise FileNotFoundError(f"fragment XYZ/PDB not found: {path}")
        text = path.read_text(encoding="utf-8")
    if text is None:
        raise ValueError(
            "QuantumChemistryOracle requires XYZ in pdb_text/pdb_path "
            "(fail loud — no fabricated geometry)"
        )
    return _parse_xyz(text)


def _parse_xyz(text: str) -> tuple[list[str], list[list[float]]]:
    lines = [ln.strip() for ln in text.strip().splitlines() if ln.strip()]
    if len(lines) < 3:
        raise ValueError("XYZ text too short")
    try:
        n = int(lines[0].split()[0])
    except ValueError as exc:
        raise ValueError("XYZ first line must be atom count") from exc
    body = lines[2 : 2 + n]
    if len(body) < n:
        raise ValueError(f"XYZ declares {n} atoms but found {len(body)} coordinate lines")
    atoms: list[str] = []
    coords: list[list[float]] = []
    for line in body:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"bad XYZ atom line: {line}")
        atoms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return atoms, coords


def h2_fragment_structure(
    *,
    bond_angstrom: float = H2_BOND_ANGSTROM,
    candidate_id: UUID | None = None,
) -> ComplexStructure:
    """Build a ComplexStructure carrying H2 geometry for quantum golden tests."""
    half = bond_angstrom / 2.0
    cand = PeptideCandidate(
        candidate_id=candidate_id or uuid4(),
        sequence="AAAAA",  # placeholder — fragment geoms live in XYZ payload
        generation_method="synthetic_qchem_h2",
        metadata={"fragment": "H2", "bond_angstrom": bond_angstrom},
    )
    coords = [[0.0, 0.0, -half], [0.0, 0.0, half]]
    xyz = (
        f"2\nH2 bond={bond_angstrom} A synthetic_qchem\n"
        f"H {coords[0][0]:.6f} {coords[0][1]:.6f} {coords[0][2]:.6f}\n"
        f"H {coords[1][0]:.6f} {coords[1][1]:.6f} {coords[1][2]:.6f}\n"
    )
    return ComplexStructure(
        candidate_id=cand.candidate_id,
        target_id="synthetic_h2",
        sequence=cand.sequence,
        pdb_text=xyz,
        confidence=1.0,
        fold_method="synthetic_fragment_xyz",
    )


def he_atom_structure(*, candidate_id: UUID | None = None) -> ComplexStructure:
    """He atom for HF/STO-3G literature comparison golden test."""
    cand = PeptideCandidate(
        candidate_id=candidate_id or uuid4(),
        sequence="AAAAA",
        generation_method="synthetic_qchem_he",
        metadata={"fragment": "He"},
    )
    xyz = "1\nHe atom synthetic_qchem\nHe 0.0 0.0 0.0\n"
    return ComplexStructure(
        candidate_id=cand.candidate_id,
        target_id="synthetic_he",
        sequence=cand.sequence,
        pdb_text=xyz,
        confidence=1.0,
        fold_method="synthetic_fragment_xyz",
    )


class QuantumChemistryOracle:
    """Oracle for tiny fragments: classical PySCF + optional simulator VQE."""

    AVAILABLE_TIERS: tuple[OracleTier, ...] = (
        OracleTier.QCHEM_CLASSICAL,
        OracleTier.QCHEM_VQE,
    )

    def __init__(self, config: QuantumOracleConfig | None = None) -> None:
        # Classical path always required (VQE needs baseline).
        require_pyscf()
        self.config = config or QuantumOracleConfig()
        if self.config.use_hardware:
            raise NotImplementedError(
                "use_hardware=True is reserved for a future backend; "
                "simulator is the only supported path (fail loud)."
            )

    def evaluate(
        self,
        complex_structure: ComplexStructure,
        *,
        tier: OracleTier | None = None,
        cost_cap: float | None = None,
    ) -> OracleResult:
        resolved = resolve_tier(tier, available=self.AVAILABLE_TIERS)
        cost = enforce_cost_cap(resolved, cost_cap)
        atoms, coords = geometry_from_complex(complex_structure)
        if len(atoms) > self.config.max_atoms:
            raise ValueError(
                f"fragment has {len(atoms)} atoms > max_atoms={self.config.max_atoms}; "
                "quantum tier is for small active-site models only"
            )

        classical_e, classical_meta = self._classical_energy(atoms, coords)

        if resolved == OracleTier.QCHEM_CLASSICAL:
            pyscf = require_pyscf()
            return OracleResult(
                candidate_id=complex_structure.candidate_id,
                complex_id=complex_structure.complex_id,
                value=classical_e,
                uncertainty=0.0,
                cost_estimate=cost,
                tier=OracleTier.QCHEM_CLASSICAL,
                unit="hartree",
                classical_baseline=classical_e,
                metadata={
                    **classical_meta,
                    "atoms": atoms,
                    "coordinates_angstrom": coords,
                },
                provenance=Provenance(
                    tool_versions={"pyscf": getattr(pyscf, "__version__", "unknown")},
                    config_hash=str(self.config.seed),
                ),
            )

        # VQE path — classical baseline ALWAYS attached
        vqe_e, vqe_meta = self._vqe_energy(atoms, coords, classical_baseline=classical_e)
        matched = abs(vqe_e - classical_e) <= HARTREE_TOLERANCE
        qml, _ = require_pennylane()
        pyscf = require_pyscf()
        return OracleResult(
            candidate_id=complex_structure.candidate_id,
            complex_id=complex_structure.complex_id,
            value=vqe_e,
            uncertainty=abs(vqe_e - classical_e),
            cost_estimate=cost,
            tier=OracleTier.QCHEM_VQE,
            unit="hartree",
            classical_baseline=classical_e,
            metadata={
                **classical_meta,
                **vqe_meta,
                "vqe_matched_classical": matched,
                "vqe_minus_classical_hartree": vqe_e - classical_e,
                "tolerance_hartree": HARTREE_TOLERANCE,
                "atoms": atoms,
                "coordinates_angstrom": coords,
                "simulator_device": self.config.simulator_device,
                "use_hardware": False,
            },
            provenance=Provenance(
                tool_versions={
                    "pyscf": getattr(pyscf, "__version__", "unknown"),
                    "pennylane": getattr(qml, "__version__", "unknown"),
                },
                config_hash=str(self.config.seed),
            ),
        )

    def _classical_energy(
        self, atoms: list[str], coords: list[list[float]]
    ) -> tuple[float, dict[str, Any]]:
        require_pyscf()
        from pyscf import cc, dft, fci, gto, scf

        mol = gto.M(
            atom=[[sym, tuple(xyz)] for sym, xyz in zip(atoms, coords, strict=True)],
            basis=self.config.basis,
            unit="Angstrom",
            verbose=0,
        )
        method = self.config.classical_method.lower()
        if method == "hf":
            mf = scf.RHF(mol).run()
            energy = float(mf.e_tot)
        elif method == "dft":
            mf = dft.RKS(mol)
            mf.xc = self.config.dft_xc
            mf = mf.run()
            energy = float(mf.e_tot)
        elif method == "ccsd":
            mf = scf.RHF(mol).run()
            mycc = cc.CCSD(mf).run()
            energy = float(mycc.e_tot)
        elif method == "fci":
            mf = scf.RHF(mol).run()
            cisolver = fci.FCI(mol, mf.mo_coeff)
            e_fci, _ = cisolver.kernel()
            # FCI returns correlation+HF electronic; include nuclear repulsion
            energy = float(e_fci)
        else:
            raise ValueError(f"unsupported classical_method={method}")

        return energy, {
            "classical_method": method,
            "basis": self.config.basis,
            "dft_xc": self.config.dft_xc if method == "dft" else None,
            "n_electrons": int(mol.nelectron),
            "n_atoms": len(atoms),
        }

    def _vqe_energy(
        self,
        atoms: list[str],
        coords: list[list[float]],
        *,
        classical_baseline: float,
    ) -> tuple[float, dict[str, Any]]:
        """Run a simulator VQE; classical_baseline is required (caller must pass it)."""
        _ = classical_baseline  # documented requirement on the call signature
        qml, qchem = require_pennylane()
        import numpy as np

        symbols = atoms
        geometry = np.array(coords, dtype=float)
        try:
            hamiltonian, qubits = qchem.molecular_hamiltonian(
                symbols,
                geometry,
                charge=0,
                mult=1,
                basis=self.config.basis,
                unit="angstrom",
            )
        except TypeError:
            # Older pennylane.qchem API without unit= — convert Å → bohr
            bohr = geometry / 0.52917721092
            hamiltonian, qubits = qchem.molecular_hamiltonian(
                symbols,
                bohr,
                charge=0,
                mult=1,
                basis=self.config.basis,
            )

        n_electrons = _valence_electrons(atoms)
        hf_state = qchem.hf_state(n_electrons, qubits)
        singles, doubles = qchem.excitations(n_electrons, qubits)

        if not singles and not doubles:
            # Fallback RY layers when no excitations (e.g. trivial systems)
            dev = qml.device(self.config.simulator_device, wires=qubits)

            @qml.qnode(dev)  # type: ignore[untyped-decorator]
            def circuit_ry(params: Any) -> Any:
                qml.BasisState(np.array(hf_state), wires=range(qubits))
                for i in range(qubits):
                    qml.RY(params[i], wires=i)
                return qml.expval(hamiltonian)

            rng = np.random.default_rng(self.config.seed)
            params = qml.numpy.array(rng.normal(0.0, 0.1, size=(qubits,)), requires_grad=True)
            opt = qml.AdamOptimizer(stepsize=self.config.vqe_stepsize)
            energy = float(circuit_ry(params))
            for _ in range(self.config.vqe_steps):
                params, energy = opt.step_and_cost(circuit_ry, params)
            return float(energy), {
                "vqe_ansatz": "hf_plus_ry",
                "n_qubits": int(qubits),
                "n_steps": self.config.vqe_steps,
            }

        wires = range(qubits)
        dev = qml.device(self.config.simulator_device, wires=qubits)

        @qml.qnode(dev)  # type: ignore[untyped-decorator]
        def circuit_uccsd(weights: Any) -> Any:
            qml.AllSinglesDoubles(
                weights,
                wires,
                hf_state=np.array(hf_state),
                singles=singles,
                doubles=doubles,
            )
            return qml.expval(hamiltonian)

        n_params = len(singles) + len(doubles)
        # Start near HF (zero amplitudes) and mark trainable for autodiff.
        weights = qml.numpy.zeros(n_params, requires_grad=True)
        opt = qml.AdamOptimizer(stepsize=self.config.vqe_stepsize)
        energy = float(circuit_uccsd(weights))
        for _ in range(self.config.vqe_steps):
            weights, energy = opt.step_and_cost(circuit_uccsd, weights)

        return float(energy), {
            "vqe_ansatz": "AllSinglesDoubles",
            "n_qubits": int(qubits),
            "n_params": n_params,
            "n_steps": self.config.vqe_steps,
            "n_singles": len(singles),
            "n_doubles": len(doubles),
        }


def _valence_electrons(atoms: list[str]) -> int:
    """Minimal electron count for H/He/C/N/O/F fragments used in golden tests."""
    table = {
        "H": 1,
        "HE": 2,
        "C": 4,
        "N": 5,
        "O": 6,
        "F": 7,
        "LI": 1,
        "BE": 2,
        "B": 3,
    }
    total = 0
    for a in atoms:
        key = a.upper()
        if key not in table:
            raise ValueError(
                f"atom {a} not in minimal valence table for VQE path; "
                "extend _valence_electrons or use classical-only tier"
            )
        total += table[key]
    return total


# Re-export errors for callers
__all__ = [
    "H2_BOND_ANGSTROM",
    "HARTREE_TOLERANCE",
    "HE_HF_STO3G_LITERATURE_HARTREE",
    "PennyLaneUnavailableError",
    "PySCFUnavailableError",
    "QuantumChemistryOracle",
    "QuantumOracleConfig",
    "geometry_from_complex",
    "h2_fragment_structure",
    "he_atom_structure",
]

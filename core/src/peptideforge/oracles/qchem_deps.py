"""Optional quantum chemistry dependency guards — fail loud, never fabricate."""

from __future__ import annotations

from typing import Any


class PySCFUnavailableError(ImportError):
    """Raised when PySCF cannot be imported."""


class PennyLaneUnavailableError(ImportError):
    """Raised when PennyLane (simulator VQE path) cannot be imported."""


def require_pyscf() -> Any:
    """Import ``pyscf`` or raise loudly."""
    try:
        import pyscf
    except ImportError as exc:
        raise PySCFUnavailableError(
            "PySCF is required for the classical quantum-chemistry oracle tier but is "
            "not installed. Install with: pip install pyscf. "
            "Refusing to fabricate energies."
        ) from exc
    return pyscf


def require_pennylane() -> tuple[Any, Any]:
    """Import PennyLane + qchem helpers or raise loudly."""
    try:
        import pennylane as qml
        from pennylane import qchem
    except ImportError as exc:
        raise PennyLaneUnavailableError(
            "PennyLane (with pennylane.qchem) is required for the VQE oracle tier but is "
            "not installed. Install with: pip install pennylane pennylane-qchem. "
            "Refusing to fabricate VQE energies."
        ) from exc
    return qml, qchem

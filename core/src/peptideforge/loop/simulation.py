"""Simulation-mode ground truth + synthetic fold/oracle for CI (P9).

All labels come from an explicit synthetic_* function — never fabricated as
real physics. Used to exercise the full DBTL loop without OpenMM/Boltz/ESM.
"""

from __future__ import annotations

from peptideforge.contracts.models import (
    ComplexStructure,
    OracleResult,
    OracleTier,
    PeptideCandidate,
    Provenance,
)
from peptideforge.surrogate.ensemble import synthetic_physics_label


class SimulationOracle:
    """Tiered-compatible synthetic oracle wrapping ``synthetic_physics_label``."""

    def __init__(self, *, noise: float = 0.05, cost: float = 1.0) -> None:
        self.noise = noise
        self.cost = cost

    def evaluate(
        self,
        complex_structure: ComplexStructure,
        *,
        tier: OracleTier | None = None,
        cost_cap: float | None = None,
    ) -> OracleResult:
        del tier  # simulation uses a single synthetic tier
        if cost_cap is not None and self.cost > cost_cap:
            from peptideforge.oracles.costs import CostCapExceededError

            raise CostCapExceededError(
                f"simulation oracle cost {self.cost} exceeds cost_cap {cost_cap}"
            )
        value = synthetic_physics_label(
            complex_structure.sequence,
            noise=self.noise,
            seed=0,
        )
        return OracleResult(
            candidate_id=complex_structure.candidate_id,
            complex_id=complex_structure.complex_id,
            value=value,
            uncertainty=self.noise,
            cost_estimate=self.cost,
            tier=OracleTier.SYNTHETIC,
            unit="synthetic_score",
            metadata={"method": "synthetic_physics_label"},
            provenance=Provenance(
                tool_versions={"simulation_oracle": "0.1.0"},
                data_version="synthetic_v1",
            ),
        )


class SyntheticStructurePredictor:
    """Minimal valid ComplexStructure for simulation mode (no fabricated physics)."""

    fold_method = "synthetic_fold"

    def fold(
        self,
        candidate: PeptideCandidate,
        *,
        target_id: str,
        target_structure: str,
        seed: int | None = None,
    ) -> ComplexStructure:
        del target_structure, seed
        # Minimal PDB text — plumbing only, not a predicted pose
        pdb_text = (
            "HEADER    SYNTHETIC PLUMBING STRUCTURE\n"
            f"TITLE     SIMULATION FOLD FOR {candidate.sequence}\n"
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 50.00           C\n"
            "END\n"
        )
        return ComplexStructure(
            candidate_id=candidate.candidate_id,
            target_id=target_id,
            sequence=candidate.sequence,
            pdb_text=pdb_text,
            confidence=0.5,
            fold_method=self.fold_method,
            provenance=Provenance(
                tool_versions={"synthetic_folder": "0.1.0"},
                data_version="synthetic_v1",
            ),
        )

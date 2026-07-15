"""Human stage-gate pauses for the agent (CURSOR_PROJECT_CONTEXT §8)."""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field


class StageGate(str, Enum):
    """Gates that require human sign-off before proceeding."""

    ORACLE_VALIDITY = "oracle_validity"
    CALIBRATION = "calibration"
    LOOP = "loop"
    QUANTUM = "quantum"
    SPEND = "spend"


class GateStatus(str, Enum):
    PENDING = "pending"
    APPROVED = "approved"
    REJECTED = "rejected"
    SKIPPED_SIMULATION = "skipped_simulation"


class GatePause(BaseModel):
    """A pause event recorded in the reasoning trace."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    gate: StageGate
    status: GateStatus
    reason: str
    evidence: dict[str, object] = Field(default_factory=dict)
    requires_human: bool = True


class StageGateManager:
    """Track gate decisions; block progress until approved or explicitly skipped."""

    def __init__(self, *, simulation_mode: bool = True) -> None:
        self.simulation_mode = simulation_mode
        self._decisions: dict[StageGate, GatePause] = {}

    @property
    def decisions(self) -> dict[StageGate, GatePause]:
        return dict(self._decisions)

    def request_pause(
        self,
        gate: StageGate,
        *,
        reason: str,
        evidence: dict[str, object] | None = None,
    ) -> GatePause:
        pause = GatePause(
            gate=gate,
            status=GateStatus.PENDING,
            reason=reason,
            evidence=evidence or {},
            requires_human=True,
        )
        self._decisions[gate] = pause
        return pause

    def approve(self, gate: StageGate, *, note: str = "human approved") -> GatePause:
        prev = self._decisions.get(gate)
        if prev is None:
            raise KeyError(f"no pending pause for gate {gate}")
        pause = GatePause(
            gate=gate,
            status=GateStatus.APPROVED,
            reason=f"{prev.reason} | {note}",
            evidence=prev.evidence,
            requires_human=False,
        )
        self._decisions[gate] = pause
        return pause

    def reject(self, gate: StageGate, *, note: str = "human rejected") -> GatePause:
        prev = self._decisions.get(gate)
        if prev is None:
            raise KeyError(f"no pending pause for gate {gate}")
        pause = GatePause(
            gate=gate,
            status=GateStatus.REJECTED,
            reason=f"{prev.reason} | {note}",
            evidence=prev.evidence,
            requires_human=False,
        )
        self._decisions[gate] = pause
        return pause

    def skip_for_simulation(self, gate: StageGate, *, reason: str) -> GatePause:
        """Allow simulation plumbing while recording that the real gate is not passed."""
        if not self.simulation_mode:
            raise RuntimeError(
                f"cannot skip_for_simulation({gate}) when simulation_mode=False"
            )
        pause = GatePause(
            gate=gate,
            status=GateStatus.SKIPPED_SIMULATION,
            reason=reason,
            evidence={"simulation_mode": True, "real_gate_passed": False},
            requires_human=False,
        )
        self._decisions[gate] = pause
        return pause

    def is_cleared(self, gate: StageGate) -> bool:
        d = self._decisions.get(gate)
        if d is None:
            return False
        return d.status in {GateStatus.APPROVED, GateStatus.SKIPPED_SIMULATION}

    def assert_cleared(self, gate: StageGate) -> None:
        if not self.is_cleared(gate):
            status = self._decisions[gate].status if gate in self._decisions else "missing"
            raise RuntimeError(
                f"stage gate {gate.value} not cleared (status={status}). "
                "Human sign-off required before proceeding."
            )

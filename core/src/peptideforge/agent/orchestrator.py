"""PeptideForge LLM agent — tool-calling orchestrator with stage gates (P10).

HARD RULES:
- Reports only real tool outputs
- Every number attributed to the tool that produced it
- Never invents values
- Pauses for human sign-off at stage gates
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field

from peptideforge.agent.gates import GateStatus, StageGate, StageGateManager
from peptideforge.agent.llm import (
    LLMClient,
    LLMMessage,
    MockLLMClient,
    default_synthetic_campaign_script,
)
from peptideforge.agent.tools import AttributionLedger, ToolRegistry, ToolResult
from peptideforge.agent.trace import (
    ReasoningTrace,
    validate_summary_against_ledger,
)
from peptideforge.loop.report import ToolAttributedNumber


class AgentReport(BaseModel):
    """Final agent report — numbers subset of the attribution ledger only."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    session_id: str
    summary: str
    numbers: tuple[ToolAttributedNumber, ...]
    propose_iterate: bool
    gate_statuses: dict[str, str]
    halted_on_gate: str | None = None
    notes: str | None = None


class PeptideForgeAgent:
    """LLM-driven campaign agent over PeptideForge tools."""

    def __init__(
        self,
        *,
        llm: LLMClient | None = None,
        tools: ToolRegistry | None = None,
        simulation_mode: bool = True,
        auto_skip_simulation_gates: bool = True,
        max_steps: int = 20,
    ) -> None:
        self.simulation_mode = simulation_mode
        self.auto_skip_simulation_gates = auto_skip_simulation_gates
        self.max_steps = max_steps
        self.tools = tools or ToolRegistry()
        self.llm = llm or MockLLMClient(default_synthetic_campaign_script())
        self.gates = StageGateManager(simulation_mode=simulation_mode)
        self.ledger = AttributionLedger()
        self.session_id = str(uuid4())
        self.trace = ReasoningTrace(session_id=self.session_id)
        self._human_approvals: dict[StageGate, bool] = {}

    def set_human_decision(self, gate: StageGate, *, approved: bool) -> None:
        """Inject human stage-gate sign-off (tests / UI)."""
        self._human_approvals[gate] = approved
        if approved:
            self.gates.approve(gate, note="human approved via set_human_decision")
        else:
            self.gates.reject(gate, note="human rejected via set_human_decision")

    def run(self, goal: str, *, trace_path: Path | None = None) -> AgentReport:
        """Execute tool-calling loop until done, rejected gate, or max_steps."""
        messages: list[LLMMessage] = [
            LLMMessage(
                role="system",
                content=(
                    "You are PeptideForge's campaign agent. "
                    "You may ONLY report numbers returned by tools. "
                    "Never invent affinities, Spearman, or costs. "
                    "Pause at stage gates for human sign-off."
                ),
            ),
            LLMMessage(role="user", content=goal),
        ]
        self.trace.thought(f"Received goal: {goal}")

        # Gate 1: oracle-validity — simulation may skip with explicit honesty
        self._handle_oracle_validity_gate()

        halted_on: str | None = None
        for step in range(self.max_steps):
            self.trace.thought(f"step={step}")
            response = self.llm.complete(messages)

            if response.tool_calls:
                for call in response.tool_calls:
                    # Spend gate before claiming spend-gate tool results as campaign go
                    if call.name == "check_spend_gate":
                        pause = self.gates.request_pause(
                            StageGate.SPEND,
                            reason="Spend gate validation requested",
                        )
                        self.trace.gate_pause(pause)
                        if not self._resolve_gate(StageGate.SPEND):
                            halted_on = StageGate.SPEND.value
                            return self._finalize(
                                propose_iterate=False,
                                halted_on_gate=halted_on,
                                notes="Halted: spend gate not approved",
                                trace_path=trace_path,
                            )

                    self.trace.tool_call(call.name, call.arguments)
                    result = self.tools.call(call.name, **call.arguments)
                    self.ledger.record(result)
                    self.trace.tool_result(result)
                    messages.append(
                        LLMMessage(
                            role="tool",
                            content=result.model_dump_json(),
                        )
                    )
                continue

            if response.done or response.content.strip().upper() == "STOP":
                break

            if response.content:
                self.trace.thought(response.content)
                messages.append(LLMMessage(role="assistant", content=response.content))

        return self._finalize(
            propose_iterate=self._should_iterate(),
            halted_on_gate=halted_on,
            notes=None,
            trace_path=trace_path,
        )

    def build_summary(self) -> tuple[str, tuple[ToolAttributedNumber, ...]]:
        """Compose a summary using ONLY ledger numbers (anti-hallucination)."""
        lines = ["## Campaign summary (tool-attributed only)", ""]
        for entry in self.ledger.entries:
            unit = f" {entry.unit}" if entry.unit else ""
            lines.append(
                f"- **{entry.name}** = {entry.value}{unit} _(tool: {entry.tool})_"
            )
        if not self.ledger.entries:
            lines.append("_No tool numbers recorded._")
        text = "\n".join(lines)
        bad = validate_summary_against_ledger(text, self.ledger)
        if bad:
            raise ValueError(f"summary contains unattributed numbers: {bad}")
        return text, self.ledger.entries

    def _finalize(
        self,
        *,
        propose_iterate: bool,
        halted_on_gate: str | None,
        notes: str | None,
        trace_path: Path | None,
    ) -> AgentReport:
        summary, numbers = self.build_summary()
        self.trace.summary(summary, numbers)
        if trace_path is not None:
            self.trace.save(trace_path)
        return AgentReport(
            session_id=self.session_id,
            summary=summary,
            numbers=numbers,
            propose_iterate=propose_iterate,
            gate_statuses={g.value: p.status.value for g, p in self.gates.decisions.items()},
            halted_on_gate=halted_on_gate,
            notes=notes,
        )

    def _handle_oracle_validity_gate(self) -> None:
        pause = self.gates.request_pause(
            StageGate.ORACLE_VALIDITY,
            reason=(
                "Oracle-validity gate: real MM-GBSA Spearman must pass before "
                "trusting binding campaigns (see ACCEPTANCE.md)."
            ),
            evidence={"simulation_mode": self.simulation_mode},
        )
        self.trace.gate_pause(pause)
        if self.simulation_mode:
            # Always record that the *real* gate is not passed when simulating.
            skipped = self.gates.skip_for_simulation(
                StageGate.ORACLE_VALIDITY,
                reason=(
                    "Simulation mode: proceeding with synthetic_* oracle only. "
                    "Real oracle-validity gate remains NOT passed."
                ),
            )
            self.trace.gate_pause(skipped)
            return
        if not self._resolve_gate(StageGate.ORACLE_VALIDITY):
            raise RuntimeError("oracle_validity gate blocked campaign start")

    def _resolve_gate(self, gate: StageGate) -> bool:
        """Apply pending human decision or simulation auto-approve when enabled."""
        if gate in self._human_approvals:
            return self._human_approvals[gate]

        decision = self.gates.decisions.get(gate)
        if decision and decision.status == GateStatus.PENDING:
            if (
                self.simulation_mode
                and self.auto_skip_simulation_gates
                and gate == StageGate.SPEND
            ):
                self.gates.approve(
                    gate,
                    note=(
                        "Simulation auto-approve for CI: still requires human for "
                        "production campaigns"
                    ),
                )
                self.trace.gate_pause(self.gates.decisions[gate])
                return True
            return False
        return self.gates.is_cleared(gate)

    def _should_iterate(self) -> bool:
        for e in self.ledger.entries:
            if e.name == "reached_target" and e.value < 0.5:
                return True
            if e.name == "spend_gate_passed" and e.value < 0.5:
                return True
        return False

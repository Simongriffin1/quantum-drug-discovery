"""Reasoning-trace log for agent runs (auditable, no invented numbers)."""

from __future__ import annotations

import json
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from peptideforge.agent.gates import GatePause
from peptideforge.agent.tools import AttributionLedger, ToolResult
from peptideforge.loop.report import ToolAttributedNumber


class TraceEvent(BaseModel):
    """One step in the agent reasoning / tool-call log."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    timestamp: datetime = Field(default_factory=lambda: datetime.now(UTC))
    kind: str  # thought | tool_call | tool_result | gate_pause | summary | error
    content: str
    tool_name: str | None = None
    call_id: str | None = None
    data: dict[str, Any] = Field(default_factory=dict)


class ReasoningTrace(BaseModel):
    """Full session trace — persistable for UI / audit."""

    model_config = ConfigDict(extra="forbid")

    session_id: str
    events: list[TraceEvent] = Field(default_factory=list)
    summary_numbers: tuple[ToolAttributedNumber, ...] = ()

    def add(self, event: TraceEvent) -> None:
        self.events.append(event)

    def thought(self, text: str) -> None:
        self.add(TraceEvent(kind="thought", content=text))

    def tool_call(self, name: str, args: dict[str, Any]) -> None:
        self.add(
            TraceEvent(
                kind="tool_call",
                content=f"call {name}",
                tool_name=name,
                data={"args": args},
            )
        )

    def tool_result(self, result: ToolResult) -> None:
        self.add(
            TraceEvent(
                kind="tool_result",
                content=result.message,
                tool_name=result.tool_name,
                call_id=result.call_id,
                data={"ok": result.ok, "data": result.data},
            )
        )

    def gate_pause(self, pause: GatePause) -> None:
        self.add(
            TraceEvent(
                kind="gate_pause",
                content=f"{pause.gate.value}: {pause.status.value} — {pause.reason}",
                data=pause.model_dump(),
            )
        )

    def summary(self, text: str, numbers: tuple[ToolAttributedNumber, ...]) -> None:
        self.summary_numbers = numbers
        self.add(
            TraceEvent(
                kind="summary",
                content=text,
                data={"n_numbers": len(numbers)},
            )
        )

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2), encoding="utf-8")


_FLOAT_RE = re.compile(
    r"(?<![A-Za-z_])([+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?)(?![A-Za-z_])"
)


def extract_floats(text: str) -> list[float]:
    """Extract numeric literals from agent prose for hallucination checks."""
    out: list[float] = []
    for match in _FLOAT_RE.finditer(text):
        raw = match.group(1)
        # Skip lone integers used as list indices / years if needed — keep all for strictness
        try:
            out.append(float(raw))
        except ValueError:
            continue
    return out


def validate_summary_against_ledger(
    text: str,
    ledger: AttributionLedger,
    *,
    ignore_values: frozenset[float] | None = None,
) -> list[float]:
    """Return list of floats in ``text`` that are NOT in the ledger (hallucinations).

    ``ignore_values`` may include benign constants (0, 1, iteration counts already
    registered). Empty list means the summary is clean.
    """
    ignore = ignore_values or frozenset({0.0, 1.0})
    bad: list[float] = []
    for value in extract_floats(text):
        if value in ignore:
            continue
        if ledger.has_exact(value, tol=1e-6):
            continue
        # Also allow values that match after rounding to 4–6 decimals as printed
        if any(
            abs(value - e.value) < 1e-4 or abs(value - round(e.value, 4)) < 1e-9
            for e in ledger.entries
        ):
            continue
        bad.append(value)
    return bad

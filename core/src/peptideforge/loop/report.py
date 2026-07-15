"""Iteration report stub — attribution of numbers to tools (expanded in P10)."""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict

from peptideforge.contracts.models import LoopState


class ToolAttributedNumber(BaseModel):
    """A numeric claim that must cite its producing tool (anti-hallucination)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    value: float
    tool: str
    unit: str | None = None


class IterationReport(BaseModel):
    """Human/agent-readable summary of one DBTL iteration."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    iteration: int
    oracle_calls: int
    total_cost: float
    best_oracle_value: float | None
    batch_size: int
    acquisition_method: str
    numbers: tuple[ToolAttributedNumber, ...] = ()
    notes: str | None = None


def build_iteration_report(
    state: LoopState,
    *,
    best_oracle_value: float | None,
    batch_size: int,
    acquisition_method: str,
) -> IterationReport:
    numbers: list[ToolAttributedNumber] = [
        ToolAttributedNumber(
            name="oracle_calls",
            value=float(state.oracle_calls),
            tool="ClosedLoopOrchestrator",
            unit="count",
        ),
        ToolAttributedNumber(
            name="total_cost",
            value=state.total_cost,
            tool="ClosedLoopOrchestrator",
            unit="cost_units",
        ),
    ]
    if best_oracle_value is not None:
        numbers.append(
            ToolAttributedNumber(
                name="best_oracle_value",
                value=best_oracle_value,
                tool="Oracle.evaluate",
                unit="oracle_unit",
            )
        )
    return IterationReport(
        iteration=state.iteration,
        oracle_calls=state.oracle_calls,
        total_cost=state.total_cost,
        best_oracle_value=best_oracle_value,
        batch_size=batch_size,
        acquisition_method=acquisition_method,
        numbers=tuple(numbers),
        notes=state.notes,
    )

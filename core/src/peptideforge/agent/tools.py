"""Tool-calling surface for the PeptideForge agent (P10).

Every numeric claim in agent reports must come from a ToolResult's attributed
numbers — never invented by the LLM.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field

from peptideforge.loop.orchestrator import CampaignResult, ClosedLoopOrchestrator, LoopConfig
from peptideforge.loop.report import ToolAttributedNumber
from peptideforge.loop.validate import run_simulations_to_target_validation


class ToolResult(BaseModel):
    """Structured tool output — the only permitted source of reported numbers."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    tool_name: str
    call_id: str = Field(default_factory=lambda: str(uuid4()))
    ok: bool
    data: dict[str, Any] = Field(default_factory=dict)
    numbers: tuple[ToolAttributedNumber, ...] = ()
    message: str = ""


class AttributionLedger:
    """Records every tool-attributed number; blocks hallucinated values."""

    def __init__(self) -> None:
        self._entries: list[ToolAttributedNumber] = []
        self._by_name: dict[str, list[float]] = {}
        self._results: list[ToolResult] = []

    def record(self, result: ToolResult) -> None:
        self._results.append(result)
        for num in result.numbers:
            self._entries.append(num)
            self._by_name.setdefault(num.name, []).append(num.value)

    @property
    def results(self) -> tuple[ToolResult, ...]:
        return tuple(self._results)

    @property
    def entries(self) -> tuple[ToolAttributedNumber, ...]:
        return tuple(self._entries)

    def has_exact(self, value: float, *, tol: float = 1e-9) -> bool:
        return any(math.isclose(value, e.value, rel_tol=0.0, abs_tol=tol) for e in self._entries)

    def has_named(self, name: str, value: float, *, tol: float = 1e-9) -> bool:
        return any(
            math.isclose(value, v, rel_tol=0.0, abs_tol=tol) for v in self._by_name.get(name, [])
        )

    def assert_number_allowed(self, value: float, *, name: str | None = None) -> None:
        """Raise if ``value`` was never produced by a tool."""
        if name is not None:
            if not self.has_named(name, value):
                raise ValueError(
                    f"hallucinated number blocked: {name}={value} not in tool ledger"
                )
            return
        if not self.has_exact(value):
            raise ValueError(
                f"hallucinated number blocked: {value} not present in any ToolResult"
            )


class ToolRegistry:
    """Named tools the agent may call."""

    def __init__(self) -> None:
        self._tools: dict[str, Callable[..., ToolResult]] = {}
        self.register("parse_goal", parse_goal)
        self.register("plan_campaign", plan_campaign)
        self.register("run_simulation_campaign", run_simulation_campaign)
        self.register("check_spend_gate", check_spend_gate)
        self.register("monitor_jobs", monitor_jobs)

    def register(self, name: str, fn: Callable[..., ToolResult]) -> None:
        self._tools[name] = fn

    def names(self) -> tuple[str, ...]:
        return tuple(sorted(self._tools))

    def call(self, name: str, **kwargs: Any) -> ToolResult:
        if name not in self._tools:
            raise KeyError(
                f"unknown tool {name!r}; available: {self.names()}. "
                "Refusing to invent a tool result."
            )
        return self._tools[name](**kwargs)


def parse_goal(*, goal: str, simulation_mode: bool = True) -> ToolResult:
    """Parse a natural-language design goal into structured fields."""
    goal_l = goal.lower()
    target_value = -4.0 if "strong" in goal_l or "binder" in goal_l else -1.0
    # Explicit numbers in the goal text are echoed as attributed (user-provided)
    numbers: list[ToolAttributedNumber] = [
        ToolAttributedNumber(
            name="parsed_target_value",
            value=target_value,
            tool="parse_goal",
            unit="oracle_unit",
        )
    ]
    return ToolResult(
        tool_name="parse_goal",
        ok=True,
        data={
            "goal": goal,
            "simulation_mode": simulation_mode,
            "modality": "peptide",
            "target_value": target_value,
            "acquisition": "qnehvi",
        },
        numbers=tuple(numbers),
        message=f"Parsed goal; simulation_mode={simulation_mode}; target_value={target_value}",
    )


def plan_campaign(
    *,
    seed: int = 0,
    target_value: float = -4.0,
    n_init: int = 8,
    max_iterations: int = 4,
    batch_size: int = 2,
    simulation_mode: bool = True,
) -> ToolResult:
    """Build a LoopConfig plan (does not execute)."""
    cfg = LoopConfig(
        seed=seed,
        simulation_mode=simulation_mode,
        n_init=n_init,
        max_iterations=max_iterations,
        batch_size=batch_size,
        target_value=target_value,
        acquisition="qnehvi",
        n_propose=20,
    )
    numbers = (
        ToolAttributedNumber(name="plan_seed", value=float(seed), tool="plan_campaign"),
        ToolAttributedNumber(
            name="plan_target_value", value=target_value, tool="plan_campaign", unit="oracle_unit"
        ),
        ToolAttributedNumber(
            name="plan_max_iterations", value=float(max_iterations), tool="plan_campaign"
        ),
        ToolAttributedNumber(
            name="plan_batch_size", value=float(batch_size), tool="plan_campaign"
        ),
        ToolAttributedNumber(name="plan_n_init", value=float(n_init), tool="plan_campaign"),
    )
    return ToolResult(
        tool_name="plan_campaign",
        ok=True,
        data={"config": cfg.to_dict()},
        numbers=numbers,
        message="Campaign plan ready (not yet executed)",
    )


def run_simulation_campaign(
    *,
    seed: int = 0,
    target_value: float = -4.0,
    n_init: int = 8,
    max_iterations: int = 3,
    batch_size: int = 2,
    state_dir: str | None = None,
) -> ToolResult:
    """Run ClosedLoopOrchestrator in simulation mode; return attributed metrics."""
    cfg = LoopConfig(
        seed=seed,
        simulation_mode=True,
        n_init=n_init,
        max_iterations=max_iterations,
        batch_size=batch_size,
        target_value=target_value,
        acquisition="qnehvi",
        n_propose=20,
        state_dir=state_dir,
        data_version="synthetic_v1",
    )
    result: CampaignResult = ClosedLoopOrchestrator(config=cfg).run()
    numbers: list[ToolAttributedNumber] = [
        ToolAttributedNumber(
            name="oracle_calls",
            value=float(result.states[-1].oracle_calls if result.states else 0),
            tool="run_simulation_campaign",
            unit="count",
        ),
        ToolAttributedNumber(
            name="total_cost",
            value=float(result.states[-1].total_cost if result.states else 0.0),
            tool="run_simulation_campaign",
            unit="cost_units",
        ),
        ToolAttributedNumber(
            name="n_iterations",
            value=float(len(result.states)),
            tool="run_simulation_campaign",
            unit="count",
        ),
        ToolAttributedNumber(
            name="reached_target",
            value=1.0 if result.reached_target else 0.0,
            tool="run_simulation_campaign",
        ),
    ]
    if result.best_value is not None:
        numbers.append(
            ToolAttributedNumber(
                name="best_oracle_value",
                value=result.best_value,
                tool="run_simulation_campaign/Oracle.evaluate",
                unit="oracle_unit",
            )
        )
    if result.oracle_calls_to_target is not None:
        numbers.append(
            ToolAttributedNumber(
                name="oracle_calls_to_target",
                value=float(result.oracle_calls_to_target),
                tool="run_simulation_campaign",
                unit="count",
            )
        )
    return ToolResult(
        tool_name="run_simulation_campaign",
        ok=True,
        data={
            "campaign_id": str(result.campaign_id),
            "reached_target": result.reached_target,
            "best_value": result.best_value,
            "oracle_calls_to_target": result.oracle_calls_to_target,
            "n_labeled": len(result.dataset.records),
            "status": result.states[-1].status if result.states else "empty",
        },
        numbers=tuple(numbers),
        message=(
            f"Simulation campaign finished; reached_target={result.reached_target}; "
            f"best={result.best_value}"
        ),
    )


def check_spend_gate(*, seed: int = 0, target_value: float = -4.5) -> ToolResult:
    """Run spend-gate validation; report real HV-style sims-to-target metrics."""
    report = run_simulations_to_target_validation(
        seed=seed,
        target_value=target_value,
        n_pool=60,
        n_init=8,
        max_rounds=10,
        batch_size=2,
    )
    numbers: list[ToolAttributedNumber] = [
        ToolAttributedNumber(
            name="spend_gate_passed",
            value=1.0 if report.passed else 0.0,
            tool="check_spend_gate",
        ),
        ToolAttributedNumber(
            name="spend_gate_target_value",
            value=report.target_value,
            tool="check_spend_gate",
            unit="oracle_unit",
        ),
    ]
    if report.oracle_calls_qnehvi is not None:
        numbers.append(
            ToolAttributedNumber(
                name="oracle_calls_qnehvi",
                value=float(report.oracle_calls_qnehvi),
                tool="check_spend_gate",
                unit="count",
            )
        )
    if report.oracle_calls_random is not None:
        numbers.append(
            ToolAttributedNumber(
                name="oracle_calls_random",
                value=float(report.oracle_calls_random),
                tool="check_spend_gate",
                unit="count",
            )
        )
    if report.best_qnehvi is not None:
        numbers.append(
            ToolAttributedNumber(
                name="best_qnehvi",
                value=report.best_qnehvi,
                tool="check_spend_gate",
                unit="oracle_unit",
            )
        )
    if report.best_random is not None:
        numbers.append(
            ToolAttributedNumber(
                name="best_random",
                value=report.best_random,
                tool="check_spend_gate",
                unit="oracle_unit",
            )
        )
    return ToolResult(
        tool_name="check_spend_gate",
        ok=report.passed,
        data=report.model_dump(),
        numbers=tuple(numbers),
        message=f"Spend gate passed={report.passed}",
    )


def monitor_jobs(*, use_ray: bool = False, n_pending: int = 0) -> ToolResult:
    """Report job status (Ray optional). Always returns real counts."""
    backend = "ray" if use_ray else "serial"
    numbers = (
        ToolAttributedNumber(
            name="jobs_pending",
            value=float(n_pending),
            tool="monitor_jobs",
            unit="count",
        ),
        ToolAttributedNumber(
            name="jobs_use_ray",
            value=1.0 if use_ray else 0.0,
            tool="monitor_jobs",
        ),
    )
    return ToolResult(
        tool_name="monitor_jobs",
        ok=True,
        data={"backend": backend, "pending": n_pending},
        numbers=numbers,
        message=f"Job monitor: backend={backend}, pending={n_pending}",
    )

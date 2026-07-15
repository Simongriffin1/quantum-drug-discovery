"""Tests for P10 LLM agent — attribution, stage gates, synthetic campaign."""

from __future__ import annotations

from pathlib import Path

import pytest

from peptideforge.agent import (
    AttributionLedger,
    GateStatus,
    LLMResponse,
    LLMUnavailableError,
    MockLLMClient,
    PeptideForgeAgent,
    StageGate,
    ToolCallRequest,
    ToolRegistry,
    ToolResult,
    extract_floats,
    require_openai_client,
    validate_summary_against_ledger,
)
from peptideforge.loop.report import ToolAttributedNumber


def test_ledger_blocks_hallucinated_number() -> None:
    ledger = AttributionLedger()
    ledger.record(
        ToolResult(
            tool_name="t",
            ok=True,
            numbers=(
                ToolAttributedNumber(name="oracle_calls", value=10.0, tool="t"),
            ),
        )
    )
    ledger.assert_number_allowed(10.0, name="oracle_calls")
    with pytest.raises(ValueError, match="hallucinated"):
        ledger.assert_number_allowed(0.99, name="spearman")


def test_summary_validation_flags_invented_values() -> None:
    ledger = AttributionLedger()
    ledger.record(
        ToolResult(
            tool_name="run_simulation_campaign",
            ok=True,
            numbers=(
                ToolAttributedNumber(
                    name="best_oracle_value", value=-4.2, tool="run_simulation_campaign"
                ),
            ),
        )
    )
    clean = "best_oracle_value = -4.2 (tool: run_simulation_campaign)"
    assert validate_summary_against_ledger(clean, ledger) == []
    dirty = "Spearman rho = 0.87 on PDBbind"
    bad = validate_summary_against_ledger(dirty, ledger)
    assert 0.87 in bad


def test_oracle_validity_gate_skipped_in_simulation() -> None:
    agent = PeptideForgeAgent(
        llm=MockLLMClient([LLMResponse(content="STOP", done=True)]),
        simulation_mode=True,
        auto_skip_simulation_gates=True,
    )
    report = agent.run("dry run")
    assert report.gate_statuses[StageGate.ORACLE_VALIDITY.value] == GateStatus.SKIPPED_SIMULATION.value
    assert any(e.kind == "gate_pause" for e in agent.trace.events)


def test_oracle_validity_gate_blocks_without_approval() -> None:
    agent = PeptideForgeAgent(
        llm=MockLLMClient([LLMResponse(content="STOP", done=True)]),
        simulation_mode=False,
        auto_skip_simulation_gates=False,
    )
    with pytest.raises(RuntimeError, match="oracle_validity"):
        agent.run("live campaign")


def test_spend_gate_pause_fires_and_can_reject(tmp_path: Path) -> None:
    script = [
        LLMResponse(
            tool_calls=(
                ToolCallRequest(
                    name="check_spend_gate",
                    arguments={"seed": 0, "target_value": -4.5},
                ),
            )
        ),
        LLMResponse(content="STOP", done=True),
    ]
    agent = PeptideForgeAgent(
        llm=MockLLMClient(script),
        simulation_mode=True,
        auto_skip_simulation_gates=False,  # force human path
    )
    # Reject spend when pause is pending — agent halts before tool runs
    # First request_pause happens inside run; we need to reject via callback.
    # Without auto-skip, _resolve_gate returns False → halt
    report = agent.run("check spend", trace_path=tmp_path / "trace.json")
    assert report.halted_on_gate == StageGate.SPEND.value
    assert (tmp_path / "trace.json").is_file()


@pytest.mark.eval
def test_agent_synthetic_campaign_end_to_end(tmp_path: Path) -> None:
    agent = PeptideForgeAgent(
        simulation_mode=True,
        auto_skip_simulation_gates=True,
    )
    report = agent.run(
        "Find a strong binder in simulation mode",
        trace_path=tmp_path / "agent_trace.json",
    )
    assert report.numbers, "expected tool-attributed numbers"
    assert report.summary
    # Every reported number is in the ledger
    for num in report.numbers:
        agent.ledger.assert_number_allowed(num.value, name=num.name)
    # Summary must not contain floats outside the ledger
    bad = validate_summary_against_ledger(report.summary, agent.ledger)
    assert bad == [], f"hallucinated floats in summary: {bad}"
    assert StageGate.ORACLE_VALIDITY.value in report.gate_statuses
    assert (tmp_path / "agent_trace.json").is_file()
    # Tool results present in trace
    kinds = {e.kind for e in agent.trace.events}
    assert "tool_result" in kinds
    assert "gate_pause" in kinds


def test_agent_rejects_summary_with_planted_hallucination() -> None:
    agent = PeptideForgeAgent(
        llm=MockLLMClient([LLMResponse(content="STOP", done=True)]),
        simulation_mode=True,
    )
    agent.run("noop")
    dirty = "## Campaign summary\n\n- fake_spearman = 0.99 _(tool: nowhere)_"
    bad = validate_summary_against_ledger(dirty, agent.ledger)
    assert 0.99 in bad


def test_openai_backend_fails_loud_without_key(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    try:
        import openai  # noqa: F401
    except ImportError:
        with pytest.raises(LLMUnavailableError):
            require_openai_client()
    else:
        with pytest.raises(LLMUnavailableError, match="OPENAI_API_KEY"):
            require_openai_client()


def test_unknown_tool_fails_loud() -> None:
    reg = ToolRegistry()
    with pytest.raises(KeyError, match="unknown tool"):
        reg.call("invent_affinity")


def test_extract_floats() -> None:
    assert 0.4 in extract_floats("Spearman ρ ≥ 0.40")
    assert -4.5 in extract_floats("target=-4.5 kcal")

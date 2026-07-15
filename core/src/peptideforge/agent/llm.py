"""LLM backends for the agent — mock for CI; real APIs fail loud if missing."""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable

from pydantic import BaseModel, ConfigDict, Field


class LLMMessage(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    role: str  # system | user | assistant | tool
    content: str


class ToolCallRequest(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    arguments: dict[str, Any] = Field(default_factory=dict)


class LLMResponse(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    content: str = ""
    tool_calls: tuple[ToolCallRequest, ...] = ()
    done: bool = False


class LLMUnavailableError(ImportError):
    """Raised when a real LLM backend is requested but credentials/deps are missing."""


@runtime_checkable
class LLMClient(Protocol):
    def complete(self, messages: list[LLMMessage]) -> LLMResponse:
        """Return next assistant turn (tool calls and/or final content)."""
        ...


class MockLLMClient:
    """Deterministic scripted LLM for CI — never invents tool results."""

    def __init__(self, script: list[LLMResponse] | None = None) -> None:
        self._script = list(script or [])
        self._idx = 0

    def complete(self, messages: list[LLMMessage]) -> LLMResponse:
        del messages
        if self._idx >= len(self._script):
            return LLMResponse(
                content="No further actions scripted.",
                done=True,
            )
        resp = self._script[self._idx]
        self._idx += 1
        return resp


def default_synthetic_campaign_script(
    *,
    seed: int = 0,
    target_value: float = -4.0,
) -> list[LLMResponse]:
    """Script: parse → plan → run campaign → monitor → spend gate → stop."""
    return [
        LLMResponse(
            tool_calls=(
                ToolCallRequest(
                    name="parse_goal",
                    arguments={
                        "goal": "Find a strong binder in simulation mode",
                        "simulation_mode": True,
                    },
                ),
            )
        ),
        LLMResponse(
            tool_calls=(
                ToolCallRequest(
                    name="plan_campaign",
                    arguments={
                        "seed": seed,
                        "target_value": target_value,
                        "n_init": 8,
                        "max_iterations": 3,
                        "batch_size": 2,
                        "simulation_mode": True,
                    },
                ),
            )
        ),
        LLMResponse(
            tool_calls=(
                ToolCallRequest(
                    name="run_simulation_campaign",
                    arguments={
                        "seed": seed,
                        "target_value": target_value,
                        "n_init": 8,
                        "max_iterations": 3,
                        "batch_size": 2,
                    },
                ),
            )
        ),
        LLMResponse(
            tool_calls=(
                ToolCallRequest(
                    name="monitor_jobs",
                    arguments={"use_ray": False, "n_pending": 0},
                ),
            )
        ),
        LLMResponse(
            tool_calls=(
                ToolCallRequest(
                    name="check_spend_gate",
                    arguments={"seed": seed, "target_value": -4.5},
                ),
            )
        ),
        LLMResponse(
            content="STOP",
            done=True,
        ),
    ]


def require_openai_client(*, model: str = "gpt-4o-mini") -> LLMClient:
    """Fail loud if OpenAI SDK / API key are unavailable — never fake completions."""
    try:
        from openai import OpenAI
    except ImportError as exc:
        raise LLMUnavailableError(
            "OpenAI SDK required for live LLM agent but is not installed. "
            "Install with: pip install openai. Use MockLLMClient for CI."
        ) from exc
    import os

    if not os.environ.get("OPENAI_API_KEY"):
        raise LLMUnavailableError(
            "OPENAI_API_KEY is not set. Refusing to invent LLM responses. "
            "Use MockLLMClient for CI / offline runs."
        )

    class _OpenAIClient:
        def __init__(self) -> None:
            self._client = OpenAI()
            self._model = model

        def complete(self, messages: list[LLMMessage]) -> LLMResponse:
            # Minimal live path — tool calling full schema left for product wiring.
            # Still never invents numeric claims; tools remain the source of truth.
            resp = self._client.chat.completions.create(
                model=self._model,
                messages=[{"role": m.role, "content": m.content} for m in messages],
            )
            text = resp.choices[0].message.content or ""
            return LLMResponse(content=text, done=True)

    return _OpenAIClient()

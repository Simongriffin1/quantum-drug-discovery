"""LLM agent orchestrator — P10.

Tool-calling agent with attribution ledger, stage-gate pauses, and reasoning
trace. Reports only real tool outputs — never invents numbers.
"""

from peptideforge.agent.gates import (
    GatePause,
    GateStatus,
    StageGate,
    StageGateManager,
)
from peptideforge.agent.llm import (
    LLMClient,
    LLMMessage,
    LLMResponse,
    LLMUnavailableError,
    MockLLMClient,
    ToolCallRequest,
    default_synthetic_campaign_script,
    require_openai_client,
)
from peptideforge.agent.orchestrator import AgentReport, PeptideForgeAgent
from peptideforge.agent.tools import (
    AttributionLedger,
    ToolRegistry,
    ToolResult,
)
from peptideforge.agent.trace import (
    ReasoningTrace,
    extract_floats,
    validate_summary_against_ledger,
)

__all__ = [
    "AgentReport",
    "AttributionLedger",
    "GatePause",
    "GateStatus",
    "LLMClient",
    "LLMMessage",
    "LLMResponse",
    "LLMUnavailableError",
    "MockLLMClient",
    "PeptideForgeAgent",
    "ReasoningTrace",
    "StageGate",
    "StageGateManager",
    "ToolCallRequest",
    "ToolRegistry",
    "ToolResult",
    "default_synthetic_campaign_script",
    "extract_floats",
    "require_openai_client",
    "validate_summary_against_ledger",
]

"""Optional ESM-2 dependency guards — fail loud, never fabricate sequences."""

from __future__ import annotations

from typing import Any


class ESM2UnavailableError(ImportError):
    """Raised when fair-esm / ESM-2 cannot be imported or loaded."""


def require_esm2() -> tuple[Any, Any]:
    """Import ``esm`` and ``torch`` or raise loudly."""
    try:
        import esm
        import torch
    except ImportError as exc:
        raise ESM2UnavailableError(
            "ESM-2 (fair-esm) is required for masked LM sampling but is not installed. "
            "Install with: pip install fair-esm torch. "
            "Use MutationGenerator for dependency-free generation in CI."
        ) from exc
    return esm, torch

"""Optional Boltz-2 dependency guards."""

from __future__ import annotations

import shutil


class Boltz2UnavailableError(ImportError):
    """Raised when Boltz-2 CLI / package is not available."""


def require_boltz_cli() -> str:
    """Return path to ``boltz`` executable or raise loudly."""
    path = shutil.which("boltz")
    if path is None:
        raise Boltz2UnavailableError(
            "Boltz-2 CLI is required for structure prediction but `boltz` was not found "
            "on PATH. Install with: pip install boltz. "
            "Use FixtureStructurePredictor for CI without large model weights."
        )
    return path

"""Provenance helpers — git SHA, data version, tool versions on every artifact."""

from __future__ import annotations

import subprocess
from functools import lru_cache

from peptideforge import __version__ as core_version
from peptideforge.contracts.models import Provenance


@lru_cache(maxsize=1)
def git_sha() -> str | None:
    """Return current HEAD SHA or None if not in a git repo / git missing."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:
        return None
    sha = out.stdout.strip()
    return sha or None


def campaign_provenance(*, data_version: str) -> Provenance:
    """Build provenance for a campaign artifact."""
    return Provenance(
        git_sha=git_sha(),
        data_version=data_version,
        tool_versions={
            "peptideforge": core_version,
            "peptideforge_loop": "0.1.0",
            "peptideforge_agent": "0.1.0",
            "platform": "p11",
        },
    )

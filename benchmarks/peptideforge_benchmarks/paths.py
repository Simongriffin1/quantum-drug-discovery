"""Paths and fixture discovery (no network)."""

from __future__ import annotations

from pathlib import Path


def fixtures_dir() -> Path:
    """Return ``benchmarks/fixtures`` next to this package."""
    return Path(__file__).resolve().parent.parent / "fixtures"


def require_fixture(name: str) -> Path:
    """Resolve a fixture file or raise FileNotFoundError (fail loud)."""
    path = fixtures_dir() / name
    if not path.is_file():
        raise FileNotFoundError(
            f"benchmark fixture missing: {path}. "
            "CI must ship fixtures under benchmarks/fixtures/ (no live download)."
        )
    return path

"""Optional Ray parallel executor — serial fallback for CI (fail loud if Ray required)."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")


class RayUnavailableError(ImportError):
    """Raised when Ray is required but not installed."""


def require_ray() -> object:
    try:
        import ray
    except ImportError as exc:
        raise RayUnavailableError(
            "Ray is required for parallel fold/oracle jobs but is not installed. "
            "Install with: pip install ray. Use require_ray=False for serial CI."
        ) from exc
    return ray


def map_parallel(
    fn: Callable[[T], R],
    items: Sequence[T],
    *,
    use_ray: bool = False,
) -> list[R]:
    """Map ``fn`` over items; Ray if requested and available, else serial."""
    if not items:
        return []
    if not use_ray:
        return [fn(item) for item in items]
    ray = require_ray()
    if not ray.is_initialized():  # type: ignore[attr-defined]
        ray.init(ignore_reinit_error=True, logging_level="ERROR")  # type: ignore[attr-defined]
    remote_fn = ray.remote(fn)  # type: ignore[attr-defined]
    futures = [remote_fn.remote(item) for item in items]
    return list(ray.get(futures))  # type: ignore[attr-defined]

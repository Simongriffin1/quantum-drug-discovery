"""Unit tests for cost caps and tier resolution."""

from __future__ import annotations

import pytest

from peptideforge.contracts.models import OracleTier
from peptideforge.oracles.costs import (
    CostCapExceededError,
    enforce_cost_cap,
    resolve_tier,
)


def test_default_tier_is_cheapest() -> None:
    assert resolve_tier(None) == OracleTier.DOCKING


def test_cost_cap_blocks_expensive_tier() -> None:
    with pytest.raises(CostCapExceededError):
        enforce_cost_cap(OracleTier.MM_GBSA, cost_cap=0.1)


def test_cost_cap_allows_docking() -> None:
    cost = enforce_cost_cap(OracleTier.DOCKING, cost_cap=0.1)
    assert cost == pytest.approx(0.05)


def test_unknown_tier_fails_loud() -> None:
    with pytest.raises(ValueError):
        resolve_tier(OracleTier.FEP, available=(OracleTier.DOCKING, OracleTier.MM_GBSA))

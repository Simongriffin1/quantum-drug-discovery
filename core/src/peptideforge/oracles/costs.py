"""Cost ladder and cost-cap enforcement for physics oracles."""

from __future__ import annotations

from peptideforge.contracts.models import OracleTier


class CostCapExceededError(RuntimeError):
    """Raised when a requested tier's estimated cost exceeds the job cost cap."""


# Relative compute-cost units (GPU-minutes-ish). Default selection uses cheapest.
TIER_COST_ESTIMATE: dict[OracleTier, float] = {
    OracleTier.DOCKING: 0.05,
    OracleTier.MM_GBSA: 1.0,
    OracleTier.MD: 10.0,
    OracleTier.FEP: 100.0,
    OracleTier.QCHEM_CLASSICAL: 50.0,
    OracleTier.QCHEM_VQE: 80.0,
    OracleTier.SYNTHETIC: 0.0,
}

# Prefer cheapest → most expensive when tier is omitted.
DEFAULT_TIER_ORDER: tuple[OracleTier, ...] = (
    OracleTier.DOCKING,
    OracleTier.MM_GBSA,
    OracleTier.MD,
)


def estimate_cost(tier: OracleTier) -> float:
    if tier not in TIER_COST_ESTIMATE:
        raise KeyError(f"no cost estimate for tier={tier}")
    return TIER_COST_ESTIMATE[tier]


def enforce_cost_cap(tier: OracleTier, cost_cap: float | None) -> float:
    """Return estimated cost, or raise if it exceeds ``cost_cap``."""
    cost = estimate_cost(tier)
    if cost_cap is not None and cost > cost_cap:
        raise CostCapExceededError(
            f"tier={tier.value} cost_estimate={cost} exceeds cost_cap={cost_cap}"
        )
    return cost


def resolve_tier(
    requested: OracleTier | None,
    *,
    available: tuple[OracleTier, ...] = DEFAULT_TIER_ORDER,
) -> OracleTier:
    """Default to the cheapest available tier. Unknown tiers fail loud."""
    if requested is None:
        return available[0]
    if requested not in available and requested not in TIER_COST_ESTIMATE:
        raise ValueError(f"unsupported oracle tier: {requested}")
    if requested not in available:
        raise ValueError(
            f"tier={requested.value} not enabled on this oracle; "
            f"available={[t.value for t in available]}"
        )
    return requested

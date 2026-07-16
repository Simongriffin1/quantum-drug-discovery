"""Claim-language guards reflecting validated correlation ceiling (Step 4C).

A validated ρ≈0.38 supports RANKING/ENRICHMENT, not confident per-peptide
affinity point calls. Live predicted-fold campaigns use the (possibly degraded)
predicted ρ when available.
"""

from __future__ import annotations

from dataclasses import dataclass

from peptideforge.authorization import AuthorizationRecord, InputType, TaskType


class PointAffinityClaimError(ValueError):
    """Raised when a point-affinity claim is emitted without authorization."""


@dataclass(frozen=True)
class ClaimCeiling:
    """How strongly we may speak about within-target results."""

    validated_rho: float
    framing: str  # ranking_enrichment | blocked
    max_claim: str
    allow_point_affinity: bool = False


def claim_ceiling_from_authorization(rec: AuthorizationRecord) -> ClaimCeiling:
    if not rec.authorized:
        return ClaimCeiling(
            validated_rho=0.0,
            framing="blocked",
            max_claim="No binding affinity claims authorized for this campaign class.",
            allow_point_affinity=False,
        )
    rho = float(rec.validated_rho or 0.0)
    # ρ≈0.38 → ~15% variance — enrichment language only
    return ClaimCeiling(
        validated_rho=rho,
        framing="ranking_enrichment",
        max_claim=(
            f"Prioritized candidates with expected enrichment over random "
            f"(validated within-target Spearman ρ≈{rho:.2f}; not a point ΔΔG call)."
        ),
        allow_point_affinity=False,
    )


def assert_no_point_affinity_claim(
    *,
    text: str,
    task_type: TaskType,
    input_type: InputType,
    authorization: AuthorizationRecord | None,
) -> None:
    """Reject copy that asserts unsupported point affinities for within-target results."""
    if task_type != TaskType.WITHIN_TARGET:
        return
    if input_type == InputType.SIMULATION:
        return
    lowered = text.lower()
    banned = (
        "predicted affinity is",
        "kd =",
        "δg =",
        "ddg =",
        "binding free energy of",
        "will bind with",
    )
    if any(b in lowered for b in banned):
        ceiling = (
            claim_ceiling_from_authorization(authorization)
            if authorization is not None
            else None
        )
        raise PointAffinityClaimError(
            "Point-affinity claim not authorized for within-target results. "
            f"Use ranking/enrichment framing. Ceiling={ceiling}"
        )


def enrich_interval_width(base_uncertainty: float, validated_rho: float) -> float:
    """Inflate displayed uncertainty when validated ρ is modest.

    Maps ρ∈[0,1] → scale ∈[~1, ~4]; ρ=0.38 → ~2.6× wider intervals.
    """
    rho = max(0.0, min(1.0, validated_rho))
    scale = 1.0 + 3.0 * (1.0 - rho)
    return base_uncertainty * scale

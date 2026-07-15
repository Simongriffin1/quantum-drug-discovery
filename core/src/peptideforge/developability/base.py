"""Shared helpers for developability predictors."""

from __future__ import annotations

from peptideforge.contracts.models import (
    DevelopabilityProperty,
    DevelopabilityScores,
    PeptideCandidate,
    PropertyScore,
)


def single_score_result(
    candidate: PeptideCandidate,
    *,
    property_name: DevelopabilityProperty | str,
    value: float,
    uncertainty: float,
    higher_is_better: bool,
    method: str,
) -> DevelopabilityScores:
    """Build a DevelopabilityScores payload with one PropertyScore."""
    return DevelopabilityScores(
        candidate_id=candidate.candidate_id,
        scores=(
            PropertyScore(
                property_name=property_name,
                value=value,
                uncertainty=uncertainty,
                higher_is_better=higher_is_better,
                method=method,
            ),
        ),
    )

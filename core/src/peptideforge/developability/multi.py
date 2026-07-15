"""Multi-objective developability evaluator — per-property vector, no scalar collapse."""

from __future__ import annotations

from peptideforge.contracts.models import (
    DevelopabilityProperty,
    DevelopabilityScores,
    PeptideCandidate,
    PropertyScore,
    Provenance,
)
from peptideforge.contracts.protocols import DevelopabilityPredictor
from peptideforge.developability.aggregation import AggregationPredictor
from peptideforge.developability.half_life import HalfLifePredictor
from peptideforge.developability.immunogenicity import ImmunogenicityPredictor
from peptideforge.developability.solubility import SolubilityPredictor
from peptideforge.developability.synthesizability import SynthesizabilityPredictor

DEFAULT_PREDICTORS: tuple[DevelopabilityPredictor, ...] = (
    AggregationPredictor(),
    SolubilityPredictor(),
    ImmunogenicityPredictor(),
    SynthesizabilityPredictor(),
    HalfLifePredictor(),
)

EXPECTED_PROPERTIES = frozenset(p.value for p in DevelopabilityProperty)


class PeptideDevelopabilityEvaluator:
    """Aggregate all developability axes into one DevelopabilityScores vector.

    Biological rationale: developability failures (aggregation, poor solubility,
    immunogenicity, synthesis, rapid clearance) are independent axes — collapsing
    to a single score hides trade-offs the multi-objective loop must see.
    """

    def __init__(
        self,
        predictors: tuple[DevelopabilityPredictor, ...] | None = None,
    ) -> None:
        self.predictors = predictors if predictors is not None else DEFAULT_PREDICTORS
        if not self.predictors:
            raise ValueError("PeptideDevelopabilityEvaluator requires ≥1 predictor")

    def evaluate(self, candidate: PeptideCandidate) -> DevelopabilityScores:
        """Return per-property (value, uncertainty) for every registered predictor."""
        scores: list[PropertyScore] = []
        for predictor in self.predictors:
            result = predictor.predict(candidate)
            if len(result.scores) != 1:
                raise ValueError(
                    f"predictor {predictor.property_name} returned "
                    f"{len(result.scores)} scores; expected 1"
                )
            scores.append(result.scores[0])
        return DevelopabilityScores(
            candidate_id=candidate.candidate_id,
            scores=tuple(scores),
            provenance=Provenance(
                tool_versions={
                    "peptideforge_developability": "0.1.0",
                    "predictors": ",".join(p.property_name for p in self.predictors),
                },
            ),
        )

    @property
    def property_names(self) -> tuple[str, ...]:
        return tuple(p.property_name for p in self.predictors)

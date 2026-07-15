"""Aggregation propensity (AGGRESCAN-style hydrophobic patch score).

Method: sliding-window average of Kyte-Doolittle hydropathy; the maximum window
score approximates the strongest hydrophobic patch (Conchillo-Solé et al., 2007,
Bioinformatics — AGGRESCAN concept). Higher score → higher aggregation risk.
No training data; purely algorithmic.
"""

from __future__ import annotations

from peptideforge.contracts.models import DevelopabilityProperty, DevelopabilityScores, PeptideCandidate
from peptideforge.developability.amino_acids import KYTE_DOOLITTLE, validate_sequence
from peptideforge.developability.base import single_score_result

DEFAULT_WINDOW = 5


class AggregationPredictor:
    """AGGRESCAN-style hydrophobic patch aggregation score."""

    property_name = DevelopabilityProperty.AGGREGATION.value

    def __init__(self, window: int = DEFAULT_WINDOW) -> None:
        if window < 2:
            raise ValueError("aggregation window must be ≥ 2")
        self.window = window

    def predict(self, candidate: PeptideCandidate) -> DevelopabilityScores:
        seq = validate_sequence(candidate.sequence)
        if len(seq) < self.window:
            # Short peptides: use whole-sequence average hydropathy as proxy
            patch_max = sum(KYTE_DOOLITTLE[aa] for aa in seq) / len(seq)
            uncertainty = 0.15
        else:
            windows = [
                sum(KYTE_DOOLITTLE[aa] for aa in seq[i : i + self.window]) / self.window
                for i in range(len(seq) - self.window + 1)
            ]
            patch_max = max(windows)
            uncertainty = 0.08
        # Normalize roughly to [0, 1] for interpretability (typical patches ~0–4)
        normalized = max(0.0, min(1.0, (patch_max + 1.0) / 5.0))
        return single_score_result(
            candidate,
            property_name=DevelopabilityProperty.AGGREGATION,
            value=normalized,
            uncertainty=uncertainty,
            higher_is_better=False,
            method=f"aggrescan_style_w{self.window}_kytedoolittle",
        )

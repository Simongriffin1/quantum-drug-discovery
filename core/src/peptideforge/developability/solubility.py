"""Solubility heuristics: GRAVY + net charge at pH 7.4.

Method: GRAVY (Kyte-Doolittle grand average of hydropathy) correlates inversely with
aqueous solubility; net positive charge at pH 7.4 improves colloidal stability for
peptides. Composite solubility score = −GRAVY + 0.15×|charge| (signed charge bonus
for |Q|>0). Higher score → better predicted solubility. No training data.
"""

from __future__ import annotations

from peptideforge.contracts.models import DevelopabilityProperty, DevelopabilityScores, PeptideCandidate
from peptideforge.developability.amino_acids import gravy, net_charge_ph74
from peptideforge.developability.base import single_score_result


class SolubilityPredictor:
    """GRAVY + net-charge composite solubility score."""

    property_name = DevelopabilityProperty.SOLUBILITY.value

    def __init__(self, charge_weight: float = 0.15) -> None:
        self.charge_weight = charge_weight

    def predict(self, candidate: PeptideCandidate) -> DevelopabilityScores:
        seq = candidate.sequence
        g = gravy(seq)
        q = net_charge_ph74(seq)
        # Higher is better: negative GRAVY helps; moderate |charge| helps solubility
        score = -g + self.charge_weight * abs(q)
        # Uncertainty grows with length (charge estimate coarser)
        uncertainty = 0.05 + 0.002 * len(seq)
        return single_score_result(
            candidate,
            property_name=DevelopabilityProperty.SOLUBILITY,
            value=score,
            uncertainty=uncertainty,
            higher_is_better=True,
            method="gravy_net_charge_ph74",
        )

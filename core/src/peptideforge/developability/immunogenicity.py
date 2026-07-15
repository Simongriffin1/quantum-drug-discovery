"""MHC-II immunogenicity heuristic (algorithmic — not NetMHCIIpan).

Method: scan all 9-mer windows; score each core using simplified anchor propensities
(P1, P4, P6, P9 positions in the 9-mer, Sturniolo-inspired heuristic). Return the
maximum window score as immunogenicity risk. Higher → more likely MHC-II presentation.

LICENSE NOTE: NetMHCIIpan / NetMHC family require commercial licenses for commercial
use. This module uses **no NetMHC weights** — only published-style propensity tables.
For production immunogenicity, plug in a licensed predictor behind the same interface.
"""

from __future__ import annotations

from peptideforge.contracts.models import DevelopabilityProperty, DevelopabilityScores, PeptideCandidate
from peptideforge.developability.amino_acids import MHC2_ANCHOR_PROPENSITY, validate_sequence
from peptideforge.developability.base import single_score_result

MHC2_CORE_LEN = 9
# Anchor indices within 9-mer (0-based): P1, P4, P6, P9
ANCHOR_IDX = (0, 3, 5, 8)


class ImmunogenicityPredictor:
    """MHC-II 9-mer propensity heuristic (no proprietary NetMHC weights)."""

    property_name = DevelopabilityProperty.IMMUNOGENICITY.value

    def predict(self, candidate: PeptideCandidate) -> DevelopabilityScores:
        seq = validate_sequence(candidate.sequence)
        if len(seq) < MHC2_CORE_LEN:
            # Very short: average propensity
            score = sum(MHC2_ANCHOR_PROPENSITY[aa] for aa in seq) / len(seq)
            uncertainty = 0.2
        else:
            window_scores: list[float] = []
            for i in range(len(seq) - MHC2_CORE_LEN + 1):
                core = seq[i : i + MHC2_CORE_LEN]
                anchors = [core[j] for j in ANCHOR_IDX]
                window_scores.append(sum(MHC2_ANCHOR_PROPENSITY[a] for a in anchors) / len(anchors))
            score = max(window_scores)
            uncertainty = 0.1
        # Normalize to ~[0,1]
        normalized = max(0.0, min(1.0, (score + 0.5) / 1.5))
        return single_score_result(
            candidate,
            property_name=DevelopabilityProperty.IMMUNOGENICITY,
            value=normalized,
            uncertainty=uncertainty,
            higher_is_better=False,
            method="mhc2_9mer_anchor_heuristic_not_netmhc",
        )

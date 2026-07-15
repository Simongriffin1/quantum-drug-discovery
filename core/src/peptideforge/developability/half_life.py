"""Half-life heuristics: N-end rule + protease cleavage-site density.

Method: mammalian N-end rule (Varshavsky/Bachmair) — destabilizing N-terminal residues
shorten cytosolic half-life. Protease cleavage motifs (trypsin-like K/R, chymotrypsin-
like F/W/Y/L) increase degradation susceptibility. Higher score → longer predicted
half-life. No training data.
"""

from __future__ import annotations

import re

from peptideforge.contracts.models import DevelopabilityProperty, DevelopabilityScores, PeptideCandidate
from peptideforge.developability.amino_acids import N_END_DESTABILIZING, N_END_STABILIZING, validate_sequence
from peptideforge.developability.base import single_score_result

# Trypsin: K or R not followed by P; chymotrypsin-like: after F/W/Y/L
TRYPSIN_PATTERN = re.compile(r"[KR](?!P)")
CHYMOTRYPSIN_PATTERN = re.compile(r"[FWYL]")


class HalfLifePredictor:
    """N-end rule + cleavage-site heuristic half-life score."""

    property_name = DevelopabilityProperty.HALF_LIFE.value

    def predict(self, candidate: PeptideCandidate) -> DevelopabilityScores:
        seq = validate_sequence(candidate.sequence)
        score = 0.5  # neutral baseline
        n_term = seq[0]
        if n_term in N_END_STABILIZING:
            score += 0.25
        elif n_term in N_END_DESTABILIZING:
            score -= 0.25
        # Cleavage sites reduce effective half-life
        n_trypsin = len(TRYPSIN_PATTERN.findall(seq))
        n_chymo = len(CHYMOTRYPSIN_PATTERN.findall(seq))
        cleavage_density = (n_trypsin + 0.5 * n_chymo) / len(seq)
        score -= min(0.4, cleavage_density * 2.0)
        score = max(0.0, min(1.0, score))
        uncertainty = 0.12
        return single_score_result(
            candidate,
            property_name=DevelopabilityProperty.HALF_LIFE,
            value=score,
            uncertainty=uncertainty,
            higher_is_better=True,
            method="n_end_rule_cleavage_heuristic",
        )

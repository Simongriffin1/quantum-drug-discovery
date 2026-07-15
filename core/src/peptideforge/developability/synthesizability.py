"""Synthesizability heuristics for solid-phase peptide synthesis (SPPS).

Method: penalize length, difficult residues (C/M/W), N-terminal Q/N (cyclization/
deamidation risk), high hydrophobic runs (on-resin aggregation), and Pro/Gly
stretches that complicate coupling. Higher score → easier / more reliable synthesis.
No training data.
"""

from __future__ import annotations

from peptideforge.contracts.models import DevelopabilityProperty, DevelopabilityScores, PeptideCandidate
from peptideforge.developability.amino_acids import (
    DIFFICULT_SYNTHESIS,
    KYTE_DOOLITTLE,
    N_TERM_RISK,
    validate_sequence,
)
from peptideforge.developability.base import single_score_result

MAX_LEN_COMFORT = 30
HYDROPHOBIC = frozenset("FILMVWY")


class SynthesizabilityPredictor:
    """Length + difficult-residue + aggregation-on-resin heuristics."""

    property_name = DevelopabilityProperty.SYNTHESIZABILITY.value

    def predict(self, candidate: PeptideCandidate) -> DevelopabilityScores:
        seq = validate_sequence(candidate.sequence)
        n = len(seq)
        score = 1.0
        # Length penalty beyond comfort zone
        if n > MAX_LEN_COMFORT:
            score -= 0.02 * (n - MAX_LEN_COMFORT)
        # Difficult residues
        difficult_frac = sum(1 for aa in seq if aa in DIFFICULT_SYNTHESIS) / n
        score -= 0.4 * difficult_frac
        # N-terminal risk
        if seq[0] in N_TERM_RISK:
            score -= 0.15
        # Long hydrophobic stretches (≥4) — resin aggregation during SPPS
        max_hydro_run = _max_run(seq, HYDROPHOBIC)
        if max_hydro_run >= 4:
            score -= 0.1 * (max_hydro_run - 3)
        # Very high average hydrophobicity
        avg_hydro = sum(KYTE_DOOLITTLE[aa] for aa in seq) / n
        if avg_hydro > 1.5:
            score -= 0.1 * (avg_hydro - 1.5)
        score = max(0.0, min(1.0, score))
        uncertainty = 0.05 + 0.001 * n
        return single_score_result(
            candidate,
            property_name=DevelopabilityProperty.SYNTHESIZABILITY,
            value=score,
            uncertainty=uncertainty,
            higher_is_better=True,
            method="spps_length_difficult_residue_heuristic",
        )


def _max_run(seq: str, chars: frozenset[str]) -> int:
    best = cur = 0
    for aa in seq:
        if aa in chars:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best

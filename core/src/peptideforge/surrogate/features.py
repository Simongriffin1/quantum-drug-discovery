"""Sequence features for surrogate models — published scales only, no training data.

Biological rationale: amino-acid composition + hydropathy + net charge are cheap
physicochemical descriptors that weakly predict affinity/developability until
physics labels arrive. They never replace the oracle.
"""

from __future__ import annotations

from peptideforge.contracts.models import CANONICAL_AA, PeptideCandidate
from peptideforge.developability.amino_acids import gravy, net_charge_ph74

# Fixed feature order for reproducibility across fit/predict.
AA_ORDER: tuple[str, ...] = tuple(sorted(CANONICAL_AA))
FEATURE_NAMES: tuple[str, ...] = (
    *(f"frac_{aa}" for aa in AA_ORDER),
    "length_norm",
    "gravy",
    "net_charge_ph74",
)


def sequence_features(sequence: str) -> tuple[float, ...]:
    """Return fixed-length feature vector for a peptide sequence."""
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        raise ValueError("empty sequence")
    counts = {aa: 0 for aa in AA_ORDER}
    for aa in seq:
        if aa not in counts:
            raise ValueError(f"non-canonical residue: {aa}")
        counts[aa] += 1
    fracs = tuple(counts[aa] / n for aa in AA_ORDER)
    # Normalize length to [0,1] over peptide range 5–50
    length_norm = (n - 5) / 45.0
    return (*fracs, length_norm, gravy(seq), net_charge_ph74(seq))


def candidate_features(candidate: PeptideCandidate) -> tuple[float, ...]:
    return sequence_features(candidate.sequence)

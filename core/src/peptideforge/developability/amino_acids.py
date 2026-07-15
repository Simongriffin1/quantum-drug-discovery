"""Amino-acid biophysical constants for developability scoring.

All scales are published, algorithmic lookups — no proprietary training data.
"""

from __future__ import annotations

# Kyte-Doolittle hydropathy (1982). Positive = hydrophobic.
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

# Side-chain pKa values (approx.) for net charge at pH 7.4.
PKA_SIDE: dict[str, float] = {
    "D": 3.9,
    "E": 4.3,
    "H": 6.0,
    "C": 8.3,
    "Y": 10.1,
    "K": 10.5,
    "R": 12.5,
}
PKA_N_TERM = 9.6
PKA_C_TERM = 2.3

# Residues flagged as difficult for SPPS (synthesis aggregation / coupling).
DIFFICULT_SYNTHESIS = frozenset("CMW")
# N-terminal glutamine cyclizes to pyroglutamate; asparagine deamidation-prone.
N_TERM_RISK = frozenset("QN")

# Mammalian N-end rule (Varshavsky): destabilizing N-terminal residues shorten half-life.
N_END_DESTABILIZING = frozenset("RKHFYWILM")  # broad destabilizing set for peptides
N_END_STABILIZING = frozenset("GSTAVCNDEQP")

# Simplified MHC-II pocket propensity (9-mer core positions P1,P4,P6,P9 anchors).
# Algorithmic heuristic inspired by Sturniolo et al.; not NetMHCIIpan weights.
MHC2_ANCHOR_PROPENSITY: dict[str, float] = {
    "A": 0.2,
    "C": 0.1,
    "D": -0.3,
    "E": -0.2,
    "F": 0.8,
    "G": 0.0,
    "H": 0.1,
    "I": 0.7,
    "K": -0.4,
    "L": 0.9,
    "M": 0.6,
    "N": -0.2,
    "P": -0.5,
    "Q": -0.1,
    "R": -0.3,
    "S": 0.0,
    "T": 0.1,
    "V": 0.5,
    "W": 0.7,
    "Y": 0.4,
}


def validate_sequence(seq: str) -> str:
    """Uppercase canonical sequence or raise."""
    s = seq.upper()
    invalid = {c for c in s if c not in KYTE_DOOLITTLE}
    if invalid:
        raise ValueError(f"non-canonical residues: {sorted(invalid)}")
    return s


def gravy(seq: str) -> float:
    """Grand average of hydropathy (Kyte-Doolittle). Lower → more soluble."""
    s = validate_sequence(seq)
    return sum(KYTE_DOOLITTLE[aa] for aa in s) / len(s)


def net_charge_ph74(seq: str) -> float:
    """Estimated net charge at pH 7.4 (Henderson-Hasselbalch, side chains + termini)."""
    s = validate_sequence(seq)
    charge = 0.0
    # N-terminus
    charge += 1.0 / (1.0 + 10 ** (7.4 - PKA_N_TERM))
    # C-terminus
    charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERM - 7.4))
    for aa in s:
        if aa not in PKA_SIDE:
            continue
        pka = PKA_SIDE[aa]
        if aa in "DE":
            charge -= 1.0 / (1.0 + 10 ** (PKA_SIDE[aa] - 7.4))
        elif aa in "KRH":
            charge += 1.0 / (1.0 + 10 ** (7.4 - PKA_SIDE[aa]))
        elif aa in "CY":
            # deprotonated thiol / phenol contributions (simplified)
            if aa == "C":
                charge -= 1.0 / (1.0 + 10 ** (PKA_SIDE[aa] - 7.4))
            else:
                charge -= 1.0 / (1.0 + 10 ** (PKA_SIDE[aa] - 7.4))
    return charge

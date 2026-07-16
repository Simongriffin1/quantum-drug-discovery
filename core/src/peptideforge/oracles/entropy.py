"""Configurational entropy estimates for MM-GBSA end-point scoring.

Supports a cheap truncated-NMA / fluctuation proxy on a minimized structure
(Interaction-Entropy-like −TΔS from energy variance when MD samples exist).
Flags are config-level; A/B on train/dev only — never invent a favorable ΔS.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class EntropyEstimate:
    """Entropy correction in kcal/mol (negative = favorable binding contribution)."""

    tds_kcal_mol: float
    method: str
    n_samples: int
    metadata: dict[str, Any]


def interaction_entropy_from_energies(
    energy_samples_kcal: list[float],
    *,
    temperature_K: float = 300.0,
) -> EntropyEstimate:
    """Interaction Entropy: −TΔS ≈ kT ln⟨exp(βΔE)⟩ − ⟨ΔE⟩ (Duan et al.).

    ``energy_samples_kcal`` are per-frame interaction energies (E_c − E_r − E_l).
    Requires ≥2 frames; raises if empty.
    """
    if len(energy_samples_kcal) < 2:
        raise ValueError(
            "interaction_entropy_from_energies needs ≥2 MD frames; "
            "refuse silent zero entropy."
        )
    import math

    kB = 0.0019872041  # kcal/(mol·K)
    beta = 1.0 / (kB * temperature_K)
    mean_e = sum(energy_samples_kcal) / len(energy_samples_kcal)
    # Shift for numerical stability
    shifted = [e - mean_e for e in energy_samples_kcal]
    exp_mean = sum(math.exp(beta * s) for s in shifted) / len(shifted)
    tds = kB * temperature_K * math.log(max(exp_mean, 1e-300))
    # IE convention: −TΔS ≈ kT ln⟨e^{βΔE}⟩ − ⟨ΔE⟩; with shift, −⟨ΔE⟩ cancels mean
    # Effective correction added to ΔG is +TΔS_config (unfavorable if fluctuations large)
    return EntropyEstimate(
        tds_kcal_mol=tds,
        method="interaction_entropy",
        n_samples=len(energy_samples_kcal),
        metadata={"temperature_K": temperature_K, "mean_E": mean_e},
    )


def truncated_nma_entropy_proxy(
    *,
    n_heavy_atoms_interface: int,
    temperature_K: float = 300.0,
) -> EntropyEstimate:
    """Cheap minimized-structure entropy proxy (not full NMA).

    Rough −TΔS ∝ √N_interface heavy atoms (order-of-magnitude placeholder for
    A/B gating). Documented as crude — prefer Interaction Entropy when MD runs.
    """
    if n_heavy_atoms_interface < 1:
        raise ValueError("n_heavy_atoms_interface must be ≥1")
    # Scale so ~50 heavy atoms → ~3 kcal/mol at 300 K (literature ballpark)
    tds = 0.4 * (n_heavy_atoms_interface**0.5) * (temperature_K / 300.0)
    return EntropyEstimate(
        tds_kcal_mol=tds,
        method="truncated_nma_proxy",
        n_samples=1,
        metadata={"n_heavy_atoms_interface": n_heavy_atoms_interface},
    )


def apply_entropy_correction(dg_kcal: float, entropy: EntropyEstimate) -> float:
    """ΔG_corrected = ΔG_MMGBSA + (−TΔS) with IE sign: add tds_kcal_mol."""
    return dg_kcal + entropy.tds_kcal_mol

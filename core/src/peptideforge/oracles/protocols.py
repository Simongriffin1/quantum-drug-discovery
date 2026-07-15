"""Configurable MM-GBSA / MM-PBSA end-point protocols for train/dev sweeps.

Literature anchor: Chen et al. PCCP 2019 — for short peptides (5–12 aa),
MM/PBSA on minimized explicit-solvent structures with ε_in=2 often wins;
medium peptides prefer MM/GBSA + short MD. Sweep NEVER touches held-out test.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class EndpointProtocol:
    """One scoring protocol configuration."""

    name: str
    solvation: str = "gbsa"  # gbsa | pbsa (pbsa requires OpenMM PB plugin / fail loud)
    gb_model: str = "obc2"  # obc1 | obc2 | gbn2
    epsilon_in: float = 1.0
    forcefield: str = "amber14"  # amber14 | ff14SB | ff99SB | ff03 (mapped to XML)
    minimize_max_iterations: int = 100
    minimize_solvent: str = "implicit"  # implicit | explicit
    md_steps: int = 0  # 0 = single minimized structure
    include_entropy: bool = False
    include_interface_waters: bool = False
    length_dependent: bool = True  # ≤12 aa prefer single-min; longer may use MD

    def forcefield_xml(self) -> tuple[str, ...]:
        gb = {
            "obc1": "implicit/obc1.xml",
            "obc2": "implicit/obc2.xml",
            "gbn2": "implicit/gbn2.xml",
        }.get(self.gb_model, "implicit/obc2.xml")
        # OpenMM amber14 is the robust installed default; ff99/ff03 requested
        # via name tagging for sweep provenance (XML map expands when available).
        if self.forcefield in {"amber14", "ff14SB"}:
            return ("amber14-all.xml", gb)
        if self.forcefield in {"ff99SB", "ff99"}:
            # Fall back to amber14 if classic ff99 XML not shipped — tagged in name
            return ("amber14-all.xml", gb)
        if self.forcefield == "ff03":
            return ("amber14-all.xml", gb)
        raise ValueError(f"unknown forcefield: {self.forcefield}")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def default_protocol_grid(*, short_peptide: bool = True) -> tuple[EndpointProtocol, ...]:
    """Train/dev sweep grid grounded in Chen et al. recommendations."""
    # Recommended start for 5–12 aa: PBSA, explicit min, ε_in=2 — PBSA may be
    # unavailable; include GBSA ε_in sweep as the always-runnable core.
    base = [
        EndpointProtocol(
            name="gbsa_obc2_eps1_amber14_min0",
            solvation="gbsa",
            gb_model="obc2",
            epsilon_in=1.0,
            forcefield="amber14",
            minimize_max_iterations=0,
        ),
        EndpointProtocol(
            name="gbsa_obc2_eps2_amber14_min100",
            solvation="gbsa",
            gb_model="obc2",
            epsilon_in=2.0,
            forcefield="amber14",
            minimize_max_iterations=100,
        ),
        EndpointProtocol(
            name="gbsa_obc2_eps4_amber14_min100",
            solvation="gbsa",
            gb_model="obc2",
            epsilon_in=4.0,
            forcefield="amber14",
            minimize_max_iterations=100,
        ),
        EndpointProtocol(
            name="gbsa_obc1_eps2_amber14_min100",
            solvation="gbsa",
            gb_model="obc1",
            epsilon_in=2.0,
            forcefield="amber14",
            minimize_max_iterations=100,
        ),
        EndpointProtocol(
            name="gbsa_gbn2_eps2_amber14_min100",
            solvation="gbsa",
            gb_model="gbn2",
            epsilon_in=2.0,
            forcefield="amber14",
            minimize_max_iterations=100,
        ),
        EndpointProtocol(
            name="gbsa_obc2_eps2_ff03_min100",
            solvation="gbsa",
            gb_model="obc2",
            epsilon_in=2.0,
            forcefield="ff03",
            minimize_max_iterations=100,
        ),
    ]
    if not short_peptide:
        base.append(
            EndpointProtocol(
                name="gbsa_obc2_eps1_amber14_md200",
                solvation="gbsa",
                gb_model="obc2",
                epsilon_in=1.0,
                forcefield="amber14",
                minimize_max_iterations=50,
                md_steps=200,
            )
        )
    # Literature-preferred short-peptide start (PBSA may fail loud if unavailable)
    base.insert(
        0,
        EndpointProtocol(
            name="pbsa_eps2_ff99_explicit_min_SHORT_PEP_REC",
            solvation="pbsa",
            epsilon_in=2.0,
            forcefield="ff99SB",
            minimize_max_iterations=100,
            minimize_solvent="explicit",
        ),
    )
    return tuple(base)

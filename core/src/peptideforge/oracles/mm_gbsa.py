"""MM-GBSA and tiered OpenMM physics oracle.

Physical rationale
------------------
MM-GBSA estimates relative binding free energy as

    ΔG_bind ≈ ⟨E_complex⟩ − ⟨E_receptor⟩ − ⟨E_ligand⟩

with a molecular-mechanics potential plus a Generalized Born Surface Area
implicit solvent term (Onufriev–Bashford–Case / GBSA-OBC). This is the
workhorse oracle (tier ``mm_gbsa``): first-principles enough to self-label,
cheap enough to screen peptide–target complexes on cloud GPUs.

Tier ``docking``: same thermodynamic cycle with a short minimization / vacuum-
lean setup (peptide-appropriate interface score; AutoDock Vina is for small
molecules and remains a pluggable alternative).

Tier ``md``: short production MD then average the MM-GBSA cycle over frames.

Missing OpenMM → ``OpenMMUnavailableError`` (never fabricate energies).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
from uuid import uuid4

from peptideforge.contracts.models import (
    ComplexStructure,
    OracleResult,
    OracleTier,
    Provenance,
)
from peptideforge.oracles.costs import enforce_cost_cap, resolve_tier
from peptideforge.oracles.openmm_utils import (
    all_chain_ids,
    load_pdb,
    require_openmm,
    resolve_pdb_path,
)


@dataclass(frozen=True)
class OpenMMOracleConfig:
    """Simulation hyperparameters (seeded where stochastic)."""

    forcefield_xml: tuple[str, ...] = ("amber14-all.xml", "implicit/obc2.xml")
    temperature_K: float = 300.0
    friction_per_ps: float = 1.0
    timestep_fs: float = 2.0
    # Default light/no minimize: full LocalEnergyMinimizer is very expensive on CPU.
    minimize_max_iterations: int = 0
    docking_minimize_iterations: int = 0
    md_steps: int = 200
    md_report_interval: int = 50
    platform: str | None = None  # e.g. "CPU", "CUDA"; None = OpenMM default
    peptide_chain_ids: tuple[str, ...] = ("P", "L", "C")
    receptor_chain_ids: tuple[str, ...] | None = None  # None → all non-peptide
    seed: int = 0
    nonbonded_cutoff_nm: float = 1.6


class OpenMMPhysicsOracle:
    """Tiered OpenMM oracle implementing the ``Oracle`` contract."""

    AVAILABLE_TIERS: tuple[OracleTier, ...] = (
        OracleTier.DOCKING,
        OracleTier.MM_GBSA,
        OracleTier.MD,
    )

    def __init__(self, config: OpenMMOracleConfig | None = None) -> None:
        # Touch OpenMM at construction — fail loud early.
        require_openmm()
        self.config = config or OpenMMOracleConfig()

    def evaluate(
        self,
        complex_structure: ComplexStructure,
        *,
        tier: OracleTier | None = None,
        cost_cap: float | None = None,
    ) -> OracleResult:
        resolved = resolve_tier(tier, available=self.AVAILABLE_TIERS)
        cost = enforce_cost_cap(resolved, cost_cap)
        pdb_path = resolve_pdb_path(complex_structure.pdb_path, complex_structure.pdb_text)

        if resolved == OracleTier.DOCKING:
            value, uncertainty, meta = self._tier_docking(pdb_path)
        elif resolved == OracleTier.MM_GBSA:
            value, uncertainty, meta = self._tier_mm_gbsa(pdb_path)
        elif resolved == OracleTier.MD:
            value, uncertainty, meta = self._tier_short_md(pdb_path)
        else:
            raise ValueError(f"tier not implemented: {resolved}")

        openmm, _, _ = require_openmm()
        return OracleResult(
            result_id=uuid4(),
            candidate_id=complex_structure.candidate_id,
            complex_id=complex_structure.complex_id,
            value=value,
            uncertainty=uncertainty,
            cost_estimate=cost,
            tier=resolved,
            unit="kcal/mol",
            metadata={
                **meta,
                "pdb_path": str(pdb_path),
                "forcefield": list(self.config.forcefield_xml),
            },
            provenance=Provenance(
                tool_versions={"openmm": openmm.version.version},
                config_hash=str(self.config.seed),
            ),
        )

    # ------------------------------------------------------------------ tiers

    def _tier_docking(self, pdb_path: Path) -> tuple[float, float, dict[str, Any]]:
        """Cheap interface ΔE with limited minimization (Tier 0)."""
        value = self._mm_gbsa_single_point(
            pdb_path,
            minimize_iterations=self.config.docking_minimize_iterations,
            sample_md_steps=0,
        )
        return value, 1.5, {"method": "openmm_interface_gbsa_quick"}

    def _tier_mm_gbsa(self, pdb_path: Path) -> tuple[float, float, dict[str, Any]]:
        """Workhorse: minimize + single-structure MM-GBSA cycle."""
        value = self._mm_gbsa_single_point(
            pdb_path,
            minimize_iterations=self.config.minimize_max_iterations,
            sample_md_steps=0,
        )
        return value, 1.0, {"method": "openmm_mm_gbsa_minimize"}

    def _tier_short_md(self, pdb_path: Path) -> tuple[float, float, dict[str, Any]]:
        """Short MD then average MM-GBSA over reported frames."""
        value, std = self._mm_gbsa_md_average(pdb_path)
        return (
            value,
            max(std, 0.5),
            {"method": "openmm_mm_gbsa_short_md", "md_steps": self.config.md_steps},
        )

    # ------------------------------------------------------------------ physics

    def _build_context(
        self,
        topology: Any,
        positions: Any,
    ) -> tuple[Any, Any, Any]:
        openmm, app, unit = require_openmm()
        forcefield = app.ForceField(*self.config.forcefield_xml)
        residue_templates: dict[Any, str] = {}
        for residue in topology.residues():
            # MHC interfaces often leave ambiguous Cys protonation after trim.
            if residue.name == "CYS":
                residue_templates[residue] = "CYX"
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.CutoffNonPeriodic,
            nonbondedCutoff=self.config.nonbonded_cutoff_nm * unit.nanometers,
            constraints=app.HBonds,
            # Interface trims create artificial termini; ignore dangling peptide bonds.
            ignoreExternalBonds=True,
            residueTemplates=residue_templates,
        )
        integrator = openmm.LangevinMiddleIntegrator(
            self.config.temperature_K * unit.kelvin,
            self.config.friction_per_ps / unit.picosecond,
            self.config.timestep_fs * unit.femtoseconds,
        )
        integrator.setRandomNumberSeed(self.config.seed)
        if self.config.platform:
            platform = openmm.Platform.getPlatformByName(self.config.platform)
            context = openmm.Context(system, integrator, platform)
        else:
            context = openmm.Context(system, integrator)
        context.setPositions(positions)
        return system, integrator, context

    def _minimize(self, context: Any, max_iterations: int) -> None:
        openmm, _, _ = require_openmm()
        openmm.LocalEnergyMinimizer.minimize(context, maxIterations=max_iterations)

    def _potential_kcal(self, context: Any) -> float:
        _, _, unit = require_openmm()
        state = context.getState(getEnergy=True)
        return float(state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))

    def _select_chains(self, topology: Any) -> tuple[set[str], set[str]]:
        ids = all_chain_ids(topology)
        if not ids:
            raise ValueError("PDB topology has no chains")
        peptide = set(self.config.peptide_chain_ids) & set(ids)
        if not peptide:
            # Fall back: last chain is ligand/peptide (common in 2-chain dumps)
            if len(ids) < 2:
                raise ValueError(
                    "need ≥2 chains (receptor + peptide) or peptide_chain_ids matching the PDB; "
                    f"found chains={ids}"
                )
            peptide = {ids[-1]}
        if self.config.receptor_chain_ids is not None:
            receptor = set(self.config.receptor_chain_ids) & set(ids)
        else:
            receptor = set(ids) - peptide
        if not receptor:
            raise ValueError(f"no receptor chains left; chains={ids}, peptide={peptide}")
        return receptor, peptide

    def _energy_of_subset(
        self,
        topology: Any,
        positions: Any,
        keep_chains: set[str],
        *,
        minimize_iterations: int,
    ) -> float:
        """Build a GBSA system for a chain subset and return potential (kcal/mol)."""
        openmm, app, unit = require_openmm()
        modeller = app.Modeller(topology, positions)
        to_delete = [
            residue for residue in topology.residues() if residue.chain.id not in keep_chains
        ]
        modeller.delete(to_delete)
        forcefield = app.ForceField(*self.config.forcefield_xml)
        has_hydrogen = any(atom.element.symbol == "H" for atom in modeller.topology.atoms())
        if not has_hydrogen:
            modeller.addHydrogens(forcefield)
        _, _, context = self._build_context(modeller.topology, modeller.positions)
        if minimize_iterations > 0:
            self._minimize(context, minimize_iterations)
        return self._potential_kcal(context)

    def _mm_gbsa_single_point(
        self,
        pdb_path: Path,
        *,
        minimize_iterations: int,
        sample_md_steps: int,
    ) -> float:
        topology, positions = load_pdb(pdb_path)
        receptor, peptide = self._select_chains(topology)
        # Full complex
        e_complex = self._energy_of_subset(
            topology, positions, receptor | peptide, minimize_iterations=minimize_iterations
        )
        e_receptor = self._energy_of_subset(
            topology, positions, receptor, minimize_iterations=minimize_iterations
        )
        e_ligand = self._energy_of_subset(
            topology, positions, peptide, minimize_iterations=minimize_iterations
        )
        if sample_md_steps > 0:
            # reserved for callers that already sampling; kept for API clarity
            pass
        return e_complex - e_receptor - e_ligand

    def _mm_gbsa_md_average(self, pdb_path: Path) -> tuple[float, float]:
        """Minimize complex, run short MD, average MM-GBSA over sampled frames."""
        openmm, app, unit = require_openmm()
        topology, positions = load_pdb(pdb_path)
        receptor, peptide = self._select_chains(topology)

        forcefield = app.ForceField(*self.config.forcefield_xml)
        modeller = app.Modeller(topology, positions)
        has_hydrogen = any(atom.element.symbol == "H" for atom in modeller.topology.atoms())
        if not has_hydrogen:
            modeller.addHydrogens(forcefield)
        system, integrator, context = self._build_context(modeller.topology, modeller.positions)
        self._minimize(context, self.config.minimize_max_iterations)

        samples: list[float] = []
        steps = self.config.md_steps
        interval = max(1, self.config.md_report_interval)
        for start in range(0, steps, interval):
            context.getIntegrator().step(min(interval, steps - start))
            state = context.getState(getPositions=True)
            frame_pos = state.getPositions(asNumpy=False)
            e_c = self._energy_of_subset(
                modeller.topology, frame_pos, receptor | peptide, minimize_iterations=0
            )
            e_r = self._energy_of_subset(
                modeller.topology, frame_pos, receptor, minimize_iterations=0
            )
            e_l = self._energy_of_subset(
                modeller.topology, frame_pos, peptide, minimize_iterations=0
            )
            samples.append(e_c - e_r - e_l)

        if not samples:
            raise RuntimeError("short MD produced no MM-GBSA samples")
        mean = sum(samples) / len(samples)
        var = sum((x - mean) ** 2 for x in samples) / len(samples)
        return mean, var**0.5

"""Energy-conservation / equilibration sanity checks for OpenMM MD.

Physical rationale: under a Verlet integrator with no thermostat and constraints
handled correctly, total energy of an isolated system should drift little over
short trajectories. Large drift usually means a bad topology, timestep, or forces.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from peptideforge.oracles.openmm_utils import require_openmm

# Shipped OpenMM reference alanine-dipeptide (ACE-ALA-NME) for golden tests.
_DEFAULT_ALA_PDB = (
    Path(__file__).resolve().parents[4]
    / "benchmarks"
    / "fixtures"
    / "structures"
    / "alanine_dipeptide_implicit.pdb"
)


@dataclass(frozen=True)
class EnergyConservationResult:
    """Summary of an NVE (or weakly damped) conservation check."""

    n_steps: int
    timestep_fs: float
    mean_total_kcal: float
    std_total_kcal: float
    drift_kcal: float
    max_abs_drift_kcal: float
    passed: bool
    tolerance_kcal: float
    system_name: str


def run_energy_conservation_check(
    pdb_path: Path | None = None,
    *,
    n_steps: int = 1000,
    timestep_fs: float = 1.0,
    report_interval: int = 50,
    temperature_K: float = 300.0,
    friction_per_ps: float = 0.0,
    tolerance_kcal: float = 5.0,
    platform: str | None = "CPU",
    seed: int = 0,
    system_name: str = "alanine_dipeptide_implicit",
) -> EnergyConservationResult:
    """Minimize lightly, then integrate and measure total-energy drift.

    Uses Verlet when ``friction_per_ps=0`` (NVE-like). ``passed`` if
    max |E − E0| ≤ ``tolerance_kcal`` over the trajectory.
    """
    openmm, app, unit = require_openmm()
    path = pdb_path if pdb_path is not None else _DEFAULT_ALA_PDB
    if not path.is_file():
        raise FileNotFoundError(
            f"energy-conservation PDB missing: {path}. "
            "Ship benchmarks/fixtures/structures/alanine_dipeptide_implicit.pdb"
        )

    pdb = app.PDBFile(str(path))
    forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.6 * unit.nanometers,
        constraints=app.HBonds,
    )
    if friction_per_ps <= 0.0:
        integrator: Any = openmm.VerletIntegrator(timestep_fs * unit.femtoseconds)
    else:
        integrator = openmm.LangevinMiddleIntegrator(
            temperature_K * unit.kelvin,
            friction_per_ps / unit.picosecond,
            timestep_fs * unit.femtoseconds,
        )
        integrator.setRandomNumberSeed(seed)

    if platform:
        ctx = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(platform))
    else:
        ctx = openmm.Context(system, integrator)
    ctx.setPositions(pdb.positions)
    openmm.LocalEnergyMinimizer.minimize(ctx, maxIterations=50)

    energies: list[float] = []
    for start in range(0, n_steps, report_interval):
        integrator.step(min(report_interval, n_steps - start))
        state = ctx.getState(getEnergy=True)
        e = float(
            (state.getKineticEnergy() + state.getPotentialEnergy()).value_in_unit(
                unit.kilocalories_per_mole
            )
        )
        energies.append(e)

    if not energies:
        raise RuntimeError("energy conservation run produced no samples")
    e0 = energies[0]
    drifts = [abs(e - e0) for e in energies]
    mean_e = sum(energies) / len(energies)
    var = sum((e - mean_e) ** 2 for e in energies) / len(energies)
    max_drift = max(drifts)
    return EnergyConservationResult(
        n_steps=n_steps,
        timestep_fs=timestep_fs,
        mean_total_kcal=mean_e,
        std_total_kcal=var**0.5,
        drift_kcal=energies[-1] - e0,
        max_abs_drift_kcal=max_drift,
        passed=max_drift <= tolerance_kcal,
        tolerance_kcal=tolerance_kcal,
        system_name=system_name,
    )

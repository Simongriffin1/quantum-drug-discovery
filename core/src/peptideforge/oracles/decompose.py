"""Per-component MM-GBSA energy decomposition (OpenMM force-group / scaling).

Method
------
For complex / receptor / peptide subsets we evaluate single-point energies under
amber14 + OBC2, then form Δ = E_c − E_r − E_l for each component:

1. **total** — full system potential (current oracle).
2. **vdW** — NonbondedForce with all charges set to 0; CustomGBForce disabled.
3. **elec** — NonbondedForce with LJ ε set to 0; CustomGBForce disabled
   (vacuum Coulomb).
4. **polar_solv** — CustomGBForce only (GB polar solvation).
5. **nonpolar** — residual: total − (bonded + vdW + elec + polar_solv), which
   captures SASA / other leftover terms when present.

Bonded (bonds/angles/torsions) cancel approximately in ΔG for rigid poses and
are reported for diagnostics only.

This uses OpenMM directly — not AmberTools MMPBSA.py / gmx_MMPBSA.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from peptideforge.eval.affinity_validity import charge_from_sequence
from peptideforge.eval.metrics import bootstrap_ci, spearman_rho
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle
from peptideforge.oracles.openmm_utils import all_chain_ids, load_pdb, require_openmm


@dataclass(frozen=True)
class ComponentEnergies:
    record_id: str
    total: float
    vdw: float
    elec: float
    polar_solv: float
    nonpolar: float
    elec_plus_polar: float
    vdw_plus_nonpolar: float
    bonded: float
    experimental_pkd: float
    net_charge: float


def _prepare_modeller(
    topology: Any,
    positions: Any,
    keep_chains: set[str],
    forcefield: Any,
) -> Any:
    _, app, _ = require_openmm()
    modeller = app.Modeller(topology, positions)
    to_delete = [r for r in topology.residues() if r.chain.id not in keep_chains]
    modeller.delete(to_delete)
    has_h = any(a.element is not None and a.element.symbol == "H" for a in modeller.topology.atoms())
    if not has_h:
        modeller.addHydrogens(forcefield)
    return modeller


def _build_system(topology: Any, config: OpenMMOracleConfig) -> Any:
    openmm, app, unit = require_openmm()
    forcefield = app.ForceField(*config.forcefield_xml)
    residue_templates: dict[Any, str] = {}
    for residue in topology.residues():
        if residue.name == "CYX":
            residue_templates[residue] = "CYX"
        elif residue.name == "CYS":
            sg = [a for a in residue.atoms() if a.name.strip() == "SG"]
            if sg:
                for bond in topology.bonds():
                    atoms = {bond.atom1, bond.atom2}
                    if sg[0] in atoms:
                        other = bond.atom2 if bond.atom1 == sg[0] else bond.atom1
                        if other.name.strip() == "SG":
                            residue_templates[residue] = "CYX"
                            break
    return forcefield.createSystem(
        topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=config.nonbonded_cutoff_nm * unit.nanometers,
        constraints=app.HBonds,
        ignoreExternalBonds=True,
        residueTemplates=residue_templates,
    ), forcefield


def _energy_kcal(system: Any, topology: Any, positions: Any, *, platform: str | None) -> float:
    openmm, _, unit = require_openmm()
    integrator = openmm.VerletIntegrator(0.001)
    if platform:
        ctx = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(platform))
    else:
        ctx = openmm.Context(system, integrator)
    ctx.setPositions(positions)
    return float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))


def _clone_system(system: Any) -> Any:
    openmm, _, _ = require_openmm()
    return openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system))


def _disable_force_types(system: Any, *, disable_gb: bool = False, disable_nonbonded: bool = False) -> None:
    openmm, _, _ = require_openmm()
    for force in system.getForces():
        name = type(force).__name__
        if disable_gb and "CustomGB" in name:
            force.setForceGroup(31)
            # Zero contribution by removing particles is heavy; instead scale via group mask
        if disable_nonbonded and name == "NonbondedForce":
            force.setForceGroup(31)


def _potential_groups(system: Any, topology: Any, positions: Any, groups: int, platform: str | None) -> float:
    """Energy from selected force groups (bitmask)."""
    openmm, _, unit = require_openmm()
    integrator = openmm.VerletIntegrator(0.001)
    if platform:
        ctx = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(platform))
    else:
        ctx = openmm.Context(system, integrator)
    ctx.setPositions(positions)
    return float(
        ctx.getState(getEnergy=True, groups=groups).getPotentialEnergy().value_in_unit(
            unit.kilocalories_per_mole
        )
    )


def _assign_groups(system: Any) -> dict[str, int]:
    """Map component → force group id."""
    openmm, _, _ = require_openmm()
    mapping = {"bonded": 0, "nonbonded": 1, "gb": 2, "other": 3}
    for force in system.getForces():
        name = type(force).__name__
        if name in {"HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce", "CMMotionRemover"}:
            force.setForceGroup(mapping["bonded"])
        elif name == "NonbondedForce":
            force.setForceGroup(mapping["nonbonded"])
        elif "CustomGB" in name or "GBSA" in name:
            force.setForceGroup(mapping["gb"])
        else:
            force.setForceGroup(mapping["other"])
    return mapping


def _mutate_nonbonded(
    system: Any,
    *,
    zero_charge: bool = False,
    zero_lj: bool = False,
) -> None:
    openmm, _, unit = require_openmm()
    for force in system.getForces():
        if type(force).__name__ != "NonbondedForce":
            continue
        for i in range(force.getNumParticles()):
            charge, sigma, eps = force.getParticleParameters(i)
            if zero_charge:
                charge = 0.0 * unit.elementary_charge
            if zero_lj:
                eps = 0.0 * unit.kilojoules_per_mole
                sigma = 1.0 * unit.nanometer
            force.setParticleParameters(i, charge, sigma, eps)
        for i in range(force.getNumExceptions()):
            a, b, charge_prod, sigma, eps = force.getExceptionParameters(i)
            if zero_charge:
                charge_prod = 0.0 * unit.elementary_charge**2
            if zero_lj:
                eps = 0.0 * unit.kilojoules_per_mole
            force.setExceptionParameters(i, a, b, charge_prod, sigma, eps)


def component_energies_for_structure(
    pdb_path: Path,
    *,
    peptide_chain: str,
    config: OpenMMOracleConfig | None = None,
    platform: str | None = "CPU",
) -> dict[str, float]:
    """Return Δ component energies (complex − receptor − peptide) in kcal/mol."""
    config = config or OpenMMOracleConfig(minimize_max_iterations=0, platform=platform)
    openmm, app, unit = require_openmm()
    topology, positions = load_pdb(pdb_path)
    oracle = OpenMMPhysicsOracle(config)
    receptor, peptide = oracle._select_chains(topology)  # noqa: SLF001 — shared chain logic
    # Ensure peptide_chain preferred if present
    if peptide_chain in all_chain_ids(topology):
        peptide = {peptide_chain}
        receptor = set(all_chain_ids(topology)) - peptide

    forcefield = app.ForceField(*config.forcefield_xml)

    def subset_components(keep: set[str]) -> dict[str, float]:
        modeller = _prepare_modeller(topology, positions, keep, forcefield)
        system, _ = _build_system(modeller.topology, config)
        _assign_groups(system)
        pos = modeller.positions
        total = _energy_kcal(system, modeller.topology, pos, platform=platform)

        sys_vdw = _clone_system(system)
        _assign_groups(sys_vdw)
        _mutate_nonbonded(sys_vdw, zero_charge=True)
        # disable GB for pure vdW
        for force in sys_vdw.getForces():
            if "CustomGB" in type(force).__name__:
                force.setForceGroup(31)
        vdw = _potential_groups(sys_vdw, modeller.topology, pos, groups=1 << 1, platform=platform)

        sys_elec = _clone_system(system)
        _assign_groups(sys_elec)
        _mutate_nonbonded(sys_elec, zero_lj=True)
        for force in sys_elec.getForces():
            if "CustomGB" in type(force).__name__:
                force.setForceGroup(31)
        elec = _potential_groups(sys_elec, modeller.topology, pos, groups=1 << 1, platform=platform)

        polar = _potential_groups(system, modeller.topology, pos, groups=1 << 2, platform=platform)
        bonded = _potential_groups(system, modeller.topology, pos, groups=1 << 0, platform=platform)
        other = _potential_groups(system, modeller.topology, pos, groups=1 << 3, platform=platform)
        nonpolar = total - (bonded + vdw + elec + polar + other)
        return {
            "total": total,
            "vdw": vdw,
            "elec": elec,
            "polar_solv": polar,
            "nonpolar": nonpolar,
            "bonded": bonded,
            "other": other,
        }

    c = subset_components(receptor | peptide)
    r = subset_components(receptor)
    l = subset_components(peptide)
    out = {k: c[k] - r[k] - l[k] for k in c}
    out["elec_plus_polar"] = out["elec"] + out["polar_solv"]
    out["vdw_plus_nonpolar"] = out["vdw"] + out["nonpolar"]
    return out


def decompose_partition(
    record_ids: list[str],
    *,
    catalog: dict[str, dict[str, str]],
    pdb_by_id: dict[str, str],
    platform: str = "CPU",
) -> list[ComponentEnergies]:
    rows: list[ComponentEnergies] = []
    for rid in record_ids:
        if rid not in catalog or rid not in pdb_by_id:
            continue
        entry = catalog[rid]
        try:
            comps = component_energies_for_structure(
                Path(pdb_by_id[rid]),
                peptide_chain=entry.get("peptide_chain") or "C",
                platform=platform,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"decompose skip {rid}: {exc}")
            continue
        seq = entry["peptide_seq"]
        rows.append(
            ComponentEnergies(
                record_id=rid,
                total=comps["total"],
                vdw=comps["vdw"],
                elec=comps["elec"],
                polar_solv=comps["polar_solv"],
                nonpolar=comps["nonpolar"],
                elec_plus_polar=comps["elec_plus_polar"],
                vdw_plus_nonpolar=comps["vdw_plus_nonpolar"],
                bonded=comps["bonded"],
                experimental_pkd=float(entry["pKd"]),
                net_charge=charge_from_sequence(seq),
            )
        )
    return rows


def correlate_components(
    rows: list[ComponentEnergies],
    *,
    n_bootstrap: int = 500,
    seed: int = 0,
) -> dict[str, dict[str, float]]:
    if len(rows) < 3:
        raise ValueError(f"need ≥3 decomposed complexes, got {len(rows)}")
    y = [r.experimental_pkd for r in rows]
    # Predictions: more negative energy → tighter bind → invert for Spearman vs pKd
    terms = {
        "neg_total": [-r.total for r in rows],
        "neg_vdw": [-r.vdw for r in rows],
        "neg_elec": [-r.elec for r in rows],
        "neg_polar_solv": [-r.polar_solv for r in rows],
        "neg_nonpolar": [-r.nonpolar for r in rows],
        "neg_elec_plus_polar": [-r.elec_plus_polar for r in rows],
        "neg_vdw_plus_nonpolar": [-r.vdw_plus_nonpolar for r in rows],
        "net_charge": [r.net_charge for r in rows],
    }
    out: dict[str, dict[str, float]] = {}
    for name, x in terms.items():
        try:
            rho, lo, hi = bootstrap_ci(x, y, statistic="spearman", n_resamples=n_bootstrap, seed=seed)
            r, rlo, rhi = bootstrap_ci(x, y, statistic="pearson", n_resamples=n_bootstrap, seed=seed + 1)
            out[name] = {
                "spearman": rho,
                "spearman_ci_low": lo,
                "spearman_ci_high": hi,
                "pearson": r,
                "pearson_ci_low": rlo,
                "pearson_ci_high": rhi,
            }
        except Exception as exc:  # noqa: BLE001
            out[name] = {"error": str(exc)}  # type: ignore[dict-item]
        seed += 2
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--prep-manifest", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/peptide_affinity/data/decomposition_traindev.json"),
    )
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--n-bootstrap", type=int, default=500)
    args = parser.parse_args()

    splits = json.loads(args.splits.read_text(encoding="utf-8"))
    cold = splits["cold_start"]
    traindev = list(cold["train_ids"]) + list(cold.get("val_ids") or [])
    test_ids = set(cold["test_ids"])
    if set(traindev) & test_ids:
        raise SystemExit("REFUSING: train/dev overlaps test")

    catalog = {
        r["record_id"]: r
        for r in csv.DictReader(args.catalog.open(encoding="utf-8"), delimiter="\t")
    }
    pdb_by_id: dict[str, str] = {}
    for row in csv.DictReader(args.prep_manifest.open(encoding="utf-8"), delimiter="\t"):
        if row.get("ok") in {"1", "true", "True"}:
            pdb_by_id[row["record_id"]] = row.get("scoreable_path") or row["complex_path"]

    rows = decompose_partition(
        traindev, catalog=catalog, pdb_by_id=pdb_by_id, platform=args.platform
    )
    corr = correlate_components(rows, n_bootstrap=args.n_bootstrap)
    payload = {
        "method": "OpenMM force-group + charge/LJ scaling (not AmberTools MMPBSA.py)",
        "partition": "train+val",
        "n": len(rows),
        "test_touched": False,
        "correlations": corr,
        "rows": [asdict(r) for r in rows],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {args.out} N={len(rows)}")
    for name, stats in corr.items():
        if "spearman" in stats:
            print(
                f"  {name}: ρ={stats['spearman']:.3f} "
                f"CI=[{stats['spearman_ci_low']:.3f},{stats['spearman_ci_high']:.3f}]"
            )


if __name__ == "__main__":
    main()

"""Structure preparation for experimental peptide–protein complexes.

Pipeline: PDBFixer (missing atoms/hydrogens) → optional PDB2PQR/PROPKA at pH 7.4
→ strip crystallographic waters except flagged interface waters → split
receptor / peptide / complex. Failures are logged and excluded (never silently patched).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

from peptide_affinity.load import DATA_DIR, load_peptide_affinity_catalog

AA3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


@dataclass(frozen=True)
class PrepResult:
    record_id: str
    ok: bool
    complex_path: str | None
    receptor_path: str | None
    peptide_path: str | None
    interface_water_path: str | None
    error: str | None
    n_interface_waters: int = 0


def _require_openmm_pdbfixer() -> tuple[object, object, object, type]:
    try:
        import openmm
        from openmm import app, unit
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise ImportError(
            "Structure prep requires openmm + pdbfixer. Install with:\n"
            "  poetry run pip install openmm pdbfixer\n"
            "Refusing to fabricate prepared coordinates."
        ) from exc
    return openmm, app, unit, PDBFixer


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def prepare_complex(
    pdb_path: Path,
    *,
    out_dir: Path,
    record_id: str,
    peptide_chain: str,
    receptor_chains: tuple[str, ...],
    ph: float = 7.4,
    interface_water_cutoff_A: float = 3.5,
    use_pdb2pqr: bool = False,
    trim_cutoff_nm: float = 1.2,
) -> PrepResult:
    """Prepare one experimental complex; on failure return ok=False with error logged."""
    openmm, app, unit, PDBFixer = _require_openmm_pdbfixer()
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        fixer = PDBFixer(filename=str(pdb_path))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=True)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)

        # Optional PDB2PQR — fail loud if requested but unavailable
        work_pdb = out_dir / f"{record_id}_fixed.pdb"
        app.PDBFile.writeFile(fixer.topology, fixer.positions, open(work_pdb, "w"))
        if use_pdb2pqr:
            try:
                from pdb2pqr.main import main as pdb2pqr_main  # type: ignore
            except ImportError as exc:
                raise ImportError(
                    "use_pdb2pqr=True but pdb2pqr is not installed. "
                    "pip install pdb2pqr  OR re-run with --no-pdb2pqr."
                ) from exc
            pqr_out = out_dir / f"{record_id}.pqr"
            pdb2pqr_main(
                [
                    f"--ff=AMBER",
                    f"--with-ph={ph}",
                    str(work_pdb),
                    str(pqr_out),
                ]
            )

        pdb = app.PDBFile(str(work_pdb))
        topology = pdb.topology
        positions = pdb.positions

        # Identify interface waters (O of HOH within cutoff of peptide heavy atoms)
        pep_coords: list[tuple[float, float, float]] = []
        for atom in topology.atoms():
            if atom.residue.chain.id != peptide_chain:
                continue
            if atom.element is not None and atom.element.symbol == "H":
                continue
            p = positions[atom.index]
            pep_coords.append((p.x, p.y, p.z))
        if not pep_coords:
            raise ValueError(f"no peptide atoms on chain {peptide_chain}")

        interface_water_res: set[int] = set()
        for residue in topology.residues():
            if residue.name not in {"HOH", "WAT", "H2O"}:
                continue
            for atom in residue.atoms():
                if atom.name.strip() not in {"O", "OW"}:
                    continue
                p = positions[atom.index]
                pt = (p.x, p.y, p.z)
                if any(_dist(pt, pc) * 10.0 <= interface_water_cutoff_A for pc in pep_coords):
                    interface_water_res.add(residue.index)
                break

        keep_chains = set(receptor_chains) | {peptide_chain}
        # Trim receptor: keep residues with COM within trim_cutoff of peptide
        from openmm import Vec3

        pos_nm = [Vec3(float(p.x), float(p.y), float(p.z)) for p in positions]
        pep_coms = [
            pos_nm[a.index]
            for a in topology.atoms()
            if a.residue.chain.id == peptide_chain and a.element and a.element.symbol != "H"
        ]

        def res_com(residue: object) -> Vec3 | None:
            coords = [pos_nm[a.index] for a in residue.atoms()]  # type: ignore[attr-defined]
            if not coords:
                return None
            n = len(coords)
            return Vec3(
                sum(c.x for c in coords) / n,
                sum(c.y for c in coords) / n,
                sum(c.z for c in coords) / n,
            )

        keep_res: set[int] = set()
        for residue in topology.residues():
            cid = residue.chain.id
            if residue.name in {"HOH", "WAT", "H2O"}:
                if residue.index in interface_water_res:
                    keep_res.add(residue.index)
                continue
            if cid == peptide_chain:
                keep_res.add(residue.index)
                continue
            if cid not in receptor_chains:
                continue
            com = res_com(residue)
            if com is None:
                continue
            if any(
                math.sqrt((com.x - p.x) ** 2 + (com.y - p.y) ** 2 + (com.z - p.z) ** 2)
                <= trim_cutoff_nm
                for p in pep_coms
            ):
                keep_res.add(residue.index)

        # Expand to contiguous receptor window per chain
        for cid in receptor_chains:
            residues = [r for r in topology.residues() if r.chain.id == cid]
            idxs = [i for i, r in enumerate(residues) if r.index in keep_res]
            if not idxs:
                continue
            lo, hi = min(idxs), max(idxs)
            for i in range(lo, hi + 1):
                keep_res.add(residues[i].index)

        modeller = app.Modeller(topology, positions)
        to_delete = [r for r in modeller.topology.residues() if r.index not in keep_res]
        modeller.delete(to_delete)

        # QC: broken peptide chain or empty receptor
        pep_res = [r for r in modeller.topology.residues() if r.chain.id == peptide_chain]
        rec_res = [
            r
            for r in modeller.topology.residues()
            if r.chain.id in receptor_chains and r.name not in {"HOH", "WAT", "H2O"}
        ]
        if len(pep_res) < 5:
            raise ValueError(f"peptide too short after prep: {len(pep_res)} residues")
        if not rec_res:
            raise ValueError("receptor empty after trim — QC fail")

        complex_path = out_dir / f"{record_id}_complex.pdb"
        with complex_path.open("w", encoding="utf-8") as handle:
            handle.write(f"REMARK   1 PeptideForge prep record_id={record_id} ph={ph}\n")
            handle.write(
                f"REMARK   1 interface_waters={len(interface_water_res)} "
                f"cutoff_A={interface_water_cutoff_A}\n"
            )
            app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)

        # Split receptor / peptide (no waters in peptide file)
        def write_subset(path: Path, chain_pred: object) -> None:
            subset = app.Modeller(modeller.topology, modeller.positions)
            delete = [r for r in subset.topology.residues() if not chain_pred(r)]  # type: ignore[operator]
            subset.delete(delete)
            with path.open("w", encoding="utf-8") as handle:
                app.PDBFile.writeFile(subset.topology, subset.positions, handle)

        receptor_path = out_dir / f"{record_id}_receptor.pdb"
        peptide_path = out_dir / f"{record_id}_peptide.pdb"
        water_path = out_dir / f"{record_id}_interface_waters.pdb"
        write_subset(
            receptor_path,
            lambda r: r.chain.id in receptor_chains and r.name not in {"HOH", "WAT", "H2O"},
        )
        write_subset(
            peptide_path,
            lambda r: r.chain.id == peptide_chain and r.name not in {"HOH", "WAT", "H2O"},
        )
        write_subset(water_path, lambda r: r.index in interface_water_res)

        return PrepResult(
            record_id=record_id,
            ok=True,
            complex_path=str(complex_path),
            receptor_path=str(receptor_path),
            peptide_path=str(peptide_path),
            interface_water_path=str(water_path),
            error=None,
            n_interface_waters=len(interface_water_res),
        )
    except Exception as exc:  # noqa: BLE001 — log & exclude
        return PrepResult(
            record_id=record_id,
            ok=False,
            complex_path=None,
            receptor_path=None,
            peptide_path=None,
            interface_water_path=None,
            error=str(exc),
        )


def prepare_catalog(
    catalog_path: Path,
    *,
    structure_root: Path,
    out_dir: Path,
    use_pdb2pqr: bool = False,
) -> list[PrepResult]:
    entries = load_peptide_affinity_catalog(catalog_path, require_structure=False)
    results: list[PrepResult] = []
    for entry in entries:
        if not entry.structure_path or not entry.peptide_chain or not entry.receptor_chains:
            results.append(
                PrepResult(
                    record_id=entry.record_id,
                    ok=False,
                    complex_path=None,
                    receptor_path=None,
                    peptide_path=None,
                    interface_water_path=None,
                    error="missing structure_path or chain annotations",
                )
            )
            continue
        pdb_path = structure_root / entry.structure_path
        if not pdb_path.is_file():
            # also try relative to data dir
            alt = DATA_DIR / entry.structure_path
            pdb_path = alt if alt.is_file() else pdb_path
        if not pdb_path.is_file():
            results.append(
                PrepResult(
                    record_id=entry.record_id,
                    ok=False,
                    complex_path=None,
                    receptor_path=None,
                    peptide_path=None,
                    interface_water_path=None,
                    error=f"missing PDB: {entry.structure_path}",
                )
            )
            continue
        rec_chains = tuple(c.strip() for c in entry.receptor_chains.split(",") if c.strip())
        results.append(
            prepare_complex(
                pdb_path,
                out_dir=out_dir,
                record_id=entry.record_id,
                peptide_chain=entry.peptide_chain,
                receptor_chains=rec_chains,
                use_pdb2pqr=use_pdb2pqr,
            )
        )
    return results


def write_prep_manifest(results: list[PrepResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "record_id",
                "ok",
                "complex_path",
                "receptor_path",
                "peptide_path",
                "interface_water_path",
                "n_interface_waters",
                "error",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for r in results:
            writer.writerow(
                {
                    "record_id": r.record_id,
                    "ok": int(r.ok),
                    "complex_path": r.complex_path or "",
                    "receptor_path": r.receptor_path or "",
                    "peptide_path": r.peptide_path or "",
                    "interface_water_path": r.interface_water_path or "",
                    "n_interface_waters": r.n_interface_waters,
                    "error": r.error or "",
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog",
        type=Path,
        default=DATA_DIR / "peptide_affinity_catalog_v2.tsv",
    )
    parser.add_argument("--structure-root", type=Path, default=DATA_DIR)
    parser.add_argument("--out-dir", type=Path, default=DATA_DIR / "prepared")
    parser.add_argument("--pdb2pqr", action="store_true")
    parser.add_argument(
        "--log-json",
        type=Path,
        default=DATA_DIR / "prepared" / "prep_log.json",
    )
    args = parser.parse_args()
    results = prepare_catalog(
        args.catalog,
        structure_root=args.structure_root,
        out_dir=args.out_dir,
        use_pdb2pqr=args.pdb2pqr,
    )
    write_prep_manifest(results, args.out_dir / "prep_manifest.tsv")
    payload = [
        {
            "record_id": r.record_id,
            "ok": r.ok,
            "error": r.error,
            "complex_path": r.complex_path,
            "n_interface_waters": r.n_interface_waters,
        }
        for r in results
    ]
    args.log_json.parent.mkdir(parents=True, exist_ok=True)
    args.log_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    n_ok = sum(1 for r in results if r.ok)
    n_fail = len(results) - n_ok
    print(f"prep done: ok={n_ok} failed={n_fail} (failures logged, not patched)")
    if n_ok == 0:
        raise SystemExit("all structures failed QC — see prep_log.json")


if __name__ == "__main__":
    main()

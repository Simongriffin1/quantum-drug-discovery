#!/usr/bin/env python3
"""Prepare local peptide–target complex PDB fixtures for oracle-validity.

No live network by default. Builds **deterministic** two-chain (receptor + peptide)
complexes from the affinity fixture sequences. Geometry is generated from a seeded
hash of the record id — **not** optimized to maximize Spearman (would be cheating).

Physical note: these are simplified computational models for oracle plumbing and an
honest first Spearman report. Gate status may FAIL; report truthfully.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

# Backbone atoms for a residue in a short α-helix (approx. rise/rotation).
AA_HEAVY = {
    "A": [("N", 0.0, 0.0, 0.0), ("CA", 1.46, 0.0, 0.0), ("C", 2.0, 1.4, 0.0), ("O", 1.4, 2.4, 0.0), ("CB", 1.9, -1.2, 1.0)],
    "G": [("N", 0.0, 0.0, 0.0), ("CA", 1.46, 0.0, 0.0), ("C", 2.0, 1.4, 0.0), ("O", 1.4, 2.4, 0.0)],
}


def _atoms_for_aa(aa: str) -> list[tuple[str, float, float, float]]:
    aa = aa.upper()
    if aa not in AA_HEAVY:
        # Map non-A/G to alanine-like sidechain stub for fixture loading.
        return list(AA_HEAVY["A"])
    return list(AA_HEAVY[aa])


def _helix_transform(i: int, radius: float = 2.3) -> tuple[float, float, float, float]:
    """Return (dx, dy, dz, yaw) for residue i on a crude α-helix."""
    rise = 1.5
    twist = math.radians(100.0)
    ang = i * twist
    return radius * math.cos(ang), radius * math.sin(ang), i * rise, ang


def _seed_offset(record_id: str) -> tuple[float, float, float]:
    digest = hashlib.sha256(record_id.encode()).digest()
    # Offsets in Å — vary binding geometry across the subset without target-hacking.
    ox = (digest[0] / 255.0 - 0.5) * 6.0
    oy = (digest[1] / 255.0 - 0.5) * 6.0
    oz = (digest[2] / 255.0) * 4.0
    return ox, oy, oz


def write_complex_pdb(
    path: Path,
    *,
    receptor_seq: str,
    peptide_seq: str,
    record_id: str,
) -> None:
    ox, oy, oz = _seed_offset(record_id)
    lines: list[str] = [
        f"HEADER    peptideforge_fixture {record_id}",
        "REMARK    deterministic helix+peptide model for OpenMM MM-GBSA plumbing",
        "REMARK    NOT an experimental structure; experimental affinity is separate",
    ]
    atom_i = 1

    def emit(chain: str, res_i: int, resname: str, name: str, x: float, y: float, z: float) -> None:
        nonlocal atom_i
        # PDB atom name field is columns 13-16
        nm = f"{name:>4}" if len(name) < 4 else name[:4]
        lines.append(
            f"ATOM  {atom_i:5d} {nm} {resname:>3} {chain}{res_i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {name[0]:>2}"
        )
        atom_i += 1

    # Receptor chain A — poly-ish helix from receptor_seq (or polyA fallback)
    rec = (receptor_seq or "AAAAAAAAAAAA").upper()
    rec = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "A" for c in rec)[:16]
    if len(rec) < 8:
        rec = (rec + "AAAAAAAA")[:12]

    for i, aa in enumerate(rec):
        dx, dy, dz, yaw = _helix_transform(i)
        for name, lx, ly, lz in _atoms_for_aa(aa):
            # local rotate about z by yaw (approx)
            rx = lx * math.cos(yaw) - ly * math.sin(yaw)
            ry = lx * math.sin(yaw) + ly * math.cos(yaw)
            emit("A", i + 1, "ALA" if aa != "G" else "GLY", name, rx + dx, ry + dy, lz + dz)

    # Peptide chain P — short helix translated by seeded offset toward receptor
    pep = peptide_seq.upper()
    pep = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "A" for c in pep)
    if len(pep) < 5:
        pep = (pep + "AAAAA")[:5]
    pep = pep[:12]

    for i, aa in enumerate(pep):
        dx, dy, dz, yaw = _helix_transform(i, radius=2.0)
        for name, lx, ly, lz in _atoms_for_aa(aa):
            rx = lx * math.cos(yaw) - ly * math.sin(yaw)
            ry = lx * math.sin(yaw) + ly * math.cos(yaw)
            emit(
                "P",
                i + 1,
                "ALA" if aa != "G" else "GLY",
                name,
                rx + dx + 8.0 + ox,
                ry + dy + oy,
                lz + dz + oz,
            )

    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_affinity_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--affinity-tsv",
        type=Path,
        default=Path(__file__).resolve().parents[2]
        / "fixtures"
        / "pdbbind_peptide_affinity_v1.tsv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "fixtures" / "structures",
    )
    parser.add_argument("--limit", type=int, default=6, help="How many records to materialize")
    args = parser.parse_args()

    rows = [r for r in load_affinity_tsv(args.affinity_tsv) if r.get("record_id")]
    # One per cluster for diversity
    seen_cluster: set[str] = set()
    selected: list[dict[str, str]] = []
    for row in rows:
        c = row.get("cluster_id") or row["record_id"]
        if c in seen_cluster:
            continue
        seen_cluster.add(c)
        selected.append(row)
        if len(selected) >= args.limit:
            break

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out_dir / "structure_manifest_v1.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "record_id",
                "pdb_id",
                "pdb_file",
                "peptide_chain",
                "peptide_sequence",
                "experimental_pk",
                "source_note",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in selected:
            pdb_name = f"{row['record_id']}_complex.pdb"
            out_pdb = args.out_dir / pdb_name
            write_complex_pdb(
                out_pdb,
                receptor_seq=row.get("receptor_sequence") or "AAAKAAAAKAAAA",
                peptide_seq=row["peptide_sequence"],
                record_id=row["record_id"],
            )
            writer.writerow(
                {
                    "record_id": row["record_id"],
                    "pdb_id": row["pdb_id"],
                    "pdb_file": pdb_name,
                    "peptide_chain": "P",
                    "peptide_sequence": row["peptide_sequence"],
                    "experimental_pk": row["pk"],
                    "source_note": (
                        "deterministic_model_geometry; experimental_pk from "
                        "pdbbind_peptide_affinity_v1 fixture; not an experimental pose"
                    ),
                }
            )
    print(f"wrote {len(selected)} complexes → {args.out_dir}")
    print(f"manifest → {manifest_path}")


if __name__ == "__main__":
    main()

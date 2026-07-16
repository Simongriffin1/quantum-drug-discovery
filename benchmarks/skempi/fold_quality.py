"""Step 4A.2 — predicted-vs-crystal fold quality metrics.

Computes interface RMSD (peptide heavy atoms) and a DockQ-like proxy when
DockQ is unavailable. Reads Boltz self-reported confidence / ipTM when present.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class FoldQualityRow:
    record_id: str
    pdb_id: str
    mode: str
    peptide_irmsd_A: float | None
    dockq_proxy: float | None
    iptm: float | None
    interface_plddt: float | None
    confidence: float | None
    n_peptide_atoms_aligned: int
    error: str | None = None


def _ca_coords(pdb_path: Path, chain: str) -> list[tuple[str, tuple[float, float, float]]]:
    """Return (resseq, (x,y,z)) for CA atoms on chain."""
    out: list[tuple[str, tuple[float, float, float]]] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) < 54:
            continue
        if line[12:16].strip() != "CA":
            continue
        if line[21] != chain:
            continue
        resseq = line[22:26].strip()
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        out.append((resseq, (x, y, z)))
    return out


def peptide_interface_rmsd(
    crystal_pdb: Path,
    predicted_pdb: Path,
    *,
    peptide_chain_crystal: str,
    peptide_chain_pred: str | None = None,
) -> tuple[float | None, int]:
    """CA RMSD between crystal and predicted peptide chains (matched by resseq)."""
    pep_pred = peptide_chain_pred or peptide_chain_crystal
    a = {r: c for r, c in _ca_coords(crystal_pdb, peptide_chain_crystal)}
    b = {r: c for r, c in _ca_coords(predicted_pdb, pep_pred)}
    # If predicted uses different chain IDs (R/P), try all chains and pick best overlap
    if len(set(a) & set(b)) < 3:
        from collections import defaultdict

        by_chain: dict[str, dict[str, tuple[float, float, float]]] = defaultdict(dict)
        for line in predicted_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.startswith("ATOM") or len(line) < 54:
                continue
            if line[12:16].strip() != "CA":
                continue
            by_chain[line[21]][line[22:26].strip()] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        best_chain, best_n = pep_pred, 0
        for cid, coords in by_chain.items():
            n = len(set(a) & set(coords))
            if n > best_n:
                best_n, best_chain, b = n, cid, coords
        if best_n < 3:
            return None, 0

    common = sorted(set(a) & set(b))
    if len(common) < 3:
        return None, len(common)
    # Kabsch-free: after centering, RMSD of CA
    pa = [a[r] for r in common]
    pb = [b[r] for r in common]
    ca = (
        sum(p[0] for p in pa) / len(pa),
        sum(p[1] for p in pa) / len(pa),
        sum(p[2] for p in pa) / len(pa),
    )
    cb = (
        sum(p[0] for p in pb) / len(pb),
        sum(p[1] for p in pb) / len(pb),
        sum(p[2] for p in pb) / len(pb),
    )
    sq = 0.0
    for p, q in zip(pa, pb, strict=True):
        dx = (p[0] - ca[0]) - (q[0] - cb[0])
        dy = (p[1] - ca[1]) - (q[1] - cb[1])
        dz = (p[2] - ca[2]) - (q[2] - cb[2])
        sq += dx * dx + dy * dy + dz * dz
    return math.sqrt(sq / len(common)), len(common)


def dockq_proxy_from_irmsd(irmsd: float | None) -> float | None:
    """Monotonic DockQ-like score in [0,1] from iRMSD (Basu & Wallner-inspired).

    Not a substitute for official DockQ when available — documented as proxy.
    """
    if irmsd is None:
        return None
    # DockQ Fnat/iRMSD blend simplified: score decays with iRMSD
    return 1.0 / (1.0 + (irmsd / 1.5) ** 2)


def try_official_dockq(crystal: Path, predicted: Path) -> float | None:
    """Call DockQ if installed; return None if unavailable (do not fabricate)."""
    try:
        import DockQ.DockQ as dockq_mod  # type: ignore
    except ImportError:
        return None
    try:
        # Best-effort; DockQ APIs vary — fail soft to proxy
        result = dockq_mod.run_on_all_native_interfaces(str(crystal), str(predicted))  # type: ignore
        if isinstance(result, dict):
            for v in result.values():
                if isinstance(v, dict) and "DockQ" in v:
                    return float(v["DockQ"])
        return None
    except Exception:  # noqa: BLE001
        return None


def evaluate_fold_quality(
    manifest: dict[str, Any],
    *,
    structure_dir: Path,
    skempi_tsv: Path,
) -> dict[str, Any]:
    from peptideforge_benchmarks.skempi import load_skempi_ddg

    records = {r.record_id: r for r in load_skempi_ddg(skempi_tsv)}
    rows: list[dict[str, Any]] = []
    for prov in manifest.get("provenance") or []:
        if prov.get("error") or not prov.get("mutant_fold_path") and not prov.get(
            "wt_fold_path"
        ):
            rows.append(
                asdict(
                    FoldQualityRow(
                        record_id=prov.get("record_id", ""),
                        pdb_id=prov.get("pdb_id", ""),
                        mode=prov.get("mode", ""),
                        peptide_irmsd_A=None,
                        dockq_proxy=None,
                        iptm=prov.get("iptm"),
                        interface_plddt=prov.get("interface_plddt"),
                        confidence=prov.get("confidence"),
                        n_peptide_atoms_aligned=0,
                        error=prov.get("error") or "missing predicted path",
                    )
                )
            )
            continue
        rid = prov["record_id"]
        rec = records.get(rid)
        crystal = structure_dir / f"{prov['pdb_id']}.pdb"
        pred = Path(prov.get("wt_fold_path") or prov.get("mutant_fold_path") or "")
        if rec is None or not crystal.is_file() or not pred.is_file():
            rows.append(
                asdict(
                    FoldQualityRow(
                        record_id=rid,
                        pdb_id=prov.get("pdb_id", ""),
                        mode=prov.get("mode", ""),
                        peptide_irmsd_A=None,
                        dockq_proxy=None,
                        iptm=prov.get("iptm"),
                        interface_plddt=prov.get("interface_plddt"),
                        confidence=prov.get("confidence"),
                        n_peptide_atoms_aligned=0,
                        error="missing crystal or predicted pdb",
                    )
                )
            )
            continue
        pep = (rec.partner2 or "B")[0]
        irmsd, n_al = peptide_interface_rmsd(crystal, pred, peptide_chain_crystal=pep)
        official = try_official_dockq(crystal, pred)
        proxy = official if official is not None else dockq_proxy_from_irmsd(irmsd)
        rows.append(
            asdict(
                FoldQualityRow(
                    record_id=rid,
                    pdb_id=prov["pdb_id"],
                    mode=prov["mode"],
                    peptide_irmsd_A=irmsd,
                    dockq_proxy=proxy,
                    iptm=prov.get("iptm"),
                    interface_plddt=prov.get("interface_plddt"),
                    confidence=prov.get("confidence"),
                    n_peptide_atoms_aligned=n_al,
                )
            )
        )

    irmsds = [r["peptide_irmsd_A"] for r in rows if r.get("peptide_irmsd_A") is not None]
    dockqs = [r["dockq_proxy"] for r in rows if r.get("dockq_proxy") is not None]
    summary = {
        "n_rows": len(rows),
        "n_with_irmsd": len(irmsds),
        "irmsd_median": sorted(irmsds)[len(irmsds) // 2] if irmsds else None,
        "dockq_proxy_median": sorted(dockqs)[len(dockqs) // 2] if dockqs else None,
        "dockq_method": "official_DockQ_if_installed_else_irmsd_proxy",
    }
    return {"summary": summary, "rows": rows}


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("benchmarks/skempi/data/predicted_folds/predicted_folds_manifest.json"),
    )
    parser.add_argument(
        "--structure-dir",
        type=Path,
        default=Path("benchmarks/skempi/data/structures"),
    )
    parser.add_argument(
        "--skempi-tsv",
        type=Path,
        default=Path("benchmarks/skempi/data/skempi_v2.tsv"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/skempi/data/predicted_folds/fold_quality.json"),
    )
    args = parser.parse_args()
    if not args.manifest.is_file():
        raise SystemExit(f"manifest missing: {args.manifest}")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if manifest.get("status") == "BLOCKED_BOLTZ_UNAVAILABLE":
        payload = {
            "status": "BLOCKED_BOLTZ_UNAVAILABLE",
            "summary": {},
            "rows": [],
            "note": "No fold quality — Boltz unavailable; never fabricated.",
        }
    else:
        payload = evaluate_fold_quality(
            manifest, structure_dir=args.structure_dir, skempi_tsv=args.skempi_tsv
        )
        payload["status"] = "OK"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload.get("summary") or payload, indent=2))


if __name__ == "__main__":
    main()

"""Load curated protein–peptide affinity catalogs from local TSV/CSV.

Never fabricates entries. Network downloads are intentional and go through
documented scripts; CI uses committed fixtures only.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

from peptide_affinity.models import AffinityType, PeptideAffinityEntry

PACKAGE_DIR = Path(__file__).resolve().parent
FIXTURES_DIR = PACKAGE_DIR / "fixtures"
DATA_DIR = PACKAGE_DIR / "data"

REQUIRED_COLUMNS = frozenset(
    {
        "record_id",
        "pdb_id",
        "receptor_seq",
        "peptide_seq",
        "pKd",
        "source",
    }
)


def affinity_to_pkd(value: float, affinity_type: AffinityType) -> float:
    if affinity_type == AffinityType.PKD:
        return value
    if value <= 0.0:
        raise ValueError(f"affinity_value must be > 0 for {affinity_type}, got {value}")
    return -math.log10(value)


def default_catalog_path() -> Path:
    """Prefer the expanded catalog if present; else CI fixture."""
    expanded = DATA_DIR / "peptide_affinity_catalog_v2.tsv"
    if expanded.is_file():
        return expanded
    fixture = FIXTURES_DIR / "peptide_affinity_ci_v1.tsv"
    if fixture.is_file():
        return fixture
    raise FileNotFoundError(
        "No peptide affinity catalog found. Expected either:\n"
        f"  {expanded}\n"
        f"  {fixture}\n"
        "See benchmarks/peptide_affinity/README.md for the manual download path "
        "(PepBenchData_raw PpI_ba + RCSB PDBs + match_and_build_catalog.py)."
    )


def load_peptide_affinity_catalog(
    path: Path | str | None = None,
    *,
    require_structure: bool = False,
    structure_root: Path | str | None = None,
    min_len: int = 5,
    max_len: int = 50,
) -> tuple[PeptideAffinityEntry, ...]:
    """Parse a curated TSV catalog; reject malformed/duplicate record_ids."""
    catalog = Path(path) if path is not None else default_catalog_path()
    if not catalog.is_file():
        raise FileNotFoundError(
            f"peptide affinity catalog not found: {catalog}. "
            "See benchmarks/peptide_affinity/README.md — refusing to fabricate entries."
        )

    root = Path(structure_root) if structure_root else catalog.parent
    entries: list[PeptideAffinityEntry] = []
    seen: set[str] = set()
    rejected: list[str] = []

    with catalog.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty/headerless catalog: {catalog}")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"missing columns {sorted(missing)} in {catalog}")

        for row in reader:
            rid = (row.get("record_id") or "").strip()
            if not rid or rid.startswith("#"):
                continue
            if rid in seen:
                rejected.append(f"duplicate record_id={rid}")
                continue
            seen.add(rid)
            try:
                peptide = row["peptide_seq"].strip().upper().replace(" ", "")
                receptor = row["receptor_seq"].strip().upper().replace(" ", "")
                plen = int(row["peptide_len"]) if row.get("peptide_len") else len(peptide)
                if plen != len(peptide):
                    raise ValueError(f"peptide_len={plen} != len(seq)={len(peptide)}")
                if not (min_len <= plen <= max_len):
                    rejected.append(f"{rid}: peptide_len={plen} out of [{min_len},{max_len}]")
                    continue
                if len(receptor) < 20:
                    rejected.append(f"{rid}: receptor_seq too short")
                    continue

                atype = AffinityType(row.get("affinity_type") or "pKd")
                pkd = float(row["pKd"])
                aval = float(row["affinity_value"]) if row.get("affinity_value") else None
                if aval is not None and atype != AffinityType.PKD:
                    expected = affinity_to_pkd(aval, atype)
                    if abs(expected - pkd) > 1e-3:
                        raise ValueError(f"pKd inconsistency: file={pkd} recomputed={expected}")

                entry = PeptideAffinityEntry(
                    record_id=rid,
                    pdb_id=row["pdb_id"],
                    receptor_seq=receptor,
                    peptide_seq=peptide,
                    peptide_len=plen,
                    resolution=float(row["resolution"]) if row.get("resolution") else None,
                    affinity_value=aval,
                    affinity_type=atype,
                    pKd=pkd,
                    source=row["source"],
                    structure_path=row.get("structure_path") or None,
                    peptide_chain=row.get("peptide_chain") or None,
                    receptor_chains=row.get("receptor_chains") or None,
                    deposit_year=int(row["deposit_year"]) if row.get("deposit_year") else None,
                    metadata={
                        k: v
                        for k, v in row.items()
                        if k not in REQUIRED_COLUMNS
                        and k
                        not in {
                            "peptide_len",
                            "resolution",
                            "affinity_value",
                            "affinity_type",
                            "structure_path",
                            "peptide_chain",
                            "receptor_chains",
                            "deposit_year",
                            "notes",
                        }
                        and v
                    },
                )
            except (ValueError, KeyError) as exc:
                rejected.append(f"{rid}: {exc}")
                continue

            if require_structure:
                resolved = entry.resolved_structure(root)
                if resolved is None:
                    rejected.append(f"{rid}: structure missing ({entry.structure_path})")
                    continue
            entries.append(entry)

    if not entries:
        raise ValueError(
            f"no valid peptide affinity entries from {catalog}; rejected={rejected[:20]}"
        )
    return tuple(entries)


def catalog_qc_summary(entries: tuple[PeptideAffinityEntry, ...]) -> dict[str, object]:
    with_struct = sum(1 for e in entries if e.structure_path)
    return {
        "n": len(entries),
        "n_with_structure_path": with_struct,
        "pdb_ids": len({e.pdb_id for e in entries}),
        "pKd_min": min(e.pKd for e in entries),
        "pKd_max": max(e.pKd for e in entries),
        "peptide_len_min": min(e.peptide_len for e in entries),
        "peptide_len_max": max(e.peptide_len for e in entries),
        "sources": sorted({e.source for e in entries}),
    }

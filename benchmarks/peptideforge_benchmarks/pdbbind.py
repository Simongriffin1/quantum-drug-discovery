"""PDBbind-derived protein–peptide affinity loader.

Parses a TSV fixture with columns matching a documented subset schema.
Does NOT download PDBbind — redistribution of full PDBbind may be
restricted; see ``fixtures/LICENSE.md``.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

from peptideforge_benchmarks.models import AffinityRecord, AffinityUnit
from peptideforge_benchmarks.paths import require_fixture

REQUIRED_COLUMNS = frozenset(
    {
        "record_id",
        "pdb_id",
        "peptide_sequence",
        "affinity_value",
        "affinity_unit",
        "pk",
    }
)


def affinity_to_pk(value: float, unit: AffinityUnit) -> float:
    """Convert experimental affinity to pK = -log10(M).

    Physical rationale: ranking oracles against pK keeps Spearman comparable
    across Kd / Ki / IC50 reports. IC50 is treated as an approximate proxy
    for affinity (known limitation — document when correlating).
    """
    if unit == AffinityUnit.PK:
        return value
    if value <= 0.0:
        raise ValueError(f"affinity_value must be > 0 for unit={unit}, got {value}")
    # Fixtures store molar (or already converted); if value looks like nM/uM
    # the fixture MUST already convert to M before writing affinity_value.
    return -math.log10(value)


def load_pdbbind_peptide_affinity(
    path: Path | str | None = None,
    *,
    require_clusters: bool = False,
) -> tuple[AffinityRecord, ...]:
    """Load protein–peptide affinity records from a local TSV fixture.

    Default path: ``benchmarks/fixtures/pdbbind_peptide_affinity_v1.tsv``.
    """
    fixture = Path(path) if path is not None else require_fixture("pdbbind_peptide_affinity_v1.tsv")
    if not fixture.is_file():
        raise FileNotFoundError(f"PDBbind peptide affinity fixture not found: {fixture}")

    with fixture.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty or headerless TSV: {fixture}")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"missing columns {sorted(missing)} in {fixture}")

        records: list[AffinityRecord] = []
        for row in reader:
            if not row.get("record_id") or row["record_id"].startswith("#"):
                continue
            unit = AffinityUnit(row["affinity_unit"])
            pk = float(row["pk"])
            # Consistency check — fail loud if fixture pk disagrees with value
            expected_pk = affinity_to_pk(float(row["affinity_value"]), unit)
            if abs(expected_pk - pk) > 1e-3:
                raise ValueError(
                    f"pk inconsistency for {row['record_id']}: "
                    f"file pk={pk}, recomputed={expected_pk}"
                )
            cluster = row.get("cluster_id") or None
            if require_clusters and not cluster:
                raise ValueError(f"cluster_id required but missing for {row['record_id']}")
            records.append(
                AffinityRecord(
                    record_id=row["record_id"],
                    pdb_id=row["pdb_id"],
                    peptide_sequence=row["peptide_sequence"],
                    receptor_sequence=row.get("receptor_sequence") or None,
                    cluster_id=cluster,
                    affinity_value=float(row["affinity_value"]),
                    affinity_unit=unit,
                    pk=pk,
                    ligand_name=row.get("ligand_name") or None,
                    resolution_A=(float(row["resolution_A"]) if row.get("resolution_A") else None),
                    source=row.get("source") or "pdbbind_peptide_subset",
                    metadata={
                        k: v
                        for k, v in row.items()
                        if k not in REQUIRED_COLUMNS
                        and k
                        not in {
                            "receptor_sequence",
                            "cluster_id",
                            "ligand_name",
                            "resolution_A",
                            "source",
                        }
                        and v
                    },
                )
            )

    if not records:
        raise ValueError(f"no affinity records parsed from {fixture}")
    return tuple(records)

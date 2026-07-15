"""SKEMPI mutation ΔΔG loader.

Parses a TSV fixture with a SKEMPI-like column subset. Does not download
SKEMPI over the network — ship fixtures under ``benchmarks/fixtures/``.
"""

from __future__ import annotations

import csv
from pathlib import Path

from peptideforge_benchmarks.models import MutationRecord
from peptideforge_benchmarks.paths import require_fixture

REQUIRED_COLUMNS = frozenset(
    {
        "record_id",
        "pdb_id",
        "mutant",
        "ddg_kcal_mol",
    }
)


def load_skempi_ddg(
    path: Path | str | None = None,
    *,
    require_clusters: bool = False,
) -> tuple[MutationRecord, ...]:
    """Load mutation ΔΔG records from a local TSV fixture.

    Default path: ``benchmarks/fixtures/skempi_ddg_v1.tsv``.
    """
    fixture = Path(path) if path is not None else require_fixture("skempi_ddg_v1.tsv")
    if not fixture.is_file():
        raise FileNotFoundError(f"SKEMPI ddG fixture not found: {fixture}")

    with fixture.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty or headerless TSV: {fixture}")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"missing columns {sorted(missing)} in {fixture}")

        records: list[MutationRecord] = []
        for row in reader:
            if not row.get("record_id") or row["record_id"].startswith("#"):
                continue
            cluster = row.get("cluster_id") or None
            if require_clusters and not cluster:
                raise ValueError(f"cluster_id required but missing for {row['record_id']}")
            records.append(
                MutationRecord(
                    record_id=row["record_id"],
                    pdb_id=row["pdb_id"],
                    mutant=row["mutant"],
                    partner1=row.get("partner1") or None,
                    partner2=row.get("partner2") or None,
                    ddg_kcal_mol=float(row["ddg_kcal_mol"]),
                    temperature_K=(
                        float(row["temperature_K"]) if row.get("temperature_K") else None
                    ),
                    cluster_id=cluster,
                    wildtype_sequence=row.get("wildtype_sequence") or None,
                    mutant_sequence=row.get("mutant_sequence") or None,
                    source=row.get("source") or "skempi",
                    metadata={
                        k: v
                        for k, v in row.items()
                        if k not in REQUIRED_COLUMNS
                        and k
                        not in {
                            "partner1",
                            "partner2",
                            "temperature_K",
                            "cluster_id",
                            "wildtype_sequence",
                            "mutant_sequence",
                            "source",
                        }
                        and v
                    },
                )
            )

    if not records:
        raise ValueError(f"no mutation records parsed from {fixture}")
    return tuple(records)

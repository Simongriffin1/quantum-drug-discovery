"""Public benchmark loaders for PeptideForge (PDBbind-derived, SKEMPI).

Fixtures live under ``benchmarks/fixtures/``. No live network in CI.
License notes: see fixtures/LICENSE.md — check before redistributing.
"""

from peptideforge_benchmarks.models import (
    AffinityRecord,
    BenchmarkSplit,
    MutationRecord,
)
from peptideforge_benchmarks.pdbbind import load_pdbbind_peptide_affinity
from peptideforge_benchmarks.skempi import load_skempi_ddg
from peptideforge_benchmarks.splits import (
    homology_aware_split,
    load_cluster_assignments,
    sequence_identity,
)

__all__ = [
    "AffinityRecord",
    "BenchmarkSplit",
    "MutationRecord",
    "homology_aware_split",
    "load_cluster_assignments",
    "load_pdbbind_peptide_affinity",
    "load_skempi_ddg",
    "sequence_identity",
]

__version__ = "0.1.0"

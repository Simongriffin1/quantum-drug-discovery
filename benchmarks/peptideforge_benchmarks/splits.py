"""Homology-aware train/test splits for surrogate and oracle benchmarks.

Prefer precomputed cluster assignments (MMseqs2 / CD-HIT offline). For small
fixtures, a pure-Python identity greedy clustering is available. External
cluster tools are optional and fail loud when required but missing.
"""

from __future__ import annotations

import csv
import random
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from peptideforge_benchmarks.models import BenchmarkSplit
from peptideforge_benchmarks.paths import require_fixture


def sequence_identity(a: str, b: str) -> float:
    """Global length-normalized identity for equal-length sequences.

    For unequal lengths, uses the shorter length in the denominator after a
    simple contiguous alignment approximation (min-length window on the longer
    sequence via best identity) — sufficient for peptide-scale clustering in
    fixtures. Raises if either sequence is empty.
    """
    sa, sb = a.upper(), b.upper()
    if not sa or not sb:
        raise ValueError("sequences must be non-empty")
    if len(sa) == len(sb):
        matches = sum(x == y for x, y in zip(sa, sb, strict=True))
        return matches / len(sa)
    # Slide shorter over longer; take best identity over shorter length.
    short, long = (sa, sb) if len(sa) <= len(sb) else (sb, sa)
    best = 0
    for start in range(len(long) - len(short) + 1):
        window = long[start : start + len(short)]
        matches = sum(x == y for x, y in zip(short, window, strict=True))
        best = max(best, matches)
    return best / len(short)


def load_cluster_assignments(path: Path | str | None = None) -> dict[str, str]:
    """Load ``record_id → cluster_id`` from a TSV (columns: record_id, cluster_id)."""
    fixture = Path(path) if path is not None else require_fixture("clusters_v1.tsv")
    if not fixture.is_file():
        raise FileNotFoundError(f"cluster assignment file not found: {fixture}")
    assignments: dict[str, str] = {}
    with fixture.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"record_id", "cluster_id"} <= set(reader.fieldnames):
            raise ValueError(f"clusters TSV needs record_id, cluster_id columns: {fixture}")
        for row in reader:
            rid = row["record_id"]
            if not rid or rid.startswith("#"):
                continue
            assignments[rid] = row["cluster_id"]
    if not assignments:
        raise ValueError(f"no cluster assignments in {fixture}")
    return assignments


def greedy_identity_clusters(
    sequences: dict[str, str],
    *,
    identity_threshold: float = 0.4,
) -> dict[str, str]:
    """Greedy clustering: assign each sequence to first cluster with identity ≥ threshold.

    Used when MMseqs2/CD-HIT are unavailable. Deterministic given insertion order
    of ``sequences``.
    """
    if not 0.0 < identity_threshold <= 1.0:
        raise ValueError("identity_threshold must be in (0, 1]")
    representatives: list[tuple[str, str]] = []  # (cluster_id, sequence)
    assignments: dict[str, str] = {}
    for record_id, seq in sequences.items():
        placed = False
        for cluster_id, rep_seq in representatives:
            if sequence_identity(seq, rep_seq) >= identity_threshold:
                assignments[record_id] = cluster_id
                placed = True
                break
        if not placed:
            cluster_id = f"c{len(representatives)}"
            representatives.append((cluster_id, seq))
            assignments[record_id] = cluster_id
    return assignments


def homology_aware_split(
    record_ids: tuple[str, ...] | list[str],
    cluster_of: dict[str, str],
    *,
    test_fraction: float = 0.25,
    val_fraction: float = 0.0,
    seed: int = 0,
    split_name: str = "homology_v1",
) -> BenchmarkSplit:
    """Split by holding out entire homology clusters (no cluster shared across splits).

    Physical/eval rationale: random splits leak homologous sequences and inflate
    surrogate Spearman. Cluster hold-out is the minimum credible protocol.
    """
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("test_fraction must be in (0, 1)")
    if val_fraction < 0.0 or test_fraction + val_fraction >= 1.0:
        raise ValueError("val_fraction invalid relative to test_fraction")

    missing = [rid for rid in record_ids if rid not in cluster_of]
    if missing:
        raise KeyError(f"missing cluster assignments for records: {missing[:5]}")

    clusters: dict[str, list[str]] = defaultdict(list)
    for rid in record_ids:
        clusters[cluster_of[rid]].append(rid)

    cluster_ids = sorted(clusters.keys())
    rng = random.Random(seed)
    rng.shuffle(cluster_ids)

    n_test = max(1, int(round(len(cluster_ids) * test_fraction)))
    n_val = int(round(len(cluster_ids) * val_fraction)) if val_fraction > 0 else 0
    if n_test + n_val >= len(cluster_ids):
        raise ValueError("not enough clusters for requested fractions")

    test_clusters = set(cluster_ids[:n_test])
    val_clusters = set(cluster_ids[n_test : n_test + n_val])
    train_clusters = set(cluster_ids[n_test + n_val :])

    def ids_for(cset: set[str]) -> tuple[str, ...]:
        out: list[str] = []
        for cid in sorted(cset):
            out.extend(sorted(clusters[cid]))
        return tuple(out)

    return BenchmarkSplit(
        split_name=split_name,
        train_ids=ids_for(train_clusters),
        test_ids=ids_for(test_clusters),
        val_ids=ids_for(val_clusters),
        method="cluster_holdout",
        seed=seed,
        max_train_test_identity=None,
        notes=f"clusters_train={len(train_clusters)} test={len(test_clusters)} val={len(val_clusters)}",
    )


def run_mmseqs_easy_cluster(
    fasta_path: Path,
    output_prefix: Path,
    *,
    min_seq_id: float = 0.4,
) -> Path:
    """Optional MMseqs2 clustering — raises if ``mmseqs`` is not on PATH."""
    mmseqs = shutil.which("mmseqs")
    if mmseqs is None:
        raise RuntimeError(
            "mmseqs not found on PATH. Install MMseqs2 or use precomputed "
            "cluster fixtures / greedy_identity_clusters (fail loud, no silent fallback)."
        )
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        mmseqs,
        "easy-cluster",
        str(fasta_path),
        str(output_prefix),
        str(output_prefix.parent / "tmp"),
        "--min-seq-id",
        str(min_seq_id),
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"mmseqs failed ({result.returncode}): {result.stderr}")
    tsv = Path(str(output_prefix) + "_cluster.tsv")
    if not tsv.is_file():
        raise FileNotFoundError(f"mmseqs did not produce {tsv}")
    return tsv

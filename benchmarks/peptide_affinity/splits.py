"""Homology-aware, leakage-controlled splits for the peptide-affinity catalog.

Cold-start: cluster receptors (and peptides) at 30% identity — no cluster shared
across train/dev/test. Time-based split: deposit_year ≤ 2019 vs > 2019 as an
independent leakage check (LP-PDBBind style).
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from peptideforge_benchmarks.splits import (
    greedy_identity_clusters,
    homology_aware_split,
    sequence_identity,
)

from peptide_affinity.load import DATA_DIR, load_peptide_affinity_catalog
from peptide_affinity.models import PeptideAffinityEntry


def cluster_entries(
    entries: tuple[PeptideAffinityEntry, ...],
    *,
    identity_threshold: float = 0.30,
) -> dict[str, str]:
    """Union-find clusters on receptor OR peptide identity (≥ threshold).

    Receptor-only clustering previously allowed peptide-sequence leakage across
    the train/test cut (red-team stop-the-line). Two records share a cluster if
    receptors OR peptides are ≥ ``identity_threshold`` identical (transitive).
    """
    ids = [e.record_id for e in entries]
    parent = {i: i for i in ids}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    by_id = {e.record_id: e for e in entries}
    # Pairwise identity edges (N is small: O(N²) peptides/receptors is fine for ~100)
    for i, a in enumerate(ids):
        ea = by_id[a]
        for b in ids[i + 1 :]:
            eb = by_id[b]
            if sequence_identity(ea.receptor_seq, eb.receptor_seq) >= identity_threshold:
                union(a, b)
            elif sequence_identity(ea.peptide_seq, eb.peptide_seq) >= identity_threshold:
                union(a, b)

    # Stable cluster labels
    roots = {find(i) for i in ids}
    root_to_cid = {r: f"j{k}" for k, r in enumerate(sorted(roots))}
    return {i: root_to_cid[find(i)] for i in ids}


def assert_no_sequence_leakage(
    train_ids: tuple[str, ...],
    test_ids: tuple[str, ...],
    entries: tuple[PeptideAffinityEntry, ...],
    *,
    val_ids: tuple[str, ...] = (),
    max_identity: float = 0.30,
) -> None:
    """Fail loud if any train/val peptide is ≥ max_identity to a test peptide."""
    by_id = {e.record_id: e for e in entries}
    train_pool = tuple(train_ids) + tuple(val_ids)
    train_seqs = {i: by_id[i].peptide_seq for i in train_pool if i in by_id}
    test_seqs = {i: by_id[i].peptide_seq for i in test_ids if i in by_id}
    # identical
    shared = set(train_seqs.values()) & set(test_seqs.values())
    if shared:
        raise AssertionError(
            f"identical peptide sequence leakage across split ({len(shared)} sequences)"
        )
    close: list[str] = []
    for tid, tseq in test_seqs.items():
        for sid, sseq in train_seqs.items():
            if sequence_identity(tseq, sseq) > max_identity:
                close.append(f"{tid}~{sid}")
                if len(close) >= 5:
                    break
        if len(close) >= 5:
            break
    if close:
        raise AssertionError(
            f"peptide identity >{max_identity:.0%} leakage: {close}"
        )


def peptide_clusters(
    entries: tuple[PeptideAffinityEntry, ...],
    *,
    identity_threshold: float = 0.30,
) -> dict[str, str]:
    return greedy_identity_clusters(
        {e.record_id: e.peptide_seq for e in entries},
        identity_threshold=identity_threshold,
    )


def assert_zero_cluster_overlap(
    train_ids: tuple[str, ...],
    test_ids: tuple[str, ...],
    cluster_of: dict[str, str],
    *,
    val_ids: tuple[str, ...] = (),
) -> None:
    parts = [("train", train_ids), ("test", test_ids)]
    if val_ids:
        parts.append(("val", val_ids))
    cluster_sets = {
        name: {cluster_of[i] for i in ids if i in cluster_of} for name, ids in parts
    }
    names = list(cluster_sets)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            overlap = cluster_sets[a] & cluster_sets[b]
            if overlap:
                raise AssertionError(f"cluster overlap {a}∩{b}: {sorted(overlap)[:10]}")


def cold_start_split(
    entries: tuple[PeptideAffinityEntry, ...],
    *,
    test_fraction: float = 0.25,
    val_fraction: float = 0.15,
    seed: int = 0,
    identity_threshold: float = 0.30,
    min_test_clusters: int = 10,
) -> tuple[object, dict[str, str]]:
    cluster_of = cluster_entries(entries, identity_threshold=identity_threshold)
    ids = tuple(e.record_id for e in entries)
    split = homology_aware_split(
        ids,
        cluster_of,
        test_fraction=test_fraction,
        val_fraction=val_fraction,
        seed=seed,
        split_name="cold_start_30pct_receptor",
    )
    n_test_clusters = len({cluster_of[i] for i in split.test_ids})
    if n_test_clusters < min_test_clusters:
        raise ValueError(
            f"test clusters={n_test_clusters} < min_test_clusters={min_test_clusters}; "
            "expand the catalog before claiming a measurable gate."
        )
    assert_zero_cluster_overlap(
        split.train_ids, split.test_ids, cluster_of, val_ids=split.val_ids
    )
    assert_no_sequence_leakage(
        split.train_ids,
        split.test_ids,
        entries,
        val_ids=split.val_ids,
        max_identity=identity_threshold,
    )
    return split, cluster_of


def time_based_split(
    entries: tuple[PeptideAffinityEntry, ...],
    *,
    cutoff_year: int = 2019,
) -> dict[str, tuple[str, ...]]:
    """Structures deposited ≤ cutoff → train/dev pool; > cutoff → test."""
    early: list[str] = []
    late: list[str] = []
    unknown: list[str] = []
    for e in entries:
        if e.deposit_year is None:
            unknown.append(e.record_id)
        elif e.deposit_year <= cutoff_year:
            early.append(e.record_id)
        else:
            late.append(e.record_id)
    if not late:
        raise ValueError(
            f"time-based split has empty test (no deposit_year > {cutoff_year}). "
            "Populate deposit_year or expand the catalog."
        )
    return {
        "train_dev": tuple(sorted(early)),
        "test": tuple(sorted(late)),
        "unknown_year": tuple(sorted(unknown)),
    }


def max_cross_identity(
    train_seqs: dict[str, str],
    test_seqs: dict[str, str],
) -> float:
    best = 0.0
    for tid, tseq in test_seqs.items():
        for sid, sseq in train_seqs.items():
            best = max(best, sequence_identity(tseq, sseq))
            _ = tid, sid
    return best


def write_split_artifact(
    path: Path,
    *,
    cold_split: object,
    cluster_of: dict[str, str],
    time_split: dict[str, tuple[str, ...]] | None,
    counts: dict[str, int],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "cold_start": {
            "split_name": cold_split.split_name,  # type: ignore[attr-defined]
            "train_ids": list(cold_split.train_ids),  # type: ignore[attr-defined]
            "val_ids": list(cold_split.val_ids),  # type: ignore[attr-defined]
            "test_ids": list(cold_split.test_ids),  # type: ignore[attr-defined]
            "method": cold_split.method,  # type: ignore[attr-defined]
            "seed": cold_split.seed,  # type: ignore[attr-defined]
            "notes": cold_split.notes,  # type: ignore[attr-defined]
        },
        "clusters": cluster_of,
        "time_based": {k: list(v) for k, v in (time_split or {}).items()},
        "counts": counts,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog",
        type=Path,
        default=DATA_DIR / "peptide_affinity_catalog_v2.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DATA_DIR / "splits_v2.json",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--min-test-clusters", type=int, default=10)
    args = parser.parse_args()

    entries = load_peptide_affinity_catalog(args.catalog)
    cluster_probe = cluster_entries(entries)
    n_clusters = len(set(cluster_probe.values()))
    min_tc = min(args.min_test_clusters, max(2, n_clusters // 4))
    split, cluster_of = cold_start_split(
        entries, seed=args.seed, min_test_clusters=min_tc
    )

    time_split = None
    try:
        time_split = time_based_split(entries)
    except ValueError as exc:
        print(f"time-based split skipped: {exc}")

    counts = {
        "n": len(entries),
        "n_clusters": n_clusters,
        "train": len(split.train_ids),
        "val": len(split.val_ids),
        "test": len(split.test_ids),
        "test_clusters": len({cluster_of[i] for i in split.test_ids}),
    }
    write_split_artifact(
        args.out,
        cold_split=split,
        cluster_of=cluster_of,
        time_split=time_split,
        counts=counts,
    )
    print(
        f"cold-start partitions: train={counts['train']} val={counts['val']} "
        f"test={counts['test']} (test_clusters={counts['test_clusters']}) → {args.out}"
    )


if __name__ == "__main__":
    main()

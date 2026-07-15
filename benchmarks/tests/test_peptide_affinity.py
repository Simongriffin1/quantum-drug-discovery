"""Tests for expanded peptide-affinity catalog loader (offline fixtures)."""

from __future__ import annotations

from pathlib import Path

import pytest

from peptide_affinity.load import catalog_qc_summary, load_peptide_affinity_catalog
from peptide_affinity.splits import (
    assert_zero_cluster_overlap,
    cluster_entries,
    cold_start_split,
)

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "peptide_affinity"
    / "fixtures"
    / "peptide_affinity_ci_v1.tsv"
)


def test_loader_parses_ci_fixture_and_reports_n() -> None:
    entries = load_peptide_affinity_catalog(FIXTURE)
    summary = catalog_qc_summary(entries)
    assert summary["n"] == len(entries) >= 5
    assert all(5 <= e.peptide_len <= 50 for e in entries)
    assert all(e.pKd > 0 for e in entries)


def test_loader_rejects_duplicates_idempotently(tmp_path: Path) -> None:
    text = FIXTURE.read_text(encoding="utf-8")
    dup = tmp_path / "dup.tsv"
    lines = text.strip().splitlines()
    dup.write_text("\n".join(lines + [lines[-1]]) + "\n", encoding="utf-8")
    entries = load_peptide_affinity_catalog(dup)
    ids = [e.record_id for e in entries]
    assert len(ids) == len(set(ids))


def test_loader_rejects_bad_pkd_consistency(tmp_path: Path) -> None:
    bad = tmp_path / "bad.tsv"
    bad.write_text(
        "record_id\tpdb_id\treceptor_seq\tpeptide_seq\tpeptide_len\tresolution\t"
        "affinity_value\taffinity_type\tpKd\tsource\tstructure_path\tpeptide_chain\t"
        "receptor_chains\tdeposit_year\tnotes\n"
        "BAD1\t1BD2\t"
        + ("A" * 40)
        + "\tLLFGYPVYV\t9\t2.5\t1e-6\tKd\t9.0\tci\t\tC\tA\t2000\tbad\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="no valid peptide affinity"):
        load_peptide_affinity_catalog(bad)


def test_cold_start_split_zero_cluster_overlap() -> None:
    entries = load_peptide_affinity_catalog(FIXTURE)
    # Small fixture: relax min_test_clusters
    cluster_of = cluster_entries(entries, identity_threshold=0.30)
    # Force distinct clusters by using high identity threshold for this tiny set
    split, cluster_of = cold_start_split(
        entries,
        test_fraction=0.34,
        val_fraction=0.0,
        seed=1,
        identity_threshold=0.99,
        min_test_clusters=1,
    )
    assert_zero_cluster_overlap(split.train_ids, split.test_ids, cluster_of)
    assert len(split.train_ids) >= 1
    assert len(split.test_ids) >= 1
    print(
        f"partitions train={len(split.train_ids)} test={len(split.test_ids)} "
        f"clusters={len(set(cluster_of.values()))}"
    )

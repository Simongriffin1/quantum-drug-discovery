"""Benchmark loader + split tests (fixture-only, no network)."""

from __future__ import annotations

import math

import pytest

from peptideforge_benchmarks.pdbbind import affinity_to_pk, load_pdbbind_peptide_affinity
from peptideforge_benchmarks.models import AffinityUnit
from peptideforge_benchmarks.skempi import load_skempi_ddg
from peptideforge_benchmarks.splits import (
    greedy_identity_clusters,
    homology_aware_split,
    load_cluster_assignments,
    sequence_identity,
)


def test_load_pdbbind_peptide_affinity_fixture() -> None:
    records = load_pdbbind_peptide_affinity()
    assert len(records) == 12
    assert records[0].pdb_id == "1A1N"
    assert records[0].pk == pytest.approx(8.0)
    # No network — loader is local TSV only
    assert all(r.source == "pdbbind_peptide_subset" for r in records)


def test_affinity_to_pk_consistency() -> None:
    assert affinity_to_pk(1e-8, AffinityUnit.KD) == pytest.approx(8.0)
    assert affinity_to_pk(8.0, AffinityUnit.PK) == pytest.approx(8.0)
    with pytest.raises(ValueError):
        affinity_to_pk(0.0, AffinityUnit.KD)


def test_load_skempi_ddg_fixture() -> None:
    records = load_skempi_ddg()
    assert len(records) == 8
    assert records[0].ddg_kcal_mol == pytest.approx(0.80)
    assert records[0].mutant == "AI38A"


def test_missing_fixture_fails_loud(tmp_path) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(FileNotFoundError):
        load_pdbbind_peptide_affinity(tmp_path / "does_not_exist.tsv")


def test_sequence_identity_identical() -> None:
    assert sequence_identity("ACDEF", "ACDEF") == pytest.approx(1.0)
    assert sequence_identity("ACDEF", "AAAAA") == pytest.approx(0.2)


def test_homology_aware_split_no_cluster_leakage() -> None:
    records = load_pdbbind_peptide_affinity()
    clusters = {r.record_id: r.cluster_id for r in records if r.cluster_id}
    assert len(clusters) == 12
    split = homology_aware_split(
        tuple(clusters.keys()),
        clusters,  # type: ignore[arg-type]
        test_fraction=0.25,
        seed=0,
    )
    train_c = {clusters[i] for i in split.train_ids}
    test_c = {clusters[i] for i in split.test_ids}
    assert train_c.isdisjoint(test_c)
    assert len(split.train_ids) + len(split.test_ids) == 12


def test_load_cluster_assignments() -> None:
    assignments = load_cluster_assignments()
    assert assignments["PP001"] == "c0"
    assert assignments["SK001"] == "c0"


def test_greedy_identity_clusters() -> None:
    seqs = {
        "a": "ACDEFGHIKL",
        "b": "ACDEFGHIKL",  # identical → same cluster
        "c": "YYYYYYYYYY",  # distant
    }
    clusters = greedy_identity_clusters(seqs, identity_threshold=0.9)
    assert clusters["a"] == clusters["b"]
    assert clusters["c"] != clusters["a"]


def test_pk_values_match_recompute() -> None:
    for rec in load_pdbbind_peptide_affinity():
        recomputed = affinity_to_pk(rec.affinity_value, rec.affinity_unit)
        assert recomputed == pytest.approx(rec.pk, abs=1e-3)
        assert math.isfinite(recomputed)

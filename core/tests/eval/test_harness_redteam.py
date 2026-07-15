"""Oracle-validity harness + red-team tests."""

from __future__ import annotations

import pytest
from peptideforge_benchmarks.pdbbind import load_pdbbind_peptide_affinity
from peptideforge_benchmarks.splits import homology_aware_split

from peptideforge.contracts.models import (
    ComplexStructure,
    OracleResult,
    OracleTier,
    PeptideCandidate,
)
from peptideforge.eval.harness import (
    DEFAULT_THRESHOLDS,
    PredictionLabelPair,
    evaluate_predictions,
    run_oracle_on_complexes,
)
from peptideforge.eval.metrics import rmse, spearman_rho
from peptideforge.eval.redteam import (
    label_shuffle_control,
    leakage_audit,
    run_red_team,
    trivial_baseline_check,
)


@pytest.mark.eval
def test_spearman_perfect_and_anticorrelated() -> None:
    assert spearman_rho([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]) == pytest.approx(1.0)
    assert spearman_rho([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0]) == pytest.approx(-1.0)


@pytest.mark.eval
def test_rmse() -> None:
    assert rmse([1.0, 2.0], [1.0, 2.0]) == pytest.approx(0.0)
    # rms of errors [3,4] = sqrt((9+16)/2) = sqrt(12.5)
    assert rmse([0.0, 0.0], [3.0, 4.0]) == pytest.approx((12.5) ** 0.5)


@pytest.mark.eval
def test_evaluate_predictions_report() -> None:
    pairs = [
        PredictionLabelPair(record_id="a", predicted=1.0, experimental=1.0),
        PredictionLabelPair(record_id="b", predicted=2.0, experimental=2.0),
        PredictionLabelPair(record_id="c", predicted=3.0, experimental=2.5),
    ]
    report = evaluate_predictions(
        pairs, spearman_threshold=DEFAULT_THRESHOLDS.oracle_affinity_spearman
    )
    assert report.n == 3
    assert report.spearman > 0.9
    assert report.passed_spearman_threshold


@pytest.mark.eval
def test_run_oracle_on_complexes_invert_sign() -> None:
    """Toy oracle: returns -pK as ΔG so -value recovers pK."""

    class _ToyOracle:
        def __init__(self, pk_by_id: dict) -> None:  # type: ignore[type-arg]
            self._pk = pk_by_id

        def evaluate(self, complex_structure, *, tier=None, cost_cap=None):  # type: ignore[no-untyped-def]
            pk = self._pk[complex_structure.candidate_id]
            return OracleResult(
                candidate_id=complex_structure.candidate_id,
                complex_id=complex_structure.complex_id,
                value=-pk,
                uncertainty=0.1,
                cost_estimate=0.01,
                tier=tier or OracleTier.SYNTHETIC,
            )

    labels: dict = {}
    complexes = []
    for pk in (5.0, 6.0, 7.0, 8.0):
        cand = PeptideCandidate(sequence="ACDEFGHIKLM", generation_method="synthetic_eval")
        labels[cand.candidate_id] = pk
        complexes.append(
            ComplexStructure(
                candidate_id=cand.candidate_id,
                target_id="synthetic_target",
                sequence=cand.sequence,
                pdb_text="HEADER synthetic\nEND\n",
                confidence=0.9,
                fold_method="synthetic_folder",
            )
        )

    report, results = run_oracle_on_complexes(
        _ToyOracle(labels),
        complexes,
        labels,
        spearman_threshold=0.4,
        invert_sign=True,
    )
    assert len(results) == 4
    assert report.spearman == pytest.approx(1.0)
    assert report.passed_spearman_threshold


@pytest.mark.eval
def test_label_shuffle_collapses_perfect_predictor() -> None:
    pairs = [
        PredictionLabelPair(record_id=str(i), predicted=float(i), experimental=float(i))
        for i in range(10)
    ]
    true_rho = spearman_rho([p.predicted for p in pairs], [p.experimental for p in pairs])
    passed, shuffled_rho = label_shuffle_control(pairs, seed=1)
    assert passed
    # Either |ρ| collapses or drops by the pre-registered margin (ACCEPTANCE.md)
    assert abs(shuffled_rho) < 0.20 or (true_rho - shuffled_rho) >= 0.40


@pytest.mark.eval
def test_trivial_baseline_beaten_by_real_signal() -> None:
    pairs = [
        PredictionLabelPair(record_id=str(i), predicted=float(i), experimental=float(i))
        for i in range(8)
    ]
    passed, model_rho, baseline_rho = trivial_baseline_check(pairs)
    assert passed
    assert model_rho > baseline_rho


@pytest.mark.eval
def test_leakage_audit_flags_overlapping_ids() -> None:
    ok, findings = leakage_audit(["a", "b", "c"], ["c", "d"])
    assert not ok
    assert any("record_id overlap" in f for f in findings)


@pytest.mark.eval
def test_leakage_audit_flags_cluster_and_sequence_leak() -> None:
    ok, findings = leakage_audit(
        ["a", "b"],
        ["c"],
        train_clusters={"a": "c0", "b": "c1"},
        test_clusters={"c": "c0"},
        train_sequences={"a": "ACDEF", "b": "GGGGG"},
        test_sequences={"c": "ACDEF"},
    )
    assert not ok
    assert any("cluster overlap" in f for f in findings)
    assert any("identical sequence" in f for f in findings)


@pytest.mark.eval
def test_red_team_flags_deliberately_leaked_toy_set() -> None:
    """Deliberately leaky: train/test share record_ids.

    Red-team MUST fail (leakage_passed=False). synthetic_* plumbing only.
    """
    records = load_pdbbind_peptide_affinity()
    ids = [r.record_id for r in records[:6]]
    sequ = {r.record_id: r.peptide_sequence for r in records[:6]}
    clusters = {r.record_id: r.cluster_id or "x" for r in records[:6]}

    pairs = [
        PredictionLabelPair(
            record_id=r.record_id,
            predicted=r.pk + 0.1,
            experimental=r.pk,
        )
        for r in records[:6]
    ]

    report = run_red_team(
        pairs,
        train_ids=ids,
        test_ids=ids,  # deliberate ID leakage
        train_clusters=clusters,
        test_clusters=clusters,
        train_sequences=sequ,
        test_sequences=sequ,
        seed=0,
    )
    assert not report.passed
    assert not report.leakage_passed
    assert len(report.leakage_findings) >= 1


@pytest.mark.eval
def test_red_team_passes_clean_homology_split() -> None:
    records = load_pdbbind_peptide_affinity()
    clusters = {r.record_id: r.cluster_id for r in records if r.cluster_id}
    split = homology_aware_split(
        tuple(clusters.keys()),
        {k: v for k, v in clusters.items() if v is not None},  # type: ignore[misc]
        test_fraction=0.25,
        seed=2,
    )
    pairs = [
        PredictionLabelPair(record_id=r.record_id, predicted=r.pk, experimental=r.pk)
        for r in records
    ]
    test_pairs = [p for p in pairs if p.record_id in set(split.test_ids)]
    if len(test_pairs) < 3:
        test_pairs = pairs[:6]

    report = run_red_team(
        test_pairs,
        train_ids=split.train_ids,
        test_ids=split.test_ids,
        train_clusters={i: clusters[i] for i in split.train_ids},  # type: ignore[misc]
        test_clusters={i: clusters[i] for i in split.test_ids},  # type: ignore[misc]
        seed=0,
    )
    assert report.leakage_passed
    assert report.passed


@pytest.mark.eval
def test_thresholds_match_acceptance_doc() -> None:
    assert DEFAULT_THRESHOLDS.oracle_affinity_spearman == 0.40
    assert DEFAULT_THRESHOLDS.oracle_ddg_spearman == 0.30
    assert DEFAULT_THRESHOLDS.surrogate_ece == 0.10

"""Tests for P7 surrogate + calibrated UQ."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from peptideforge.contracts.models import Candidates, PeptideCandidate
from peptideforge.contracts.protocols import Surrogate
from peptideforge.eval.harness import DEFAULT_THRESHOLDS
from peptideforge.surrogate import (
    DeepEnsembleSurrogate,
    SurrogateNotFittedError,
    evaluate_calibration,
    make_oracle_labels,
    run_surrogate_acceptance,
    sequence_features,
    synthetic_physics_label,
    write_acceptance_report,
)


def _synthetic_pool(n: int, *, seed: int = 0, prefix: str = "A") -> list[PeptideCandidate]:
    """Diverse synthetic peptides with gravy/charge variation (plumbing only)."""
    rng = random.Random(seed)
    # Residue pools for hydrophobic vs polar/charged
    hydro = list("FLIVMWA")
    polar = list("DESTNQKR")
    items: list[PeptideCandidate] = []
    for i in range(n):
        length = rng.randint(7, 14)
        if i % 2 == 0:
            seq = "".join(rng.choice(hydro) for _ in range(length))
        else:
            seq = "".join(rng.choice(polar) for _ in range(length))
        # Inject a distinctive residue from prefix cluster tag to avoid identical seqs
        tag = prefix[0] if prefix[0] in "ACDEFGHIKLMNPQRSTVWY" else "A"
        seq = (tag + seq)[:length] if length >= 5 else seq
        if len(seq) < 5:
            seq = seq + "AAAAA"
            seq = seq[:5]
        items.append(
            PeptideCandidate(
                sequence=seq.upper(),
                generation_method="synthetic_surrogate_pool",
                metadata={"cluster": prefix, "idx": i},
            )
        )
    # Deduplicate sequences
    seen: set[str] = set()
    unique: list[PeptideCandidate] = []
    for c in items:
        if c.sequence in seen:
            # mutate one position
            chars = list(c.sequence)
            chars[1] = "G" if chars[1] != "G" else "P"
            c = PeptideCandidate(
                sequence="".join(chars),
                generation_method="synthetic_surrogate_pool",
                metadata=c.metadata,
            )
        seen.add(c.sequence)
        unique.append(c)
    return unique


@pytest.mark.eval
def test_sequence_features_fixed_dim() -> None:
    feats = sequence_features("LLFGYPVYV")
    assert len(feats) == 23  # 20 AA + length + gravy + charge


@pytest.mark.eval
def test_surrogate_protocol() -> None:
    s = DeepEnsembleSurrogate(n_ensemble=4)
    assert isinstance(s, Surrogate)


@pytest.mark.eval
def test_predict_before_fit_raises() -> None:
    s = DeepEnsembleSurrogate()
    cands = Candidates(
        items=(
            PeptideCandidate(sequence="AAAAA", generation_method="synthetic_test"),
        )
    )
    with pytest.raises(SurrogateNotFittedError):
        s.predict(cands)


@pytest.mark.eval
def test_surrogate_learns_synthetic_labels_and_meets_calibration_gate(
    tmp_path: Path,
) -> None:
    """On a learnable synthetic_* oracle, ECE < 0.10 and red-team pass."""
    train = _synthetic_pool(40, seed=1, prefix="T")
    test = _synthetic_pool(20, seed=2, prefix="E")
    # Ensure no identical sequence leakage
    train_seqs = {c.sequence for c in train}
    test = [c for c in test if c.sequence not in train_seqs]
    assert len(test) >= 12

    train_batch = Candidates(items=tuple(train), seed=1)
    test_batch = Candidates(items=tuple(test), seed=2)

    train_vals = {
        c.candidate_id: synthetic_physics_label(c.sequence, noise=0.05, seed=11)
        for c in train
    }
    test_vals = {
        c.candidate_id: synthetic_physics_label(c.sequence, noise=0.05, seed=11)
        for c in test
    }
    labels = make_oracle_labels(train_batch, train_vals)

    train_clusters = {str(c.candidate_id): "cluster_train" for c in train}
    test_clusters = {str(c.candidate_id): "cluster_test" for c in test}

    surrogate = DeepEnsembleSurrogate(
        n_ensemble=6,
        l2=0.5,
        coverage_target=0.90,
        calib_fraction=0.25,
        objective_name="synthetic_binding",
    )
    report = run_surrogate_acceptance(
        surrogate,
        train_batch,
        labels,
        test_batch,
        test_vals,
        train_ids=tuple(str(c.candidate_id) for c in train),
        test_ids=tuple(str(c.candidate_id) for c in test),
        train_clusters=train_clusters,
        test_clusters=test_clusters,
        train_sequences={str(c.candidate_id): c.sequence for c in train},
        test_sequences={str(c.candidate_id): c.sequence for c in test},
        seed=42,
        thresholds=DEFAULT_THRESHOLDS,
        data_version="synthetic_surrogate_v1",
        notes="synthetic_learnable_oracle_not_physics",
    )
    out = tmp_path / "surrogate_acceptance.json"
    write_acceptance_report(report, out)
    assert out.is_file()

    assert report.calibration.ece < DEFAULT_THRESHOLDS.surrogate_ece
    assert report.calibration.passed
    assert report.red_team.leakage_passed
    assert report.red_team.trivial_baseline_passed
    assert report.red_team.label_shuffle_passed
    assert report.accepted
    assert report.red_team.model_rho > 0.3


@pytest.mark.eval
def test_calibrated_intervals_ordered() -> None:
    train = _synthetic_pool(24, seed=3, prefix="X")
    batch = Candidates(items=tuple(train), seed=3)
    vals = {c.candidate_id: synthetic_physics_label(c.sequence) for c in train}
    labels = make_oracle_labels(batch, vals)
    s = DeepEnsembleSurrogate(n_ensemble=4)
    s.register(batch)
    s.fit(tuple(c.candidate_id for c in train), labels, seed=0)
    # Hold out a few for predict shape check
    hold = Candidates(items=tuple(train[:5]), seed=0)
    vectors = s.predict(hold)
    assert len(vectors) == 5
    for vec in vectors:
        pred = vec.predictions[0]
        assert pred.lower <= pred.mean <= pred.upper
        assert pred.epistemic_std >= 0.0


@pytest.mark.eval
def test_ece_constant_width() -> None:
    # Perfectly calibrated: all observations inside intervals at exact coverage
    lowers = [-1.0] * 20
    uppers = [1.0] * 20
    # 18/20 = 0.9 coverage
    obs = [-0.5] * 18 + [2.0, 2.0]
    report = evaluate_calibration(
        lowers, uppers, obs, coverage_target=0.90, ece_threshold=0.10
    )
    assert report.empirical_coverage == pytest.approx(0.9)
    assert report.ece == pytest.approx(0.0)
    assert report.passed

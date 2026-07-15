"""Red-team controls: label shuffle, trivial baseline, leakage audit.

Every claimed correlation must pass these before it is believed
(see ACCEPTANCE.md and CURSOR_PROJECT_CONTEXT.md).
"""

from __future__ import annotations

import random
from collections.abc import Sequence

from pydantic import BaseModel, ConfigDict

from peptideforge.eval.harness import DEFAULT_THRESHOLDS, PredictionLabelPair, Thresholds
from peptideforge.eval.metrics import spearman_rho


class RedTeamReport(BaseModel):
    """Pass/fail summary for red-team battery."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    passed: bool
    label_shuffle_passed: bool
    label_shuffle_rho: float
    trivial_baseline_passed: bool
    model_rho: float
    baseline_rho: float
    leakage_passed: bool
    leakage_findings: tuple[str, ...] = ()
    notes: str | None = None


def label_shuffle_control(
    pairs: Sequence[PredictionLabelPair],
    *,
    seed: int = 0,
    thresholds: Thresholds = DEFAULT_THRESHOLDS,
) -> tuple[bool, float]:
    """Permute experimental labels; Spearman must collapse if predictions are real.

    Returns (passed, shuffled_spearman).
    """
    if len(pairs) < 3:
        raise ValueError("label_shuffle_control needs ≥ 3 pairs")
    predicted = [p.predicted for p in pairs]
    experimental = [p.experimental for p in pairs]
    true_rho = spearman_rho(predicted, experimental)

    rng = random.Random(seed)
    shuffled = list(experimental)
    # Ensure a real shuffle (not identity) when possible
    for _ in range(100):
        rng.shuffle(shuffled)
        if shuffled != experimental:
            break
    else:
        raise RuntimeError("could not produce a non-identity label shuffle")

    shuffled_rho = spearman_rho(predicted, shuffled)
    near_zero = abs(shuffled_rho) < thresholds.label_shuffle_max_abs_rho
    dropped = (true_rho - shuffled_rho) >= thresholds.label_shuffle_min_drop
    # If true_rho itself is weak, require near-zero shuffled correlation
    passed = near_zero or dropped
    return passed, shuffled_rho


def trivial_baseline_predictions(pairs: Sequence[PredictionLabelPair]) -> list[float]:
    """Mean-constant baseline (no sequence features) — trivial predictor."""
    mean_y = sum(p.experimental for p in pairs) / len(pairs)
    return [mean_y for _ in pairs]


def length_baseline_from_meta(
    pairs: Sequence[PredictionLabelPair],
    lengths: Sequence[int],
) -> list[float]:
    """Rank by sequence length (often a weakly informative trivial heuristic)."""
    if len(lengths) != len(pairs):
        raise ValueError("lengths must align with pairs")
    return [float(length) for length in lengths]


def trivial_baseline_check(
    pairs: Sequence[PredictionLabelPair],
    *,
    baseline_predicted: Sequence[float] | None = None,
    thresholds: Thresholds = DEFAULT_THRESHOLDS,
) -> tuple[bool, float, float]:
    """Model Spearman must beat trivial baseline by Δρ ≥ threshold.

    Returns (passed, model_rho, baseline_rho).
    """
    predicted = [p.predicted for p in pairs]
    experimental = [p.experimental for p in pairs]
    model_rho = spearman_rho(predicted, experimental)
    baseline = (
        list(baseline_predicted)
        if baseline_predicted is not None
        else trivial_baseline_predictions(pairs)
    )
    # Constant baseline → Spearman undefined; treat as rho = 0.0
    try:
        baseline_rho = spearman_rho(baseline, experimental)
    except ValueError:
        baseline_rho = 0.0
    passed = (model_rho - baseline_rho) >= thresholds.trivial_baseline_min_delta_rho
    return passed, model_rho, baseline_rho


def leakage_audit(
    train_ids: Sequence[str],
    test_ids: Sequence[str],
    *,
    train_clusters: dict[str, str] | None = None,
    test_clusters: dict[str, str] | None = None,
    train_sequences: dict[str, str] | None = None,
    test_sequences: dict[str, str] | None = None,
) -> tuple[bool, tuple[str, ...]]:
    """Detect train/test identity leakage (IDs, clusters, identical sequences)."""
    findings: list[str] = []
    train_set = set(train_ids)
    test_set = set(test_ids)
    overlap_ids = train_set & test_set
    if overlap_ids:
        findings.append(f"record_id overlap: {sorted(overlap_ids)[:10]}")

    if train_clusters is not None and test_clusters is not None:
        train_c = {train_clusters[i] for i in train_ids if i in train_clusters}
        test_c = {test_clusters[i] for i in test_ids if i in test_clusters}
        overlap_c = train_c & test_c
        if overlap_c:
            findings.append(f"cluster overlap: {sorted(overlap_c)[:10]}")

    if train_sequences is not None and test_sequences is not None:
        train_seq_set = {train_sequences[i] for i in train_ids if i in train_sequences}
        leaked = [
            tid
            for tid in test_ids
            if tid in test_sequences and test_sequences[tid] in train_seq_set
        ]
        if leaked:
            findings.append(f"identical sequence leakage: {leaked[:10]}")

    return (len(findings) == 0, tuple(findings))


def run_red_team(
    pairs: Sequence[PredictionLabelPair],
    *,
    train_ids: Sequence[str] | None = None,
    test_ids: Sequence[str] | None = None,
    train_clusters: dict[str, str] | None = None,
    test_clusters: dict[str, str] | None = None,
    train_sequences: dict[str, str] | None = None,
    test_sequences: dict[str, str] | None = None,
    baseline_predicted: Sequence[float] | None = None,
    seed: int = 0,
    thresholds: Thresholds = DEFAULT_THRESHOLDS,
) -> RedTeamReport:
    """Run full red-team battery; ``passed`` is True only if all controls pass."""
    shuffle_ok, shuffle_rho = label_shuffle_control(pairs, seed=seed, thresholds=thresholds)
    triv_ok, model_rho, baseline_rho = trivial_baseline_check(
        pairs, baseline_predicted=baseline_predicted, thresholds=thresholds
    )

    if train_ids is not None and test_ids is not None:
        leak_ok, findings = leakage_audit(
            train_ids,
            test_ids,
            train_clusters=train_clusters,
            test_clusters=test_clusters,
            train_sequences=train_sequences,
            test_sequences=test_sequences,
        )
    else:
        leak_ok, findings = True, ()

    return RedTeamReport(
        passed=shuffle_ok and triv_ok and leak_ok,
        label_shuffle_passed=shuffle_ok,
        label_shuffle_rho=shuffle_rho,
        trivial_baseline_passed=triv_ok,
        model_rho=model_rho,
        baseline_rho=baseline_rho,
        leakage_passed=leak_ok,
        leakage_findings=findings,
    )

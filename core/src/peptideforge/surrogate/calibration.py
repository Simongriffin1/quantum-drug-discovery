"""Calibration metrics: coverage, ECE, reliability diagram bins."""

from __future__ import annotations

from collections.abc import Sequence

from pydantic import BaseModel, ConfigDict, Field


class ReliabilityBin(BaseModel):
    """One bin of a reliability diagram for interval coverage."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    bin_index: int
    n: int
    predicted_coverage: float
    empirical_coverage: float
    mean_interval_width: float


class CalibrationReport(BaseModel):
    """Surrogate calibration summary — real metrics only, never fabricated.

    Gate: ECE < ``ece_threshold`` at ``coverage_target`` (see ACCEPTANCE.md).
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    n: int
    coverage_target: float
    empirical_coverage: float
    ece: float
    ece_threshold: float
    passed: bool
    mean_interval_width: float
    reliability_bins: tuple[ReliabilityBin, ...]
    notes: str | None = None


def interval_covered(lower: float, upper: float, observed: float) -> bool:
    return lower <= observed <= upper


def empirical_coverage(
    lowers: Sequence[float],
    uppers: Sequence[float],
    observed: Sequence[float],
) -> float:
    if not (len(lowers) == len(uppers) == len(observed)):
        raise ValueError("interval arrays must align")
    if not observed:
        raise ValueError("need ≥ 1 observation")
    hits = sum(
        1 for lo, hi, y in zip(lowers, uppers, observed, strict=True) if interval_covered(lo, hi, y)
    )
    return hits / len(observed)


def expected_calibration_error(
    lowers: Sequence[float],
    uppers: Sequence[float],
    observed: Sequence[float],
    *,
    coverage_target: float,
    n_bins: int = 10,
) -> tuple[float, tuple[ReliabilityBin, ...]]:
    """ECE over bins of predicted interval width (proxy for confidence).

    For conformal intervals with a single global width, width is constant and ECE
    collapses to |empirical_coverage − coverage_target|. With heterogeneous
    widths (ensemble std absorbed into conformal expansion), we bin by width.
    """
    if n_bins < 1:
        raise ValueError("n_bins must be ≥ 1")
    if not (len(lowers) == len(uppers) == len(observed)):
        raise ValueError("interval arrays must align")
    n = len(observed)
    if n == 0:
        raise ValueError("need ≥ 1 observation")

    widths = [hi - lo for lo, hi in zip(lowers, uppers, strict=True)]
    # If all widths equal, single-bin ECE
    if max(widths) - min(widths) < 1e-12:
        cov = empirical_coverage(lowers, uppers, observed)
        ece = abs(cov - coverage_target)
        bin0 = ReliabilityBin(
            bin_index=0,
            n=n,
            predicted_coverage=coverage_target,
            empirical_coverage=cov,
            mean_interval_width=widths[0],
        )
        return ece, (bin0,)

    w_min, w_max = min(widths), max(widths)
    edges = [w_min + (w_max - w_min) * i / n_bins for i in range(n_bins + 1)]
    bins: list[ReliabilityBin] = []
    ece = 0.0
    for b in range(n_bins):
        lo_e, hi_e = edges[b], edges[b + 1]
        idxs = [
            i
            for i, w in enumerate(widths)
            if (w >= lo_e and w < hi_e) or (b == n_bins - 1 and w <= hi_e)
        ]
        if not idxs:
            continue
        sub_lo = [lowers[i] for i in idxs]
        sub_hi = [uppers[i] for i in idxs]
        sub_y = [observed[i] for i in idxs]
        cov = empirical_coverage(sub_lo, sub_hi, sub_y)
        weight = len(idxs) / n
        ece += weight * abs(cov - coverage_target)
        bins.append(
            ReliabilityBin(
                bin_index=b,
                n=len(idxs),
                predicted_coverage=coverage_target,
                empirical_coverage=cov,
                mean_interval_width=sum(widths[i] for i in idxs) / len(idxs),
            )
        )
    return ece, tuple(bins)


def evaluate_calibration(
    lowers: Sequence[float],
    uppers: Sequence[float],
    observed: Sequence[float],
    *,
    coverage_target: float = 0.90,
    ece_threshold: float = 0.10,
    n_bins: int = 10,
    notes: str | None = None,
) -> CalibrationReport:
    """Compute coverage + ECE and compare against the pre-registered gate."""
    cov = empirical_coverage(lowers, uppers, observed)
    ece, bins = expected_calibration_error(
        lowers, uppers, observed, coverage_target=coverage_target, n_bins=n_bins
    )
    widths = [hi - lo for lo, hi in zip(lowers, uppers, strict=True)]
    mean_w = sum(widths) / len(widths)
    return CalibrationReport(
        n=len(observed),
        coverage_target=coverage_target,
        empirical_coverage=cov,
        ece=ece,
        ece_threshold=ece_threshold,
        passed=ece < ece_threshold,
        mean_interval_width=mean_w,
        reliability_bins=bins,
        notes=notes,
    )

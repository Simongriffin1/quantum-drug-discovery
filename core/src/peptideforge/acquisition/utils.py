"""Constraint filtering and prediction utilities for acquisition."""

from __future__ import annotations

from collections.abc import Iterable
from uuid import UUID

from peptideforge.contracts.models import Candidates, ObjectiveVector, RankedCandidate


def align_predictions(
    pool: Candidates,
    predictions: tuple[ObjectiveVector, ...],
) -> dict[UUID, ObjectiveVector]:
    """Map candidate_id → ObjectiveVector; raise if pool members lack predictions."""
    by_id = {v.candidate_id: v for v in predictions}
    missing = [c.candidate_id for c in pool.items if c.candidate_id not in by_id]
    if missing:
        raise KeyError(f"predictions missing for candidates: {missing[:5]}")
    return by_id


def objective_means(vec: ObjectiveVector) -> tuple[float, ...]:
    return tuple(p.mean for p in vec.predictions)


def objective_stds(vec: ObjectiveVector) -> tuple[float, ...]:
    """Epistemic std; fall back to half conformal width if std is zero."""
    out: list[float] = []
    for p in vec.predictions:
        if p.epistemic_std > 0.0:
            out.append(p.epistemic_std)
        else:
            half = max(0.0, (p.upper - p.lower) / 2.0)
            out.append(half if half > 0.0 else 1e-6)
    return tuple(out)


def apply_constraints(
    pool: Candidates,
    pred_map: dict[UUID, ObjectiveVector],
    constraints: dict[str, object] | None,
) -> tuple[list[UUID], dict[UUID, str]]:
    """Return (eligible_ids, constrained_out reasons).

    Supported constraints:
    - ``max_<objective_name>``: float upper bound on predicted mean
    - ``min_<objective_name>``: float lower bound on predicted mean
    - ``exclude_ids``: iterable of UUID / str
    """
    constrained: dict[UUID, str] = {}
    exclude: set[UUID] = set()
    if constraints:
        raw_ex = constraints.get("exclude_ids")
        if isinstance(raw_ex, Iterable) and not isinstance(raw_ex, (str, bytes)):
            for item in raw_ex:
                exclude.add(item if isinstance(item, UUID) else UUID(str(item)))

    eligible: list[UUID] = []
    for cand in pool.items:
        cid = cand.candidate_id
        if cid in exclude:
            constrained[cid] = "exclude_ids"
            continue
        reason = _check_bounds(pred_map[cid], constraints)
        if reason is not None:
            constrained[cid] = reason
            continue
        eligible.append(cid)
    return eligible, constrained


def _check_bounds(
    vec: ObjectiveVector,
    constraints: dict[str, object] | None,
) -> str | None:
    if not constraints:
        return None
    for pred in vec.predictions:
        max_key = f"max_{pred.objective_name}"
        min_key = f"min_{pred.objective_name}"
        if max_key in constraints:
            bound = _as_float(constraints[max_key])
            if pred.mean > bound:
                return f"{max_key}={bound}"
        if min_key in constraints:
            bound = _as_float(constraints[min_key])
            if pred.mean < bound:
                return f"{min_key}={bound}"
    return None


def _as_float(value: object) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    raise TypeError(f"expected numeric constraint, got {type(value).__name__}")


def cost_per_candidate(constraints: dict[str, object] | None) -> float:
    if constraints is None:
        return 1.0
    raw = constraints.get("cost_per_candidate", 1.0)
    cost = _as_float(raw)
    if cost <= 0.0:
        raise ValueError("cost_per_candidate must be > 0")
    return cost


def pack_batch(
    *,
    selected: list[tuple[UUID, float]],
    constrained: dict[UUID, str],
    batch_size: int,
    budget_remaining: float,
    method: str,
    seed: int | None,
) -> tuple[RankedCandidate, ...]:
    """Build RankedCandidate list: selected first (by rank), then constrained-out."""
    _ = budget_remaining, seed  # reserved for provenance in later versions
    ranked: list[RankedCandidate] = []
    for i, (cid, score) in enumerate(selected):
        ranked.append(
            RankedCandidate(
                candidate_id=cid,
                acquisition_score=score,
                rank=i,
                constrained_out=False,
            )
        )
    offset = len(ranked)
    for j, (cid, reason) in enumerate(sorted(constrained.items(), key=lambda kv: str(kv[0]))):
        ranked.append(
            RankedCandidate(
                candidate_id=cid,
                acquisition_score=float("-inf"),
                rank=offset + j,
                constrained_out=True,
                reason=reason,
            )
        )
    if len([r for r in ranked if not r.constrained_out]) > batch_size:
        raise ValueError("internal error: more selected than batch_size")
    return tuple(ranked)

"""Persist and restore LoopState (reproducible from seed + config)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from uuid import UUID

from peptideforge.contracts.models import LoopState, Provenance


def save_loop_state(state: LoopState, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(state.model_dump_json(indent=2), encoding="utf-8")


def load_loop_state(path: Path) -> LoopState:
    if not path.is_file():
        raise FileNotFoundError(f"LoopState not found: {path}")
    return LoopState.model_validate_json(path.read_text(encoding="utf-8"))


def write_state_history(states: list[LoopState], directory: Path) -> None:
    """Write ``iteration_NNNN.json`` plus ``latest.json``."""
    directory.mkdir(parents=True, exist_ok=True)
    for state in states:
        save_loop_state(state, directory / f"iteration_{state.iteration:04d}.json")
    if states:
        save_loop_state(states[-1], directory / "latest.json")
        (directory / "history_index.json").write_text(
            json.dumps(
                {
                    "n_iterations": len(states),
                    "campaign_id": str(states[-1].campaign_id),
                    "seed": states[-1].seed,
                },
                indent=2,
            ),
            encoding="utf-8",
        )


def make_loop_state(
    *,
    campaign_id: UUID,
    iteration: int,
    seed: int,
    config: dict[str, Any],
    candidate_ids: tuple[UUID, ...],
    labeled_ids: tuple[UUID, ...],
    oracle_calls: int,
    total_cost: float,
    pareto_front_ids: tuple[UUID, ...] = (),
    surrogate_version: str | None = None,
    data_version: str | None = None,
    status: str = "running",
    notes: str | None = None,
) -> LoopState:
    return LoopState(
        campaign_id=campaign_id,
        iteration=iteration,
        seed=seed,
        config=config,
        candidate_ids=candidate_ids,
        labeled_ids=labeled_ids,
        oracle_calls=oracle_calls,
        total_cost=total_cost,
        pareto_front_ids=pareto_front_ids,
        surrogate_version=surrogate_version,
        data_version=data_version,
        status=status,
        notes=notes,
        provenance=Provenance(
            data_version=data_version,
            tool_versions={"peptideforge_loop": "0.1.0"},
        ),
    )

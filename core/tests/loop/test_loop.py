"""Tests for P9 closed-loop orchestrator."""

from __future__ import annotations

from pathlib import Path

import pytest

from peptideforge.contracts.models import LoopState
from peptideforge.loop import (
    ClosedLoopOrchestrator,
    LoopConfig,
    RayUnavailableError,
    load_loop_state,
    require_ray,
    run_public_sequence_space_validation,
    run_simulations_to_target_validation,
)
from peptideforge.loop.parallel import map_parallel


def test_serial_map_parallel() -> None:
    assert map_parallel(lambda x: x * 2, [1, 2, 3], use_ray=False) == [2, 4, 6]


def test_ray_required_fails_loud() -> None:
    try:
        require_ray()
    except RayUnavailableError:
        with pytest.raises(RayUnavailableError):
            map_parallel(lambda x: x, [1], use_ray=True)
    else:
        pytest.skip("Ray installed; missing-dep path skipped")


def test_simulation_campaign_persists_state(tmp_path: Path) -> None:
    cfg = LoopConfig(
        seed=0,
        simulation_mode=True,
        n_init=8,
        n_propose=20,
        batch_size=2,
        max_iterations=3,
        target_value=-4.0,
        acquisition="qnehvi",
        state_dir=str(tmp_path / "campaign"),
    )
    result = ClosedLoopOrchestrator(config=cfg).run()
    assert result.states
    assert (tmp_path / "campaign" / "latest.json").is_file()
    assert (tmp_path / "campaign" / "dataset.json").is_file()
    loaded = load_loop_state(tmp_path / "campaign" / "latest.json")
    assert isinstance(loaded, LoopState)
    assert loaded.seed == 0
    assert loaded.oracle_calls >= 8
    assert result.reports[0].numbers
    # Every reported number cites a tool
    for report in result.reports:
        for num in report.numbers:
            assert num.tool


@pytest.mark.eval
def test_loop_beats_random_simulations_to_target(tmp_path: Path) -> None:
    report = run_simulations_to_target_validation(
        seed=0,
        target_value=-4.5,
        n_pool=80,
        n_init=8,
        max_rounds=12,
        batch_size=2,
        state_dir=tmp_path / "spend",
    )
    assert report.passed, (
        f"qNEHVI calls={report.oracle_calls_qnehvi} best={report.best_qnehvi}; "
        f"random calls={report.oracle_calls_random} best={report.best_random}"
    )


@pytest.mark.eval
def test_loop_public_sequence_space() -> None:
    report = run_public_sequence_space_validation(
        seed=1,
        target_value=-3.5,
        n_pool=60,
        max_rounds=10,
    )
    assert report.passed, (
        f"public seq space: q={report.oracle_calls_qnehvi}/{report.best_qnehvi} "
        f"r={report.oracle_calls_random}/{report.best_random}"
    )
    assert report.mode.endswith("synthetic_oracle")
    assert report.best_qnehvi is not None
    assert report.best_random is not None


def test_reproducible_seed_config() -> None:
    cfg = LoopConfig(
        seed=42,
        n_init=8,
        max_iterations=2,
        batch_size=2,
        n_propose=16,
        target_value=-10.0,  # unreachable → full budget
        acquisition="random",
    )
    a = ClosedLoopOrchestrator(config=cfg).run()
    b = ClosedLoopOrchestrator(config=cfg).run()
    assert a.best_value == b.best_value
    assert a.states[-1].oracle_calls == b.states[-1].oracle_calls
    assert [r.candidate.sequence for r in a.dataset.records] == [
        r.candidate.sequence for r in b.dataset.records
    ]

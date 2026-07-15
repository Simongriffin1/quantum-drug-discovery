"""Closed-loop DBTL orchestrator — P9.

Simulation mode + spend-gate validation (simulations-to-target vs random).
"""

from peptideforge.loop.dataset import LabeledRecord, VersionedDataset
from peptideforge.loop.orchestrator import CampaignResult, ClosedLoopOrchestrator, LoopConfig
from peptideforge.loop.parallel import RayUnavailableError, map_parallel, require_ray
from peptideforge.loop.report import IterationReport, ToolAttributedNumber, build_iteration_report
from peptideforge.loop.simulation import SimulationOracle, SyntheticStructurePredictor
from peptideforge.loop.state import load_loop_state, save_loop_state, write_state_history
from peptideforge.loop.validate import (
    LoopValidationReport,
    run_pool_campaign,
    run_public_sequence_space_validation,
    run_simulations_to_target_validation,
)

__all__ = [
    "CampaignResult",
    "ClosedLoopOrchestrator",
    "IterationReport",
    "LabeledRecord",
    "LoopConfig",
    "LoopValidationReport",
    "RayUnavailableError",
    "SimulationOracle",
    "SyntheticStructurePredictor",
    "ToolAttributedNumber",
    "VersionedDataset",
    "build_iteration_report",
    "load_loop_state",
    "map_parallel",
    "require_ray",
    "run_pool_campaign",
    "run_public_sequence_space_validation",
    "run_simulations_to_target_validation",
    "save_loop_state",
    "write_state_history",
]

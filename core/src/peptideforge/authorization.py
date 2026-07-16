"""Machine-readable campaign authorization (Step 4B).

Binds validation runs to runtime scope gates. Cross-target is never authorized.
Predicted-fold within-target requires a measured, pre-registered PASS.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from enum import Enum
from pathlib import Path
from typing import Any
from uuid import UUID, uuid4

from pydantic import BaseModel, ConfigDict, Field


class TaskType(str, Enum):
    WITHIN_TARGET = "within_target"
    CROSS_TARGET = "cross_target"


class InputType(str, Enum):
    EXPERIMENTAL = "experimental"
    PREDICTED = "predicted"
    SIMULATION = "simulation"


class AuthorizationRecord(BaseModel):
    """Versioned license for a (task_type, input_type) campaign class."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    record_id: UUID = Field(default_factory=uuid4)
    task_type: TaskType
    input_type: InputType
    structure_source: str
    validated_rho: float | None = None
    ci_low: float | None = None
    ci_high: float | None = None
    n: int | None = None
    split_id: str | None = None
    skempi_run_id: str | None = None
    fold_confidence_metric: str | None = None
    fold_confidence_threshold: float | None = None
    authorized: bool
    reason: str
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    git_sha: str | None = None
    pre_registration_ref: str = "ACCEPTANCE.md#predicted-fold-within-target"

    def to_dict(self) -> dict[str, Any]:
        payload = self.model_dump(mode="json")
        payload["record_id"] = str(self.record_id)
        payload["created_at"] = self.created_at.isoformat()
        return payload

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> AuthorizationRecord:
        return cls.model_validate(data)


class AuthorizationDenied(PermissionError):
    """Raised when a campaign is outside the validated authorization boundary."""


def _git_sha() -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


def default_authorization_dir() -> Path:
    return (
        Path(__file__).resolve().parents[3]
        / "benchmarks"
        / "authorization"
    )


def build_authorization_bundle(
    *,
    experimental_skempi: dict[str, Any],
    predicted_degradation: dict[str, Any] | None,
    split_id: str,
) -> list[AuthorizationRecord]:
    """Populate records from measured SKEMPI (experimental) + degradation (predicted)."""
    sha = _git_sha()
    records: list[AuthorizationRecord] = []

    exp_ok = bool(experimental_skempi.get("gate_pass"))
    records.append(
        AuthorizationRecord(
            task_type=TaskType.WITHIN_TARGET,
            input_type=InputType.EXPERIMENTAL,
            structure_source="experimental_crystal_WT",
            validated_rho=experimental_skempi.get("spearman"),
            ci_low=experimental_skempi.get("spearman_ci_low"),
            ci_high=experimental_skempi.get("spearman_ci_high"),
            n=experimental_skempi.get("n"),
            split_id=split_id,
            skempi_run_id=experimental_skempi.get("artifact")
            or "skempi_ddg_powered_last_run.json",
            authorized=exp_ok,
            reason=(
                "Within-target experimental SKEMPI gate PASS"
                if exp_ok
                else "Within-target experimental SKEMPI gate not passed"
            ),
            git_sha=sha,
        )
    )

    # Cross-target always blocked
    records.append(
        AuthorizationRecord(
            task_type=TaskType.CROSS_TARGET,
            input_type=InputType.EXPERIMENTAL,
            structure_source="any",
            authorized=False,
            reason="Cross-target absolute affinity is not authorized (ACCEPTANCE.md)",
            git_sha=sha,
            split_id=split_id,
        )
    )
    records.append(
        AuthorizationRecord(
            task_type=TaskType.CROSS_TARGET,
            input_type=InputType.PREDICTED,
            structure_source="boltz2",
            authorized=False,
            reason="Cross-target absolute affinity is not authorized (ACCEPTANCE.md)",
            git_sha=sha,
            split_id=split_id,
        )
    )

    # Predicted within-target per pre-registered 4A.5 rule
    if predicted_degradation is None:
        records.append(
            AuthorizationRecord(
                task_type=TaskType.WITHIN_TARGET,
                input_type=InputType.PREDICTED,
                structure_source="boltz2",
                authorized=False,
                reason=(
                    "Predicted-fold degradation not measured (Boltz-2 unavailable or "
                    "run incomplete). Live predicted-fold campaigns remain BLOCKED "
                    "per ACCEPTANCE.md Step 4A.5 pre-registration."
                ),
                git_sha=sha,
                split_id=split_id,
            )
        )
        return records

    mode_a = predicted_degradation.get("mode_a") or {}
    stratified = predicted_degradation.get("stratified_authorization") or {}
    unconditional = bool(mode_a.get("gate_pass")) and int(mode_a.get("n") or 0) >= 30
    stratified_ok = bool(stratified.get("gate_pass"))
    thr = stratified.get("threshold")
    metric = stratified.get("metric")

    if unconditional:
        records.append(
            AuthorizationRecord(
                task_type=TaskType.WITHIN_TARGET,
                input_type=InputType.PREDICTED,
                structure_source="boltz2_mutate_in_place",
                validated_rho=mode_a.get("spearman"),
                ci_low=mode_a.get("spearman_ci_low"),
                ci_high=mode_a.get("spearman_ci_high"),
                n=mode_a.get("n"),
                split_id=split_id,
                skempi_run_id=predicted_degradation.get("artifact"),
                authorized=True,
                reason="Predicted-fold MODE A unconditional gate PASS (4A.5)",
                git_sha=sha,
            )
        )
    elif stratified_ok and thr is not None:
        records.append(
            AuthorizationRecord(
                task_type=TaskType.WITHIN_TARGET,
                input_type=InputType.PREDICTED,
                structure_source="boltz2_mutate_in_place",
                validated_rho=stratified.get("spearman"),
                ci_low=stratified.get("spearman_ci_low"),
                ci_high=stratified.get("spearman_ci_high"),
                n=stratified.get("n"),
                split_id=split_id,
                skempi_run_id=predicted_degradation.get("artifact"),
                fold_confidence_metric=str(metric) if metric else None,
                fold_confidence_threshold=float(thr),
                authorized=True,
                reason=(
                    f"Predicted-fold MODE A PASS restricted to {metric}>={thr} "
                    f"(N={stratified.get('n')}) per 4A.5"
                ),
                git_sha=sha,
            )
        )
    else:
        records.append(
            AuthorizationRecord(
                task_type=TaskType.WITHIN_TARGET,
                input_type=InputType.PREDICTED,
                structure_source="boltz2_mutate_in_place",
                validated_rho=mode_a.get("spearman"),
                ci_low=mode_a.get("spearman_ci_low"),
                ci_high=mode_a.get("spearman_ci_high"),
                n=mode_a.get("n"),
                split_id=split_id,
                skempi_run_id=predicted_degradation.get("artifact"),
                authorized=False,
                reason=(
                    predicted_degradation.get("block_reason")
                    or "Predicted-fold within-target gate FAIL under 4A.5 pre-registration"
                ),
                git_sha=sha,
            )
        )
    return records


def write_authorization_bundle(
    records: list[AuthorizationRecord],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "records": [r.to_dict() for r in records],
        "written_at": datetime.now(UTC).isoformat(),
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def load_authorization_bundle(path: Path) -> list[AuthorizationRecord]:
    if not path.is_file():
        raise FileNotFoundError(
            f"Authorization bundle missing: {path}. "
            "Run Step 4B build_authorization before non-simulation campaigns."
        )
    data = json.loads(path.read_text(encoding="utf-8"))
    return [AuthorizationRecord.from_dict(r) for r in data["records"]]


def find_authorization(
    records: list[AuthorizationRecord],
    *,
    task_type: TaskType,
    input_type: InputType,
) -> AuthorizationRecord | None:
    for rec in records:
        if rec.task_type == task_type and rec.input_type == input_type:
            return rec
    return None


def assert_campaign_authorized(
    records: list[AuthorizationRecord],
    *,
    task_type: TaskType,
    input_type: InputType,
    simulation_mode: bool = False,
    fold_confidence: float | None = None,
) -> AuthorizationRecord:
    """Enforce the authorization boundary; raise AuthorizationDenied if blocked.

    Simulation mode is always allowed (plumbing only — synthetic_*).
    """
    if simulation_mode or input_type == InputType.SIMULATION:
        return AuthorizationRecord(
            task_type=task_type,
            input_type=InputType.SIMULATION,
            structure_source="synthetic",
            authorized=True,
            reason="simulation_mode — synthetic plumbing only; not a physics claim",
            git_sha=_git_sha(),
        )

    rec = find_authorization(records, task_type=task_type, input_type=input_type)
    if rec is None:
        raise AuthorizationDenied(
            f"No authorization record for task_type={task_type.value} "
            f"input_type={input_type.value}. Refusing campaign."
        )
    if not rec.authorized:
        raise AuthorizationDenied(
            f"Campaign BLOCKED by authorization record {rec.record_id}: {rec.reason} "
            f"(rho={rec.validated_rho}, CI=[{rec.ci_low},{rec.ci_high}], n={rec.n}, "
            f"skempi_run_id={rec.skempi_run_id}). "
            "See ACCEPTANCE.md predicted-fold rule."
        )
    if (
        rec.fold_confidence_threshold is not None
        and fold_confidence is not None
        and fold_confidence < rec.fold_confidence_threshold
    ):
        raise AuthorizationDenied(
            f"Fold confidence {fold_confidence:.3f} < authorized threshold "
            f"{rec.fold_confidence_threshold:.3f} "
            f"(metric={rec.fold_confidence_metric}) for record {rec.record_id}."
        )
    return rec


def filter_by_fold_confidence(
    confidences: dict[str, float],
    *,
    threshold: float,
) -> tuple[list[str], list[str]]:
    """Return (kept_ids, excluded_ids) for confidence-scoped authorization."""
    kept = [i for i, c in confidences.items() if c >= threshold]
    excluded = [i for i, c in confidences.items() if c < threshold]
    return kept, excluded

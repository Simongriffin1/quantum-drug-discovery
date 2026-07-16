"""Collect benchmark numbers from local JSON artifacts and/or live seeded runs.

Never fabricates metrics. Missing physics/affinity runs are reported as NOT_RUN /
FAIL with an honest pointer to the artifact path.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from peptideforge.contracts.models import Candidates, PeptideCandidate
from peptideforge.generators.mutation import MutationGenerator
from peptideforge.acquisition.validate import run_branin_currin_validation
from peptideforge.eval.harness import DEFAULT_THRESHOLDS
from peptideforge.loop.validate import run_simulations_to_target_validation
from peptideforge.reports.models import BenchmarkReport, ReportSection, TraceableNumber
from peptideforge.surrogate.ensemble import DeepEnsembleSurrogate, synthetic_physics_label
from peptideforge.surrogate.report import (
    SurrogateAcceptanceReport,
    make_oracle_labels,
    run_surrogate_acceptance,
    write_acceptance_report,
)


def _repo_root() -> Path:
    # core/src/peptideforge/reports -> repo root
    return Path(__file__).resolve().parents[4]


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
    if out.returncode != 0:
        return None
    return out.stdout.strip() or None


def default_oracle_artifact() -> Path:
    v3 = (
        _repo_root()
        / "benchmarks"
        / "peptide_affinity"
        / "data"
        / "oracle_validity_v3_oneshot_test.json"
    )
    if v3.is_file():
        return v3
    v2 = (
        _repo_root()
        / "benchmarks"
        / "peptide_affinity"
        / "data"
        / "oracle_validity_v2_last_run.json"
    )
    if v2.is_file():
        return v2
    return (
        _repo_root()
        / "benchmarks"
        / "fixtures"
        / "structures"
        / "oracle_validity_last_run.json"
    )


def load_json_artifact(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            f"benchmark artifact missing: {path}. "
            "Run the corresponding eval and write the JSON before reporting."
        )
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}")
    return data


def collect_oracle_validity(artifact: Path | None = None) -> ReportSection:
    path = artifact or default_oracle_artifact()
    if not path.is_file():
        return ReportSection(
            title="Oracle validity (affinity)",
            status="NOT_RUN",
            summary=f"No oracle-validity artifact at {path}",
            numbers=(),
            details={"artifact": str(path)},
        )
    data = load_json_artifact(path)
    spearman = float(data["spearman"])
    threshold = float(data.get("spearman_threshold", DEFAULT_THRESHOLDS.oracle_affinity_spearman))
    # Prefer explicit gate boolean from v2 artifact (CI + min_n + red-team aware)
    if "passed_threshold" in data:
        passed = bool(data["passed_threshold"])
    else:
        passed = spearman >= threshold
    n = int(data["n"])
    measurable = bool(data.get("measurable", n >= 30))
    ci_lo = data.get("spearman_ci_low")
    ci_hi = data.get("spearman_ci_high")
    mlflow = data.get("mlflow") or {}
    try:
        display_path = path.resolve().relative_to(_repo_root())
    except ValueError:
        display_path = path
    source = (
        f"file:{display_path}"
        + (f"; mlflow_run_id={mlflow['run_id']}" if mlflow.get("run_id") else "")
    )
    ci_txt = (
        f" CI95%=[{float(ci_lo):.3f},{float(ci_hi):.3f}]"
        if ci_lo is not None and ci_hi is not None
        else " (no CI — insufficient for gate)"
    )
    status = "PASS" if passed else ("FAIL" if measurable else "NOT_MEASURABLE")
    numbers = [
        TraceableNumber(
            name="spearman",
            value=spearman,
            source=source,
            data_version=data.get("subset_name"),
        ),
        TraceableNumber(
            name="rmse",
            value=float(data["rmse"]),
            source=source,
            data_version=data.get("subset_name"),
        ),
        TraceableNumber(
            name="n",
            value=float(n),
            unit="count",
            source=source,
            data_version=data.get("subset_name"),
        ),
        TraceableNumber(
            name="spearman_threshold",
            value=threshold,
            source="ACCEPTANCE.md",
        ),
        TraceableNumber(
            name="passed_threshold",
            value=1.0 if passed else 0.0,
            source=source,
        ),
    ]
    if ci_lo is not None and ci_hi is not None:
        numbers.extend(
            [
                TraceableNumber(
                    name="spearman_ci_low",
                    value=float(ci_lo),
                    source=source,
                    data_version=data.get("subset_name"),
                ),
                TraceableNumber(
                    name="spearman_ci_high",
                    value=float(ci_hi),
                    source=source,
                    data_version=data.get("subset_name"),
                ),
                TraceableNumber(
                    name="pearson",
                    value=float(data.get("pearson", float("nan"))),
                    source=source,
                    data_version=data.get("subset_name"),
                ),
            ]
        )
    return ReportSection(
        title="Oracle validity (affinity)",
        status=status,
        summary=(
            f"Subset `{data.get('subset_name', 'unknown')}` N={n} "
            f"partition={data.get('partition', 'n/a')}: "
            f"Spearman={spearman:.4f}{ci_txt} (threshold ≥ {threshold}, "
            f"require CI_low>0, N≥30, red-team). "
            f"Gate {'PASSED' if passed else 'FAILED'} — reported honestly."
        ),
        numbers=tuple(numbers),
        details={
            "subset_name": data.get("subset_name"),
            "protocol": data.get("protocol"),
            "mlflow": mlflow,
            "red_team": data.get("red_team"),
            "measurable": measurable,
            "record_ids": [d.get("record_id") for d in data.get("details", [])],
            "artifact": str(path),
        },
    )


def collect_ddg_status(artifact: Path | None = None) -> ReportSection:
    """SKEMPI within-target ΔΔG co-primary gate (ACCEPTANCE.md)."""
    path = artifact or (
        _repo_root() / "benchmarks" / "skempi" / "data" / "skempi_ddg_last_run.json"
    )
    if not path.is_file():
        return ReportSection(
            title="Stability / ddG (SKEMPI within-target)",
            status="NOT_RUN",
            summary=(
                "No logged SKEMPI ΔΔG run. Threshold pre-registered in ACCEPTANCE.md "
                "(Spearman ≥ 0.30, CI_low > 0). Refusing to invent a correlation."
            ),
            numbers=(
                TraceableNumber(
                    name="ddg_spearman_threshold",
                    value=DEFAULT_THRESHOLDS.oracle_ddg_spearman,
                    source="ACCEPTANCE.md",
                ),
            ),
            details={"fixture": "benchmarks/fixtures/skempi_ddg_v1.tsv"},
        )

    data = load_json_artifact(path)
    rho = data.get("spearman")
    lo = data.get("spearman_ci_low")
    hi = data.get("spearman_ci_high")
    n = data.get("n")
    gate_pass = bool(data.get("gate_pass"))
    status = "PASS" if gate_pass else "FAIL"
    numbers = [
        TraceableNumber(
            name="skempi_ddg_spearman",
            value=float(rho) if rho is not None else None,
            source=f"file:{path}",
            notes=f"N={n}; CI=[{lo}, {hi}]",
        ),
        TraceableNumber(
            name="skempi_ddg_spearman_ci_low",
            value=float(lo) if lo is not None else None,
            source=f"file:{path}",
        ),
        TraceableNumber(
            name="skempi_ddg_gate_threshold",
            value=float(data.get("gate_threshold") or 0.30),
            source="ACCEPTANCE.md",
        ),
    ]
    return ReportSection(
        title="Stability / ddG (SKEMPI within-target)",
        status=status,
        summary=(
            f"Within-target held-out SKEMPI ΔΔG: N={n}, Spearman ρ={rho} "
            f"(95% CI [{lo}, {hi}]). Pre-registered gate ρ≥0.30 with CI_low>0: "
            f"{'PASSED' if gate_pass else 'FAILED'}. "
            f"Red-team={'pass' if (data.get('red_team') or {}).get('passed') else 'fail'}. "
            "Cross-target affinity gate is reported separately and is not replaced."
        ),
        numbers=tuple(numbers),
        details={
            "artifact": str(path),
            "protocol": data.get("protocol"),
            "red_team": data.get("red_team"),
            "gate_pass": gate_pass,
        },
    )


def _synthetic_pool(n: int, *, seed: int) -> list[PeptideCandidate]:
    gen = MutationGenerator(generation_method="synthetic_report_pool")
    seeds = (
        "FLIVVFLIV",
        "DDEEGGSSS",
        "AAAAKKKKK",
        "WWYYLLIIV",
        "GILGFVFTL",
        "LLFGYPVYV",
    )
    batch = gen.propose(n=n, seed_sequences=seeds, seed=seed)
    return list(batch.items)


def collect_surrogate_calibration(
    artifacts_dir: Path,
    *,
    seed: int = 0,
) -> ReportSection:
    """Run (or reload) synthetic_* surrogate acceptance — plumbing/UQ only."""
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    out_path = artifacts_dir / "surrogate_acceptance_last_run.json"
    if out_path.is_file():
        data = load_json_artifact(out_path)
        report = SurrogateAcceptanceReport.model_validate(data)
        source = f"file:{out_path}"
    else:
        train = _synthetic_pool(40, seed=seed)
        test = _synthetic_pool(20, seed=seed + 1)
        train_seqs = {c.sequence for c in train}
        test = [c for c in test if c.sequence not in train_seqs][:16]
        train_batch = Candidates(items=tuple(train), seed=seed)
        test_batch = Candidates(items=tuple(test), seed=seed + 1)
        train_vals = {
            c.candidate_id: synthetic_physics_label(c.sequence, noise=0.05, seed=11)
            for c in train
        }
        test_vals = {
            c.candidate_id: synthetic_physics_label(c.sequence, noise=0.05, seed=11)
            for c in test
        }
        labels = make_oracle_labels(train_batch, train_vals)
        surrogate = DeepEnsembleSurrogate(
            n_ensemble=6,
            l2=0.5,
            coverage_target=0.90,
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
            train_clusters={str(c.candidate_id): "cluster_train" for c in train},
            test_clusters={str(c.candidate_id): "cluster_test" for c in test},
            train_sequences={str(c.candidate_id): c.sequence for c in train},
            test_sequences={str(c.candidate_id): c.sequence for c in test},
            seed=seed,
            thresholds=DEFAULT_THRESHOLDS,
            data_version="synthetic_surrogate_v1",
            notes="synthetic_learnable_oracle_not_physics",
        )
        write_acceptance_report(report, out_path)
        source = f"live:run_surrogate_acceptance;file:{out_path}"

    cal = report.calibration
    return ReportSection(
        title="Surrogate calibration (synthetic plumbing)",
        status="PASS" if report.accepted else "FAIL",
        summary=(
            f"ECE={cal.ece:.4f} (threshold < {cal.ece_threshold}); "
            f"red-team passed={report.red_team.passed}. "
            "This validates UQ plumbing on synthetic_* labels — NOT physics affinity."
        ),
        numbers=(
            TraceableNumber(
                name="ece",
                value=cal.ece,
                source=source,
                data_version=report.data_version,
            ),
            TraceableNumber(
                name="empirical_coverage",
                value=cal.empirical_coverage,
                source=source,
                data_version=report.data_version,
            ),
            TraceableNumber(
                name="ece_threshold",
                value=cal.ece_threshold,
                source="ACCEPTANCE.md",
            ),
            TraceableNumber(
                name="red_team_passed",
                value=1.0 if report.red_team.passed else 0.0,
                source=source,
            ),
            TraceableNumber(
                name="model_rho",
                value=report.red_team.model_rho,
                source=source,
                data_version=report.data_version,
            ),
        ),
        details={
            "data_version": report.data_version,
            "notes": report.notes,
            "artifact": str(out_path),
        },
    )


def collect_loop_efficiency(
    artifacts_dir: Path,
    *,
    seed: int = 0,
) -> ReportSection:
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    out_path = artifacts_dir / "loop_spend_gate_last_run.json"
    if out_path.is_file():
        data = load_json_artifact(out_path)
        from peptideforge.loop.validate import LoopValidationReport

        report = LoopValidationReport.model_validate(data)
        source = f"file:{out_path}"
    else:
        report = run_simulations_to_target_validation(
            seed=seed,
            target_value=-4.5,
            n_pool=60,
            n_init=8,
            max_rounds=10,
            batch_size=2,
            state_dir=artifacts_dir / "loop_runs",
        )
        out_path.write_text(report.model_dump_json(indent=2), encoding="utf-8")
        source = f"live:run_simulations_to_target_validation;file:{out_path}"

    return ReportSection(
        title="Loop efficiency (simulations-to-target)",
        status="PASS" if report.passed else "FAIL",
        summary=(
            f"qNEHVI calls={report.oracle_calls_qnehvi} best={report.best_qnehvi}; "
            f"random calls={report.oracle_calls_random} best={report.best_random}; "
            f"target={report.target_value}. Mode={report.mode}."
        ),
        numbers=(
            TraceableNumber(
                name="oracle_calls_qnehvi",
                value=float(report.oracle_calls_qnehvi)
                if report.oracle_calls_qnehvi is not None
                else None,
                unit="count",
                source=source,
                data_version=report.mode,
            ),
            TraceableNumber(
                name="oracle_calls_random",
                value=float(report.oracle_calls_random)
                if report.oracle_calls_random is not None
                else None,
                unit="count",
                source=source,
                data_version=report.mode,
            ),
            TraceableNumber(
                name="best_qnehvi",
                value=report.best_qnehvi,
                source=source,
                data_version=report.mode,
            ),
            TraceableNumber(
                name="best_random",
                value=report.best_random,
                source=source,
                data_version=report.mode,
            ),
            TraceableNumber(
                name="target_value",
                value=report.target_value,
                source=source,
            ),
            TraceableNumber(
                name="spend_gate_passed",
                value=1.0 if report.passed else 0.0,
                source=source,
            ),
        ),
        details={"notes": report.notes, "artifact": str(out_path)},
    )


def collect_acquisition(
    artifacts_dir: Path,
    *,
    seed: int = 7,
) -> ReportSection:
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    out_path = artifacts_dir / "acquisition_branin_currin_last_run.json"
    if out_path.is_file():
        data = load_json_artifact(out_path)
        from peptideforge.acquisition.validate import AcquisitionValidationReport

        report = AcquisitionValidationReport.model_validate(data)
        source = f"file:{out_path}"
    else:
        report = run_branin_currin_validation(
            n_pool=60,
            n_init=6,
            n_rounds=5,
            batch_size=2,
            seed=seed,
            epistemic_std=0.8,
        )
        out_path.write_text(report.model_dump_json(indent=2), encoding="utf-8")
        source = f"live:run_branin_currin_validation;file:{out_path}"

    return ReportSection(
        title="Acquisition (Branin–Currin hypervolume)",
        status="PASS" if report.passed else "FAIL",
        summary=(
            f"qNEHVI HV={report.hv_qnehvi:.4f} vs random HV={report.hv_random:.4f} "
            f"(same budget). Validates acquisition logic, not biology."
        ),
        numbers=(
            TraceableNumber(
                name="hv_qnehvi",
                value=report.hv_qnehvi,
                source=source,
                notes=report.notes,
            ),
            TraceableNumber(
                name="hv_random",
                value=report.hv_random,
                source=source,
                notes=report.notes,
            ),
            TraceableNumber(
                name="acquisition_passed",
                value=1.0 if report.passed else 0.0,
                source=source,
            ),
        ),
        details={
            "n_pool": report.n_pool,
            "n_init": report.n_init,
            "n_rounds": report.n_rounds,
            "batch_size": report.batch_size,
            "reference": list(report.reference),
            "artifact": str(out_path),
        },
    )


def try_mlflow_latest_oracle(
    tracking_uri: str,
    *,
    experiment: str = "peptideforge-oracle-validity",
) -> dict[str, Any] | None:
    """Fetch latest MLflow run metrics if mlflow is installed; else None."""
    try:
        import mlflow
        from mlflow.tracking import MlflowClient
    except ImportError:
        return None
    mlflow.set_tracking_uri(tracking_uri)
    client = MlflowClient()
    exp = client.get_experiment_by_name(experiment)
    if exp is None:
        return None
    runs = client.search_runs(
        experiment_ids=[exp.experiment_id],
        order_by=["attributes.start_time DESC"],
        max_results=1,
    )
    if not runs:
        return None
    run = runs[0]
    return {
        "run_id": run.info.run_id,
        "experiment_id": run.info.experiment_id,
        "tracking_uri": tracking_uri,
        "metrics": dict(run.data.metrics),
        "params": dict(run.data.params),
    }


def build_benchmark_report(
    *,
    artifacts_dir: Path | None = None,
    oracle_artifact: Path | None = None,
    mlflow_uri: str | None = None,
    seed: int = 0,
) -> BenchmarkReport:
    root = _repo_root()
    art = artifacts_dir or (root / "benchmarks" / "reports" / "artifacts")
    art.mkdir(parents=True, exist_ok=True)

    oracle_section = collect_oracle_validity(oracle_artifact)
    if mlflow_uri:
        latest = try_mlflow_latest_oracle(mlflow_uri)
        if latest is not None:
            details = dict(oracle_section.details)
            details["mlflow_latest"] = latest
            oracle_section = ReportSection(
                title=oracle_section.title,
                status=oracle_section.status,
                summary=oracle_section.summary
                + f" MLflow latest run_id={latest['run_id']}.",
                numbers=oracle_section.numbers,
                details=details,
            )

    sections = (
        oracle_section,
        collect_ddg_status(),
        collect_surrogate_calibration(art, seed=seed),
        collect_acquisition(art, seed=seed + 7),
        collect_loop_efficiency(art, seed=seed),
    )
    caveats = (
        "Oracle-validity FAIL means binding campaigns are not authorized "
        "(CURSOR_PROJECT_CONTEXT §8).",
        "Surrogate / loop / acquisition sections marked synthetic_* validate "
        "algorithmic plumbing only — not physics fidelity.",
        "All numbers below are traced to JSON artifacts and/or MLflow run IDs; "
        "none are invented.",
    )
    return BenchmarkReport(
        generated_at=datetime.now(UTC).isoformat(),
        git_sha=_git_sha(),
        sections=sections,
        caveats=caveats,
    )

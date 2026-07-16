"""Extended affinity oracle-validity with bootstrap CIs and red-team battery.

Validates MM-GBSA on **experimental** prepared structures (oracle error isolated
from folding error). Gate (ACCEPTANCE.md):

  Spearman ρ ≥ 0.40 on held-out test (N≥30) with lower CI bound > 0,
  AND all red-team controls pass.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field

from peptideforge.contracts.models import ComplexStructure, OracleTier, PeptideCandidate
from peptideforge.eval.harness import PredictionLabelPair
from peptideforge.eval.metrics import bootstrap_ci, pearson_r, rmse, spearman_rho
from peptideforge.eval.redteam import RedTeamReport, run_red_team
from peptideforge.oracles.mm_gbsa import OpenMMOracleConfig, OpenMMPhysicsOracle

SUBSET_NAME_V2 = "peptide_affinity_v2_experimental_openmm"


class AffinityValidityReport(BaseModel):
    """Honest oracle-validity with CIs — never declare PASS on a point estimate alone."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    n: int
    spearman: float
    spearman_ci_low: float
    spearman_ci_high: float
    pearson: float
    pearson_ci_low: float
    pearson_ci_high: float
    rmse: float
    spearman_threshold: float
    min_n: int
    measurable: bool
    passed: bool
    pairs: tuple[PredictionLabelPair, ...]
    red_team: RedTeamReport | None = None
    notes: str | None = None


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


def evaluate_affinity_with_ci(
    pairs: list[PredictionLabelPair],
    *,
    thresholds: ExtendedThresholds | None = None,
    min_n: int = 30,
    n_bootstrap: int = 1000,
    seed: int = 0,
    red_team: RedTeamReport | None = None,
) -> AffinityValidityReport:
    thr = thresholds or EXTENDED_THRESHOLDS
    if len(pairs) < 3:
        raise ValueError("evaluate_affinity_with_ci needs ≥ 3 pairs")
    pred = [p.predicted for p in pairs]
    exp = [p.experimental for p in pairs]
    rho, rho_lo, rho_hi = bootstrap_ci(
        pred, exp, statistic="spearman", n_resamples=n_bootstrap, seed=seed
    )
    r, r_lo, r_hi = bootstrap_ci(
        pred, exp, statistic="pearson", n_resamples=n_bootstrap, seed=seed + 1
    )
    err = rmse(pred, exp)
    measurable = len(pairs) >= min_n
    # Gate: measurable AND ρ≥threshold AND lower CI > 0 AND red-team (if provided)
    passed = (
        measurable
        and rho >= thr.oracle_affinity_spearman
        and rho_lo > 0.0
        and (red_team is None or red_team.passed)
    )
    notes = None
    if not measurable:
        notes = f"N={len(pairs)} < min_n={min_n}: point estimate not sufficient for PASS/FAIL"
    return AffinityValidityReport(
        n=len(pairs),
        spearman=rho,
        spearman_ci_low=rho_lo,
        spearman_ci_high=rho_hi,
        pearson=r,
        pearson_ci_low=r_lo,
        pearson_ci_high=r_hi,
        rmse=err,
        spearman_threshold=thr.oracle_affinity_spearman,
        min_n=min_n,
        measurable=measurable,
        passed=passed,
        pairs=tuple(pairs),
        red_team=red_team,
        notes=notes,
    )


def trivial_baselines_for_structures(
    pairs: list[PredictionLabelPair],
    *,
    peptide_lengths: dict[str, int],
    net_charges: dict[str, float] | None = None,
    interface_contacts: dict[str, float] | None = None,
    buried_sasa: dict[str, float] | None = None,
) -> dict[str, list[float]]:
    """Naive predictors the oracle must beat (fail loud if a baseline wins)."""
    baselines: dict[str, list[float]] = {
        "peptide_length": [float(peptide_lengths[p.record_id]) for p in pairs],
    }
    if net_charges is not None:
        baselines["net_charge"] = [float(net_charges[p.record_id]) for p in pairs]
    if interface_contacts is not None:
        baselines["interface_contacts"] = [
            float(interface_contacts[p.record_id]) for p in pairs
        ]
    if buried_sasa is not None:
        baselines["buried_sasa"] = [float(buried_sasa[p.record_id]) for p in pairs]
    return baselines


def charge_from_sequence(seq: str) -> float:
    return float(sum(1 for c in seq.upper() if c in "KR") - sum(1 for c in seq.upper() if c in "DE"))


def run_trivial_baseline_battery(
    pairs: list[PredictionLabelPair],
    baselines: dict[str, list[float]],
    *,
    min_delta_rho: float = 0.05,
) -> tuple[bool, dict[str, float], str | None]:
    """Oracle must beat ALL trivial baselines; otherwise halt and report."""
    model_rho = spearman_rho(
        [p.predicted for p in pairs], [p.experimental for p in pairs]
    )
    baseline_rhos: dict[str, float] = {}
    exp = [p.experimental for p in pairs]
    for name, preds in baselines.items():
        try:
            baseline_rhos[name] = spearman_rho(preds, exp)
        except ValueError:
            baseline_rhos[name] = 0.0
    winners = [
        name
        for name, br in baseline_rhos.items()
        if (model_rho - br) < min_delta_rho
    ]
    if winners:
        return (
            False,
            baseline_rhos,
            f"trivial baseline(s) win or tie within Δρ: {winners}; "
            f"model_rho={model_rho:.4f} baselines={baseline_rhos}",
        )
    return True, baseline_rhos, None


def score_prepared_structures(
    manifest_rows: list[dict[str, str]],
    *,
    pk_by_id: dict[str, float],
    seq_by_id: dict[str, str],
    minimize_max_iterations: int = 0,
    seed: int = 0,
    platform: str | None = "CPU",
    forcefield_xml: tuple[str, ...] = ("amber14-all.xml", "implicit/gbn2.xml"),
    solute_dielectric: float = 1.0,
    salt_conc_M: float = 0.0,
) -> tuple[list[PredictionLabelPair], list[dict[str, Any]], list[dict[str, str]]]:
    """Run OpenMM MM-GBSA on prepared experimental complexes."""
    pairs: list[PredictionLabelPair] = []
    details: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    for row in manifest_rows:
        rid = row["record_id"]
        pdb_file = Path(row["pdb_path"])
        if not pdb_file.is_file():
            failures.append({"record_id": rid, "error": f"missing {pdb_file}"})
            continue
        try:
            experimental_pk = float(row.get("experimental_pk") or pk_by_id[rid])
            peptide_seq = (row.get("peptide_sequence") or seq_by_id[rid]).upper()
        except KeyError as exc:
            failures.append({"record_id": rid, "error": str(exc)})
            continue
        peptide_chain = row.get("peptide_chain") or "C"
        cand = PeptideCandidate(
            candidate_id=uuid4(),
            sequence=peptide_seq[:50] if len(peptide_seq) >= 5 else (peptide_seq + "AAAAA")[:5],
            generation_method="benchmark_experimental",
            metadata={"record_id": rid, "peptide_chain": peptide_chain},
        )
        complex_structure = ComplexStructure(
            candidate_id=cand.candidate_id,
            target_id=row.get("pdb_id") or rid,
            sequence=cand.sequence,
            pdb_path=str(pdb_file.resolve()),
            confidence=1.0,
            fold_method="experimental_prepared",
        )
        try:
            oracle = OpenMMPhysicsOracle(
                OpenMMOracleConfig(
                    forcefield_xml=forcefield_xml,
                    minimize_max_iterations=minimize_max_iterations,
                    docking_minimize_iterations=0,
                    md_steps=200,
                    seed=seed,
                    platform=platform,
                    peptide_chain_ids=(peptide_chain, "P", "L", "C"),
                    solute_dielectric=solute_dielectric,
                    salt_conc_M=salt_conc_M,
                )
            )
            result = oracle.evaluate(complex_structure, tier=OracleTier.MM_GBSA)
        except Exception as exc:  # noqa: BLE001
            failures.append({"record_id": rid, "error": str(exc)})
            continue
        predicted = -float(result.value)
        pairs.append(
            PredictionLabelPair(
                record_id=rid,
                predicted=predicted,
                experimental=experimental_pk,
                unit="pK_proxy_neg_dG",
            )
        )
        details.append(
            {
                "record_id": rid,
                "pdb_path": str(pdb_file),
                "experimental_pk": experimental_pk,
                "mm_gbsa_kcal_mol": result.value,
                "predicted_neg_dG": predicted,
            }
        )
    return pairs, details, failures


class ExtendedThresholds(BaseModel):
    """Extended M2 affinity gate (ACCEPTANCE.md)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    oracle_affinity_spearman: float = 0.40
    oracle_ddg_spearman: float = 0.30
    surrogate_ece: float = 0.10
    label_shuffle_max_abs_rho: float = 0.20
    label_shuffle_min_drop: float = 0.40
    trivial_baseline_min_delta_rho: float = 0.05
    min_n: int = Field(30, description="Minimum test N for measurable gate")
    require_ci_lower_positive: bool = True


EXTENDED_THRESHOLDS = ExtendedThresholds()


def report_to_dict(report: AffinityValidityReport, extras: dict[str, Any]) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "subset_name": extras.get("subset_name", SUBSET_NAME_V2),
        "n": report.n,
        "spearman": report.spearman,
        "spearman_ci_low": report.spearman_ci_low,
        "spearman_ci_high": report.spearman_ci_high,
        "pearson": report.pearson,
        "pearson_ci_low": report.pearson_ci_low,
        "pearson_ci_high": report.pearson_ci_high,
        "rmse": report.rmse,
        "spearman_threshold": report.spearman_threshold,
        "min_n": report.min_n,
        "measurable": report.measurable,
        "passed_threshold": report.passed,
        "git_sha": _git_sha(),
        "data_version": extras.get("data_version"),
        "notes": report.notes,
        "pairs": [p.model_dump() for p in report.pairs],
    }
    if report.red_team is not None:
        payload["red_team"] = report.red_team.model_dump()
    payload.update({k: v for k, v in extras.items() if k not in payload})
    return payload


def write_validity_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

"""Deep-ensemble ridge surrogate + split conformal UQ (P7).

Epistemic uncertainty: bootstrap ensemble of ridge regressors (Lakatos et al.-
style deep ensemble idea without requiring PyTorch). Aleatoric/calibration:
inductive conformal prediction intervals around the ensemble mean.

Evaluate ONLY via the eval harness with homology-aware splits. The calibration
gate (ECE < 0.10) is declared in ACCEPTANCE.md — report honestly.
"""

from __future__ import annotations

import random
from pathlib import Path
from uuid import UUID

from peptideforge.contracts.models import (
    CalibratedPrediction,
    Candidates,
    ObjectiveVector,
    OracleResult,
    PeptideCandidate,
)
from peptideforge.surrogate.calibration import CalibrationReport, evaluate_calibration
from peptideforge.surrogate.conformal import absolute_residuals, conformal_quantile
from peptideforge.surrogate.features import candidate_features, sequence_features
from peptideforge.surrogate.linalg import fit_ridge, predict_linear


class SurrogateNotFittedError(RuntimeError):
    """Raised when predict is called before fit."""


class DeepEnsembleSurrogate:
    """Bootstrap ridge ensemble + conformal calibrated intervals.

    Sequence features must be registered (via ``register`` or ``Candidates`` seen
    at predict time) before ``fit`` so labels can be joined by ``candidate_id``.
    """

    def __init__(
        self,
        *,
        n_ensemble: int = 8,
        l2: float = 1.0,
        coverage_target: float = 0.90,
        calib_fraction: float = 0.25,
        objective_name: str = "binding",
        higher_is_better: bool = False,
    ) -> None:
        if n_ensemble < 2:
            raise ValueError("n_ensemble must be ≥ 2 for epistemic std")
        if not 0.0 < calib_fraction < 0.5:
            raise ValueError("calib_fraction must be in (0, 0.5)")
        if not 0.0 < coverage_target < 1.0:
            raise ValueError("coverage_target must be in (0, 1)")
        self.n_ensemble = n_ensemble
        self.l2 = l2
        self.coverage_target = coverage_target
        self.calib_fraction = calib_fraction
        self.objective_name = objective_name
        self.higher_is_better = higher_is_better
        self._feature_store: dict[UUID, tuple[float, ...]] = {}
        self._sequence_store: dict[UUID, str] = {}
        self._members: list[list[float]] = []
        self._conformal_radius: float | None = None
        self._fitted = False
        self.last_calibration: CalibrationReport | None = None

    @property
    def is_fitted(self) -> bool:
        return self._fitted

    def register(self, candidates: Candidates | PeptideCandidate) -> None:
        """Index sequence features by candidate_id for later fit."""
        items: tuple[PeptideCandidate, ...]
        if isinstance(candidates, PeptideCandidate):
            items = (candidates,)
        else:
            items = candidates.items
        for cand in items:
            feats = candidate_features(cand)
            self._feature_store[cand.candidate_id] = feats
            self._sequence_store[cand.candidate_id] = cand.sequence

    def fit(
        self,
        candidate_ids: tuple[UUID, ...],
        labels: tuple[OracleResult, ...],
        *,
        seed: int | None = None,
    ) -> None:
        """Fit ensemble on physics-oracle labels; calibrate conformal radius."""
        if len(candidate_ids) != len(labels):
            raise ValueError("candidate_ids and labels must align")
        if len(candidate_ids) < 6:
            raise ValueError("fit requires ≥ 6 labeled examples (train + calib)")
        for cid, lab in zip(candidate_ids, labels, strict=True):
            if lab.candidate_id != cid:
                raise ValueError(
                    f"label candidate_id {lab.candidate_id} != fit key {cid}"
                )
            if cid not in self._feature_store:
                raise KeyError(
                    f"no features for candidate_id={cid}; call register() first"
                )

        rng = random.Random(seed)
        pairs = list(zip(candidate_ids, labels, strict=True))
        rng.shuffle(pairs)
        n_calib = max(2, int(round(len(pairs) * self.calib_fraction)))
        n_train = len(pairs) - n_calib
        if n_train < 3:
            raise ValueError("insufficient train size after conformal split")
        train = pairs[:n_train]
        calib = pairs[n_train:]

        x_train = [list(self._feature_store[cid]) for cid, _ in train]
        # Add bias column
        x_train_b = [[1.0, *row] for row in x_train]
        y_train = [lab.value for _, lab in train]

        self._members = []
        for m in range(self.n_ensemble):
            # Bootstrap resample
            idxs = [rng.randrange(n_train) for _ in range(n_train)]
            xb = [x_train_b[i] for i in idxs]
            yb = [y_train[i] for i in idxs]
            weights = fit_ridge(xb, yb, l2=self.l2)
            self._members.append(weights)

        # Calibration residuals on ensemble mean
        calib_preds: list[float] = []
        calib_obs: list[float] = []
        for cid, lab in calib:
            mean, _ = self._ensemble_predict(self._feature_store[cid])
            calib_preds.append(mean)
            calib_obs.append(lab.value)
        scores = absolute_residuals(calib_preds, calib_obs)
        self._conformal_radius = conformal_quantile(scores, coverage=self.coverage_target)
        self._fitted = True

        # In-sample calibration report on calib fold (honest, not test metrics)
        lowers = [p - self._conformal_radius for p in calib_preds]
        uppers = [p + self._conformal_radius for p in calib_preds]
        self.last_calibration = evaluate_calibration(
            lowers,
            uppers,
            calib_obs,
            coverage_target=self.coverage_target,
            notes="conformal_calib_fold",
        )

    def predict(self, candidates: Candidates) -> tuple[ObjectiveVector, ...]:
        """Predict mean + conformal interval per candidate."""
        if not self._fitted or self._conformal_radius is None:
            raise SurrogateNotFittedError("call fit() before predict()")
        self.register(candidates)
        radius = self._conformal_radius
        out: list[ObjectiveVector] = []
        for cand in candidates.items:
            feats = self._feature_store[cand.candidate_id]
            mean, epi_std = self._ensemble_predict(feats)
            # Pure conformal half-width (finite-sample coverage guarantee).
            # Epistemic std is reported separately — do not inflate the interval
            # or ECE becomes a different quantity than the calibrated coverage.
            half = radius
            out.append(
                ObjectiveVector(
                    candidate_id=cand.candidate_id,
                    predictions=(
                        CalibratedPrediction(
                            candidate_id=cand.candidate_id,
                            objective_name=self.objective_name,
                            mean=mean,
                            lower=mean - half,
                            upper=mean + half,
                            epistemic_std=epi_std,
                            coverage_target=self.coverage_target,
                        ),
                    ),
                )
            )
        return tuple(out)

    def _ensemble_predict(self, feats: tuple[float, ...]) -> tuple[float, float]:
        x = [1.0, *feats]
        preds = [predict_linear(w, x) for w in self._members]
        mean = sum(preds) / len(preds)
        var = sum((p - mean) ** 2 for p in preds) / max(1, len(preds) - 1)
        return mean, var**0.5

    def write_calibration_report(self, path: Path) -> CalibrationReport:
        """Persist last calibration report as JSON artifact."""
        if self.last_calibration is None:
            raise SurrogateNotFittedError("no calibration report — fit first")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.last_calibration.model_dump_json(indent=2), encoding="utf-8")
        return self.last_calibration


def evaluate_surrogate_on_pairs(
    surrogate: DeepEnsembleSurrogate,
    candidates: Candidates,
    experimental: dict[UUID, float],
    *,
    coverage_target: float = 0.90,
    ece_threshold: float = 0.10,
    notes: str | None = None,
) -> CalibrationReport:
    """Score calibrated intervals vs held-out experimental / oracle labels."""
    vectors = surrogate.predict(candidates)
    lowers: list[float] = []
    uppers: list[float] = []
    obs: list[float] = []
    for vec in vectors:
        if vec.candidate_id not in experimental:
            raise KeyError(f"no label for {vec.candidate_id}")
        pred = vec.predictions[0]
        lowers.append(pred.lower)
        uppers.append(pred.upper)
        obs.append(experimental[vec.candidate_id])
    return evaluate_calibration(
        lowers,
        uppers,
        obs,
        coverage_target=coverage_target,
        ece_threshold=ece_threshold,
        notes=notes,
    )


def synthetic_physics_label(sequence: str, *, noise: float = 0.0, seed: int = 0) -> float:
    """Deterministic synthetic oracle for plumbing — MUST be named synthetic_*.

    Linear function of gravy + charge features (not real physics).
    """
    import hashlib

    feats = sequence_features(sequence)
    # FEATURE_NAMES: AA frac…, length_norm, gravy, net_charge
    gravy_v = feats[-2]
    charge_v = feats[-1]
    length_v = feats[-3]
    base = -1.5 * gravy_v - 0.3 * abs(charge_v) + 0.2 * length_v
    if noise == 0.0:
        return base
    digest = hashlib.sha256(f"{seed}:{sequence}".encode()).hexdigest()
    rng = random.Random(int(digest[:16], 16))
    return base + rng.gauss(0.0, noise)

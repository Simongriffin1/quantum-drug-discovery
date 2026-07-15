"""Contract tests: validate accepted shapes and reject invalid payloads.

No business logic — only schema / Protocol surface checks.
Synthetic fixtures use the synthetic_* naming convention.
"""

from __future__ import annotations

from uuid import uuid4

import pytest
from pydantic import ValidationError

from peptideforge.contracts.models import (
    AcquisitionBatch,
    CalibratedPrediction,
    Candidates,
    ComplexStructure,
    DevelopabilityProperty,
    DevelopabilityScores,
    LoopState,
    ObjectiveVector,
    OracleResult,
    OracleTier,
    PeptideCandidate,
    PropertyScore,
    RankedCandidate,
)
from peptideforge.contracts.protocols import (
    AcquisitionFunction,
    Generator,
    Oracle,
    StructurePredictor,
    Surrogate,
)


@pytest.fixture
def synthetic_peptide() -> PeptideCandidate:
    return PeptideCandidate(
        sequence="ACDEFGHIKLM",
        generation_method="synthetic_fixture",
        metadata={"fixture": "synthetic_shape_test"},
    )


@pytest.mark.contract
def test_peptide_candidate_rejects_short_sequence() -> None:
    with pytest.raises(ValidationError):
        PeptideCandidate(sequence="ACDE", generation_method="synthetic_fixture")


@pytest.mark.contract
def test_peptide_candidate_rejects_non_canonical() -> None:
    with pytest.raises(ValidationError):
        PeptideCandidate(sequence="ACDEFGHIKLX", generation_method="synthetic_fixture")


@pytest.mark.contract
def test_candidates_reject_empty(synthetic_peptide: PeptideCandidate) -> None:
    with pytest.raises(ValidationError):
        Candidates(items=())


@pytest.mark.contract
def test_candidates_accept_batch(synthetic_peptide: PeptideCandidate) -> None:
    batch = Candidates(items=(synthetic_peptide,), seed=0)
    assert len(batch.items) == 1
    assert batch.seed == 0


@pytest.mark.contract
def test_complex_requires_structure_payload(synthetic_peptide: PeptideCandidate) -> None:
    with pytest.raises(ValidationError):
        ComplexStructure(
            candidate_id=synthetic_peptide.candidate_id,
            target_id="synthetic_target",
            sequence=synthetic_peptide.sequence,
            confidence=0.9,
            fold_method="synthetic_folder",
        )


@pytest.mark.contract
def test_complex_accepts_pdb_text(synthetic_peptide: PeptideCandidate) -> None:
    complex_ = ComplexStructure(
        candidate_id=synthetic_peptide.candidate_id,
        target_id="synthetic_target",
        sequence=synthetic_peptide.sequence,
        pdb_text="HEADER    synthetic_plumbing\nEND\n",
        confidence=0.85,
        fold_method="synthetic_folder",
    )
    assert complex_.confidence == 0.85


@pytest.mark.contract
def test_oracle_result_qchem_vqe_requires_classical_baseline(
    synthetic_peptide: PeptideCandidate,
) -> None:
    with pytest.raises(ValidationError):
        OracleResult(
            candidate_id=synthetic_peptide.candidate_id,
            value=-1.0,
            uncertainty=0.1,
            cost_estimate=10.0,
            tier=OracleTier.QCHEM_VQE,
        )


@pytest.mark.contract
def test_oracle_result_qchem_vqe_with_baseline(
    synthetic_peptide: PeptideCandidate,
) -> None:
    result = OracleResult(
        candidate_id=synthetic_peptide.candidate_id,
        value=-1.174,
        uncertainty=0.01,
        cost_estimate=5.0,
        tier=OracleTier.QCHEM_VQE,
        classical_baseline=-1.174,
        unit="hartree",
    )
    assert result.classical_baseline == -1.174


@pytest.mark.contract
def test_calibrated_prediction_interval_order(
    synthetic_peptide: PeptideCandidate,
) -> None:
    with pytest.raises(ValidationError):
        CalibratedPrediction(
            candidate_id=synthetic_peptide.candidate_id,
            objective_name="binding",
            mean=0.0,
            lower=1.0,
            upper=-1.0,
            epistemic_std=0.1,
        )


@pytest.mark.contract
def test_developability_scores_vector(synthetic_peptide: PeptideCandidate) -> None:
    scores = DevelopabilityScores(
        candidate_id=synthetic_peptide.candidate_id,
        scores=(
            PropertyScore(
                property_name=DevelopabilityProperty.AGGREGATION,
                value=0.2,
                uncertainty=0.05,
                higher_is_better=False,
                method="synthetic_aggrescan_style",
            ),
            PropertyScore(
                property_name=DevelopabilityProperty.SOLUBILITY,
                value=-0.5,
                higher_is_better=True,
                method="synthetic_gravy",
            ),
        ),
    )
    assert len(scores.scores) == 2


@pytest.mark.contract
def test_acquisition_batch_respects_batch_size(synthetic_peptide: PeptideCandidate) -> None:
    ranked = (
        RankedCandidate(candidate_id=synthetic_peptide.candidate_id, acquisition_score=1.0, rank=0),
        RankedCandidate(candidate_id=uuid4(), acquisition_score=0.5, rank=1),
    )
    with pytest.raises(ValidationError):
        AcquisitionBatch(
            ranked=ranked,
            batch_size=1,
            budget_remaining=10.0,
            method="synthetic_random",
        )


@pytest.mark.contract
def test_loop_state_shape() -> None:
    state = LoopState(
        campaign_id=uuid4(),
        iteration=0,
        seed=42,
        config={"mode": "synthetic"},
        status="running",
    )
    assert state.oracle_calls == 0
    assert state.iteration == 0


@pytest.mark.contract
def test_generator_protocol_runtime_checkable() -> None:
    class _ToyGen:
        def propose(
            self,
            *,
            n: int,
            seed_sequences: tuple[str, ...] | None = None,
            seed: int | None = None,
            constraints: dict[str, object] | None = None,
        ) -> Candidates:
            item = PeptideCandidate(
                sequence="ACDEFGHIKLM",
                generation_method="synthetic_toy",
            )
            return Candidates(items=(item,) * n, seed=seed)

    assert isinstance(_ToyGen(), Generator)


@pytest.mark.contract
def test_oracle_protocol_runtime_checkable(synthetic_peptide: PeptideCandidate) -> None:
    class _ToyOracle:
        def evaluate(
            self,
            complex_structure: ComplexStructure,
            *,
            tier: OracleTier | None = None,
            cost_cap: float | None = None,
        ) -> OracleResult:
            return OracleResult(
                candidate_id=complex_structure.candidate_id,
                complex_id=complex_structure.complex_id,
                value=-5.0,
                uncertainty=1.0,
                cost_estimate=0.1,
                tier=tier or OracleTier.SYNTHETIC,
            )

    assert isinstance(_ToyOracle(), Oracle)


@pytest.mark.contract
def test_structure_predictor_and_surrogate_and_acquisition_protocols() -> None:
    assert hasattr(StructurePredictor, "fold")
    assert hasattr(Surrogate, "predict")
    assert hasattr(AcquisitionFunction, "rank")


@pytest.mark.contract
def test_objective_vector_roundtrip(synthetic_peptide: PeptideCandidate) -> None:
    pred = CalibratedPrediction(
        candidate_id=synthetic_peptide.candidate_id,
        objective_name="binding",
        mean=-6.0,
        lower=-8.0,
        upper=-4.0,
        epistemic_std=0.5,
    )
    vec = ObjectiveVector(candidate_id=synthetic_peptide.candidate_id, predictions=(pred,))
    assert vec.model_dump()["predictions"][0]["mean"] == -6.0

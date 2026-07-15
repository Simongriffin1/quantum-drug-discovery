"""Tests for P5 developability predictors and multi-objective evaluator."""

from __future__ import annotations

import pytest

from peptideforge.contracts.models import DevelopabilityProperty, PeptideCandidate
from peptideforge.contracts.protocols import DevelopabilityPredictor, MultiObjectiveEvaluator
from peptideforge.developability import (
    AggregationPredictor,
    HalfLifePredictor,
    ImmunogenicityPredictor,
    PeptideDevelopabilityEvaluator,
    SolubilityPredictor,
    SynthesizabilityPredictor,
)
from peptideforge.developability.amino_acids import gravy, net_charge_ph74


@pytest.fixture
def hydrophobic_peptide() -> PeptideCandidate:
    return PeptideCandidate(
        sequence="FLIVVFLIVVFLIVV",
        generation_method="synthetic_hydrophobic",
    )


@pytest.fixture
def charged_peptide() -> PeptideCandidate:
    return PeptideCandidate(
        sequence="KKKKDDDDEEEE",
        generation_method="synthetic_charged",
    )


@pytest.fixture
def destabilizing_n_term() -> PeptideCandidate:
    return PeptideCandidate(
        sequence="RAAAAAA",
        generation_method="synthetic_n_end",
    )


@pytest.fixture
def stabilizing_n_term() -> PeptideCandidate:
    return PeptideCandidate(
        sequence="GAAAAAA",
        generation_method="synthetic_n_end",
    )


@pytest.mark.parametrize(
    "predictor_cls",
    [
        AggregationPredictor,
        SolubilityPredictor,
        ImmunogenicityPredictor,
        SynthesizabilityPredictor,
        HalfLifePredictor,
    ],
)
def test_predictor_satisfies_protocol(predictor_cls: type) -> None:
    predictor = predictor_cls()
    assert isinstance(predictor, DevelopabilityPredictor)
    assert predictor.property_name in {p.value for p in DevelopabilityProperty}


def test_aggregation_higher_for_hydrophobic_patch(hydrophobic_peptide: PeptideCandidate) -> None:
    agg = AggregationPredictor()
    hydro = agg.predict(hydrophobic_peptide).scores[0].value
    polar = agg.predict(
        PeptideCandidate(sequence="DDEEGGSS", generation_method="synthetic_polar")
    ).scores[0].value
    assert hydro > polar
    assert agg.predict(hydrophobic_peptide).scores[0].higher_is_better is False


def test_solubility_charged_beats_hydrophobic(
    hydrophobic_peptide: PeptideCandidate,
    charged_peptide: PeptideCandidate,
) -> None:
    sol = SolubilityPredictor()
    assert sol.predict(charged_peptide).scores[0].value > sol.predict(hydrophobic_peptide).scores[0].value
    assert sol.predict(charged_peptide).scores[0].higher_is_better is True


def test_gravy_and_charge_helpers() -> None:
    assert gravy("AAAA") == pytest.approx(1.8)
    # Poly-E should be negative charge at pH 7.4
    assert net_charge_ph74("EEEE") < 0


def test_synthesizability_penalizes_long_hydrophobic() -> None:
    syn = SynthesizabilityPredictor()
    short = syn.predict(
        PeptideCandidate(sequence="AAAAAA", generation_method="synthetic_short")
    ).scores[0].value
    long_hydro = syn.predict(
        PeptideCandidate(
            sequence="FLIVVFLIVVFLIVVFLIVVFLIVVFLIVVFLIVV",
            generation_method="synthetic_long",
        )
    ).scores[0].value
    assert short > long_hydro


def test_half_life_n_end_rule(
    destabilizing_n_term: PeptideCandidate,
    stabilizing_n_term: PeptideCandidate,
) -> None:
    hl = HalfLifePredictor()
    assert (
        hl.predict(stabilizing_n_term).scores[0].value
        > hl.predict(destabilizing_n_term).scores[0].value
    )


def test_immunogenicity_returns_bounded_score() -> None:
    imm = ImmunogenicityPredictor()
    result = imm.predict(
        PeptideCandidate(sequence="FLIVVFLIVV", generation_method="synthetic_immuno")
    )
    score = result.scores[0]
    assert 0.0 <= score.value <= 1.0
    assert score.higher_is_better is False
    assert "not_netmhc" in score.method


def test_multi_evaluator_returns_all_axes() -> None:
    evaluator = PeptideDevelopabilityEvaluator()
    assert isinstance(evaluator, MultiObjectiveEvaluator)
    candidate = PeptideCandidate(sequence="ACDEFGHIKLM", generation_method="synthetic_multi")
    result = evaluator.evaluate(candidate)
    assert result.candidate_id == candidate.candidate_id
    assert len(result.scores) == 5
    names = {s.property_name for s in result.scores}
    assert names == set(DevelopabilityProperty)
    for score in result.scores:
        assert score.uncertainty >= 0.0
        assert score.method


def test_multi_evaluator_no_scalar_collapse() -> None:
    """Vector must retain distinct per-property values — no merged objective."""
    evaluator = PeptideDevelopabilityEvaluator()
    candidate = PeptideCandidate(sequence="FLIVVFLIVVFLIVV", generation_method="synthetic_multi")
    result = evaluator.evaluate(candidate)
    values = [s.value for s in result.scores]
    assert len(set(values)) > 1  # not all identical

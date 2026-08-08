from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import evaluate_candidate_evidence
import summarize_multimodal_evidence
from contracts import ContractError


def test_signed_protein_uses_nonnegative_detection_source() -> None:
    corrected = np.array([[-2.0], [2.0], [-1.0], [3.0]])
    logcounts = np.array([[1.0], [0.0], [0.0], [2.0]])
    result = summarize_multimodal_evidence.summarize_binary_features(
        expression=corrected,
        detection=logcounts,
        groups=pd.Series(["0", "0", "1", "1"]),
        features=pd.Index(["P1"]),
        assay="protein",
        assay_kind="protein",
        feature_space="antibody",
        representation="corrected",
        detection_semantics="protein logcounts > 0",
    ).iloc[0]
    assert result["target_mean"] == 0.0
    assert result["target_fraction_detected"] == 0.5
    assert result["detection_semantics"] == "protein logcounts > 0"


def test_detection_can_be_unavailable_but_signed_detection_is_rejected() -> None:
    expression = np.array([[-1.0], [1.0]])
    groups = pd.Series(["0", "1"])
    result = summarize_multimodal_evidence.summarize_binary_features(
        expression=expression,
        detection=None,
        groups=groups,
        features=pd.Index(["peak"]),
        assay="chromatin",
        assay_kind="atac",
        feature_space="peak",
        representation="tfidf",
    ).iloc[0]
    assert not bool(result["detection_available"])
    assert np.isnan(result["target_fraction_detected"])
    with pytest.raises(ContractError, match="nonnegative"):
        summarize_multimodal_evidence.summarize_binary_features(
            expression=expression,
            detection=expression,
            groups=groups,
            features=pd.Index(["peak"]),
            assay="chromatin",
            assay_kind="atac",
            feature_space="peak",
            representation="tfidf",
        )


def _candidate(**updates) -> dict:
    row = {
        "scope": "s",
        "candidate_key": "c1",
        "parent_id": "0",
        "parent_uid": "s::0",
        "candidate_child_uid": "s::0::candidate_1",
        "n_cells": 2,
        "n_donors": 2,
        "n_captures": 2,
        "max_donor_fraction": 0.25,
        "max_capture_fraction": 0.25,
        "technical_effect_size": 0.01,
        "doublet_fraction": 0.01,
        "incompatible_coexpression_fraction": 0.01,
        "stability_pass": True,
        "candidate_vs_parent_remainder_pass": True,
        "candidate_vs_closest_sibling_pass": True,
        "positive_program_pass": True,
        "expected_negative_program_pass": True,
        "cell_level_separation_pass": True,
        "source_corroboration_pass": True,
        "lane__transcriptome__logcounts__evidence_pass": True,
        "biological_label_L1": "",
        "cl_id": "",
    }
    row.update(updates)
    return row


def _audit() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "scope": ["s"] * 4,
            "parent_id": ["0", "1", "2", "3"],
            "parent_uid": ["s::0", "s::1", "s::2", "s::3"],
            "audit_flagged": [True] * 4,
            "doublet_fraction": [0.01, 0.01, 0.40, 0.01],
            "incompatible_program_fraction": [0.01, 0.01, 0.50, 0.01],
            "donor_dominated": [False, True, False, False],
            "capture_dominated": [False] * 4,
            "qc_outlier": [False] * 4,
        }
    )


def test_configured_assay_gates_and_parent_outcomes(config: dict) -> None:
    evidence = pd.DataFrame(
        [
            _candidate(),
            _candidate(candidate_key="c2", parent_id="1", parent_uid="s::1", candidate_child_uid="s::1::candidate_1", max_donor_fraction=0.8),
            _candidate(candidate_key="c3", parent_id="2", parent_uid="s::2", candidate_child_uid="s::2::candidate_1", incompatible_coexpression_fraction=0.6),
            _candidate(candidate_key="c4", parent_id="3", parent_uid="s::3", candidate_child_uid="s::3::candidate_1", expected_negative_program_pass=False),
        ]
    )
    children = evaluate_candidate_evidence.evaluate_candidates(config, evidence)
    indexed = children.set_index("candidate_key")
    decisions = evaluate_candidate_evidence.decide_parents(
        _audit(), children, config
    ).set_index("parent_id")
    assert bool(indexed.loc["c1", "accepted_child_evidence_gate"])
    assert "donor_nondominance" in indexed.loc["c2", "failed_gates"]
    assert "incompatible_coexpression" in indexed.loc["c3", "failed_gates"]
    assert "expected_negative_program_pass" in indexed.loc["c4", "failed_gates"]
    assert decisions.loc["0", "recommended_outcome"] == "ACCEPT SPLIT"
    assert decisions.loc["1", "recommended_outcome"] == "TECHNICAL/UNRESOLVED"
    assert decisions.loc["2", "recommended_outcome"] == "DOUBLET/MIXED"
    assert decisions.loc["3", "recommended_outcome"] == "RETAIN PARENT"


def test_unconfigured_assays_do_not_create_gates(config: dict) -> None:
    config["adaptive"]["eligibility"].update({"min_donors": None, "min_captures": None})
    config["adaptive"]["audit"].update(
        {
            "donor_max_fraction": None,
            "capture_max_fraction": None,
            "doublet_fraction_max": None,
            "technical_effect_size_max": None,
            "incompatible_program_fraction_min": None,
        }
    )
    evidence = pd.DataFrame([_candidate()]).drop(
        columns=[
            "n_donors",
            "n_captures",
            "max_donor_fraction",
            "max_capture_fraction",
            "technical_effect_size",
            "doublet_fraction",
            "incompatible_coexpression_fraction",
        ]
    )
    assert bool(
        evaluate_candidate_evidence.evaluate_candidates(config, evidence).iloc[0][
            "accepted_child_evidence_gate"
        ]
    )

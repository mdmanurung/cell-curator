from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import cell_curator
from cell_curator.benchmark import evaluate_predictions
from cell_curator.capabilities import build_capability_registry
from cell_curator.config import CellCuratorConfig
from cell_curator.contracts import ContractError, sha256
from cell_curator.evidence import score_marker_programs
from cell_curator.knowledge import LocalMarkerProvider, collect_programs
from cell_curator.models import AssayKind, Capability, CapabilityLevel, MarkerProgram
from cell_curator.review import (
    import_critic_findings,
    import_reconciliation,
    validate_reconciliation,
)


def test_strict_config_and_capability_levels(config: dict) -> None:
    with pytest.raises(ValueError, match="3"):
        MarkerProgram(schema_version=2, label="invalid")
    config["unexpected"] = True
    with pytest.raises(ValueError, match="extra"):
        CellCuratorConfig.model_validate(config)
    config.pop("unexpected")
    config["review"]["require_reconciliation"] = False
    with pytest.raises(ValueError, match="require_reconciliation"):
        CellCuratorConfig.model_validate(config)
    config["review"]["require_reconciliation"] = True
    config["review"]["require_explicit_human_approval"] = False
    with pytest.raises(ValueError, match="explicit_human_approval"):
        CellCuratorConfig.model_validate(config)
    config["review"]["require_explicit_human_approval"] = True
    config["review"]["type_and_state_separate"] = False
    with pytest.raises(ValueError, match="type_and_state_separate"):
        CellCuratorConfig.model_validate(config)
    config["review"]["type_and_state_separate"] = True
    config["knowledge"]["require_local_grounding_for_labels"] = False
    with pytest.raises(ValueError, match="require_local_grounding_for_labels"):
        CellCuratorConfig.model_validate(config)
    config["knowledge"]["require_local_grounding_for_labels"] = True
    config["adaptive"]["eligibility"]["rare_exceptions_require_approval"] = False
    with pytest.raises(ValueError, match="rare_exceptions_require_approval"):
        CellCuratorConfig.model_validate(config)
    config["adaptive"]["eligibility"]["rare_exceptions_require_approval"] = True
    config["adaptive"]["parent_local"]["reject_one_cluster"] = False
    with pytest.raises(ValueError, match="reject_one_cluster"):
        CellCuratorConfig.model_validate(config)
    config["adaptive"]["parent_local"]["reject_one_cluster"] = True
    config["comparators"]["closest_sibling_strategy"] = "arbitrary"
    with pytest.raises(ValueError, match="embedding_centroid"):
        CellCuratorConfig.model_validate(config)
    config["comparators"]["closest_sibling_strategy"] = "embedding_centroid"
    config["input"]["assays"]["transcriptome"]["representations"]["logcounts"][
        "ranking_backend"
    ] = "scipy_mannwhitney"
    model = CellCuratorConfig.model_validate(config)
    capability = build_capability_registry(model)[0]
    assert capability.level is CapabilityLevel.ANNOTATION
    assert not capability.requires_gpu


def test_public_registration_api_returns_strict_config_copies(config: dict) -> None:
    model = CellCuratorConfig.model_validate(config)
    extended = cell_curator.register_assay_representation(
        model,
        assay="transcriptome",
        name="secondary",
        representation={
            "source": "layer",
            "key": "logcounts",
            "ranking_backend": "scipy_mannwhitney",
            "ontology_compatible": True,
        },
    )
    assert "secondary" in extended.input.assays["transcriptome"].representations
    assert "secondary" not in model.input.assays["transcriptome"].representations
    grounded = cell_curator.register_knowledge_provider(
        extended,
        {"kind": "local_markers", "path": "markers.json", "version": "fixture-v1"},
    )
    assert len(grounded.knowledge.providers) == len(model.knowledge.providers) + 1


def test_annotation_mode_requires_local_grounding(config: dict) -> None:
    config["run"]["mode"] = "annotation"
    config["knowledge"]["providers"] = []
    with pytest.raises(ValueError, match="local marker, labeled reference, or cached-table"):
        CellCuratorConfig.model_validate(config)


def test_local_marker_json_supports_signed_programs(tmp_path: Path) -> None:
    path = tmp_path / "markers.json"
    path.write_text(
        json.dumps(
            {
                "Type A": {
                    "pos": ["A1", "A2", "A1"],
                    "neg": ["B1"],
                    "cl_id": "CL:0000001",
                }
            }
        )
    )
    programs = collect_programs([LocalMarkerProvider(path, version="fixture-v1")])
    assert programs == [
        MarkerProgram(
            label="Type A",
            positive=("A1", "A2"),
            negative=("B1",),
            cl_id="CL:0000001",
            source=f"local-markers:{path.resolve()}:fixture-v1",
            source_version="fixture-v1",
        )
    ]


def test_signed_program_scoring_penalizes_expected_negative_markers() -> None:
    matrix = np.array(
        [
            [6.0, 5.0, 0.0],
            [5.0, 6.0, 0.0],
            [0.0, 0.0, 7.0],
            [0.0, 0.0, 6.0],
        ]
    )
    capability = Capability(
        assay="rna",
        representation="logcounts",
        assay_kind=AssayKind.RNA,
        feature_space="gene",
        level=CapabilityLevel.ANNOTATION,
        claims=("cell_identity",),
        backend="scipy_mannwhitney",
    )
    programs = [
        MarkerProgram(label="Type A", positive=("A1", "A2"), negative=("B1",)),
        MarkerProgram(label="Type B", positive=("B1",), negative=("A1", "A2")),
    ]
    result = score_marker_programs(
        matrix,
        pd.Series(["0", "0", "1", "1"]),
        pd.Index(["A1", "A2", "B1"]),
        programs,
        detection=matrix,
        run_id="fixture",
        assay="rna",
        representation="logcounts",
        capability=capability,
    )
    winners = (
        result.sort_values(["cluster_id", "score"], ascending=[True, False])
        .groupby("cluster_id", observed=True)
        .first()["label"]
        .to_dict()
    )
    assert winners == {"0": "Type A", "1": "Type B"}
    wrong = result[(result["cluster_id"] == "0") & (result["label"] == "Type B")].iloc[0]
    assert wrong["negative_violation"] == 1.0


def test_neutral_benchmark_measures_abstention_and_calibration() -> None:
    frame = pd.DataFrame(
        {
            "barcode": ["a", "b", "c", "d"],
            "true_label": ["A", "B", "Unknown", "A"],
            "predicted_label": ["A", "Unknown", "Unknown", "B"],
            "confidence": [0.9, 0.4, 0.8, 0.7],
            "true_path": ["root/A", "root/B", "root", "root/A"],
            "predicted_path": ["root/A", "root", "root", "root/B"],
            "true_mixed": [False, False, True, False],
            "predicted_mixed": [False, False, True, False],
            "true_split_decision": ["retain", "split", "retain", "split"],
            "predicted_split_decision": ["retain", "split", "retain", "retain"],
            "repeat_predicted_label": ["A", "Unknown", "Unknown", "B"],
            "input_invariant": [True] * 4,
            "unsplit_parent_invariant": [True] * 4,
            "runtime_seconds": [1.5] * 4,
            "peak_memory_mb": [64.0] * 4,
            "artifact_complete": [True] * 4,
        }
    )
    result = evaluate_predictions(frame, system="fixture", dataset="synthetic")
    assert result.n == 4
    assert result.coverage == 0.5
    assert result.unknown_precision == 0.5
    assert result.brier_score is not None
    assert result.hierarchy_error is not None
    assert result.mixed_detection_f1 == 1.0
    assert result.split_decision_accuracy == 0.75
    assert result.repeat_reproducibility == 1.0
    assert result.artifact_completeness == 1.0


def test_critic_findings_and_reconciliation_import_atomically(tmp_path: Path) -> None:
    run_root = tmp_path / "run"
    evidence_path = run_root / "evidence/summary/multimodal_label_scores.tsv"
    evidence_path.parent.mkdir(parents=True)
    evidence_path.write_text("cluster_id\tlabel\tscore\n0\tA\t1.0\n")
    findings_input = tmp_path / "findings.tsv"
    reconciliation_input = tmp_path / "reconciliation.tsv"
    pd.DataFrame(
        [
            {
                "finding_id": "f-1",
                "target_kind": "cluster",
                "target_id": "0",
                "critic_pass": "biology",
                "severity": "warning",
                "claim": "Alternative identity remains plausible.",
                "evidence_path": "evidence/summary/multimodal_label_scores.tsv",
                "proposed_resolution": "Inspect the expected-negative program.",
            }
        ]
    ).to_csv(findings_input, sep="\t", index=False)
    finding_receipt = import_critic_findings(run_root, findings_input)
    assert finding_receipt["n_findings"] == 1
    assert (
        pd.read_csv(run_root / "review/critic_findings.tsv", sep="\t").loc[0, "severity"]
        == "WARNING"
    )
    imported_findings = run_root / "review/critic_findings.tsv"
    imported_sha256 = sha256(imported_findings)
    outside = tmp_path / "outside.tsv"
    outside.write_text("not run evidence\n")
    escaped = pd.read_csv(findings_input, sep="\t", keep_default_na=False)
    escaped.loc[0, "evidence_path"] = "../outside.tsv"
    escaped.to_csv(findings_input, sep="\t", index=False)
    with pytest.raises(ContractError, match="escapes"):
        import_critic_findings(run_root, findings_input)
    assert sha256(imported_findings) == imported_sha256

    pd.DataFrame(
        [
            {
                "finding_id": "f-1",
                "status": "RESOLVED",
                "resolution": "Expected-negative evidence excludes the alternative.",
                "reviewer": "fixture-reviewer",
            }
        ]
    ).to_csv(reconciliation_input, sep="\t", index=False)
    receipt = import_reconciliation(run_root, reconciliation_input)
    assert receipt["all_reconciled"] is True

    tampered = pd.read_csv(
        run_root / "review/critic_reconciliation.tsv", sep="\t", keep_default_na=False
    )
    tampered.loc[0, "reviewer"] = ""
    tampered.to_csv(run_root / "review/critic_reconciliation.tsv", sep="\t", index=False)
    with pytest.raises(ContractError, match="non-empty reviewer"):
        validate_reconciliation(run_root)

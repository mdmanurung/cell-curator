from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

from cell_curator.config import CellCuratorConfig, load_config
from cell_curator.contracts import ContractError, sha256
from cell_curator.guided import approve_assumptions, construct_plan
from cell_curator.hierarchy import profile_refined_clusters
from cell_curator.ontology import CellOntology
from cell_curator.pipeline import run_annotation, run_root
from cell_curator.review import (
    build_critic_packet,
    build_evidence_cards,
    build_html_report,
    check_report_completeness,
    resolve_durable_mapping_values,
)
from cell_curator.writeback import validate_writeback, write_back


def _synthetic_input(path: Path) -> None:
    rng = np.random.default_rng(7)
    counts = np.zeros((60, 6), dtype=np.float32)
    counts[:30, :2] = rng.poisson(8, size=(30, 2)) + 1
    counts[30:, 2:4] = rng.poisson(8, size=(30, 2)) + 1
    obs = pd.DataFrame(
        {
            "canonical_cluster": ["0"] * 30 + ["1"] * 30,
            "lineage": ["immune"] * 60,
        },
        index=[f"cell-{index:03d}" for index in range(60)],
    )
    var = pd.DataFrame(
        {"gene_symbol": ["A1", "A2", "B1", "B2", "N1", "N2"]},
        index=["A1", "A2", "B1", "B2", "N1", "N2"],
    )
    adata = ad.AnnData(X=np.log1p(counts), obs=obs, var=var)
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = np.vstack(
        [rng.normal(-3, 0.2, size=(30, 3)), rng.normal(3, 0.2, size=(30, 3))]
    ).astype(np.float32)
    adata.write_h5ad(path)


def _config(tmp_path: Path, source: Path, markers: Path) -> Path:
    root = Path(__file__).resolve().parents[1]
    value = yaml.safe_load((root / "assets" / "config.template.yaml").read_text())
    value["run"].update(
        {
            "run_id": "synthetic-v1",
            "output_root": str(tmp_path / "results"),
            "mode": "annotation",
        }
    )
    value["input"]["path"] = str(source)
    value["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    value["adaptive"]["eligibility"].update(
        {"min_cells": 5, "min_donors": None, "min_captures": None}
    )
    value["adaptive"]["parent_local"]["enabled"] = False
    value["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": str(markers), "version": "synthetic-v1"}
    ]
    value["guidance"] = {
        "interactive": False,
        "preauthorized": True,
        "reviewer": "synthetic-fixture-reviewer",
        "pause_after_each_level": False,
        "assumptions": ["Synthetic marker programs and context are pre-authorized."],
        "hierarchy_levels": ["L1"],
        "context": {
            "organism": "human",
            "tissue": "synthetic tissue",
            "developmental_stage": "adult",
            "condition": "healthy",
            "technology": "dissociated",
            "experimental_context": "synthetic executable integration fixture",
            "composition_covariates": ["lineage"],
            "vocabulary_contract": "local synthetic marker programs",
            "reference_contract": "local_only",
            "visualization_key": "X_pca_harmony",
            "writeback_prefix": "cell_curator",
        },
    }
    value["evidence"].update(
        {
            "min_positive_coverage": 1.0,
            "max_negative_violation": 0.0,
            "min_confidence": 0.6,
            "min_margin": 0.2,
        }
    )
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(value, sort_keys=False))
    return path


def _completed_synthetic_run(
    tmp_path: Path,
    *,
    existing_writeback_contract: bool = False,
) -> tuple[Path, CellCuratorConfig, Path]:
    source = tmp_path / "source.h5ad"
    _synthetic_input(source)
    if existing_writeback_contract:
        adata = ad.read_h5ad(source)
        adata.obs["cell_curator_type"] = pd.Categorical(
            ["Type A"] * 30 + ["Type B"] * 30,
            categories=["Type B", "Type A", "Legacy"],
            ordered=True,
        )
        adata.uns["cell_curator_type_colors"] = np.asarray(["#222222", "#AAAAAA", "#777777"])
        adata.obs["cell_curator_state"] = pd.Categorical(
            ["legacy"] * adata.n_obs,
            categories=["not_assessed", "legacy"],
            ordered=True,
        )
        adata.uns["cell_curator_state_colors"] = np.asarray(["#008800", "#880000"])
        adata.write_h5ad(source)
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            {
                "Type A": {"pos": ["A1", "A2"], "neg": ["B1", "B2"]},
                "Type B": {"pos": ["B1", "B2"], "neg": ["A1", "A2"]},
            }
        )
    )
    config_path = _config(tmp_path, source, markers)
    run_annotation(config_path)
    config = load_config(config_path)
    return source, config, run_root(config)


def _rebuild_strict_review(config: CellCuratorConfig, root: Path) -> None:
    build_evidence_cards(config, root)
    build_critic_packet(config, root)
    build_html_report(config, root, strict=True, level="L1")
    build_html_report(config, root, strict=True)


def _writeback_approval(preview: dict[str, object], *, reviewer: str) -> dict[str, object]:
    return {
        "approved": True,
        "assumptions_approved": True,
        "critics_reconciled": True,
        "reviewer": reviewer,
        "run_id": preview["run_id"],
        "input_sha256": preview["input_sha256"],
        "mapping_sha256": preview["mapping_sha256"],
        "labels_sha256": preview["labels_sha256"],
        "diff_sha256": preview["diff_sha256"],
    }


def test_durable_mapping_resolution_uses_deepest_parent_scoped_review() -> None:
    mapping = pd.DataFrame(
        {
            "barcode": ["a", "b"],
            "canonical_parent": ["0", "1"],
            "leaf_id": ["global::0::candidate::child-a", "global::1"],
        }
    )
    labels = pd.DataFrame(
        [
            {
                "level": "L1",
                "parent_scope": "ROOT",
                "cluster_id": "0",
                "label": "Parent A",
                "cell_state": "normal",
            },
            {
                "level": "L1",
                "parent_scope": "ROOT",
                "cluster_id": "1",
                "label": "Parent B",
                "cell_state": "normal",
            },
            {
                "level": "L2",
                "parent_scope": "Parent A",
                "cluster_id": "child-a",
                "label": "Reviewed child",
                "cell_state": "cycling",
            },
        ]
    )
    l1 = resolve_durable_mapping_values(mapping, labels, max_level="L1")
    full = resolve_durable_mapping_values(mapping, labels)
    assert list(l1["cell_type"]) == ["Parent A", "Parent B"]
    assert list(full["cell_type"]) == ["Reviewed child", "Parent B"]
    assert list(full["cell_state"]) == ["cycling", "normal"]


def test_independent_end_to_end_and_guarded_writeback(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    _synthetic_input(source)
    source_before = sha256(source)
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            {
                "Type A": {"pos": ["A1", "A2"], "neg": ["B1", "B2"]},
                "Type B": {"pos": ["B1", "B2"], "neg": ["A1", "A2"]},
            }
        )
    )
    config_path = _config(tmp_path, source, markers)
    packet = run_annotation(config_path)
    config = load_config(config_path)
    root = run_root(config)
    assert packet.is_dir()
    assert (packet / "review_packet.manifest.json").is_file()
    assert sha256(source) == source_before
    assert (root / "adaptive/technical_summary.tsv").is_file()
    assert (root / "adaptive/doublet_summary.tsv").is_file()
    assert (root / "evidence/summary/program_cooccurrence.tsv").is_file()
    assert (root / "evidence/summary/source_corroboration.tsv").is_file()
    report = (root / "review/report.html").read_text()
    assert "Program co-occurrence" in report
    assert "Independent-source corroboration" in report

    proposals = pd.read_csv(
        root / "evidence" / "summary" / "annotation_proposals.tsv",
        sep="\t",
        keep_default_na=False,
    ).set_index("cluster_id")
    assert proposals.loc[0, "cell_type"] == "Type A"
    assert proposals.loc[1, "cell_type"] == "Type B"
    assert not proposals["abstained"].astype(bool).any()

    preview = validate_writeback(config, root)
    invalid_approval = tmp_path / "invalid-approval.json"
    invalid_approval.write_text(json.dumps({"approved": True, "reviewer": "fixture-reviewer"}))
    with pytest.raises(ContractError, match="assumptions_approved"):
        write_back(
            config,
            root,
            output_path=tmp_path / "must-not-exist.h5ad",
            approval_path=invalid_approval,
        )
    assert not (tmp_path / "must-not-exist.h5ad").exists()
    approval = tmp_path / "approval.json"
    approval.write_text(
        json.dumps(
            {
                "approved": True,
                "assumptions_approved": True,
                "critics_reconciled": True,
                "reviewer": "fixture-reviewer",
                "run_id": config.run.run_id,
                "input_sha256": preview["input_sha256"],
                "mapping_sha256": preview["mapping_sha256"],
                "labels_sha256": preview["labels_sha256"],
                "diff_sha256": preview["diff_sha256"],
            }
        )
    )
    output = tmp_path / "annotated.h5ad"
    write_back(
        config,
        root,
        output_path=output,
        approval_path=approval,
    )
    assert output.is_file()
    assert sha256(source) == source_before
    annotated = ad.read_h5ad(output)
    assert set(annotated.obs["cell_curator_type"].astype(str)) == {"Type A", "Type B"}


def test_review_rebuild_preserves_manual_reasoning_and_detects_tampering(
    tmp_path: Path,
) -> None:
    _, config, root = _completed_synthetic_run(tmp_path)
    labels_path = root / "labels/L1.tsv"
    labels = pd.read_csv(labels_path, sep="\t", keep_default_na=False)
    labels.loc[0, "rationale"] = "Manual reviewer rationale retained verbatim."
    labels.loc[0, "reviewer"] = "reviewer-one"
    labels.to_csv(labels_path, sep="\t", index=False)

    with pytest.raises(ContractError, match="stale"):
        check_report_completeness(config, root, require_report=True)
    _rebuild_strict_review(config, root)
    card_path = root / "review/evidence_cards/L1__ROOT__cluster_0.md"
    assert "Manual reviewer rationale retained verbatim." in card_path.read_text()
    assert (
        "Manual reviewer rationale retained verbatim."
        in (root / "review/report.html").read_text()
    )
    assert check_report_completeness(config, root)["status"] == "complete"

    card_path.write_text(card_path.read_text() + "\nTAMPERED\n")
    with pytest.raises(ContractError, match="stale"):
        check_report_completeness(config, root)

    _rebuild_strict_review(config, root)
    report_path = root / "review/report.html"
    report_path.write_text(
        report_path.read_text() + '<img src="https://example.invalid/pixel.png">'
    )
    manifest_path = root / "review/report.manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["report_sha256"] = sha256(report_path)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ContractError, match="self-contained"):
        check_report_completeness(config, root)


def test_writeback_uses_durable_labels_and_previews_categories_colors_and_cells(
    tmp_path: Path,
) -> None:
    source, config, root = _completed_synthetic_run(tmp_path, existing_writeback_contract=True)
    source_before = sha256(source)
    labels_path = root / "labels/L1.tsv"
    labels = pd.read_csv(labels_path, sep="\t", keep_default_na=False)
    labels.loc[labels["cluster_id"].astype(str).eq("0"), "label"] = "Reviewed A"
    labels.loc[labels["cluster_id"].astype(str).eq("0"), "rationale"] = (
        "Reviewer-selected durable identity after inspecting raw evidence."
    )
    labels.loc[labels["cluster_id"].astype(str).eq("0"), "cl_id"] = "CL:0000001"
    labels.loc[labels["cluster_id"].astype(str).eq("1"), "cl_id"] = "CL:0000002"
    labels.to_csv(labels_path, sep="\t", index=False)
    ontology_path = tmp_path / "writeback-cl.obo"
    ontology_path.write_text(
        """format-version: 1.2

[Term]
id: CL:0000000
name: cell

[Term]
id: CL:0000001
name: Reviewed A
is_a: CL:0000000 ! cell

[Term]
id: CL:0000002
name: Type B
is_a: CL:0000000 ! cell
"""
    )
    ontology = CellOntology.load(ontology_path)
    _rebuild_strict_review(config, root)

    invalid_labels = labels.copy()
    invalid_labels.loc[invalid_labels["cluster_id"].astype(str).eq("0"), "cl_id"] = "CL:0000002"
    invalid_labels.to_csv(labels_path, sep="\t", index=False)
    _rebuild_strict_review(config, root)
    with pytest.raises(ContractError, match="inconsistent"):
        validate_writeback(config, root, ontology=ontology)
    labels.to_csv(labels_path, sep="\t", index=False)
    _rebuild_strict_review(config, root)

    preview = validate_writeback(config, root, ontology=ontology)
    assert sha256(source) == source_before
    assert preview["n_cell_types"] == 2
    assert preview["labels_sha256"] == {"labels/L1.tsv": sha256(labels_path)}
    assert preview["ontology_validation"]["ontology_status"] == "validated"
    assert preview["hierarchy_validation"]["invariants"]["one_cell_one_leaf"] is True
    diff = pd.read_csv(preview["diff_path"], sep="\t", keep_default_na=False)
    metadata = diff[diff["barcode"].eq("<METADATA>")]
    assert len(diff) == 3 * 60 + 9
    assert len(metadata) == 9
    assert {
        "cell_curator_type.__categories__",
        "cell_curator_type.__colors__",
        "cell_curator_type.__ordered__",
        "cell_curator_state.__categories__",
        "cell_curator_state.__colors__",
        "cell_curator_state.__ordered__",
        "cell_curator_leaf.__categories__",
        "cell_curator_leaf.__colors__",
        "cell_curator_leaf.__ordered__",
    } == set(metadata["column"])
    type_contract = preview["category_and_color_contract"]["cell_curator_type"]
    assert type_contract["categories"] == ["Type B", "Reviewed A"]
    assert type_contract["colors"][0] == "#222222"
    state_contract = preview["category_and_color_contract"]["cell_curator_state"]
    assert state_contract == {
        "categories": ["not_assessed", "legacy"],
        "colors": ["#008800", "#880000"],
    }

    stale_approval = tmp_path / "stale-approval.json"
    stale_approval.write_text(json.dumps(_writeback_approval(preview, reviewer="reviewer-one")))
    labels.loc[0, "rationale"] = "Final manually reviewed rationale."
    labels.to_csv(labels_path, sep="\t", index=False)
    _rebuild_strict_review(config, root)
    with pytest.raises(ContractError, match="labels_sha256"):
        write_back(
            config,
            root,
            output_path=tmp_path / "stale-must-not-exist.h5ad",
            approval_path=stale_approval,
            ontology=ontology,
        )
    assert not (tmp_path / "stale-must-not-exist.h5ad").exists()

    current_preview = validate_writeback(config, root, ontology=ontology)
    approval_path = tmp_path / "approval.json"
    approval_path.write_text(
        json.dumps(_writeback_approval(current_preview, reviewer="reviewer-one"))
    )
    output = tmp_path / "reviewed.h5ad"
    write_back(
        config,
        root,
        output_path=output,
        approval_path=approval_path,
        ontology=ontology,
    )
    assert sha256(source) == source_before
    written = ad.read_h5ad(output)
    assert set(written.obs["cell_curator_type"].astype(str)) == {
        "Reviewed A",
        "Type B",
    }
    assert list(written.obs["cell_curator_type"].cat.categories) == type_contract["categories"]
    assert list(map(str, written.uns["cell_curator_type_colors"])) == type_contract["colors"]
    assert (
        list(written.obs["cell_curator_state"].cat.categories) == state_contract["categories"]
    )
    assert list(map(str, written.uns["cell_curator_state_colors"])) == state_contract["colors"]


def test_interactive_equivalent_end_to_end_requires_named_gate_approval(
    tmp_path: Path,
) -> None:
    source = tmp_path / "interactive.h5ad"
    _synthetic_input(source)
    markers = tmp_path / "interactive-markers.json"
    markers.write_text(
        json.dumps(
            {
                "Type A": {"pos": ["A1", "A2"], "neg": ["B1", "B2"]},
                "Type B": {"pos": ["B1", "B2"], "neg": ["A1", "A2"]},
            }
        )
    )
    config_path = _config(tmp_path, source, markers)
    value = yaml.safe_load(config_path.read_text())
    value["run"]["run_id"] = "interactive-v1"
    value["guidance"].update({"interactive": True, "preauthorized": False, "reviewer": ""})
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))

    with pytest.raises(ContractError, match="paused after scene setting"):
        run_annotation(config_path)
    construct_plan(config_path)
    with pytest.raises(ContractError, match="Gate A blocks"):
        run_annotation(config_path)
    approval = approve_assumptions(config_path, reviewer="interactive-reviewer")
    assert approval["approval_mode"] == "interactive"
    packet = run_annotation(config_path)
    root = tmp_path / "results" / "interactive-v1"
    assert packet.is_dir()
    assert (root / "labels/L1.tsv").is_file()
    assert (root / "checkpoints/L1/ROOT.json").is_file()


def test_stored_refinement_promotes_only_evidence_backed_children(
    tmp_path: Path,
) -> None:
    rng = np.random.default_rng(11)
    n = 40
    barcodes = [f"split-{index:02d}" for index in range(n)]
    child = np.array(["a"] * 20 + ["b"] * 20)
    counts = np.ones((n, 3), dtype=np.float32)
    counts[:20, 0] = rng.poisson(8, size=20) + 1
    counts[:20, 1] = 0
    counts[20:, 0] = 0
    counts[20:, 1] = rng.poisson(8, size=20) + 1
    obs = pd.DataFrame(
        {
            "canonical_cluster": "P",
            "lineage": "immune",
            "stored_low": child,
            "stored_high": child,
        },
        index=barcodes,
    )
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=obs,
        var=pd.DataFrame({"gene_symbol": ["A", "B", "HK"]}, index=["A", "B", "HK"]),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca"] = np.column_stack(
        [
            np.r_[rng.normal(-4, 0.1, 20), rng.normal(4, 0.1, 20)],
            rng.normal(size=n),
        ]
    )
    source = tmp_path / "split.h5ad"
    adata.write_h5ad(source)
    markers = tmp_path / "split_markers.json"
    markers.write_text(
        json.dumps(
            {
                "Type A": {"positive": ["A"], "negative": ["B"]},
                "Type B": {"positive": ["B"], "negative": ["A"]},
            }
        )
    )
    config_path = _config(tmp_path, source, markers)
    config = yaml.safe_load(config_path.read_text())
    config["run"]["run_id"] = "accepted-split-v1"
    config["input"]["keys"].update(
        {"embedding": "X_pca", "stored_resolutions": ["stored_low", "stored_high"]}
    )
    config["guidance"]["context"]["visualization_key"] = "X_pca"
    config["input"]["assays"]["transcriptome"]["graph_representations"] = [
        {"source": "obsm", "key": "X_pca"}
    ]
    config["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 1}
    }
    config["adaptive"]["parent_local"]["enabled"] = True
    config["adaptive"]["stored_candidates"].update(
        {"min_parent_fraction": 0.1, "max_parent_fraction": 0.9}
    )
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    packet = run_annotation(config_path)
    root = tmp_path / "results" / "accepted-split-v1"
    assert packet.is_dir()
    profile = profile_refined_clusters(config_path)
    assert {"size", "qc", "composition"}.issubset(profile["panels"])
    assert Path(profile["profile_path"]).is_file()
    manifest = pd.read_csv(root / "comparators" / "comparator_manifest.tsv", sep="\t")
    assert len(manifest) == 4
    assert manifest["lane_id"].nunique() == 4
    assert all(
        (root / "evidence" / "candidate_lanes" / lane_id / "receipt.json").is_file()
        for lane_id in manifest["lane_id"].astype(str)
    )
    children = pd.read_csv(root / "decisions" / "child_evidence_summary.tsv", sep="\t")
    assert children["accepted_child_evidence_gate"].astype(bool).all()
    assert set(children["cell_type"]) == {"Type A", "Type B"}
    parents = pd.read_csv(root / "decisions" / "parent_decision_summary.tsv", sep="\t")
    assert parents.loc[0, "recommended_outcome"] == "ACCEPT SPLIT"
    mapping = pd.read_csv(root / "mapping" / "parent_child_mapping.tsv", sep="\t")
    assert mapping["barcode"].nunique() == n
    assert mapping["leaf_id"].nunique() == 2
    assert set(mapping["membership_role"]) == {"accepted_child"}
    assert run_annotation(config_path) == packet


def test_zero_applicable_lanes_complete_with_unknown_calls(tmp_path: Path) -> None:
    source = tmp_path / "zero-lane.h5ad"
    _synthetic_input(source)
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps({"Type A": {"positive": ["A1"]}, "Type B": {"positive": ["B1"]}})
    )
    config_path = _config(tmp_path, source, markers)
    value = yaml.safe_load(config_path.read_text())
    value["run"]["run_id"] = "zero-lane-v1"
    value["evidence"]["lanes"] = []
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))

    packet = run_annotation(config_path)
    root = tmp_path / "results" / "zero-lane-v1"
    receipt_manifest = json.loads((root / "evidence/receipts.validated.json").read_text())
    proposals = pd.read_csv(root / "evidence/summary/annotation_proposals.tsv", sep="\t")
    assert packet.is_dir()
    assert receipt_manifest["zero_lane_run"] is True
    assert proposals["cell_type"].eq("Unknown").all()
    assert proposals["abstained"].astype(bool).all()

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cell_curator.config import CellCuratorConfig
from cell_curator.contracts import ContractError, sha256
from cell_curator.models import MarkerProgram
from cell_curator.ontology import (
    CellOntology,
    ontology_from_config,
    validate_durable_label_ontology,
    validate_program_ontology,
)
from cell_curator.refinement import (
    adjusted_rand_index,
    evaluate_local_partitions,
    mutual_best_matches,
    partition_jaccard_index,
)


def test_local_obo_resolution_obsolete_and_hierarchy(tmp_path: Path) -> None:
    path = tmp_path / "cl.obo"
    path.write_text(
        """format-version: 1.2

[Term]
id: CL:0000000
name: cell

[Term]
id: CL:0000001
name: lymphocyte
synonym: "lymphoid cell" EXACT []
is_a: CL:0000000 ! cell

[Term]
id: CL:9999999
name: old lymphocyte
is_obsolete: true
replaced_by: CL:0000001
"""
    )
    ontology = CellOntology.load(path)
    assert ontology.resolve("lymphoid cell").accession == "CL:0000001"
    assert ontology.is_descendant("CL:0000001", "CL:0000000")
    assert ontology.validate_term("CL:0000001", expected_name="lymphocyte").name == "lymphocyte"
    validate_program_ontology(
        [
            MarkerProgram(
                label="lymphocyte",
                positive=("PTPRC",),
                cl_id="CL:0000001",
                parent_cl_id="CL:0000000",
            )
        ],
        ontology,
    )
    with pytest.raises(ContractError, match="inconsistent"):
        validate_program_ontology(
            [
                MarkerProgram(
                    label="monocyte",
                    positive=("LYZ",),
                    cl_id="CL:0000001",
                )
            ],
            ontology,
        )


def test_durable_labels_enforce_name_obsolete_and_parent_semantics(
    tmp_path: Path,
) -> None:
    path = tmp_path / "durable-cl.obo"
    path.write_text(
        """format-version: 1.2

[Term]
id: CL:0000000
name: cell

[Term]
id: CL:0000001
name: lymphocyte
is_a: CL:0000000 ! cell

[Term]
id: CL:0000002
name: T cell
synonym: "T lymphocyte" EXACT []
is_a: CL:0000001 ! lymphocyte

[Term]
id: CL:0000003
name: epithelial cell
is_a: CL:0000000 ! cell

[Term]
id: CL:9999999
name: old T cell
is_obsolete: true
replaced_by: CL:0000002
"""
    )
    ontology = CellOntology.load(path)
    labels = pd.DataFrame(
        [
            {
                "level": "L1",
                "parent_scope": "ROOT",
                "cluster_id": "immune",
                "label": "lymphocyte",
                "cl_id": "CL:0000001",
                "parent_cl_id": "",
            },
            {
                "level": "L2",
                "parent_scope": "lymphocyte",
                "cluster_id": "t",
                "label": "T lymphocyte",
                "cl_id": "CL:0000002",
                "parent_cl_id": "CL:0000001",
            },
        ]
    )
    result = validate_durable_label_ontology(labels, ontology)
    assert result["ontology_status"] == "validated"
    assert result["n_supplied_cl_ids"] == 2
    assert result["n_validated_parent_links"] == 1

    cl_free = labels.assign(cl_id="", parent_cl_id="")
    assert validate_durable_label_ontology(cl_free, None)["ontology_status"] == (
        "unavailable-cl-free"
    )

    mismatched = labels.copy()
    mismatched.loc[1, "label"] = "epithelial cell"
    with pytest.raises(ContractError, match="inconsistent"):
        validate_durable_label_ontology(mismatched, ontology)

    incompatible = labels.copy()
    incompatible.loc[1, ["label", "cl_id"]] = ["epithelial cell", "CL:0000003"]
    with pytest.raises(ContractError, match="not compatible"):
        validate_durable_label_ontology(incompatible, ontology)

    obsolete = labels.copy()
    obsolete.loc[1, ["label", "cl_id"]] = ["old T cell", "CL:9999999"]
    with pytest.raises(ContractError, match="obsolete.*replaced by CL:0000002"):
        validate_durable_label_ontology(obsolete, ontology)

    unknown = labels.iloc[[0]].assign(label="Unknown", cl_id="CL:0000001")
    with pytest.raises(ContractError, match="Unknown.*cannot carry"):
        validate_durable_label_ontology(unknown, ontology)


def test_checksum_pinned_cached_remote_ontology_and_unavailable_receipt(
    tmp_path: Path, config: dict
) -> None:
    cache = tmp_path / "cached-cl.obo"
    cache.write_text(
        """format-version: 1.2

[Term]
id: CL:0000000
name: cell
"""
    )
    config["knowledge"]["providers"] = [
        {
            "kind": "remote_ontology",
            "url": "https://example.invalid/cl.obo",
            "cache_path": str(cache),
            "sha256": sha256(cache),
            "version": "fixture-v1",
            "timeout_seconds": 0.1,
        }
    ]
    model = CellCuratorConfig.model_validate(config)
    ontology = ontology_from_config(model)
    assert ontology is not None
    assert ontology.resolve("cell").accession == "CL:0000000"
    status_path = (
        Path(model.run.output_root)
        / model.run.run_id
        / "knowledge/ontology_provider_status.json"
    )
    assert json.loads(status_path.read_text())["status"] == "complete"

    config["run"]["run_id"] = "missing-ontology-v1"
    config["knowledge"]["providers"][0].update(
        {
            "cache_path": str(tmp_path / "missing.obo"),
            "url": "file:///definitely/missing/cl.obo",
        }
    )
    unavailable_model = CellCuratorConfig.model_validate(config)
    assert ontology_from_config(unavailable_model) is None
    unavailable = json.loads(
        (
            Path(unavailable_model.run.output_root)
            / unavailable_model.run.run_id
            / "knowledge/ontology_provider_status.json"
        ).read_text()
    )
    assert unavailable["status"] == "unavailable"
    assert "CL-free" in unavailable["remediation"]


def test_refinement_metrics_and_graph_partition_are_computed(config: dict) -> None:
    assert adjusted_rand_index(np.array([0, 0, 1, 1]), np.array([1, 1, 0, 0])) == 1.0
    assert partition_jaccard_index(
        np.array([0, 0, 1, 1]), np.array([1, 1, 0, 0])
    ) == pytest.approx(1.0)
    matches = mutual_best_matches(
        {"a": {"1", "2"}, "b": {"3", "4"}},
        {"x": {"1", "2"}, "y": {"3", "4"}},
        threshold=0.8,
    )
    assert matches == [("a", "x", 1.0), ("b", "y", 1.0)]

    rng = np.random.default_rng(3)
    embedding = np.vstack(
        [rng.normal(-3, 0.15, size=(20, 3)), rng.normal(3, 0.15, size=(20, 3))]
    )
    config["input"]["keys"].update({"donor": "", "capture": "", "batch": ""})
    config["adaptive"]["eligibility"].update(
        {"min_cells": 5, "min_donors": None, "min_captures": None}
    )
    config["adaptive"]["parent_local"].update(
        {
            "resolutions": [0.5, 1.0],
            "stability_threshold_strictly_greater_than": 0.0,
            "reference_seeds": [0, 1, 2],
            "knn_k": 8,
        }
    )
    config["input"]["assays"]["transcriptome"]["representations"]["logcounts"][
        "ranking_backend"
    ] = "scipy_mannwhitney"
    model = CellCuratorConfig.model_validate(config)
    evaluations = evaluate_local_partitions(
        embedding,
        pd.DataFrame(index=[f"c{i}" for i in range(40)]),
        model,
    )
    assert any(item.n_clusters >= 2 for item in evaluations)
    assert all(0 <= item.stability <= 1 for item in evaluations)
    assert all(item.separation >= 0 for item in evaluations)

    config["adaptive"]["parent_local"]["stability_metric"] = "jaccard"
    jaccard_model = CellCuratorConfig.model_validate(config)
    jaccard_evaluations = evaluate_local_partitions(
        embedding,
        pd.DataFrame(index=[f"c{i}" for i in range(40)]),
        jaccard_model,
    )
    assert jaccard_evaluations
    assert {item.stability_metric for item in jaccard_evaluations} == {"jaccard"}
    assert all(0 <= item.stability <= 1 for item in jaccard_evaluations)

    config["adaptive"]["parent_local"]["max_passes_per_parent"] = 0
    disabled_model = CellCuratorConfig.model_validate(config)
    assert (
        evaluate_local_partitions(
            embedding,
            pd.DataFrame(index=[f"c{i}" for i in range(40)]),
            disabled_model,
        )
        == []
    )

from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from cell_curator.config import CellCuratorConfig
from cell_curator.contracts import ContractError, sha256
from cell_curator.models import MarkerProgram
from cell_curator.ontology import (
    CellOntology,
    OntologyTerm,
    corroborate_terms_live,
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


def test_depth_and_roots_describe_the_supplied_graph_not_the_cell_ontology(
    tmp_path: Path,
) -> None:
    """A subtree export has a cut at the top, and the cut is what a root means.

    Depth is the number of `is_a` edges to a root of the file that was loaded, so
    a term whose parent was not exported counts as a root and a term that hangs
    off nothing measurable resolves to 0. A distance built on these two
    primitives is only meaningful against the snapshot it names, which is why
    they report the graph rather than a hard-coded Cell Ontology root. The
    dangling parent also stays in `ancestors`, so a consumer comparing two
    accessions has to check that a term is present before measuring it.
    """

    path = tmp_path / "subtree-cl.tsv"
    path.write_text(
        """id\tname\tparents
CL:0000000\tcell\t
CL:0000738\tleukocyte\tCL:0000000
CL:0000542\tlymphocyte\tCL:0000738
CL:0000084\tT cell\tCL:0000542
CL:0000786\tplasma cell\tCL:0009000
"""
    )
    ontology = CellOntology.load(path)
    assert ontology.roots() == frozenset({"CL:0000000", "CL:0000786"})
    assert ontology.depth("CL:0000000") == 0
    assert ontology.depth("CL:0000738") == 1
    assert ontology.depth("CL:0000542") == 2
    assert ontology.depth("CL:0000084") == 3
    assert ontology.depth("CL:0000084") == ontology.depth("CL:0000084")
    with pytest.raises(ContractError, match="unknown"):
        ontology.resolve("CL:0009000")
    assert ontology.ancestors("CL:0000786") == {"CL:0009000"}


def test_a_parent_cycle_resolves_to_depth_zero_rather_than_hanging() -> None:
    """A malformed snapshot must not be able to produce a plausible depth.

    A cycle has no top, so there is no distance to measure in it. The walk is
    iterative and memoized precisely so a file with a cycle terminates and
    reports 0 for the terms involved, instead of exhausting the stack or
    inventing a depth a reader would take for a measurement.
    """

    ontology = CellOntology(
        {
            "CL:0000000": OntologyTerm("CL:0000000", "a", parents=("CL:0000001",)),
            "CL:0000001": OntologyTerm("CL:0000001", "b", parents=("CL:0000000",)),
            "CL:0000002": OntologyTerm("CL:0000002", "c", parents=("CL:0000001",)),
        },
        source="cycle",
    )
    assert ontology.roots() == frozenset()
    assert ontology.depth("CL:0000000") == 0
    assert ontology.depth("CL:0000001") == 1
    assert ontology.depth("CL:0000002") == 2


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


def _cl_obo(tmp_path: Path) -> Path:
    path = tmp_path / "local-cl.obo"
    path.write_text(
        """format-version: 1.2

[Term]
id: CL:0000084
name: T cell
synonym: "T-lymphocyte" EXACT []
"""
    )
    return path


def _serve_ols(monkeypatch: pytest.MonkeyPatch, payload: object) -> list[str]:
    """Answer one OLS-shaped query without reaching the network."""

    requested: list[str] = []

    class _Response:
        def __init__(self, body: bytes) -> None:
            self._buffer = io.BytesIO(body)

        def read(self, size: int = -1) -> bytes:
            return self._buffer.read(size)

        def __enter__(self) -> Any:
            return self

        def __exit__(self, *exc: object) -> bool:
            return False

    def fake_urlopen(request: Any, timeout: float = 0.0) -> Any:
        requested.append(request.full_url)
        return _Response(json.dumps(payload).encode())

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    return requested


def _with_live_provider(config: dict, tmp_path: Path, **overrides: Any) -> CellCuratorConfig:
    provider = {
        "kind": "live_ontology",
        "url": "https://service.invalid/api/ontologies/cl/terms",
        "version": "fixture-live-v1",
        "timeout_seconds": 0.5,
    }
    provider.update(overrides)
    config["knowledge"]["providers"] = [
        {"kind": "ontology", "path": str(_cl_obo(tmp_path)), "version": "local-v1"},
        provider,
    ]
    return CellCuratorConfig.model_validate(config)


def test_live_ontology_provider_takes_only_an_https_endpoint(config: dict) -> None:
    """A live provider has no digest to pin, so it must not pretend to carry one."""

    config["knowledge"]["providers"] = [
        {"kind": "live_ontology", "url": "http://service.invalid/terms"}
    ]
    with pytest.raises(ValueError, match="requires an https url"):
        CellCuratorConfig.model_validate(config)

    config["knowledge"]["providers"] = [
        {
            "kind": "live_ontology",
            "url": "https://service.invalid/terms",
            "sha256": "a" * 64,
        }
    ]
    with pytest.raises(ValueError, match="only url, version, and timeout_seconds"):
        CellCuratorConfig.model_validate(config)


def test_live_corroboration_without_a_provider_never_attempts_the_network(
    config: dict,
) -> None:
    model = CellCuratorConfig.model_validate(config)
    status = corroborate_terms_live(model, ["CL:0000084"])
    assert status["status"] == "unavailable"
    assert status["network_attempted"] is False
    assert "no live ontology provider" in status["reason"]


def test_live_corroboration_disabled_by_policy_never_attempts_the_network(
    config: dict, tmp_path: Path
) -> None:
    config["knowledge"]["allow_optional_remote_enhancement"] = False
    model = _with_live_provider(config, tmp_path)
    status = corroborate_terms_live(model, ["CL:0000084"])
    assert status["status"] == "unavailable"
    assert status["network_attempted"] is False
    assert "CL-free annotation" in status["remediation"]


def test_live_corroboration_records_agreement_with_a_hashed_answer(
    monkeypatch: pytest.MonkeyPatch, config: dict, tmp_path: Path
) -> None:
    """Reproducibility comes from the stored answer, not from the live endpoint."""

    model = _with_live_provider(config, tmp_path)
    requested = _serve_ols(
        monkeypatch,
        {
            "_embedded": {
                "terms": [
                    {
                        "obo_id": "CL:0000084",
                        "label": "T cell",
                        "is_obsolete": False,
                        "synonyms": ["T-lymphocyte"],
                    }
                ]
            }
        },
    )

    status = corroborate_terms_live(
        model,
        ["CL:0000084"],
        expected_names=["T cell"],
        ontology=ontology_from_config(model),
    )

    assert requested == [
        "https://service.invalid/api/ontologies/cl/terms?obo_id=CL%3A0000084"
    ]
    assert status["status"] == "complete"
    assert status["n_agreement"] == 1
    assert status["n_disagreement"] == 0
    row = status["corroborated"][0]
    assert row["status"] == "agreement"
    assert row["live_name"] == "T cell"
    assert row["local_name"] == "T cell"
    stored = Path(row["artifact"])
    assert sha256(stored) == row["sha256"]
    assert json.loads(stored.read_text())["_embedded"]["terms"][0]["label"] == "T cell"


def test_live_corroboration_reports_obsolescence_and_name_disagreement(
    monkeypatch: pytest.MonkeyPatch, config: dict, tmp_path: Path
) -> None:
    model = _with_live_provider(config, tmp_path)
    _serve_ols(
        monkeypatch,
        {
            "_embedded": {
                "terms": [
                    {
                        "obo_id": "CL:0000084",
                        "label": "alpha-beta T cell",
                        "is_obsolete": True,
                        "term_replaced_by": "CL:0000789",
                    }
                ]
            }
        },
    )

    status = corroborate_terms_live(
        model,
        ["CL:0000084"],
        expected_names=["T cell"],
        ontology=ontology_from_config(model),
    )

    row = status["corroborated"][0]
    assert status["n_disagreement"] == 1
    assert row["status"] == "disagreement"
    assert row["live_replaced_by"] == "CL:0000789"
    joined = " | ".join(row["findings"])
    assert "obsolete" in joined
    assert "disagrees with live name" in joined
    assert "not the live name or a live synonym" in joined


def test_an_unreachable_ontology_service_is_a_receipt_not_an_exception(
    monkeypatch: pytest.MonkeyPatch, config: dict, tmp_path: Path
) -> None:
    def refuse(request: Any, timeout: float = 0.0) -> Any:
        raise OSError("network is unreachable")

    monkeypatch.setattr("urllib.request.urlopen", refuse)
    model = _with_live_provider(config, tmp_path)

    status = corroborate_terms_live(model, ["CL:0000084"])

    assert status["status"] == "complete"
    assert status["n_unavailable"] == 1
    assert "network is unreachable" in " ".join(status["corroborated"][0]["findings"])


def test_a_live_provider_never_becomes_required_for_a_label(
    monkeypatch: pytest.MonkeyPatch, config: dict, tmp_path: Path
) -> None:
    """The whole point of the optional tier: offline still reaches validated labels.

    `require_local_grounding_for_labels` is pinned true, so a configured live
    service may corroborate and nothing more. With the network refused, the local
    ontology must still resolve terms and still validate durable labels.
    """

    def refuse(request: Any, timeout: float = 0.0) -> Any:
        raise OSError("network is unreachable")

    monkeypatch.setattr("urllib.request.urlopen", refuse)
    model = _with_live_provider(config, tmp_path)

    ontology = ontology_from_config(model)
    assert ontology is not None
    assert ontology.resolve("CL:0000084").name == "T cell"
    assert ontology.validate_term("CL:0000084", expected_name="T cell").name == "T cell"

    labels = pd.DataFrame(
        {
            "level": ["L1"],
            "parent_scope": ["ROOT"],
            "cluster_id": ["0"],
            "label": ["T cell"],
            "cl_id": ["CL:0000084"],
        }
    )
    assert validate_durable_label_ontology(labels, ontology)["status"] == "complete"

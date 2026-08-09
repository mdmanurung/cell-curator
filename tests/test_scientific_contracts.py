from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

from cell_curator.capabilities import build_capability_registry
from cell_curator.config import CellCuratorConfig
from cell_curator.contracts import (
    ContractError,
    require_safe_path_component,
    safe_path_component,
)
from cell_curator.data import feature_ids_for
from cell_curator.evidence import (
    aggregate_multimodal_evidence,
    compute_bottom_up_markers,
    compute_gpu_bottom_up_markers,
    score_marker_programs,
)
from cell_curator.models import MarkerProgram
from cell_curator.provenance import atomic_publish
from cell_curator.refinement import _stored_evaluation


def _model(config: dict) -> CellCuratorConfig:
    value = yaml.safe_load(yaml.safe_dump(config))
    value["run"]["mode"] = "annotation"
    value["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": "markers.json", "version": "fixture"}
    ]
    value["input"]["assays"]["transcriptome"]["representations"]["logcounts"][
        "ranking_backend"
    ] = "scipy_mannwhitney"
    value["markers"].update(
        {
            "rapids_singlecell_version": "",
            "rapids_distribution": "",
            "cuda_major": None,
            "allowed_gpu_contracts": [],
        }
    )
    return CellCuratorConfig.model_validate(value)


def test_bottom_up_bh_correction_covers_all_features(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pvalues = iter([0.01, *([0.5] * 99)] * 2)

    class Result:
        def __init__(self, pvalue: float) -> None:
            self.pvalue = pvalue

    monkeypatch.setattr(
        "cell_curator.evidence.stats.mannwhitneyu",
        lambda *_args, **_kwargs: Result(next(pvalues)),
    )
    matrix = np.zeros((8, 100), dtype=float)
    matrix[:4, 0] = 10
    result = compute_bottom_up_markers(
        matrix,
        pd.Series(["target"] * 4 + ["reference"] * 4),
        pd.Index([f"g{index:03d}" for index in range(100)]),
        detection=matrix,
        top_n=1,
    )
    target = result[result["cluster_id"].eq("target")].iloc[0]
    assert target["feature"] == "g000"
    assert target["pvalue"] == pytest.approx(0.01)
    assert target["adjusted_pvalue"] == pytest.approx(0.5)


def test_missing_detection_and_unmeasured_program_features_force_abstention(
    config: dict,
) -> None:
    model = _model(config)
    capabilities = build_capability_registry(model)
    matrix = np.array([[5.0], [4.0], [0.0], [0.0]])
    groups = pd.Series(["0", "0", "1", "1"])
    program = MarkerProgram(
        label="Target",
        positive=("measured", *tuple(f"missing-{index}" for index in range(99))),
    )
    evidence = score_marker_programs(
        matrix,
        groups,
        pd.Index(["measured"]),
        [program],
        detection=None,
        run_id=model.run.run_id,
        assay="transcriptome",
        representation="logcounts",
        capability=capabilities[0],
    )
    assert evidence["detection_valid"].eq(False).all()
    assert evidence["positive_coverage"].isna().all()
    assert evidence["negative_violation"].isna().all()
    assert evidence["program_coverage"].eq(0.01).all()
    _, decisions = aggregate_multimodal_evidence(
        evidence,
        capabilities,
        model,
        scope="global",
        cluster_ids=["0", "1"],
    )
    assert all(item.abstained and item.cell_type == "Unknown" for item in decisions)


def test_unmeasured_expected_negative_forces_abstention(config: dict) -> None:
    model = _model(config)
    capabilities = build_capability_registry(model)
    matrix = np.array([[5.0], [4.0], [0.0], [0.0]])
    groups = pd.Series(["0", "0", "1", "1"])
    program = MarkerProgram(
        label="Target",
        positive=("measured",),
        negative=("unmeasured-negative",),
    )
    evidence = score_marker_programs(
        matrix,
        groups,
        pd.Index(["measured"]),
        [program],
        detection=matrix,
        run_id=model.run.run_id,
        assay="transcriptome",
        representation="logcounts",
        capability=capabilities[0],
    )
    assert evidence["negative_violation"].isna().all()
    _, decisions = aggregate_multimodal_evidence(
        evidence,
        capabilities,
        model,
        scope="global",
        cluster_ids=["0", "1"],
    )
    assert all(item.abstained and item.cell_type == "Unknown" for item in decisions)


def test_config_rejects_missing_comparator_and_duplicate_lane_pair(config: dict) -> None:
    missing_comparator = yaml.safe_load(yaml.safe_dump(config))
    missing_comparator["comparators"]["roles"] = ["parent_remainder"]
    with pytest.raises(ValueError, match="closest_sibling"):
        CellCuratorConfig.model_validate(missing_comparator)

    duplicated = yaml.safe_load(yaml.safe_dump(config))
    duplicated["evidence"]["lanes"].append(
        {
            "assay": "transcriptome",
            "representation": "logcounts",
            "lane_key": "duplicate-key",
        }
    )
    with pytest.raises(ValueError, match="pair is duplicated"):
        CellCuratorConfig.model_validate(duplicated)


def test_unsafe_parent_identifiers_have_distinct_reversible_path_components() -> None:
    left = safe_path_component("parent/a")
    right = safe_path_component("parent?a")
    assert left != right
    assert safe_path_component("/") != safe_path_component("u-Lw")
    assert left != safe_path_component("u-cGFyZW50L2E")
    assert "/" not in left
    assert "?" not in right
    assert safe_path_component("already-safe") == "already-safe"
    with pytest.raises(ContractError, match="safe path component"):
        require_safe_path_component("../../outside", field="hierarchy level")


def test_raw_representation_uses_raw_feature_axis() -> None:
    query = ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(
            {"gene_symbol": ["current-a", "current-b"]},
            index=["current-a", "current-b"],
        ),
    )
    query.raw = ad.AnnData(
        X=np.ones((2, 3)),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(
            {"gene_symbol": ["raw-a", "raw-b", "raw-c"]},
            index=["raw-a", "raw-b", "raw-c"],
        ),
    )
    assert list(
        feature_ids_for(query, "gene_symbol", declaration={"source": "raw"})
    ) == ["raw-a", "raw-b", "raw-c"]
    assert list(feature_ids_for(query, "gene_symbol", declaration={"source": "X"})) == [
        "current-a",
        "current-b",
    ]


def test_evidence_feature_ids_apply_declared_gene_mapping(tmp_path) -> None:
    query = ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(index=["ENSG000001.1", "ENSG000002.4"]),
    )
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text("ensembl_id\tgene_symbol\nENSG000001\tA\nENSG000002\tB\n")
    assert list(
        feature_ids_for(
            query,
            "",
            declaration={"source": "X"},
            gene_mapping_path=mapping,
        )
    ) == ["A", "B"]


def test_gpu_marker_lane_rejects_one_group_before_claiming_methods() -> None:
    with pytest.raises(ContractError, match="no marker methods were executed"):
        compute_gpu_bottom_up_markers(
            np.ones((2, 1)),
            pd.Series(["only", "only"]),
            pd.Index(["gene"]),
            detection=None,
            rsc=object(),
            marker_config={},
        )


def test_stored_resolution_excludes_transient_children(config: dict) -> None:
    value = yaml.safe_load(yaml.safe_dump(config))
    value["input"]["keys"]["stored_resolutions"] = ["stored_low", "stored_high"]
    model = _model(value)
    parent_obs = pd.DataFrame(
        {
            "stored_low": ["A"] * 6 + ["B"] * 6,
            "stored_high": ["A"] * 5 + ["transient"] + ["B"] * 6,
        },
        index=[f"cell-{index}" for index in range(12)],
    )
    embedding = np.vstack(
        [np.column_stack((np.arange(6), np.zeros(6))),
         np.column_stack((np.arange(6), np.ones(6) * 10))]
    )
    evaluation = _stored_evaluation(parent_obs, embedding, model)
    assert evaluation is not None
    assert evaluation.n_clusters == 2
    assert evaluation.eligible_children is not None
    selected_names = {
        str(pd.Categorical(parent_obs["stored_high"]).categories[index])
        for index in evaluation.eligible_children
    }
    assert selected_names == {"A", "B"}


def test_immutable_atomic_publisher_resumes_identical_and_refuses_overwrite(
    tmp_path,
) -> None:
    path = tmp_path / "artifact.txt"
    stale = tmp_path / ".artifact.txt.partial.previous-job"
    stale.write_text("incomplete\n")
    first = atomic_publish(path, b"complete\n")
    assert atomic_publish(path, b"complete\n") == first
    with pytest.raises(ContractError, match="non-identical"):
        atomic_publish(path, b"replacement\n")
    assert path.read_text() == "complete\n"
    assert stale.read_text() == "incomplete\n"

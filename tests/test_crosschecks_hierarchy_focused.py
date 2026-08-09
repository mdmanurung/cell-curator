from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import cell_curator.crosschecks as crosschecks_module
import cell_curator.environment as environment_module
from cell_curator.contracts import ContractError
from cell_curator.crosschecks import (
    discover_celltypist_models,
    existing_prediction_hypothesis,
    reference_knn_hypothesis,
    scimilarity_hypothesis,
)
from cell_curator.environment import GPUValidationError
from cell_curator.hierarchy import (
    advisory_confidence,
    assign_parent_scoped_level,
    classify_parallel_state,
    recursive_parent_scoped_assignment,
    validate_hierarchy_nesting,
)


def _query_and_reference() -> tuple[ad.AnnData, ad.AnnData]:
    genes = ["G1", "G2", "G3"]
    reference = ad.AnnData(
        X=np.asarray(
            [
                [5.0, 0.0, 0.1],
                [4.5, 0.5, 0.2],
                [0.0, 5.0, 0.1],
                [0.5, 4.5, 0.2],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(
            {"cell_type": ["TypeA", "TypeA", "TypeB", "TypeB"]},
            index=["r0", "r1", "r2", "r3"],
        ),
        var=pd.DataFrame(index=genes),
    )
    query = ad.AnnData(
        X=np.asarray(
            [
                [4.8, 0.2, 0.1],
                [2.5, 2.5, 0.15],
                [0.2, 4.8, 0.1],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(
            {"cluster": ["left", "middle", "right"]},
            index=["q0", "q1", "q2"],
        ),
        var=pd.DataFrame(index=genes),
    )
    return query, reference


def _high_evidence_scores(
    cluster_id: str,
    winner: str,
    runner_up: str,
    *,
    winner_score: float = 0.95,
    runner_score: float = 0.10,
) -> list[dict[str, object]]:
    shared = {
        "specificity": 0.95,
        "agreement": 0.90,
        "coverage": 0.95,
    }
    return [
        {"cluster_id": cluster_id, "label": winner, "score": winner_score, **shared},
        {"cluster_id": cluster_id, "label": runner_up, "score": runner_score, **shared},
    ]


def test_existing_predictions_remain_cluster_hypotheses_with_uncertainty(
    tmp_path,
) -> None:
    adata = ad.AnnData(
        X=np.ones((4, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {
                "cluster": ["0", "0", "1", "1"],
                "prior_prediction": ["T cell", "T cell", "B cell", "Myeloid"],
            },
            index=["c0", "c1", "c2", "c3"],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )

    cells, clusters, status = existing_prediction_hypothesis(
        adata,
        prediction_column="prior_prediction",
        cluster_column="cluster",
        output_dir=tmp_path / "existing",
    )

    assert cells["predicted_label"].tolist() == ["T cell", "T cell", "B cell", "Myeloid"]
    cluster_zero = clusters.set_index("cluster_id").loc["0"]
    cluster_one = clusters.set_index("cluster_id").loc["1"]
    assert cluster_zero["support"] == pytest.approx(1.0)
    assert cluster_zero["uncertainty"] == pytest.approx(0.0)
    assert cluster_one["support"] == pytest.approx(0.5)
    assert cluster_one["uncertainty"] == pytest.approx(1.0)
    assert set(clusters["role"]) == {"hypothesis_not_verdict"}
    assert status["status"] == "complete"
    assert status["automated_output_role"] == "hypothesis_not_verdict"
    receipt = json.loads((tmp_path / "existing" / "provider_status.json").read_text())
    assert receipt["provider"] == "existing_prediction"


def test_cpu_reference_knn_exports_cell_and_cluster_uncertainty() -> None:
    query, reference = _query_and_reference()

    cells, clusters, status = reference_knn_hypothesis(
        query,
        reference,
        reference_label_column="cell_type",
        query_cluster_column="cluster",
        backend="cpu",
        k=3,
        components=2,
    )

    assert cells["barcode"].tolist() == ["q0", "q1", "q2"]
    assert set(cells["backend"]) == {"cpu"}
    assert set(cells["role"]) == {"hypothesis_not_verdict"}
    assert cells["uncertainty"].between(0.0, 1.0).all()
    assert cells.loc[cells["barcode"].eq("q1"), "uncertainty"].iloc[0] > 0.0
    assert set(clusters["cluster_id"]) == {"left", "middle", "right"}
    assert clusters["uncertainty"].between(0.0, 1.0).all()
    assert status["status"] == "complete"
    assert status["backend"] == "cpu"
    assert status["n_shared_genes"] == 3
    assert status["uncertainty_semantics"] == "one_minus_inverse_distance_vote_fraction"


def test_optional_providers_return_actionable_unavailable_status(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    monkeypatch.setattr(crosschecks_module.importlib.util, "find_spec", lambda _name: None)

    celltypist = discover_celltypist_models(output_dir=tmp_path / "celltypist")
    cells, scimilarity = scimilarity_hypothesis(
        ad.AnnData(
            X=np.ones((1, 2), dtype=np.float32),
            obs=pd.DataFrame(index=["c0"]),
            var=pd.DataFrame(index=["G1", "G2"]),
        ),
        model_path=tmp_path / "missing-model",
        output_dir=tmp_path / "scimilarity",
    )

    assert cells.empty
    for provider_status, package in (
        (celltypist, "celltypist"),
        (scimilarity, "scimilarity"),
    ):
        assert provider_status["status"] == "unavailable"
        assert package in provider_status["reason"]
        assert "Install cell-curator[" in provider_status["remediation"]
        assert provider_status["automated_output_role"] == "hypothesis_not_verdict"


def test_gpu_reference_knn_fails_closed_before_cpu_substitution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def reject_gpu() -> None:
        raise GPUValidationError(
            "fixture has no validated real CUDA GPU; CPU substitution is disabled"
        )

    monkeypatch.setattr(environment_module, "require_real_gpu", reject_gpu)
    monkeypatch.setattr(
        crosschecks_module.importlib.util,
        "find_spec",
        lambda _name: pytest.fail("cuML discovery ran before real-GPU validation"),
    )
    monkeypatch.setattr(
        crosschecks_module,
        "harmonize_reference_genes",
        lambda *_args, **_kwargs: pytest.fail("GPU failure read or transformed input data"),
    )
    monkeypatch.setattr(
        crosschecks_module,
        "cKDTree",
        lambda *_args, **_kwargs: pytest.fail("GPU route silently invoked the CPU kNN mapper"),
    )

    with pytest.raises(GPUValidationError, match="CPU substitution is disabled"):
        reference_knn_hypothesis(
            object(),
            object(),
            reference_label_column="cell_type",
            backend="gpu",
            k=3,
            components=2,
        )


def test_gpu_reference_knn_requires_cuml_before_data_access(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(environment_module, "require_real_gpu", lambda: None)
    monkeypatch.setattr(crosschecks_module.importlib.util, "find_spec", lambda _name: None)
    monkeypatch.setattr(
        crosschecks_module,
        "harmonize_reference_genes",
        lambda *_args, **_kwargs: pytest.fail("missing cuML route accessed input data"),
    )

    with pytest.raises(ContractError) as error:
        reference_knn_hypothesis(
            object(),
            object(),
            reference_label_column="cell_type",
            backend="gpu",
        )
    assert str(error.value) == (
        "GPU kNN was requested and a real GPU was detected, but cuML is unavailable; "
        "install a CUDA-compatible cuML build in the declared GPU environment. "
        "CPU substitution is prohibited"
    )


@pytest.mark.parametrize(
    ("labels", "message"),
    [
        ([], "reference kNN requires at least one labeled reference cell"),
        (
            ["TypeA", None],
            "reference label column 'cell_type' contains 1 missing or blank labels",
        ),
        (
            ["TypeA", "  "],
            "reference label column 'cell_type' contains 1 missing or blank labels",
        ),
    ],
)
def test_reference_knn_rejects_empty_or_incomplete_reference_labels_before_matrix_work(
    labels: list[str | None],
    message: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    query, _ = _query_and_reference()
    reference = ad.AnnData(
        X=np.ones((len(labels), 3), dtype=np.float32),
        obs=pd.DataFrame(
            {"cell_type": labels},
            index=[f"r{index}" for index in range(len(labels))],
        ),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    monkeypatch.setattr(
        crosschecks_module,
        "harmonize_reference_genes",
        lambda *_args, **_kwargs: pytest.fail("invalid labels reached matrix harmonization"),
    )

    with pytest.raises(ContractError) as error:
        reference_knn_hypothesis(
            query,
            reference,
            reference_label_column="cell_type",
            backend="cpu",
        )
    assert str(error.value) == message


def test_parent_scope_unknown_confidence_and_cell_state_are_independent() -> None:
    vocabulary = pd.DataFrame(
        [
            {"level": "L2", "parent_scope": "Immune", "label": "T", "cl_id": "CL:0000084"},
            {"level": "L2", "parent_scope": "Immune", "label": "B", "cl_id": "CL:0000236"},
            {
                "level": "L2",
                "parent_scope": "Stromal",
                "label": "Endothelial",
                "cl_id": "CL:0000115",
            },
        ]
    )
    scores = pd.DataFrame(
        [
            *_high_evidence_scores("cycling", "T", "B"),
            *_high_evidence_scores("low_quality", "T", "B"),
            {
                "cluster_id": "uncertain",
                "label": "T",
                "score": 0.51,
                "specificity": 0.8,
                "agreement": 0.8,
                "coverage": 0.8,
            },
            {
                "cluster_id": "uncertain",
                "label": "B",
                "score": 0.50,
                "specificity": 0.8,
                "agreement": 0.8,
                "coverage": 0.8,
            },
            {
                "cluster_id": "cycling",
                "label": "Endothelial",
                "score": 10.0,
                "specificity": 1.0,
                "agreement": 1.0,
                "coverage": 1.0,
            },
        ]
    )
    states = {
        "cycling": classify_parallel_state(cycling_fraction=0.9),
        "low_quality": classify_parallel_state(depth_robust_z=-4.0),
    }

    assigned = assign_parent_scoped_level(
        scores,
        vocabulary,
        level="L2",
        parent_scope="Immune",
        state_by_cluster=states,
    ).set_index("cluster_id")

    assert assigned.loc["cycling", "label"] == "T"
    assert assigned.loc["cycling", "cell_state"] == "cycling"
    assert "Endothelial" not in assigned.loc["cycling", "rejected_alternatives"]
    assert assigned.loc["low_quality", "label"] == "Unknown"
    assert assigned.loc["low_quality", "cell_state"] == "low-quality"
    assert assigned.loc["uncertain", "label"] == "Unknown"
    assert "evidence margin below threshold" in assigned.loc["uncertain", "rationale"]
    assert assigned["confidence"].between(0.0, 1.0).all()
    assert advisory_confidence(
        margin=0.8, specificity=0.9, agreement=0.9, coverage=0.9
    ) > advisory_confidence(margin=0.01, specificity=0.2, agreement=0.2, coverage=0.2)


def test_recursive_assignment_retains_parent_without_child_checkpoint() -> None:
    vocabulary = pd.DataFrame(
        [
            {"level": "L1", "parent_scope": "ROOT", "label": "Immune"},
            {"level": "L1", "parent_scope": "ROOT", "label": "Stromal"},
            {"level": "L2", "parent_scope": "Immune", "label": "T"},
            {"level": "L2", "parent_scope": "Immune", "label": "B"},
            {"level": "L2", "parent_scope": "Stromal", "label": "Fibroblast"},
            {"level": "L2", "parent_scope": "Stromal", "label": "Endothelial"},
        ]
    )
    score_tables = {
        ("L1", "ROOT"): pd.DataFrame(
            [
                *_high_evidence_scores("root_immune", "Immune", "Stromal"),
                *_high_evidence_scores("root_stromal", "Stromal", "Immune"),
            ]
        ),
        ("L2", "Immune"): pd.DataFrame(_high_evidence_scores("immune_child", "T", "B")),
    }

    outputs = recursive_parent_scoped_assignment(
        vocabulary=vocabulary,
        score_tables=score_tables,
    )

    assert set(outputs) == {("L1", "ROOT"), ("L2", "Immune")}
    assert outputs[("L2", "Immune")].iloc[0]["parent_scope"] == "Immune"
    assert outputs[("L2", "Immune")].iloc[0]["label"] == "T"
    assert "Stromal" in set(outputs[("L1", "ROOT")]["label"])


def test_strict_nesting_allows_unknown_but_rejects_reused_child_and_duplicate_leaf() -> None:
    immune = pd.DataFrame(
        [{"level": "L2", "parent_scope": "Immune", "cluster_id": "i0", "label": "Shared"}]
    )
    stromal = pd.DataFrame(
        [{"level": "L2", "parent_scope": "Stromal", "cluster_id": "s0", "label": "Shared"}]
    )
    with pytest.raises(ContractError, match="strict hierarchy nesting failed"):
        validate_hierarchy_nesting([immune, stromal])

    stromal.loc[0, "label"] = "Unknown"
    membership = pd.DataFrame({"barcode": ["c0", "c1"], "leaf_id": ["Shared", "Unknown"]})
    result = validate_hierarchy_nesting([immune, stromal], cell_membership=membership)
    assert result["invariants"] == {
        "one_parent_per_label": True,
        "unknown_permitted": True,
        "one_cell_one_leaf": True,
    }

    duplicated = pd.DataFrame({"barcode": ["c0", "c0"], "leaf_id": ["Shared", "Unknown"]})
    with pytest.raises(ContractError, match="duplicate barcodes"):
        validate_hierarchy_nesting([immune, stromal], cell_membership=duplicated)

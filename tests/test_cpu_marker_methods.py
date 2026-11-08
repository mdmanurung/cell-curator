from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import optimize, sparse

from cell_curator.config import load_config
from cell_curator.contracts import ContractError
from cell_curator.evidence import (
    compute_bottom_up_markers,
    run_lane,
    validate_backend_receipt,
)
from cell_curator.models import BackendReceipt, MarkerProgram


@pytest.mark.parametrize("sparse_input", [False, True])
def test_balanced_binary_logistic_coefficients_match_the_independent_scalar_optimum(
    sparse_input: bool,
) -> None:
    # Unequal group sizes test balanced weights; the constant column tests that
    # the unpenalized intercept, rather than a gene, absorbs the baseline.
    matrix = np.column_stack([np.r_[np.ones(2), -np.ones(18)], np.full(20, 7.0)])
    groups = pd.Series(["rare"] * 2 + ["common"] * 18)
    result = compute_bottom_up_markers(
        sparse.csr_matrix(matrix) if sparse_input else matrix,
        groups,
        pd.Index(["signal", "constant"]),
        detection=np.abs(matrix),
        marker_config={"logreg": {"tol": 1e-8}},
    )
    expected = optimize.brentq(
        lambda coefficient: coefficient - 20 / (1 + np.exp(coefficient)), 0, 20
    )
    logistic = result[result["method"].eq("logreg")].set_index(["cluster_id", "feature"])
    assert logistic.loc[("rare", "signal"), "effect"] == pytest.approx(expected, rel=1e-5)
    assert logistic.loc[("common", "signal"), "effect"] == pytest.approx(-expected, rel=1e-5)
    assert logistic.loc[("rare", "constant"), "effect"] == pytest.approx(0, abs=1e-6)
    assert logistic["pvalue"].isna().all()
    assert logistic["adjusted_pvalue"].isna().all()
    assert set(logistic["effect_semantics"]) == {"logistic_coefficient"}
    assert result.attrs["ranking_group_coverage"] == {
        "wilcoxon": ["common", "rare"],
        "logreg": ["common", "rare"],
    }
    assert result.attrs["runtime"]["logreg"]["fits"]["rare"]["intercept"] == pytest.approx(
        0, abs=1e-6
    )


def test_multiclass_markers_have_independent_binary_fits_and_repeat_exactly() -> None:
    matrix = np.repeat(np.eye(3) * 5, 8, axis=0)
    groups = pd.Series(np.repeat(["A", "B", "C"], 8))
    arguments = (matrix, groups, pd.Index(["a", "b", "c"]))
    first = compute_bottom_up_markers(*arguments, detection=matrix, top_n=1)
    second = compute_bottom_up_markers(*arguments, detection=matrix, top_n=1)
    pd.testing.assert_frame_equal(first, second)
    assert len(first) == 6
    assert first["method_concordant"].all()
    for label, feature in zip(["A", "B", "C"], ["a", "b", "c"], strict=True):
        selected = first[first["cluster_id"].eq(label)]
        assert set(selected["method"]) == {"wilcoxon", "logreg"}
        assert set(selected["feature"]) == {feature}
        assert selected["effect"].gt(0).all()


def test_constant_panel_does_not_create_markers() -> None:
    matrix = np.full((12, 2), 3.0)
    result = compute_bottom_up_markers(
        matrix, pd.Series(["A"] * 6 + ["B"] * 6), pd.Index(["x", "y"]), detection=matrix
    )
    assert result["effect"].eq(0).all()
    assert not result["method_concordant"].any()
    assert result.loc[result["method"].eq("wilcoxon"), "adjusted_pvalue"].eq(1).all()


@pytest.mark.parametrize(
    ("matrix", "groups", "detection", "options", "message"),
    [
        (np.ones((3, 2)), ["A"] * 3, None, {}, "at least two"),
        (np.array([[np.nan, 1], [0, 0]]), ["A", "B"], None, {}, "finite"),
        (np.ones((2, 2)), ["A", None], None, {}, "missing or blank"),
        (np.ones((2, 2)), ["A", "B"], -np.ones((2, 2)), {}, "nonnegative"),
        (np.ones((2, 2)), ["A", "B"], None, {"logreg": {"C": 0}}, "C must"),
        (np.ones((2, 2)), ["A", "B"], None, {"logreg": {"class_weight": "none"}}, "balanced"),
        (
            np.ones((2, 2)),
            ["A", "B"],
            None,
            {"wilcoxon": {"tie_correct": False}},
            "tie correction",
        ),
    ],
)
def test_invalid_inputs_fail_before_publishing_method_evidence(
    matrix: np.ndarray, groups: list, detection: np.ndarray | None, options: dict, message: str
) -> None:
    with pytest.raises(ContractError, match=message):
        compute_bottom_up_markers(
            matrix,
            pd.Series(groups),
            pd.Index(["x", "y"]),
            detection=detection,
            marker_config=options,
        )


def test_failed_logistic_fit_cannot_publish_a_wilcoxon_only_success(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def failed(_objective, initial, **_kwargs):
        return optimize.OptimizeResult(
            success=False, x=initial, fun=1.0, nit=1, message="iteration limit"
        )

    monkeypatch.setattr("cell_curator.evidence.optimize.minimize", failed)
    matrix = np.eye(2)
    with pytest.raises(ContractError, match="logistic regression did not converge"):
        compute_bottom_up_markers(
            matrix, pd.Series(["A", "B"]), pd.Index(["x", "y"]), detection=matrix
        )


@pytest.mark.parametrize("single_group", [False, True])
def test_cpu_lane_receipt_binds_executed_methods_and_rejects_incomplete_resume(
    tmp_path: Path,
    config: dict,
    single_group: bool,
) -> None:
    matrix = np.repeat(np.eye(2) * 5, 8, axis=0)
    source = tmp_path / "input.h5ad"
    adata = ad.AnnData(
        matrix,
        obs=pd.DataFrame(
            {
                "canonical_cluster": ["A"] * 16 if single_group else ["A"] * 8 + ["B"] * 8,
                "lineage": ["immune"] * 16,
            },
            index=[f"c{i}" for i in range(16)],
        ),
        var=pd.DataFrame({"gene_symbol": ["x", "y"]}, index=["x", "y"]),
    )
    adata.layers["logcounts"] = matrix.copy()
    adata.layers["counts"] = matrix.copy()
    adata.write_h5ad(source)
    config["input"]["path"] = str(source)
    representation = config["input"]["assays"]["transcriptome"]["representations"]["logcounts"]
    representation["ranking_backend"] = "scipy_mannwhitney"
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    model = load_config(config_path)
    root = Path(model.run.output_root) / model.run.run_id
    result = run_lane(
        config=model,
        config_path=config_path,
        run_root=root,
        assay_id="transcriptome",
        representation_name="logcounts",
        programs=[MarkerProgram(label="A", positive=("x",), negative=("y",))],
    )
    receipt = BackendReceipt.model_validate_json(result.receipt_path.read_text())
    if single_group:
        assert result.bottom_up.empty
        assert receipt.marker_methods == ()
        assert receipt.convergence_validated is None
        assert receipt.runtime["marker_status"] == "unavailable_no_comparator"
    else:
        assert receipt.marker_methods == ("wilcoxon", "logreg")
        assert receipt.convergence_validated is True
        assert receipt.runtime["ranking_group_coverage"] == {
            "wilcoxon": ["A", "B"],
            "logreg": ["A", "B"],
        }
        assert receipt.runtime["logreg"]["fit_intercept"] is True
    value = json.loads(result.receipt_path.read_text())
    value["marker_methods"] = ["wilcoxon"]
    result.receipt_path.write_text(json.dumps(value))
    message = "cannot claim" if single_group else "Wilcoxon and balanced logreg"
    with pytest.raises(ContractError, match=message):
        run_lane(
            config=model,
            config_path=config_path,
            run_root=root,
            assay_id="transcriptome",
            representation_name="logcounts",
            programs=[],
        )


def test_cpu_receipt_requires_successful_logistic_convergence(config: dict) -> None:
    from cell_curator.capabilities import build_capability_registry
    from cell_curator.config import CellCuratorConfig

    config["input"]["assays"]["transcriptome"]["representations"]["logcounts"][
        "ranking_backend"
    ] = "scipy_mannwhitney"
    capability = build_capability_registry(CellCuratorConfig.model_validate(config))[0]
    receipt = BackendReceipt(
        run_id="fixture",
        lane_key="cpu",
        backend="scipy_mannwhitney",
        capability_level=capability.level,
        status="complete",
        input_sha256="a" * 64,
        config_sha256="b" * 64,
        marker_methods=("wilcoxon", "logreg"),
        convergence_validated=False,
    )
    with pytest.raises(ContractError, match="convergence"):
        validate_backend_receipt(receipt, capability)


def test_standalone_marker_artifacts_preserve_configured_method_runtime(tmp_path: Path) -> None:
    from cell_curator.markers import run_marker_analysis, write_marker_artifacts

    matrix = np.repeat(np.eye(2) * 5, 8, axis=0)
    analysis = run_marker_analysis(
        matrix,
        pd.Series(["A"] * 8 + ["B"] * 8),
        pd.Index(["a", "b"]),
        [MarkerProgram(label="A", positive=("a",), negative=("b",))],
        detection=matrix,
        marker_config={"logreg": {"C": 0.3}},
    )
    artifacts = write_marker_artifacts(analysis, tmp_path / "markers")
    manifest = json.loads(artifacts["manifest"].read_text())
    assert manifest["marker_methods"] == ["logreg", "wilcoxon"]
    assert manifest["marker_runtime"]["logreg"]["C"] == 0.3
    assert manifest["marker_runtime"]["logreg"]["fits"]["A"]["converged"] is True

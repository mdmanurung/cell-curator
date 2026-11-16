from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import optimize, sparse

from cell_curator.contracts import ContractError
from cell_curator.evidence import compute_bottom_up_markers
from cell_curator.models import MarkerProgram

BOTH = ("wilcoxon", "logreg")  # decision row 18: logistic regression is opt-in


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
        methods=BOTH,
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
    first = compute_bottom_up_markers(*arguments, detection=matrix, top_n=1, methods=BOTH)
    second = compute_bottom_up_markers(*arguments, detection=matrix, top_n=1, methods=BOTH)
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
        matrix, pd.Series(["A"] * 6 + ["B"] * 6), pd.Index(["x", "y"]), detection=matrix,
        methods=BOTH,
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
            methods=BOTH,
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
            matrix, pd.Series(["A", "B"]), pd.Index(["x", "y"]), detection=matrix, methods=BOTH
        )






def test_standalone_marker_artifacts_record_the_default_wilcoxon_runtime(tmp_path: Path) -> None:
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
    assert manifest["marker_methods"] == ["wilcoxon"]
    assert "logreg" not in manifest["marker_runtime"]


def test_configured_logistic_runtime_is_kept_when_the_method_is_requested() -> None:
    matrix = np.repeat(np.eye(2) * 5, 8, axis=0)
    result = compute_bottom_up_markers(
        matrix, pd.Series(["A"] * 8 + ["B"] * 8), pd.Index(["a", "b"]), detection=matrix,
        methods=BOTH, marker_config={"logreg": {"C": 0.3}},
    )
    assert result.attrs["runtime"]["logreg"]["C"] == 0.3
    assert result.attrs["runtime"]["logreg"]["fits"]["A"]["converged"] is True


def test_wilcoxon_is_the_default_and_fits_no_logistic_model() -> None:
    matrix = np.repeat(np.eye(2) * 5, 8, axis=0)
    result = compute_bottom_up_markers(
        matrix, pd.Series(["A"] * 8 + ["B"] * 8), pd.Index(["a", "b"]), detection=matrix
    )
    assert set(result["method"]) == {"wilcoxon"}
    assert "logreg" not in result.attrs["runtime"]

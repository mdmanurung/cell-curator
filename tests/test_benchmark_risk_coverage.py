"""A risk-coverage curve, because one operating point cannot carry a
selective-prediction claim.

`evaluate_predictions` reports coverage, selective accuracy, and risk at exactly
one point, so the package's own headline claim -- comparable accuracy at matched
coverage, with abstention doing real work -- had no artifact behind it. The
single operating point it did report is the `Unknown`-label coverage, which
measures a different thing entirely: a system can rank its own confidence well
while abstaining arbitrarily, or abstain well while ranking noise. Reading one
number off that axis and calling it the selective-prediction result is the same
class of error as publishing an alias.

These tests pin the boundary semantics (the sweep is worthless if its endpoints
do not mean "everything" and "only maximal confidence"), the refusal to
fabricate a point where nothing is called, and the failure mode for a table that
carries no confidence at all.
"""

from __future__ import annotations

import pandas as pd
import pytest
from conftest import REPO_ROOT

import cell_curator
from cell_curator.api import sweep_benchmark_risk_coverage
from cell_curator.benchmark import (
    DEFAULT_RISK_COVERAGE_THRESHOLDS,
    evaluate_predictions,
    risk_coverage_sweep,
)
from cell_curator.contracts import ContractError

FIXTURE = REPO_ROOT / "benchmarks" / "fixtures" / "synthetic_predictions.tsv"


def _fixture_frame() -> pd.DataFrame:
    return pd.read_csv(FIXTURE, sep="\t", keep_default_na=False)


def test_the_default_grid_spans_exactly_zero_to_one() -> None:
    """Endpoints have to be exact, because their meaning is a promise.

    A grid that started at 0.05 would silently refuse to report full coverage, so
    a table whose lowest confidence is 0.03 would be scored as though it had
    abstained on those cells. A grid that stopped at 0.95 would never report the
    maximal-confidence point. Both are the kind of off-by-one that only shows up
    as a slightly flattering curve.
    """

    grid = DEFAULT_RISK_COVERAGE_THRESHOLDS
    assert grid[0] == 0.0
    assert grid[-1] == 1.0
    assert all(0.0 <= value <= 1.0 for value in grid)
    assert all(
        left < right for left, right in zip(grid[:-1], grid[1:], strict=True)
    )
    assert len(grid) == 21


def test_threshold_zero_calls_everything_and_reproduces_plain_accuracy() -> None:
    """The two endpoints of the instrument have to agree with the scalar metrics.

    At 0.0 the sweep's selective accuracy is the table's plain accuracy over
    every row, and its risk is one minus that. A sweep whose full-coverage point
    disagreed with `accuracy` would be scoring a different notion of correctness,
    which is the fastest way to make a curve that flatters one system and
    penalises another.
    """

    frame = _fixture_frame()
    measured = evaluate_predictions(frame, system="fixture", dataset="synthetic")
    points = risk_coverage_sweep(frame)
    full = points[0]
    assert full.threshold == 0.0
    assert full.n_called == len(frame) == 8
    assert full.coverage == 1.0
    assert full.selective_accuracy == measured.accuracy == 0.875
    assert full.risk == pytest.approx(1 - measured.accuracy)
    assert full.n_correct == 7


def test_threshold_one_calls_only_cells_claiming_maximal_confidence() -> None:
    """Inclusive at both ends: a cell at 1.0 is called at threshold 1.0.

    Exclusive at the top would make the last point unreachable for any table that
    does not emit exact 1.0, which is most of them, and the curve would end one
    step early. The fixture's highest confidence is 0.95, so its 1.0 point calls
    nothing -- and that has to read as unmeasured, not as a system that called
    every cell wrong.
    """

    frame = _fixture_frame()
    assert frame["confidence"].max() == 0.95
    empty = risk_coverage_sweep(frame, thresholds=[1.0])[0]
    assert empty.threshold == 1.0
    assert empty.n_called == 0
    assert empty.coverage == 0.0
    assert empty.selective_accuracy is None
    assert empty.risk is None

    at_max = _frame(["A", "B", "A"], [0.5, 1.0, 1.0], ["A", "B", "B"])
    called = risk_coverage_sweep(at_max, thresholds=[1.0])[0]
    assert called.n_called == 2
    assert called.n_correct == 1
    assert called.coverage == pytest.approx(2 / 3)
    assert called.selective_accuracy == 0.5
    assert called.risk == 0.5

    below_max = risk_coverage_sweep(at_max, thresholds=[0.99])[0]
    assert below_max.n_called == 2
    everything = risk_coverage_sweep(at_max, thresholds=[0.0])[0]
    assert everything.n_called == 3
    assert everything.selective_accuracy == pytest.approx(2 / 3)


def test_a_point_where_nothing_is_called_reports_no_risk() -> None:
    """An empty call set has no measured risk, and must not be given one.

    Reporting 0.0 accuracy would read as a system that abstained perfectly and
    reporting 1.0 risk would read as a system that abstained disastrously. Both
    are fabrications, and a curve that contains either would make a comparison
    at matched coverage meaningless exactly where comparisons are hardest.
    """

    frame = _frame(["A", "A"], [0.1, 0.2], ["A", "A"])
    points = {point.threshold: point for point in risk_coverage_sweep(frame)}
    assert points[0.5].n_called == 0
    assert points[0.5].coverage == 0.0
    assert points[0.5].selective_accuracy is None
    assert points[0.5].risk is None
    assert points[0.1].n_called == 2
    assert points[0.1].risk == 0.0


def test_the_curve_separates_a_system_that_ranks_its_own_uncertainty() -> None:
    """The scientific content: selective accuracy above plain accuracy.

    The fixture's single error is its lowest-confidence row, so raising the
    threshold buys perfect accuracy at 87.5% coverage. This is the whole reason
    to publish a curve, and it is invisible in `BenchmarkResult`, whose coverage
    counts `Unknown` labels instead and therefore reports a different operating
    point on a different axis.
    """

    frame = _fixture_frame()
    measured = evaluate_predictions(frame, system="fixture", dataset="synthetic")
    points = {point.threshold: point for point in risk_coverage_sweep(frame)}
    selective = points[0.5]
    assert selective.coverage == 0.875
    assert selective.selective_accuracy == 1.0
    assert selective.risk == 0.0
    assert selective.selective_accuracy > measured.accuracy
    assert selective.coverage != measured.coverage
    coverages = [point.coverage for point in risk_coverage_sweep(frame)]
    assert coverages == sorted(coverages, reverse=True)
    assert coverages[0] == 1.0
    assert coverages[-1] == 0.0


def test_a_table_without_confidence_is_refused_rather_than_explained() -> None:
    """A missing confidence column is a broken table, not a result to interpret.

    Returning a reason string or an empty curve was the alternative, and it is
    the wrong one: at the call site an empty curve is indistinguishable from a
    table that abstains on every cell, and a reason string in a return type that
    is otherwise a list of points is a union nobody will type-check. Raising
    keeps the return type a measurement or nothing.
    """

    frame = pd.DataFrame(
        {
            "barcode": ["a", "b"],
            "true_label": ["A", "B"],
            "predicted_label": ["A", "A"],
        }
    )
    with pytest.raises(ContractError, match="confidence"):
        risk_coverage_sweep(frame)


def test_the_sweep_shares_the_calibration_guard_on_confidence() -> None:
    """One definition of a usable confidence, or two instruments to reconcile.

    A value outside [0, 1] is refused with the message the calibration figures
    already use, so a table cannot be scored by one and rejected by the other.
    """

    with pytest.raises(ContractError, match=r"within \[0, 1\]"):
        risk_coverage_sweep(_frame(["A"], [1.5], ["A"]))
    with pytest.raises(ContractError, match=r"within \[0, 1\]"):
        risk_coverage_sweep(_frame(["A"], [float("nan")], ["A"]))


def test_the_sweep_needs_the_two_label_columns_and_nothing_else() -> None:
    """A sweep that demanded a barcode would be scoring a table, not a system.

    The three columns it reads are the three it uses, so a table carrying a
    duplicated barcode or an unexpected extra column is neither rejected nor
    silently reinterpreted: the one-row-per-barcode guarantee belongs to
    `evaluate_predictions`, which is where a malformed export is reported.
    """

    frame = pd.DataFrame(
        {
            "true_label": ["A", "A"],
            "predicted_label": ["A", "B"],
            "confidence": [0.9, 0.8],
        }
    )
    point = risk_coverage_sweep(frame, thresholds=[0.0])[0]
    assert point.n_total == 2
    assert point.selective_accuracy == 0.5
    with pytest.raises(ContractError, match="true_label"):
        risk_coverage_sweep(frame.drop(columns=["true_label"]))


@pytest.mark.parametrize(
    "thresholds",
    [
        [1.0, 0.5],
        [0.0, 0.0],
        [0.0, 1.5],
        [-0.1],
        [float("nan")],
        [],
    ],
)
def test_a_supplied_grid_must_have_usable_endpoints(thresholds: list[float]) -> None:
    """Order and range are part of what a curve means, not a formatting detail.

    Silently sorting an unsorted grid would make two reports of the same table
    look identical when they are not, and accepting an out-of-range threshold
    would produce a point whose "coverage" is a constant 0.0 or 1.0 that a reader
    could mistake for a measurement.
    """

    with pytest.raises(ContractError, match="risk-coverage thresholds"):
        risk_coverage_sweep(_frame(["A"], [0.9], ["A"]), thresholds=thresholds)


def test_the_sweep_is_repeatable_for_one_table() -> None:
    """Two reads of the same table must not produce two different curves.

    A curve that moved between runs would be quotable evidence, so the tie-break
    on the deepest shared ontology ancestor and the fixed grid order matter as
    much as the arithmetic.
    """

    frame = _fixture_frame()
    first = risk_coverage_sweep(frame)
    second = risk_coverage_sweep(frame)
    assert [point.model_dump(mode="json") for point in first] == [
        point.model_dump(mode="json") for point in second
    ]


def test_a_two_point_grid_is_honoured_exactly() -> None:
    """A caller-supplied grid is used as given, not merged with the default.

    Merging would silently insert points a caller did not ask for and make the
    endpoints of their curve lie.
    """

    frame = _fixture_frame()
    points = risk_coverage_sweep(frame, thresholds=[0.0, 1.0])
    assert [point.threshold for point in points] == [0.0, 1.0]
    assert len(points) == 2


def test_the_public_entry_point_reads_a_table_from_disk() -> None:
    """The curve is public because a comparative claim is checked from outside.

    Whoever holds the competing export is the only person who can compare two
    systems at matched coverage, so the exported function takes a table path and
    returns the curve as records. It publishes nothing, because a published
    result is a commitment and a curve is a reading.
    """

    curve = sweep_benchmark_risk_coverage(FIXTURE, thresholds=[0.0, 0.5])
    assert [point.threshold for point in curve] == [0.0, 0.5]
    assert [point.coverage for point in curve] == [1.0, 0.875]
    assert [point.model_dump(mode="json") for point in curve] == [
        point.model_dump(mode="json")
        for point in risk_coverage_sweep(_fixture_frame(), thresholds=[0.0, 0.5])
    ]
    assert "sweep_benchmark_risk_coverage" in cell_curator.__all__


def _frame(
    truth: list[str], confidence: list[float], prediction: list[str]
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "true_label": truth,
            "predicted_label": prediction,
            "confidence": confidence,
        }
    )


def test_abstain_unknown_never_counts_an_abstention_as_a_call() -> None:
    """The fixture has three Unknown predictions; with abstention on, none is called."""

    plain, abstaining = (
        sweep_benchmark_risk_coverage(FIXTURE, thresholds=[0.0], abstain_unknown=flag)[0]
        for flag in (False, True)
    )
    assert plain.n_called == 8
    assert abstaining.n_called == 5
    assert abstaining.coverage == 0.625  # matches BenchmarkResult.coverage for the table
    assert abstaining.selective_accuracy == 1.0


def test_the_sweep_command_publishes_one_immutable_curve(tmp_path) -> None:
    import json

    from cell_curator.cli import build_parser, dispatch

    output = tmp_path / "curve.json"
    arguments = [
        "benchmark", "sweep", "--predictions", str(FIXTURE), "--output", str(output),
        "--thresholds", "0", "0.9", "--abstain-unknown",
    ]  # fmt: skip
    dispatch(build_parser().parse_args(arguments))
    curve = json.loads(output.read_text())
    assert curve["artifact_kind"] == "cell_curator_benchmark_risk_coverage"
    assert curve["prediction_path"] == FIXTURE.name
    assert curve["abstain_unknown"] is True
    assert [point["n_called"] for point in curve["points"]] == [5, 3]
    with pytest.raises(ContractError, match="immutable"):
        dispatch(build_parser().parse_args([*arguments[:-1]]))

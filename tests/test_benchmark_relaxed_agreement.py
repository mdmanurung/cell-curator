"""A correct call at the wrong Cell Ontology depth is not an error.

`accuracy`, `macro_f1`, `micro_f1`, and the risk-coverage sweep all compare
exported label strings for exact equality. The intended benchmark tables carry
ground truth at the annotator's own granularity -- truth is "central memory
CD4-positive, alpha-beta T cell" where a system that found the cell is expected
to answer "T cell" -- so exact equality reads a correct, more general call as a
complete error. Every published number made that way is misleading about the
thing it claims to measure, and a selective-prediction curve built on it could
never approach 1.0 however good the system was, which removes the only form in
which two systems can be compared at matched coverage.

`relaxed_accuracy` credits a prediction whose Cell Ontology term is the true term
or an `is_a` ancestor or descendant of it, and is reported beside `accuracy`
rather than in place of it. That symmetry is the metric's one real caveat and it
is pinned from both sides here: being too general and being too specific are the
same kind of mistake to this number, so a system that answers one blunt term for
everything scores 1.0 and only the strict field beside it says what it gave up.
Collapsing every call onto a shared ancestor would let a system buy apparent
accuracy by declining to be specific, which is the failure the strict numbers
exist to expose, so `accuracy` and its F1 pair are asserted here to be untouched
by any of this.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from conftest import REPO_ROOT

from cell_curator import sweep_benchmark_risk_coverage
from cell_curator.benchmark import evaluate_file, evaluate_predictions, risk_coverage_sweep
from cell_curator.contracts import ContractError
from cell_curator.models import BenchmarkResult
from cell_curator.ontology import CellOntology, OntologyTerm

FIXTURE = REPO_ROOT / "benchmarks" / "fixtures" / "synthetic_predictions.tsv"

CELL = "CL:0000000"
NATIVE_CELL = "CL:0000003"
LEUKOCYTE = "CL:0000738"
LYMPHOCYTE = "CL:0000542"
T_CELL = "CL:0000084"
CENTRAL_MEMORY_CD4 = "CL:0000904"
MONOCYTE = "CL:0000576"
B_CELL = "CL:0000236"
CD4_POSITIVE_T = "CL:0000624"
FIBROBLAST = "CL:0000057"
ABSENT_TERM = "CL:1234567"


def _ontology() -> CellOntology:
    """A miniature Cell Ontology, shaped so a granularity difference exists.

    `CENTRAL_MEMORY_CD4` hangs directly off `T_CELL` where the real ontology
    reaches it through memory T cell, because these tests are about the shape of
    the relaxation and not about the depth of any particular term.
    """

    terms = {
        CELL: OntologyTerm(CELL, "cell"),
        NATIVE_CELL: OntologyTerm(NATIVE_CELL, "native cell", (CELL,)),
        LEUKOCYTE: OntologyTerm(LEUKOCYTE, "leukocyte", (NATIVE_CELL,)),
        LYMPHOCYTE: OntologyTerm(LYMPHOCYTE, "lymphocyte", (LEUKOCYTE,)),
        T_CELL: OntologyTerm(T_CELL, "T cell", (LYMPHOCYTE,)),
        CENTRAL_MEMORY_CD4: OntologyTerm(
            CENTRAL_MEMORY_CD4, "central memory CD4-positive, alpha-beta T cell", (T_CELL,)
        ),
        MONOCYTE: OntologyTerm(MONOCYTE, "monocyte", (LEUKOCYTE,)),
        B_CELL: OntologyTerm(B_CELL, "B cell", (LYMPHOCYTE,)),
        FIBROBLAST: OntologyTerm(FIBROBLAST, "fibroblast", (NATIVE_CELL,)),
    }
    return CellOntology(terms, source="fixture")


def _frame(
    true_labels: list[str],
    predicted_labels: list[str],
    true_ids: list[str],
    predicted_ids: list[str],
    confidence: list[float] | None = None,
) -> pd.DataFrame:
    count = len(true_labels)
    values: dict[str, object] = {
        "barcode": [f"b{index}" for index in range(count)],
        "true_label": true_labels,
        "predicted_label": predicted_labels,
        "true_cl_id": true_ids,
        "predicted_cl_id": predicted_ids,
    }
    if confidence is not None:
        values["confidence"] = confidence
    return pd.DataFrame(values)


def _evaluate(frame: pd.DataFrame) -> BenchmarkResult:
    return evaluate_predictions(
        frame, system="fixture", dataset="relaxed-agreement", ontology=_ontology()
    )


def test_a_correct_but_coarser_call_is_a_complete_error_under_exact_equality() -> None:
    """The defect this metric exists to correct, stated as one row.

    The cell really is a central memory CD4-positive, alpha-beta T cell, the
    system's answer "T cell" is true of it, and `accuracy` scores the row wrong
    along with the F1 pair. If the relaxed field did not read 1.0 here, the only
    thing it could be reading is something other than Cell Ontology ancestry.
    """

    frame = _frame(
        ["central memory CD4-positive, alpha-beta T cell"],
        ["T cell"],
        [CENTRAL_MEMORY_CD4],
        [T_CELL],
    )
    result = _evaluate(frame)
    assert result.accuracy == 0.0
    assert result.macro_f1 == 0.0
    assert result.micro_f1 == 0.0
    assert result.selective_accuracy == 0.0
    assert result.relaxed_accuracy == 1.0


def test_the_relaxation_credits_specificity_and_bluntness_alike() -> None:
    """The asymmetry, pinned from both sides, because it is the whole caveat.

    Reading "central memory CD4-positive, alpha-beta T cell" as truth, a blunter
    `T cell` and a more specific descendant of it are both credited, and so is
    `leukocyte`, which is blunter still. A system that answered "leukocyte" for
    every cell would therefore score 1.0 here. That is the price of the
    relaxation and the reason the strict field beside it is not replaced: the
    gap between the two numbers is the only evidence in the result about how much
    of a system's accuracy is identification and how much is specificity.
    """

    blunt = _frame(["T cell"], ["leukocyte"], [T_CELL], [LEUKOCYTE])
    assert _evaluate(blunt).relaxed_accuracy == 1.0
    assert _evaluate(blunt).accuracy == 0.0

    over_specific = _frame(
        ["T cell"],
        ["central memory CD4-positive, alpha-beta T cell"],
        [T_CELL],
        [CENTRAL_MEMORY_CD4],
    )
    assert _evaluate(over_specific).relaxed_accuracy == 1.0
    assert _evaluate(over_specific).accuracy == 0.0

    identical = _frame(
        ["T cell"],
        ["T cell"],
        [T_CELL],
        [T_CELL],
    )
    assert _evaluate(identical).relaxed_accuracy == 1.0
    assert _evaluate(identical).accuracy == 1.0


def test_a_prediction_in_another_branch_is_never_credited() -> None:
    """Sharing a high ancestor is not agreement, and this is where that bites.

    Every pair of real Cell Ontology terms shares `cell` and `native cell`, so a
    rule that accepted any shared ancestry would credit a monocyte called a T cell
    and then a fibroblast called a T cell. Requiring the shared term to be one of
    the two themselves is what keeps the metric a statement about granularity
    rather than a statement about being a cell at all.
    """

    wrong_branch = _frame(["T cell"], ["monocyte"], [T_CELL], [MONOCYTE])
    assert _evaluate(wrong_branch).relaxed_accuracy == 0.0
    other_tissue = _frame(["T cell"], ["fibroblast"], [T_CELL], [FIBROBLAST])
    assert _evaluate(other_tissue).relaxed_accuracy == 0.0
    fine_against_coarse_branch = _frame(
        ["central memory CD4-positive, alpha-beta T cell"],
        ["B cell"],
        [CENTRAL_MEMORY_CD4],
        [B_CELL],
    )
    assert _evaluate(fine_against_coarse_branch).relaxed_accuracy == 0.0


def test_the_strict_metrics_do_not_move_when_a_relaxed_number_is_published() -> None:
    """The relaxation is additive, and the pin is on every strict field at once.

    Comparing a table with its CL identifier columns against the same table
    without them isolates what the new evidence can touch. A change to
    `accuracy`, either F1, the coverage, or the coverage risk would mean the
    relaxed reading had leaked into a number that other results are already
    compared on, and a change to `relaxed_accuracy` itself would mean the metric
    is not reading the identifiers it is named after.
    """

    frame = _frame(
        ["central memory CD4-positive, alpha-beta T cell", "T cell", "monocyte"],
        ["T cell", "Unknown", "monocyte"],
        [CENTRAL_MEMORY_CD4, T_CELL, MONOCYTE],
        [T_CELL, "", MONOCYTE],
    )
    bare = frame.drop(columns=["true_cl_id", "predicted_cl_id"])
    with_identifiers = _evaluate(frame)
    without_identifiers = _evaluate(bare)
    for field in (
        "accuracy",
        "macro_f1",
        "micro_f1",
        "coverage",
        "selective_accuracy",
        "coverage_risk",
        "hierarchical_accuracy",
        "unknown_precision",
        "unknown_recall",
        "unknown_f1",
        "hierarchy_error",
    ):
        assert getattr(with_identifiers, field) == getattr(without_identifiers, field), field
    assert with_identifiers.relaxed_accuracy == 1.0
    assert without_identifiers.relaxed_accuracy is None


def test_a_table_that_makes_no_ontology_claim_reports_no_relaxed_number() -> None:
    """The synthetic fixture's case: no CL columns, so nothing to relax.

    The strict fields are all real measurements of that table and
    `relaxed_accuracy` is the one field with no evidence behind it. Substituting
    the strict number, or the slash-path `hierarchical_accuracy`, would publish a
    measurement of a different quantity under a name that promises an ontology
    one, which is the defect the sibling `ontology_distance_error` field shipped
    once already.
    """

    frame = pd.read_csv(FIXTURE, sep="\t", keep_default_na=False)
    result = evaluate_predictions(frame, system="fixture", dataset="synthetic")
    assert result.accuracy == 0.875
    assert result.hierarchical_accuracy == pytest.approx(0.9375)
    assert result.relaxed_accuracy is None
    supplied = evaluate_predictions(
        frame, system="fixture", dataset="synthetic", ontology=_ontology()
    )
    assert supplied.relaxed_accuracy is None


def test_ontology_identifiers_without_a_snapshot_fail_closed() -> None:
    """A table claiming to be scored against CL cannot be scored unreported.

    Returning `None` would be indistinguishable in the artifact from a table that
    made no claim, so a caller could ship a relaxed comparison that was never
    performed. The message names the metric and the columns, as the distance
    metric's does, so the reader knows which claim went unchecked.
    """

    frame = _frame(["T cell"], ["T cell"], [T_CELL], [T_CELL])
    with pytest.raises(ContractError, match="relaxed_accuracy"):
        evaluate_predictions(frame, system="fixture", dataset="relaxed-unbacked")
    with pytest.raises(ContractError, match="true_cl_id"):
        evaluate_predictions(frame, system="fixture", dataset="relaxed-unbacked")


def test_an_identifier_the_snapshot_does_not_know_is_refused() -> None:
    """An unresolvable CL id must not be scored as maximally far.

    The alternative reading -- "not in the ontology" means "as wrong as possible"
    -- would let a typo become the worst row of a published relaxed accuracy
    without anything in the result saying the table was never scored.
    """

    frame = _frame(["T cell"], ["T cell"], [T_CELL], [ABSENT_TERM])
    with pytest.raises(ContractError, match="ontology term is unknown"):
        _evaluate(frame)


def test_an_abstention_is_left_out_of_the_relaxed_mean() -> None:
    """A blank predicted term is not a maximally wrong term, in either metric.

    An `Unknown` prediction has no CL term to be near or far from, and this
    instrument scores abstention through the Unknown precision and recall.
    Charging the row as a relaxed error as well would count one mistake twice, so
    the identifiers of such rows are excluded and the mean is taken over the rows
    that carry a pair. A table whose identifiers are all blank therefore has
    nothing to measure and says so.
    """

    frame = _frame(
        ["T cell"] * 3,
        ["lymphocyte", "Unknown", "Unknown"],
        [T_CELL] * 3,
        [LYMPHOCYTE, "", "  "],
    )
    assert _evaluate(frame).relaxed_accuracy == 1.0
    blank = _frame(["T cell"] * 2, ["Unknown", "Unknown"], [T_CELL] * 2, ["", ""])
    assert _evaluate(blank).relaxed_accuracy is None


def test_the_relaxed_curve_measures_granularity_and_the_strict_one_cannot() -> None:
    """The scientific content: confidence ranking against a fine ground truth.

    The most confident cell here is the one answered at a coarser granularity,
    and it is right. Under exact equality that row is an error, so the strict
    curve falls to 0.0 selective accuracy at the top of its range while the three
    exact calls are all less confident: the curve is anti-correlated with quality
    and unusable for the comparison the package claims to support. The relaxed
    curve reads the same table in the intended order.
    """

    frame = _frame(
        [
            "central memory CD4-positive, alpha-beta T cell",
            "T cell",
            "T cell",
            "T cell",
        ],
        ["T cell", "T cell", "T cell", "T cell"],
        [CENTRAL_MEMORY_CD4, T_CELL, T_CELL, T_CELL],
        [T_CELL, T_CELL, T_CELL, T_CELL],
        [1.0, 0.4, 0.35, 0.3],
    )
    strict = {point.threshold: point for point in risk_coverage_sweep(frame)}
    relaxed = {
        point.threshold: point
        for point in risk_coverage_sweep(frame, agreement="relaxed", ontology=_ontology())
    }
    assert strict[1.0].selective_accuracy == 0.0
    assert relaxed[1.0].selective_accuracy == 1.0
    assert strict[0.0].selective_accuracy == 0.75
    assert relaxed[0.0].selective_accuracy == 1.0
    assert relaxed[0.0].selective_accuracy == _evaluate(frame).relaxed_accuracy
    assert strict[0.0].selective_accuracy == _evaluate(frame).accuracy
    assert [point.agreement for point in risk_coverage_sweep(frame)] == ["strict"] * 21
    assert relaxed[0.5].n_called == 1
    assert relaxed[0.5].coverage == 0.25


def test_a_relaxed_curve_keeps_the_strict_boundary_semantics() -> None:
    """Inclusive at both ends, and no point where nothing is called.

    The relaxed variant is a different notion of correctness and nothing else: a
    cell claiming maximal confidence is still called at threshold 1.0, a threshold
    that calls nothing still reports no measured risk, and coverage still falls
    monotonically. Changing the boundary along with the correctness would make
    the two curves incomparable at matched coverage, which is the only comparison
    either of them exists for.
    """

    frame = _frame(
        ["central memory CD4-positive, alpha-beta T cell", "T cell", "T cell"],
        ["T cell", "T cell", "T cell"],
        [CENTRAL_MEMORY_CD4, T_CELL, T_CELL],
        [T_CELL, T_CELL, T_CELL],
        [1.0, 0.9, 0.2],
    )
    points = {
        point.threshold: point
        for point in risk_coverage_sweep(frame, agreement="relaxed", ontology=_ontology())
    }
    assert points[1.0].n_called == 1
    assert points[1.0].selective_accuracy == 1.0
    assert points[0.9].n_called == 2
    assert points[0.9].coverage == pytest.approx(2 / 3)
    assert points[0.3].n_called == 2
    assert points[0.2].n_called == 3
    below_max = _frame(
        ["central memory CD4-positive, alpha-beta T cell", "T cell"],
        ["T cell", "T cell"],
        [CENTRAL_MEMORY_CD4, T_CELL],
        [T_CELL, T_CELL],
        [0.9, 0.8],
    )
    empty = risk_coverage_sweep(
        below_max, thresholds=[1.0], agreement="relaxed", ontology=_ontology()
    )[0]
    assert empty.n_called == 0
    assert empty.coverage == 0.0
    assert empty.selective_accuracy is None
    assert empty.risk is None
    coverages = [
        point.coverage
        for point in risk_coverage_sweep(
            frame, agreement="relaxed", ontology=_ontology()
        )
    ]
    assert coverages == sorted(coverages, reverse=True)


def test_an_unmeasurable_row_is_charged_by_the_curve_and_excluded_from_the_mean() -> None:
    """The one place the two numbers have different denominators, pinned.

    `relaxed_accuracy` drops a row with no predicted term, because an abstention
    is scored by the Unknown precision and recall. A curve cannot: its denominator
    is the set of cells called at a threshold, so dropping the row would report a
    coverage the table does not have. The full-coverage point therefore sits one
    step below the scalar on a table containing an abstention, and saying so here
    is what stops the two from being read as a disagreement later.
    """

    frame = _frame(
        ["T cell"] * 2,
        ["T cell", "Unknown"],
        [T_CELL] * 2,
        [T_CELL, ""],
        [0.9, 0.8],
    )
    result = _evaluate(frame)
    full = risk_coverage_sweep(frame, agreement="relaxed", ontology=_ontology())[0]
    assert result.relaxed_accuracy == 1.0
    assert full.selective_accuracy == 0.5
    assert full.n_correct == 1
    assert full.n_called == 2
    assert result.accuracy == 0.5


def test_a_relaxed_curve_needs_identifiers_and_a_snapshot() -> None:
    """No CL claim, no relaxed curve: the return type has nowhere to say so.

    A curve of points with `None` accuracy would be indistinguishable from a
    table whose confidence abstains everywhere, and that state already means
    something specific. Refusing keeps the return type a measurement or nothing,
    the same choice the sweep already makes about a missing confidence column.
    """

    frame = _frame(
        ["T cell", "T cell"], ["T cell", "T cell"], [T_CELL, T_CELL], [T_CELL, T_CELL]
    )
    frame["confidence"] = [0.9, 0.8]
    without_identifiers = frame.drop(columns=["true_cl_id", "predicted_cl_id"])
    with pytest.raises(ContractError, match="Cell Ontology identifier"):
        risk_coverage_sweep(without_identifiers, agreement="relaxed")
    with pytest.raises(ContractError, match="local ontology snapshot"):
        risk_coverage_sweep(frame, agreement="relaxed")
    with pytest.raises(ContractError, match="risk-coverage sweep requires"):
        risk_coverage_sweep(without_identifiers.drop(columns=["confidence"]))


def test_an_unrecognised_agreement_is_refused_and_a_snapshot_changes_nothing() -> None:
    """A misspelled mode must not fall back to the strict curve.

    `agreement` is the whole content of the request, so a typo that quietly
    produced the strict curve would hand the caller the numbers this metric was
    added to stop them from having, with nothing to say they were the wrong ones.
    Supplying a snapshot alongside the default strict agreement is a caller
    preparing for the relaxed variant, and it must not change what the strict
    curve measures.
    """

    frame = _frame(
        ["central memory CD4-positive, alpha-beta T cell"],
        ["T cell"],
        [CENTRAL_MEMORY_CD4],
        [T_CELL],
        [0.95],
    )
    with pytest.raises(ContractError, match="'strict', 'descendant', or 'relaxed'"):
        risk_coverage_sweep(frame, agreement="granularity", ontology=_ontology())
    strict = risk_coverage_sweep(frame)
    with_snapshot = risk_coverage_sweep(frame, ontology=_ontology())
    assert [point.model_dump(mode="json") for point in strict] == [
        point.model_dump(mode="json") for point in with_snapshot
    ]
    assert strict[0].selective_accuracy == 0.0


def test_a_published_result_carries_both_numbers(tmp_path: Path) -> None:
    """The reader of the artifact has to be able to see the two of them.

    A relaxed accuracy that exists only on the returned model would leave every
    published comparison reporting the strict number alone, which is the state
    this change exists to end. The ontology table on disk is the same code path
    `benchmark evaluate --ontology` takes, so the published provenance and the
    relaxed field come from one call.
    """

    ontology_path = tmp_path / "mini_cl.tsv"
    ontology_path.write_text(
        "id\tname\tparents\n"
        f"{NATIVE_CELL}\tnative cell\t{CELL}\n"
        f"{LEUKOCYTE}\tleukocyte\t{NATIVE_CELL}\n"
        f"{LYMPHOCYTE}\tlymphocyte\t{LEUKOCYTE}\n"
        f"{T_CELL}\tT cell\t{LYMPHOCYTE}\n"
        f"{CENTRAL_MEMORY_CD4}\tcentral memory CD4-positive, alpha-beta T cell\t{T_CELL}\n"
    )
    predictions = tmp_path / "predictions.tsv"
    predictions.write_text(
        "barcode\ttrue_label\tpredicted_label\tconfidence\ttrue_cl_id\tpredicted_cl_id\n"
        f"b0\tcentral memory CD4-positive, alpha-beta T cell\tT cell\t0.9"
        f"\t{CENTRAL_MEMORY_CD4}\t{T_CELL}\n"
        f"b1\tT cell\tT cell\t0.8\t{T_CELL}\t{T_CELL}\n"
    )
    result = evaluate_file(
        predictions,
        system="fixture",
        dataset="relaxed-publish",
        output_path=tmp_path / "result.json",
        ontology=ontology_path,
    )
    assert result.accuracy == 0.5
    assert result.relaxed_accuracy == 1.0
    payload = json.loads((tmp_path / "result.json").read_text())
    assert payload["accuracy"] == 0.5
    assert payload["relaxed_accuracy"] == 1.0
    assert payload["ontology_distance_error"] is not None


def test_the_public_sweep_reads_a_relaxed_curve_from_disk(tmp_path: Path) -> None:
    """A comparative claim is checked from outside, so the curve is reachable.

    Whoever holds the competing export is the only person who can compare two
    systems at matched coverage, and with a ground truth annotated more finely
    than either system's answer they can only do it against the relaxed
    criterion. The public entry point therefore takes the same snapshot argument
    the scalar evaluation does, as a path or as a loaded ontology.
    """

    predictions = tmp_path / "predictions.tsv"
    predictions.write_text(
        "barcode\ttrue_label\tpredicted_label\tconfidence"
        f"\ttrue_cl_id\tpredicted_cl_id\n"
        f"b0\tcentral memory CD4-positive, alpha-beta T cell\tT cell\t0.99"
        f"\t{CENTRAL_MEMORY_CD4}\t{T_CELL}\n"
        f"b1\tT cell\tmonocyte\t0.8\t{T_CELL}\t{MONOCYTE}\n"
    )
    frame = pd.read_csv(predictions, sep="\t", keep_default_na=False)
    expected = risk_coverage_sweep(
        frame, thresholds=[0.0, 0.95], agreement="relaxed", ontology=_ontology()
    )
    from_path = sweep_benchmark_risk_coverage(
        predictions, thresholds=[0.0, 0.95], agreement="relaxed", ontology=_ontology()
    )
    assert [point.model_dump(mode="json") for point in from_path] == [
        point.model_dump(mode="json") for point in expected
    ]
    assert [point.agreement for point in from_path] == ["relaxed", "relaxed"]
    assert [point.selective_accuracy for point in from_path] == [0.5, 1.0]
    assert [
        point.selective_accuracy
        for point in sweep_benchmark_risk_coverage(
            predictions, thresholds=[0.0, 0.95]
        )
    ] == [0.0, 0.0]


def _directional_ontology() -> CellOntology:
    """A chain deep enough to have a child, a parent, and a grandparent of a term.

    Built here rather than reused from the module fixture, because that one has
    no term with a child and these tests are about which side of the truth a
    prediction falls on.
    """

    return CellOntology(
        {
            CELL: OntologyTerm(CELL, "cell"),
            NATIVE_CELL: OntologyTerm(NATIVE_CELL, "native cell", (CELL,)),
            LEUKOCYTE: OntologyTerm(LEUKOCYTE, "leukocyte", (NATIVE_CELL,)),
            LYMPHOCYTE: OntologyTerm(LYMPHOCYTE, "lymphocyte", (LEUKOCYTE,)),
            T_CELL: OntologyTerm(T_CELL, "T cell", (LYMPHOCYTE,)),
            CD4_POSITIVE_T: OntologyTerm(CD4_POSITIVE_T, "CD4-positive T cell", (T_CELL,)),
            CENTRAL_MEMORY_CD4: OntologyTerm(
                CENTRAL_MEMORY_CD4, "central memory CD4-positive T cell", (CD4_POSITIVE_T,)
            ),
            B_CELL: OntologyTerm(B_CELL, "B cell", (LYMPHOCYTE,)),
        },
        source="directional",
    )


def _directional_frame(truth: str, predicted: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "barcode": ["a"],
            "true_label": ["truth"],
            "predicted_label": ["prediction"],
            "confidence": [0.9],
            "true_cl_id": [truth],
            "predicted_cl_id": [predicted],
        }
    )


def test_the_two_relaxations_differ_in_exactly_one_direction() -> None:
    """A finer answer is credited by both; a coarser one only by the symmetric rule.

    The one-directional rule is the one the largest published cross-tissue
    benchmark uses, and it states its reason: a valid subtype remains consistent
    with the original label, whereas a coarser label asserts less than the truth
    and is counted incorrect. That asymmetry is what stops a system buying
    apparent accuracy by answering vaguely.

    Both are reported beside the strict number, because publishing only the
    flattering direction would mean choosing it after seeing the result. The gap
    between them is the evidence about how much of an accuracy is identification
    and how much is specificity.
    """

    ontology = _directional_ontology()
    assert ontology.is_descendant(CENTRAL_MEMORY_CD4, CD4_POSITIVE_T)
    assert ontology.is_descendant(CD4_POSITIVE_T, T_CELL)

    for predicted, descendant, relaxed in (
        (CD4_POSITIVE_T, 1.0, 1.0),
        (CENTRAL_MEMORY_CD4, 1.0, 1.0),
        (T_CELL, 0.0, 1.0),
        (LYMPHOCYTE, 0.0, 1.0),
        (B_CELL, 0.0, 0.0),
    ):
        measured = evaluate_predictions(
            _directional_frame(CD4_POSITIVE_T, predicted),
            system="fixture",
            dataset="cl-direction",
            ontology=ontology,
        )
        assert measured.descendant_accuracy == pytest.approx(descendant), predicted
        assert measured.relaxed_accuracy == pytest.approx(relaxed), predicted

    coarse = evaluate_predictions(
        _directional_frame(CD4_POSITIVE_T, T_CELL),
        system="fixture",
        dataset="cl-direction",
        ontology=ontology,
    )
    assert coarse.accuracy == 0.0
    assert coarse.descendant_accuracy == 0.0
    assert coarse.relaxed_accuracy == 1.0
    assert coarse.ontology_distance_error is not None


def test_the_sweep_under_each_relaxation_is_a_different_measurement() -> None:
    """Three curves over one table are three measurements and each records its own.

    A curve that credits a coarser call answers a different question from one
    that does not, so a point has to say which question it answered, or two curves
    could be compared as though they measured the same thing.
    """

    ontology = _directional_ontology()
    frame = pd.DataFrame(
        {
            "barcode": ["a", "b"],
            "true_label": ["CD4 T", "CD4 T"],
            "predicted_label": ["T cell", "central memory CD4 T"],
            "confidence": [0.9, 0.9],
            "true_cl_id": [CD4_POSITIVE_T, CD4_POSITIVE_T],
            "predicted_cl_id": [T_CELL, CENTRAL_MEMORY_CD4],
        }
    )
    curves = {
        mode: risk_coverage_sweep(frame, agreement=mode, ontology=ontology)
        for mode in ("strict", "descendant", "relaxed")
    }
    assert {points[0].coverage for points in curves.values()} == {1.0}
    at_zero = {mode: points[0].selective_accuracy for mode, points in curves.items()}
    assert at_zero["strict"] == pytest.approx(0.0)
    assert at_zero["descendant"] == pytest.approx(0.5)
    assert at_zero["relaxed"] == pytest.approx(1.0)
    for mode, points in curves.items():
        assert {point.agreement for point in points} == {mode}

    with pytest.raises(ContractError, match="descendant"):
        risk_coverage_sweep(frame, agreement="strictest", ontology=ontology)

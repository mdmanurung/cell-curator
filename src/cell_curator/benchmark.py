"""Neutral, dependency-free evaluation of exported annotation predictions."""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pydantic import ValidationError

from .contracts import ContractError, sha256
from .models import (
    Agreement,
    BenchmarkDatasetManifest,
    BenchmarkResult,
    RiskCoveragePoint,
)
from .ontology import CellOntology
from .provenance import publish_json

REQUIRED_COLUMNS = {"barcode", "true_label", "predicted_label"}

ONTOLOGY_COLUMN_PAIRS: tuple[tuple[str, str], ...] = (
    ("true_cl_id", "predicted_cl_id"),
    ("true_cell_type_ontology_term_id", "predicted_cell_type_ontology_term_id"),
)

DEFAULT_RISK_COVERAGE_THRESHOLDS: tuple[float, ...] = tuple(
    round(step / 20, 2) for step in range(21)
)


def _class_metrics(truth: pd.Series, prediction: pd.Series) -> tuple[float, float, float]:
    labels = sorted(set(truth.astype(str)) | set(prediction.astype(str)))
    f1s = []
    total_tp = total_fp = total_fn = 0
    for label in labels:
        true = truth.astype(str).eq(label)
        predicted = prediction.astype(str).eq(label)
        tp = int((true & predicted).sum())
        fp = int((~true & predicted).sum())
        fn = int((true & ~predicted).sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1s.append(2 * precision * recall / (precision + recall) if precision + recall else 0.0)
        total_tp += tp
        total_fp += fp
        total_fn += fn
    micro_precision = total_tp / (total_tp + total_fp) if total_tp + total_fp else 0.0
    micro_recall = total_tp / (total_tp + total_fn) if total_tp + total_fn else 0.0
    micro_f1 = (
        2 * micro_precision * micro_recall / (micro_precision + micro_recall)
        if micro_precision + micro_recall
        else 0.0
    )
    accuracy = float(truth.astype(str).eq(prediction.astype(str)).mean())
    return accuracy, float(np.mean(f1s)) if f1s else 0.0, micro_f1


def _unknown_metrics(
    truth: pd.Series, prediction: pd.Series
) -> tuple[float | None, float | None]:
    true_unknown = truth.astype(str).eq("Unknown")
    predicted_unknown = prediction.astype(str).eq("Unknown")
    if not true_unknown.any() and not predicted_unknown.any():
        return None, None
    tp = int((true_unknown & predicted_unknown).sum())
    precision = tp / int(predicted_unknown.sum()) if predicted_unknown.any() else 0.0
    recall = tp / int(true_unknown.sum()) if true_unknown.any() else 0.0
    return precision, recall


def _confidence_values(frame: pd.DataFrame) -> np.ndarray:
    """Read the per-cell confidence column, refusing anything unusable.

    Every claim the instrument makes about calibration or selective prediction
    rests on this one column, and a value outside ``[0, 1]`` or a NaN would let
    a table assert a confidence it cannot have. Reading it in one place keeps
    the calibration figures and the risk-coverage curve from drifting into two
    different notions of confidence.
    """

    confidence = pd.to_numeric(frame["confidence"], errors="coerce").to_numpy()
    if not np.isfinite(confidence).all() or ((confidence < 0) | (confidence > 1)).any():
        raise ContractError("benchmark confidence must be finite and within [0, 1]")
    return confidence


def _calibration(frame: pd.DataFrame, bins: int = 10) -> tuple[float | None, float | None]:
    if "confidence" not in frame:
        return None, None
    confidence = _confidence_values(frame)
    correct = (
        frame["true_label"]
        .astype(str)
        .eq(frame["predicted_label"].astype(str))
        .astype(float)
        .to_numpy()
    )
    brier = float(np.mean((confidence - correct) ** 2))
    edges = np.linspace(0, 1, bins + 1)
    ece = 0.0
    for left, right in zip(edges[:-1], edges[1:], strict=True):
        selected = (confidence >= left) & (
            confidence <= right if right == 1 else confidence < right
        )
        if selected.any():
            ece += selected.mean() * abs(
                float(confidence[selected].mean() - correct[selected].mean())
            )
    return brier, float(ece)


def _hierarchy_error(frame: pd.DataFrame) -> float | None:
    if not {"true_path", "predicted_path"}.issubset(frame.columns):
        return None
    distances = []
    for true_path, predicted_path in zip(
        frame["true_path"], frame["predicted_path"], strict=True
    ):
        truth = [item for item in str(true_path).split("/") if item]
        prediction = [item for item in str(predicted_path).split("/") if item]
        # A blank side is unmeasured: no truth path, or an abstention. Scoring
        # blank against blank as an exact match made the metric equal coverage.
        # Abstentions are left out, as in the ontology metrics.
        if not truth or not prediction:
            continue
        common = 0
        for left, right in zip(truth, prediction, strict=False):
            if left != right:
                break
            common += 1
        distance = (len(truth) - common) + (len(prediction) - common)
        maximum = max(1, len(truth) + len(prediction))
        distances.append(distance / maximum)
    return float(np.mean(distances)) if distances else None


def _term_distance(ontology: CellOntology, truth: str, prediction: str) -> float:
    """Distance between two Cell Ontology terms over the supplied `is_a` graph.

    For a true term ``T`` and a predicted term ``P`` write ``A(x)`` for ``x``
    together with all of its ancestors, ``d`` for :meth:`CellOntology.depth`, and
    let ``S`` be the *deepest* term of ``A(T) ∩ A(P)``, breaking ties by
    accession so that one table always yields one number. With ``a = d(T) -
    d(S)`` and ``b = d(P) - d(S)`` the distance is::

        (a + b) / (a + b + d(S))

    The numerator is the number of `is_a` edges on which the two terms
    disagree: the ``a`` edges the truth climbs from ``T`` to the agreement
    point, plus the ``b`` edges the prediction invents below it. The denominator
    adds the ``d(S)`` edges they genuinely share above the agreement point, so
    agreeing high in the hierarchy is worth more than agreeing between two
    leaves. Three properties follow from that shape and are the reason it is
    used rather than a bare path edit: an exact match is 0; a prediction that
    descends from the truth is charged only for its extra depth; and a
    prediction that shares nothing but the roots of the supplied graph is
    charged the whole path, which is ``d(S) = 0`` and therefore exactly 1.0.

    Three limits are deliberate and a caller reading a result needs all three.
    The measure scores *how far* a prediction sits from the truth rather than
    *which branch* it took, so a wrong sibling is distinguished from a right one
    only by its extra depth; the label-level accuracy this instrument already
    publishes is what separates the branches. A pair with no shared ancestry
    below the graph roots cannot be placed on a path, so it takes the
    specification's 1.0 rather than a fabricated intermediate value.

    The third is a property of the ontology rather than the formula, and it
    changes how a number should be read. On the complete, connected Cell
    Ontology every pair of real terms shares `cell` and `native cell`, and
    those are not roots of the graph, so the agreement point always sits at
    least one edge below a root and the score never reaches 1.0. A B cell
    predicted as a fibroblast measures roughly 0.8, not 1.0, because the two
    really are both native cells. Only a subtree export, a broken snapshot, or a
    term in another component produces the 1.0 above. A published value should
    therefore be read as an ordering among systems on the same snapshot, not as
    a fraction of errors made, and the label-level accuracy is what answers the
    latter question.
    """

    if truth == prediction:
        return 0.0
    shared = (
        (ontology.ancestors(truth) | {truth})
        & (ontology.ancestors(prediction) | {prediction})
    ) - ontology.roots()
    shared = {accession for accession in shared if accession in ontology.terms}
    if not shared:
        return 1.0
    agreement = max(shared, key=lambda accession: (ontology.depth(accession), accession))
    depth_agreement = ontology.depth(agreement)
    apart = (ontology.depth(truth) - depth_agreement) + (
        ontology.depth(prediction) - depth_agreement
    )
    return apart / (apart + depth_agreement)


def _ontology_identifier_claim(
    frame: pd.DataFrame, ontology: CellOntology | None, *, metric: str
) -> tuple[tuple[str, str], CellOntology] | None:
    """The CL columns a table offers and the snapshot they were checked against.

    A table may present its terms under either spelling in ``ONTOLOGY_COLUMN_PAIRS``
    and the first pair present wins, so a table carrying both is scored under the
    package's own naming and the Census spelling is read only when it is the one
    there. Presenting neither means the table makes no ontology claim at all, which
    is reported as ``None`` by the caller rather than as a measurement. Presenting
    identifiers with no snapshot behind them raises here, because an ontology claim
    scored against nothing is a number with no evidence in it, and the message names
    both the metric and the columns at fault so a caller reading it can tell which
    of the two ontology metrics it was asking for.
    """

    pair = next(
        (pair for pair in ONTOLOGY_COLUMN_PAIRS if set(pair).issubset(frame.columns)), None
    )
    if pair is None:
        return None
    if ontology is None:
        raise ContractError(
            "benchmark ontology identifiers require a local ontology snapshot; "
            f"supply one to measure {metric} from {list(pair)}"
        )
    return pair, ontology


def _resolved_identifier_pairs(
    frame: pd.DataFrame, pair: tuple[str, str], ontology: CellOntology
) -> list[tuple[int, str, str]]:
    """Resolve every row's identifier pair, keeping the rows that carry one.

    Both identifiers go through the same :meth:`CellOntology.resolve` the distance
    metric uses, so every ontology metric in this instrument accepts and refuses
    exactly the same tables: a term name is resolved to its accession when the
    snapshot is unambiguous, and an identifier the snapshot does not know raises
    rather than being scored as maximally far away. Rows where either side is blank
    are left out, because an abstention has no term to be near or far from. The row
    index travels with each pair so a caller that needs a per-row vector can say
    which rows it was unable to measure.
    """

    resolved: list[tuple[int, str, str]] = []
    for row, (truth, prediction) in enumerate(
        zip(frame[pair[0]].astype(str), frame[pair[1]].astype(str), strict=True)
    ):
        true_id = truth.strip()
        predicted_id = prediction.strip()
        if not true_id or not predicted_id:
            continue
        resolved.append(
            (
                row,
                ontology.resolve(true_id).accession,
                ontology.resolve(predicted_id).accession,
            )
        )
    return resolved


def _ontology_distance_error(
    frame: pd.DataFrame, ontology: CellOntology | None
) -> float | None:
    """Measure CL distance, or report nothing when the table carries no CL claim.

    A table may present a true and a predicted Cell Ontology identifier under
    either the package's ``true_cl_id``/``predicted_cl_id`` naming or the Census
    ``*_cell_type_ontology_term_id`` spelling; the first pair present wins. The
    three outcomes are kept distinct on purpose, because a measurement that
    cannot be taken is not the same as a measurement of zero:

    * no pair of identifier columns, so the table makes no ontology claim, and
      the metric is ``None`` -- never a copy of the slash-path hierarchy error,
      which is a different quantity over different evidence;
    * identifier columns present but no ontology supplied, which raises rather
      than reporting ``None``: the table claims to be scored against the Cell
      Ontology and nothing in the result would say that the claim went
      unchecked;
    * an accession absent from the supplied snapshot, which raises for the same
      reason a label inconsistent with its CL id raises elsewhere.

    Rows where either identifier is blank carry no pair, because an abstention
    has no term to be near or far from and is already accounted for by the
    Unknown precision and recall this instrument reports. They are excluded
    rather than charged full error, and the reported value is the mean over the
    rows that do carry a pair; a table whose identifiers are all blank is
    therefore unmeasured and reports ``None``.
    """

    claim = _ontology_identifier_claim(frame, ontology, metric="ontology_distance_error")
    if claim is None:
        return None
    pair, snapshot = claim
    distances = [
        _term_distance(snapshot, truth, prediction)
        for _, truth, prediction in _resolved_identifier_pairs(frame, pair, snapshot)
    ]
    return float(np.mean(distances)) if distances else None


def _relaxed_agreement(
    ontology: CellOntology, truth: str, prediction: str, *, mode: str = "relaxed"
) -> bool:
    """Whether two resolved Cell Ontology terms agree, at one of two relaxations.

    :meth:`CellOntology.is_descendant` read over the same ``ancestors`` walk the
    distance metric uses and never a second copy of the graph. Both terms arrive
    already resolved, so an identifier the snapshot does not know has raised
    before this is reached and a term it considers obsolete is treated exactly as
    the distance metric treats it.

    The two relaxations differ in which direction they forgive, and the difference
    is not cosmetic:

    ``descendant`` credits a prediction that is the truth or a more specific
    descendant of it, and does *not* credit a coarser answer. This is the
    one-directional rule used by the largest published cross-tissue benchmark,
    which states the reason: a valid subtype remains consistent with the original
    label, whereas a coarser label asserts less than the truth and is counted
    incorrect. It is the rule that cannot be satisfied by answering vaguely.

    ``relaxed`` credits both directions, because a table whose vocabulary is
    genuinely coarser than its ground truth is measuring a real difference of
    design rather than a mistake. It is a ceiling on what the system was asked,
    and it is why a run that abstains from specificity scores well here.

    Publishing only one of them would be choosing the flattering direction, so
    both are reported beside the strict number and the gap between them is the
    evidence in a result about how much of an accuracy is identification and how
    much is specificity.
    """

    if truth == prediction:
        return True
    # `truth in ancestors(prediction)` means the prediction is a descendant of the
    # truth, that is, a more specific answer. `prediction in ancestors(truth)` is
    # the opposite: the prediction is a coarser answer, asserting less than the
    # truth. The one-directional rule credits the first and not the second.
    if mode in {"descendant", "relaxed"} and truth in ontology.ancestors(prediction):
        return True
    return mode == "relaxed" and prediction in ontology.ancestors(truth)


def _granularity_accuracy(
    frame: pd.DataFrame, ontology: CellOntology | None, *, mode: str
) -> float | None:
    """Share of calls that survive a difference of Cell Ontology granularity.

    A prediction agrees when its term is the true term, or under ``relaxed`` when
    one is an ``is_a`` ancestor of the other, so a call of "T cell" scored
    against a truth of "central memory CD4-positive, alpha-beta T cell" is
    credited where ``accuracy`` reads it as a complete error. That gap is what
    this metric exists to close. The benchmark ground truth is annotated at a
    finer Cell Ontology granularity than any system is asked to reproduce, so
    exact string equality measures how nearly a system imitates the annotator's
    chosen depth rather than whether it found the cell, and a system that is
    right in every cell publishes a number near zero. A selective-prediction
    curve built on the same rule inherits the defect: its selective accuracy
    could never approach one, so two systems could not be
    compared at matched coverage on it at all.

    The relaxation is deliberately symmetric, and that is the one caveat a reader
    has to carry into the number. Being too specific and being too general are the
    same kind of mistake here, so a system that answers "leukocyte" for every cell
    is credited on every row, and the number is a ceiling on what the system was
    asked rather than a score for what it did. It is a relaxation and not a
    replacement for that reason: ``accuracy``, ``macro_f1``, ``micro_f1``,
    ``coverage``, ``selective_accuracy``, and ``coverage_risk`` are computed exactly
    as they were, and the result carries both numbers side by side so a reader can
    always see the strict one. Collapsing every prediction onto a shared ancestor
    would let a system buy apparent accuracy by declining to be specific, which is
    the failure the strict numbers exist to expose, and it is why this is a second
    field rather than a corrected ``accuracy``.

    Two further limits belong to the evidence rather than the definition. The
    measurement exists only where a table supplies CL identifiers, and it says
    ``None`` rather than borrowing the strict number when it does not; the
    three-way refusal, including a table that claims identifiers with no snapshot
    to check them against, is that of ``ontology_distance_error`` and is documented
    there. And where the strict metric is a fraction of all rows, this is a fraction
    of the rows that carry a pair, because an abstention is scored through the
    Unknown precision and recall rather than charged as a term that is maximally
    far away; a table whose identifiers are all blank is therefore unmeasured
    rather than perfectly wrong.
    """

    claim = _ontology_identifier_claim(frame, ontology, metric="relaxed_accuracy")
    if claim is None:
        return None
    pair, snapshot = claim
    agreements = [
        _relaxed_agreement(snapshot, truth, prediction, mode=mode)
        for _, truth, prediction in _resolved_identifier_pairs(frame, pair, snapshot)
    ]
    return float(np.mean(agreements)) if agreements else None


def _relaxed_correctness(
    frame: pd.DataFrame, ontology: CellOntology | None, *, mode: str = "relaxed"
) -> np.ndarray:
    """Per-cell relaxed correctness, for the selective-prediction curve.

    A cell is relaxed-correct on exactly the condition ``relaxed_accuracy`` counts
    rows, so a full-coverage point of the relaxed curve reproduces that field on a
    table whose identifiers are all present. One rule differs, and the denominator
    is why: a row whose identifiers are blank is dropped from ``relaxed_accuracy``
    and charged as incorrect here. A curve's denominator is the number of cells
    called at a threshold, so a row that cannot be shown to agree has to lower the
    selective accuracy of every threshold that calls it; removing the row instead
    would publish a coverage the table does not have and an accuracy no system
    produced. A table with no identifiers at all has no relaxed curve to draw, and
    that is a refusal rather than an empty curve, because a ``None`` selective
    accuracy on a returned point already means something else -- a threshold that
    called nothing -- and reusing it here would make an unmeasured curve
    indistinguishable from an abstaining system.
    """

    claim = _ontology_identifier_claim(frame, ontology, metric="relaxed risk-coverage agreement")
    if claim is None:
        raise ContractError(
            "benchmark relaxed risk-coverage agreement requires a true and a predicted "
            f"Cell Ontology identifier on each row; add one of {list(ONTOLOGY_COLUMN_PAIRS)} "
            "to the table or ask for the strict sweep"
        )
    pair, snapshot = claim
    correct = np.zeros(len(frame), dtype=bool)
    for row, truth, prediction in _resolved_identifier_pairs(frame, pair, snapshot):
        correct[row] = _relaxed_agreement(snapshot, truth, prediction, mode=mode)
    return correct


def _binary_f1(frame: pd.DataFrame, truth_column: str, prediction_column: str) -> float | None:
    if not {truth_column, prediction_column}.issubset(frame.columns):
        return None

    def boolean(value: Any) -> bool:
        return str(value).strip().lower() in {"1", "true", "yes", "mixed", "doublet"}

    truth = frame[truth_column].map(boolean)
    predicted = frame[prediction_column].map(boolean)
    tp = int((truth & predicted).sum())
    fp = int((~truth & predicted).sum())
    fn = int((truth & ~predicted).sum())
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0


def _optional_agreement(frame: pd.DataFrame, left: str, right: str) -> float | None:
    if not {left, right}.issubset(frame.columns):
        return None
    return float(frame[left].astype(str).eq(frame[right].astype(str)).mean())


def _optional_boolean_fraction(frame: pd.DataFrame, column: str) -> float | None:
    if column not in frame:
        return None
    values = frame[column].astype(str).str.lower().isin({"1", "true", "yes", "pass"})
    return float(values.mean())


def _optional_numeric(
    frame: pd.DataFrame, column: str, *, aggregation: str = "max"
) -> float | None:
    if column not in frame:
        return None
    values = pd.to_numeric(frame[column], errors="coerce")
    if values.isna().all():
        return None
    return float(values.max() if aggregation == "max" else values.mean())


def risk_coverage_sweep(
    frame: pd.DataFrame,
    *,
    thresholds: Sequence[float] | None = None,
    agreement: Agreement = "strict",
    ontology: CellOntology | None = None,
) -> tuple[RiskCoveragePoint, ...]:
    """Report coverage, selective accuracy, and risk across a confidence grid.

    The published claim for this package is a selective-prediction claim:
    comparable accuracy at matched coverage, with abstention doing real work.
    That claim is a statement about a curve, not about the single operating
    point ``evaluate_predictions`` reports, and a curve is the only form in
    which two systems can be compared at equal coverage. This function is that
    curve.

    It is a pure function over the same neutral prediction frame and reads only
    ``true_label``, ``predicted_label``, and ``confidence``, so it scores any
    system's export, not cell-curator's own, and it cannot be steered by a
    column it does not understand. A cell is *called* at threshold ``t`` when
    ``confidence >= t``, which makes the grid endpoints exact and inclusive: at
    ``0.0`` every cell is called and the point reproduces the table's plain
    accuracy, and at ``1.0`` only cells claiming maximal confidence are called.
    Exact-match correctness is reused rather than reimplemented so a point can
    never disagree with ``accuracy`` at full coverage.

    Note that this is a different operating-point axis from the ``coverage``
    published in ``BenchmarkResult``, which counts cells that were not labelled
    ``Unknown``. Both are legitimate; a system whose label-based abstention
    works and whose confidence ranking does not will show one curve and not the
    other.

    A threshold at which nothing is called reports ``coverage`` 0.0 with a
    ``None`` accuracy and a ``None`` risk, because an empty set of called cells
    has no measured risk; reporting 0.0 accuracy or 1.0 risk there would
    fabricate the worst and best possible system at the same point. A frame
    without a confidence column raises rather than degrading to a reason
    string: the sweep's entire content is a function of that column, and an
    empty or explanatory return would be indistinguishable at the call site from
    a table that abstains everywhere. The same guard and the same message as
    the calibration figures apply to values outside ``[0, 1]``.

    The default grid is twenty-one thresholds from 0.0 to 1.0 in steps of 0.05.
    A supplied grid must be finite, within ``[0, 1]``, non-empty, and strictly
    increasing, so that the endpoints of a reported curve are unambiguous and
    two runs over the same table cannot disagree about the order of the points.

    ``agreement="relaxed"`` counts a cell correct when its predicted Cell
    Ontology term equals the true term or is an ``is_a`` ancestor or descendant
    of it, which is the relaxation ``relaxed_accuracy`` reports and is
    documented there, asymmetry included. It is opt-in because the strict curve
    is the one the published numbers are, and it changes what a curve is a claim
    about: a curve that credits a coarser or finer call than the truth is asking
    whether the system found the cell, not whether it named it the way the
    annotator did. Everything else is unchanged, so the two variants are
    comparable to each other -- same grid, same inclusive ``confidence >=
    threshold`` boundary, same refusal to report a point where nothing is
    called, same ``None`` rather than a fabricated accuracy. A relaxed curve
    requires a true and a predicted CL identifier on each row and a local
    ontology snapshot to check them against, and raises without them for the
    reason ``_relaxed_correctness`` gives. The two variants differ in exactly
    what they read: the strict curve counts equality of ``true_label`` and
    ``predicted_label`` and the relaxed curve counts agreement of the two CL
    identifiers, so each demands the evidence it scores and neither is run on the
    other's, and the strict curve ignores ``ontology`` entirely, so supplying a
    snapshot cannot silently change what an existing call measures. Every point
    records the agreement it was computed under, because two curves over one
    table are two different measurements and neither is the other.
    """

    if agreement not in {"strict", "descendant", "relaxed"}:
        raise ContractError(
            "risk-coverage agreement must be 'strict', 'descendant', or "
            f"'relaxed', not {agreement!r}"
        )
    if "confidence" not in frame:
        raise ContractError(
            "benchmark risk-coverage sweep requires a per-cell confidence column"
        )
    if agreement == "strict":
        missing = {"true_label", "predicted_label"}.difference(frame.columns)
        if missing:
            raise ContractError(
                f"neutral prediction table misses columns: {sorted(missing)}"
            )
    confidence = _confidence_values(frame)
    grid = _validated_thresholds(thresholds)
    if agreement in {"descendant", "relaxed"}:
        correct = _relaxed_correctness(frame, ontology, mode=agreement)
    else:
        correct = (
            frame["true_label"]
            .astype(str)
            .eq(frame["predicted_label"].astype(str))
            .to_numpy()
        )
    points: list[RiskCoveragePoint] = []
    for threshold in grid:
        called = confidence >= threshold
        n_called = int(called.sum())
        n_correct = int(correct[called].sum())
        accuracy = float(correct[called].mean()) if n_called else None
        points.append(
            RiskCoveragePoint(
                threshold=threshold,
                n_total=len(frame),
                n_called=n_called,
                n_correct=n_correct,
                coverage=n_called / len(frame) if len(frame) else 0.0,
                selective_accuracy=accuracy,
                risk=None if accuracy is None else 1.0 - accuracy,
                agreement=agreement,
            )
        )
    return tuple(points)


def _validated_thresholds(thresholds: Sequence[float] | None) -> tuple[float, ...]:
    """Accept a caller-supplied grid only if its endpoints mean what they say."""

    if thresholds is None:
        return DEFAULT_RISK_COVERAGE_THRESHOLDS
    grid = tuple(float(value) for value in thresholds)
    if not grid:
        raise ContractError("risk-coverage thresholds must not be empty")
    for value in grid:
        if not np.isfinite(value) or value < 0.0 or value > 1.0:
            raise ContractError("risk-coverage thresholds must be finite and within [0, 1]")
    if any(left >= right for left, right in zip(grid[:-1], grid[1:], strict=True)):
        raise ContractError("risk-coverage thresholds must be strictly increasing")
    return grid


def evaluate_predictions(
    frame: pd.DataFrame,
    *,
    system: str,
    dataset: str,
    ontology: CellOntology | None = None,
) -> BenchmarkResult:
    """Score a neutral prediction table at the operating point it was run at.

    The optional ``ontology`` is a local Cell Ontology snapshot, supplied by the
    caller rather than discovered here, for three reasons: the instrument is
    dependency-free and neutral, so it must not fetch or pin an ontology on its
    own; the snapshot is part of the claim, so whoever publishes the result has
    to name the file it measured against; and a plain evaluation of a table that
    carries no CL identifiers must not be made to depend on an ontology download
    that changes nothing about the numbers. When the table does carry CL
    identifiers, the snapshot becomes mandatory, because an unchecked ontology
    claim is exactly the unmeasured claim this instrument exists to refuse.

    ``ontology_distance_error`` is a measurement over the ``is_a`` graph and is
    reported as ``None`` when the table makes no ontology claim. It is not, and
    must not be, a restatement of ``hierarchy_error``: the latter is a distance
    over slash-separated label paths supplied by the exporting system, so
    copying one into the other would publish a number under an ontology field
    name while no ontology was consulted.

    ``relaxed_accuracy`` is the granularity-relaxed counterpart of ``accuracy``
    and is reported beside it rather than in place of it. A benchmark table whose
    truth is annotated more finely than any system is asked to answer would
    otherwise read every correct call as a complete error, and the strict
    accuracy, the F1 pair, ``coverage``, and ``selective_accuracy`` are left
    exactly as strict as they were for that reason: a reader has to be able to
    see both the exact-label number and the relaxed one, because the gap between
    them is the only evidence in the result about how much of a system's
    accuracy is specificity rather than identification. The two are not
    interchangeable and neither is derived from the other.
    """

    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise ContractError(f"neutral prediction table misses columns: {sorted(missing)}")
    if frame["barcode"].astype(str).duplicated().any():
        raise ContractError("neutral prediction table must contain one row per barcode")
    accuracy, macro_f1, micro_f1 = _class_metrics(frame["true_label"], frame["predicted_label"])
    called = ~frame["predicted_label"].astype(str).eq("Unknown")
    coverage = float(called.mean())
    selective_accuracy = (
        float(
            frame.loc[called, "true_label"]
            .astype(str)
            .eq(frame.loc[called, "predicted_label"].astype(str))
            .mean()
        )
        if called.any()
        else 0.0
    )
    unknown_precision, unknown_recall = _unknown_metrics(
        frame["true_label"], frame["predicted_label"]
    )
    brier, ece = _calibration(frame)
    hierarchy_error = _hierarchy_error(frame)
    unknown_f1 = (
        2 * unknown_precision * unknown_recall / (unknown_precision + unknown_recall)
        if unknown_precision is not None
        and unknown_recall is not None
        and unknown_precision + unknown_recall
        else None
    )
    return BenchmarkResult(
        system=system,
        dataset=dataset,
        n=len(frame),
        accuracy=accuracy,
        macro_f1=macro_f1,
        micro_f1=micro_f1,
        coverage=coverage,
        selective_accuracy=selective_accuracy,
        coverage_risk=1.0 - selective_accuracy,
        hierarchical_accuracy=(1.0 - hierarchy_error) if hierarchy_error is not None else None,
        descendant_accuracy=_granularity_accuracy(frame, ontology, mode="descendant"),
        relaxed_accuracy=_granularity_accuracy(frame, ontology, mode="relaxed"),
        unknown_precision=unknown_precision,
        unknown_recall=unknown_recall,
        unknown_f1=unknown_f1,
        brier_score=brier,
        expected_calibration_error=ece,
        hierarchy_error=hierarchy_error,
        ontology_distance_error=_ontology_distance_error(frame, ontology),
        mixed_detection_f1=_binary_f1(frame, "true_mixed", "predicted_mixed"),
        split_decision_accuracy=_optional_agreement(
            frame, "true_split_decision", "predicted_split_decision"
        ),
        repeat_reproducibility=_optional_agreement(
            frame, "predicted_label", "repeat_predicted_label"
        ),
        input_invariance=_optional_boolean_fraction(frame, "input_invariant"),
        unsplit_parent_invariance=_optional_boolean_fraction(frame, "unsplit_parent_invariant"),
        runtime_seconds=_optional_numeric(frame, "runtime_seconds"),
        peak_memory_mb=_optional_numeric(frame, "peak_memory_mb"),
        artifact_completeness=_optional_boolean_fraction(frame, "artifact_complete"),
    )


def read_predictions(prediction_path: str | Path) -> pd.DataFrame:
    """Read a neutral prediction table, choosing the separator from the suffix.

    A comma-separated table named ``.tsv`` otherwise arrives as one column and
    fails later as a missing required column, which reports a formatting slip as
    a contract violation.
    """

    path = Path(prediction_path)
    separator = "\t" if path.suffix.lower() == ".tsv" else ","
    return pd.read_csv(path, sep=separator, keep_default_na=False)


def evaluate_file(
    prediction_path: Path,
    *,
    system: str,
    dataset: str,
    output_path: Path,
    ontology: str | Path | CellOntology | None = None,
) -> BenchmarkResult:
    """Score a prediction table from disk and publish the result immutably.

    An ``ontology`` given as a path is loaded here, so the published result can
    record the snapshot's resolved path and digest next to the prediction
    table's, in the same way it records what it scored. A measurement that
    cannot be traced to the ontology it was taken against is not reproducible,
    and a benchmark artifact that claims to be is worse than one that declines
    to measure.
    """

    if isinstance(ontology, CellOntology):
        resolved: CellOntology | None = ontology
    else:
        resolved = None if ontology is None else CellOntology.load(ontology)
    frame = read_predictions(prediction_path)
    result = evaluate_predictions(frame, system=system, dataset=dataset, ontology=resolved)
    source = None if resolved is None else Path(resolved.source)
    payload = {
        **result.model_dump(mode="json"),
        "prediction_path": str(prediction_path.resolve()),
        "prediction_sha256": sha256(prediction_path),
        "neutral_format": True,
        "ontology_source": None if source is None else str(source.resolve()),
        "ontology_sha256": None if source is None or not source.is_file() else sha256(source),
    }
    publish_json(output_path, payload)
    return result


def validate_dataset_manifest(path: Path) -> dict[str, Any]:
    try:
        manifest = BenchmarkDatasetManifest.model_validate_json(path.read_text())
    except (OSError, ValidationError, json.JSONDecodeError) as exc:
        raise ContractError(f"invalid benchmark manifest {path}: {exc}") from exc
    value = manifest.model_dump(mode="json", exclude_none=True)
    if manifest.status == "available":
        local = Path(manifest.local_path or "")
        if not local.is_absolute() and not local.is_file():
            # Manifests ship under <package-data>/benchmarks/manifests while
            # local_path is rooted at <package-data> (and at the repository
            # root during development). This keeps the same manifest portable.
            packaged = path.resolve().parents[2] / local
            if packaged.is_file():
                local = packaged
        if not local.is_file() or sha256(local) != manifest.sha256:
            raise ContractError("available benchmark data are absent or fail checksum")
    return value


def write_pending_status(
    output_path: Path,
    *,
    dataset: str,
    modality: str,
    reason: str,
    reproduce: list[str],
) -> None:
    publish_json(
        output_path,
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_benchmark_status",
            "dataset": dataset,
            "modality": modality,
            "status": "pending-external-resource",
            "reason": reason,
            "reproduce": reproduce,
            "superiority_claimed": False,
        },
    )

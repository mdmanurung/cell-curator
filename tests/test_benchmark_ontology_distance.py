"""`ontology_distance_error` has to be an ontology measurement or nothing at all.

The field shipped as a second name for the slash-path hierarchy error: one
number was written into both, so every published result carried two apparently
independent ontology metrics of which exactly one existed. In a package whose
README says comparative claims require measured results, that is the failure
mode the instrument is built to prevent -- a claim under a name that no code
ever checked.

The value in the committed synthetic result, 0.0625, therefore came from the
alias and is not reproducible now that the alias is gone. It cannot be made
reproducible again without minting a number for a measurement nobody performed,
so the artifact has to be superseded rather than re-derived. That decision is
recorded in
`test_the_committed_ontology_value_came_from_the_removed_alias` below, which is
the one place in the suite that knows about it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from conftest import REPO_ROOT

from cell_curator.benchmark import evaluate_file, evaluate_predictions
from cell_curator.cli import build_parser, dispatch
from cell_curator.contracts import ContractError, sha256
from cell_curator.ontology import CellOntology, OntologyTerm

COMMITTED = REPO_ROOT / "benchmarks" / "results" / "synthetic.json"
FIXTURE = REPO_ROOT / "benchmarks" / "fixtures" / "synthetic_predictions.tsv"

B_CELL = "CL:0000236"
PLASMA_CELL = "CL:0000786"
LYMPHOCYTE = "CL:0000542"
T_CELL = "CL:0000084"
LEUKOCYTE = "CL:0000738"
FIBROBLAST = "CL:0000057"
DETACHED = "CL:9000001"
ABSENT_PARENT = "CL:9000000"
ABSENT_TERM = "CL:1234567"


def _ontology() -> CellOntology:
    """A miniature Cell Ontology: one immune branch and one somatic branch.

    The parent of `DETACHED` is deliberately not in the snapshot, which is what
    a subtree export of an ontology looks like and what makes a root of the
    supplied graph a term rather than a fact about the Cell Ontology.
    """

    terms = {
        "CL:0000000": OntologyTerm("CL:0000000", "cell"),
        LEUKOCYTE: OntologyTerm(LEUKOCYTE, "leukocyte", parents=("CL:0000000",)),
        LYMPHOCYTE: OntologyTerm(LYMPHOCYTE, "lymphocyte", parents=(LEUKOCYTE,)),
        T_CELL: OntologyTerm(T_CELL, "T cell", parents=(LYMPHOCYTE,)),
        B_CELL: OntologyTerm(B_CELL, "B cell", parents=(LYMPHOCYTE,)),
        PLASMA_CELL: OntologyTerm(PLASMA_CELL, "plasma cell", parents=(B_CELL,)),
        FIBROBLAST: OntologyTerm(FIBROBLAST, "fibroblast", parents=("CL:0000000",)),
        DETACHED: OntologyTerm(DETACHED, "detached term", parents=(ABSENT_PARENT,)),
    }
    return CellOntology(terms, source="fixture")


def _cl_frame(true_ids: list[str], predicted_ids: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "barcode": [f"b{index}" for index in range(len(true_ids))],
            "true_label": ["B cell"] * len(true_ids),
            "predicted_label": ["B cell"] * len(true_ids),
            "true_cl_id": true_ids,
            "predicted_cl_id": predicted_ids,
        }
    )


def _distance(true_id: str, predicted_id: str) -> float | None:
    return evaluate_predictions(
        _cl_frame([true_id], [predicted_id]),
        system="fixture",
        dataset="cl-distance",
        ontology=_ontology(),
    ).ontology_distance_error


def test_a_cl_descendant_scores_better_than_an_unrelated_term() -> None:
    """The point of an ontology distance: being near is not the same as being wrong.

    Reading the ladder for a B cell as truth: the exact term is free, a
    descendant is charged only for the depth it adds, an ancestor is charged for
    the depth it drops, a sibling branch is charged for taking the wrong branch,
    and a term placed outside the immune branch is charged in full. A metric that
    could not order these would not be measuring semantics, and one that scored
    them all alike would be the alias again wearing a different name.

    The last rung is full only because this miniature hangs the somatic branch
    directly off the root. `test_the_metric_does_not_saturate_on_a_connected_ontology`
    pins what the same pair measures on a graph shaped like the real Cell
    Ontology, which is not 1.0.
    """

    frame = _cl_frame(
        [B_CELL] * 5,
        [B_CELL, PLASMA_CELL, LYMPHOCYTE, T_CELL, FIBROBLAST],
    )
    result = evaluate_predictions(
        frame, system="fixture", dataset="cl-distance", ontology=_ontology()
    )
    ladder = [
        _distance(B_CELL, predicted)
        for predicted in (B_CELL, PLASMA_CELL, LYMPHOCYTE, T_CELL, FIBROBLAST)
    ]
    assert ladder == [
        pytest.approx(0.0),
        pytest.approx(1 / 4),
        pytest.approx(1 / 3),
        pytest.approx(0.5),
        1.0,
    ]
    assert ladder == sorted(ladder)
    assert result.ontology_distance_error == pytest.approx(sum(ladder) / 5)


def test_agreement_near_the_root_is_worth_more_than_agreement_near_the_leaves() -> None:
    """The shared depth belongs in the denominator, and this is why.

    Being too coarse by one level and being too coarse by two are both one step
    away from the truth, yet the second is a blunter claim. Charging the first at
    a third and the second at two thirds is what stops a system from buying
    apparent accuracy by collapsing everything to a high-level ancestor.
    """

    assert _distance(B_CELL, LYMPHOCYTE) == pytest.approx(1 / 3)
    assert _distance(B_CELL, LEUKOCYTE) == pytest.approx(2 / 3)
    assert _distance(B_CELL, LEUKOCYTE) > _distance(B_CELL, LYMPHOCYTE)


def test_a_term_sharing_no_ancestor_scores_the_full_error() -> None:
    """No shared ancestry means the is_a graph supplies no evidence of relation.

    Two ways of arriving at that: a term whose only common ancestry with the
    truth is a root of the supplied graph, and a term in a different component
    of it. Neither can be placed on a path between them, so the specification's
    1.0 is the honest answer and anything between 0 and 1 would be invented.

    Only the second of the two survives contact with a real ontology; see
    `test_the_metric_does_not_saturate_on_a_connected_ontology`.
    """

    ontology = _ontology()
    frame = _cl_frame([B_CELL, B_CELL], [FIBROBLAST, DETACHED])
    result = evaluate_predictions(
        frame, system="fixture", dataset="cl-unrelated", ontology=ontology
    )
    assert result.ontology_distance_error == 1.0
    assert ontology.roots() == frozenset({"CL:0000000", DETACHED})


def test_the_metric_does_not_saturate_on_a_connected_ontology() -> None:
    """The bound a real Cell Ontology result has to be read against.

    Every pair of real CL terms shares `cell` and `native cell`, and neither is a
    root of the graph, so the deepest shared term always sits at least one edge
    down and the score stays below 1.0 however wrong the prediction is. A B cell
    called a fibroblast is therefore charged about 0.8, not 1.0, because the two
    really are both native cells.

    This is why `ontology_distance_error` orders systems against each other on
    one snapshot and is not a fraction of errors made. An implementation that
    reached 1.0 here would be discarding genuine shared ancestry to make a
    number look finished.
    """

    native_cell = "CL:0000003"
    realistic = CellOntology(
        {
            "CL:0000000": OntologyTerm("CL:0000000", "cell"),
            native_cell: OntologyTerm(native_cell, "native cell", ("CL:0000000",)),
            LEUKOCYTE: OntologyTerm(LEUKOCYTE, "leukocyte", (native_cell,)),
            LYMPHOCYTE: OntologyTerm(LYMPHOCYTE, "lymphocyte", (LEUKOCYTE,)),
            B_CELL: OntologyTerm(B_CELL, "B cell", (LYMPHOCYTE,)),
            FIBROBLAST: OntologyTerm(FIBROBLAST, "fibroblast", (native_cell,)),
        },
        source="realistic",
    )
    shared = (
        (realistic.ancestors(B_CELL) | {B_CELL})
        & (realistic.ancestors(FIBROBLAST) | {FIBROBLAST})
    ) - realistic.roots()
    assert shared == {native_cell}

    distance = evaluate_predictions(
        _cl_frame([B_CELL], [FIBROBLAST]),
        system="fixture",
        dataset="cl-connected",
        ontology=realistic,
    ).ontology_distance_error
    assert distance == pytest.approx(4 / 5)
    assert distance < 1.0
    assert realistic.roots() == frozenset({"CL:0000000"})


def test_a_table_without_ontology_ids_reports_no_measurement() -> None:
    """The synthetic fixture's case, and the one the alias used to paper over.

    `hierarchy_error` is a real measurement of slash-separated label paths, so it
    is present and non-zero here. `ontology_distance_error` has no evidence
    behind it in this table at all, so it must be absent. A copied hierarchy
    number would be indistinguishable from a measurement to every reader of the
    published artifact, which is the entire defect.
    """

    frame = pd.read_csv(FIXTURE, sep="\t", keep_default_na=False)
    assert "true_cl_id" not in frame.columns
    assert "true_cell_type_ontology_term_id" not in frame.columns
    result = evaluate_predictions(frame, system="fixture", dataset="synthetic")
    assert result.hierarchy_error == 0.0625
    assert result.hierarchical_accuracy == pytest.approx(0.9375)
    assert result.ontology_distance_error is None


def test_supplying_an_ontology_does_not_invent_a_measurement() -> None:
    """Constructing the snapshot is the caller's choice, not the metric's.

    The ontology is optional so that a plain evaluation does not depend on a
    download, but supplying one is not itself permission to report a number. If
    the table makes no ontology claim, the result says the same thing either way.
    """

    frame = pd.read_csv(FIXTURE, sep="\t", keep_default_na=False)
    without = evaluate_predictions(frame, system="fixture", dataset="synthetic")
    with_ontology = evaluate_predictions(
        frame, system="fixture", dataset="synthetic", ontology=_ontology()
    )
    assert with_ontology.ontology_distance_error is None
    assert with_ontology.model_dump(mode="json") == without.model_dump(mode="json")


def test_ontology_ids_without_a_snapshot_fail_closed() -> None:
    """A table that claims to be scored against CL cannot be scored unreported.

    Returning `None` here would be indistinguishable in the published artifact
    from a table that made no claim, so a caller could ship an ontology
    comparison that was never performed. The same reasoning makes
    `validate_durable_label_ontology` refuse to write a CL id without a
    validated snapshot, and it is why the error names the columns at fault.
    """

    frame = _cl_frame([B_CELL], [T_CELL])
    with pytest.raises(ContractError, match="true_cl_id"):
        evaluate_predictions(frame, system="fixture", dataset="cl-unbacked")


def test_the_census_column_spelling_is_accepted_too() -> None:
    """A Census export names the same field differently and means the same thing.

    Rejecting the alternative spelling would push an external table to be
    renamed before it could be scored, and the renamed copy would no longer be
    the file the system published.
    """

    frame = pd.DataFrame(
        {
            "barcode": ["b0"],
            "true_label": ["B cell"],
            "predicted_label": ["lymphocyte"],
            "true_cell_type_ontology_term_id": [B_CELL],
            "predicted_cell_type_ontology_term_id": [LYMPHOCYTE],
        }
    )
    result = evaluate_predictions(
        frame, system="fixture", dataset="cl-census-naming", ontology=_ontology()
    )
    assert result.ontology_distance_error == pytest.approx(1 / 3)


def test_an_abstention_carries_no_pair_and_is_not_charged_full_error() -> None:
    """A blank predicted term is not a maximally wrong term.

    An `Unknown` prediction has no CL term to be near or far from, and this
    instrument already scores abstention through the Unknown precision and
    recall. Charging a row the maximum distance as well would count one mistake
    twice, so blank identifiers are excluded and the reported value is the mean
    over the rows that carry a pair. A table whose identifiers are all blank
    therefore has nothing to measure and says so.
    """

    ontology = _ontology()
    frame = _cl_frame([B_CELL, B_CELL, B_CELL], [T_CELL, "", "  "])
    result = evaluate_predictions(
        frame, system="fixture", dataset="cl-abstention", ontology=ontology
    )
    assert result.ontology_distance_error == 0.5
    blank = _cl_frame([B_CELL, B_CELL], ["", ""])
    assert (
        evaluate_predictions(
            blank, system="fixture", dataset="cl-abstention", ontology=ontology
        ).ontology_distance_error
        is None
    )


def test_an_identifier_the_snapshot_does_not_know_is_refused() -> None:
    """An unresolvable CL id must not be scored as maximally distant.

    The alternative reading -- "not in the ontology" means "as far as possible"
    -- would let a typo silently become the worst score in a published result.
    """

    frame = _cl_frame([B_CELL], [ABSENT_TERM])
    with pytest.raises(ContractError, match="ontology term is unknown"):
        evaluate_predictions(
            frame, system="fixture", dataset="cl-unknown", ontology=_ontology()
        )


def test_a_published_result_records_the_snapshot_it_measured_against(
    tmp_path: Path,
) -> None:
    """A measurement nobody can trace to a file is not reproducible.

    The result already pins the prediction table's path and digest. An ontology
    distance is only meaningful against a named snapshot, and the Cell Ontology
    is a moving target, so the digest is what lets a later reader tell whether
    the number can be reproduced at all. The snapshot is loaded from a table
    here, which is the same code path the `--ontology` flag takes.
    """

    ontology_path = _write_ontology_table(tmp_path)
    predictions = _write_cl_predictions(tmp_path)
    output = tmp_path / "result.json"
    result = evaluate_file(
        predictions,
        system="fixture",
        dataset="cl-provenance",
        output_path=output,
        ontology=ontology_path,
    )
    assert result.ontology_distance_error == pytest.approx((0.5 + 1.0) / 2)
    payload = json.loads(output.read_text())
    assert payload["ontology_source"] == str(ontology_path.resolve())
    assert payload["ontology_sha256"] == sha256(ontology_path)
    assert payload["prediction_sha256"] == sha256(predictions)
    assert payload["neutral_format"] is True


def test_a_result_without_an_ontology_publishes_no_ontology_provenance(
    tmp_path: Path,
) -> None:
    """The absent snapshot is recorded as absent rather than omitted.

    A reader comparing two results needs to tell "measured against no ontology"
    apart from "measured against a snapshot that has since been deleted", and
    the second is the one that invalidates a claim.
    """

    output = tmp_path / "result.json"
    predictions = tmp_path / "plain_predictions.tsv"
    predictions.write_text(
        "barcode\ttrue_label\tpredicted_label\nb0\tB cell\tB cell\nb1\tB cell\tT cell\n"
    )
    evaluate_file(
        predictions,
        system="fixture",
        dataset="cl-provenance",
        output_path=output,
    )
    payload = json.loads(output.read_text())
    assert payload["ontology_source"] is None
    assert payload["ontology_sha256"] is None
    assert payload["ontology_distance_error"] is None


def test_the_committed_ontology_value_came_from_the_removed_alias() -> None:
    """The one value in the suite that cannot be re-derived, and why.

    `benchmarks/results/synthetic.json` is immutable and records
    `ontology_distance_error` 0.0625, which is `hierarchy_error` under a second
    name. No code path produces that number any more, so
    `test_the_committed_result_is_reproducible_from_its_fixture` fails on this
    key alone while the other 26 metrics still re-derive exactly. Re-deriving it
    would require restoring the alias, and the only way to make the artifact
    true is to supersede it: re-publish the result with the metric absent, or
    ship a fixture that carries CL identifiers and a matching snapshot. Either
    way the change is a deliberate, committed act, which is exactly what
    `publish_json` immutability is for.
    """

    committed = json.loads(COMMITTED.read_text())
    assert committed["hierarchy_error"] == 0.0625
    assert committed["ontology_distance_error"] == 0.0625
    frame = pd.read_csv(FIXTURE, sep="\t", keep_default_na=False)
    assert (
        evaluate_predictions(frame, system="fixture", dataset="synthetic").ontology_distance_error
        is None
    )


def test_the_cli_flag_reaches_the_evaluator(tmp_path: Path) -> None:
    """The snapshot is reachable from the command line, not only from Python.

    A comparative claim is usually checked by someone running the published
    command, and the flag is the only route by which a Cell Ontology can reach it.
    Omitting it would leave the metric reachable only from a notebook, which is
    how a claim quietly stops being reproducible for everyone else.
    """

    arguments = build_parser().parse_args(
        [
            "benchmark",
            "evaluate",
            "--predictions",
            str(_write_cl_predictions(tmp_path)),
            "--system",
            "fixture",
            "--dataset",
            "cl-cli",
            "--output",
            str(tmp_path / "cli-result.json"),
            "--ontology",
            str(_write_ontology_table(tmp_path)),
        ]
    )
    reported = dispatch(arguments)
    assert reported["ontology_distance_error"] == pytest.approx((0.5 + 1.0) / 2)
    published = json.loads((tmp_path / "cli-result.json").read_text())
    assert published["ontology_source"] == str(
        _write_ontology_table(tmp_path).resolve()
    )
    assert published["ontology_sha256"] == sha256(_write_ontology_table(tmp_path))


def _write_ontology_table(directory: Path) -> Path:
    path = directory / "mini_cl.tsv"
    path.write_text(
        "id\tname\tparents\n"
        "CL:0000000\tcell\t\n"
        f"{LEUKOCYTE}\tleukocyte\tCL:0000000\n"
        f"{LYMPHOCYTE}\tlymphocyte\t{LEUKOCYTE}\n"
        f"{T_CELL}\tT cell\t{LYMPHOCYTE}\n"
        f"{B_CELL}\tB cell\t{LYMPHOCYTE}\n"
        f"{FIBROBLAST}\tfibroblast\tCL:0000000\n"
    )
    return path


def _write_cl_predictions(directory: Path) -> Path:
    path = directory / "predictions.tsv"
    path.write_text(
        "barcode\ttrue_label\tpredicted_label\ttrue_cl_id\tpredicted_cl_id\n"
        f"b0\tB cell\tT cell\t{B_CELL}\t{T_CELL}\n"
        f"b1\tB cell\tfibroblast\t{B_CELL}\t{FIBROBLAST}\n"
    )
    return path

"""The published neutral-format contract has to name the columns the code reads.

`benchmarks/neutral-format.schema.json` described a table's required columns and
then said only that some optional columns might exist, without naming one. Every
manifest in `benchmarks/` tells an exporter to follow that file, so an exporter
following it in good faith had no way to know that `confidence` was the column
the risk-coverage sweep is computed from, or that a CL identifier column turns
on a metric the evaluator will refuse to report without a local ontology
snapshot. The header of the shipped fixture carried fourteen optional columns and
the contract accounted for none of them.

Documentation of a wire format drifts silently, so the enumeration is pinned here
in both directions against the two sources that actually decide it: the fixture
header, which shows what an export looks like, and `benchmark.py`, which decides
what is scored. A new scored column that nobody documents, or a documented column
the evaluator no longer reads, both fail.
"""

from __future__ import annotations

import json

from conftest import REPO_ROOT

SCHEMA = REPO_ROOT / "benchmarks" / "neutral-format.schema.json"
FIXTURE = REPO_ROOT / "benchmarks" / "fixtures" / "synthetic_predictions.tsv"
BENCHMARK_SOURCE = REPO_ROOT / "src" / "cell_curator" / "benchmark.py"

REQUIRED = ("barcode", "true_label", "predicted_label")


def _schema() -> dict:
    return json.loads(SCHEMA.read_text())


def _fixture_columns() -> list[str]:
    return FIXTURE.read_text().splitlines()[0].split("\t")


def _scored() -> list[str]:
    return _schema()["properties"]["scored_optional_columns"]["items"]["enum"]


def _unscored() -> list[str]:
    return _schema()["properties"]["carried_unscored_columns"]["items"]["enum"]


def test_every_optional_column_in_the_shipped_export_is_accounted_for() -> None:
    """A column the contract does not name is a column nobody can rely on.

    The fixture is the only worked example of the format, so its header is the
    closest thing to the contract's subject matter. Anything optional in it must
    be either scored or explicitly declared carried-and-unscored, because a
    reader who finds an undeclared column has to assume it is a measurement and
    it is not.
    """

    optional = [column for column in _fixture_columns() if column not in REQUIRED]
    declared = set(_scored()) | set(_unscored())
    assert [column for column in optional if column not in declared] == []
    assert set(_scored()) | set(_unscored()) == set(optional) | {
        "true_cl_id",
        "predicted_cl_id",
        "true_cell_type_ontology_term_id",
        "predicted_cell_type_ontology_term_id",
    }


def test_every_scored_column_is_one_the_evaluator_actually_reads() -> None:
    """A promised column that is silently ignored is worse than an undocumented one.

    An exporter following the contract would reasonably conclude that a column
    it omits costs it a metric. If the evaluator in fact stopped reading that
    column, the contract would be charging for nothing.
    """

    source = BENCHMARK_SOURCE.read_text()
    assert [column for column in _scored() if f'"{column}"' not in source] == []


def test_every_scored_column_says_which_metric_it_moves() -> None:
    """Naming a column without naming its metric is half an enumeration.

    The mapping is also the only place a reader learns that `confidence` is what
    makes the risk-coverage curve computable, and that a CL identifier column
    commits the result to an ontology snapshot.
    """

    metrics = _schema()["properties"]["scored_optional_column_metrics"]["properties"]
    assert list(metrics) == _scored()
    assert metrics["confidence"]["const"] == [
        "brier_score",
        "expected_calibration_error",
    ]
    assert all(
        value["const"] for value in metrics.values() if isinstance(value, dict)
    )
    ontology = {
        column: value["const"]
        for column, value in metrics.items()
        if column in {"true_cl_id", "predicted_cl_id"}
        or column.endswith("_cell_type_ontology_term_id")
    }
    assert {tuple(value) for value in ontology.values()} == {
        ("ontology_distance_error", "descendant_accuracy", "relaxed_accuracy")
    }
    assert len(ontology) == 4


def test_the_required_columns_still_match_what_the_evaluator_demands() -> None:
    """The contract's floor is the evaluator's floor, restated rather than copied.

    The evaluator's own guard lives in code, and this file is what an exporter
    reads before writing a table. A divergence between the two would fail a
    correct export at the worst moment, so the two are compared rather than
    trusted.
    """

    from cell_curator.benchmark import REQUIRED_COLUMNS

    assert list(_schema()["properties"]["required_columns"]["const"]) == list(REQUIRED)
    assert set(REQUIRED_COLUMNS) == set(REQUIRED)


def test_the_document_shape_itself_is_unchanged() -> None:
    """Additive documentation only: a consumer's document still validates.

    The description object gained keys but kept the three it required before, so
    an exporter that wrote the old shape still validates and a reader's tooling
    that reads `required_columns` still finds it.
    """

    schema = _schema()
    assert schema["type"] == "object"
    assert schema["required"] == ["schema_version", "required_columns", "optional_columns"]
    assert schema["properties"]["schema_version"] == {"const": 3}
    assert schema["properties"]["optional_columns"]["items"] == {"type": "string"}


def test_the_ontology_pair_is_required_to_be_supplied_as_a_pair() -> None:
    """Half a CL column is a claim the evaluator cannot check.

    The contract states the pairing for the evaluator's sake: a table carrying
    only `true_cl_id` makes no measurement and a table carrying only
    `predicted_cl_id` does not either, so pairing them in the schema is how an
    exporter learns it before shipping an unscoreable table.
    """

    rules = _schema()["properties"]["scored_optional_columns"]["allOf"]
    assert any("true_path" in json.dumps(rule) for rule in rules)
    assert any("predicted_path" in json.dumps(rule) for rule in rules)
    assert any("true_cl_id" in json.dumps(rule) for rule in rules)
    assert any("predicted_cl_id" in json.dumps(rule) for rule in rules)


def test_every_public_manifest_points_an_exporter_at_the_enumerated_contract() -> None:
    """The five public manifests are how a third party learns the format.

    Each names this schema in its instructions, so the enumerated columns reach
    an external system only while those instructions survive. The synthetic
    manifest is the one exception and is not one: it describes a fixture that
    already exists on disk rather than data an external system must export, so
    it has nothing to instruct.
    """

    manifests = REPO_ROOT / "benchmarks" / "manifests"
    pointing = {
        path.stem
        for path in manifests.glob("*.json")
        if "benchmarks/neutral-format.schema.json" in path.read_text()
    }
    assert pointing == {"atac", "cite-seq", "multiome", "rna", "spatial"}
    fixture_manifest = json.loads((manifests / "synthetic.json").read_text())
    assert fixture_manifest["status"] == "available"

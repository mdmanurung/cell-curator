"""The benchmark instrument is itself held to account before any number is read.

`benchmarks/results/synthetic.json` is a committed, immutable, claim-bearing
artifact and `test_clean_wheel.py` only asserts the packaging *path*. Nothing
recomputed it, so an edit to `benchmark.py` or to the fixture it scores would
leave CI green while the shipped number silently became false. That is the
worst failure mode available to a package whose README states that comparative
claims require measured results: the measurement is unfalsifiable.

These tests pin the instrument, so a later real-data result is a measurement
rather than a reading of whatever the code happened to emit that day.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cell_curator.benchmark import (
    evaluate_file,
)
from cell_curator.contracts import ContractError, sha256

REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "benchmarks" / "fixtures" / "synthetic_predictions.tsv"
COMMITTED = REPO / "benchmarks" / "results" / "synthetic.v2.json"
SUPERSEDED = REPO / "benchmarks" / "results" / "synthetic.json"
SUPERSESSION = REPO / "benchmarks" / "results" / "synthetic.supersession.json"
MANIFEST = REPO / "benchmarks" / "manifests" / "synthetic.json"


def test_the_committed_result_is_reproducible_from_its_fixture(tmp_path: Path) -> None:
    """Recomputing the shipped result must reproduce every metric in it, exactly.

    The artifact stores file names rather than resolved paths, so it is the same
    on every machine and the comparison covers every key the scorer returns.
    """

    committed = json.loads(COMMITTED.read_text())
    recomputed = evaluate_file(
        FIXTURE,
        system=committed["system"],
        dataset=committed["dataset"],
        output_path=tmp_path / "recomputed.json",
    ).model_dump(mode="json")

    for key, value in recomputed.items():
        assert committed[key] == value, f"{key} drifted from the committed result"
    assert committed["prediction_path"] == FIXTURE.name
    assert committed["neutral_format"] is True


def test_the_superseded_result_is_kept_and_the_supersession_names_both() -> None:
    """Superseding is a recorded act: the old bytes stay, and a record says why."""

    record = json.loads(SUPERSESSION.read_text())
    assert record["supersedes"] == "benchmarks/results/synthetic.json"
    assert record["superseded_by"] == "benchmarks/results/synthetic.v2.json"
    assert SUPERSEDED.is_file()
    assert json.loads(SUPERSEDED.read_text())["prediction_sha256"] == record["prediction_sha256"]


def test_the_result_manifest_and_fixture_all_name_the_same_bytes() -> None:
    """Three files carry the same digest; nothing detected a divergence before."""

    digest = sha256(FIXTURE)
    assert json.loads(COMMITTED.read_text())["prediction_sha256"] == digest
    assert json.loads(MANIFEST.read_text())["sha256"] == digest


def test_the_committed_metrics_are_the_edge_cases_the_fixture_documents() -> None:
    """Pin the values themselves, not just their reproducibility.

    Re-derivation proves the code still produces this file. It does not prove
    the file is *right*. The fixture's one wrong answer is a deliberate
    abstention at confidence 0.42, so selective accuracy must be perfect while
    plain accuracy is not; a harness that lost that distinction would still
    reproduce the bytes.
    """

    committed = json.loads(COMMITTED.read_text())
    assert committed["n"] == 8
    assert committed["accuracy"] == 0.875
    assert committed["coverage"] == 0.625
    assert committed["selective_accuracy"] == 1.0
    assert committed["coverage_risk"] == 0.0
    assert committed["unknown_recall"] == 1.0
    assert committed["unknown_precision"] == pytest.approx(2 / 3)
    assert committed["hierarchy_error"] == 0.0625
    # Calibration over calls only: five called rows, abstentions excluded.
    assert committed["brier_score"] == pytest.approx(
        (0.05**2 + 0.08**2 + 0.1**2 + 0.14**2 + 0.11**2) / 5
    )
    assert committed["ontology_distance_error"] is None
    assert committed["per_type"]["T cell"] == {
        "precision": 1.0,
        "recall": 0.5,
        "f1": pytest.approx(2 / 3),
        "support": 2,
    }


def test_evaluate_file_chooses_its_separator_from_the_suffix(tmp_path: Path) -> None:
    """`.tsv` and everything else are read on different separators.

    A comma-separated table named `.tsv` reads as a single column and fails on
    a missing required column, which is a confusing error for a real mistake.
    """

    row = "barcode\ttrue_label\tpredicted_label\nx\tT cell\tT cell\n"
    tsv = tmp_path / "predictions.tsv"
    tsv.write_text(row)
    assert evaluate_file(tsv, system="s", dataset="d", output_path=tmp_path / "a.json").n == 1

    csv = tmp_path / "predictions.csv"
    csv.write_text(row.replace("\t", ","))
    assert evaluate_file(csv, system="s", dataset="d", output_path=tmp_path / "b.json").n == 1


def test_evaluate_file_refuses_to_rewrite_a_result_with_new_numbers(tmp_path: Path) -> None:
    """A published result is a commitment, so the second disagreement fails."""

    output = tmp_path / "result.json"
    evaluate_file(FIXTURE, system="s", dataset="d", output_path=output)

    evaluate_file(FIXTURE, system="s", dataset="d", output_path=output)
    assert json.loads(output.read_text())["n"] == 8

    other = tmp_path / "other.tsv"
    other.write_text("barcode\ttrue_label\tpredicted_label\nz\tB cell\tT cell\n")
    with pytest.raises(ContractError, match="immutable"):
        evaluate_file(other, system="s", dataset="d", output_path=output)


def test_the_scrambled_8000_gene_result_stays_retracted() -> None:
    """The v10 result ran on a row-scrambled input; its retraction cannot vanish."""

    retraction = json.loads((COMMITTED.parent / "aida-pbmc-v10.retraction.json").read_text())
    assert retraction["input_sha256"].startswith("a8c5e332")
    for path in retraction["retracts"]:
        assert (COMMITTED.parents[2] / path).is_file(), path

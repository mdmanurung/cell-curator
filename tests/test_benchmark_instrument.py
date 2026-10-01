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

import pandas as pd
import pytest

from cell_curator.benchmark import (
    evaluate_file,
    evaluate_predictions,
    write_pending_status,
)
from cell_curator.contracts import ContractError, sha256

REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "benchmarks" / "fixtures" / "synthetic_predictions.tsv"
COMMITTED = REPO / "benchmarks" / "results" / "synthetic.json"
MANIFEST = REPO / "benchmarks" / "manifests" / "synthetic.json"


def test_the_committed_result_is_reproducible_from_its_fixture(tmp_path: Path) -> None:
    """Recomputing the shipped result must reproduce every metric in it.

    Only `prediction_path` is excluded. `evaluate_file` bakes an absolute
    resolved machine path into the payload, so that field is by construction
    different on every machine and cannot participate in a cross-machine pin.
    Everything that is a measurement is asserted.

    `ontology_distance_error`, `descendant_accuracy`, and `relaxed_accuracy` are
    pinned to their pending values rather than to what the code returns. `ontology_distance_error`
    recorded 0.0625, which was `hierarchy_error` under a second name; the alias
    is gone and a table with no Cell Ontology ids makes no ontology claim, so
    the honest answer is `None`. `relaxed_accuracy` did not exist when the
    artifact was published, so the artifact predates it. Superseding the
    immutable artifact is deferred until the first real benchmark result lands,
    and `test_the_committed_artifact_still_awaiting_supersession` is the strict
    marker that fails the moment that deferral is resolved either way.
    """

    committed = json.loads(COMMITTED.read_text())
    recomputed = evaluate_file(
        FIXTURE,
        system=committed["system"],
        dataset=committed["dataset"],
        output_path=tmp_path / "recomputed.json",
    ).model_dump(mode="json")

    for key, value in recomputed.items():
        if key in {"ontology_distance_error", "descendant_accuracy", "relaxed_accuracy"}:
            continue
        assert committed[key] == value, f"{key} drifted from the committed result"
    assert committed["prediction_path"] == str(FIXTURE.resolve())
    assert committed["neutral_format"] is True
    assert committed["ontology_distance_error"] == 0.0625
    assert recomputed["ontology_distance_error"] is None
    for key in ("descendant_accuracy", "relaxed_accuracy"):
        assert key not in committed
        assert recomputed[key] is None


@pytest.mark.xfail(
    strict=True,
    reason=(
        "the committed synthetic result records ontology_distance_error from the "
        "alias that was removed once the metric became a real is_a-graph "
        "measurement; superseding the immutable artifact is deferred until the "
        "first real benchmark result lands"
    ),
)
def test_the_committed_artifact_still_awaiting_supersession(tmp_path: Path) -> None:
    """Fail the day the deferral is resolved, so it cannot be forgotten.

    This asserts what the artifact will have to look like once it is superseded:
    every metric re-derives from the fixture, `ontology_distance_error`
    included. It is `strict`, so when the artifact is replaced with a result
    carrying a measured ontology distance this turns into a failure and has to
    be removed deliberately, at which point the pin above is tightened back to
    comparing every key.
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
    assert committed["brier_score"] == pytest.approx(0.0337875)
    assert committed["expected_calibration_error"] == pytest.approx(0.14875)


def test_evaluate_file_chooses_its_separator_from_the_suffix(tmp_path: Path) -> None:
    """`.tsv` and everything else are read on different separators.

    A comma-separated table named `.tsv` reads as a single column and fails on
    a missing required column, which is a confusing error for a real mistake.
    """

    row = "barcode\ttrue_label\tpredicted_label\nx\tT cell\tT cell\n"
    tsv = tmp_path / "predictions.tsv"
    tsv.write_text(row)
    assert evaluate_file(
        tsv, system="s", dataset="d", output_path=tmp_path / "a.json"
    ).n == 1

    csv = tmp_path / "predictions.csv"
    csv.write_text(row.replace("\t", ","))
    assert evaluate_file(
        csv, system="s", dataset="d", output_path=tmp_path / "b.json"
    ).n == 1


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


def test_a_pending_status_is_not_a_measurement_and_never_claims_superiority(
    tmp_path: Path,
) -> None:
    """Every real run writes pending markers for all five public modalities.

    `pipeline.py` emits one per run, so a reader scanning `benchmarks/` for a
    score can land on one of these. It must be structurally impossible to
    mistake for a `BenchmarkResult`, and it must never assert a comparison.
    """

    output = tmp_path / "public_rna.pending.json"
    reason = "Public dataset and external comparator exports were not supplied."
    reproduce = ["curl -fL --output matrix.h5 https://example.invalid/matrix.h5"]
    write_pending_status(
        output, dataset="10x-pbmc", modality="rna", reason=reason, reproduce=reproduce
    )

    payload = json.loads(output.read_text())
    assert payload["status"] == "pending-external-resource"
    assert payload["superiority_claimed"] is False
    assert payload["reason"] == reason
    assert payload["reproduce"] == reproduce

    measured = evaluate_predictions(_frame(["T cell"]), system="s", dataset="d")
    envelope = {"schema_version", "system", "dataset", "n", "artifact_path"}
    reported = {
        key
        for key, value in measured.model_dump(mode="json").items()
        if key not in envelope and isinstance(value, (int, float))
    }
    assert reported, "the measurement surface must not be empty"
    assert not reported & payload.keys(), "a pending marker must report no measurement"

    with pytest.raises(ContractError, match="immutable"):
        write_pending_status(
            output,
            dataset="10x-pbmc",
            modality="rna",
            reason="a different reason",
            reproduce=reproduce,
        )


def _frame(labels: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "barcode": [f"b{i}" for i in range(len(labels))],
            "true_label": labels,
            "predicted_label": labels,
        }
    )


def test_the_scrambled_8000_gene_result_stays_retracted() -> None:
    """The v10 result ran on a row-scrambled input; its retraction cannot vanish."""

    retraction = json.loads(
        (COMMITTED.parent / "aida-pbmc-v10.retraction.json").read_text()
    )
    assert retraction["input_sha256"].startswith("a8c5e332")
    for path in retraction["retracts"]:
        assert (COMMITTED.parents[2] / path).is_file(), path

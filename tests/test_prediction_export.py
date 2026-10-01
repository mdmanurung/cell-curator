"""The run-to-neutral-table join is the one step nobody owned, so it is pinned here.

`benchmark.py` is engine-blind on purpose, so nothing in the package could read
a run and write a prediction table. The only table in `benchmarks/fixtures/` was
typed by hand, which means the scored cell set was not a function of anything
hashable: a cell dropped between the run and the score left no trace, and the
artifact could not be re-derived.

The contract this module establishes is fail-closed, and each refusal below
corresponds to a way the join could quietly produce a number nobody chose:

- truth that does not cover exactly the run's barcodes, in either direction,
  because the intersection of two sets is not a decision anyone made;
- a truth row repeated, because the harness scores one row per barcode and
  would otherwise average the duplicates in;
- a label carrying two different Cell Ontology identifiers, where the winner
  would be decided by row order and the metric would inherit that;
- a confidence outside [0, 1], which the calibration metric refuses and which
  would otherwise reach a published result.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from conftest import SKILL_ROOT

from cell_curator.config import CellCuratorConfig, load_config
from cell_curator.contracts import ContractError, sha256
from cell_curator.export import export_benchmark_predictions
from cell_curator.pipeline import run_annotation, run_root


def _run_artifacts(
    tmp_path: Path,
    *,
    mapping: pd.DataFrame,
    labels: pd.DataFrame | dict[str, pd.DataFrame],
) -> CellCuratorConfig:
    """Lay down the two artifacts the exporter reads, without running a pipeline.

    `labels` is a single frame for a single-level run, or a mapping of filename
    to frame for a recursive one, which is how the two cases are told apart.
    """

    root = tmp_path / "results" / "export-fixture-v1"
    (root / "mapping").mkdir(parents=True)
    (root / "labels").mkdir(parents=True)
    mapping.to_csv(
        root / "mapping" / "parent_child_mapping.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    tables = {name: frame for name, frame in labels.items()} if isinstance(labels, dict) else {
        "L1.tsv": labels
    }
    for name, frame in tables.items():
        frame.to_csv(root / "labels" / name, sep="\t", index=False, lineterminator="\n")
    value = yaml.safe_load((SKILL_ROOT / "assets" / "config.template.yaml").read_text())
    value["run"].update(
        {
            "run_id": "export-fixture-v1",
            "output_root": str(tmp_path / "results"),
            "mode": "evidence-only",
        }
    )
    value["input"]["path"] = str(tmp_path / "source.h5ad")
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(value, sort_keys=False))
    return load_config(path)


def _mapping(
    barcodes: list[str], types: list[str], clusters: list[str] | None = None
) -> pd.DataFrame:
    """A durable mapping. `clusters` assigns each cell its parent; all share one by default."""
    parents = clusters or ["0"] * len(barcodes)
    assert len(parents) == len(barcodes)
    return pd.DataFrame(
        {
            "schema_version": [3] * len(barcodes),
            "run_id": ["export-fixture-v1"] * len(barcodes),
            "barcode": barcodes,
            "scope": ["ROOT"] * len(barcodes),
            "canonical_parent": parents,
            "parent_uid": [f"ROOT::{value}" for value in parents],
            "leaf_id": [f"ROOT::{parent}::{value}" for parent, value in zip(parents, types, strict=True)],
            "cell_type": types,
            "cell_state": ["normal"] * len(barcodes),
            "membership_role": ["retained_parent"] * len(barcodes),
            "decision_outcome": ["RETAIN PARENT"] * len(barcodes),
        }
    )


def _labels(entries: list[tuple[str, str, str, float]], *, level: str = "L1") -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "level": level,
                "parent_scope": "ROOT",
                "cluster_id": cluster,
                "label": label,
                "cell_state": "normal",
                "cl_id": accession,
                "confidence": confidence,
                "margin": 0.4,
                "specificity": 0.5,
                "agreement": 0.5,
                "coverage": 0.5,
                "rationale": "fixture",
                "rejected_alternatives": "none measured",
                "reviewer": "",
            }
            for cluster, label, accession, confidence in entries
        ]
    )


def _truth(barcodes: list[str], labels: list[str], **extra: str) -> Path:
    path = Path(extra.pop("__path", "truth.tsv"))
    frame = pd.DataFrame({"barcode": barcodes, "true_label": labels, **extra})
    path.write_text(frame.to_csv(sep="\t", index=False, lineterminator="\n"))
    return path


def test_a_completed_export_carries_both_sides_of_every_comparable_column(
    tmp_path: Path,
) -> None:
    """The table has to be scoreable, and the receipt has to be attributable.

    Both halves matter. Without the CL identifiers the ontology metrics are
    absent, and without a receipt naming the truth and the mapping digests the
    published result cannot be traced to the artifacts it came from.
    """

    config = _run_artifacts(
        tmp_path,
        mapping=_mapping(
            ["c0", "c1", "c2", "c3", "c4"],
            ["T cell", "B cell", "T cell", "Unknown", "T cell"],
            ["0", "1", "0", "2", "0"],
        ),
        # Cluster 2 was decided against the T cell program but abstained.
        labels=_labels(
            [
                ("0", "T cell", "CL:0000084", 0.9),
                ("1", "B cell", "CL:0000236", 0.7),
                ("2", "Unknown", "CL:0000084", 0.4),
            ]
        ),
    )
    truth = _truth(
        ["c0", "c1", "c2", "c3", "c4"],
        ["T cell", "B cell", "Unknown", "T cell", "T cell"],
        true_cl_id=["CL:0000084", "CL:0000236", "", "CL:0000084", "CL:0000084"],
        true_path=[
            "immune/T cell",
            "immune/B cell",
            "immune",
            "immune/T cell",
            "immune/T cell",
        ],
        __path=tmp_path / "truth.tsv",
    )
    output = tmp_path / "predictions.tsv"

    receipt = export_benchmark_predictions(config, truth_path=truth, output_path=output)

    table = pd.read_csv(output, sep="\t", keep_default_na=False)
    assert list(table.columns) == [
        "barcode",
        "true_label",
        "predicted_label",
        "confidence",
        "cluster_id",
        "true_path",
        "predicted_path",
        "true_cl_id",
        "predicted_cl_id",
    ]
    assert table["barcode"].tolist() == ["c0", "c1", "c2", "c3", "c4"]
    assert table["predicted_label"].tolist() == [
        "T cell",
        "B cell",
        "T cell",
        "Unknown",
        "T cell",
    ]
    assert table["predicted_cl_id"].tolist() == [
        "CL:0000084",
        "CL:0000236",
        "CL:0000084",
        "",
        "CL:0000084",
    ]
    assert table["true_cl_id"].tolist() == [
        "CL:0000084",
        "CL:0000236",
        "",
        "CL:0000084",
        "CL:0000084",
    ]
    # The abstained cell sits in a cluster that was decided against a named
    # program, and it reports neither that program's accession nor a path: what
    # the cell received was `Unknown`, and anything else would put a term into
    # the ontology metrics that the system never asserted.
    assert table["predicted_path"].tolist() == [
        "T cell",
        "B cell",
        "T cell",
        "",
        "T cell",
    ]

    assert receipt["n_cells"] == 5
    assert receipt["n_predicted_unknown"] == 1
    assert receipt["coverage"] == pytest.approx(0.8)
    assert receipt["truth_sha256"] == sha256(truth)
    assert receipt["mapping_sha256"]
    assert receipt["label_tables_sha256"] == {"L1.tsv": receipt["label_tables_sha256"]["L1.tsv"]}
    assert receipt["prediction_sha256"] == sha256(output)
    assert receipt["accuracy_claimed"] is False
    assert receipt["hierarchy_paths_emitted"] is True
    assert receipt["granularity"] == "per-cluster, carried to each cell of the cluster"


def test_confidence_is_a_cluster_property_and_says_so(tmp_path: Path) -> None:
    """Annotation is per cluster, so a per-cell confidence column would overclaim.

    Both cells of cluster 0 share one confidence, and the receipt names the
    granularity so a reader comparing against a system that varies within a
    cluster knows the two are not the same quantity.
    """

    config = _run_artifacts(
        tmp_path,
        mapping=_mapping(["c0", "c1"], ["T cell", "T cell"]),
        labels=_labels([("0", "T cell", "CL:0000084", 0.42)]),
    )
    truth = _truth(["c0", "c1"], ["T cell", "T cell"], __path=tmp_path / "truth.tsv")
    output = tmp_path / "predictions.tsv"

    receipt = export_benchmark_predictions(config, truth_path=truth, output_path=output)

    table = pd.read_csv(output, sep="\t", keep_default_na=False)
    assert table["confidence"].tolist() == [0.42, 0.42]
    assert table["cluster_id"].tolist() == [0, 0]
    assert "per-cluster" in receipt["granularity"]


def test_truth_must_cover_exactly_the_run_in_either_direction(tmp_path: Path) -> None:
    """The intersection of two barcode sets is not a decision anyone made.

    A benchmark scoring whatever fraction happened to overlap reports a number
    with no stated denominator, and the two ways to be wrong -- truth covering
    fewer cells than the run, or more -- are indistinguishable from the metric
    side.
    """

    config = _run_artifacts(
        tmp_path,
        mapping=_mapping(["c0", "c1", "c2"], ["T cell"] * 3),
        labels=_labels([("0", "T cell", "CL:0000084", 0.9)]),
    )

    short = _truth(["c0", "c1"], ["T cell", "T cell"], __path=tmp_path / "short.tsv")
    with pytest.raises(ContractError, match="1 run cells have no truth row"):
        export_benchmark_predictions(config, truth_path=short, output_path=tmp_path / "a.tsv")

    long = _truth(
        ["c0", "c1", "c2", "c3"],
        ["T cell"] * 4,
        __path=tmp_path / "long.tsv",
    )
    with pytest.raises(ContractError, match="1 truth rows name a cell the run did not label"):
        export_benchmark_predictions(config, truth_path=long, output_path=tmp_path / "b.tsv")


def test_a_repeated_truth_barcode_is_refused(tmp_path: Path) -> None:
    """The harness scores one row per barcode, so duplicates would be averaged in."""

    config = _run_artifacts(
        tmp_path,
        mapping=_mapping(["c0", "c1"], ["T cell"] * 2),
        labels=_labels([("0", "T cell", "CL:0000084", 0.9)]),
    )
    truth = _truth(["c0", "c0", "c1"], ["T cell"] * 3, __path=tmp_path / "dup.tsv")
    with pytest.raises(ContractError, match="one row per barcode"):
        export_benchmark_predictions(config, truth_path=truth, output_path=tmp_path / "a.tsv")


def test_a_label_with_two_cell_ontology_ids_is_refused(tmp_path: Path) -> None:
    """Row order would silently decide which identifier the metric inherits."""

    first = _labels([("0", "T cell", "CL:0000084", 0.9)])
    # One label, two clusters, two accessions: each cluster decided its own, and
    # row order would otherwise choose which id reaches the published metric.
    conflict = pd.concat(
        [first, _labels([("1", "T cell", "CL:0000625", 0.9)])], ignore_index=True
    )
    config = _run_artifacts(
        tmp_path, mapping=_mapping(["c0"], ["T cell"]), labels=conflict
    )
    truth = _truth(["c0"], ["T cell"], __path=tmp_path / "truth.tsv")
    with pytest.raises(ContractError, match="conflicting Cell Ontology ids"):
        export_benchmark_predictions(config, truth_path=truth, output_path=tmp_path / "a.tsv")


def test_confidence_outside_the_unit_interval_is_refused(tmp_path: Path) -> None:
    """The calibration metric bounds confidence, so the export must too.

    Catching it here names the offending label, where the alternative is a
    published result failing later with no indication which run produced it.
    """

    config = _run_artifacts(
        tmp_path,
        mapping=_mapping(["c0"], ["T cell"]),
        labels=_labels([("0", "T cell", "CL:0000084", 1.4)]),
    )
    truth = _truth(["c0"], ["T cell"], __path=tmp_path / "truth.tsv")
    with pytest.raises(ContractError, match=r"confidence 1.4 outside \[0, 1\]"):
        export_benchmark_predictions(config, truth_path=truth, output_path=tmp_path / "a.tsv")


def test_a_cluster_described_at_two_levels_is_refused_rather_than_guessed(
    tmp_path: Path,
) -> None:
    """The two artifacts share the cluster, not the scope, and only that is checked.

    A recursive run writes one labels file per level, and a cluster id can appear
    at more than one level. `mapping.scope` is the input's scope column and
    `parent_scope` is the vocabulary's parent, and on a configured run they are
    different things, so neither disambiguates. Attaching the wrong level's
    confidence to a cell would be silent.
    """

    config = _run_artifacts(
        tmp_path,
        mapping=_mapping(["c0"], ["T cell"]),
        labels={
            # Distinct labels, so the cluster-id collision is the only ambiguity
            # present. With the same label the accession check fires first, which
            # is the more specific diagnosis and is pinned separately.
            "L1.tsv": _labels([("0", "T cell", "CL:0000084", 0.9)]),
            "L2.tsv": _labels([("0", "naive T cell", "CL:0000895", 0.3)], level="L2"),
        },
    )
    truth = _truth(["c0"], ["T cell"], __path=tmp_path / "truth.tsv")
    with pytest.raises(ContractError, match="more than one level or parent scope"):
        export_benchmark_predictions(
            config, truth_path=truth, output_path=tmp_path / "a.tsv"
        )


def test_a_recursive_run_gets_no_path_rather_than_a_wrong_one(tmp_path: Path) -> None:
    """`hierarchy_error` is the one metric that notices a hierarchy at all.

    The durable mapping records a cell's leaf and the scope its parent was
    decided in, but not the ancestor chain, so a multi-level run cannot supply
    a path without re-walking the checkpoints. Inventing one would corrupt the
    metric silently, so the column is dropped and the receipt records why.
    """

    config = _run_artifacts(
        tmp_path,
        mapping=_mapping(["c0", "c1"], ["T cell", "T cell"]),
        labels={
            "L1.tsv": _labels([("0", "T cell", "CL:0000084", 0.9)]),
            "L2.tsv": _labels(
                [("1", "CD4-positive T cell", "CL:0000624", 0.8)], level="L2"
            ),
        },
    )
    truth = _truth(["c0", "c1"], ["T cell"] * 2, __path=tmp_path / "truth.tsv")
    output = tmp_path / "predictions.tsv"

    receipt = export_benchmark_predictions(config, truth_path=truth, output_path=output)

    table = pd.read_csv(output, sep="\t", keep_default_na=False)
    assert "predicted_path" not in table.columns
    assert receipt["hierarchy_paths_emitted"] is False


def test_the_receipt_is_immutable_against_its_own_inputs(tmp_path: Path) -> None:
    """Re-exporting identical inputs must reproduce; changing truth must be visible.

    `publish_json` refuses a non-identical rewrite, so an export that varies
    with the same run and the same truth fails instead of quietly replacing an
    earlier claim. A changed truth is a different evaluation and belongs at a
    different output path.
    """

    config = _run_artifacts(
        tmp_path,
        mapping=_mapping(["c0"], ["T cell"]),
        labels=_labels([("0", "T cell", "CL:0000084", 0.9)]),
    )
    truth = _truth(["c0"], ["T cell"], __path=tmp_path / "truth.tsv")
    output = tmp_path / "predictions.tsv"
    first = export_benchmark_predictions(config, truth_path=truth, output_path=output)
    again = export_benchmark_predictions(config, truth_path=truth, output_path=output)
    assert again == first

    other_truth = _truth(["c0"], ["B cell"], __path=tmp_path / "other.tsv")
    with pytest.raises(ContractError, match="immutable"):
        export_benchmark_predictions(config, truth_path=other_truth, output_path=output)


def test_a_run_without_a_mapping_or_labels_is_refused_before_any_write(
    tmp_path: Path,
) -> None:
    """A missing run artifact is named, not turned into an empty table."""

    value = yaml.safe_load((SKILL_ROOT / "assets" / "config.template.yaml").read_text())
    value["run"].update(
        {
            "run_id": "absent-v1",
            "output_root": str(tmp_path / "results"),
            "mode": "evidence-only",
        }
    )
    value["input"]["path"] = str(tmp_path / "source.h5ad")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))
    config = load_config(config_path)
    truth = _truth(["c0"], ["T cell"], __path=tmp_path / "truth.tsv")
    output = tmp_path / "predictions.tsv"

    with pytest.raises(ContractError, match="no durable mapping to export"):
        export_benchmark_predictions(config, truth_path=truth, output_path=output)
    assert not output.exists()


def test_the_exported_table_of_a_real_run_scores_end_to_end(tmp_path: Path) -> None:
    """The join has to work against artifacts a pipeline actually produced.

    Every other test here builds its run directory by hand, which would pass
    even if the real artifact names had moved. This one runs the pipeline, so it
    fails if the mapping or the label tables stop being where the exporter
    expects, and it then scores the export through the same public entry point
    an operator would use.
    """

    from cell_curator import evaluate_benchmark_predictions

    # The exported table makes an ontology claim whenever the labels carry CL
    # identifiers, and the harness refuses to score such a table without a
    # snapshot. That refusal is the point: an ontology metric that quietly
    # abstains is the failure mode this work removed.
    ontology = tmp_path / "cl.tsv"
    pd.DataFrame(
        [
            {"id": "CL:0000000", "name": "cell", "parents": ""},
            {"id": "CL:0000003", "name": "native cell", "parents": "CL:0000000"},
            {"id": "CL:0000738", "name": "leukocyte", "parents": "CL:0000003"},
            {"id": "CL:0000542", "name": "lymphocyte", "parents": "CL:0000738"},
            {"id": "CL:0000084", "name": "T cell", "parents": "CL:0000542"},
            {"id": "CL:0000236", "name": "B cell", "parents": "CL:0000542"},
            {"id": "CL:0000624", "name": "CD4-positive T cell", "parents": "CL:0000084"},
        ]
    ).to_csv(ontology, sep="\t", index=False, lineterminator="\n")

    source = tmp_path / "source.h5ad"
    rng = np.random.default_rng(11)
    counts = np.zeros((40, 4), dtype=np.float32)
    counts[:20, :2] = rng.poisson(9, size=(20, 2)) + 1
    counts[20:, 2:4] = rng.poisson(9, size=(20, 2)) + 1
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=pd.DataFrame(
            {"canonical_cluster": ["0"] * 20 + ["1"] * 20, "lineage": ["immune"] * 40},
            index=[f"cell-{index:03d}" for index in range(40)],
        ),
        var=pd.DataFrame({"gene_symbol": ["A1", "A2", "B1", "B2"]}, index=["A1", "A2", "B1", "B2"]),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = np.vstack(
        [rng.normal(-3, 0.2, size=(20, 3)), rng.normal(3, 0.2, size=(20, 3))]
    ).astype(np.float32)
    adata.write_h5ad(source)

    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            {
                "T cell": {"pos": ["A1", "A2"], "neg": ["B1", "B2"]},
                "B cell": {"pos": ["B1", "B2"], "neg": ["A1", "A2"]},
            }
        )
    )
    value = yaml.safe_load((SKILL_ROOT / "assets" / "config.template.yaml").read_text())
    value["run"].update(
        {
            "run_id": "export-real-v1",
            "output_root": str(tmp_path / "results"),
            "mode": "annotation",
        }
    )
    value["input"]["path"] = str(source)
    value["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    value["adaptive"]["eligibility"].update(
        {"min_cells": 5, "min_donors": None, "min_captures": None}
    )
    value["adaptive"]["parent_local"]["enabled"] = False
    value["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": str(markers), "version": "export-real-v1"}
    ]
    value["guidance"] = {
        "interactive": False,
        "preauthorized": True,
        "reviewer": "export-fixture-reviewer",
        "pause_after_each_level": False,
        "assumptions": ["Synthetic markers and context are pre-authorized."],
        "hierarchy_levels": ["L1"],
        "context": {
            "organism": "human",
            "tissue": "synthetic tissue",
            "developmental_stage": "adult",
            "condition": "healthy",
            "technology": "dissociated",
            "experimental_context": "export integration fixture",
            "composition_covariates": [],
            "vocabulary_contract": "local synthetic marker programs",
            "reference_contract": "local_only",
            "visualization_key": "X_pca_harmony",
            "writeback_prefix": "cell_curator",
        },
    }
    value["evidence"].update(
        {
            "min_positive_coverage": 1.0,
            "max_negative_violation": 0.0,
            "min_confidence": 0.6,
            "min_margin": 0.2,
        }
    )
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))
    run_annotation(config_path)
    config = load_config(config_path)

    root = run_root(config)
    mapping = pd.read_csv(root / "mapping" / "parent_child_mapping.tsv", sep="\t")
    truth_path = tmp_path / "truth.tsv"
    pd.DataFrame(
        {
            "barcode": mapping["barcode"].astype(str),
            "true_label": ["T cell"] * 20 + ["B cell"] * 20,
            "true_cl_id": ["CL:0000084"] * 20 + ["CL:0000236"] * 20,
        }
    ).to_csv(truth_path, sep="\t", index=False, lineterminator="\n")

    predictions = tmp_path / "predictions.tsv"
    receipt = export_benchmark_predictions(config, truth_path=truth_path, output_path=predictions)
    assert receipt["n_cells"] == 40

    result = evaluate_benchmark_predictions(
        predictions,
        system="cell-curator-fixture",
        dataset="export-integration-v1",
        output_path=tmp_path / "result.json",
        ontology=ontology,
    )
    assert result.n == 40
    assert result.coverage == 1.0
    assert result.brier_score is not None

    # The markers declare no Cell Ontology identifiers, so the run has no
    # validated accessions to carry and the export records none. The ontology
    # metrics are then absent rather than zero, and absent rather than borrowed
    # from the slash-path metric, which is the behaviour this work established.
    # A table that *does* make an ontology claim is scored by
    # test_benchmark_ontology_distance.py and test_benchmark_relaxed_agreement.py.
    assert result.ontology_distance_error is None
    assert result.relaxed_accuracy is None

    exported = pd.read_csv(predictions, sep="\t", keep_default_na=False, dtype=str)
    assert "None" not in set(exported["predicted_cl_id"])
    assert set(exported["predicted_cl_id"]) == {""}

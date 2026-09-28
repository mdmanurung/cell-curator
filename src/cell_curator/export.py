"""Turn a completed run into a neutral benchmark prediction table.

`benchmark.py` is deliberately blind to the annotation engine so that any system
can be scored on equal footing, which is also why it cannot read a run: there is
nothing in it that knows what a run directory looks like. That leaves the join
unowned, and the only table in `benchmarks/fixtures/` was written by hand. A
hand-assembled table is not an artifact -- it cannot be re-derived, its barcode
set cannot be audited, and a cell silently dropped between the run and the
score is invisible at exactly the point where it would matter. This module owns
that join so the table is a function of two inputs that are both hashed.

The contract is fail-closed on the join. The truth table must name exactly the
barcodes the run labelled: not a subset, not a superset. A benchmark that scores
whatever fraction happened to intersect reports a number nobody chose, and the
two ways it can be wrong -- truth covering fewer cells than the run, or more --
look identical from the metric side. Both are refused here with the counts, and
subsampling is done by re-running on a subsetted input rather than by trimming
the truth after the fact.

Annotation is cluster-level, so per-cell evidence is a lookup, not a
measurement: every cell in a cluster inherits that cluster's label, confidence,
and Cell Ontology identifier. That is a real property of the engine and the
table is honest about it by carrying the cluster key alongside, so a reader can
see that a per-cell confidence column is really a per-cluster one. A system
that varies its call within a cluster is being measured on a different quantity
and must not be compared to this one without saying so.

Label hierarchies are not reconstructed here. The durable mapping records each
cell's leaf and the scope its parent was decided in, but not the chain of
ancestors between them, and a single-level run therefore yields a real
`predicted_path` while a recursive run does not. Inventing a chain would
silently corrupt `hierarchy_error`, which is the one metric that notices a
hierarchy at all, so the path columns are emitted only when they are derivable
and the receipt records which case applied.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .config import CellCuratorConfig
from .contracts import ContractError, require_columns, sha256
from .provenance import atomic_publish, publish_json

# The truth side is a third-party artifact, so its contract is declared here
# rather than inherited: a table that lacks these cannot be scored, and one that
# carries extra columns is fine because a curator's export usually has many.
TRUTH_REQUIRED_COLUMNS = {"barcode", "true_label"}
TRUTH_OPTIONAL_COLUMNS = ("true_cl_id", "true_path")

# Deliberately narrow. The neutral format accepts many more optional columns, and
# emitting a column whose value is invented would corrupt the metric it feeds.
EMITTED_COLUMNS = (
    "barcode",
    "true_label",
    "predicted_label",
    "confidence",
    "cluster_id",
    "true_path",
    "predicted_path",
    "true_cl_id",
    "predicted_cl_id",
)


def _label_lookups(
    frames: list[pd.DataFrame],
) -> tuple[dict[str, str], dict[str, float], dict[str, str]]:
    """Index the label tables by label text, refusing any label that disagrees with itself.

    A label reaching two levels or two parent scopes is already a contract
    violation elsewhere, but it would corrupt the export quietly here: one of the
    two Cell Ontology identifiers would win by row order. Refusing is the only
    outcome that names the problem.
    """

    cl_id: dict[str, str] = {}
    confidence: dict[str, float] = {}
    path: dict[str, str] = {}
    for frame in frames:
        for row in frame.to_dict("records"):
            label = str(row.get("label", ""))
            if not label:
                continue
            accession = str(row.get("cl_id", "") or "")
            previous = cl_id.get(label)
            if previous is not None and previous != accession:
                raise ContractError(
                    f"label {label!r} carries conflicting Cell Ontology ids: "
                    f"{previous!r} and {accession!r}"
                )
            cl_id[label] = accession
            raw_confidence = str(row.get("confidence", "") or "").strip()
            if raw_confidence:
                value = float(raw_confidence)
                if not 0.0 <= value <= 1.0:
                    raise ContractError(
                        f"label {label!r} carries confidence {value} outside [0, 1]"
                    )
                recorded = confidence.get(label)
                if recorded is not None and abs(recorded - value) > 1e-9:
                    raise ContractError(
                        f"label {label!r} carries conflicting confidence: {recorded} and {value}"
                    )
                confidence[label] = value
            rendered = "/".join(
                str(row.get(column, "")).strip()
                for column in ("level", "parent_scope", "label")
            )
            if not path.setdefault(label, rendered) == rendered:
                path[label] = ""
    return cl_id, confidence, path


def _load_truth(truth_path: Path) -> pd.DataFrame:
    """Read a third-party truth table under this module's declared contract."""

    if not truth_path.is_file():
        raise ContractError(f"truth table is absent: {truth_path}")
    separator = "\t" if truth_path.suffix.lower() in {".tsv", ".txt"} else ","
    frame = pd.read_csv(truth_path, sep=separator, keep_default_na=False, dtype=str)
    require_columns(frame, TRUTH_REQUIRED_COLUMNS, artifact=f"truth table {truth_path}")
    frame = frame.fillna("")
    frame["barcode"] = frame["barcode"].astype(str)
    frame["true_label"] = frame["true_label"].astype(str)
    if frame["barcode"].duplicated().any():
        duplicated = sorted(frame.loc[frame["barcode"].duplicated(), "barcode"].unique())
        raise ContractError(
            f"truth table must contain one row per barcode; repeated: {duplicated[:5]}"
        )
    if (frame["true_label"].str.strip() == "").any():
        raise ContractError("truth table contains a row with a blank true_label")
    for column in TRUTH_OPTIONAL_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
        frame[column] = frame[column].astype(str)
    return frame


def export_benchmark_predictions(
    config: CellCuratorConfig,
    *,
    truth_path: str | Path,
    output_path: str | Path,
) -> dict[str, Any]:
    """Write a run's predictions as a neutral-format table plus a binding receipt.

    The receipt is immutable, as is the table beside it: re-exporting the same
    run against the same truth reproduces both byte for byte, and an export that
    would change either is refused rather than replacing an earlier claim. A
    genuine change has to move a hashed input rather than edit the output.
    """

    truth_file = Path(truth_path)
    destination = Path(output_path)
    root = Path(config.run.output_root) / config.run.run_id

    mapping_path = root / "mapping" / "parent_child_mapping.tsv"
    if not mapping_path.is_file():
        raise ContractError(
            f"run has no durable mapping to export: {mapping_path}; complete a run first"
        )
    mapping = pd.read_csv(mapping_path, sep="\t", keep_default_na=False, dtype=str)
    require_columns(
        mapping,
        {"barcode", "cell_type", "cell_state", "canonical_parent", "scope"},
        artifact="run mapping",
    )
    mapping["barcode"] = mapping["barcode"].astype(str)
    if mapping["barcode"].duplicated().any():
        raise ContractError("run mapping must contain one row per barcode")

    label_paths = sorted((root / "labels").glob("*.tsv"))
    if not label_paths:
        raise ContractError(f"run has no durable label tables to export: {root / 'labels'}")
    frames = [pd.read_csv(path, sep="\t", keep_default_na=False) for path in label_paths]
    cl_lookup, confidence_lookup, path_lookup = _label_lookups(frames)

    single_level = len({str(frame["level"].iloc[0]) for frame in frames if len(frame)}) == 1
    truth = _load_truth(truth_file)

    labelled = set(mapping["barcode"])
    annotated = set(truth["barcode"])
    missing_truth = sorted(labelled - annotated)
    unknown_truth = sorted(annotated - labelled)
    if missing_truth or unknown_truth:
        raise ContractError(
            "truth table and run must cover exactly the same barcodes; "
            f"{len(missing_truth)} run cells have no truth row, "
            f"{len(unknown_truth)} truth rows name a cell the run did not label. "
            f"Run cells without truth: {missing_truth[:5]}; "
            f"unknown truth barcodes: {unknown_truth[:5]}. "
            "Subsample the input and re-run rather than trimming the truth after "
            "the fact, or the scored cell set is one nobody chose."
        )

    joined = truth.merge(
        mapping[["barcode", "cell_type", "cell_state", "canonical_parent", "scope"]],
        on="barcode",
        how="inner",
        validate="one_to_one",
    )
    joined["predicted_label"] = joined["cell_type"].astype(str)
    joined["predicted_cl_id"] = joined["predicted_label"].map(cl_lookup).fillna("")
    joined["confidence"] = [
        confidence_lookup.get(label, "") for label in joined["predicted_label"]
    ]
    joined["cluster_id"] = joined["canonical_parent"].astype(str)
    for column in ("true_cl_id", "true_path"):
        joined[column] = joined[column].astype(str)
    if single_level:
        joined["predicted_path"] = joined["predicted_label"].map(path_lookup).fillna("")
    else:
        joined["predicted_path"] = ""

    emitted = [column for column in EMITTED_COLUMNS if column in joined.columns]
    if not joined["predicted_path"].astype(str).str.strip().any():
        emitted = [column for column in emitted if column != "predicted_path"]
        for column in ("true_path",):
            if column in emitted and not joined[column].astype(str).str.strip().any():
                emitted = [item for item in emitted if item != column]
    table = joined[emitted].sort_values("barcode", ignore_index=True)

    confidence_values = pd.to_numeric(table["confidence"], errors="coerce").dropna()
    if len(confidence_values) and (
        (confidence_values < 0.0).any() or (confidence_values > 1.0).any()
    ):
        raise ContractError("exported confidence must lie within [0, 1]")

    atomic_publish(
        destination,
        table.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    abstained = int(table["predicted_label"].astype(str).eq("Unknown").sum())
    receipt = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_prediction_export",
        "run_id": config.run.run_id,
        "run_root": str(root),
        "n_cells": int(len(table)),
        "n_predicted_unknown": abstained,
        "coverage": float(1.0 - abstained / len(table)) if len(table) else None,
        "columns": emitted,
        "truth_path": str(truth_file.resolve()),
        "truth_sha256": sha256(truth_file),
        "mapping_sha256": sha256(mapping_path),
        "label_tables_sha256": {path.name: sha256(path) for path in label_paths},
        "prediction_sha256": sha256(destination),
        "hierarchy_paths_emitted": bool(single_level),
        "granularity": "per-cluster, carried to each cell of the cluster",
        "neutral_format": True,
        "accuracy_claimed": False,
    }
    publish_json(destination.with_suffix(destination.suffix + ".receipt.json"), receipt)
    return receipt


__all__ = ["export_benchmark_predictions", "EMITTED_COLUMNS", "TRUTH_REQUIRED_COLUMNS"]

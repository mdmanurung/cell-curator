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
from .contracts import ContractError, identifier_text, measured_float, require_columns, sha256
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


def _cluster_rows(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Index the label tables by the parent scope and cluster they describe.

    Confidence and the Cell Ontology identifier are properties of a *cluster's*
    decision, not of a label's text. A label reaches many clusters, and a run
    whose clusters mostly abstain produces a different confidence for every one
    of them under the same text, so looking these up by label conflates distinct
    measurements. The first benchmark run refused to export for exactly that
    reason, having compared two clusters' confidences under one `Unknown`.

    Joining on the cluster is also the only join that can detect a label
    appearing at two levels or in two parent scopes: a cluster that the mapping
    places in a scope the label table does not describe is a hole in the
    vocabulary, and refusing names it instead of letting a neighbouring cluster's
    value fill the gap.
    """

    rows: list[dict[str, object]] = []
    for frame in frames:
        for row in frame.to_dict("records"):
            confidence = measured_float(row.get("confidence"))
            if confidence is not None and not 0.0 <= confidence <= 1.0:
                raise ContractError(
                    f"cluster {row.get('cluster_id')!r} carries confidence {confidence} "
                    "outside [0, 1]"
                )
            parent_scope = str(row.get("parent_scope", ""))
            cluster_id = str(row.get("cluster_id", ""))
            rows.append(
                {
                    "parent_scope": parent_scope,
                    "cluster_id": cluster_id,
                    "level": str(row.get("level", "")),
                    "label": str(row.get("label", "")),
                    "cl_id": identifier_text(row.get("cl_id")),
                    "confidence": confidence,
                    # Biological ancestry only: level and ROOT tokens never match
                    # a truth path, which made every correct call a full error.
                    "path": "/".join(
                        item
                        for item in (parent_scope.strip(), str(row.get("label", "")).strip())
                        if item and item != "ROOT"
                    ),
                }
            )
    if not rows:
        raise ContractError("run label tables contain no rows to index")
    table = pd.DataFrame(rows)
    # A label is one biological identity, so one accession. Two clusters reaching
    # the same label with different accessions means one of them is wrong, and
    # picking between them by row order would put an unearned id in a published
    # metric.
    for label, block in table[table["cl_id"].astype(str).ne("")].groupby("label"):
        accessions = sorted(set(block["cl_id"].astype(str)))
        if len(accessions) > 1:
            raise ContractError(
                f"label {label!r} carries conflicting Cell Ontology ids across clusters: "
                f"{accessions}"
            )
    return table


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
    clusters = _cluster_rows(frames)
    from .review import resolve_durable_mapping_values

    # Score what write-back writes: the reviewed labels resolved down each parent
    # path, not the pre-review mapping column, which ignored reviewer edits and
    # accepted children.
    resolved = resolve_durable_mapping_values(mapping, pd.concat(frames, ignore_index=True))
    mapping["cell_type"] = resolved["cell_type"]
    mapping["cell_state"] = resolved["cell_state"]

    single_level = len({str(frame["level"].iloc[0]) for frame in frames if len(frame)}) == 1
    # `mapping.scope` is the input's scope column and `parent_scope` is the
    # vocabulary's parent, and on a configured run those are different things --
    # the AIDA run set its scope column to "ROOT" so they happened to agree, and a
    # run that stratifies by lineage does not. `canonical_parent` is the
    # cluster, and that is the only column the two artifacts genuinely share.
    described = set(clusters["cluster_id"].astype(str))
    duplicated = sorted(
        clusters.loc[clusters["cluster_id"].astype(str).duplicated(keep=False), "cluster_id"]
        .astype(str)
        .unique()
    )
    if duplicated:
        raise ContractError(
            "the same cluster is described at more than one level or parent scope, so "
            f"its evidence cannot be attached to it unambiguously: {duplicated[:5]}. "
            "Export one level at a time, or use a run whose cluster ids are unique."
        )
    undescribed = sorted(
        set(mapping["canonical_parent"].astype(str)).difference(described)
    )
    if undescribed:
        raise ContractError(
            "the run mapping places cells in parents the label tables do not describe, "
            f"so no evidence can be attached to them: {sorted(set(undescribed))[:5]}. "
            "This means the vocabulary and the executed scopes disagree."
        )
    clusters = clusters.drop_duplicates(subset=["cluster_id"], keep="last")
    mapping = mapping.merge(
        clusters[["cluster_id", "cl_id", "confidence", "path"]].rename(
            columns={
                "cluster_id": "canonical_parent",
                "cl_id": "_cl_id",
                "confidence": "_confidence",
                "path": "_path",
            }
        ),
        on="canonical_parent",
        how="left",
        validate="many_to_one",
    )
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
        mapping[
            [
                "barcode",
                "cell_type",
                "cell_state",
                "canonical_parent",
                "scope",
                "_cl_id",
                "_confidence",
                "_path",
            ]
        ],
        on="barcode",
        how="inner",
        validate="one_to_one",
    )
    joined["predicted_label"] = joined["cell_type"].astype(str)
    # The accession, confidence, and path are the canonical cluster's. A cell whose
    # label came from an accepted child carries none of them rather than its parent's.
    parent_label = clusters.set_index(clusters["cluster_id"].astype(str))["label"].astype(str)
    refined = joined["predicted_label"].ne(
        joined["canonical_parent"].astype(str).map(parent_label).fillna("")
    )
    joined.loc[refined, ["_cl_id", "_path"]] = ""
    joined.loc[refined, "_confidence"] = None
    # An abstention has no term, so it carries no accession even when its cluster
    # was decided against a named program. The cluster's own decision is kept in
    # the label table; what the cell received is `Unknown`, and reporting an
    # accession for it would put a term in the ontology metrics that the system
    # never asserted. This is the same rule the label builder uses.
    abstained = joined["predicted_label"].eq("Unknown")
    joined["predicted_cl_id"] = [
        "" if is_abstention else str(accession)
        for is_abstention, accession in zip(abstained, joined["_cl_id"], strict=True)
    ]
    joined["confidence"] = [
        "" if value is None or (isinstance(value, float) and pd.isna(value)) else value
        for value in joined["_confidence"]
    ]
    joined["cluster_id"] = joined["canonical_parent"].astype(str)
    for column in ("true_cl_id", "true_path"):
        joined[column] = joined[column].astype(str)
    if single_level:
        joined["predicted_path"] = [
            ""
            if is_abstention
            or not str(value).strip()
            or str(value).endswith("/")
            else str(value)
            for is_abstention, value in zip(abstained, joined["_path"], strict=True)
        ]
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

"""Approval-gated dry-run and atomic AnnData/MuData annotation write-back."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import CellCuratorConfig
from .contracts import ContractError, sha256
from .hierarchy import validate_hierarchy_nesting
from .ontology import (
    CellOntology,
    ontology_from_config,
    validate_durable_label_ontology,
)
from .provenance import atomic_replace
from .review import (
    check_report_completeness,
    hierarchy_level_sort_key,
    resolve_durable_mapping_values,
    validate_reconciliation,
)


def _load_durable_labels(run_root: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    paths = sorted((run_root / "labels").glob("*.tsv"))
    if not paths:
        raise ContractError("write-back requires durable labels tables")
    labels = pd.concat(
        [pd.read_csv(path, sep="\t", keep_default_na=False) for path in paths],
        ignore_index=True,
    )
    required = {
        "level",
        "parent_scope",
        "cluster_id",
        "label",
        "cell_state",
        "cl_id",
    }
    missing = required.difference(labels.columns)
    if missing:
        raise ContractError(f"durable labels miss write-back columns: {sorted(missing)}")
    if labels.duplicated(["level", "parent_scope", "cluster_id"]).any():
        raise ContractError("durable labels duplicate a level/parent/cluster key")
    return labels, {str(path.relative_to(run_root)): sha256(path) for path in paths}


def _validate_durable_parent_paths(labels: pd.DataFrame) -> None:
    available_parents = {"ROOT"}
    for index, level in enumerate(
        sorted(labels["level"].astype(str).unique(), key=hierarchy_level_sort_key)
    ):
        block = labels[labels["level"].astype(str).eq(level)]
        parents = set(block["parent_scope"].astype(str))
        if index == 0 and parents != {"ROOT"}:
            raise ContractError("the first durable hierarchy level must be ROOT-scoped")
        missing = parents.difference(available_parents)
        if missing:
            raise ContractError(
                f"durable hierarchy level {level} has unresolved parent scopes: {sorted(missing)}"
            )
        available_parents.update(
            block.loc[~block["label"].astype(str).eq("Unknown"), "label"].astype(str)
        )


def _category_color_plan(
    source_object: Any,
    column: str,
    values: np.ndarray,
    report_contract: dict[str, Any],
) -> dict[str, list[str]]:
    observed = list(dict.fromkeys(map(str, values)))
    report_categories = list(map(str, report_contract.get("categories", [])))
    report_colors = list(map(str, report_contract.get("colors", [])))
    if report_colors and len(report_colors) != len(report_categories):
        raise ContractError(f"report color contract is misaligned for {column}")
    if len(set(report_categories)) != len(report_categories):
        raise ContractError(f"report category contract contains duplicates for {column}")

    existing_categories: list[str] = []
    if column in source_object.obs and isinstance(
        source_object.obs[column].dtype, pd.CategoricalDtype
    ):
        existing_categories = list(map(str, source_object.obs[column].cat.categories))
    categories = list(report_categories or existing_categories)
    categories.extend(item for item in observed if item not in categories)

    color_lookup: dict[str, str] = {}
    if report_colors:
        color_lookup.update(zip(report_categories, report_colors, strict=True))
    elif existing_categories:
        existing_colors = source_object.uns.get(f"{column}_colors")
        if existing_colors is not None and len(existing_colors) == len(existing_categories):
            color_lookup.update(
                zip(existing_categories, map(str, existing_colors), strict=True)
            )
    colors = [color_lookup.get(category, "#BBBBBB") for category in categories]
    if not color_lookup:
        colors = []
    return {"categories": categories, "colors": colors}


def _align_mapping(mapping: pd.DataFrame, obs_names: pd.Index) -> pd.DataFrame:
    indexed = mapping.set_index("barcode", verify_integrity=True)
    aligned = indexed.loc[list(map(str, obs_names))].copy()
    aligned.insert(0, "barcode", aligned.index.astype(str))
    return aligned.reset_index(drop=True)


def _approval(
    config: CellCuratorConfig,
    run_root: Path,
    approval_path: Path,
    preview: dict[str, Any],
) -> dict[str, Any]:
    if not approval_path.is_file():
        raise ContractError("write-back requires an explicit approval JSON file")
    value = json.loads(approval_path.read_text())
    required_flags = ("approved", "assumptions_approved", "critics_reconciled")
    missing_flags = [key for key in required_flags if value.get(key) is not True]
    if missing_flags or not str(value.get("reviewer", "")).strip():
        raise ContractError(
            "write-back approval must name a reviewer and set approved, "
            f"assumptions_approved, and critics_reconciled to true; missing {missing_flags}"
        )
    expected = {
        "run_id": config.run.run_id,
        "input_sha256": preview["input_sha256"],
        "mapping_sha256": preview["mapping_sha256"],
        "labels_sha256": preview["labels_sha256"],
        "diff_sha256": preview["diff_sha256"],
        "assumptions_sha256": preview["assumptions_sha256"],
        "report_manifest_sha256": preview["report_manifest_sha256"],
        "report_sha256": preview["report_sha256"],
        "critic_packet_sha256": preview["critic_packet_sha256"],
        "critic_findings_sha256": preview["critic_findings_sha256"],
        "critic_reconciliation_sha256": preview["critic_reconciliation_sha256"],
        "output_path": preview["output_path"],
        "in_place": preview["in_place"],
    }
    for key, expected_value in expected.items():
        if value.get(key) != expected_value:
            raise ContractError(f"write-back approval {key} does not match the reviewed run")
    return value


def validate_writeback(
    config: CellCuratorConfig,
    run_root: Path,
    *,
    ontology: CellOntology | None = None,
    output_path: Path | None = None,
    in_place: bool = False,
) -> dict[str, Any]:
    ontology = ontology or ontology_from_config(config)
    source_path = Path(config.input.path)
    target = output_path or (
        source_path
        if in_place
        else run_root / "writeback" / f"annotated.{config.input.kind}"
    )
    source_resolved = source_path.resolve()
    target_resolved = target.resolve()
    if in_place != (target_resolved == source_resolved):
        raise ContractError(
            "write-back intent mismatch: --in-place must target the source and a preview "
            "without --in-place must target a distinct output"
        )
    source_sha256 = sha256(source_path)
    mapping_path = run_root / "mapping" / "parent_child_mapping.tsv"
    validation_path = run_root / "mapping" / "mapping.validation.json"
    if not mapping_path.is_file() or not validation_path.is_file():
        raise ContractError("write-back requires the immutable mapping and validation")
    validation = json.loads(validation_path.read_text())
    mapping_sha256 = sha256(mapping_path)
    invariants = validation.get("invariants")
    if (
        validation.get("status") != "complete"
        or validation.get("mapping_sha256") != mapping_sha256
        or not isinstance(invariants, dict)
        or not invariants
        or not all(value is True for value in invariants.values())
    ):
        raise ContractError("mapping validation is incomplete or hash-stale")
    reconciliation = validate_reconciliation(run_root)
    state_path = run_root / "run_state.json"
    if not state_path.is_file():
        raise ContractError("write-back requires the guided run state")
    from .guided import load_state

    state = load_state(run_root)
    if state.input_sha256 != source_sha256:
        raise ContractError("source input changed after the guided run was frozen")
    assumptions_path = run_root / "assumptions.json"
    if (
        not state.assumptions_approved
        or not state.approved_assumptions_sha256
        or not assumptions_path.is_file()
        or sha256(assumptions_path) != state.approved_assumptions_sha256
    ):
        raise ContractError("write-back requires the recorded assumptions gate")
    check_report_completeness(config, run_root, require_report=True)
    mapping = pd.read_csv(mapping_path, sep="\t", keep_default_na=False)
    required_mapping_columns = {"barcode", "leaf_id", "canonical_parent"}
    missing_mapping_columns = required_mapping_columns.difference(mapping.columns)
    if missing_mapping_columns:
        raise ContractError(
            f"write-back mapping misses columns: {sorted(missing_mapping_columns)}"
        )
    labels, label_hashes = _load_durable_labels(run_root)
    _validate_durable_parent_paths(labels)
    ontology_validation = validate_durable_label_ontology(labels, ontology)
    hierarchy_validation = validate_hierarchy_nesting(
        [labels],
        cell_membership=mapping[["barcode", "leaf_id"]],
    )
    if config.run.mode.value != "annotation":
        raise ContractError("evidence-only runs cannot write biological labels")
    context_path = run_root / "context.json"
    prefix = config.guidance.context.writeback_prefix
    if context_path.is_file():
        context = json.loads(context_path.read_text())
        prefix = str(context.get("biological_context", {}).get("writeback_prefix", prefix))
    if not prefix or any(
        character not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_"
        for character in prefix
    ):
        raise ContractError(
            "write-back prefix must contain only letters, digits, and underscore"
        )
    output_columns = [f"{prefix}_type", f"{prefix}_state", f"{prefix}_leaf"]
    colors_path = run_root / "review" / "colors.json"
    report_colors = (
        json.loads(colors_path.read_text()).get("columns", {}) if colors_path.is_file() else {}
    )
    if config.input.kind == "h5ad":
        import anndata as ad

        source_object = ad.read_h5ad(config.input.path, backed="r")
    else:
        import mudata as md

        source_object = md.read_h5mu(config.input.path, backed=True)
    try:
        obs_names = source_object.obs_names.astype(str)
        if set(obs_names) != set(mapping["barcode"].astype(str)):
            raise ContractError("write-back mapping barcodes differ from the source object")
        aligned = _align_mapping(mapping, obs_names)
        durable_values = resolve_durable_mapping_values(aligned, labels)
        new_values = {
            output_columns[0]: durable_values["cell_type"],
            output_columns[1]: durable_values["cell_state"],
            output_columns[2]: durable_values["leaf_id"],
        }
        diff_rows: list[pd.DataFrame] = []
        changed_by_column: dict[str, int] = {}
        metadata_changed_by_column: dict[str, int] = {}
        planned_contract: dict[str, dict[str, list[str]]] = {}
        for column, values in new_values.items():
            old = (
                source_object.obs[column].astype(str).to_numpy()
                if column in source_object.obs
                else np.full(len(obs_names), "<ABSENT>", dtype=object)
            )
            changed = old != values
            changed_by_column[column] = int(np.sum(changed))
            diff_rows.append(
                pd.DataFrame(
                    {
                        "barcode": obs_names,
                        "column": column,
                        "old_value": old,
                        "new_value": values,
                        "changed": changed,
                    }
                )
            )
            contract = _category_color_plan(
                source_object,
                column,
                values,
                report_colors.get(column, {}),
            )
            planned_contract[column] = contract
            if column in source_object.obs and isinstance(
                source_object.obs[column].dtype, pd.CategoricalDtype
            ):
                old_categories = list(map(str, source_object.obs[column].cat.categories))
                old_ordered: bool | str = bool(source_object.obs[column].cat.ordered)
            else:
                old_categories = []
                old_ordered = "<ABSENT>"
            old_color_values = source_object.uns.get(f"{column}_colors")
            old_colors = (
                list(map(str, old_color_values)) if old_color_values is not None else []
            )
            metadata_rows = [
                {
                    "barcode": "<METADATA>",
                    "column": f"{column}.__categories__",
                    "old_value": json.dumps(old_categories),
                    "new_value": json.dumps(contract["categories"]),
                    "changed": old_categories != contract["categories"],
                },
                {
                    "barcode": "<METADATA>",
                    "column": f"{column}.__colors__",
                    "old_value": json.dumps(old_colors),
                    "new_value": json.dumps(contract["colors"]),
                    "changed": old_colors != contract["colors"],
                },
                {
                    "barcode": "<METADATA>",
                    "column": f"{column}.__ordered__",
                    "old_value": json.dumps(old_ordered),
                    "new_value": json.dumps(True),
                    "changed": old_ordered is not True,
                },
            ]
            metadata_changed_by_column[column] = int(
                sum(bool(row["changed"]) for row in metadata_rows)
            )
            diff_rows.append(pd.DataFrame(metadata_rows))
    finally:
        file = getattr(source_object, "file", None)
        if file is not None:
            file.close()
    if sha256(source_path) != source_sha256:
        raise ContractError("read-only write-back validation mutated the source input")
    diff = pd.concat(diff_rows, ignore_index=True)
    diff_path = run_root / "writeback" / "dry_run_diff.tsv"
    atomic_replace(diff_path, diff.to_csv(sep="\t", index=False).encode())
    resolved_types = pd.Series(durable_values["cell_type"], dtype="string")
    resolved_states = pd.Series(durable_values["cell_state"], dtype="string")
    review_paths = {
        "report_manifest_sha256": run_root / "review" / "report.manifest.json",
        "report_sha256": run_root / "review" / "report.html",
        "critic_packet_sha256": run_root / "review" / "critic_packet.json",
        "critic_findings_sha256": run_root / "review" / "critic_findings.tsv",
        "critic_reconciliation_sha256": run_root
        / "review"
        / "critic_reconciliation.tsv",
    }
    missing_review = [str(path) for path in review_paths.values() if not path.is_file()]
    if missing_review:
        raise ContractError(f"write-back review artifacts are absent: {missing_review}")
    return {
        "schema_version": 3,
        "artifact_kind": "cell_curator_writeback_preview",
        "status": "complete",
        "run_id": config.run.run_id,
        "n_cells": len(mapping),
        "n_cell_types": int(resolved_types.nunique()),
        "n_states": int(resolved_states.nunique()),
        "n_unknown": int(resolved_types.eq("Unknown").sum()),
        "input_sha256": source_sha256,
        "mapping_sha256": mapping_sha256,
        "labels_sha256": label_hashes,
        "ontology_validation": ontology_validation,
        "hierarchy_validation": hierarchy_validation,
        "critics_reconciled": reconciliation["all_reconciled"],
        "assumptions_sha256": state.approved_assumptions_sha256,
        **{name: sha256(path) for name, path in review_paths.items()},
        "output_path": str(target_resolved),
        "in_place": in_place,
        "output_columns": output_columns,
        "changed_by_column": changed_by_column,
        "metadata_changed_by_column": metadata_changed_by_column,
        "diff_rows": len(diff),
        "diff_path": str(diff_path),
        "diff_sha256": sha256(diff_path),
        "category_and_color_contract": planned_contract,
    }


def write_back(
    config: CellCuratorConfig,
    run_root: Path,
    *,
    output_path: Path,
    approval_path: Path,
    ontology: CellOntology | None = None,
    in_place: bool = False,
) -> Path:
    preview = validate_writeback(
        config,
        run_root,
        ontology=ontology,
        output_path=output_path,
        in_place=in_place,
    )
    _approval(config, run_root, approval_path, preview)
    source = Path(config.input.path)
    source_sha256 = sha256(source)
    receipt_path = run_root / "writeback" / "receipt.json"
    if receipt_path.exists():
        prior = json.loads(receipt_path.read_text())
        if (
            prior.get("output_path") == str(output_path.resolve())
            and prior.get("in_place") is in_place
            and output_path.is_file()
            and prior.get("output_sha256") == sha256(output_path)
            and prior.get("approval_sha256") == sha256(approval_path)
            and all(prior.get(key) == preview.get(key) for key in preview)
        ):
            return output_path
        raise ContractError(
            "an immutable write-back receipt already exists for a different or drifted output"
        )
    mapping = pd.read_csv(
        run_root / "mapping" / "parent_child_mapping.tsv",
        sep="\t",
        keep_default_na=False,
    )
    labels, _ = _load_durable_labels(run_root)
    output_columns = list(preview["output_columns"])
    color_contract = preview.get("category_and_color_contract", {})
    partial = output_path.with_name(f".{output_path.name}.partial.{os.getpid()}")
    receipt_partial = receipt_path.with_name(
        f".{receipt_path.name}.partial.{os.getpid()}"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and output_path.resolve() != source.resolve():
        raise ContractError(f"refusing to overwrite write-back output: {output_path}")
    try:
        if config.input.kind == "h5ad":
            import anndata as ad

            obj = ad.read_h5ad(source)
            aligned = _align_mapping(mapping, obj.obs_names.astype(str))
            durable_values = resolve_durable_mapping_values(aligned, labels)
            value_columns = ("cell_type", "cell_state", "leaf_id")
            for output_column, value_column in zip(output_columns, value_columns, strict=True):
                values = durable_values[value_column]
                contract = color_contract.get(output_column, {})
                categories = list(map(str, contract.get("categories", [])))
                categories.extend(
                    item for item in dict.fromkeys(values) if item not in categories
                )
                obj.obs[output_column] = pd.Categorical(
                    values, categories=categories, ordered=True
                )
                colors = list(map(str, contract.get("colors", [])))
                if colors:
                    if len(colors) != len(categories):
                        raise ContractError(
                            f"write-back color contract is misaligned for {output_column}"
                        )
                    obj.uns[f"{output_column}_colors"] = colors
                else:
                    obj.uns.pop(f"{output_column}_colors", None)
            obj.write_h5ad(partial)
        else:
            import mudata as md

            obj = md.read_h5mu(source)
            aligned = _align_mapping(mapping, obj.obs_names.astype(str))
            durable_values = resolve_durable_mapping_values(aligned, labels)
            value_columns = ("cell_type", "cell_state", "leaf_id")
            for output_column, value_column in zip(output_columns, value_columns, strict=True):
                values = durable_values[value_column]
                contract = color_contract.get(output_column, {})
                categories = list(map(str, contract.get("categories", [])))
                categories.extend(
                    item for item in dict.fromkeys(values) if item not in categories
                )
                obj.obs[output_column] = pd.Categorical(
                    values, categories=categories, ordered=True
                )
                colors = list(map(str, contract.get("colors", [])))
                if colors:
                    if len(colors) != len(categories):
                        raise ContractError(
                            f"write-back color contract is misaligned for {output_column}"
                        )
                    obj.uns[f"{output_column}_colors"] = colors
                else:
                    obj.uns.pop(f"{output_column}_colors", None)
            obj.write_h5mu(partial)
        receipt = {
            **preview,
            "artifact_kind": "cell_curator_writeback_receipt",
            "output_path": str(output_path.resolve()),
            "output_sha256": sha256(partial),
            "approval_sha256": sha256(approval_path),
        }
        receipt_partial.write_text(
            json.dumps(receipt, indent=2, sort_keys=True, default=str) + "\n"
        )
        os.replace(partial, output_path)
        os.replace(receipt_partial, receipt_path)
    finally:
        partial.unlink(missing_ok=True)
        receipt_partial.unlink(missing_ok=True)
    if not in_place and sha256(source) != source_sha256:
        raise ContractError("write-back changed the immutable source input")
    return output_path

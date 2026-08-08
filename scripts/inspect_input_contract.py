#!/usr/bin/env python3
"""Inspect and freeze a capability-aware AnnData or MuData input contract."""

from __future__ import annotations

import argparse
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import sparse

from contracts import (
    ContractError,
    assay_registry,
    atomic_write_json,
    atomic_write_text,
    evidence_lanes,
    read_config,
    representation_decl,
    sha256,
    values_sha256,
)


def _minimum(matrix) -> float:
    if sparse.issparse(matrix):
        stored = float(matrix.data.min()) if matrix.nnz else 0.0
        return min(stored, 0.0)
    return float(np.asarray(matrix).min())


def _matrix(adata, declaration: dict):
    source = str(declaration["source"]).lower()
    if source == "x":
        return adata.X
    key = str(declaration["key"])
    if key not in adata.layers:
        raise ContractError(f"declared layer is absent: {key}")
    return adata.layers[key]


def _inspect_obs(adata, config: dict, *, is_mudata: bool) -> dict:
    keys = config["input"]["keys"]
    required_names = list(config["input"].get("required_obs_keys", []))
    required_names.extend(["canonical_parent", "scope"])
    required_obs = set()
    for logical_name in required_names:
        key = str(keys.get(logical_name, "")).strip()
        if not key:
            raise ContractError(f"required obs key {logical_name!r} is blank")
        required_obs.add(key)
    required_obs.update(config["input"].get("guide_columns", []))
    required_obs.update(config["input"].get("technical_metrics", {}).values())
    missing = required_obs.difference(adata.obs.columns)
    if missing:
        raise ContractError(f"input observations miss keys: {sorted(missing)}")
    stored = list(keys.get("stored_resolutions", []))
    missing_stored = set(stored).difference(adata.obs.columns)
    if missing_stored:
        raise ContractError(f"stored resolution keys are absent: {sorted(missing_stored)}")
    if not adata.obs_names.is_unique:
        raise ContractError("input barcodes are not unique")
    parent_key = keys["canonical_parent"]
    if adata.obs[parent_key].isna().any():
        raise ContractError("canonical parent key contains missing values")
    return {
        "n_cells": int(adata.n_obs),
        "n_features_root": int(adata.n_vars),
        "observation_columns": sorted(map(str, adata.obs.columns)),
        "canonical_parent_count": int(adata.obs[parent_key].astype(str).nunique()),
        "canonical_parent_counts": {
            str(key): int(value)
            for key, value in adata.obs[parent_key]
            .astype(str)
            .value_counts()
            .sort_index()
            .items()
        },
        "stored_resolution_keys": stored,
        "ordered_barcode_sha256": values_sha256(list(adata.obs_names)),
        "sorted_barcode_sha256": values_sha256(
            list(adata.obs_names), sort_values=True
        ),
        "is_mudata": is_mudata,
    }


def _inspect_assay(config: dict, assay_id: str, adata) -> tuple[dict, list[str]]:
    assay = assay_registry(config)[assay_id]
    feature_id_column = str(assay.get("feature_id_column", "")).strip()
    if feature_id_column and feature_id_column not in adata.var.columns:
        raise ContractError(
            f"assay {assay_id!r} feature ID column is absent: {feature_id_column}"
        )
    graph_summaries = []
    for graph in assay.get("graph_representations", []):
        source = str(graph.get("source", "")).lower()
        key = str(graph.get("key", ""))
        if source == "obsm":
            if key not in adata.obsm:
                raise ContractError(
                    f"assay {assay_id!r} graph representation is absent: {key}"
                )
            shape = adata.obsm[key].shape
        elif source in {"x", "layer"}:
            shape = _matrix(adata, graph).shape
        else:
            raise ContractError(
                f"assay {assay_id!r} graph source must be X, layer, or obsm"
            )
        graph_summaries.append({"source": source, "key": key, "shape": list(shape)})
    representations = {}
    blockers = []
    for name in assay["representations"]:
        declaration = representation_decl(config, assay_id, str(name))
        expression = _matrix(adata, declaration)
        detection = declaration.get("detection")
        detection_summary = {"available": False, "required": False}
        if detection is not None:
            detection_matrix = _matrix(adata, detection)
            minimum = _minimum(detection_matrix)
            valid = minimum >= 0
            required = bool(detection.get("required", False))
            detection_summary = {
                "available": True,
                "required": required,
                "source": str(detection["source"]),
                "key": str(detection.get("key", "")),
                "minimum": minimum,
                "nonnegative": valid,
                "semantics": str(detection.get("semantics", "> 0")),
            }
            if required and not valid:
                blockers.append(
                    f"{assay_id}/{name} requires detection but its declared source is signed"
                )
        representations[str(name)] = {
            "source": str(declaration["source"]),
            "key": str(declaration.get("key", "")),
            "shape": list(map(int, expression.shape)),
            "ranking_backend": str(declaration["ranking_backend"]),
            "feature_space": str(
                declaration.get("feature_space", assay["feature_space"])
            ),
            "ontology_compatible": bool(
                declaration.get("ontology_compatible", False)
            ),
            "detection": detection_summary,
        }
    return (
        {
            "kind": str(assay["kind"]),
            "feature_space": str(assay["feature_space"]),
            "n_features": int(adata.n_vars),
            "feature_id_column": feature_id_column,
            "graph_representations": graph_summaries,
            "layers": sorted(map(str, adata.layers.keys())),
            "obsm": sorted(map(str, adata.obsm.keys())),
            "representations": representations,
        },
        blockers,
    )


def inspect_input(config: dict) -> dict:
    declaration = config["input"]
    path = Path(declaration["path"])
    if not path.is_file():
        raise ContractError(f"input object does not exist: {path}")
    observed_sha = sha256(path)
    expected_sha = str(declaration.get("expected_sha256", "")).strip()
    if expected_sha and expected_sha != observed_sha:
        raise ContractError("input SHA-256 does not match the frozen declaration")

    kind = str(declaration["kind"]).lower()
    assays = assay_registry(config)
    assay_summaries: dict[str, dict] = {}
    blockers: list[str] = []
    if kind == "h5ad":
        if len(assays) != 1:
            raise ContractError("h5ad input must declare exactly one root assay")
        import anndata as ad

        root = ad.read_h5ad(path, backed="r")
        try:
            summary = _inspect_obs(root, config, is_mudata=False)
            assay_id, assay = next(iter(assays.items()))
            modality = str(assay.get("modality", "root"))
            if modality not in {"", "root"}:
                raise ContractError("h5ad assay modality must be root")
            assay_summaries[assay_id], found = _inspect_assay(config, assay_id, root)
            blockers.extend(found)
        finally:
            root.file.close()
    elif kind == "h5mu":
        import mudata as md

        mdata = md.read_h5mu(path, backed=True)
        try:
            summary = _inspect_obs(mdata, config, is_mudata=True)
            root_barcodes = list(mdata.obs_names)
            for assay_id, assay in assays.items():
                modality = str(assay.get("modality", "")).strip()
                if not modality or modality not in mdata.mod:
                    raise ContractError(
                        f"MuData modality for assay {assay_id!r} is absent: {modality!r}"
                    )
                adata = mdata.mod[modality]
                if list(adata.obs_names) != root_barcodes:
                    raise ContractError(
                        f"assay {assay_id!r} barcode order differs from MuData root"
                    )
                assay_summaries[assay_id], found = _inspect_assay(
                    config, assay_id, adata
                )
                assay_summaries[assay_id]["modality"] = modality
                blockers.extend(found)
        finally:
            mdata.file.close()
    else:
        raise ContractError("input.kind must be h5ad or h5mu")

    configured_lanes = evidence_lanes(config)
    reference = config["evidence"].get("validated_reference", {})
    validated_reference = bool(reference.get("available", False)) and bool(
        reference.get("ontology_validated", False)
    )
    if not any(lane["ontology_compatible"] for lane in configured_lanes) and not validated_reference:
        blockers.append(
            "no ontology-compatible evidence lane is declared: biological labels and CL IDs are blocked; label-free audit and subclustering remain available"
        )
    summary.update(
        {
            "input_path": str(path.resolve()),
            "input_sha256": observed_sha,
            "input_kind": kind,
            "assays": assay_summaries,
            "evidence_lanes": configured_lanes,
            "blockers": blockers,
        }
    )
    return summary


def validate_mode(config: dict) -> tuple[str, list[str]]:
    mode = str(config["run"]["mode"])
    if mode not in {"annotation", "evidence-only"}:
        raise ContractError("run.mode must be annotation or evidence-only")
    available = bool(config.get("knowledgebase", {}).get("available", False))
    if mode == "annotation" and not available:
        raise ContractError("annotation mode requires a callable Knowledgebase MCP")
    blockers = []
    if not available:
        blockers.append(
            "Knowledgebase MCP unavailable: labels, CL IDs, finalization, and write-back blocked"
        )
    return mode, blockers


def run(config_path: Path, output_root: Path) -> dict:
    config = read_config(config_path)
    mode, blockers = validate_mode(config)
    summary = inspect_input(config)
    blockers.extend(summary.get("blockers", []))
    expected = sum(
        int(scope["expected_parent_count"])
        for scope in config["canonical_partitions"].values()
    )
    if expected != int(summary["canonical_parent_count"]):
        raise ContractError("declared canonical parent total does not match input")
    now = datetime.now(timezone.utc).isoformat()
    manifest = {
        "schema_version": 2,
        "artifact_kind": "adaptive_assay_input_contract",
        "status": "complete",
        "run_id": config["run"]["run_id"],
        "mode": mode,
        "created_utc": now,
        "input": summary,
        "config": {"path": str(config_path.resolve()), "sha256": sha256(config_path)},
        "blockers": blockers,
        "software": {"python": platform.python_version()},
        "label_fields_allowed": mode == "annotation" and not blockers,
    }
    context = (
        "# Input context\n\n"
        f"- Run ID: `{config['run']['run_id']}`\n"
        f"- Mode: `{mode}`\n"
        f"- Input: `{summary['input_path']}`\n"
        f"- Input SHA-256: `{summary['input_sha256']}`\n"
        f"- Cells: {summary['n_cells']}\n"
        f"- Canonical parents: {summary['canonical_parent_count']}\n"
        f"- Assays: {', '.join(summary['assays'])}\n"
        f"- Blockers: {'; '.join(blockers) if blockers else 'none'}\n"
    )
    state = {
        "schema_version": 2,
        "run_id": config["run"]["run_id"],
        "mode": mode,
        "phase": "INPUT_FROZEN",
        "blockers": blockers,
        "labels_permitted": False,
        "writeback_permitted": False,
        "updated_utc": now,
    }
    atomic_write_text(output_root / "context.md", context)
    atomic_write_json(output_root / "run_state.json", state)
    atomic_write_json(output_root / "input_contract.manifest.json", manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        run(args.config, args.output_root)
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()

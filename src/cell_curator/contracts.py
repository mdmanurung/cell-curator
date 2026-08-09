#!/usr/bin/env python3
"""Shared deterministic contracts for the bundled skill utilities."""

from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any

import yaml

OUTCOMES = (
    "ACCEPT SPLIT",
    "RETAIN PARENT",
    "DOUBLET/MIXED",
    "TECHNICAL/UNRESOLVED",
)
COMPARATOR_ROLES = ("parent_remainder", "closest_sibling")
RAPIDS_BACKEND = "rapids_singlecell"


class ContractError(RuntimeError):
    """Raised when an artifact violates a reusable workflow contract."""


def read_config(path: Path) -> dict[str, Any]:
    # Import lazily to keep lightweight hash/artifact utilities importable without
    # constructing the full Pydantic model.
    from .config import load_config

    try:
        model = load_config(path)
    except (OSError, ValueError, yaml.YAMLError) as exc:
        raise ContractError(str(exc)) from exc
    value = model.model_dump(mode="json")
    evidence_lanes(value)
    return value


def assay_registry(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    assays = config["input"]["assays"]
    for assay_id, assay in assays.items():
        if not isinstance(assay, dict):
            raise ContractError(f"assay {assay_id!r} must be a mapping")
        if not str(assay.get("kind", "")).strip():
            raise ContractError(f"assay {assay_id!r} lacks kind")
        if not str(assay.get("feature_space", "")).strip():
            raise ContractError(f"assay {assay_id!r} lacks feature_space")
        representations = assay.get("representations")
        if not isinstance(representations, dict) or not representations:
            raise ContractError(f"assay {assay_id!r} lacks representations")
    return assays


def representation_decl(
    config: dict[str, Any], assay_id: str, representation: str
) -> dict[str, Any]:
    assays = assay_registry(config)
    if assay_id not in assays:
        raise ContractError(f"evidence lane refers to unknown assay {assay_id!r}")
    representations = assays[assay_id]["representations"]
    if representation not in representations:
        raise ContractError(f"assay {assay_id!r} lacks representation {representation!r}")
    declaration = representations[representation]
    source = str(declaration.get("source", "")).lower()
    if source not in {"x", "layer"}:
        raise ContractError(f"{assay_id}/{representation} source must be X or layer")
    if source == "layer" and not str(declaration.get("key", "")).strip():
        raise ContractError(f"{assay_id}/{representation} layer key is blank")
    backend = str(declaration.get("ranking_backend", "")).strip()
    if not backend:
        raise ContractError(f"{assay_id}/{representation} lacks ranking_backend")
    detection = declaration.get("detection")
    if detection is not None:
        if not isinstance(detection, dict):
            raise ContractError(
                f"{assay_id}/{representation} detection must be a mapping or null"
            )
        detection_source = str(detection.get("source", "")).lower()
        if detection_source not in {"x", "layer"}:
            raise ContractError(
                f"{assay_id}/{representation} detection source must be X or layer"
            )
        if detection_source == "layer" and not str(detection.get("key", "")).strip():
            raise ContractError(f"{assay_id}/{representation} detection layer key is blank")
    return declaration


def evidence_lanes(config: dict[str, Any]) -> list[dict[str, Any]]:
    lanes: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw in config["evidence"]["lanes"]:
        if not isinstance(raw, dict):
            raise ContractError("each evidence lane must be a mapping")
        assay_id = str(raw.get("assay", "")).strip()
        representation = str(raw.get("representation", "")).strip()
        declaration = representation_decl(config, assay_id, representation)
        lane_key = str(raw.get("lane_key") or f"{assay_id}__{representation}")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", lane_key):
            raise ContractError(f"invalid evidence lane_key {lane_key!r}")
        if lane_key in seen:
            raise ContractError(f"duplicate evidence lane_key {lane_key!r}")
        seen.add(lane_key)
        assay = config["input"]["assays"][assay_id]
        lanes.append(
            {
                "lane_key": lane_key,
                "assay": assay_id,
                "assay_kind": str(assay["kind"]),
                "feature_space": str(declaration.get("feature_space", assay["feature_space"])),
                "representation": representation,
                "backend": str(declaration["ranking_backend"]),
                "required_for_acceptance": as_bool(raw.get("required_for_acceptance", True)),
                "ontology_compatible": as_bool(declaration.get("ontology_compatible", False)),
                "detection_required": bool(
                    declaration.get("detection", {}).get("required", False)
                    if declaration.get("detection") is not None
                    else False
                ),
            }
        )
    return lanes


def lane_gate_column(lane_key: str) -> str:
    return f"lane__{lane_key}__evidence_pass"


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def values_sha256(values: list[str], *, sort_values: bool = False) -> str:
    normalized = sorted(map(str, values)) if sort_values else list(map(str, values))
    digest = hashlib.sha256()
    for value in normalized:
        digest.update(value.encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def natural_key(value: object) -> tuple[object, ...]:
    return tuple(
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", str(value))
    )


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no", ""}:
        return False
    raise ContractError(f"cannot interpret Boolean value {value!r}")


def require_blank_labels(frame, *, prefix: str = "biological_label_") -> None:
    columns = [column for column in frame.columns if str(column).startswith(prefix)]
    columns += [column for column in ("cell_type", "cl_id") if column in frame]
    for column in columns:
        if frame[column].fillna("").astype(str).str.strip().ne("").any():
            raise ContractError(f"label field {column!r} must remain blank")


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ContractError(f"refusing to overwrite immutable output: {path}")
    partial = path.with_name(f".{path.name}.partial.{os.getpid()}")
    try:
        partial.write_text(text)
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_write_tsv(path: Path, frame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ContractError(f"refusing to overwrite immutable output: {path}")
    partial = path.with_name(f".{path.name}.partial.{os.getpid()}")
    try:
        frame.to_csv(partial, sep="\t", index=False, lineterminator="\n")
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def require_columns(frame, columns: set[str], *, artifact: str) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ContractError(f"{artifact} misses columns: {sorted(missing)}")

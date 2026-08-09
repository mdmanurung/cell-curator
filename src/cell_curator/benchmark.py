"""Neutral, dependency-free evaluation of exported annotation predictions."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pydantic import ValidationError

from .contracts import ContractError, sha256
from .models import BenchmarkDatasetManifest, BenchmarkResult
from .provenance import publish_json

REQUIRED_COLUMNS = {"barcode", "true_label", "predicted_label"}


def _class_metrics(truth: pd.Series, prediction: pd.Series) -> tuple[float, float, float]:
    labels = sorted(set(truth.astype(str)) | set(prediction.astype(str)))
    f1s = []
    total_tp = total_fp = total_fn = 0
    for label in labels:
        true = truth.astype(str).eq(label)
        predicted = prediction.astype(str).eq(label)
        tp = int((true & predicted).sum())
        fp = int((~true & predicted).sum())
        fn = int((true & ~predicted).sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1s.append(2 * precision * recall / (precision + recall) if precision + recall else 0.0)
        total_tp += tp
        total_fp += fp
        total_fn += fn
    micro_precision = total_tp / (total_tp + total_fp) if total_tp + total_fp else 0.0
    micro_recall = total_tp / (total_tp + total_fn) if total_tp + total_fn else 0.0
    micro_f1 = (
        2 * micro_precision * micro_recall / (micro_precision + micro_recall)
        if micro_precision + micro_recall
        else 0.0
    )
    accuracy = float(truth.astype(str).eq(prediction.astype(str)).mean())
    return accuracy, float(np.mean(f1s)) if f1s else 0.0, micro_f1


def _unknown_metrics(
    truth: pd.Series, prediction: pd.Series
) -> tuple[float | None, float | None]:
    true_unknown = truth.astype(str).eq("Unknown")
    predicted_unknown = prediction.astype(str).eq("Unknown")
    if not true_unknown.any() and not predicted_unknown.any():
        return None, None
    tp = int((true_unknown & predicted_unknown).sum())
    precision = tp / int(predicted_unknown.sum()) if predicted_unknown.any() else 0.0
    recall = tp / int(true_unknown.sum()) if true_unknown.any() else 0.0
    return precision, recall


def _calibration(frame: pd.DataFrame, bins: int = 10) -> tuple[float | None, float | None]:
    if "confidence" not in frame:
        return None, None
    confidence = pd.to_numeric(frame["confidence"], errors="coerce").to_numpy()
    if not np.isfinite(confidence).all() or ((confidence < 0) | (confidence > 1)).any():
        raise ContractError("benchmark confidence must be finite and within [0, 1]")
    correct = (
        frame["true_label"]
        .astype(str)
        .eq(frame["predicted_label"].astype(str))
        .astype(float)
        .to_numpy()
    )
    brier = float(np.mean((confidence - correct) ** 2))
    edges = np.linspace(0, 1, bins + 1)
    ece = 0.0
    for left, right in zip(edges[:-1], edges[1:], strict=True):
        selected = (confidence >= left) & (
            confidence <= right if right == 1 else confidence < right
        )
        if selected.any():
            ece += selected.mean() * abs(
                float(confidence[selected].mean() - correct[selected].mean())
            )
    return brier, float(ece)


def _hierarchy_error(frame: pd.DataFrame) -> float | None:
    if not {"true_path", "predicted_path"}.issubset(frame.columns):
        return None
    distances = []
    for true_path, predicted_path in zip(
        frame["true_path"], frame["predicted_path"], strict=True
    ):
        truth = [item for item in str(true_path).split("/") if item]
        prediction = [item for item in str(predicted_path).split("/") if item]
        common = 0
        for left, right in zip(truth, prediction, strict=False):
            if left != right:
                break
            common += 1
        distance = (len(truth) - common) + (len(prediction) - common)
        maximum = max(1, len(truth) + len(prediction))
        distances.append(distance / maximum)
    return float(np.mean(distances)) if distances else None


def _binary_f1(frame: pd.DataFrame, truth_column: str, prediction_column: str) -> float | None:
    if not {truth_column, prediction_column}.issubset(frame.columns):
        return None

    def boolean(value: Any) -> bool:
        return str(value).strip().lower() in {"1", "true", "yes", "mixed", "doublet"}

    truth = frame[truth_column].map(boolean)
    predicted = frame[prediction_column].map(boolean)
    tp = int((truth & predicted).sum())
    fp = int((~truth & predicted).sum())
    fn = int((truth & ~predicted).sum())
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0


def _optional_agreement(frame: pd.DataFrame, left: str, right: str) -> float | None:
    if not {left, right}.issubset(frame.columns):
        return None
    return float(frame[left].astype(str).eq(frame[right].astype(str)).mean())


def _optional_boolean_fraction(frame: pd.DataFrame, column: str) -> float | None:
    if column not in frame:
        return None
    values = frame[column].astype(str).str.lower().isin({"1", "true", "yes", "pass"})
    return float(values.mean())


def _optional_numeric(
    frame: pd.DataFrame, column: str, *, aggregation: str = "max"
) -> float | None:
    if column not in frame:
        return None
    values = pd.to_numeric(frame[column], errors="coerce")
    if values.isna().all():
        return None
    return float(values.max() if aggregation == "max" else values.mean())


def evaluate_predictions(
    frame: pd.DataFrame,
    *,
    system: str,
    dataset: str,
) -> BenchmarkResult:
    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise ContractError(f"neutral prediction table misses columns: {sorted(missing)}")
    if frame["barcode"].astype(str).duplicated().any():
        raise ContractError("neutral prediction table must contain one row per barcode")
    accuracy, macro_f1, micro_f1 = _class_metrics(frame["true_label"], frame["predicted_label"])
    called = ~frame["predicted_label"].astype(str).eq("Unknown")
    coverage = float(called.mean())
    selective_accuracy = (
        float(
            frame.loc[called, "true_label"]
            .astype(str)
            .eq(frame.loc[called, "predicted_label"].astype(str))
            .mean()
        )
        if called.any()
        else 0.0
    )
    unknown_precision, unknown_recall = _unknown_metrics(
        frame["true_label"], frame["predicted_label"]
    )
    brier, ece = _calibration(frame)
    hierarchy_error = _hierarchy_error(frame)
    unknown_f1 = (
        2 * unknown_precision * unknown_recall / (unknown_precision + unknown_recall)
        if unknown_precision is not None
        and unknown_recall is not None
        and unknown_precision + unknown_recall
        else None
    )
    return BenchmarkResult(
        system=system,
        dataset=dataset,
        n=len(frame),
        accuracy=accuracy,
        macro_f1=macro_f1,
        micro_f1=micro_f1,
        coverage=coverage,
        selective_accuracy=selective_accuracy,
        coverage_risk=1.0 - selective_accuracy,
        hierarchical_accuracy=(1.0 - hierarchy_error) if hierarchy_error is not None else None,
        unknown_precision=unknown_precision,
        unknown_recall=unknown_recall,
        unknown_f1=unknown_f1,
        brier_score=brier,
        expected_calibration_error=ece,
        hierarchy_error=hierarchy_error,
        ontology_distance_error=hierarchy_error,
        mixed_detection_f1=_binary_f1(frame, "true_mixed", "predicted_mixed"),
        split_decision_accuracy=_optional_agreement(
            frame, "true_split_decision", "predicted_split_decision"
        ),
        repeat_reproducibility=_optional_agreement(
            frame, "predicted_label", "repeat_predicted_label"
        ),
        input_invariance=_optional_boolean_fraction(frame, "input_invariant"),
        unsplit_parent_invariance=_optional_boolean_fraction(frame, "unsplit_parent_invariant"),
        runtime_seconds=_optional_numeric(frame, "runtime_seconds"),
        peak_memory_mb=_optional_numeric(frame, "peak_memory_mb"),
        artifact_completeness=_optional_boolean_fraction(frame, "artifact_complete"),
    )


def evaluate_file(
    prediction_path: Path,
    *,
    system: str,
    dataset: str,
    output_path: Path,
) -> BenchmarkResult:
    separator = "\t" if prediction_path.suffix.lower() == ".tsv" else ","
    frame = pd.read_csv(prediction_path, sep=separator, keep_default_na=False)
    result = evaluate_predictions(frame, system=system, dataset=dataset)
    payload = {
        **result.model_dump(mode="json"),
        "prediction_path": str(prediction_path.resolve()),
        "prediction_sha256": sha256(prediction_path),
        "neutral_format": True,
    }
    publish_json(output_path, payload)
    return result


def validate_dataset_manifest(path: Path) -> dict[str, Any]:
    try:
        manifest = BenchmarkDatasetManifest.model_validate_json(path.read_text())
    except (OSError, ValidationError, json.JSONDecodeError) as exc:
        raise ContractError(f"invalid benchmark manifest {path}: {exc}") from exc
    value = manifest.model_dump(mode="json", exclude_none=True)
    if manifest.status == "available":
        local = Path(manifest.local_path or "")
        if not local.is_absolute() and not local.is_file():
            # Manifests ship under <package-data>/benchmarks/manifests while
            # local_path is rooted at <package-data> (and at the repository
            # root during development). This keeps the same manifest portable.
            packaged = path.resolve().parents[2] / local
            if packaged.is_file():
                local = packaged
        if not local.is_file() or sha256(local) != manifest.sha256:
            raise ContractError("available benchmark data are absent or fail checksum")
    return value


def write_pending_status(
    output_path: Path,
    *,
    dataset: str,
    modality: str,
    reason: str,
    reproduce: list[str],
) -> None:
    publish_json(
        output_path,
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_benchmark_status",
            "dataset": dataset,
            "modality": modality,
            "status": "pending-external-resource",
            "reason": reason,
            "reproduce": reproduce,
            "superiority_claimed": False,
        },
    )

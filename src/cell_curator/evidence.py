"""Deterministic marker, program, multimodal, and abstention evidence engine."""

from __future__ import annotations

import json
import logging
import math
import warnings
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

from .capabilities import build_capability_registry
from .config import CellCuratorConfig
from .contracts import ContractError, as_bool, sha256
from .data import (
    LoadedInput,
    close_input,
    feature_ids_for,
    load_input,
    matrix_for,
    resolve_obsm,
)
from .models import AnnotationDecision, Capability, CapabilityLevel, MarkerProgram
from .provenance import make_receipt, publish_json


@dataclass(frozen=True)
class LaneResult:
    lane_key: str
    evidence: pd.DataFrame
    bottom_up: pd.DataFrame
    receipt_path: Path


LOG = logging.getLogger(__name__)


def validate_backend_receipt(receipt: Any, capability: Capability) -> None:
    """Enforce the common no-fallback contract and strict RAPIDS method receipt."""

    if receipt.fallback_used or receipt.cpu_fallback_available or receipt.cpu_fallback_used:
        raise ContractError("backend receipt records a prohibited CPU fallback")
    if not capability.requires_gpu:
        return
    if set(receipt.marker_methods) != {"wilcoxon", "logreg"}:
        raise ContractError("RAPIDS receipt must record native Wilcoxon and balanced logreg")
    if receipt.convergence_validated is not True:
        raise ContractError("RAPIDS receipt lacks successful logreg convergence validation")
    required_runtime = {
        "device_count",
        "device_name",
        "compute_capability",
        "validated_gpu_contract",
        "rapids_singlecell",
        "rapids_distribution",
        "cuda_major",
        "cuda_runtime_version",
        "cuda_driver_version",
        "native_api",
    }
    missing = required_runtime.difference(receipt.runtime)
    if missing or int(receipt.runtime.get("device_count", 0)) != 1:
        raise ContractError(f"RAPIDS receipt runtime contract is incomplete: {sorted(missing)}")


def lane_key(assay: str, representation: str, configured: str | None = None) -> str:
    return configured or f"{assay}__{representation}"


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    if values.size == 0:
        return values
    values = np.nan_to_num(values, nan=1.0, posinf=1.0, neginf=1.0)
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def _matrix_rows(matrix: Any, mask: np.ndarray) -> np.ndarray:
    selected = matrix[mask]
    return selected.toarray() if sparse.issparse(selected) else np.asarray(selected)


def compute_bottom_up_markers(
    matrix: Any,
    groups: pd.Series,
    features: pd.Index,
    *,
    detection: Any | None,
    top_n: int = 100,
) -> pd.DataFrame:
    """Compute one-vs-remainder Mann-Whitney evidence with deterministic BH correction."""

    if matrix.shape[0] != len(groups) or matrix.shape[1] != len(features):
        raise ContractError("matrix, groups, and feature dimensions disagree")
    rows: list[dict[str, Any]] = []
    group_values = groups.astype(str).to_numpy()
    for cluster_id in sorted(pd.unique(group_values)):
        target_mask = group_values == cluster_id
        reference_mask = ~target_mask
        if not target_mask.any() or not reference_mask.any():
            continue
        target = _matrix_rows(matrix, target_mask)
        reference = _matrix_rows(matrix, reference_mask)
        target_mean = np.asarray(target.mean(axis=0)).ravel()
        reference_mean = np.asarray(reference.mean(axis=0)).ravel()
        effect = target_mean - reference_mean
        pvalues = []
        for column in range(len(features)):
            try:
                pvalue = stats.mannwhitneyu(
                    target[:, column],
                    reference[:, column],
                    alternative="two-sided",
                    method="asymptotic",
                ).pvalue
            except ValueError:
                pvalue = 1.0
            pvalues.append(float(pvalue))
        adjusted = benjamini_hochberg(np.asarray(pvalues))
        order = np.lexsort((features.astype(str).to_numpy(), -effect))
        selected = order[: min(top_n, len(order))]
        if detection is not None:
            target_detect = (_matrix_rows(detection, target_mask)[:, selected] > 0).mean(axis=0)
            reference_detect = (_matrix_rows(detection, reference_mask)[:, selected] > 0).mean(
                axis=0
            )
        else:
            target_detect = np.full(len(selected), np.nan)
            reference_detect = np.full(len(selected), np.nan)
        for rank, (column, pvalue, padj, target_fraction, reference_fraction) in enumerate(
            zip(
                selected,
                np.asarray(pvalues)[selected],
                adjusted[selected],
                target_detect,
                reference_detect,
                strict=True,
            ),
            1,
        ):
            rows.append(
                {
                    "cluster_id": cluster_id,
                    "feature": str(features[column]),
                    "rank": rank,
                    "target_mean": float(target_mean[column]),
                    "reference_mean": float(reference_mean[column]),
                    "effect": float(effect[column]),
                    "pvalue": pvalue,
                    "adjusted_pvalue": float(padj),
                    "target_fraction_detected": float(target_fraction),
                    "reference_fraction_detected": float(reference_fraction),
                }
            )
    columns = [
        "cluster_id",
        "feature",
        "rank",
        "target_mean",
        "reference_mean",
        "effect",
        "pvalue",
        "adjusted_pvalue",
        "target_fraction_detected",
        "reference_fraction_detected",
    ]
    return pd.DataFrame(rows, columns=columns)


def compute_gpu_bottom_up_markers(
    matrix: Any,
    groups: pd.Series,
    features: pd.Index,
    *,
    detection: Any | None,
    rsc: Any,
    marker_config: dict[str, Any],
    top_n: int = 100,
) -> pd.DataFrame:
    """Run native RAPIDS Wilcoxon and balanced logreg without a CPU fallback."""

    import anndata as ad

    categories = sorted(groups.astype(str).unique())
    if len(categories) < 2:
        raise ContractError(
            "RAPIDS marker evidence requires at least two observed groups; "
            "no marker methods were executed"
        )
    gpu_input = ad.AnnData(
        X=matrix.copy(),
        obs=pd.DataFrame(
            {"cell_curator_group": pd.Categorical(groups.astype(str), categories=categories)},
            index=groups.index.astype(str),
        ),
        var=pd.DataFrame(index=features.astype(str)),
    )
    rsc.get.anndata_to_GPU(gpu_input)
    wilcoxon_key = "cell_curator_native_rapids_wilcoxon"
    rsc.tl.rank_genes_groups(
        gpu_input,
        groupby="cell_curator_group",
        groups=categories,
        reference="rest",
        n_genes=min(top_n, len(features)),
        rankby_abs=False,
        pts=True,
        key_added=wilcoxon_key,
        method="wilcoxon",
        corr_method="benjamini-hochberg",
        tie_correct=bool(marker_config["wilcoxon"]["tie_correct"]),
        use_continuity=bool(marker_config["wilcoxon"]["continuity_correct"]),
        use_raw=False,
        layer=None,
        chunk_size=int(marker_config["wilcoxon"]["chunk_size"]),
        multi_gpu=False,
    )
    logreg_key = "cell_curator_native_rapids_logreg"
    with warnings.catch_warnings(record=True) as observed_warnings:
        warnings.simplefilter("always")
        rsc.tl.rank_genes_groups(
            gpu_input,
            groupby="cell_curator_group",
            groups=categories,
            reference="rest",
            n_genes=min(top_n, len(features)),
            rankby_abs=False,
            pts=True,
            key_added=logreg_key,
            method="logreg",
            use_raw=False,
            layer=None,
            class_weight="balanced",
            penalty="l2",
            C=float(marker_config["logreg"]["C"]),
            fit_intercept=True,
            max_iter=int(marker_config["logreg"]["retry_max_iter"]),
            tol=float(marker_config["logreg"]["tol"]),
            solver="qn",
        )
    convergence_warnings = [
        str(item.message)
        for item in observed_warnings
        if "converg" in str(item.message).casefold()
    ]
    if convergence_warnings:
        raise ContractError(
            "native RAPIDS logistic regression did not converge: "
            + "; ".join(convergence_warnings)
        )
    positions = {str(feature): index for index, feature in enumerate(features)}
    group_values = groups.astype(str).to_numpy()
    rows: list[dict[str, Any]] = []
    for cluster_id in categories:
        target_mask = group_values == cluster_id
        reference_mask = ~target_mask
        target = _matrix_rows(matrix, target_mask)
        reference = _matrix_rows(matrix, reference_mask)
        target_mean = np.asarray(target.mean(axis=0)).ravel()
        reference_mean = np.asarray(reference.mean(axis=0)).ravel()
        for method, key in (("wilcoxon", wilcoxon_key), ("logreg", logreg_key)):
            result = gpu_input.uns[key]
            names = result["names"][cluster_id]
            for rank, feature in enumerate(names, 1):
                feature = str(feature)
                column = positions[feature]

                def scalar(
                    field: str,
                    default: float = float("nan"),
                    *,
                    group: str = cluster_id,
                    row_rank: int = rank,
                    method_result: dict[str, Any] = result,
                ) -> float:
                    value = method_result.get(field)
                    if value is None or group not in (value.dtype.names or ()):
                        return default
                    return float(value[group][row_rank - 1])

                if detection is None:
                    target_fraction = reference_fraction = float("nan")
                else:
                    target_fraction = float(
                        (_matrix_rows(detection, target_mask)[:, column] > 0).mean()
                    )
                    reference_fraction = float(
                        (_matrix_rows(detection, reference_mask)[:, column] > 0).mean()
                    )
                rows.append(
                    {
                        "cluster_id": cluster_id,
                        "method": method,
                        "feature": feature,
                        "rank": rank,
                        "target_mean": float(target_mean[column]),
                        "reference_mean": float(reference_mean[column]),
                        "effect": (
                            scalar("scores")
                            if method == "logreg"
                            else float(target_mean[column] - reference_mean[column])
                        ),
                        "pvalue": scalar("pvals") if method == "wilcoxon" else float("nan"),
                        "adjusted_pvalue": (
                            scalar("pvals_adj") if method == "wilcoxon" else float("nan")
                        ),
                        "target_fraction_detected": target_fraction,
                        "reference_fraction_detected": reference_fraction,
                    }
                )
    result = pd.DataFrame(rows)
    positive_methods = {
        (str(row.cluster_id), str(row.feature), str(row.method))
        for row in result.itertuples(index=False)
        if float(row.effect) > 0
        and (
            str(row.method) == "logreg"
            or (
                np.isfinite(float(row.adjusted_pvalue))
                and float(row.adjusted_pvalue) <= 1.0
            )
        )
    }
    result["method_concordant"] = [
        all(
            (str(row.cluster_id), str(row.feature), method) in positive_methods
            for method in ("wilcoxon", "logreg")
        )
        for row in result.itertuples(index=False)
    ]
    return result


def _group_summaries(
    matrix: Any,
    detection: Any | None,
    groups: pd.Series,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    labels = sorted(groups.astype(str).unique())
    group_values = groups.astype(str).to_numpy()
    means = []
    fractions = []
    for label in labels:
        mask = group_values == label
        values = _matrix_rows(matrix, mask)
        means.append(np.asarray(values.mean(axis=0)).ravel())
        if detection is None:
            fractions.append(np.full(matrix.shape[1], np.nan))
        else:
            fractions.append(
                np.asarray((_matrix_rows(detection, mask) > 0).mean(axis=0)).ravel()
            )
    return labels, np.vstack(means), np.vstack(fractions)


def _zscore_columns(values: np.ndarray) -> np.ndarray:
    center = np.nanmean(values, axis=0)
    spread = np.nanstd(values, axis=0)
    spread = np.where(spread <= 1e-12, 1.0, spread)
    return (values - center) / spread


def score_marker_programs(
    matrix: Any,
    groups: pd.Series,
    features: pd.Index,
    programs: list[MarkerProgram],
    *,
    detection: Any | None,
    run_id: str,
    assay: str,
    representation: str,
    capability: Capability,
    positive_detection_threshold: float = 0.1,
) -> pd.DataFrame:
    """Score signed programs without treating unmeasured features as negatives."""

    cluster_ids, means, fractions = _group_summaries(matrix, detection, groups)
    standardized = _zscore_columns(means)
    positions = {str(feature): index for index, feature in enumerate(features)}
    rows: list[dict[str, Any]] = []
    for program in programs:
        compatible_spaces = {capability.feature_space}
        if capability.assay_kind.value == "gene_activity":
            compatible_spaces.add("gene")
        if program.feature_space not in compatible_spaces:
            continue
        positive = [positions[item] for item in program.positive if item in positions]
        negative = [positions[item] for item in program.negative if item in positions]
        missing_positive = tuple(item for item in program.positive if item not in positions)
        missing_negative = tuple(item for item in program.negative if item not in positions)
        if not positive:
            continue
        program_coverage = len(positive) / len(program.positive)
        for row_index, cluster_id in enumerate(cluster_ids):
            positive_signal = float(np.mean(standardized[row_index, positive]))
            negative_signal = (
                float(np.mean(np.maximum(standardized[row_index, negative], 0.0)))
                if negative
                else 0.0
            )
            if detection is not None:
                positive_coverage = float(
                    np.mean(fractions[row_index, positive] >= positive_detection_threshold)
                )
                if missing_negative:
                    negative_violation = float("nan")
                elif negative:
                    negative_violation = float(
                        np.mean(fractions[row_index, negative] >= positive_detection_threshold)
                    )
                else:
                    negative_violation = 0.0
            else:
                positive_coverage = float("nan")
                negative_violation = float("nan")
            score = positive_signal - negative_signal - (
                negative_violation if np.isfinite(negative_violation) else 0.0
            )
            rows.append(
                {
                    "schema_version": 3,
                    "run_id": run_id,
                    "cluster_id": cluster_id,
                    "label": program.label,
                    "axis": program.axis,
                    "parent_label": program.parent or "",
                    "cl_id": program.cl_id or "",
                    "assay": assay,
                    "representation": representation,
                    "capability_level": capability.level.value,
                    "feature_space": capability.feature_space,
                    "score": score,
                    "positive_signal": positive_signal,
                    "negative_signal": negative_signal,
                    "positive_coverage": positive_coverage,
                    "negative_violation": negative_violation,
                    "program_coverage": program_coverage,
                    "measured_positive": len(positive),
                    "total_positive": len(program.positive),
                    "measured_negative": len(negative),
                    "total_negative": len(program.negative),
                    "missing_positive": ";".join(missing_positive),
                    "missing_negative": ";".join(missing_negative),
                    "detection_valid": detection is not None,
                    "source": program.source,
                }
            )
    columns = [
        "schema_version",
        "run_id",
        "cluster_id",
        "label",
        "axis",
        "parent_label",
        "cl_id",
        "assay",
        "representation",
        "capability_level",
        "feature_space",
        "score",
        "positive_signal",
        "negative_signal",
        "positive_coverage",
        "negative_violation",
        "program_coverage",
        "measured_positive",
        "total_positive",
        "measured_negative",
        "total_negative",
        "missing_positive",
        "missing_negative",
        "detection_valid",
        "source",
    ]
    return pd.DataFrame(rows, columns=columns)


def _softmax(values: np.ndarray) -> np.ndarray:
    shifted = values - np.max(values)
    exponent = np.exp(np.clip(shifted, -700, 700))
    return exponent / exponent.sum()


def _normalized_entropy(probabilities: np.ndarray) -> float:
    if len(probabilities) <= 1:
        return 0.0
    values = probabilities[probabilities > 0]
    return float(-(values * np.log(values)).sum() / math.log(len(probabilities)))


def _call_cell_state(
    evidence: pd.DataFrame,
    capability_by_key: dict[tuple[str, str], Capability],
    config: CellCuratorConfig,
) -> str:
    if evidence.empty:
        return "not_assessed"
    state_scores = (
        evidence.groupby("label", observed=True)
        .agg(
            score=("score", "mean"),
            positive_coverage=("positive_coverage", "mean"),
            negative_violation=("negative_violation", "max"),
            program_coverage=("program_coverage", "min"),
            detection_valid=("detection_valid", "all"),
        )
        .reset_index()
        .sort_values(["score", "label"], ascending=[False, True])
    )
    winner = state_scores.iloc[0]
    direct = any(
        capability_by_key[(str(item["assay"]), str(item["representation"]))].level
        is CapabilityLevel.ANNOTATION
        for item in evidence[evidence["label"].astype(str).eq(str(winner["label"]))].to_dict(
            "records"
        )
    )
    passes = bool(
        direct
        and bool(winner["detection_valid"])
        and float(winner["program_coverage"]) >= config.evidence.min_positive_coverage
        and float(winner["positive_coverage"]) >= config.evidence.min_positive_coverage
        and float(winner["negative_violation"]) <= config.evidence.max_negative_violation
    )
    return str(winner["label"]) if passes else "Unknown"


def aggregate_multimodal_evidence(
    evidence: pd.DataFrame,
    capabilities: list[Capability],
    config: CellCuratorConfig,
    *,
    scope: str,
    cluster_ids: list[str],
) -> tuple[pd.DataFrame, list[AnnotationDecision]]:
    """Aggregate independent lanes while preserving disagreement and missing gates."""

    capability_by_key = {(item.assay, item.representation): item for item in capabilities}
    required_lanes = {
        lane_key(lane.assay, lane.representation, lane.lane_key)
        for lane in config.evidence.lanes
        if lane.required_for_acceptance
    }
    rows: list[dict[str, Any]] = []
    decisions: list[AnnotationDecision] = []
    for cluster_id in cluster_ids:
        if "cluster_id" in evidence.columns:
            block = evidence[evidence["cluster_id"].astype(str).eq(str(cluster_id))].copy()
        else:
            block = pd.DataFrame(columns=["cluster_id", "axis", "parent_label"])
        if "axis" not in block:
            block["axis"] = "cell_type"
        parent_name = config.evidence.parent_vocabulary.get(str(cluster_id), str(cluster_id))
        if not block.empty:
            compatible_parent = (
                block["parent_label"].astype(str).isin({"", str(cluster_id), parent_name})
            )
            block = block[compatible_parent].copy()
        state_block = block[block["axis"].astype(str).eq("cell_state")].copy()
        block = block[block["axis"].astype(str).eq("cell_type")].copy()
        state_label = _call_cell_state(state_block, capability_by_key, config)
        if block.empty:
            decisions.append(
                AnnotationDecision(
                    run_id=config.run.run_id,
                    scope=scope,
                    cluster_id=str(cluster_id),
                    cell_type="Unknown",
                    cell_state=state_label,
                    confidence=0.0,
                    margin=0.0,
                    normalized_entropy=1.0,
                    abstained=True,
                    rationale="No locally grounded marker program produced measurable evidence.",
                    required_capabilities_complete=False,
                )
            )
            continue
        if "lane_key" not in block:
            block["lane_key"] = (
                block["assay"].astype(str) + "__" + block["representation"].astype(str)
            )
        if not state_block.empty and "lane_key" not in state_block:
            state_block["lane_key"] = (
                state_block["assay"].astype(str)
                + "__"
                + state_block["representation"].astype(str)
            )
        lane_winners = {
            str(lane): str(
                lane_block.sort_values(["score", "label"], ascending=[False, True]).iloc[0][
                    "label"
                ]
            )
            for lane, lane_block in block.groupby("lane_key", observed=True)
        }
        labels = sorted(block["label"].astype(str).unique())
        label_rows: list[dict[str, Any]] = []
        for label in labels:
            label_block = block[block["label"].astype(str).eq(label)]
            weighted_scores = []
            weights = []
            passing_lanes: set[str] = set()
            annotation_lanes: set[str] = set()
            for item in label_block.to_dict("records"):
                capability = capability_by_key[
                    (str(item["assay"]), str(item["representation"]))
                ]
                weight = next(
                    config.input.assays[capability.assay]
                    .representations[capability.representation]
                    .weight
                    for _ in [0]
                )
                if capability.level is CapabilityLevel.CORROBORATION:
                    weight *= 0.5
                elif capability.level in {CapabilityLevel.QC, CapabilityLevel.UNAVAILABLE}:
                    weight = 0.0
                configured_key = str(
                    item.get("lane_key")
                    or lane_key(str(item["assay"]), str(item["representation"]))
                )
                detection_valid = as_bool(item["detection_valid"])
                lane_passes = bool(
                    detection_valid
                    and np.isfinite(float(item["positive_coverage"]))
                    and np.isfinite(float(item["negative_violation"]))
                    and float(item["program_coverage"])
                    >= config.evidence.min_positive_coverage
                    and float(item["positive_coverage"])
                    >= config.evidence.min_positive_coverage
                    and float(item["negative_violation"])
                    <= config.evidence.max_negative_violation
                )
                if lane_passes:
                    passing_lanes.add(configured_key)
                    if capability.level is CapabilityLevel.ANNOTATION:
                        annotation_lanes.add(configured_key)
                weighted_scores.append(float(item["score"]) * weight)
                weights.append(weight)
            score = sum(weighted_scores) / sum(weights) if sum(weights) else -math.inf
            label_rows.append(
                {
                    "label": label,
                    "axis": "cell_type",
                    "score": score,
                    "positive_coverage": float(label_block["positive_coverage"].mean()),
                    "negative_violation": float(label_block["negative_violation"].max()),
                    "program_coverage": float(label_block["program_coverage"].min()),
                    "annotation_lanes": len(annotation_lanes),
                    "required_program_evidence_complete": required_lanes.issubset(
                        passing_lanes
                    ),
                    "cl_id": next(
                        (str(value) for value in label_block["cl_id"] if str(value).strip()),
                        "",
                    ),
                }
            )
        label_rows.sort(key=lambda item: (-item["score"], item["label"]))
        finite = [item for item in label_rows if np.isfinite(item["score"])]
        if not finite:
            probabilities = np.array([])
        else:
            probabilities = _softmax(np.asarray([item["score"] for item in finite]))
        for item, probability in zip(finite, probabilities, strict=True):
            item["probability"] = float(probability)
            rows.append(
                {
                    "schema_version": 3,
                    "run_id": config.run.run_id,
                    "scope": scope,
                    "cluster_id": str(cluster_id),
                    **item,
                    "lane_winners": ";".join(
                        f"{key}={value}" for key, value in sorted(lane_winners.items())
                    ),
                    "multimodal_disagreement": len(set(lane_winners.values())) > 1,
                }
            )
        if not finite:
            winner = None
            confidence = margin = 0.0
            entropy = 1.0
        else:
            winner = finite[0]
            confidence = float(probabilities[0])
            margin = (
                float(probabilities[0] - probabilities[1])
                if len(probabilities) > 1
                else confidence
            )
            entropy = _normalized_entropy(probabilities)
        observed_required = set(block["lane_key"].astype(str))
        required_complete = required_lanes.issubset(observed_required)
        thresholds = config.evidence
        passes = bool(
            winner
            and config.run.mode.value == "annotation"
            and winner["annotation_lanes"] >= 1
            and winner["required_program_evidence_complete"]
            and winner["program_coverage"] >= thresholds.min_positive_coverage
            and winner["positive_coverage"] >= thresholds.min_positive_coverage
            and winner["negative_violation"] <= thresholds.max_negative_violation
            and confidence >= thresholds.min_confidence
            and margin >= thresholds.min_margin
            and entropy <= thresholds.max_normalized_entropy
            and required_complete
        )
        reasons = []
        if winner is None:
            reasons.append("no finite annotation-capable score")
        else:
            if config.run.mode.value != "annotation":
                reasons.append("run is explicitly evidence-only")
            if winner["annotation_lanes"] < 1:
                reasons.append("no direct annotation-capable lane supported the winner")
            if not winner["required_program_evidence_complete"]:
                reasons.append("required lanes lacked valid detection-dependent program evidence")
            if winner["program_coverage"] < thresholds.min_positive_coverage:
                reasons.append("too little of the signed marker program was measured")
            if not np.isfinite(winner["positive_coverage"]):
                reasons.append("positive-marker detection was unavailable")
            elif winner["positive_coverage"] < thresholds.min_positive_coverage:
                reasons.append("positive-program coverage was insufficient")
            if not np.isfinite(winner["negative_violation"]):
                reasons.append("expected-negative detection was unavailable")
            elif winner["negative_violation"] > thresholds.max_negative_violation:
                reasons.append("expected-negative markers were violated")
            if confidence < thresholds.min_confidence or margin < thresholds.min_margin:
                reasons.append("confidence or winner margin was insufficient")
            if entropy > thresholds.max_normalized_entropy:
                reasons.append("candidate evidence remained too entropic")
        if not required_complete:
            reasons.append("one or more required evidence lanes were missing")
        rejected = tuple(item["label"] for item in finite[1:4])
        decisions.append(
            AnnotationDecision(
                run_id=config.run.run_id,
                scope=scope,
                cluster_id=str(cluster_id),
                cell_type=str(winner["label"]) if passes and winner else "Unknown",
                cell_state=state_label,
                cl_id=(str(winner["cl_id"]) or None) if passes and winner else None,
                confidence=confidence,
                margin=margin,
                normalized_entropy=entropy,
                abstained=not passes,
                rationale=(
                    "Independent signed-program evidence passed coverage, negative-marker, "
                    "uncertainty, and capability gates."
                    if passes
                    else "; ".join(reasons)
                ),
                rejected_alternatives=rejected,
                required_capabilities_complete=required_complete,
            )
        )
    aggregate_columns = [
        "schema_version",
        "run_id",
        "scope",
        "cluster_id",
        "label",
        "axis",
        "score",
        "positive_coverage",
        "negative_violation",
        "program_coverage",
        "annotation_lanes",
        "required_program_evidence_complete",
        "cl_id",
        "probability",
        "lane_winners",
        "multimodal_disagreement",
    ]
    return pd.DataFrame(rows, columns=aggregate_columns), decisions


def spatial_context_evidence(
    coordinates: np.ndarray,
    groups: pd.Series,
    *,
    k: int = 8,
) -> pd.DataFrame:
    if coordinates.ndim != 2 or len(coordinates) != len(groups):
        raise ContractError("spatial coordinates must be a cells-by-dimensions matrix")
    effective_k = min(k + 1, len(coordinates))
    _, indices = cKDTree(coordinates).query(coordinates, k=effective_k)
    if indices.ndim == 1:
        indices = indices[:, None]
    neighbors = indices[:, 1:]
    values = groups.astype(str).to_numpy()
    coherence = (
        np.mean(values[neighbors] == values[:, None], axis=1)
        if neighbors.size
        else np.ones(len(values))
    )
    return (
        pd.DataFrame(
            {
                "cluster_id": values,
                "spatial_neighbor_coherence": coherence,
            }
        )
        .groupby("cluster_id", observed=True, as_index=False)
        .mean()
    )


def _feature_mapping(
    features: pd.Index,
    mapping_path: str,
    *,
    target_column: str,
) -> tuple[Any, pd.Index]:
    mapping = pd.read_csv(mapping_path, sep=None, engine="python")
    required = {"peak", target_column}
    if not required.issubset(mapping.columns):
        raise ContractError(f"feature mapping requires peak and {target_column} columns")
    weights = (
        mapping["weight"].astype(float)
        if "weight" in mapping
        else pd.Series(1.0, index=mapping.index)
    )
    feature_pos = {str(value): index for index, value in enumerate(features)}
    targets = sorted(mapping[target_column].astype(str).unique())
    target_pos = {target: index for index, target in enumerate(targets)}
    rows, columns, values = [], [], []
    for row, weight in zip(mapping.to_dict("records"), weights, strict=True):
        if str(row["peak"]) not in feature_pos:
            continue
        rows.append(feature_pos[str(row["peak"])])
        columns.append(target_pos[str(row[target_column])])
        values.append(float(weight))
    if not rows:
        raise ContractError(f"feature mapping contains no measured peaks: {mapping_path}")
    transform = sparse.csr_matrix(
        (values, (rows, columns)), shape=(len(features), len(targets))
    )
    return transform, pd.Index(targets)


def _lane_materials(
    loaded: LoadedInput,
    config: CellCuratorConfig,
    assay_id: str,
    representation_name: str,
    *,
    key: str,
) -> tuple[Any, Any | None, pd.Index, Any, Any]:
    adata = loaded.assays[assay_id]
    assay = config.input.assays[assay_id]
    representation = assay.representations[representation_name]
    matrix = matrix_for(adata, representation)
    source_features = feature_ids_for(
        adata,
        assay.feature_id_column,
        declaration=representation,
        gene_mapping_path=assay.gene_mapping_path,
    )
    features = source_features
    transform = None
    if representation.peak_to_gene_path:
        transform, features = _feature_mapping(
            features, representation.peak_to_gene_path, target_column="gene"
        )
        matrix = matrix @ transform
    elif representation.motif_mapping_path:
        transform, features = _feature_mapping(
            features, representation.motif_mapping_path, target_column="motif"
        )
        matrix = matrix @ transform
    detection = None
    if representation.detection is not None:
        detection = matrix_for(adata, representation.detection.model_dump())
        detection_features = feature_ids_for(
            adata,
            assay.feature_id_column,
            declaration=representation.detection,
            gene_mapping_path=assay.gene_mapping_path,
        )
        if not detection_features.equals(source_features):
            raise ContractError(
                f"representation and detection feature axes differ for {key}"
            )
        if transform is not None:
            if transform.data.size and float(transform.data.min()) < 0:
                raise ContractError(
                    f"mapped detection for {key} requires nonnegative mapping weights"
                )
            detection = detection @ transform
        if sparse.issparse(detection):
            minimum = min(float(detection.data.min()), 0.0) if detection.nnz else 0.0
        else:
            minimum = float(np.asarray(detection).min())
        if minimum < 0:
            raise ContractError(f"detection source for {key} contains negative values")
    return matrix, detection, features, assay, representation


def _binary_separation(matrix: Any, groups: pd.Series) -> float:
    target = groups.astype(str).eq("candidate").to_numpy()
    reference = groups.astype(str).eq("reference").to_numpy()
    if not target.any() or not reference.any():
        raise ContractError("candidate comparison must contain candidate and reference cells")
    target_values = matrix[target]
    reference_values = matrix[reference]
    difference = (
        np.asarray(target_values.mean(axis=0)).ravel()
        - np.asarray(reference_values.mean(axis=0)).ravel()
    )
    if sparse.issparse(matrix):
        scale = np.sqrt(np.asarray(matrix.power(2).mean(axis=0)).ravel() + 1e-8)
    else:
        scale = np.sqrt(np.mean(np.asarray(matrix) ** 2, axis=0) + 1e-8)
    return float(np.sqrt(np.mean((difference / scale) ** 2)))


def run_lane(
    *,
    config: CellCuratorConfig,
    config_path: Path,
    run_root: Path,
    assay_id: str,
    representation_name: str,
    programs: list[MarkerProgram],
) -> LaneResult:
    capability = next(
        item
        for item in build_capability_registry(config)
        if item.assay == assay_id and item.representation == representation_name
    )
    lane = next(
        item
        for item in config.evidence.lanes
        if item.assay == assay_id and item.representation == representation_name
    )
    key = lane_key(assay_id, representation_name, lane.lane_key)
    lane_dir = run_root / "evidence" / "lanes" / key
    evidence_path = lane_dir / "program_evidence.tsv"
    bottom_up_path = lane_dir / "bottom_up_markers.tsv"
    receipt_path = lane_dir / "receipt.json"
    completed = (receipt_path.is_file(), evidence_path.is_file(), bottom_up_path.is_file())
    if all(completed):
        from .models import BackendReceipt

        receipt = BackendReceipt.model_validate_json(receipt_path.read_text())
        expected_artifacts = {
            "program_evidence": sha256(evidence_path),
            "bottom_up_markers": sha256(bottom_up_path),
        }
        if (
            receipt.status != "complete"
            or receipt.run_id != config.run.run_id
            or receipt.lane_key != key
            or receipt.backend != capability.backend
            or receipt.input_sha256 != sha256(Path(config.input.path))
            or receipt.config_sha256 != sha256(config_path)
            or receipt.output_sha256 != expected_artifacts["program_evidence"]
            or receipt.artifact_sha256 != expected_artifacts
            or receipt.fallback_used
        ):
            raise ContractError(f"evidence lane receipt or artifact hash drifted: {key}")
        validate_backend_receipt(receipt, capability)
        return LaneResult(
            key,
            pd.read_csv(evidence_path, sep="\t", keep_default_na=False),
            pd.read_csv(bottom_up_path, sep="\t", keep_default_na=False),
            receipt_path,
        )
    if any(completed):
        raise ContractError(f"evidence lane is partially materialized and cannot resume: {key}")
    LOG.info(
        "evidence lane starting: key=%s assay=%s representation=%s",
        key,
        assay_id,
        representation_name,
    )
    gpu_runtime: dict[str, Any] | None = None
    gpu_rsc: Any | None = None
    if capability.requires_gpu:
        # Validate hardware/software before opening the potentially large input.
        from .rank_markers_gpu import inspect_runtime

        _, gpu_rsc, gpu_runtime = inspect_runtime(config.model_dump(mode="json"))
    loaded = load_input(config, backed=False)
    try:
        matrix, detection, features, assay, _ = _lane_materials(
            loaded, config, assay_id, representation_name, key=key
        )
        groups = loaded.root.obs[config.input.keys.canonical_parent].astype(str)
        if capability.requires_gpu:
            bottom_up = compute_gpu_bottom_up_markers(
                matrix,
                groups,
                features,
                detection=detection,
                rsc=gpu_rsc,
                marker_config=config.markers.model_dump(mode="json"),
            )
        else:
            bottom_up = compute_bottom_up_markers(
                matrix,
                groups,
                features,
                detection=detection,
            )
        evidence = score_marker_programs(
            matrix,
            groups,
            features,
            programs,
            detection=detection,
            run_id=config.run.run_id,
            assay=assay_id,
            representation=representation_name,
            capability=capability,
        )
        if assay.kind == "spatial" and config.input.keys.spatial:
            coordinates, _, _ = resolve_obsm(
                loaded, config.input.keys.spatial, preferred_assay=assay_id
            )
            spatial = spatial_context_evidence(
                coordinates, groups
            )
            evidence = evidence.merge(spatial, on="cluster_id", how="left")
        lane_dir.mkdir(parents=True, exist_ok=True)
        for path, frame in ((evidence_path, evidence), (bottom_up_path, bottom_up)):
            payload = frame.to_csv(sep="\t", index=False, lineterminator="\n").encode()
            from .provenance import atomic_publish

            atomic_publish(path, payload)
        receipt = make_receipt(
            run_id=config.run.run_id,
            lane_key=key,
            capability=capability,
            input_path=Path(config.input.path),
            config_path=config_path,
            output_path=evidence_path,
            artifacts={
                "program_evidence": evidence_path,
                "bottom_up_markers": bottom_up_path,
            },
            feature_ids=list(map(str, features)),
            groups=list(groups.astype(str)),
            limitations=capability.reasons,
            runtime=gpu_runtime,
            marker_methods=("wilcoxon", "logreg") if capability.requires_gpu else (),
            convergence_validated=True if capability.requires_gpu else None,
        )
        publish_json(receipt_path, receipt.model_dump(mode="json"))
        LOG.info("evidence lane complete: key=%s n_rows=%d", key, len(evidence))
        return LaneResult(key, evidence, bottom_up, receipt_path)
    finally:
        close_input(loaded)


def run_candidate_lane(
    *,
    config: CellCuratorConfig,
    config_path: Path,
    run_root: Path,
    comparison_lane_id: str,
    programs: list[MarkerProgram],
) -> LaneResult:
    """Execute one manifest-selected binary candidate comparison lane."""

    manifest_path = run_root / "comparators" / "comparator_manifest.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t", keep_default_na=False)
    selected = manifest[manifest["lane_id"].astype(str).eq(comparison_lane_id)]
    if len(selected) != 1:
        raise ContractError(f"candidate lane is absent or duplicated: {comparison_lane_id}")
    row = selected.iloc[0]
    assay_id = str(row["assay"])
    representation_name = str(row["representation"])
    capability = next(
        item
        for item in build_capability_registry(config)
        if item.assay == assay_id and item.representation == representation_name
    )
    lane_dir = run_root / "evidence" / "candidate_lanes" / comparison_lane_id
    receipt_path = lane_dir / "receipt.json"
    evidence_path = lane_dir / "program_evidence.tsv"
    bottom_up_path = lane_dir / "bottom_up_markers.tsv"
    comparison_path = lane_dir / "comparison.json"
    completed = (
        receipt_path.is_file(),
        evidence_path.is_file(),
        bottom_up_path.is_file(),
        comparison_path.is_file(),
    )
    if all(completed):
        from .models import BackendReceipt

        receipt = BackendReceipt.model_validate_json(receipt_path.read_text())
        if (
            receipt.status != "complete"
            or receipt.run_id != config.run.run_id
            or receipt.lane_key != comparison_lane_id
            or receipt.backend != capability.backend
            or receipt.input_sha256 != sha256(Path(config.input.path))
            or receipt.config_sha256 != sha256(config_path)
            or receipt.output_sha256 != sha256(evidence_path)
            or receipt.artifact_sha256
            != {
                "program_evidence": sha256(evidence_path),
                "bottom_up_markers": sha256(bottom_up_path),
                "comparison": sha256(comparison_path),
            }
            or receipt.fallback_used
        ):
            raise ContractError(f"candidate lane receipt drifted: {comparison_lane_id}")
        validate_backend_receipt(receipt, capability)
        comparison = json.loads(comparison_path.read_text())
        groups_path = Path(str(row["groups_path"]))
        if not groups_path.is_file() or comparison.get("groups_sha256") != sha256(groups_path):
            raise ContractError(
                f"candidate comparator group artifact drifted: {comparison_lane_id}"
            )
        return LaneResult(
            comparison_lane_id,
            pd.read_csv(evidence_path, sep="\t", keep_default_na=False),
            pd.read_csv(bottom_up_path, sep="\t", keep_default_na=False),
            receipt_path,
        )
    if any(completed):
        raise ContractError(
            f"candidate lane is partially materialized and cannot resume: {comparison_lane_id}"
        )
    gpu_runtime: dict[str, Any] | None = None
    gpu_rsc: Any | None = None
    if capability.requires_gpu:
        from .rank_markers_gpu import inspect_runtime

        _, gpu_rsc, gpu_runtime = inspect_runtime(config.model_dump(mode="json"))
    groups_path = Path(str(row["groups_path"]))
    groups_frame = pd.read_csv(groups_path, sep="\t", keep_default_na=False)
    if groups_frame["barcode"].astype(str).duplicated().any():
        raise ContractError("candidate comparator group table duplicates barcodes")
    target = sorted(
        groups_frame.loc[
            groups_frame["analysis_group"].astype(str).eq("candidate"), "barcode"
        ].astype(str)
    )
    reference = sorted(
        groups_frame.loc[
            groups_frame["analysis_group"].astype(str).eq("reference"), "barcode"
        ].astype(str)
    )
    from .contracts import values_sha256

    if values_sha256(target) != str(row["target_barcode_sha256"]) or values_sha256(
        reference
    ) != str(row["comparator_barcode_sha256"]):
        raise ContractError("candidate comparator membership hash mismatch")
    loaded = load_input(config, backed=False)
    try:
        root_positions = {
            str(barcode): index
            for index, barcode in enumerate(loaded.root.obs_names.astype(str))
        }
        try:
            positions = np.asarray(
                [root_positions[str(barcode)] for barcode in groups_frame["barcode"]]
            )
        except KeyError as exc:
            raise ContractError(
                f"candidate comparator contains an unknown barcode: {exc}"
            ) from exc
        key = str(row["lane_key"])
        matrix, detection, features, assay, _ = _lane_materials(
            loaded, config, assay_id, representation_name, key=key
        )
        matrix = matrix[positions]
        detection = detection[positions] if detection is not None else None
        groups = groups_frame["analysis_group"].astype(str)
        if capability.requires_gpu:
            bottom_up = compute_gpu_bottom_up_markers(
                matrix,
                groups,
                features,
                detection=detection,
                rsc=gpu_rsc,
                marker_config=config.markers.model_dump(mode="json"),
            )
        else:
            bottom_up = compute_bottom_up_markers(matrix, groups, features, detection=detection)
        evidence = score_marker_programs(
            matrix,
            groups,
            features,
            programs,
            detection=detection,
            run_id=config.run.run_id,
            assay=assay_id,
            representation=representation_name,
            capability=capability,
        )
        if assay.kind == "spatial" and config.input.keys.spatial:
            coordinates, _, _ = resolve_obsm(
                loaded, config.input.keys.spatial, preferred_assay=assay_id
            )
            spatial = spatial_context_evidence(
                coordinates[positions], groups
            )
            evidence = evidence.merge(spatial, on="cluster_id", how="left")
        separation = _binary_separation(matrix, groups)
        lane_dir.mkdir(parents=True, exist_ok=True)
        from .provenance import atomic_publish

        atomic_publish(
            evidence_path,
            evidence.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
        )
        atomic_publish(
            bottom_up_path,
            bottom_up.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
        )
        publish_json(
            comparison_path,
            {
                "schema_version": 3,
                "artifact_kind": "cell_curator_candidate_comparison",
                "run_id": config.run.run_id,
                "lane_id": comparison_lane_id,
                "candidate_key": str(row["candidate_key"]),
                "comparator_role": str(row["comparator_role"]),
                "separation": separation,
                "groups_sha256": sha256(groups_path),
                "program_evidence_sha256": sha256(evidence_path),
                "bottom_up_sha256": sha256(bottom_up_path),
            },
        )
        receipt = make_receipt(
            run_id=config.run.run_id,
            lane_key=comparison_lane_id,
            capability=capability,
            input_path=Path(config.input.path),
            config_path=config_path,
            output_path=evidence_path,
            artifacts={
                "program_evidence": evidence_path,
                "bottom_up_markers": bottom_up_path,
                "comparison": comparison_path,
            },
            feature_ids=list(map(str, features)),
            groups=list(groups),
            limitations=capability.reasons,
            runtime=gpu_runtime,
            marker_methods=("wilcoxon", "logreg") if capability.requires_gpu else (),
            convergence_validated=True if capability.requires_gpu else None,
        )
        publish_json(receipt_path, receipt.model_dump(mode="json"))
        return LaneResult(comparison_lane_id, evidence, bottom_up, receipt_path)
    finally:
        close_input(loaded)


def summarize_program_cooccurrence(evidence: pd.DataFrame) -> pd.DataFrame:
    """Summarize cluster-level co-support without claiming cell-level detection."""

    columns = [
        "schema_version",
        "run_id",
        "cluster_id",
        "lane_key",
        "axis",
        "label_a",
        "label_b",
        "cooccurrence_scope",
        "co_support",
        "both_scores_positive",
        "maximum_negative_violation",
    ]
    if evidence.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, Any]] = []
    grouped = evidence.groupby(["cluster_id", "lane_key", "axis"], observed=True)
    for (cluster_id, evidence_lane, axis), block in grouped:
        by_label = (
            block.groupby("label", observed=True)
            .agg(
                score=("score", "mean"),
                positive_coverage=("positive_coverage", "mean"),
                negative_violation=("negative_violation", "max"),
            )
            .sort_index()
        )
        for label_a, label_b in combinations(map(str, by_label.index), 2):
            left = by_label.loc[label_a]
            right = by_label.loc[label_b]
            rows.append(
                {
                    "schema_version": 3,
                    "run_id": str(block.iloc[0]["run_id"]),
                    "cluster_id": str(cluster_id),
                    "lane_key": str(evidence_lane),
                    "axis": str(axis),
                    "label_a": label_a,
                    "label_b": label_b,
                    "cooccurrence_scope": "cluster_aggregate_program_support",
                    "co_support": min(
                        float(left["positive_coverage"]),
                        float(right["positive_coverage"]),
                    ),
                    "both_scores_positive": bool(
                        float(left["score"]) > 0 and float(right["score"]) > 0
                    ),
                    "maximum_negative_violation": max(
                        float(left["negative_violation"]),
                        float(right["negative_violation"]),
                    ),
                }
            )
    return pd.DataFrame(rows, columns=columns)


def summarize_source_corroboration(evidence: pd.DataFrame) -> pd.DataFrame:
    """Count independent knowledge sources and assay lanes supporting each program."""

    columns = [
        "schema_version",
        "run_id",
        "cluster_id",
        "label",
        "axis",
        "n_sources",
        "n_lanes",
        "n_annotation_capable_lanes",
        "mean_score",
        "mean_positive_coverage",
        "maximum_negative_violation",
        "score_sign_disagreement",
        "sources",
        "lanes",
    ]
    if evidence.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, Any]] = []
    for (cluster_id, label, axis), block in evidence.groupby(
        ["cluster_id", "label", "axis"], observed=True
    ):
        sources = sorted(value for value in set(block["source"].astype(str)) if value.strip())
        lanes = sorted(set(block["lane_key"].astype(str)))
        scores = pd.to_numeric(block["score"], errors="coerce").dropna()
        rows.append(
            {
                "schema_version": 3,
                "run_id": str(block.iloc[0]["run_id"]),
                "cluster_id": str(cluster_id),
                "label": str(label),
                "axis": str(axis),
                "n_sources": len(sources),
                "n_lanes": len(lanes),
                "n_annotation_capable_lanes": int(
                    block.loc[
                        block["capability_level"].astype(str).eq("annotation-capable"),
                        "lane_key",
                    ]
                    .astype(str)
                    .nunique()
                ),
                "mean_score": float(scores.mean()) if len(scores) else float("nan"),
                "mean_positive_coverage": float(block["positive_coverage"].mean()),
                "maximum_negative_violation": float(block["negative_violation"].max()),
                "score_sign_disagreement": bool(
                    len(scores) > 1 and (scores > 0).any() and (scores <= 0).any()
                ),
                "sources": ";".join(sources),
                "lanes": ";".join(lanes),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def consolidate_lanes(
    *,
    config: CellCuratorConfig,
    run_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    for lane in config.evidence.lanes:
        key = lane_key(lane.assay, lane.representation, lane.lane_key)
        path = run_root / "evidence" / "lanes" / key / "program_evidence.tsv"
        if path.is_file():
            frame = pd.read_csv(path, sep="\t", keep_default_na=False)
            frame["lane_key"] = key
            frames.append(frame)
        elif lane.required_for_acceptance:
            raise ContractError(f"required evidence lane is missing: {key}")
    evidence = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    loaded = load_input(config, backed=True)
    try:
        cluster_ids = sorted(
            loaded.root.obs[config.input.keys.canonical_parent].astype(str).unique()
        )
        scopes = loaded.root.obs[config.input.keys.scope].astype(str).unique()
        scope = str(scopes[0]) if len(scopes) == 1 else "multi_scope"
    finally:
        close_input(loaded)
    aggregate, decisions = aggregate_multimodal_evidence(
        evidence,
        build_capability_registry(config),
        config,
        scope=scope,
        cluster_ids=cluster_ids,
    )
    decision_frame = pd.DataFrame([decision.model_dump(mode="json") for decision in decisions])
    from .provenance import atomic_publish

    summary_dir = run_root / "evidence" / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    atomic_publish(
        summary_dir / "multimodal_label_scores.tsv",
        aggregate.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    atomic_publish(
        summary_dir / "annotation_proposals.tsv",
        decision_frame.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    state_evidence = (
        evidence[
            evidence.get("axis", pd.Series(index=evidence.index, dtype=str))
            .astype(str)
            .eq("cell_state")
        ]
        if not evidence.empty
        else pd.DataFrame(columns=["schema_version", "run_id", "cluster_id", "label", "axis"])
    )
    atomic_publish(
        summary_dir / "cell_state_evidence.tsv",
        state_evidence.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    cooccurrence = summarize_program_cooccurrence(evidence)
    corroboration = summarize_source_corroboration(evidence)
    cooccurrence_path = summary_dir / "program_cooccurrence.tsv"
    corroboration_path = summary_dir / "source_corroboration.tsv"
    atomic_publish(
        cooccurrence_path,
        cooccurrence.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    atomic_publish(
        corroboration_path,
        corroboration.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    manifest = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_multimodal_evidence",
        "run_id": config.run.run_id,
        "status": "complete",
        "n_lanes": len(frames),
        "n_clusters": len(cluster_ids),
        "n_unknown": int(decision_frame["cell_type"].eq("Unknown").sum()),
        "evidence_sha256": sha256(summary_dir / "multimodal_label_scores.tsv"),
        "proposals_sha256": sha256(summary_dir / "annotation_proposals.tsv"),
        "cell_state_evidence_sha256": sha256(summary_dir / "cell_state_evidence.tsv"),
        "program_cooccurrence_sha256": sha256(cooccurrence_path),
        "source_corroboration_sha256": sha256(corroboration_path),
    }
    publish_json(summary_dir / "evidence.manifest.json", manifest)
    return aggregate, decision_frame

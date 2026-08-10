"""Panel-aware marker assembly, profiling, rank scoring, and durable artifacts.

The functions in this module are deliberately label-neutral: they expose signed
marker evidence and disagreements but never turn a cross-check score into a
biological verdict.
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Any, Literal, Protocol, TypeAlias, cast

import numpy as np
import pandas as pd
from scipy import sparse

from .config import CellCuratorConfig, load_config
from .contracts import ContractError, safe_path_component, sha256
from .data import close_input, feature_ids_for, load_input, matrix_for
from .evidence import compute_bottom_up_markers
from .knowledge import (
    LocalMarkerProvider,
    MarkerProvider,
    collect_programs,
    providers_from_config,
)
from .models import MarkerProgram
from .provenance import atomic_publish, publish_json


class _ProgramSource(Protocol):
    def programs(self) -> list[MarkerProgram]: ...


@dataclass(frozen=True)
class _StaticProvider:
    values: tuple[MarkerProgram, ...]

    @property
    def identity(self) -> str:
        return "in-memory-marker-programs"

    def programs(self) -> list[MarkerProgram]:
        return list(self.values)


@dataclass(frozen=True)
class MarkerPanelProfile:
    """Top-down marker values and panel capability for every group and label."""

    mean_matrix: pd.DataFrame
    detection_matrix: pd.DataFrame
    feature_profile: pd.DataFrame
    panel_capability: pd.DataFrame
    program_summary: pd.DataFrame
    negative_violations: pd.DataFrame


@dataclass(frozen=True)
class DifferentialMarkerOverlap:
    """Bottom-up cluster markers and their overlap with signed label programs."""

    bottom_up: pd.DataFrame
    overlap: pd.DataFrame


@dataclass(frozen=True)
class MarkerAnalysis:
    """Complete in-memory marker analysis; all tables are ready for serialization."""

    programs: tuple[MarkerProgram, ...]
    profile: MarkerPanelProfile
    aucell: pd.DataFrame
    ucell: pd.DataFrame
    differential: DifferentialMarkerOverlap
    crosschecks: pd.DataFrame


@dataclass(frozen=True)
class MarkerRun:
    """A config-backed marker analysis and its durable artifact paths."""

    analysis: MarkerAnalysis
    artifacts: dict[str, Path]


MarkerSource: TypeAlias = MarkerProgram | MarkerProvider | _ProgramSource | str | Path


def _first_column(frame: pd.DataFrame, candidates: Sequence[str]) -> str | None:
    lowered = {str(column).lower(): str(column) for column in frame.columns}
    for candidate in candidates:
        if candidate in lowered:
            return lowered[candidate]
    return None


def _model_export_records(path: Path) -> pd.DataFrame | None:
    """Return long model-export records, or ``None`` for the standard provider format."""

    suffix = path.suffix.lower()
    if suffix == ".json":
        value = json.loads(path.read_text())
        if isinstance(value, dict) and isinstance(value.get("records"), list):
            frame = pd.DataFrame(value["records"])
        elif isinstance(value, list):
            frame = pd.DataFrame(value)
        elif (
            isinstance(value, dict)
            and value
            and all(
                isinstance(payload, dict)
                and payload
                and all(isinstance(weight, (int, float)) for weight in payload.values())
                for payload in value.values()
            )
        ):
            frame = pd.DataFrame(
                [
                    {"label": label, "feature": feature, "coefficient": coefficient}
                    for label, payload in value.items()
                    for feature, coefficient in payload.items()
                ]
            )
        else:
            return None
    elif suffix in {".tsv", ".csv"}:
        frame = pd.read_csv(path, sep="\t" if suffix == ".tsv" else ",")
    else:
        return None
    feature_column = _first_column(frame, ("feature", "gene", "gene_symbol", "marker"))
    label_column = _first_column(frame, ("label", "cell_type", "class"))
    wide_columns = {str(column).lower() for column in frame.columns}
    if feature_column is None or label_column is None:
        return None
    if wide_columns.intersection({"positive", "negative", "pos", "neg"}):
        return None
    return frame.rename(columns={feature_column: "feature", label_column: "label"})


def _consistent_value(rows: pd.DataFrame, columns: Sequence[str], default: str = "") -> str:
    column = _first_column(rows, columns)
    if column is None:
        return default
    values = tuple(dict.fromkeys(rows[column].fillna("").astype(str).str.strip()))
    nonblank = tuple(value for value in values if value)
    if len(nonblank) > 1:
        raise ContractError(f"model export has conflicting {column!r} values: {nonblank}")
    return nonblank[0] if nonblank else default


def _programs_from_model_export(
    path: Path,
    frame: pd.DataFrame,
    *,
    source_version: str,
    max_features_per_sign: int | None,
) -> list[MarkerProgram]:
    if frame.empty:
        return []
    direction_column = _first_column(frame, ("direction", "sign", "marker_sign"))
    weight_column = _first_column(frame, ("coefficient", "weight", "score", "importance"))
    rows = frame.copy()
    rows["label"] = rows["label"].fillna("").astype(str).str.strip()
    rows["feature"] = rows["feature"].fillna("").astype(str).str.strip()
    if rows["label"].eq("").any() or rows["feature"].eq("").any():
        raise ContractError(
            "model-export marker rows require nonblank label and feature values"
        )
    if weight_column is not None:
        rows["_weight"] = pd.to_numeric(rows[weight_column], errors="raise")
    else:
        rows["_weight"] = 1.0

    def marker_sign(row: pd.Series) -> str | None:
        if direction_column is None:
            if float(row["_weight"]) == 0:
                return None
            return "positive" if float(row["_weight"]) > 0 else "negative"
        direction = str(row[direction_column]).strip().lower()
        if direction in {"positive", "pos", "+", "1", "up"}:
            return "positive"
        if direction in {"negative", "neg", "-", "-1", "down"}:
            return "negative"
        if direction in {"", "0", "none", "nan"}:
            if float(row["_weight"]) == 0:
                return None
            return "positive" if float(row["_weight"]) > 0 else "negative"
        raise ContractError(f"invalid model-export marker direction: {direction!r}")

    rows["_sign"] = rows.apply(marker_sign, axis=1)
    rows = rows[rows["_sign"].notna()].copy()
    rows["_absolute_weight"] = rows["_weight"].abs()
    programs: list[MarkerProgram] = []
    for label, label_rows in rows.groupby("label", sort=True, observed=True):
        axis = _consistent_value(label_rows, ("axis",), "cell_type")
        if axis not in {"cell_type", "cell_state"}:
            raise ContractError(f"model export has invalid marker axis: {axis!r}")
        signed: dict[str, tuple[str, ...]] = {}
        for sign in ("positive", "negative"):
            selected = label_rows[label_rows["_sign"].eq(sign)].sort_values(
                ["_absolute_weight", "feature"], ascending=[False, True], kind="mergesort"
            )
            features = tuple(dict.fromkeys(selected["feature"].astype(str)))
            if max_features_per_sign is not None:
                features = features[:max_features_per_sign]
            signed[sign] = features
        programs.append(
            MarkerProgram(
                label=str(label),
                axis=cast(Literal["cell_type", "cell_state"], axis),
                positive=signed["positive"],
                negative=signed["negative"],
                parent=_consistent_value(label_rows, ("parent", "parent_label")) or None,
                parent_cl_id=_consistent_value(label_rows, ("parent_cl_id",)) or None,
                cl_id=_consistent_value(label_rows, ("cl_id", "cell_ontology_id")) or None,
                feature_space=_consistent_value(label_rows, ("feature_space",), "gene"),
                source=f"model-export:{path.resolve()}",
                source_version=source_version,
            )
        )
    return programs


def filter_marker_programs(
    programs: Iterable[MarkerProgram],
    *,
    parent_scope: str | None = None,
    parent_cl_id: str | None = None,
    include_unscoped: bool = False,
    axis: str | None = None,
    feature_space: str | None = None,
) -> list[MarkerProgram]:
    """Select programs valid within one parent without leaking sibling vocabularies."""

    selected: list[MarkerProgram] = []
    for program in programs:
        if axis is not None and program.axis != axis:
            continue
        if feature_space is not None and program.feature_space != feature_space:
            continue
        if parent_scope is not None and program.parent != parent_scope:
            if not (
                include_unscoped and program.parent is None and program.parent_cl_id is None
            ):
                continue
        if parent_cl_id is not None and program.parent_cl_id != parent_cl_id:
            if not (
                include_unscoped and program.parent is None and program.parent_cl_id is None
            ):
                continue
        selected.append(program)
    return sorted(selected, key=lambda item: (item.axis, item.parent or "", item.label))


def assemble_marker_programs_from_config(
    config_path: str | Path,
    *,
    source: Iterable[str | Path] = (),
    parent_scope: str | None = None,
    parent_cl_id: str | None = None,
    include_unscoped: bool = False,
    axis: str | None = None,
    feature_space: str | None = None,
    output_path: str | Path | None = None,
) -> tuple[list[MarkerProgram], Path]:
    """Assemble the signed marker vocabulary declared by a configuration.

    Combines the configured knowledge providers with any extra local sources and
    publishes the result under the run directory. Shared by the CLI and the
    Python API so the vocabulary is assembled one way only.
    """

    from .pipeline import run_root
    from .provenance import publish_json

    config = load_config(config_path)
    programs = assemble_marker_programs(
        [
            *providers_from_config(
                config.knowledge.providers,
                allow_remote_enhancement=config.knowledge.allow_optional_remote_enhancement,
            ),
            *source,
        ],
        parent_scope=parent_scope,
        parent_cl_id=parent_cl_id,
        include_unscoped=include_unscoped,
        axis=axis,
        feature_space=feature_space,
    )
    output = (
        Path(output_path)
        if output_path is not None
        else run_root(config) / "markers" / "assembled_programs.json"
    )
    publish_json(output, [item.model_dump(mode="json") for item in programs])
    return programs, output


def assemble_marker_programs(
    sources: Iterable[MarkerSource],
    *,
    parent_scope: str | None = None,
    parent_cl_id: str | None = None,
    include_unscoped: bool = False,
    axis: str | None = None,
    feature_space: str | None = None,
    source_version: str = "unversioned",
    max_features_per_sign: int | None = None,
) -> list[MarkerProgram]:
    """Assemble standard providers and signed long-form model exports into one vocabulary."""

    if max_features_per_sign is not None and max_features_per_sign < 1:
        raise ContractError("max_features_per_sign must be positive")
    providers: list[MarkerProvider] = []
    in_memory: list[MarkerProgram] = []
    for source in sources:
        if isinstance(source, MarkerProgram):
            in_memory.append(source)
            continue
        if isinstance(source, (str, Path)):
            path = Path(source)
            if not path.is_file():
                raise ContractError(f"marker source does not exist: {path}")
            model_export = _model_export_records(path)
            if model_export is None:
                providers.append(LocalMarkerProvider(path, version=source_version))
            else:
                providers.append(
                    _StaticProvider(
                        tuple(
                            _programs_from_model_export(
                                path,
                                model_export,
                                source_version=source_version,
                                max_features_per_sign=max_features_per_sign,
                            )
                        )
                    )
                )
            continue
        if not callable(getattr(source, "programs", None)):
            raise ContractError(f"unsupported marker source: {type(source).__name__}")
        providers.append(cast(MarkerProvider, source))
    if in_memory:
        providers.append(_StaticProvider(tuple(in_memory)))
    merged = collect_programs(providers)
    for program in merged:
        contradiction = sorted(set(program.positive).intersection(program.negative))
        if contradiction:
            raise ContractError(
                f"marker program {program.label!r} has contradictory signed features: {contradiction}"
            )
    return filter_marker_programs(
        merged,
        parent_scope=parent_scope,
        parent_cl_id=parent_cl_id,
        include_unscoped=include_unscoped,
        axis=axis,
        feature_space=feature_space,
    )


def _validate_inputs(
    matrix: Any,
    groups: pd.Series,
    features: pd.Index,
    detection: Any | None,
) -> None:
    if matrix.ndim != 2 or matrix.shape != (len(groups), len(features)):
        raise ContractError("matrix, groups, and feature dimensions disagree")
    if len(groups) == 0:
        raise ContractError("marker profiling requires at least one cell")
    if len(features) == 0:
        raise ContractError("marker profiling requires at least one feature")
    if features.has_duplicates:
        raise ContractError("marker profiling requires unique feature identifiers")
    if groups.isna().any():
        raise ContractError("marker profiling groups contain missing values")
    if detection is not None:
        if detection.shape != matrix.shape:
            raise ContractError("detection matrix shape differs from the expression matrix")
        values = detection.data if sparse.issparse(detection) else np.asarray(detection)
        if values.size and (not np.isfinite(values).all() or float(values.min()) < 0):
            raise ContractError("detection matrix must be finite and nonnegative")


def _mean_rows(matrix: Any, mask: np.ndarray) -> np.ndarray:
    result = matrix[mask].mean(axis=0)
    return np.asarray(result).ravel()


def _fraction_rows(matrix: Any, mask: np.ndarray, threshold: float) -> np.ndarray:
    selected = matrix[mask]
    if sparse.issparse(selected):
        result = (selected > threshold).mean(axis=0)
    else:
        result = (np.asarray(selected) > threshold).mean(axis=0)
    return np.asarray(result).ravel()


def _panel_capability(
    programs: Sequence[MarkerProgram],
    positions: dict[str, int],
    *,
    min_measured_positive: int,
    min_positive_fraction: float,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for program in programs:
        measured_positive = tuple(
            feature for feature in program.positive if feature in positions
        )
        measured_negative = tuple(
            feature for feature in program.negative if feature in positions
        )
        missing_positive = tuple(
            feature for feature in program.positive if feature not in positions
        )
        missing_negative = tuple(
            feature for feature in program.negative if feature not in positions
        )
        fraction = len(measured_positive) / len(program.positive) if program.positive else 0.0
        reasons: list[str] = []
        if not program.positive:
            reasons.append("no_positive_markers_declared")
        if len(measured_positive) < min_measured_positive:
            reasons.append("insufficient_measured_positive_markers")
        if fraction < min_positive_fraction:
            reasons.append("positive_panel_fraction_below_threshold")
        rows.append(
            {
                "schema_version": 3,
                "label": program.label,
                "axis": program.axis,
                "parent_label": program.parent or "",
                "parent_cl_id": program.parent_cl_id or "",
                "feature_space": program.feature_space,
                "total_positive": len(program.positive),
                "measured_positive": len(measured_positive),
                "positive_panel_fraction": fraction,
                "total_negative": len(program.negative),
                "measured_negative": len(measured_negative),
                "negative_panel_fraction": (
                    len(measured_negative) / len(program.negative) if program.negative else 1.0
                ),
                "missing_positive": ";".join(missing_positive),
                "missing_negative": ";".join(missing_negative),
                "callable": not reasons,
                "uncallable_reasons": ";".join(reasons),
                "source": program.source,
            }
        )
    columns = [
        "schema_version",
        "label",
        "axis",
        "parent_label",
        "parent_cl_id",
        "feature_space",
        "total_positive",
        "measured_positive",
        "positive_panel_fraction",
        "total_negative",
        "measured_negative",
        "negative_panel_fraction",
        "missing_positive",
        "missing_negative",
        "callable",
        "uncallable_reasons",
        "source",
    ]
    return pd.DataFrame(rows, columns=columns)


def negative_marker_violations(
    feature_profile: pd.DataFrame,
    *,
    detection_threshold: float,
) -> pd.DataFrame:
    """Classify each expected-negative marker without treating missing values as absent."""

    columns = [
        "schema_version",
        "cluster_id",
        "label",
        "axis",
        "parent_label",
        "feature",
        "measured",
        "mean_value",
        "fraction_detected",
        "detection_threshold",
        "status",
        "violated",
    ]
    negative = feature_profile[feature_profile["sign"].eq("negative")].copy()
    if negative.empty:
        return pd.DataFrame(columns=columns)
    statuses: list[str] = []
    violated: list[bool | None] = []
    for row in negative.itertuples(index=False):
        if not bool(row.measured):
            statuses.append("unmeasured")
            violated.append(None)
        elif pd.isna(row.fraction_detected):
            statuses.append("not_assessed")
            violated.append(None)
        elif float(row.fraction_detected) >= detection_threshold:
            statuses.append("violated")
            violated.append(True)
        else:
            statuses.append("not_violated")
            violated.append(False)
    result = negative[
        [
            "schema_version",
            "cluster_id",
            "label",
            "axis",
            "parent_label",
            "feature",
            "measured",
            "mean_value",
            "fraction_detected",
        ]
    ].copy()
    result["detection_threshold"] = detection_threshold
    result["status"] = statuses
    result["violated"] = pd.array(violated, dtype="boolean")
    return result[columns]


def profile_marker_panel(
    matrix: Any,
    groups: pd.Series,
    features: pd.Index,
    programs: Sequence[MarkerProgram],
    *,
    detection: Any | None,
    detection_threshold: float = 0.1,
    min_measured_positive: int = 1,
    min_positive_fraction: float = 0.5,
) -> MarkerPanelProfile:
    """Build top-down mean/detection matrices and signed program summaries."""

    if not 0 <= min_positive_fraction <= 1:
        raise ContractError("min_positive_fraction must lie in [0, 1]")
    if not 0 <= detection_threshold <= 1:
        raise ContractError("detection_threshold must lie in [0, 1]")
    if min_measured_positive < 1:
        raise ContractError("min_measured_positive must be positive")
    _validate_inputs(matrix, groups, features, detection)
    cluster_ids = sorted(groups.astype(str).unique())
    group_values = groups.astype(str).to_numpy()
    positions = {str(feature): index for index, feature in enumerate(features)}
    program_features = tuple(
        dict.fromkeys(
            feature
            for program in programs
            for feature in (*program.positive, *program.negative)
        )
    )
    means = pd.DataFrame(np.nan, index=cluster_ids, columns=program_features, dtype=float)
    fractions = pd.DataFrame(np.nan, index=cluster_ids, columns=program_features, dtype=float)
    all_means: dict[str, np.ndarray] = {}
    all_fractions: dict[str, np.ndarray] = {}
    for cluster_id in cluster_ids:
        mask = group_values == cluster_id
        all_means[cluster_id] = _mean_rows(matrix, mask)
        if detection is not None:
            all_fractions[cluster_id] = _fraction_rows(detection, mask, 0.0)
        for feature in program_features:
            position = positions.get(feature)
            if position is None:
                continue
            means.loc[cluster_id, feature] = all_means[cluster_id][position]
            if detection is not None:
                fractions.loc[cluster_id, feature] = all_fractions[cluster_id][position]
    means.index.name = "cluster_id"
    fractions.index.name = "cluster_id"
    panel = _panel_capability(
        programs,
        positions,
        min_measured_positive=min_measured_positive,
        min_positive_fraction=min_positive_fraction,
    )
    profile_rows: list[dict[str, Any]] = []
    for program in programs:
        for sign, marker_features in (
            ("positive", program.positive),
            ("negative", program.negative),
        ):
            for feature in marker_features:
                measured = feature in positions
                for cluster_id in cluster_ids:
                    profile_rows.append(
                        {
                            "schema_version": 3,
                            "cluster_id": cluster_id,
                            "label": program.label,
                            "axis": program.axis,
                            "parent_label": program.parent or "",
                            "feature_space": program.feature_space,
                            "feature": feature,
                            "sign": sign,
                            "measured": measured,
                            "mean_value": (
                                float(means.loc[cluster_id, feature])
                                if measured
                                else float("nan")
                            ),
                            "fraction_detected": (
                                float(fractions.loc[cluster_id, feature])
                                if measured and detection is not None
                                else float("nan")
                            ),
                            "source": program.source,
                        }
                    )
    feature_columns = [
        "schema_version",
        "cluster_id",
        "label",
        "axis",
        "parent_label",
        "feature_space",
        "feature",
        "sign",
        "measured",
        "mean_value",
        "fraction_detected",
        "source",
    ]
    feature_profile = pd.DataFrame(profile_rows, columns=feature_columns)
    violations = negative_marker_violations(
        feature_profile,
        detection_threshold=detection_threshold,
    )
    standardized = means.copy()
    for feature in standardized:
        values = standardized[feature].to_numpy(dtype=float)
        finite = np.isfinite(values)
        if not finite.any():
            continue
        spread = float(np.nanstd(values))
        standardized[feature] = (values - float(np.nanmean(values))) / (
            spread if spread > 1e-12 else 1.0
        )
    summary_rows: list[dict[str, Any]] = []
    for program_index, program in enumerate(programs):
        capability = panel.iloc[program_index]
        positive = [feature for feature in program.positive if feature in positions]
        negative = [feature for feature in program.negative if feature in positions]
        for cluster_id in cluster_ids:
            positive_signal = (
                float(standardized.loc[cluster_id, positive].mean())
                if positive
                else float("nan")
            )
            negative_signal = (
                float(np.maximum(standardized.loc[cluster_id, negative].to_numpy(), 0).mean())
                if negative
                else 0.0
            )
            positive_detection = (
                float(fractions.loc[cluster_id, positive].mean())
                if detection is not None and positive
                else float("nan")
            )
            negative_detection = (
                float(fractions.loc[cluster_id, negative].mean())
                if detection is not None and negative
                else (0.0 if detection is not None else float("nan"))
            )
            negative_violation = (
                float((fractions.loc[cluster_id, negative] >= detection_threshold).mean())
                if detection is not None and negative
                else (0.0 if detection is not None else float("nan"))
            )
            callable_program = bool(capability["callable"])
            score = (
                positive_signal
                - negative_signal
                - (negative_violation if detection is not None else 0.0)
                if callable_program
                else float("nan")
            )
            summary_rows.append(
                {
                    "schema_version": 3,
                    "cluster_id": cluster_id,
                    "label": program.label,
                    "axis": program.axis,
                    "parent_label": program.parent or "",
                    "feature_space": program.feature_space,
                    "callable": callable_program,
                    "uncallable_reasons": capability["uncallable_reasons"],
                    "top_down_score": score,
                    "positive_signal": positive_signal,
                    "negative_signal": negative_signal,
                    "positive_detection_mean": positive_detection,
                    "negative_detection_mean": negative_detection,
                    "negative_violation": negative_violation,
                    "detection_valid": detection is not None,
                    "measured_positive": int(capability["measured_positive"]),
                    "total_positive": int(capability["total_positive"]),
                    "measured_negative": int(capability["measured_negative"]),
                    "total_negative": int(capability["total_negative"]),
                    "missing_positive": capability["missing_positive"],
                    "missing_negative": capability["missing_negative"],
                }
            )
    return MarkerPanelProfile(
        mean_matrix=means,
        detection_matrix=fractions,
        feature_profile=feature_profile,
        panel_capability=panel.reset_index(drop=True),
        program_summary=pd.DataFrame(summary_rows),
        negative_violations=violations,
    )


def _cell_ids(count: int, cell_ids: Sequence[str] | None) -> tuple[str, ...]:
    if cell_ids is None:
        return tuple(str(index) for index in range(count))
    values = tuple(map(str, cell_ids))
    if len(values) != count or len(set(values)) != count:
        raise ContractError("cell IDs must be unique and match the matrix row count")
    return values


def _row_ranks(matrix: Any, features: pd.Index) -> Iterable[np.ndarray]:
    tie_breaker = features.astype(str).to_numpy()
    for index in range(matrix.shape[0]):
        row = (
            matrix.getrow(index).toarray().ravel()
            if sparse.issparse(matrix)
            else np.asarray(matrix[index]).ravel()
        )
        if not np.isfinite(row).all():
            raise ContractError("rank scoring requires finite matrix values")
        order = np.lexsort((tie_breaker, -row))
        ranks: np.ndarray = np.empty(len(order), dtype=int)
        ranks[order] = np.arange(1, len(order) + 1)
        yield ranks


def _signature_auc(ranks: np.ndarray, indices: Sequence[int], top_k: int) -> float:
    if not indices:
        return 0.0
    hits = sorted(int(ranks[index]) for index in indices if int(ranks[index]) <= top_k)
    if not hits:
        return 0.0
    area = sum(top_k - rank + 1 for rank in hits)
    maximum_hits = min(len(indices), top_k)
    maximum = sum(top_k - rank + 1 for rank in range(1, maximum_hits + 1))
    return float(area / maximum) if maximum else 0.0


def score_aucell(
    matrix: Any,
    features: pd.Index,
    programs: Sequence[MarkerProgram],
    *,
    cell_ids: Sequence[str] | None = None,
    max_rank_fraction: float = 0.05,
    min_measured_positive: int = 1,
) -> pd.DataFrame:
    """Compute deterministic AUCell-style recovery AUCs from per-cell feature ranks."""

    if matrix.ndim != 2 or matrix.shape[1] != len(features):
        raise ContractError("AUCell matrix and feature dimensions disagree")
    if not 0 < max_rank_fraction <= 1:
        raise ContractError("max_rank_fraction must lie in (0, 1]")
    if min_measured_positive < 1:
        raise ContractError("min_measured_positive must be positive")
    if len(features) == 0:
        raise ContractError("AUCell scoring requires at least one feature")
    positions = {str(feature): index for index, feature in enumerate(features)}
    identifiers = _cell_ids(matrix.shape[0], cell_ids)
    top_k = max(1, min(len(features), math.ceil(len(features) * max_rank_fraction)))
    rows: list[dict[str, Any]] = []
    for cell_id, ranks in zip(identifiers, _row_ranks(matrix, features), strict=True):
        for program in programs:
            positive = [positions[item] for item in program.positive if item in positions]
            negative = [positions[item] for item in program.negative if item in positions]
            callable_program = len(positive) >= min_measured_positive
            positive_auc = (
                _signature_auc(ranks, positive, top_k) if callable_program else float("nan")
            )
            negative_auc = _signature_auc(ranks, negative, top_k) if negative else 0.0
            rows.append(
                {
                    "schema_version": 3,
                    "cell_id": cell_id,
                    "label": program.label,
                    "axis": program.axis,
                    "parent_label": program.parent or "",
                    "feature_space": program.feature_space,
                    "callable": callable_program,
                    "measured_positive": len(positive),
                    "measured_negative": len(negative),
                    "top_rank": top_k,
                    "positive_aucell": positive_auc,
                    "negative_aucell": negative_auc,
                    "aucell_score": positive_auc - negative_auc
                    if callable_program
                    else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def _ucell_signature(ranks: np.ndarray, indices: Sequence[int], max_rank: int) -> float:
    if not indices:
        return 0.0
    capped = np.minimum(ranks[list(indices)], max_rank + 1)
    mann_whitney_u = float(capped.sum() - len(indices) * (len(indices) + 1) / 2)
    return float(np.clip(1.0 - mann_whitney_u / (len(indices) * max_rank), 0.0, 1.0))


def score_ucell(
    matrix: Any,
    features: pd.Index,
    programs: Sequence[MarkerProgram],
    *,
    cell_ids: Sequence[str] | None = None,
    max_rank: int = 1500,
    min_measured_positive: int = 1,
) -> pd.DataFrame:
    """Compute UCell-style Mann-Whitney rank scores with signed-marker penalties."""

    if matrix.ndim != 2 or matrix.shape[1] != len(features):
        raise ContractError("UCell matrix and feature dimensions disagree")
    if max_rank < 1:
        raise ContractError("max_rank must be positive")
    if min_measured_positive < 1:
        raise ContractError("min_measured_positive must be positive")
    if len(features) == 0:
        raise ContractError("UCell scoring requires at least one feature")
    effective_max_rank = min(max_rank, len(features))
    positions = {str(feature): index for index, feature in enumerate(features)}
    identifiers = _cell_ids(matrix.shape[0], cell_ids)
    rows: list[dict[str, Any]] = []
    for cell_id, ranks in zip(identifiers, _row_ranks(matrix, features), strict=True):
        for program in programs:
            positive = [positions[item] for item in program.positive if item in positions]
            negative = [positions[item] for item in program.negative if item in positions]
            callable_program = len(positive) >= min_measured_positive
            positive_score = (
                _ucell_signature(ranks, positive, effective_max_rank)
                if callable_program
                else float("nan")
            )
            negative_score = (
                _ucell_signature(ranks, negative, effective_max_rank) if negative else 0.0
            )
            rows.append(
                {
                    "schema_version": 3,
                    "cell_id": cell_id,
                    "label": program.label,
                    "axis": program.axis,
                    "parent_label": program.parent or "",
                    "feature_space": program.feature_space,
                    "callable": callable_program,
                    "measured_positive": len(positive),
                    "measured_negative": len(negative),
                    "max_rank": effective_max_rank,
                    "positive_ucell": positive_score,
                    "negative_ucell": negative_score,
                    "ucell_score": positive_score - negative_score
                    if callable_program
                    else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def compute_de_marker_overlap(
    matrix: Any,
    groups: pd.Series,
    features: pd.Index,
    programs: Sequence[MarkerProgram],
    *,
    detection: Any | None,
    top_n: int = 100,
    adjusted_pvalue_max: float = 0.05,
    min_effect: float = 0.0,
) -> DifferentialMarkerOverlap:
    """Compare label programs to one-cluster-vs-remainder markers, not condition-level DE."""

    if top_n < 1:
        raise ContractError("DE marker top_n must be positive")
    if not 0 <= adjusted_pvalue_max <= 1:
        raise ContractError("adjusted_pvalue_max must lie in [0, 1]")
    bottom_up = compute_bottom_up_markers(
        matrix,
        groups,
        features,
        detection=detection,
        top_n=top_n,
    )
    measured = set(features.astype(str))
    rows: list[dict[str, Any]] = []
    for cluster_id in sorted(groups.astype(str).unique()):
        cluster_markers = bottom_up[
            bottom_up["cluster_id"].astype(str).eq(cluster_id)
            & bottom_up["effect"].gt(min_effect)
            & bottom_up["adjusted_pvalue"].le(adjusted_pvalue_max)
        ].sort_values(["rank", "feature"], kind="mergesort")
        marker_set = set(cluster_markers["feature"].astype(str))
        marker_order = tuple(cluster_markers["feature"].astype(str))
        for program in programs:
            positive = tuple(feature for feature in program.positive if feature in measured)
            negative = tuple(feature for feature in program.negative if feature in measured)
            positive_overlap = tuple(feature for feature in marker_order if feature in positive)
            negative_overlap = tuple(feature for feature in marker_order if feature in negative)
            positive_fraction = (
                len(positive_overlap) / len(positive) if positive else float("nan")
            )
            negative_fraction = len(negative_overlap) / len(negative) if negative else 0.0
            rows.append(
                {
                    "schema_version": 3,
                    "cluster_id": cluster_id,
                    "label": program.label,
                    "axis": program.axis,
                    "parent_label": program.parent or "",
                    "feature_space": program.feature_space,
                    "n_de_markers": len(marker_set),
                    "measured_positive": len(positive),
                    "positive_overlap_n": len(positive_overlap),
                    "positive_overlap_fraction": positive_fraction,
                    "positive_overlap": ";".join(positive_overlap),
                    "measured_negative": len(negative),
                    "negative_overlap_n": len(negative_overlap),
                    "negative_overlap_fraction": negative_fraction,
                    "negative_overlap": ";".join(negative_overlap),
                    "de_overlap_score": (
                        positive_fraction - negative_fraction if positive else float("nan")
                    ),
                }
            )
    return DifferentialMarkerOverlap(bottom_up=bottom_up, overlap=pd.DataFrame(rows))


def _cell_score_means(
    scores: pd.DataFrame,
    groups: pd.Series,
    score_column: str,
) -> pd.DataFrame:
    membership = pd.DataFrame(
        {"cell_id": groups.index.astype(str), "cluster_id": groups.astype(str).to_numpy()}
    )
    if membership["cell_id"].duplicated().any():
        raise ContractError("cell-group membership contains duplicate cell IDs")
    merged = scores.merge(membership, on="cell_id", how="left", validate="many_to_one")
    if merged["cluster_id"].isna().any():
        raise ContractError("rank-score cell IDs do not match group membership")
    return (
        merged.groupby(
            ["cluster_id", "label", "axis", "parent_label", "feature_space"],
            observed=True,
        )[score_column]
        .mean()
        .rename(f"{score_column}_mean")
        .reset_index()
    )


def _unit_score(value: object, *, logistic: bool = False) -> float:
    if pd.isna(value):
        return float("nan")
    number = float(cast(Any, value))
    if logistic:
        return float(1.0 / (1.0 + math.exp(-float(np.clip(number, -50, 50)))))
    return float(np.clip((number + 1.0) / 2.0, 0.0, 1.0))


def consolidate_marker_crosschecks(
    profile: MarkerPanelProfile,
    aucell: pd.DataFrame,
    ucell: pd.DataFrame,
    differential_overlap: pd.DataFrame,
    *,
    groups: pd.Series,
) -> pd.DataFrame:
    """Consolidate independent marker routes while keeping every component visible."""

    keys = ["cluster_id", "label", "axis", "parent_label", "feature_space"]
    result = profile.program_summary.copy()
    result = result.merge(
        _cell_score_means(aucell, groups, "aucell_score"), on=keys, how="left"
    )
    result = result.merge(_cell_score_means(ucell, groups, "ucell_score"), on=keys, how="left")
    result = result.merge(
        differential_overlap[keys + ["de_overlap_score"]],
        on=keys,
        how="left",
        validate="one_to_one",
    )
    normalized_columns = [
        "top_down_component",
        "aucell_component",
        "ucell_component",
        "de_overlap_component",
    ]
    result["top_down_component"] = result["top_down_score"].map(
        lambda value: _unit_score(value, logistic=True)
    )
    result["aucell_component"] = result["aucell_score_mean"].map(_unit_score)
    result["ucell_component"] = result["ucell_score_mean"].map(_unit_score)
    result["de_overlap_component"] = result["de_overlap_score"].map(_unit_score)
    consolidated: list[float] = []
    spreads: list[float] = []
    routes: list[str] = []
    for row in result.itertuples(index=False):
        values = {
            "top_down": row.top_down_component,
            "aucell": row.aucell_component,
            "ucell": row.ucell_component,
            "de_overlap": row.de_overlap_component,
        }
        available = {name: float(value) for name, value in values.items() if pd.notna(value)}
        callable_program = bool(row.callable)
        consolidated.append(
            float(np.mean(list(available.values())))
            if available and callable_program
            else float("nan")
        )
        spreads.append(
            float(max(available.values()) - min(available.values()))
            if available
            else float("nan")
        )
        routes.append(";".join(available))
    result["crosscheck_score"] = consolidated
    result["crosscheck_spread"] = spreads
    result["routes_available"] = routes
    result["n_routes"] = result["routes_available"].map(
        lambda value: len(value.split(";")) if value else 0
    )
    ordered = [
        "schema_version",
        *keys,
        "callable",
        "uncallable_reasons",
        "crosscheck_score",
        "crosscheck_spread",
        "n_routes",
        "routes_available",
        "top_down_score",
        "aucell_score_mean",
        "ucell_score_mean",
        "de_overlap_score",
        *normalized_columns,
        "positive_signal",
        "negative_signal",
        "positive_detection_mean",
        "negative_detection_mean",
        "negative_violation",
        "detection_valid",
        "measured_positive",
        "total_positive",
        "measured_negative",
        "total_negative",
        "missing_positive",
        "missing_negative",
    ]
    return result[ordered].sort_values(
        ["cluster_id", "crosscheck_score", "label"],
        ascending=[True, False, True],
        na_position="last",
        kind="mergesort",
    )


def run_marker_analysis(
    matrix: Any,
    groups: pd.Series,
    features: pd.Index,
    programs: Sequence[MarkerProgram],
    *,
    detection: Any | None,
    detection_threshold: float = 0.1,
    min_measured_positive: int = 1,
    min_positive_fraction: float = 0.5,
    aucell_max_rank_fraction: float = 0.05,
    ucell_max_rank: int = 1500,
    de_top_n: int = 100,
    adjusted_pvalue_max: float = 0.05,
    min_de_effect: float = 0.0,
) -> MarkerAnalysis:
    """Run all marker routes once from shared matrices and a shared signed vocabulary."""

    if not programs:
        raise ContractError("marker analysis requires at least one marker program")
    _validate_inputs(matrix, groups, features, detection)
    identifiers = groups.index.astype(str).tolist()
    profile = profile_marker_panel(
        matrix,
        groups,
        features,
        programs,
        detection=detection,
        detection_threshold=detection_threshold,
        min_measured_positive=min_measured_positive,
        min_positive_fraction=min_positive_fraction,
    )
    aucell = score_aucell(
        matrix,
        features,
        programs,
        cell_ids=identifiers,
        max_rank_fraction=aucell_max_rank_fraction,
        min_measured_positive=min_measured_positive,
    )
    ucell = score_ucell(
        matrix,
        features,
        programs,
        cell_ids=identifiers,
        max_rank=ucell_max_rank,
        min_measured_positive=min_measured_positive,
    )
    differential = compute_de_marker_overlap(
        matrix,
        groups,
        features,
        programs,
        detection=detection,
        top_n=de_top_n,
        adjusted_pvalue_max=adjusted_pvalue_max,
        min_effect=min_de_effect,
    )
    crosschecks = consolidate_marker_crosschecks(
        profile,
        aucell,
        ucell,
        differential.overlap,
        groups=groups,
    )
    return MarkerAnalysis(
        programs=tuple(programs),
        profile=profile,
        aucell=aucell,
        ucell=ucell,
        differential=differential,
        crosschecks=crosschecks,
    )


def _interpolate(
    start: tuple[int, int, int], end: tuple[int, int, int], fraction: float
) -> str:
    values = tuple(
        round(left + (right - left) * fraction) for left, right in zip(start, end, strict=True)
    )
    return "#" + "".join(f"{value:02x}" for value in values)


def _dot_color(value: float) -> str:
    value = float(np.clip(value, 0, 1))
    if value <= 0.5:
        return _interpolate((33, 102, 172), (247, 247, 247), value * 2)
    return _interpolate((247, 247, 247), (178, 24, 43), (value - 0.5) * 2)


def marker_dotplot_svg(profile: MarkerPanelProfile) -> str:
    """Render a dependency-free, self-contained SVG from actual mean/detection values."""

    frame = profile.feature_profile.copy()
    if frame.empty:
        return (
            '<svg xmlns="http://www.w3.org/2000/svg" width="480" height="120" '
            'viewBox="0 0 480 120"><title>Marker dotplot</title>'
            '<text x="20" y="60">No marker features were declared.</text></svg>\n'
        )
    columns = tuple(
        dict.fromkeys(
            (str(row.label), str(row.feature_space), str(row.feature), str(row.sign))
            for row in frame.itertuples(index=False)
        )
    )
    clusters = tuple(sorted(frame["cluster_id"].astype(str).unique()))
    left, top, cell_width, cell_height = 150, 150, 48, 42
    width = max(540, left + cell_width * len(columns) + 130)
    height = max(300, top + cell_height * len(clusters) + 100)
    measured = frame[frame["measured"] & frame["mean_value"].notna()]["mean_value"].astype(
        float
    )
    low = float(measured.min()) if not measured.empty else 0.0
    high = float(measured.max()) if not measured.empty else 1.0

    def scaled_mean(value: float) -> float:
        return 0.5 if high - low <= 1e-12 else (value - low) / (high - low)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<title>cell-curator signed marker dotplot</title>",
        "<metadata>Circle area encodes detection fraction; color encodes observed group mean. "
        "Crosses are unmeasured features.</metadata>",
        "<style>text{font-family:system-ui,sans-serif;font-size:11px}.axis{fill:#222}.grid{stroke:#ddd;stroke-width:1}.negative{stroke:#a40000;stroke-width:2}.positive{stroke:#333;stroke-width:1}</style>",
        '<text x="20" y="24" font-size="16" font-weight="bold">Signed marker profile</text>',
    ]
    for column_index, (label, feature_space, feature, sign) in enumerate(columns):
        x = left + column_index * cell_width + cell_width / 2
        text = escape(
            f"{label}: {feature} [{feature_space}] ({'+' if sign == 'positive' else '−'})"
        )
        parts.append(
            f'<text class="axis" x="{x:.1f}" y="{top - 12}" transform="rotate(-55 {x:.1f} {top - 12})">{text}</text>'
        )
    lookup = {
        (
            str(row.cluster_id),
            str(row.label),
            str(row.feature_space),
            str(row.feature),
            str(row.sign),
        ): row
        for row in frame.itertuples(index=False)
    }
    for row_index, cluster_id in enumerate(clusters):
        y = top + row_index * cell_height + cell_height / 2
        parts.append(
            f'<text class="axis" x="{left - 12}" y="{y + 4:.1f}" text-anchor="end">{escape(cluster_id)}</text>'
        )
        for column_index, (label, feature_space, feature, sign) in enumerate(columns):
            x = left + column_index * cell_width + cell_width / 2
            parts.append(
                f'<line class="grid" x1="{x:.1f}" y1="{y - cell_height / 2:.1f}" x2="{x:.1f}" y2="{y + cell_height / 2:.1f}"/>'
            )
            row = lookup[(cluster_id, label, feature_space, feature, sign)]
            if not bool(row.measured):
                parts.append(
                    f'<path d="M{x - 4:.1f},{y - 4:.1f} L{x + 4:.1f},{y + 4:.1f} M{x + 4:.1f},{y - 4:.1f} L{x - 4:.1f},{y + 4:.1f}" stroke="#777"><title>{escape(feature)}: unmeasured</title></path>'
                )
                continue
            fraction = float(row.fraction_detected) if pd.notna(row.fraction_detected) else 0.0
            mean = float(row.mean_value)
            radius = 2.5 + 10.0 * math.sqrt(max(0.0, fraction))
            css_class = "positive" if sign == "positive" else "negative"
            title = escape(
                f"cluster={cluster_id}; label={label}; feature={feature}; sign={sign}; mean={mean:.6g}; detection={fraction:.6g}"
            )
            parts.append(
                f'<circle class="{css_class}" cx="{x:.1f}" cy="{y:.1f}" r="{radius:.2f}" '
                f'fill="{_dot_color(scaled_mean(mean))}"><title>{title}</title></circle>'
            )
    legend_y = height - 35
    parts.extend(
        [
            f'<text x="20" y="{legend_y}">Detection:</text>',
            f'<circle cx="95" cy="{legend_y - 4}" r="5" fill="#999"><title>25% detected</title></circle>',
            f'<circle cx="125" cy="{legend_y - 4}" r="10" fill="#999"><title>100% detected</title></circle>',
            f'<text x="150" y="{legend_y}">Color: low mean</text>',
            f'<rect x="240" y="{legend_y - 13}" width="22" height="14" fill="#2166ac"/>',
            f'<rect x="262" y="{legend_y - 13}" width="22" height="14" fill="#f7f7f7"/>',
            f'<rect x="284" y="{legend_y - 13}" width="22" height="14" fill="#b2182b"/>',
            f'<text x="312" y="{legend_y}">high mean</text>',
            "</svg>",
        ]
    )
    return "\n".join(parts) + "\n"


def _tsv(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(sep="\t", index=False, lineterminator="\n").encode()


def write_marker_artifacts(
    analysis: MarkerAnalysis,
    output_dir: str | Path,
    *,
    run_id: str = "standalone",
) -> dict[str, Path]:
    """Atomically publish all marker tables, JSON provenance, and the real SVG dotplot."""

    root = Path(output_dir)
    paths = {
        "programs": root / "marker_programs.json",
        "feature_profile": root / "marker_feature_profile.tsv",
        "mean_matrix": root / "marker_mean_matrix.tsv",
        "detection_matrix": root / "marker_detection_matrix.tsv",
        "panel_capability": root / "panel_capability.tsv",
        "program_summary": root / "top_down_program_summary.tsv",
        "negative_violations": root / "negative_marker_violations.tsv",
        "aucell": root / "aucell_scores.tsv",
        "ucell": root / "ucell_scores.tsv",
        "bottom_up": root / "bottom_up_markers.tsv",
        "de_overlap": root / "de_marker_overlap.tsv",
        "crosschecks": root / "consolidated_crosschecks.tsv",
        "dotplot": root / "marker_dotplot.svg",
    }
    publish_json(
        paths["programs"],
        [program.model_dump(mode="json") for program in analysis.programs],
    )
    frames = {
        "feature_profile": analysis.profile.feature_profile,
        "mean_matrix": analysis.profile.mean_matrix.reset_index(),
        "detection_matrix": analysis.profile.detection_matrix.reset_index(),
        "panel_capability": analysis.profile.panel_capability,
        "program_summary": analysis.profile.program_summary,
        "negative_violations": analysis.profile.negative_violations,
        "aucell": analysis.aucell,
        "ucell": analysis.ucell,
        "bottom_up": analysis.differential.bottom_up,
        "de_overlap": analysis.differential.overlap,
        "crosschecks": analysis.crosschecks,
    }
    for name, frame in frames.items():
        atomic_publish(paths[name], _tsv(frame))
    atomic_publish(paths["dotplot"], marker_dotplot_svg(analysis.profile).encode())
    manifest_path = root / "marker_artifacts.manifest.json"
    publish_json(
        manifest_path,
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_marker_analysis",
            "run_id": run_id,
            "n_programs": len(analysis.programs),
            "n_crosschecks": len(analysis.crosschecks),
            "artifacts": {
                name: {"path": str(path), "sha256": sha256(path)}
                for name, path in sorted(paths.items())
            },
        },
    )
    paths["manifest"] = manifest_path
    return paths


def _obs_column(loaded: Any, adata: Any, key: str) -> pd.Series:
    if key in adata.obs:
        values = adata.obs[key].copy()
    elif key in loaded.root.obs:
        values = loaded.root.obs[key].copy()
    else:
        raise ContractError(f"marker grouping column is absent: {key}")
    values.index = adata.obs_names.astype(str)
    return values


def _scope_slug(value: str | None) -> str:
    if value is None:
        return "all-parents"
    return safe_path_component(value)


def run_marker_analysis_from_config(
    config: CellCuratorConfig | str | Path,
    *,
    assay: str,
    representation: str,
    group_key: str,
    parent_scope: str | None = None,
    parent_membership_key: str | None = None,
    parent_cl_id: str | None = None,
    include_unscoped_programs: bool = False,
    additional_sources: Iterable[MarkerSource] = (),
    output_dir: str | Path | None = None,
    detection_threshold: float = 0.1,
    min_measured_positive: int = 1,
    min_positive_fraction: float = 0.5,
    adjusted_pvalue_max: float = 0.05,
) -> MarkerRun:
    """Execute the marker surface against one declared config assay and parent scope."""

    model = load_config(Path(config)) if isinstance(config, (str, Path)) else config
    if assay not in model.input.assays:
        raise ContractError(f"marker assay is not configured: {assay}")
    assay_config = model.input.assays[assay]
    if representation not in assay_config.representations:
        raise ContractError(
            f"marker representation is not configured: {assay}/{representation}"
        )
    providers: list[MarkerSource] = list(
        providers_from_config(
            model.knowledge.providers,
            allow_remote_enhancement=model.knowledge.allow_optional_remote_enhancement,
        )
    )
    providers.extend(additional_sources)
    programs = assemble_marker_programs(
        providers,
        parent_scope=parent_scope,
        parent_cl_id=parent_cl_id,
        include_unscoped=include_unscoped_programs,
        feature_space=assay_config.representations[representation].feature_space
        or assay_config.feature_space,
    )
    if not programs:
        raise ContractError(
            "no marker programs are callable in the requested parent/feature scope"
        )
    loaded = load_input(model, backed=False)
    try:
        adata = loaded.assays[assay]
        declaration = assay_config.representations[representation]
        matrix = matrix_for(adata, declaration)
        features = feature_ids_for(
            adata,
            assay_config.feature_id_column,
            declaration=declaration,
            gene_mapping_path=assay_config.gene_mapping_path,
        )
        detection = (
            matrix_for(adata, declaration.detection.model_dump())
            if declaration.detection is not None
            else None
        )
        if declaration.detection is not None:
            detection_features = feature_ids_for(
                adata,
                assay_config.feature_id_column,
                declaration=declaration.detection,
                gene_mapping_path=assay_config.gene_mapping_path,
            )
            if not detection_features.equals(features):
                raise ContractError(
                    "marker representation and detection source use different feature axes"
                )
        groups = _obs_column(loaded, adata, group_key)
        mask = np.ones(adata.n_obs, dtype=bool)
        if parent_scope is not None:
            membership_key = parent_membership_key or model.input.keys.canonical_parent
            membership = _obs_column(loaded, adata, membership_key)
            mask = membership.astype(str).eq(parent_scope).to_numpy()
            if not mask.any():
                raise ContractError(f"requested parent scope has no cells: {parent_scope}")
        scoped_groups = groups.iloc[np.flatnonzero(mask)].copy()
        scoped_groups.index = adata.obs_names[mask].astype(str)
        result = run_marker_analysis(
            matrix[mask],
            scoped_groups,
            features,
            programs,
            detection=detection[mask] if detection is not None else None,
            detection_threshold=detection_threshold,
            min_measured_positive=min_measured_positive,
            min_positive_fraction=min_positive_fraction,
            adjusted_pvalue_max=adjusted_pvalue_max,
        )
    finally:
        close_input(loaded)
    target = (
        Path(output_dir)
        if output_dir is not None
        else Path(model.run.output_root)
        / model.run.run_id
        / "markers"
        / f"{assay}__{representation}"
        / _scope_slug(parent_scope)
    )
    return MarkerRun(
        analysis=result,
        artifacts=write_marker_artifacts(result, target, run_id=model.run.run_id),
    )

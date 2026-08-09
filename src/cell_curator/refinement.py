"""Evidence-neutral adaptive nested refinement from frozen representations."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .config import CellCuratorConfig
from .contracts import ContractError, values_sha256
from .data import LoadedInput, close_input, dense, load_input
from .provenance import atomic_publish, publish_json


def adjusted_rand_index(left: np.ndarray, right: np.ndarray) -> float:
    """Adjusted Rand index without a heavyweight clustering dependency."""

    left = np.asarray(left)
    right = np.asarray(right)
    if left.shape != right.shape or left.ndim != 1:
        raise ContractError("ARI inputs must be same-length vectors")
    if len(left) < 2:
        return 1.0
    _, left_codes = np.unique(left, return_inverse=True)
    _, right_codes = np.unique(right, return_inverse=True)
    contingency = np.zeros((left_codes.max() + 1, right_codes.max() + 1), dtype=np.int64)
    np.add.at(contingency, (left_codes, right_codes), 1)

    def choose2(values: np.ndarray) -> np.int64:
        return np.sum(values * (values - 1) // 2)

    sum_comb = float(choose2(contingency))
    row_comb = float(choose2(contingency.sum(axis=1)))
    col_comb = float(choose2(contingency.sum(axis=0)))
    total_comb = len(left) * (len(left) - 1) / 2
    expected = row_comb * col_comb / total_comb if total_comb else 0.0
    maximum = 0.5 * (row_comb + col_comb)
    denominator = maximum - expected
    if denominator == 0:
        return 1.0
    return float((sum_comb - expected) / denominator)


def partition_jaccard_index(left: np.ndarray, right: np.ndarray) -> float:
    """Symmetric macro best-match Jaccard for two clustering partitions."""

    left = np.asarray(left)
    right = np.asarray(right)
    if left.shape != right.shape or left.ndim != 1:
        raise ContractError("partition Jaccard inputs must be same-length vectors")
    if len(left) == 0:
        return 1.0

    def memberships(values: np.ndarray) -> list[set[int]]:
        return [set(np.flatnonzero(values == label).tolist()) for label in np.unique(values)]

    left_groups = memberships(left)
    right_groups = memberships(right)

    def directional(source: list[set[int]], target: list[set[int]]) -> float:
        return float(
            np.mean(
                [
                    max(len(group & candidate) / len(group | candidate) for candidate in target)
                    for group in source
                ]
            )
        )

    return 0.5 * (
        directional(left_groups, right_groups) + directional(right_groups, left_groups)
    )


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def _knn_graph(embedding: np.ndarray, *, k: int) -> nx.Graph:
    if embedding.ndim != 2 or len(embedding) < 2:
        raise ContractError("parent-local embedding requires at least two cells")
    effective_k = min(k + 1, len(embedding))
    distances, indices = cKDTree(embedding).query(embedding, k=effective_k)
    if indices.ndim == 1:
        indices = indices[:, None]
        distances = distances[:, None]
    graph = nx.Graph()
    graph.add_nodes_from(range(len(embedding)))
    positive = distances[:, 1:][distances[:, 1:] > 0]
    scale = float(np.median(positive)) if positive.size else 1.0
    for source in range(len(embedding)):
        for distance, target in zip(distances[source, 1:], indices[source, 1:], strict=True):
            if source == int(target):
                continue
            weight = math.exp(-float(distance) / max(scale, 1e-12))
            if graph.has_edge(source, int(target)):
                graph[source][int(target)]["weight"] = max(
                    graph[source][int(target)]["weight"], weight
                )
            else:
                graph.add_edge(source, int(target), weight=weight)
    return graph


def graph_partition(
    embedding: np.ndarray,
    *,
    resolution: float,
    seed: int,
    k: int,
) -> np.ndarray:
    graph = _knn_graph(embedding, k=k)
    communities = nx.community.louvain_communities(
        graph,
        weight="weight",
        resolution=float(resolution),
        seed=int(seed),
    )
    ordered = sorted(
        (sorted(group) for group in communities), key=lambda group: (group[0], len(group))
    )
    labels = np.empty(len(embedding), dtype=int)
    for label, members in enumerate(ordered):
        labels[members] = label
    return labels


def resampling_stability(
    embedding: np.ndarray,
    reference: np.ndarray,
    *,
    resolution: float,
    seeds: list[int],
    fraction: float,
    k: int,
    metric: str = "ari",
) -> tuple[float, list[float]]:
    if metric not in {"ari", "jaccard"}:
        raise ContractError(f"unsupported resampling stability metric: {metric}")
    scores = []
    sample_size = max(2, min(len(embedding), int(round(len(embedding) * fraction))))
    for seed in seeds:
        rng = np.random.default_rng(seed)
        selected = np.sort(rng.choice(len(embedding), size=sample_size, replace=False))
        replicate = graph_partition(
            embedding[selected],
            resolution=resolution,
            seed=seed,
            k=min(k, sample_size - 1),
        )
        scores.append(
            adjusted_rand_index(reference[selected], replicate)
            if metric == "ari"
            else partition_jaccard_index(reference[selected], replicate)
        )
    return float(np.mean(scores)) if scores else 0.0, scores


def cramer_v(labels: np.ndarray, values: pd.Series) -> float:
    table = pd.crosstab(pd.Series(labels, dtype=str), values.astype(str))
    observed = table.to_numpy(dtype=float)
    if min(observed.shape) <= 1 or observed.sum() == 0:
        return 0.0
    expected = np.outer(observed.sum(axis=1), observed.sum(axis=0)) / observed.sum()
    valid = expected > 0
    chi2 = float(np.sum(((observed - expected) ** 2 / np.where(valid, expected, 1))[valid]))
    phi2 = chi2 / observed.sum()
    rows, columns = observed.shape
    corrected = max(0.0, phi2 - ((columns - 1) * (rows - 1)) / max(observed.sum() - 1, 1))
    corrected_rows = rows - ((rows - 1) ** 2) / max(observed.sum() - 1, 1)
    corrected_columns = columns - ((columns - 1) ** 2) / max(observed.sum() - 1, 1)
    denominator = max(min(corrected_columns - 1, corrected_rows - 1), 1e-12)
    return float(math.sqrt(corrected / denominator))


def eta_squared(labels: np.ndarray, values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(numeric)
    if valid.sum() < 3:
        return 0.0
    numeric = numeric[valid]
    labels = labels[valid]
    total = float(np.sum((numeric - numeric.mean()) ** 2))
    if total <= 0:
        return 0.0
    between = 0.0
    for label in np.unique(labels):
        group = numeric[labels == label]
        between += len(group) * float((group.mean() - numeric.mean()) ** 2)
    return between / total


def discrete_separation(embedding: np.ndarray, labels: np.ndarray) -> float:
    if len(np.unique(labels)) < 2:
        return 0.0
    center = embedding.mean(axis=0)
    total = float(np.sum((embedding - center) ** 2))
    if total <= 0:
        return 0.0
    between = 0.0
    for label in np.unique(labels):
        group = embedding[labels == label]
        between += len(group) * float(np.sum((group.mean(axis=0) - center) ** 2))
    return between / total


def mutual_best_matches(
    left: dict[str, set[str]],
    right: dict[str, set[str]],
    *,
    threshold: float,
) -> list[tuple[str, str, float]]:
    scores = {
        (a, b): jaccard(a_members, b_members)
        for a, a_members in left.items()
        for b, b_members in right.items()
    }
    best_right = {a: max(right, key=lambda b: (scores[(a, b)], b)) for a in left if right}
    best_left = {b: max(left, key=lambda a: (scores[(a, b)], a)) for b in right if left}
    return sorted(
        [
            (a, b, scores[(a, b)])
            for a, b in best_right.items()
            if best_left.get(b) == a and scores[(a, b)] >= threshold
        ],
        key=lambda item: (item[0], item[1]),
    )


@dataclass(frozen=True)
class PartitionEvaluation:
    source: str
    resolution: str
    labels: np.ndarray
    n_clusters: int
    minimum_cluster_size: int
    stability: float
    stability_metric: str
    separation: float
    technical_effect_size: float
    stable: bool
    eligible: bool
    technical_dominated: bool
    gradient_like: bool


def _technical_effect(
    labels: np.ndarray,
    metadata: pd.DataFrame,
    config: CellCuratorConfig,
) -> float:
    keys = config.input.keys
    categorical = [
        key for key in (keys.donor, keys.capture, keys.batch) if key and key in metadata
    ]
    numeric = [key for key in config.input.technical_metrics.values() if key in metadata]
    effects = [cramer_v(labels, metadata[key]) for key in categorical]
    effects.extend(eta_squared(labels, metadata[key]) for key in numeric)
    return max(effects, default=0.0)


def evaluate_local_partitions(
    embedding: np.ndarray,
    metadata: pd.DataFrame,
    config: CellCuratorConfig,
) -> list[PartitionEvaluation]:
    local = config.adaptive.parent_local
    if not local.enabled or local.max_passes_per_parent == 0:
        return []
    evaluations = []
    for resolution in sorted(local.resolutions):
        labels = graph_partition(
            embedding,
            resolution=resolution,
            seed=local.reference_seeds[0],
            k=min(local.knn_k, len(embedding) - 1),
        )
        counts = np.bincount(labels)
        stability, _ = resampling_stability(
            embedding,
            labels,
            resolution=resolution,
            seeds=local.reference_seeds[1:],
            fraction=local.bootstrap_fraction,
            k=min(local.knn_k, len(embedding) - 1),
            metric=local.stability_metric,
        )
        separation = discrete_separation(embedding, labels)
        technical = _technical_effect(labels, metadata, config)
        n_clusters = len(counts)
        stable = stability > local.stability_threshold_strictly_greater_than
        eligible = (
            n_clusters >= 2 and int(counts.min()) >= config.adaptive.eligibility.min_cells
        )
        technical_limit = config.adaptive.audit.technical_effect_size_max
        technical_dominated = technical_limit is not None and technical >= technical_limit
        gradient_like = separation < 0.1
        evaluations.append(
            PartitionEvaluation(
                source="parent_local_graph",
                resolution=f"{resolution:g}",
                labels=labels,
                n_clusters=n_clusters,
                minimum_cluster_size=int(counts.min()),
                stability=stability,
                stability_metric=local.stability_metric,
                separation=separation,
                technical_effect_size=technical,
                stable=stable,
                eligible=eligible,
                technical_dominated=technical_dominated,
                gradient_like=gradient_like,
            )
        )
    return evaluations


def _embedding(loaded: LoadedInput, config: CellCuratorConfig) -> np.ndarray:
    adata = next(iter(loaded.assays.values()))
    key = config.input.keys.embedding
    if key:
        if key not in adata.obsm:
            raise ContractError(f"declared refinement embedding is absent: {key}")
        return np.asarray(adata.obsm[key])
    matrix = dense(adata.X)
    centered = matrix - matrix.mean(axis=0)
    u, singular, _ = np.linalg.svd(centered, full_matrices=False)
    dimensions = min(20, len(singular))
    return u[:, :dimensions] * singular[:dimensions]


def _stored_evaluation(
    parent_obs: pd.DataFrame,
    embedding: np.ndarray,
    config: CellCuratorConfig,
) -> PartitionEvaluation | None:
    keys = config.input.keys.stored_resolutions
    if len(keys) < 2:
        return None
    threshold = config.adaptive.stored_candidates.adjacent_jaccard_min
    minimum_matches = config.adaptive.stored_candidates.min_adjacent_matches
    previous: dict[str, set[str]] | None = None
    persistence: dict[tuple[str, str], int] = {}
    for key in keys:
        current = {
            str(label): set(
                parent_obs.index[parent_obs[key].astype(str).eq(str(label))].astype(str)
            )
            for label in sorted(parent_obs[key].astype(str).unique())
        }
        if previous is not None:
            for left, right, _ in mutual_best_matches(previous, current, threshold=threshold):
                persistence[(key, right)] = max(
                    persistence.get((key, right), 0),
                    max(
                        (
                            value
                            for (stored_key, child), value in persistence.items()
                            if child == left
                        ),
                        default=0,
                    )
                    + 1,
                )
        previous = current
    for key in keys:
        labels = parent_obs[key].astype(str).to_numpy()
        counts = pd.Series(labels).value_counts()
        persistent = sum(
            1
            for label in counts.index
            if persistence.get((key, str(label)), 0) >= minimum_matches
        )
        if persistent < 2:
            continue
        encoded = pd.Categorical(labels).codes
        separation = discrete_separation(embedding, encoded)
        technical = _technical_effect(encoded, parent_obs, config)
        limit = config.adaptive.audit.technical_effect_size_max
        return PartitionEvaluation(
            source="stored_resolution",
            resolution=key,
            labels=encoded,
            n_clusters=int(counts.size),
            minimum_cluster_size=int(counts.min()),
            stability=1.0,
            stability_metric="adjacent_jaccard",
            separation=separation,
            technical_effect_size=technical,
            stable=True,
            eligible=int(counts.min()) >= config.adaptive.eligibility.min_cells,
            technical_dominated=limit is not None and technical >= limit,
            gradient_like=separation < 0.1,
        )
    return None


def discover_refinement_candidates(
    *,
    config: CellCuratorConfig,
    run_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    loaded = load_input(config, backed=False)
    try:
        obs = loaded.root.obs.copy()
        embedding = _embedding(loaded, config)
        keys = config.input.keys
        registry_rows: list[dict[str, Any]] = []
        membership_rows: list[dict[str, Any]] = []
        evaluation_rows: list[dict[str, Any]] = []
        for (scope, parent), indices in obs.groupby(
            [keys.scope, keys.canonical_parent], observed=True
        ).indices.items():
            positions = np.asarray(indices, dtype=int)
            parent_obs = obs.iloc[positions]
            parent_embedding = embedding[positions]
            if len(parent_obs) < 2 * config.adaptive.eligibility.min_cells:
                evaluation_rows.append(
                    {
                        "scope": str(scope),
                        "parent_id": str(parent),
                        "source": "none",
                        "resolution": "",
                        "n_clusters": 1,
                        "stability": 0.0,
                        "separation": 0.0,
                        "technical_effect_size": 0.0,
                        "action": "RETAIN_TOO_SMALL",
                    }
                )
                continue
            stored = _stored_evaluation(parent_obs, parent_embedding, config)
            evaluations = (
                [stored]
                if stored is not None
                else (
                    evaluate_local_partitions(parent_embedding, parent_obs, config)
                    if config.adaptive.parent_local.enabled
                    and config.adaptive.parent_local.max_passes_per_parent > 0
                    else []
                )
            )
            selected = next(
                (
                    item
                    for item in evaluations
                    if item is not None
                    and item.stable
                    and item.eligible
                    and not item.technical_dominated
                    and not item.gradient_like
                ),
                None,
            )
            for item in evaluations:
                if item is None:
                    continue
                evaluation_rows.append(
                    {
                        "scope": str(scope),
                        "parent_id": str(parent),
                        "source": item.source,
                        "resolution": item.resolution,
                        "n_clusters": item.n_clusters,
                        "minimum_cluster_size": item.minimum_cluster_size,
                        "stability": item.stability,
                        "stability_metric": item.stability_metric,
                        "separation": item.separation,
                        "technical_effect_size": item.technical_effect_size,
                        "stable": item.stable,
                        "eligible": item.eligible,
                        "technical_dominated": item.technical_dominated,
                        "gradient_like": item.gradient_like,
                        "action": "SELECTED_FOR_EVIDENCE" if item is selected else "REJECTED",
                    }
                )
            if selected is None:
                continue
            cluster_ids = sorted(np.unique(selected.labels))
            centroids = {
                child: parent_embedding[selected.labels == child].mean(axis=0)
                for child in cluster_ids
            }
            for child in cluster_ids:
                mask = selected.labels == child
                candidate_key = f"{scope}__p{parent}__c{child}"
                closest = min(
                    (other for other in cluster_ids if other != child),
                    key=lambda other: (
                        float(np.linalg.norm(centroids[child] - centroids[other])),
                        int(other),
                    ),
                )
                members = parent_obs.index[mask].astype(str).tolist()
                registry_rows.append(
                    {
                        "schema_version": 3,
                        "scope": str(scope),
                        "parent_id": str(parent),
                        "candidate_key": candidate_key,
                        "candidate_child_uid": f"{scope}::{parent}::candidate_{child}",
                        "source_kind": selected.source,
                        "source_resolution": selected.resolution,
                        "n_cells_in_parent": int(mask.sum()),
                        "n_parent_remainder": int((~mask).sum()),
                        "stability_score": selected.stability,
                        "separation_score": selected.separation,
                        "technical_effect_size": selected.technical_effect_size,
                        "closest_sibling": f"{scope}__p{parent}__c{closest}",
                        "membership_sha256": values_sha256(sorted(members)),
                        "triage_action": "FULL_EVIDENCE_TESTING",
                        "cell_type": "",
                        "cl_id": "",
                    }
                )
                for barcode, selected_value in zip(
                    parent_obs.index.astype(str), mask, strict=True
                ):
                    membership_rows.append(
                        {
                            "schema_version": 3,
                            "scope": str(scope),
                            "parent_id": str(parent),
                            "candidate_key": candidate_key,
                            "barcode": barcode,
                            "candidate_membership": bool(selected_value),
                        }
                    )
        registry = pd.DataFrame(
            registry_rows,
            columns=[
                "schema_version",
                "scope",
                "parent_id",
                "candidate_key",
                "candidate_child_uid",
                "source_kind",
                "source_resolution",
                "n_cells_in_parent",
                "n_parent_remainder",
                "stability_score",
                "separation_score",
                "technical_effect_size",
                "closest_sibling",
                "membership_sha256",
                "triage_action",
                "cell_type",
                "cl_id",
            ],
        )
        membership = pd.DataFrame(
            membership_rows,
            columns=[
                "schema_version",
                "scope",
                "parent_id",
                "candidate_key",
                "barcode",
                "candidate_membership",
            ],
        )
        evaluations = pd.DataFrame(
            evaluation_rows,
            columns=[
                "scope",
                "parent_id",
                "source",
                "resolution",
                "n_clusters",
                "minimum_cluster_size",
                "stability",
                "stability_metric",
                "separation",
                "technical_effect_size",
                "stable",
                "eligible",
                "technical_dominated",
                "gradient_like",
                "action",
            ],
        )
    finally:
        close_input(loaded)
    discovery_dir = run_root / "discovery"
    discovery_dir.mkdir(parents=True, exist_ok=True)
    for path, frame in (
        (discovery_dir / "candidate_registry.tsv", registry),
        (discovery_dir / "candidate_membership.tsv", membership),
        (discovery_dir / "partition_evaluations.tsv", evaluations),
    ):
        atomic_publish(path, frame.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    publish_json(
        discovery_dir / "discovery.manifest.json",
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_adaptive_refinement",
            "run_id": config.run.run_id,
            "status": "complete",
            "n_candidates": len(registry),
            "n_parent_partitions_evaluated": len(evaluations),
            "candidate_labels_blank": True,
            "decisions_derived_from_computed_stability": True,
        },
    )
    return registry, membership, evaluations

"""Optional automated annotation hypotheses and reference mapping providers.

Every function in this module produces evidence, never an authoritative label.  The
CPU reference mapper is dependency-light and fully local.  Optional providers return
a structured ``unavailable`` record when their declared dependency or resource is
missing; GPU routes validate real hardware and never substitute the CPU mapper.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import multiprocessing
import os
import queue
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Protocol, cast

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

from .contracts import ContractError, sha256
from .data import gene_symbols, load_gene_mapping
from .provenance import atomic_publish, publish_json


class _ProcessFactoryContext(Protocol):
    """Subset of concrete multiprocessing contexts omitted from ``BaseContext`` stubs."""

    Process: type[multiprocessing.Process]


def _census_value_filter(tissue: str) -> str:
    """Render the documented SOMA value-filter string literal form safely."""

    normalized = tissue.strip()
    if not normalized:
        raise ContractError("Census tissue must be nonblank")
    if any(ord(character) < 32 for character in normalized):
        raise ContractError("Census tissue must not contain control characters")
    escaped = normalized.replace("\\", "\\\\").replace("'", "\\'")
    return f"tissue_general == '{escaped}'"


def _normalized_census_labels(label_columns: Iterable[str]) -> tuple[str, ...]:
    labels = tuple(str(column).strip() for column in label_columns)
    if not labels or any(not column for column in labels):
        raise ContractError(
            "Census label_columns must contain at least one nonblank column name"
        )
    if len(labels) != len(set(labels)):
        raise ContractError("Census label_columns must not contain duplicates")
    if "soma_joinid" in labels:
        raise ContractError(
            "Census label_columns must not include the reserved soma_joinid column"
        )
    return labels


def _nonblank_label_mask(frame: pd.DataFrame, label_columns: tuple[str, ...]) -> pd.Series:
    missing = sorted(set(label_columns).difference(frame.columns))
    if missing:
        raise ContractError(f"Census response is missing requested label columns: {missing}")
    keep = pd.Series(True, index=frame.index, dtype=bool)
    for column in label_columns:
        keep &= frame[column].astype("string").str.strip().fillna("").ne("")
    return keep


def _status(
    provider: str,
    status: str,
    *,
    reason: str = "",
    remediation: str = "",
    artifacts: dict[str, str] | None = None,
) -> dict[str, Any]:
    return {
        "schema_version": 3,
        "artifact_kind": "cell_curator_optional_provider_status",
        "provider": provider,
        "status": status,
        "reason": reason,
        "remediation": remediation,
        "artifacts": artifacts or {},
        "automated_output_role": "hypothesis_not_verdict",
    }


def _publish_status(output_dir: Path | None, value: dict[str, Any]) -> dict[str, Any]:
    if output_dir is not None:
        publish_json(output_dir / "provider_status.json", value)
    return value


def _missing_dependency(
    provider: str,
    package: str,
    extra: str,
    output_dir: Path | None,
) -> dict[str, Any]:
    return _publish_status(
        output_dir,
        _status(
            provider,
            "unavailable",
            reason=f"optional Python distribution {package!r} is not installed",
            remediation=f"Install cell-curator[{extra}] in the execution environment and retry.",
        ),
    )


def _read_object(path: str | Path, *, modality: str = "") -> Any:
    source = Path(path)
    if not source.is_file():
        raise ContractError(f"single-cell object does not exist: {source}")
    if source.suffix.lower() == ".h5mu":
        import mudata as md

        root = md.read_h5mu(source)
        if not modality or modality not in root.mod:
            raise ContractError(
                f"MuData object requires an existing modality; got {modality!r}"
            )
        return root.mod[modality]
    if source.suffix.lower() != ".h5ad":
        raise ContractError("reference providers accept .h5ad or .h5mu inputs")
    import anndata as ad

    return ad.read_h5ad(source)


def read_single_cell_object(path: str | Path, *, modality: str = "") -> Any:
    """Read one AnnData assay or an explicitly selected MuData modality."""

    return _read_object(path, modality=modality)


def _symbol_index(
    adata: Any,
    *,
    symbol_column: str = "",
    symbol_map: dict[str, str] | None = None,
) -> pd.Index:
    """Resolve gene symbols without depending on a remote identifier service."""

    return gene_symbols(adata, symbol_column=symbol_column, mapping=symbol_map)


def load_symbol_map(path: str | Path | None) -> dict[str, str]:
    if path is None:
        return {}
    return load_gene_mapping(path)


def _matrix(adata: Any, layer: str = "") -> Any:
    if layer:
        if layer == "raw":
            if adata.raw is None:
                raise ContractError("requested .raw matrix is absent")
            return adata.raw.X
        if layer not in adata.layers:
            raise ContractError(f"requested matrix layer is absent: {layer}")
        return adata.layers[layer]
    return adata.X


def _as_dense_selected(matrix: Any, columns: np.ndarray) -> np.ndarray:
    selected = matrix[:, columns]
    if sparse.issparse(selected):
        selected = selected.toarray()
    values = np.asarray(selected, dtype=np.float32)
    if not np.isfinite(values).all():
        raise ContractError("reference mapping matrix contains non-finite values")
    return values


def harmonize_reference_genes(
    query: Any,
    reference: Any,
    *,
    query_symbol_column: str = "",
    reference_symbol_column: str = "",
    query_symbol_map: dict[str, str] | None = None,
    reference_symbol_map: dict[str, str] | None = None,
    max_genes: int = 2_000,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return aligned query/reference column positions with deterministic HVG capping."""

    query_symbols = _symbol_index(
        query, symbol_column=query_symbol_column, symbol_map=query_symbol_map
    )
    reference_symbols = _symbol_index(
        reference, symbol_column=reference_symbol_column, symbol_map=reference_symbol_map
    )
    qpos = {value: index for index, value in enumerate(query_symbols)}
    rpos = {value: index for index, value in enumerate(reference_symbols)}
    common = sorted(set(qpos).intersection(rpos))
    if len(common) < 2:
        raise ContractError(
            "reference mapping found fewer than two shared harmonized genes; provide symbol columns or --*-symbol-map"
        )
    q_columns = np.asarray([qpos[item] for item in common], dtype=int)
    r_columns = np.asarray([rpos[item] for item in common], dtype=int)
    if len(common) > max_genes:
        # Rank by reference variance using a bounded deterministic sample.  This avoids
        # materializing a full atlas solely to choose columns.
        sample_rows = np.linspace(0, reference.n_obs - 1, min(reference.n_obs, 5_000)).astype(
            int
        )
        sampled = reference.X[sample_rows][:, r_columns]
        sampled = sampled.toarray() if sparse.issparse(sampled) else np.asarray(sampled)
        variance = np.var(sampled, axis=0)
        selected = np.lexsort((np.asarray(common), -variance))[:max_genes]
        q_columns = q_columns[selected]
        r_columns = r_columns[selected]
        common = [common[index] for index in selected]
    return q_columns, r_columns, common


def existing_prediction_hypothesis(
    adata: Any,
    *,
    prediction_column: str,
    cluster_column: str,
    output_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Summarize an existing prediction column without promoting it to a verdict."""

    for column in (prediction_column, cluster_column):
        if column not in adata.obs:
            raise ContractError(f"existing-prediction column is absent: {column}")
    cells = pd.DataFrame(
        {
            "barcode": adata.obs_names.astype(str),
            "cluster_id": adata.obs[cluster_column].astype(str).to_numpy(),
            "predicted_label": adata.obs[prediction_column].astype(str).to_numpy(),
        }
    )
    rows: list[dict[str, Any]] = []
    for cluster_id, block in cells.groupby("cluster_id", observed=True, sort=True):
        fractions = block["predicted_label"].value_counts(normalize=True)
        winner = str(fractions.index[0])
        probability = float(fractions.iloc[0])
        entropy = -sum(
            float(value) * math.log(float(value)) for value in fractions if value > 0
        )
        normalized = entropy / math.log(len(fractions)) if len(fractions) > 1 else 0.0
        rows.append(
            {
                "cluster_id": str(cluster_id),
                "predicted_label": winner,
                "support": probability,
                "uncertainty": float(normalized),
                "n_cells": len(block),
                "source": f"obs:{prediction_column}",
                "role": "hypothesis_not_verdict",
            }
        )
    clusters = pd.DataFrame(rows)
    root = Path(output_dir) if output_dir is not None else None
    artifacts: dict[str, str] = {}
    if root is not None:
        root.mkdir(parents=True, exist_ok=True)
        cell_path = root / "cell_predictions.tsv"
        cluster_path = root / "cluster_predictions.tsv"
        atomic_publish(cell_path, cells.to_csv(sep="\t", index=False).encode())
        atomic_publish(cluster_path, clusters.to_csv(sep="\t", index=False).encode())
        artifacts = {"cells": str(cell_path), "clusters": str(cluster_path)}
    status = _publish_status(
        root, _status("existing_prediction", "complete", artifacts=artifacts)
    )
    return cells, clusters, status


def _latent_pair(
    query_values: np.ndarray,
    reference_values: np.ndarray,
    *,
    components: int,
    max_memory_mb: int,
) -> tuple[np.ndarray, np.ndarray]:
    estimated = (query_values.size + reference_values.size) * 4
    if estimated > max_memory_mb * 1024 * 1024:
        raise ContractError(
            f"CPU reference mapping needs about {estimated / 1024**2:.0f} MiB for the aligned matrix, "
            f"above the {max_memory_mb} MiB limit; reduce the reference or route to an authorized larger-memory profile"
        )
    combined = np.vstack([reference_values, query_values]).astype(np.float32, copy=False)
    center = combined.mean(axis=0, keepdims=True)
    spread = combined.std(axis=0, keepdims=True)
    spread[spread <= 1e-6] = 1.0
    combined = (combined - center) / spread
    n_components = max(1, min(components, combined.shape[0] - 1, combined.shape[1]))
    u, singular, _ = np.linalg.svd(combined, full_matrices=False)
    latent = u[:, :n_components] * singular[:n_components]
    return latent[len(reference_values) :], latent[: len(reference_values)]


def _vote_neighbors(
    indices: np.ndarray,
    distances: np.ndarray,
    labels: np.ndarray,
) -> tuple[list[str], list[float]]:
    predictions: list[str] = []
    uncertainty: list[float] = []
    for row_indices, row_distances in zip(indices, distances, strict=True):
        weights = 1.0 / np.maximum(np.asarray(row_distances, dtype=float), 1e-8)
        totals: dict[str, float] = {}
        for index, weight in zip(row_indices, weights, strict=True):
            label = str(labels[int(index)])
            totals[label] = totals.get(label, 0.0) + float(weight)
        ordered = sorted(totals.items(), key=lambda item: (-item[1], item[0]))
        total = sum(totals.values())
        predictions.append(ordered[0][0])
        uncertainty.append(1.0 - ordered[0][1] / total if total else 1.0)
    return predictions, uncertainty


def reference_knn_hypothesis(
    query: Any,
    reference: Any,
    *,
    reference_label_column: str,
    query_cluster_column: str = "",
    query_layer: str = "",
    reference_layer: str = "",
    query_symbol_column: str = "",
    reference_symbol_column: str = "",
    query_symbol_map: dict[str, str] | None = None,
    reference_symbol_map: dict[str, str] | None = None,
    backend: str = "cpu",
    k: int = 15,
    components: int = 50,
    max_memory_mb: int = 8_192,
    output_dir: str | Path | None = None,
    provider_name: str = "reference_knn",
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Map a labeled local reference by kNN and retain uncertainty/disagreement."""

    gpu_neighbors: Any | None = None
    if backend == "gpu":
        from .environment import require_real_gpu

        require_real_gpu()
        if importlib.util.find_spec("cuml") is None:
            raise ContractError(
                "GPU kNN was requested and a real GPU was detected, but cuML is unavailable; "
                "install a CUDA-compatible cuML build in the declared GPU environment. "
                "CPU substitution is prohibited"
            )
        from cuml.neighbors import NearestNeighbors

        gpu_neighbors = NearestNeighbors
    elif backend != "cpu":
        raise ContractError("reference kNN backend must be 'cpu' or 'gpu'")

    if len(reference) < 1:
        raise ContractError("reference kNN requires at least one labeled reference cell")
    if reference_label_column not in reference.obs:
        raise ContractError(f"reference label column is absent: {reference_label_column}")
    reference_labels = reference.obs[reference_label_column].astype("string").str.strip()
    invalid_labels = reference_labels.fillna("").eq("")
    if invalid_labels.any():
        raise ContractError(
            f"reference label column {reference_label_column!r} contains "
            f"{int(invalid_labels.sum())} missing or blank labels"
        )
    if query_cluster_column and query_cluster_column not in query.obs:
        raise ContractError(f"query cluster column is absent: {query_cluster_column}")
    if k < 1:
        raise ContractError("reference kNN requires k >= 1")
    q_columns, r_columns, common = harmonize_reference_genes(
        query,
        reference,
        query_symbol_column=query_symbol_column,
        reference_symbol_column=reference_symbol_column,
        query_symbol_map=query_symbol_map,
        reference_symbol_map=reference_symbol_map,
    )
    q_values = _as_dense_selected(_matrix(query, query_layer), q_columns)
    r_values = _as_dense_selected(_matrix(reference, reference_layer), r_columns)
    q_latent, r_latent = _latent_pair(
        q_values, r_values, components=components, max_memory_mb=max_memory_mb
    )
    actual_k = min(k, len(reference))
    if backend == "cpu":
        distances, indices = cKDTree(r_latent).query(q_latent, k=actual_k)
        if actual_k == 1:
            distances = np.asarray(distances)[:, None]
            indices = np.asarray(indices)[:, None]
    else:
        assert gpu_neighbors is not None
        model = gpu_neighbors(n_neighbors=actual_k)
        model.fit(r_latent)
        distances, indices = model.kneighbors(q_latent)
        distances, indices = np.asarray(distances), np.asarray(indices)
    labels = reference_labels.astype(str).to_numpy()
    predictions, uncertainty = _vote_neighbors(indices, distances, labels)
    cells = pd.DataFrame(
        {
            "barcode": query.obs_names.astype(str),
            "predicted_label": predictions,
            "uncertainty": uncertainty,
            "backend": backend,
            "provider": provider_name,
            "role": "hypothesis_not_verdict",
        }
    )
    if query_cluster_column:
        cells["cluster_id"] = query.obs[query_cluster_column].astype(str).to_numpy()
        cluster_rows: list[dict[str, Any]] = []
        for cluster_id, block in cells.groupby("cluster_id", observed=True, sort=True):
            counts = block["predicted_label"].value_counts()
            winner = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]
            cluster_rows.append(
                {
                    "cluster_id": str(cluster_id),
                    "predicted_label": str(winner),
                    "support": float(counts.max() / counts.sum()),
                    "uncertainty": float(block["uncertainty"].mean()),
                    "n_cells": len(block),
                }
            )
        clusters = pd.DataFrame(cluster_rows)
    else:
        clusters = pd.DataFrame()
    root = Path(output_dir) if output_dir is not None else None
    artifacts: dict[str, str] = {}
    if root is not None:
        root.mkdir(parents=True, exist_ok=True)
        cell_path = root / "cell_predictions.tsv"
        cluster_path = root / "cluster_predictions.tsv"
        gene_path = root / "shared_genes.txt"
        atomic_publish(cell_path, cells.to_csv(sep="\t", index=False).encode())
        atomic_publish(cluster_path, clusters.to_csv(sep="\t", index=False).encode())
        atomic_publish(gene_path, ("\n".join(common) + "\n").encode())
        artifacts = {
            "cells": str(cell_path),
            "clusters": str(cluster_path),
            "shared_genes": str(gene_path),
        }
    status = _publish_status(
        root,
        _status(provider_name, "complete", artifacts=artifacts),
    )
    status["backend"] = backend
    status["n_shared_genes"] = len(common)
    status["uncertainty_semantics"] = "one_minus_inverse_distance_vote_fraction"
    if root is not None:
        publish_json(root / "provider_status.json", status)
    return cells, clusters, status


def map_reference_files(
    *,
    query_path: str | Path,
    reference_path: str | Path,
    reference_label_column: str,
    query_cluster_column: str = "",
    query_modality: str = "",
    reference_modality: str = "",
    output_dir: str | Path | None = None,
    **kwargs: Any,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    query = _read_object(query_path, modality=query_modality)
    reference = _read_object(reference_path, modality=reference_modality)
    return reference_knn_hypothesis(
        query,
        reference,
        reference_label_column=reference_label_column,
        query_cluster_column=query_cluster_column,
        output_dir=output_dir,
        **kwargs,
    )


def existing_prediction_file(
    *,
    input_path: str | Path,
    prediction_column: str,
    cluster_column: str,
    modality: str = "",
    output_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    return existing_prediction_hypothesis(
        _read_object(input_path, modality=modality),
        prediction_column=prediction_column,
        cluster_column=cluster_column,
        output_dir=output_dir,
    )


def _celltypist_cache_roots(models_module: Any) -> tuple[Path, ...]:
    roots: list[Path] = []
    for attribute in ("models_path", "model_path", "MODEL_PATH"):
        raw = getattr(models_module, attribute, None)
        if isinstance(raw, (str, os.PathLike)) and str(raw).strip():
            path = Path(raw).expanduser().resolve()
            if path not in roots:
                roots.append(path)
    return tuple(roots)


def _celltypist_index_names(cache_roots: tuple[Path, ...]) -> set[str]:
    """Read only already-cached model indexes; never invoke a catalog/download API."""

    identifiers: set[str] = set()
    checked: set[Path] = set()
    for root in cache_roots:
        for path in (root / "models.json", root.parent / "models.json"):
            if path in checked or not path.is_file():
                continue
            checked.add(path)
            try:
                payload = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            items = payload.get("models", ()) if isinstance(payload, dict) else payload
            if not isinstance(items, list):
                continue
            for item in items:
                if isinstance(item, str):
                    identifier = item
                elif isinstance(item, dict):
                    identifier = next(
                        (
                            str(item[key])
                            for key in ("filename", "model", "name")
                            if item.get(key)
                        ),
                        "",
                    )
                else:
                    identifier = ""
                identifier = Path(identifier.strip()).name if identifier.strip() else ""
                if identifier:
                    identifiers.add(identifier)
    return identifiers


def _celltypist_catalog(models_module: Any) -> tuple[list[dict[str, Any]], tuple[Path, ...]]:
    roots = _celltypist_cache_roots(models_module)
    cached: dict[str, Path] = {}
    for root in roots:
        if root.is_dir():
            for path in root.glob("*.pkl"):
                if path.is_file():
                    cached[path.name] = path.resolve()
    identifiers = set(cached).union(_celltypist_index_names(roots))
    entries: list[dict[str, Any]] = []
    for identifier in sorted(identifiers):
        cached_path = cached.get(identifier)
        entry: dict[str, Any] = {
            "identifier": identifier,
            "cached": cached_path is not None,
            "availability": "cached" if cached_path is not None else "remote_metadata_only",
            "path": str(cached_path) if cached_path is not None else "",
            "sha256": "",
        }
        if cached_path is not None:
            try:
                entry["sha256"] = sha256(cached_path)
            except OSError as exc:
                entry["read_error"] = str(exc)
        entries.append(entry)
    return entries, roots


def _resolve_celltypist_model(
    model: str,
    models_module: Any,
) -> tuple[Path | None, str, tuple[Path, ...]]:
    identifier = model.strip()
    roots = _celltypist_cache_roots(models_module)
    if not identifier:
        return None, "", roots
    direct = Path(identifier).expanduser()
    if direct.is_file():
        return direct.resolve(), "explicit_local_path", roots
    if Path(identifier).name == identifier:
        for root in roots:
            candidate = root / identifier
            if candidate.is_file():
                return candidate.resolve(), "celltypist_cache", roots
    return None, "", roots


def _celltypist_model_status(
    *,
    identifier: str,
    resolved: Path | None,
    source: str,
    cache_roots: tuple[Path, ...],
    celltypist_version: str,
    model_sha256: str = "",
) -> dict[str, Any]:
    return {
        "model_identifier": identifier,
        "resolved_model_name": resolved.name if resolved is not None else "",
        "model_path": str(resolved) if resolved is not None else "",
        "model_source": source,
        "model_sha256": model_sha256,
        "model_version": f"sha256:{model_sha256}" if model_sha256 else "",
        "celltypist_version": celltypist_version,
        "cache_roots_checked": [str(path) for path in cache_roots],
        "network_accessed": False,
        "automatic_downloads_allowed": False,
    }


def discover_celltypist_models(*, output_dir: str | Path | None = None) -> dict[str, Any]:
    root = Path(output_dir) if output_dir is not None else None
    if importlib.util.find_spec("celltypist") is None:
        return _missing_dependency("celltypist", "celltypist", "celltypist", root)
    try:
        import celltypist
        from celltypist import models

        entries, cache_roots = _celltypist_catalog(models)
    except Exception as exc:  # local import/cache inspection failures must be explicit
        return _publish_status(
            root,
            _status(
                "celltypist",
                "failed",
                reason=f"CellTypist local-cache discovery failed: {exc}",
                remediation="Check the installed CellTypist package and local cache permissions; discovery never downloads models.",
            ),
        )
    value = _status("celltypist", "complete")
    value.update(
        {
            "models": [entry["identifier"] for entry in entries],
            "model_entries": entries,
            "cache_roots_checked": [str(path) for path in cache_roots],
            "celltypist_version": str(getattr(celltypist, "__version__", "unknown")),
            "discovery_mode": "local_cache_only",
            "remote_catalog_consulted": False,
            "network_accessed": False,
            "downloaded": False,
        }
    )
    return _publish_status(root, value)


def celltypist_hypothesis(
    adata: Any,
    *,
    model: str,
    cluster_column: str = "",
    output_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    root = Path(output_dir) if output_dir is not None else None
    if importlib.util.find_spec("celltypist") is None:
        return pd.DataFrame(), _missing_dependency(
            "celltypist", "celltypist", "celltypist", root
        )
    if cluster_column and cluster_column not in adata.obs:
        raise ContractError(f"CellTypist cluster column is absent: {cluster_column}")
    try:
        import celltypist
        from celltypist import models

        celltypist_version = str(getattr(celltypist, "__version__", "unknown"))
        resolved, model_source, cache_roots = _resolve_celltypist_model(model, models)
    except Exception as exc:
        return pd.DataFrame(), _publish_status(
            root,
            _status(
                "celltypist",
                "failed",
                reason=f"CellTypist local model resolution failed: {exc}",
                remediation="Check the installed CellTypist package and local cache permissions; automatic downloads are disabled.",
            ),
        )
    if resolved is None:
        value = _status(
            "celltypist",
            "unavailable",
            reason=f"CellTypist model {model!r} is not present as a local file or cached model",
            remediation="Cache the model explicitly during an authorized preparation step, then pass its local file path or cached CellTypist identifier; automatic downloads are disabled.",
        )
        value.update(
            _celltypist_model_status(
                identifier=model.strip(),
                resolved=None,
                source="",
                cache_roots=cache_roots,
                celltypist_version=celltypist_version,
            )
        )
        return pd.DataFrame(), _publish_status(root, value)
    model_digest = ""
    try:
        model_digest = sha256(resolved)

        result = celltypist.annotate(
            adata,
            model=str(resolved),
            majority_voting=bool(cluster_column),
            over_clustering=cluster_column or None,
        )
        frame = result.predicted_labels.copy()
        label_column = next(
            (column for column in ("majority_voting", "predicted_labels") if column in frame),
            str(frame.columns[0]),
        )
        cells = pd.DataFrame(
            {
                "barcode": adata.obs_names.astype(str),
                "predicted_label": frame[label_column].astype(str).to_numpy(),
                "uncertainty": (
                    1.0 - frame["conf_score"].astype(float).to_numpy()
                    if "conf_score" in frame
                    else np.nan
                ),
                "provider": "celltypist",
                "role": "hypothesis_not_verdict",
            }
        )
    except Exception as exc:
        value = _status(
            "celltypist",
            "failed",
            reason=f"CellTypist prediction failed: {exc}",
            remediation="Verify that the local model matches the tissue and the matrix is log1p-CP10k; automatic downloads remain disabled.",
        )
        value.update(
            _celltypist_model_status(
                identifier=model.strip(),
                resolved=resolved,
                source=model_source,
                cache_roots=cache_roots,
                celltypist_version=celltypist_version,
                model_sha256=model_digest,
            )
        )
        return pd.DataFrame(), _publish_status(root, value)
    artifacts: dict[str, str] = {}
    if root is not None:
        root.mkdir(parents=True, exist_ok=True)
        path = root / "cell_predictions.tsv"
        atomic_publish(path, cells.to_csv(sep="\t", index=False).encode())
        artifacts["cells"] = str(path)
    value = _status("celltypist", "complete", artifacts=artifacts)
    value.update(
        _celltypist_model_status(
            identifier=model.strip(),
            resolved=resolved,
            source=model_source,
            cache_roots=cache_roots,
            celltypist_version=celltypist_version,
            model_sha256=model_digest,
        )
    )
    return cells, _publish_status(root, value)


def predict_celltypist_file(
    *,
    input_path: str | Path,
    model: str,
    cluster_column: str = "",
    modality: str = "",
    output_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    return celltypist_hypothesis(
        _read_object(input_path, modality=modality),
        model=model,
        cluster_column=cluster_column,
        output_dir=output_dir,
    )


def _gpu_optional_status(
    provider: str,
    package: str,
    extra: str,
    output_dir: Path | None,
) -> dict[str, Any] | None:
    from .environment import require_real_gpu

    require_real_gpu()
    if importlib.util.find_spec(package) is None:
        return _missing_dependency(provider, package, extra, output_dir)
    return None


def train_scanvi_reference(
    reference: Any,
    *,
    labels_key: str,
    output_dir: str | Path,
    counts_layer: str = "counts",
    max_epochs: int = 200,
    seed: int = 0,
) -> dict[str, Any]:
    """Train a real scVI/scANVI reference under a validated GPU contract."""

    root = Path(output_dir)
    unavailable = _gpu_optional_status("scanvi_training", "scvi", "reference-models", root)
    if unavailable is not None:
        return unavailable
    if labels_key not in reference.obs:
        raise ContractError(f"scANVI reference label column is absent: {labels_key}")
    if counts_layer not in reference.layers:
        raise ContractError(f"scANVI requires a raw-count layer: {counts_layer}")
    try:
        import scvi

        scvi.settings.seed = seed
        scvi.model.SCVI.setup_anndata(reference, layer=counts_layer)
        vae = scvi.model.SCVI(reference)
        vae.train(max_epochs=max_epochs, accelerator="gpu", devices=1)
        labels = reference.obs[labels_key].astype("category")
        reference.obs[labels_key] = labels
        unlabeled = "__cell_curator_unlabeled__"
        if unlabeled not in labels.cat.categories:
            reference.obs[labels_key] = labels.cat.add_categories([unlabeled])
        scanvi = scvi.model.SCANVI.from_scvi_model(
            vae, labels_key=labels_key, unlabeled_category=unlabeled
        )
        scanvi.train(max_epochs=max_epochs, accelerator="gpu", devices=1)
        model_dir = root / "model"
        root.mkdir(parents=True, exist_ok=True)
        scanvi.save(model_dir, overwrite=False, save_anndata=False)
    except Exception as exc:
        return _publish_status(
            root,
            _status(
                "scanvi_training",
                "failed",
                reason=f"scVI/scANVI training failed without fallback: {exc}",
                remediation="Run in the declared compatible GPU environment and inspect the scvi-tools training log.",
            ),
        )
    value = _status("scanvi_training", "complete", artifacts={"model": str(model_dir)})
    value.update({"labels_key": labels_key, "seed": seed, "gpu_required": True})
    return _publish_status(root, value)


def train_scanvi_reference_file(
    *,
    reference_path: str | Path,
    labels_key: str,
    output_dir: str | Path,
    modality: str = "",
    counts_layer: str = "counts",
    max_epochs: int = 200,
    seed: int = 0,
) -> dict[str, Any]:
    return train_scanvi_reference(
        _read_object(reference_path, modality=modality),
        labels_key=labels_key,
        output_dir=output_dir,
        counts_layer=counts_layer,
        max_epochs=max_epochs,
        seed=seed,
    )


def map_scarches_query(
    query: Any,
    *,
    model_dir: str | Path,
    output_dir: str | Path,
    max_epochs: int = 100,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Perform scArches surgery and export probabilities plus uncertainty."""

    root = Path(output_dir)
    unavailable = _gpu_optional_status("scarches_mapping", "scvi", "reference-models", root)
    if unavailable is not None:
        return pd.DataFrame(), unavailable
    if not Path(model_dir).is_dir():
        raise ContractError(f"scArches model directory does not exist: {model_dir}")
    try:
        import scvi

        model = scvi.model.SCANVI.load_query_data(query, model_dir)
        model.train(max_epochs=max_epochs, accelerator="gpu", devices=1)
        probabilities = model.predict(soft=True)
        if not isinstance(probabilities, pd.DataFrame):
            probabilities = pd.DataFrame(probabilities, index=query.obs_names)
        best = probabilities.idxmax(axis=1).astype(str)
        confidence = probabilities.max(axis=1).astype(float)
        cells = pd.DataFrame(
            {
                "barcode": query.obs_names.astype(str),
                "predicted_label": best.to_numpy(),
                "uncertainty": (1.0 - confidence).to_numpy(),
                "provider": "scarches",
                "role": "hypothesis_not_verdict",
            }
        )
    except Exception as exc:
        return pd.DataFrame(), _publish_status(
            root,
            _status(
                "scarches_mapping",
                "failed",
                reason=f"scArches query mapping failed without fallback: {exc}",
                remediation="Verify model/query gene compatibility and the declared GPU runtime.",
            ),
        )
    root.mkdir(parents=True, exist_ok=True)
    path = root / "cell_predictions.tsv"
    atomic_publish(path, cells.to_csv(sep="\t", index=False).encode())
    return cells, _publish_status(
        root,
        _status("scarches_mapping", "complete", artifacts={"cells": str(path)}),
    )


def map_scarches_file(
    *,
    query_path: str | Path,
    model_dir: str | Path,
    output_dir: str | Path,
    modality: str = "",
    max_epochs: int = 100,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    return map_scarches_query(
        _read_object(query_path, modality=modality),
        model_dir=model_dir,
        output_dir=output_dir,
        max_epochs=max_epochs,
    )


def curated_human_reference_hypothesis(
    query: Any,
    reference: Any,
    *,
    reference_label_column: str,
    query_cluster_column: str = "",
    output_dir: str | Path | None = None,
    **kwargs: Any,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Azimuth-compatible equivalent: map a versioned curated human reference locally."""

    return reference_knn_hypothesis(
        query,
        reference,
        reference_label_column=reference_label_column,
        query_cluster_column=query_cluster_column,
        output_dir=output_dir,
        provider_name="curated_human_reference",
        **kwargs,
    )


def map_curated_human_reference_files(
    *,
    query_path: str | Path,
    reference_path: str | Path,
    reference_label_column: str,
    query_cluster_column: str = "",
    query_modality: str = "",
    reference_modality: str = "",
    output_dir: str | Path | None = None,
    **kwargs: Any,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    return curated_human_reference_hypothesis(
        _read_object(query_path, modality=query_modality),
        _read_object(reference_path, modality=reference_modality),
        reference_label_column=reference_label_column,
        query_cluster_column=query_cluster_column,
        output_dir=output_dir,
        **kwargs,
    )


def scimilarity_hypothesis(
    adata: Any,
    *,
    model_path: str | Path,
    output_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    root = Path(output_dir) if output_dir is not None else None
    if importlib.util.find_spec("scimilarity") is None:
        return pd.DataFrame(), _missing_dependency(
            "scimilarity", "scimilarity", "scimilarity", root
        )
    if not Path(model_path).is_dir():
        return pd.DataFrame(), _publish_status(
            root,
            _status(
                "scimilarity",
                "unavailable",
                reason=f"SCimilarity model bundle is absent: {model_path}",
                remediation="Provision a checksum-verified SCimilarity annotation bundle and pass its model directory.",
            ),
        )
    try:
        from scimilarity import CellAnnotation

        model = CellAnnotation(model_path=str(model_path))
        values = _matrix(adata)
        values = values.toarray() if sparse.issparse(values) else np.asarray(values)
        embeddings = model.get_embeddings(values)
        predictions = model.get_predictions_kNN(embeddings)
        if isinstance(predictions, tuple):
            labels = predictions[0]
            distances = np.asarray(predictions[1], dtype=float)
            uncertainty = distances / np.maximum(1.0, np.nanmax(distances))
        else:
            labels = predictions
            uncertainty = np.full(len(adata), np.nan)
        cells = pd.DataFrame(
            {
                "barcode": adata.obs_names.astype(str),
                "predicted_label": np.asarray(labels).astype(str),
                "uncertainty": uncertainty,
                "provider": "scimilarity",
                "role": "hypothesis_not_verdict",
            }
        )
    except Exception as exc:
        return pd.DataFrame(), _publish_status(
            root,
            _status(
                "scimilarity",
                "failed",
                reason=f"SCimilarity inference failed: {exc}",
                remediation="Verify the model bundle version, gene order, and scimilarity package compatibility.",
            ),
        )
    artifacts: dict[str, str] = {}
    if root is not None:
        root.mkdir(parents=True, exist_ok=True)
        path = root / "cell_predictions.tsv"
        atomic_publish(path, cells.to_csv(sep="\t", index=False).encode())
        artifacts["cells"] = str(path)
    return cells, _publish_status(root, _status("scimilarity", "complete", artifacts=artifacts))


def map_scimilarity_file(
    *,
    input_path: str | Path,
    model_path: str | Path,
    modality: str = "",
    output_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    return scimilarity_hypothesis(
        _read_object(input_path, modality=modality),
        model_path=model_path,
        output_dir=output_dir,
    )


def _acquire_census_worker(
    result_queue: Any,
    *,
    destination: str,
    organism: str,
    tissue: str,
    census_version: str,
    max_cells: int,
    label_columns: tuple[str, ...],
) -> None:
    """Execute optional Census I/O in an isolated timeout-capable process."""

    try:
        import cellxgene_census

        filter_text = _census_value_filter(tissue)
        with cellxgene_census.open_soma(census_version=census_version) as census:
            obs = cellxgene_census.get_obs(
                census,
                organism,
                value_filter=filter_text,
                column_names=["soma_joinid", *label_columns],
            )
            obs = obs.loc[_nonblank_label_mask(obs, label_columns)]
            if obs.empty:
                raise ContractError(
                    "Census returned no cells with complete nonblank labels for "
                    f"organism={organism!r}, tissue={tissue!r}"
                )
            joinids = obs["soma_joinid"].iloc[:max_cells].to_numpy()
            adata = cellxgene_census.get_anndata(
                census,
                organism,
                obs_coords=joinids,
                obs_column_names=list(label_columns),
            )
        if not _nonblank_label_mask(adata.obs, label_columns).all():
            raise ContractError("Census query returned blank values in requested label columns")
        adata.write_h5ad(destination)
        result_queue.put({"status": "complete", "n_cells": int(adata.n_obs)})
    except Exception as exc:  # subprocess must preserve the exact provider failure
        result_queue.put({"status": "failed", "reason": str(exc)})


def fetch_census_reference(
    *,
    organism: str,
    tissue: str,
    output_path: str | Path,
    census_version: str,
    max_cells: int = 50_000,
    label_columns: Iterable[str] = ("cell_type",),
    expected_sha256: str = "",
    timeout_seconds: float = 300.0,
) -> dict[str, Any]:
    """Acquire a bounded CELLxGENE Census slice with a versioned cache receipt."""

    destination = Path(output_path)
    root = destination.parent
    organism = organism.strip()
    tissue = tissue.strip()
    census_version = census_version.strip()
    labels = _normalized_census_labels(label_columns)
    if not organism:
        raise ContractError("Census organism must be nonblank")
    if not census_version:
        raise ContractError("Census census_version must be nonblank")
    tissue_filter = _census_value_filter(tissue)
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ContractError("Census timeout_seconds must be finite and positive")
    if expected_sha256 and (
        len(expected_sha256) != 64
        or any(character not in "0123456789abcdef" for character in expected_sha256)
    ):
        raise ContractError("Census expected_sha256 must be lowercase SHA-256 hex")
    if max_cells < 1 or max_cells > 1_000_000:
        raise ContractError("Census max_cells must be between 1 and 1,000,000")

    request_metadata: dict[str, Any] = {
        "organism": organism,
        "tissue": tissue,
        "tissue_filter": tissue_filter,
        "census_version": census_version,
        "max_cells": max_cells,
        "label_columns": list(labels),
        "timeout_seconds": timeout_seconds,
        "expected_sha256": expected_sha256,
    }

    def publish_census(
        value: dict[str, Any],
        *,
        cache_hit: bool,
        checksum_validated: bool,
    ) -> dict[str, Any]:
        value.update(request_metadata)
        value.update(
            {
                "cache_hit": cache_hit,
                "checksum_validated": checksum_validated,
            }
        )
        return _publish_status(root, value)

    if destination.is_file():
        digest = hashlib.sha256(destination.read_bytes()).hexdigest()
        if expected_sha256 and digest != expected_sha256:
            return publish_census(
                _status(
                    "cellxgene_census",
                    "failed",
                    reason="cached Census reference checksum does not match the declared SHA-256",
                    remediation="Remove or replace the untrusted cache entry and retry with the pinned reference receipt.",
                    artifacts={"untrusted_reference": str(destination), "sha256": digest},
                ),
                cache_hit=True,
                checksum_validated=False,
            )
        value = _status(
            "cellxgene_census",
            "complete",
            artifacts={"reference": str(destination), "sha256": digest},
        )
        return publish_census(
            value,
            cache_hit=True,
            checksum_validated=bool(expected_sha256),
        )
    if importlib.util.find_spec("cellxgene_census") is None:
        return publish_census(
            _missing_dependency("cellxgene_census", "cellxgene-census", "census", None),
            cache_hit=False,
            checksum_validated=False,
        )
    root.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(f".{destination.name}.partial.{os.getpid()}")
    context = multiprocessing.get_context("fork" if hasattr(os, "fork") else "spawn")
    result_queue = context.Queue()
    process = cast(_ProcessFactoryContext, context).Process(
        target=_acquire_census_worker,
        kwargs={
            "result_queue": result_queue,
            "destination": str(partial),
            "organism": organism,
            "tissue": tissue,
            "census_version": census_version,
            "max_cells": max_cells,
            "label_columns": labels,
        },
    )
    process.start()
    process.join(timeout_seconds)
    if process.is_alive():
        process.terminate()
        process.join()
        partial.unlink(missing_ok=True)
        result_queue.close()
        return publish_census(
            _status(
                "cellxgene_census",
                "failed",
                reason=f"bounded Census acquisition exceeded its {timeout_seconds:g}-second timeout",
                remediation="Increase --timeout-seconds only after checking network access and the bounded query size.",
            ),
            cache_hit=False,
            checksum_validated=False,
        )
    try:
        worker = result_queue.get(timeout=1.0)
    except queue.Empty:
        worker = {
            "status": "failed",
            "reason": f"Census worker exited with code {process.exitcode} without a receipt",
        }
    finally:
        result_queue.close()
    if worker.get("status") != "complete" or not partial.is_file():
        partial.unlink(missing_ok=True)
        return publish_census(
            _status(
                "cellxgene_census",
                "failed",
                reason=f"bounded Census acquisition failed: {worker.get('reason', 'no output')}",
                remediation="Check the pinned Census version, organism/tissue filter, network, and cache permissions.",
            ),
            cache_hit=False,
            checksum_validated=False,
        )
    digest = hashlib.sha256(partial.read_bytes()).hexdigest()
    if expected_sha256 and digest != expected_sha256:
        partial.unlink(missing_ok=True)
        return publish_census(
            _status(
                "cellxgene_census",
                "failed",
                reason="downloaded Census reference checksum does not match the declared SHA-256",
                remediation="Do not use the cache entry; verify the pinned Census version and expected checksum.",
                artifacts={"sha256": digest},
            ),
            cache_hit=False,
            checksum_validated=False,
        )
    os.replace(partial, destination)
    value = _status(
        "cellxgene_census",
        "complete",
        artifacts={"reference": str(destination), "sha256": digest},
    )
    value.update({"n_cells": int(worker["n_cells"])})
    return publish_census(
        value,
        cache_hit=False,
        checksum_validated=bool(expected_sha256),
    )


__all__ = [
    "celltypist_hypothesis",
    "curated_human_reference_hypothesis",
    "discover_celltypist_models",
    "existing_prediction_file",
    "existing_prediction_hypothesis",
    "fetch_census_reference",
    "harmonize_reference_genes",
    "load_symbol_map",
    "map_curated_human_reference_files",
    "map_reference_files",
    "map_scarches_file",
    "map_scarches_query",
    "map_scimilarity_file",
    "predict_celltypist_file",
    "read_single_cell_object",
    "reference_knn_hypothesis",
    "scimilarity_hypothesis",
    "train_scanvi_reference",
    "train_scanvi_reference_file",
]

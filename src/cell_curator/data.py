"""Input inspection, matrix contracts, and immutable canonical working copies."""

from __future__ import annotations

import json
import os
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from .config import CellCuratorConfig, RepresentationConfig
from .contracts import ContractError, sha256, values_sha256
from .provenance import publish_json

SYMBOL_COLUMN_CANDIDATES = (
    "gene_symbol",
    "gene_symbols",
    "symbol",
    "gene_name",
    "feature_name",
)
GENE_MAPPING_SOURCE_COLUMNS = (
    "ensembl_id",
    "gene_id",
    "feature_id",
    "source_id",
    "id",
)
GENE_MAPPING_SYMBOL_COLUMNS = (
    "gene_symbol",
    "symbol",
    "gene_name",
    "target_symbol",
)
GENE_LIKE_ASSAYS = frozenset({"rna", "gene_activity", "spatial"})
ENSEMBL_PATTERN = re.compile(r"^ENS[A-Z0-9]*G[0-9]+(?:\.[0-9]+)?$", re.IGNORECASE)
S_PHASE_GENES = frozenset(
    {
        "MCM5",
        "PCNA",
        "TYMS",
        "FEN1",
        "MCM2",
        "MCM4",
        "RRM1",
        "UNG",
        "GINS2",
        "MCM6",
        "CDCA7",
        "DTL",
        "PRIM1",
        "UHRF1",
        "CENPU",
        "HELLS",
        "RFC2",
        "RPA2",
        "NASP",
        "RAD51AP1",
        "GMNN",
        "WDR76",
        "SLBP",
        "CCNE2",
        "UBR7",
        "POLD3",
        "MSH2",
        "ATAD2",
        "RAD51",
        "RRM2",
        "CDC45",
        "CDC6",
        "EXO1",
        "TIPIN",
        "DSCC1",
        "BLM",
        "CASP8AP2",
        "USP1",
        "CLSPN",
        "POLA1",
        "CHAF1B",
        "BRIP1",
        "E2F8",
    }
)
G2M_PHASE_GENES = frozenset(
    {
        "HMGB2",
        "CDK1",
        "NUSAP1",
        "UBE2C",
        "BIRC5",
        "TPX2",
        "TOP2A",
        "NDC80",
        "CKS2",
        "NUF2",
        "CKS1B",
        "MKI67",
        "TMPO",
        "CENPF",
        "TACC3",
        "FAM64A",
        "SMC4",
        "CCNB2",
        "CKAP2L",
        "CKAP2",
        "AURKB",
        "BUB1",
        "KIF11",
        "ANP32E",
        "TUBB4B",
        "GTSE1",
        "KIF20B",
        "HJURP",
        "CDCA3",
        "HN1",
        "CDC20",
        "TTK",
        "CDC25C",
        "KIF2C",
        "RANGAP1",
        "NCAPD2",
        "DLGAP5",
        "CDCA2",
        "CDCA8",
        "ECT2",
        "KIF23",
        "HMMR",
        "AURKA",
        "PSRC1",
        "ANLN",
        "LBR",
        "CKAP5",
        "CENPE",
        "CTCF",
        "NEK2",
        "G2E3",
        "GAS2L3",
        "CBX5",
        "CENPA",
    }
)
QC_PHASE_COLORS = ("#8c8c8c", "#4c78a8", "#e45756")


@dataclass(frozen=True)
class LoadedInput:
    root: Any
    assays: dict[str, Any]
    path: Path
    kind: str


def load_input(config: CellCuratorConfig, *, backed: bool = False) -> LoadedInput:
    path = Path(config.input.path)
    if not path.is_file():
        raise ContractError(f"input object does not exist: {path}")
    if config.input.kind == "h5ad":
        import anndata as ad

        root = ad.read_h5ad(path, backed="r" if backed else None)
        if len(config.input.assays) != 1:
            close_input(LoadedInput(root, {}, path, "h5ad"))
            raise ContractError("h5ad input must declare exactly one root assay")
        assay_id, assay = next(iter(config.input.assays.items()))
        if assay.modality not in {"", "root"}:
            close_input(LoadedInput(root, {}, path, "h5ad"))
            raise ContractError("h5ad assay modality must be 'root'")
        return LoadedInput(root, {assay_id: root}, path, "h5ad")
    import mudata as md

    root = md.read_h5mu(path, backed=backed)
    assays: dict[str, Any] = {}
    for assay_id, assay in config.input.assays.items():
        if not assay.modality or assay.modality not in root.mod:
            close_input(LoadedInput(root, assays, path, "h5mu"))
            raise ContractError(
                f"MuData modality for assay {assay_id!r} is absent: {assay.modality!r}"
            )
        assays[assay_id] = root.mod[assay.modality]
    return LoadedInput(root, assays, path, "h5mu")


def close_input(loaded: LoadedInput) -> None:
    file = getattr(loaded.root, "file", None)
    if file is not None:
        try:
            file.close()
        except (AttributeError, OSError):
            pass


def _nonblank(values: Sequence[object] | pd.Index | pd.Series) -> pd.Index:
    return pd.Index([str(value).strip() for value in values])


def _unversioned_gene_id(value: object) -> str:
    text = str(value).strip()
    return text.split(".", 1)[0] if ENSEMBL_PATTERN.fullmatch(text) else text


def detect_gene_id_convention(values: Sequence[object] | pd.Index | pd.Series) -> str:
    """Classify a feature identifier collection without guessing a mapping.

    The result is one of ``ensembl``, ``gene_symbol``, ``mixed``, or ``unknown``.
    Versioned Ensembl identifiers are recognized. Peak coordinates, antibody tags,
    and blank identifiers remain ``unknown`` rather than being relabeled as genes.
    """

    identifiers = _nonblank(values)
    identifiers = identifiers[identifiers != ""]
    if identifiers.empty:
        return "unknown"
    ensembl = identifiers.to_series(index=range(len(identifiers))).str.fullmatch(
        ENSEMBL_PATTERN
    )
    ensembl_fraction = float(ensembl.mean())
    if ensembl_fraction >= 0.9:
        return "ensembl"
    symbol_like = identifiers.to_series(index=range(len(identifiers))).str.fullmatch(
        r"[A-Za-z][A-Za-z0-9_.-]*"
    )
    symbol_fraction = float(symbol_like.mean())
    if ensembl_fraction > 0.05 and symbol_fraction >= 0.9:
        return "mixed"
    if symbol_fraction >= 0.9 and ensembl_fraction == 0:
        return "gene_symbol"
    return "unknown"


def detect_symbol_column(adata: Any) -> str | None:
    """Return the first unambiguous conventional gene-symbol column, if present."""

    columns = {str(column).casefold(): str(column) for column in adata.var.columns}
    for candidate in SYMBOL_COLUMN_CANDIDATES:
        column = columns.get(candidate.casefold())
        if column is None:
            continue
        values = adata.var[column]
        normalized = values.astype("string").fillna("").str.strip()
        if normalized.ne("").all():
            return column
    return None


def _mapping_column(
    columns: Sequence[object],
    requested: str | None,
    candidates: Sequence[str],
    *,
    role: str,
) -> str:
    available = {str(column).casefold(): str(column) for column in columns}
    if requested:
        selected = available.get(requested.casefold())
        if selected is None:
            raise ContractError(f"gene mapping {role} column is absent: {requested}")
        return selected
    for candidate in candidates:
        if candidate.casefold() in available:
            return available[candidate.casefold()]
    raise ContractError(
        f"gene mapping cannot detect its {role} column; expected one of {list(candidates)}"
    )


def load_gene_mapping(
    mapping: str | Path | pd.DataFrame | Mapping[object, object],
    *,
    source_column: str | None = None,
    symbol_column: str | None = None,
) -> dict[str, str]:
    """Load a validated feature-ID-to-gene-symbol mapping.

    Mappings may be a dictionary, a DataFrame, or a CSV/TSV path. Conflicting
    duplicate IDs and blank symbols fail closed. Both versioned and unversioned
    Ensembl lookup keys are registered when they agree.
    """

    if isinstance(mapping, Mapping):
        pairs = [(str(source), str(symbol)) for source, symbol in mapping.items()]
    else:
        if isinstance(mapping, pd.DataFrame):
            frame = mapping.copy()
        else:
            path = Path(mapping)
            if not path.is_file():
                raise ContractError(f"gene mapping does not exist: {path}")
            separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
            frame = pd.read_csv(path, sep=separator, dtype="string")
        source = _mapping_column(
            frame.columns,
            source_column,
            GENE_MAPPING_SOURCE_COLUMNS,
            role="source",
        )
        symbol = _mapping_column(
            frame.columns,
            symbol_column,
            GENE_MAPPING_SYMBOL_COLUMNS,
            role="symbol",
        )
        pairs = list(zip(frame[source], frame[symbol], strict=True))
    result: dict[str, str] = {}
    for raw_source, raw_symbol in pairs:
        source = str(raw_source).strip()
        symbol = str(raw_symbol).strip()
        missing_tokens = {"nan", "<na>", "none"}
        if (
            not source
            or not symbol
            or source.casefold() in missing_tokens
            or symbol.casefold() in missing_tokens
        ):
            raise ContractError("gene mapping contains a blank source ID or symbol")
        for key in dict.fromkeys((source, _unversioned_gene_id(source))):
            existing = result.get(key)
            if existing is not None and existing != symbol:
                raise ContractError(
                    f"gene mapping contains conflicting symbols for {key!r}: "
                    f"{existing!r} and {symbol!r}"
                )
            result[key] = symbol
    if not result:
        raise ContractError("gene mapping contains no usable rows")
    return result


def gene_symbols(
    adata: Any,
    symbol_column: str | None = None,
    mapping: str | Path | pd.DataFrame | Mapping[object, object] | None = None,
) -> pd.Index:
    """Return harmonized gene symbols in the object's feature order.

    An explicit mapping maps ``var_names`` (including version-stripped Ensembl
    IDs). Otherwise an explicit or automatically detected symbol column is used.
    Symbol-like ``var_names`` are accepted directly. Incomplete mappings and bare
    Ensembl IDs fail with an actionable error instead of silently mixing ID spaces.
    """

    identifiers = _nonblank(adata.var_names)
    detected_column = symbol_column or detect_symbol_column(adata)
    fallback: pd.Index | None = None
    if detected_column is not None:
        if detected_column not in adata.var:
            raise ContractError(f"gene-symbol column is absent: {detected_column}")
        fallback = _nonblank(adata.var[detected_column])
        if (fallback == "").any():
            raise ContractError(f"gene-symbol column contains blank values: {detected_column}")
    if mapping is not None:
        lookup = load_gene_mapping(mapping)
        resolved = []
        missing = []
        for position, identifier in enumerate(identifiers):
            symbol = lookup.get(identifier) or lookup.get(_unversioned_gene_id(identifier))
            if symbol is None and fallback is not None:
                symbol = str(fallback[position])
            if symbol is None or not str(symbol).strip():
                missing.append(identifier)
                resolved.append("")
            else:
                resolved.append(str(symbol).strip())
        if missing:
            raise ContractError(
                "gene mapping is incomplete for "
                f"{len(missing)} features (examples: {missing[:5]})"
            )
        return pd.Index(resolved, name="gene_symbol")
    if fallback is not None:
        return pd.Index(fallback, name="gene_symbol")
    convention = detect_gene_id_convention(identifiers)
    if convention == "gene_symbol":
        return pd.Index(identifiers, name="gene_symbol")
    if convention == "ensembl":
        raise ContractError(
            "features use Ensembl identifiers but no gene-symbol column or mapping was supplied"
        )
    raise ContractError(
        "gene symbols are unavailable; supply a symbol column or an ID-to-symbol mapping"
    )


def _declaration_value(declaration: Any, name: str, default: Any = "") -> Any:
    if isinstance(declaration, Mapping):
        return declaration.get(name, default)
    return getattr(declaration, name, default)


def matrix_for(adata: Any, declaration: RepresentationConfig | dict[str, Any]) -> Any:
    source = str(_declaration_value(declaration, "source")).lower()
    key = str(_declaration_value(declaration, "key"))
    if source == "x":
        return adata.X
    if source == "raw":
        if adata.raw is None:
            raise ContractError("declared .raw representation is absent")
        return adata.raw.X
    if source != "layer":
        raise ContractError(f"matrix source must be X, raw, or layer, observed {source!r}")
    if key not in adata.layers:
        raise ContractError(f"declared layer is absent: {key}")
    return adata.layers[key]


def dense(matrix: Any) -> np.ndarray:
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


def matrix_minimum(matrix: Any) -> float:
    if sparse.issparse(matrix):
        stored = float(matrix.data.min()) if matrix.nnz else 0.0
        return min(stored, 0.0)
    values = np.asarray(matrix)
    return float(np.nanmin(values)) if values.size else 0.0


def matrix_summary(matrix: Any) -> dict[str, Any]:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    observed = np.asarray(values)
    finite = observed[np.isfinite(observed)]
    if finite.size == 0:
        raise ContractError("declared representation contains no finite values")
    sample = (
        finite
        if finite.size <= 100_000
        else finite[np.linspace(0, finite.size - 1, 100_000).astype(int)]
    )
    integer_fraction = float(np.mean(np.isclose(sample, np.round(sample), atol=1e-6)))
    total_values = int(np.prod(matrix.shape))
    nonfinite = int(observed.size - finite.size)
    implicit_zeros = total_values - observed.size if sparse.issparse(matrix) else 0
    return {
        "shape": [int(item) for item in matrix.shape],
        "minimum": float(sample.min()),
        "maximum": float(sample.max()),
        "mean": float(sample.mean()),
        "integer_fraction": integer_fraction,
        "nonnegative": bool(sample.min() >= 0),
        "finite_fraction": float((finite.size + implicit_zeros) / max(1, total_values)),
        "nonfinite_values": nonfinite,
        "zero_fraction": float(
            (implicit_zeros + np.count_nonzero(np.isclose(finite, 0.0))) / max(1, total_values)
        ),
    }


def normalization_fingerprint(
    matrix: Any,
    *,
    assay_kind: str,
    declared_name: str = "",
    role: str = "evidence",
) -> dict[str, Any]:
    """Classify matrix semantics and fail closed on incompatible declarations.

    Counts, log-normalized values, nonnegative corrected values, signed corrected
    values, and detection-safe matrices are distinct classes. A declaration name
    (for example ``logcounts`` or ``counts``) resolves otherwise unavoidable
    low-range ambiguity; declarations contradicted by their numeric fingerprint
    are invalid.
    """

    summary = matrix_summary(matrix)
    verdict = "valid"
    reasons: list[str] = []
    declared = declared_name.casefold()
    expects_log = any(token in declared for token in ("log", "cpm", "tpm", "normalized"))
    expects_counts = not expects_log and any(
        token in declared for token in ("count", "umi", "raw")
    )
    expects_corrected = any(
        token in declared
        for token in ("correct", "residual", "scaled", "integrated", "arcsinh", "tfidf")
    )
    integer_like = summary["integer_fraction"] >= 0.995
    binary_like = integer_like and summary["minimum"] >= 0 and summary["maximum"] <= 1
    if summary["finite_fraction"] < 1.0:
        verdict = "invalid"
        reasons.append("matrix contains non-finite values")
    if not summary["nonnegative"]:
        classification = "corrected_signed"
        detection_safe = False
        if role == "detection" or expects_counts or expects_log:
            verdict = "invalid"
            reasons.append("declared counts/log/detection matrix contains negative values")
    elif role == "detection":
        classification = "detection_safe"
        detection_safe = True
    elif binary_like:
        classification = "binary_or_counts"
        detection_safe = True
        if expects_log:
            # A declared transformation is sufficient to resolve a low-range
            # all-zero/one matrix that is numerically indistinguishable from counts.
            classification = "log_normalized"
    elif integer_like:
        classification = "counts"
        detection_safe = True
        if expects_log and summary["maximum"] > 20:
            verdict = "invalid"
            reasons.append("declared log-normalized representation appears to be raw counts")
        elif expects_log:
            classification = "log_normalized"
        elif expects_corrected:
            classification = "corrected_nonnegative"
    else:
        detection_safe = True
        if expects_counts:
            classification = "ambiguous_noninteger_counts"
            verdict = "ambiguous" if verdict != "invalid" else verdict
            reasons.append("declared count representation is substantially non-integer")
        elif expects_corrected:
            classification = "corrected_nonnegative"
        elif summary["maximum"] <= 100:
            classification = "log_normalized"
        else:
            classification = "ambiguous_nonnegative"
            verdict = "ambiguous" if verdict != "invalid" else verdict
            reasons.append("nonnegative representation has an ambiguous scale above 100")
    if assay_kind in GENE_LIKE_ASSAYS and expects_log and summary["maximum"] > 100:
        verdict = "invalid"
        reasons.append("gene-like log-normalized representation has an implausible scale")
    return {
        **summary,
        "classification": classification,
        "verdict": verdict,
        "role": role,
        "detection_safe": detection_safe,
        "reasons": list(dict.fromkeys(reasons)),
    }


def categorical_covariate_inventory(
    obj: Any,
    *,
    maximum_categories: int = 50,
) -> dict[str, dict[str, Any]]:
    """Discover categorical covariates and preserve their declared color palettes."""

    result: dict[str, dict[str, Any]] = {}
    n_observations = max(1, len(obj.obs))
    for column in map(str, obj.obs.columns):
        values = obj.obs[column]
        nonmissing = values.dropna()
        n_categories = int(nonmissing.nunique())
        categorical_dtype = isinstance(values.dtype, pd.CategoricalDtype)
        candidate_dtype = (
            categorical_dtype
            or pd.api.types.is_bool_dtype(values.dtype)
            or pd.api.types.is_object_dtype(values.dtype)
            or pd.api.types.is_string_dtype(values.dtype)
            or (pd.api.types.is_integer_dtype(values.dtype) and n_categories <= 20)
        )
        if (
            not candidate_dtype
            or n_categories == 0
            or n_categories > maximum_categories
            or (n_categories == n_observations and not categorical_dtype)
        ):
            continue
        if categorical_dtype:
            categories = [str(value) for value in values.cat.categories]
            ordered = bool(values.cat.ordered)
        else:
            categories = sorted(map(str, nonmissing.unique()))
            ordered = False
        counts = values.astype("string").fillna("<missing>").value_counts(dropna=False)
        color_key = f"{column}_colors"
        raw_colors = obj.uns.get(color_key)
        colors = [str(color) for color in raw_colors] if raw_colors is not None else []
        result[column] = {
            "dtype": str(values.dtype),
            "n_categories": n_categories,
            "categories": categories,
            "counts": {str(key): int(value) for key, value in counts.items()},
            "ordered": ordered,
            "color_key": color_key if raw_colors is not None else None,
            "colors": colors,
            "colors_match_categories": bool(colors) and len(colors) == len(categories),
        }
    return result


def targeted_panel_summary(
    feature_ids: Sequence[object] | pd.Index,
    *,
    assay_kind: str,
    gene_panel_threshold: int = 1_000,
    declared_targeted: bool = False,
) -> dict[str, Any]:
    """Return conservative panel-aware capability metadata."""

    identifiers = _nonblank(feature_ids)
    n_features = len(identifiers)
    gene_like = assay_kind in GENE_LIKE_ASSAYS
    antibody_panel = assay_kind == "protein"
    targeted = (
        declared_targeted
        or antibody_panel
        or (gene_like and n_features <= gene_panel_threshold)
    )
    if antibody_panel:
        panel_type = "antibody_panel"
    elif targeted and gene_like:
        panel_type = "targeted_gene_panel"
    elif gene_like:
        panel_type = "whole_transcriptome_or_broad_gene_space"
    else:
        panel_type = "feature_space"
    return {
        "is_targeted_panel": targeted,
        "detection_basis": (
            "declared_technology"
            if declared_targeted
            else "antibody_modality"
            if antibody_panel
            else "feature_count_threshold"
            if targeted
            else "broad_feature_space"
        ),
        "panel_type": panel_type,
        "n_features": n_features,
        "gene_id_convention": detect_gene_id_convention(identifiers)
        if gene_like
        else "not_gene",
        "supports_program_scoring": bool(n_features),
        "supports_genome_wide_marker_discovery": bool(gene_like and not targeted),
        "requires_missing_feature_accounting": targeted,
        "requires_program_coverage_gating": targeted,
    }


def _gene_mapping_path(declaration: Any) -> str:
    for field in (
        "gene_mapping_path",
        "gene_id_mapping_path",
        "symbol_mapping_path",
    ):
        value = str(_declaration_value(declaration, field, "") or "").strip()
        if value:
            return value
    return ""


def _feature_inventory(adata: Any, declaration: Any) -> tuple[pd.Index, dict[str, Any]]:
    configured_column = str(_declaration_value(declaration, "feature_id_column", "") or "")
    if configured_column:
        if configured_column not in adata.var:
            raise ContractError(f"feature ID column is absent: {configured_column}")
        identifiers = _nonblank(adata.var[configured_column])
        if adata.var[configured_column].isna().any() or (identifiers == "").any():
            raise ContractError("feature IDs contain missing or blank values")
        if identifiers.duplicated().any():
            raise ContractError("feature IDs are ambiguous")
    else:
        identifiers = _nonblank(adata.var_names)
        if (identifiers == "").any() or identifiers.duplicated().any():
            raise ContractError("feature names contain blank or duplicated identifiers")
    convention = detect_gene_id_convention(identifiers)
    symbol_column = detect_symbol_column(adata)
    mapping_path = _gene_mapping_path(declaration)
    symbols: pd.Index | None = None
    symbol_status = "unavailable"
    symbol_reason = ""
    try:
        symbols = gene_symbols(
            adata,
            symbol_column=symbol_column,
            mapping=mapping_path or None,
        )
        symbol_status = (
            "mapped" if mapping_path else ("column" if symbol_column else "var_names")
        )
    except ContractError as exc:
        if mapping_path:
            raise ContractError(f"declared gene mapping is invalid: {exc}") from exc
        symbol_reason = str(exc)
    return identifiers, {
        "configured_id_column": configured_column or None,
        "id_convention": convention,
        "n_unique_ids": int(identifiers.nunique()),
        "symbol_column": symbol_column,
        "symbol_status": symbol_status,
        "symbol_mapping_path": mapping_path or None,
        "symbol_mapping_available": symbols is not None,
        "symbol_mapping_unique": bool(symbols is not None and symbols.is_unique),
        "symbol_mapping_reason": symbol_reason,
        "symbols": symbols,
    }


def _obsm_inventory(obj: Any) -> dict[str, list[int]]:
    return {str(key): [int(item) for item in value.shape] for key, value in obj.obsm.items()}


def _obs_frame(loaded: LoadedInput, config: CellCuratorConfig) -> pd.DataFrame:
    obs = loaded.root.obs.copy()
    keys = config.input.keys
    required = {keys.canonical_parent, keys.scope, *config.input.required_obs_keys}
    required.update(config.input.guide_columns)
    required.update(config.input.technical_metrics.values())
    missing = {key for key in required if key and key not in obs.columns}
    if missing:
        raise ContractError(f"input observations miss keys: {sorted(missing)}")
    if not obs.index.is_unique:
        raise ContractError("input barcodes are not unique")
    if obs[keys.canonical_parent].isna().any() or obs[keys.scope].isna().any():
        raise ContractError("canonical parent or scope contains missing values")
    return obs


def inspect(config: CellCuratorConfig) -> dict[str, Any]:
    loaded = load_input(config, backed=True)
    try:
        observed_sha = sha256(loaded.path)
        if config.input.expected_sha256 and observed_sha != config.input.expected_sha256:
            raise ContractError("input SHA-256 does not match the frozen declaration")
        obs = _obs_frame(loaded, config)
        root_barcodes = list(map(str, loaded.root.obs_names))
        assays: dict[str, Any] = {}
        blockers: list[str] = []
        root_obsm = _obsm_inventory(loaded.root)
        available_obsm: dict[str, list[int]] = dict(root_obsm)
        for assay_id, adata in loaded.assays.items():
            declaration = config.input.assays[assay_id]
            if list(map(str, adata.obs_names)) != root_barcodes:
                raise ContractError(f"assay {assay_id!r} barcode order differs from root")
            try:
                feature_ids, feature_inventory = _feature_inventory(adata, declaration)
            except ContractError as exc:
                raise ContractError(f"assay {assay_id!r} {exc}") from exc
            if (
                declaration.kind in GENE_LIKE_ASSAYS
                and feature_inventory["id_convention"] == "mixed"
                and not feature_inventory["symbol_mapping_available"]
            ):
                blockers.append(
                    f"{assay_id}: mixed gene identifier conventions require a symbol column "
                    "or supplied mapping"
                )
            graph_representations: list[dict[str, Any]] = []
            for graph in declaration.graph_representations:
                source = graph.source.lower()
                if source == "obsm":
                    if not graph.key or graph.key not in adata.obsm:
                        raise ContractError(
                            f"assay {assay_id!r} graph representation is absent: {graph.key!r}"
                        )
                    graph_matrix = adata.obsm[graph.key]
                elif source == "layer":
                    if not graph.key or graph.key not in adata.layers:
                        raise ContractError(
                            f"assay {assay_id!r} graph layer is absent: {graph.key!r}"
                        )
                    graph_matrix = adata.layers[graph.key]
                elif source == "raw":
                    if adata.raw is None:
                        raise ContractError(
                            f"assay {assay_id!r} graph .raw representation is absent"
                        )
                    graph_matrix = adata.raw.X
                else:
                    graph_matrix = adata.X
                if int(graph_matrix.shape[0]) != int(adata.n_obs):
                    raise ContractError(
                        f"assay {assay_id!r} graph representation has the wrong cell axis"
                    )
                graph_representations.append(
                    {
                        "source": graph.source,
                        "key": graph.key,
                        "shape": [int(item) for item in graph_matrix.shape],
                    }
                )
            available_obsm.update(_obsm_inventory(adata))
            representations: dict[str, Any] = {}
            for name, representation in declaration.representations.items():
                expression = matrix_for(adata, representation)
                representation_source = str(
                    _declaration_value(representation, "source", "")
                ).lower()
                expected_features = (
                    int(adata.raw.n_vars)
                    if representation_source == "raw" and adata.raw is not None
                    else int(adata.n_vars)
                )
                if tuple(expression.shape) != (int(adata.n_obs), expected_features):
                    raise ContractError(
                        f"assay {assay_id!r}/{name!r} matrix does not match the assay shape"
                    )
                fingerprint = normalization_fingerprint(
                    expression,
                    assay_kind=declaration.kind,
                    declared_name=(
                        f"{name} "
                        f"{_declaration_value(representation, 'key', '')} "
                        f"{representation_source}"
                    ),
                )
                if fingerprint["verdict"] in {"invalid", "ambiguous"}:
                    blockers.extend(
                        f"{assay_id}/{name}: {reason}" for reason in fingerprint["reasons"]
                    )
                detection_summary: dict[str, Any] = {
                    "available": False,
                    "required": False,
                    "nonnegative": False,
                }
                if representation.detection is not None:
                    detection = matrix_for(adata, representation.detection.model_dump())
                    if tuple(detection.shape) != tuple(expression.shape):
                        raise ContractError(
                            f"assay {assay_id!r}/{name!r} detection matrix shape differs"
                        )
                    detection_fingerprint = normalization_fingerprint(
                        detection,
                        assay_kind=declaration.kind,
                        declared_name=(
                            f"{representation.detection.key} "
                            f"{representation.detection.source} detection"
                        ),
                        role="detection",
                    )
                    detection_summary = {
                        **detection_fingerprint,
                        "available": True,
                        "required": representation.detection.required,
                        "semantics": representation.detection.semantics,
                        "source": representation.detection.source,
                        "key": representation.detection.key,
                    }
                    if (
                        representation.detection.required
                        and detection_fingerprint["verdict"] != "valid"
                    ):
                        blockers.append(
                            f"{assay_id}/{name} requires a nonnegative detection source"
                        )
                representations[name] = {
                    "feature_space": representation.feature_space or declaration.feature_space,
                    "ranking_backend": representation.ranking_backend,
                    "ontology_compatible": representation.ontology_compatible,
                    "normalization": fingerprint,
                    "detection": detection_summary,
                }
            assays[assay_id] = {
                "kind": declaration.kind,
                "modality": declaration.modality,
                "feature_space": declaration.feature_space,
                "n_features": int(adata.n_vars),
                "feature_ids": {
                    key: value for key, value in feature_inventory.items() if key != "symbols"
                },
                "panel": targeted_panel_summary(
                    feature_ids,
                    assay_kind=declaration.kind,
                    declared_targeted=(
                        declaration.kind in GENE_LIKE_ASSAYS
                        and getattr(
                            getattr(getattr(config, "guidance", None), "context", None),
                            "technology",
                            "",
                        )
                        == "targeted_panel"
                    ),
                ),
                "representations": representations,
                "graph_representations": graph_representations,
                "obsm": sorted(map(str, adata.obsm.keys())),
                "obsm_shapes": _obsm_inventory(adata),
                "layers": sorted(map(str, adata.layers.keys())),
                "raw": {
                    "available": adata.raw is not None,
                    "shape": (
                        [int(adata.raw.n_obs), int(adata.raw.n_vars)]
                        if adata.raw is not None
                        else None
                    ),
                },
            }
        if config.input.keys.embedding and config.input.keys.embedding not in available_obsm:
            raise ContractError(
                f"declared frozen embedding is absent: {config.input.keys.embedding}"
            )
        if config.input.keys.embedding:
            embedding_shape = available_obsm[config.input.keys.embedding]
            if (
                embedding_shape[0] != len(obs)
                or len(embedding_shape) != 2
                or embedding_shape[1] < 2
            ):
                raise ContractError("declared embedding must be a cells-by-at-least-two matrix")
        if config.input.keys.spatial and config.input.keys.spatial not in available_obsm:
            raise ContractError(
                f"declared spatial coordinate representation is absent: {config.input.keys.spatial}"
            )
        if config.input.keys.spatial:
            spatial_shape = available_obsm[config.input.keys.spatial]
            if spatial_shape[0] != len(obs) or len(spatial_shape) != 2 or spatial_shape[1] < 2:
                raise ContractError(
                    "declared spatial coordinates must be cells-by-at-least-two"
                )
        guidance = getattr(config, "guidance", None)
        context = getattr(guidance, "context", None)
        visualization_key = str(
            getattr(context, "visualization_key", "")
            or getattr(config.input.keys, "visualization", "")
            or ""
        )
        if visualization_key and visualization_key not in available_obsm:
            raise ContractError(f"declared visualization is absent: {visualization_key}")
        expected = sum(
            item.expected_parent_count for item in config.canonical_partitions.values()
        )
        observed = int(obs[config.input.keys.canonical_parent].astype(str).nunique())
        if expected != observed:
            raise ContractError(
                f"declared canonical parent total ({expected}) differs from input ({observed})"
            )
        return {
            "schema_version": 3,
            "artifact_kind": "cell_curator_input_inspection",
            "input_path": str(loaded.path.resolve()),
            "input_sha256": observed_sha,
            "input_kind": config.input.kind,
            "n_cells": int(len(obs)),
            "canonical_parent_count": observed,
            "canonical_parent_counts": {
                str(key): int(value)
                for key, value in obs[config.input.keys.canonical_parent]
                .astype(str)
                .value_counts()
                .sort_index()
                .items()
            },
            "ordered_barcode_sha256": values_sha256(root_barcodes),
            "representations": {
                "obsm": available_obsm,
                "embedding": config.input.keys.embedding or None,
                "visualization": visualization_key or None,
                "spatial": config.input.keys.spatial or None,
            },
            "categorical_covariates": categorical_covariate_inventory(loaded.root),
            "assays": assays,
            "blockers": sorted(set(blockers)),
        }
    finally:
        close_input(loaded)


def _row_sum(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel().astype(float)


def _positive_features(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray((matrix > 0).sum(axis=1)).ravel().astype(int)
    return np.count_nonzero(np.asarray(matrix) > 0, axis=1).astype(int)


def _feature_sum(matrix: Any, positions: Sequence[int] | np.ndarray) -> np.ndarray:
    if len(positions) == 0:
        return np.zeros(matrix.shape[0], dtype=float)
    return np.asarray(matrix[:, list(positions)].sum(axis=1)).ravel().astype(float)


def _qc_matrix(
    adata: Any,
    declaration: Any,
) -> tuple[Any | None, Any | None, str, bool, list[str]]:
    candidates: list[tuple[int, Any, Any, str, bool]] = []
    limitations: list[str] = []
    for representation_name, representation in declaration.representations.items():
        detection = representation.detection
        if detection is None:
            continue
        try:
            matrix = matrix_for(adata, detection)
            fingerprint = normalization_fingerprint(
                matrix,
                assay_kind=declaration.kind,
                declared_name=f"{detection.key} {detection.semantics}",
                role="detection",
            )
        except ContractError as exc:
            limitations.append(f"{representation_name} detection unavailable: {exc}")
            continue
        if fingerprint["verdict"] != "valid":
            limitations.append(
                f"{representation_name} detection is invalid: "
                + "; ".join(fingerprint["reasons"])
            )
            continue
        holder = adata.raw if str(detection.source).lower() == "raw" else adata
        name = f"declared_detection:{representation_name}:{detection.source}:{detection.key}"
        count_tokens = set(
            re.findall(r"[a-z0-9]+", f"{detection.key} {detection.semantics}".casefold())
        )
        count_like = bool(count_tokens.intersection({"count", "counts", "umi", "umis", "raw"}))
        candidates.append((0 if count_like else 3, matrix, holder, name, count_like))
    for layer in sorted(map(str, adata.layers.keys())):
        if layer.casefold() not in {"counts", "count", "raw_counts", "umi", "umis"}:
            continue
        matrix = adata.layers[layer]
        if matrix_minimum(matrix) >= 0:
            candidates.append((0, matrix, adata, f"layer:{layer}", True))
    if adata.raw is not None and matrix_minimum(adata.raw.X) >= 0:
        raw_summary = matrix_summary(adata.raw.X)
        raw_count_like = raw_summary["integer_fraction"] >= 0.995
        candidates.append(
            (1 if raw_count_like else 4, adata.raw.X, adata.raw, "raw:X", raw_count_like)
        )
    if adata.X is not None and matrix_minimum(adata.X) >= 0:
        x_summary = matrix_summary(adata.X)
        x_count_like = x_summary["integer_fraction"] >= 0.995
        candidates.append((2 if x_count_like else 5, adata.X, adata, "X", x_count_like))
    if not candidates:
        limitations.append("no nonnegative matrix is available for read-only QC")
        return None, None, "unavailable", False, limitations
    _, matrix, holder, source, count_like = min(candidates, key=lambda item: (item[0], item[3]))
    return matrix, holder, source, count_like, limitations


def _cell_cycle_scores(
    matrix: Any,
    total: np.ndarray,
    symbols: pd.Index,
    *,
    count_scale: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]] | None:
    upper = np.asarray([str(symbol).upper() for symbol in symbols])
    s_positions = np.flatnonzero(np.isin(upper, list(S_PHASE_GENES)))
    g2m_positions = np.flatnonzero(np.isin(upper, list(G2M_PHASE_GENES)))
    if not len(s_positions) or not len(g2m_positions):
        return None

    def score(positions: np.ndarray) -> np.ndarray:
        selected = matrix[:, positions]
        selected = selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
        if not count_scale:
            return selected.mean(axis=1)
        scale = np.divide(10_000.0, total, out=np.zeros_like(total), where=total > 0)
        return np.log1p(selected * scale[:, None]).mean(axis=1)

    s_score = score(s_positions)
    g2m_score = score(g2m_positions)
    phase = np.full(len(total), "G1", dtype=object)
    phase[(s_score > g2m_score) & (s_score > 0)] = "S"
    phase[(g2m_score > s_score) & (g2m_score > 0)] = "G2M"
    return (
        s_score,
        g2m_score,
        phase,
        {"s_genes": int(len(s_positions)), "g2m_genes": int(len(g2m_positions))},
    )


def _assay_qc(adata: Any, declaration: Any) -> tuple[pd.DataFrame, dict[str, Any]]:
    matrix, feature_holder, source, count_like, limitations = _qc_matrix(adata, declaration)
    metrics = pd.DataFrame(index=adata.obs_names.astype(str))
    summary: dict[str, Any] = {
        "status": "unavailable" if matrix is None else "complete",
        "matrix_source": source,
        "count_scale": count_like,
        "limitations": limitations,
        "metrics": [],
    }
    if matrix is None or feature_holder is None:
        return metrics, summary
    if matrix.shape[0] != adata.n_obs:
        raise ContractError("QC matrix cell axis differs from the assay")
    total = _row_sum(matrix)
    detected = _positive_features(matrix)
    metrics["n_features_detected"] = detected
    summary["metrics"].append("n_features_detected")
    if count_like:
        metrics["total_counts"] = total
        metrics["n_features_by_counts"] = detected
        summary["metrics"].extend(["total_counts", "n_features_by_counts"])
        summary["total_counts"] = {
            "minimum": float(total.min()) if len(total) else float("nan"),
            "median": float(np.median(total)) if len(total) else float("nan"),
            "maximum": float(total.max()) if len(total) else float("nan"),
        }
    else:
        summary["limitations"].append(
            "count-scale depth and count-fraction metrics are unavailable; "
            "the selected QC matrix is nonnegative but not count-like"
        )
    if declaration.kind not in GENE_LIKE_ASSAYS:
        return metrics, summary
    try:
        symbols = gene_symbols(
            feature_holder,
            mapping=_gene_mapping_path(declaration) or None,
        )
    except ContractError as exc:
        summary["limitations"].append(f"gene-aware QC unavailable: {exc}")
        return metrics, summary
    upper = np.asarray([str(symbol).upper() for symbol in symbols])
    mitochondrial = np.flatnonzero(np.char.startswith(upper.astype(str), "MT-"))
    ribosomal = np.flatnonzero(
        np.char.startswith(upper.astype(str), "RPS")
        | np.char.startswith(upper.astype(str), "RPL")
    )
    denominator = np.where(total > 0, total, np.nan)
    if len(mitochondrial) and count_like:
        metrics["pct_counts_mt"] = 100.0 * _feature_sum(matrix, mitochondrial) / denominator
        summary["metrics"].append("pct_counts_mt")
    elif not len(mitochondrial):
        summary["limitations"].append("no mitochondrial features are measured")
    if len(ribosomal) and count_like:
        metrics["pct_counts_ribo"] = 100.0 * _feature_sum(matrix, ribosomal) / denominator
        summary["metrics"].append("pct_counts_ribo")
    elif not len(ribosomal):
        summary["limitations"].append("no ribosomal features are measured")
    cycle = _cell_cycle_scores(matrix, total, symbols, count_scale=count_like)
    if cycle is None:
        summary["limitations"].append("insufficient S and G2M markers for phase scoring")
    else:
        s_score, g2m_score, phase, cycle_counts = cycle
        metrics["s_score"] = s_score
        metrics["g2m_score"] = g2m_score
        metrics["phase"] = phase
        summary["metrics"].extend(["s_score", "g2m_score", "phase"])
        summary["cell_cycle_features"] = cycle_counts
    return metrics, summary


def _assign_column(frame: pd.DataFrame, name: str, values: Sequence[Any]) -> None:
    incoming = pd.Series(values, index=frame.index)
    if name in frame:
        existing = frame[name]
        if pd.api.types.is_numeric_dtype(incoming):
            same = np.allclose(
                pd.to_numeric(existing, errors="coerce"),
                pd.to_numeric(incoming, errors="coerce"),
                equal_nan=True,
            )
        else:
            same = (
                existing.astype("string")
                .fillna("")
                .equals(incoming.astype("string").fillna(""))
            )
        if not same:
            raise ContractError(
                f"canonical QC column already exists with different values: {name}"
            )
        return
    frame[name] = incoming


def _safe_identifier(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "assay"


def _prepare_canonical_copy(
    config: CellCuratorConfig,
    *,
    canonical_path: Path,
    run_root: Path,
    source_sha256: str,
) -> dict[str, Any]:
    if config.input.kind == "h5ad":
        import anndata as ad

        root = ad.read_h5ad(canonical_path)
        assay_id = next(iter(config.input.assays))
        assay_objects = {assay_id: root}
        writer = root.write_h5ad
    else:
        import mudata as md

        root = md.read_h5mu(canonical_path)
        assay_objects = {
            assay_id: root.mod[declaration.modality]
            for assay_id, declaration in config.input.assays.items()
        }
        writer = root.write_h5mu
    feature_summaries: dict[str, Any] = {}
    qc_summaries: dict[str, Any] = {}
    combined = pd.DataFrame({"barcode": root.obs_names.astype(str)})
    for assay_id, adata in assay_objects.items():
        declaration = config.input.assays[assay_id]
        identifiers, inventory = _feature_inventory(adata, declaration)
        adata.var["cell_curator_original_feature_id"] = adata.var_names.astype(str)
        adata.var["cell_curator_feature_id"] = identifiers.astype(str)
        symbols = inventory.pop("symbols")
        if symbols is not None:
            adata.var["cell_curator_gene_symbol"] = symbols.astype(str)
        feature_summaries[assay_id] = inventory
        metrics, qc_summary = _assay_qc(adata, declaration)
        qc_summaries[assay_id] = qc_summary
        safe_assay = _safe_identifier(assay_id)
        for metric in metrics:
            local_name = f"cell_curator_qc_{metric}"
            values = metrics[metric].to_numpy()
            if metric == "phase":
                values = pd.Categorical(values, categories=["G1", "S", "G2M"], ordered=True)
            _assign_column(adata.obs, local_name, values)
            root_name = f"cell_curator_qc__{safe_assay}__{metric}"
            aligned = metrics[metric].reindex(root.obs_names.astype(str)).to_numpy()
            _assign_column(root.obs, root_name, aligned)
            combined[f"{safe_assay}__{metric}"] = aligned
        if "phase" in metrics:
            adata.uns["cell_curator_qc_phase_colors"] = np.asarray(
                QC_PHASE_COLORS, dtype=object
            )
            root.uns[f"cell_curator_qc__{safe_assay}__phase_colors"] = np.asarray(
                QC_PHASE_COLORS,
                dtype=object,
            )
    root.uns["cell_curator_canonical"] = {
        "schema_version": 3,
        "source_sha256": source_sha256,
        "standardized_feature_columns": [
            "cell_curator_original_feature_id",
            "cell_curator_feature_id",
            "cell_curator_gene_symbol",
        ],
        "qc_column_prefix": "cell_curator_qc_",
    }
    partial = canonical_path.with_name(f".{canonical_path.name}.qc-partial.{os.getpid()}")
    try:
        writer(partial)
        os.replace(partial, canonical_path)
    finally:
        partial.unlink(missing_ok=True)
    from .provenance import atomic_publish

    qc_dir = run_root / "qc"
    table_path = qc_dir / "cell_metrics.tsv"
    atomic_publish(
        table_path,
        combined.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    manifest = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_read_only_qc",
        "status": "complete",
        "run_id": config.run.run_id,
        "source_sha256": source_sha256,
        "canonical_path": str(canonical_path),
        "canonical_sha256": sha256(canonical_path),
        "cell_metrics_path": str(table_path),
        "cell_metrics_sha256": sha256(table_path),
        "feature_standardization": feature_summaries,
        "assays": qc_summaries,
        "source_mutated": False,
    }
    manifest_path = qc_dir / "qc.manifest.json"
    publish_json(manifest_path, manifest)
    return {
        **manifest,
        "manifest_path": str(manifest_path),
        "manifest_sha256": sha256(manifest_path),
    }


def compute_qc(
    config: CellCuratorConfig,
    *,
    run_root: Path,
    canonical_path: Path | None = None,
    source_sha256: str | None = None,
) -> dict[str, Any]:
    """Compute or validate source-preserving QC on the canonical working copy.

    ``canonicalize`` calls the computational path while constructing a new run.
    Later calls are idempotent verification calls: hashes must still match the
    frozen source, canonical object, and durable per-cell QC table.
    """

    source = Path(config.input.path)
    observed_source_sha = sha256(source)
    if source_sha256 is not None and observed_source_sha != source_sha256:
        raise ContractError("source input changed before canonical QC")
    if canonical_path is not None:
        result = _prepare_canonical_copy(
            config,
            canonical_path=canonical_path,
            run_root=run_root,
            source_sha256=source_sha256 or observed_source_sha,
        )
        if sha256(source) != observed_source_sha:
            raise ContractError("source input changed during canonical QC")
        return result
    manifest_path = run_root / "qc" / "qc.manifest.json"
    if not manifest_path.is_file():
        raise ContractError("QC is not materialized; canonicalize the run first")
    manifest = json.loads(manifest_path.read_text())
    frozen_canonical = Path(str(manifest.get("canonical_path", "")))
    metrics = Path(str(manifest.get("cell_metrics_path", "")))
    if manifest.get("source_sha256") != observed_source_sha:
        raise ContractError("source input changed after canonical QC")
    if not frozen_canonical.is_file() or sha256(frozen_canonical) != manifest.get(
        "canonical_sha256"
    ):
        raise ContractError("canonical QC working copy is absent or hash-invalid")
    if not metrics.is_file() or sha256(metrics) != manifest.get("cell_metrics_sha256"):
        raise ContractError("canonical QC table is absent or hash-invalid")
    return {
        **manifest,
        "manifest_path": str(manifest_path),
        "manifest_sha256": sha256(manifest_path),
    }


def _copy_immutable(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if sha256(source) != sha256(target):
            raise ContractError(f"canonical working copy drifted: {target}")
        return
    partial = target.with_name(f".{target.name}.partial.{os.getpid()}")
    try:
        shutil.copy2(source, partial)
        os.replace(partial, target)
    finally:
        partial.unlink(missing_ok=True)


def canonicalize(
    config: CellCuratorConfig,
    *,
    config_path: Path,
    run_root: Path,
) -> dict[str, Any]:
    inspection = inspect(config)
    if inspection["blockers"]:
        raise ContractError(
            "input inspection has fatal representation blockers: "
            + "; ".join(inspection["blockers"])
        )
    source = Path(config.input.path)
    source_before = sha256(source)
    config_hash = sha256(config_path)
    manifest_path = run_root / "input_contract.manifest.json"
    if manifest_path.is_file():
        existing = json.loads(manifest_path.read_text())
        expected_fields = {
            "run_id": config.run.run_id,
            "input_sha256": source_before,
            "config_sha256": config_hash,
        }
        drifted = [
            field
            for field, expected in expected_fields.items()
            if str(existing.get(field, "")) != str(expected)
        ]
        if drifted:
            raise ContractError(
                "frozen run does not match the current input/configuration; "
                f"changed={drifted}. Use a new run_id so stale descendants cannot be reused"
            )
        canonical_path = Path(str(existing.get("canonical_path", "")))
        if not canonical_path.is_file() or sha256(canonical_path) != existing.get(
            "canonical_sha256"
        ):
            raise ContractError("frozen canonical working copy is absent or hash-invalid")
        prepared = existing.get("prepared", {})
        prepared_hashes = existing.get("prepared_sha256", {})
        for artifact in ("cells", "parent_membership", "stored_memberships"):
            artifact_path = Path(str(prepared.get(artifact, "")))
            if not artifact_path.is_file():
                raise ContractError(f"frozen prepared artifact is absent: {artifact}")
            if sha256(artifact_path) != str(prepared_hashes.get(artifact, "")):
                raise ContractError(f"frozen prepared artifact hash drifted: {artifact}")
        qc_manifest = Path(str(existing.get("qc", {}).get("manifest_path", "")))
        if not qc_manifest.is_file() or sha256(qc_manifest) != str(
            existing.get("qc", {}).get("manifest_sha256", "")
        ):
            raise ContractError("frozen QC manifest is absent or hash-invalid")
        compute_qc(config, run_root=run_root)
        return existing
    canonical = run_root / "canonical" / f"input.{config.input.kind}"
    _copy_immutable(source, canonical)
    qc = compute_qc(
        config,
        run_root=run_root,
        canonical_path=canonical,
        source_sha256=source_before,
    )
    canonical_input = config.input.model_copy(update={"path": str(canonical)})
    canonical_config = config.model_copy(update={"input": canonical_input})
    loaded = load_input(canonical_config, backed=True)
    try:
        obs = _obs_frame(loaded, config)
        keys = config.input.keys
        cells = pd.DataFrame(
            {
                "barcode": obs.index.astype(str),
                "scope": obs[keys.scope].astype(str).to_numpy(),
                "parent_id": obs[keys.canonical_parent].astype(str).to_numpy(),
            }
        )
        if keys.scope not in cells:
            cells[keys.scope] = obs[keys.scope].astype(str).to_numpy()
        if keys.canonical_parent not in cells:
            cells[keys.canonical_parent] = obs[keys.canonical_parent].astype(str).to_numpy()
        optional_keys = {
            "donor": keys.donor,
            "capture": keys.capture,
            "doublet_class": keys.doublet_class,
            "doublet_score": keys.doublet_score,
            "sex": keys.sex,
            "batch": keys.batch,
        }
        for _logical, key in optional_keys.items():
            if key and key in obs:
                cells[key] = obs[key].to_numpy()
        for key in config.input.guide_columns:
            cells[key] = obs[key].to_numpy()
        for key in config.input.technical_metrics.values():
            cells[key] = obs[key].to_numpy()
        for key in sorted(
            column for column in obs.columns if str(column).startswith("cell_curator_qc__")
        ):
            cells[str(key)] = obs[key].to_numpy()
        cells_path = run_root / "prepared" / "cells.tsv"
        cells_path.parent.mkdir(parents=True, exist_ok=True)
        payload = cells.to_csv(sep="\t", index=False, lineterminator="\n").encode()
        from .provenance import atomic_publish

        atomic_publish(cells_path, payload)
        parent_path = run_root / "prepared" / "parent_membership.tsv"
        parent = cells.copy()
        parent["closest_sibling"] = derive_closest_siblings(loaded, config, parent)
        atomic_publish(
            parent_path,
            parent.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
        )
        stored_rows: list[dict[str, str]] = []
        for resolution_key in keys.stored_resolutions:
            if resolution_key not in obs:
                raise ContractError(f"stored resolution key is absent: {resolution_key}")
            for barcode, scope, child in zip(
                obs.index.astype(str),
                obs[keys.scope].astype(str),
                obs[resolution_key].astype(str),
                strict=True,
            ):
                stored_rows.append(
                    {
                        "barcode": barcode,
                        "scope": scope,
                        "resolution": resolution_key,
                        "child_id": child,
                    }
                )
        stored = pd.DataFrame(
            stored_rows, columns=["barcode", "scope", "resolution", "child_id"]
        )
        stored_path = run_root / "prepared" / "stored_membership_long.tsv"
        atomic_publish(
            stored_path,
            stored.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
        )
    finally:
        close_input(loaded)
    source_after = sha256(source)
    if source_before != source_after:
        raise ContractError("source input changed during canonicalization")
    membership_hash = values_sha256(
        [
            f"{row.barcode}\t{row.scope}\t{row.parent_id}"
            for row in cells.itertuples(index=False)
        ]
    )
    manifest = {
        **inspection,
        "artifact_kind": "cell_curator_canonical_input",
        "run_id": config.run.run_id,
        "created_utc": datetime.now(UTC).isoformat(),
        "source_unchanged": True,
        "canonical_path": str(canonical),
        "canonical_sha256": sha256(canonical),
        "canonical_membership_sha256": membership_hash,
        "config_sha256": config_hash,
        "qc": {
            "manifest_path": qc["manifest_path"],
            "manifest_sha256": qc["manifest_sha256"],
            "cell_metrics_path": qc["cell_metrics_path"],
            "cell_metrics_sha256": qc["cell_metrics_sha256"],
        },
        "prepared": {
            "cells": str(cells_path),
            "parent_membership": str(parent_path),
            "stored_memberships": str(stored_path),
        },
        "prepared_sha256": {
            "cells": sha256(cells_path),
            "parent_membership": sha256(parent_path),
            "stored_memberships": sha256(stored_path),
        },
    }
    publish_json(manifest_path, manifest)
    return manifest


def derive_closest_siblings(
    loaded: LoadedInput,
    config: CellCuratorConfig,
    cells: pd.DataFrame,
) -> list[str]:
    keys = config.input.keys
    embedding_key = keys.embedding
    if embedding_key:
        candidate = next(iter(loaded.assays.values()))
        if embedding_key not in candidate.obsm:
            raise ContractError(f"declared embedding is absent: {embedding_key}")
        embedding = np.asarray(candidate.obsm[embedding_key])
    else:
        candidate = next(iter(loaded.assays.values()))
        expression = dense(candidate.X)
        embedding = expression[:, : min(20, expression.shape[1])]
    frame = cells[["scope", "parent_id"]].copy()
    frame["_row"] = np.arange(len(frame))
    result = np.empty(len(frame), dtype=object)
    for _scope, scoped in frame.groupby("scope", observed=True):
        centroids: dict[str, np.ndarray] = {}
        for parent, block in scoped.groupby("parent_id", observed=True):
            centroids[str(parent)] = embedding[block["_row"].to_numpy()].mean(axis=0)
        if len(centroids) < 2:
            only = next(iter(centroids))
            result[scoped["_row"].to_numpy()] = f"{only}__no_sibling"
            continue
        for parent, block in scoped.groupby("parent_id", observed=True):
            parent = str(parent)
            sibling = min(
                (candidate for candidate in centroids if candidate != parent),
                key=lambda candidate: (
                    float(np.linalg.norm(centroids[parent] - centroids[candidate])),
                    candidate,
                ),
            )
            result[block["_row"].to_numpy()] = sibling
    return list(map(str, result))

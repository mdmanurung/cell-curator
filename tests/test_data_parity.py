from __future__ import annotations

from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
import yaml

from cell_curator.config import CellCuratorConfig
from cell_curator.contracts import ContractError, sha256
from cell_curator.data import (
    canonicalize,
    compute_qc,
    detect_gene_id_convention,
    detect_symbol_column,
    gene_symbols,
    inspect,
    load_gene_mapping,
    matrix_for,
    normalization_fingerprint,
    targeted_panel_summary,
)

GENES = ["MT-CO1", "RPS3", "PCNA", "TYMS", "MKI67", "TOP2A", "CD3D"]


def _adata() -> ad.AnnData:
    counts = np.asarray(
        [
            [4, 2, 8, 7, 0, 0, 3],
            [2, 1, 7, 5, 0, 0, 4],
            [1, 2, 0, 0, 9, 8, 3],
            [3, 1, 0, 0, 7, 6, 2],
            [2, 1, 1, 0, 0, 0, 8],
            [1, 2, 0, 1, 0, 0, 7],
            [3, 2, 4, 3, 0, 0, 2],
            [2, 1, 0, 0, 5, 4, 2],
        ],
        dtype=np.float32,
    )
    obs = pd.DataFrame(
        {
            "canonical_cluster": ["0"] * 4 + ["1"] * 4,
            "lineage": ["immune"] * 8,
            "donor_id": ["d1", "d2"] * 4,
            "capture_id": ["c1", "c2"] * 4,
            "doublet_class": ["singlet"] * 8,
            "batch_category": pd.Categorical(
                ["b1", "b1", "b2", "b2"] * 2,
                categories=["b2", "b1"],
                ordered=True,
            ),
        },
        index=[f"cell-{index}" for index in range(8)],
    )
    var = pd.DataFrame({"gene_symbol": GENES}, index=GENES)
    result = ad.AnnData(X=np.log1p(counts), obs=obs, var=var)
    result.layers["counts"] = counts
    result.layers["logcounts"] = np.log1p(counts)
    result.obsm["X_pca_harmony"] = np.column_stack(
        [np.arange(8), np.arange(8) % 2, np.linspace(0, 1, 8)]
    ).astype(np.float32)
    result.obsm["spatial"] = np.column_stack([np.arange(8), np.arange(8)[::-1]]).astype(
        np.float32
    )
    result.uns["batch_category_colors"] = np.asarray(["#112233", "#445566"])
    result.raw = result
    return result


def _model_and_path(
    config: dict,
    tmp_path: Path,
    source: Path,
    *,
    kind: str,
    modality: str,
) -> tuple[CellCuratorConfig, Path]:
    config["input"].update({"path": str(source), "kind": kind})
    config["input"]["assays"]["transcriptome"]["modality"] = modality
    config["input"]["keys"]["spatial"] = "spatial"
    config["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    config.setdefault("guidance", {}).setdefault("context", {})["visualization_key"] = (
        "X_pca_harmony"
    )
    path = tmp_path / f"config-{kind}.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return CellCuratorConfig.model_validate(config), path


def test_gene_id_detection_mapping_and_raw_matrix(tmp_path: Path) -> None:
    identifiers = pd.Index(["ENSG000001.2", "ENSG000002.1"])
    query = ad.AnnData(
        X=np.ones((2, 2), dtype=np.float32),
        var=pd.DataFrame(index=identifiers),
    )
    query.raw = query
    mapping_path = tmp_path / "mapping.tsv"
    pd.DataFrame(
        {
            "ensembl_id": ["ENSG000001", "ENSG000002"],
            "gene_symbol": ["A", "B"],
        }
    ).to_csv(mapping_path, sep="\t", index=False)

    loaded = load_gene_mapping(mapping_path)
    assert loaded["ENSG000001"] == "A"
    assert detect_gene_id_convention(identifiers) == "ensembl"
    assert detect_gene_id_convention(["CD3D", "MS4A1"]) == "gene_symbol"
    assert detect_gene_id_convention(["chr1:1-100", "chr2:2-200"]) == "unknown"
    assert gene_symbols(query, mapping=mapping_path).tolist() == ["A", "B"]
    np.testing.assert_array_equal(matrix_for(query, {"source": "raw"}), query.raw.X)

    with_symbols = query.copy()
    with_symbols.var["gene_name"] = ["A", "B"]
    assert detect_symbol_column(with_symbols) == "gene_name"
    assert gene_symbols(with_symbols).tolist() == ["A", "B"]
    with pytest.raises(ContractError, match="incomplete"):
        gene_symbols(query, mapping={"ENSG000001": "A"})
    with pytest.raises(ContractError, match="blank source ID or symbol"):
        load_gene_mapping({"ENSG000001": pd.NA})


def test_configured_raw_representation_and_gene_mapping_are_executable(
    tmp_path: Path,
    config: dict,
) -> None:
    source = tmp_path / "ensembl-raw.h5ad"
    mapping_path = tmp_path / "ensembl-to-symbol.tsv"
    identifiers = [f"ENSG{index:011d}.1" for index in range(1, len(GENES) + 1)]
    pd.DataFrame(
        {"ensembl_id": [item.removesuffix(".1") for item in identifiers], "symbol": GENES}
    ).to_csv(mapping_path, sep="\t", index=False)
    adata = _adata()
    adata.var_names = identifiers
    adata.var = pd.DataFrame(index=identifiers)
    adata.raw = adata
    adata.write_h5ad(source)

    assay = config["input"]["assays"]["transcriptome"]
    assay["feature_id_column"] = ""
    assay["gene_mapping_path"] = str(mapping_path)
    representation = assay["representations"]["logcounts"]
    representation.update({"source": "raw", "key": ""})
    representation["detection"] = {
        "source": "raw",
        "key": "",
        "required": True,
        "semantics": "nonnegative raw matrix",
    }
    model, config_path = _model_and_path(config, tmp_path, source, kind="h5ad", modality="root")
    summary = inspect(model)
    inventory = summary["assays"]["transcriptome"]["feature_ids"]
    assert inventory["symbol_status"] == "mapped"
    assert inventory["symbol_mapping_available"] is True
    assert not summary["blockers"]

    manifest = canonicalize(model, config_path=config_path, run_root=tmp_path / "mapped-run")
    canonical = ad.read_h5ad(manifest["canonical_path"])
    assert canonical.var["cell_curator_gene_symbol"].astype(str).tolist() == GENES


def test_normalization_fingerprints_are_representation_aware() -> None:
    counts = np.asarray([[0, 2, 10], [1, 4, 3]], dtype=float)
    assert (
        normalization_fingerprint(counts, assay_kind="rna", declared_name="counts")[
            "classification"
        ]
        == "counts"
    )
    log_fingerprint = normalization_fingerprint(
        np.log1p(counts), assay_kind="rna", declared_name="logcounts"
    )
    assert log_fingerprint["classification"] == "log_normalized"
    assert log_fingerprint["verdict"] == "valid"

    signed = normalization_fingerprint(
        np.asarray([[-1.0, 1.0], [0.0, 2.0]]),
        assay_kind="protein",
        declared_name="corrected",
    )
    assert signed["classification"] == "corrected_signed"
    assert signed["detection_safe"] is False
    assert signed["verdict"] == "valid"
    assert (
        normalization_fingerprint(
            np.asarray([[-1.0, 1.0]]),
            assay_kind="protein",
            declared_name="detection",
            role="detection",
        )["verdict"]
        == "invalid"
    )
    assert (
        normalization_fingerprint(np.asarray([[0.2, 150.3]]), assay_kind="rna")["verdict"]
        == "ambiguous"
    )
    assert (
        normalization_fingerprint(
            np.asarray([[np.nan, 1.0]]), assay_kind="rna", declared_name="logcounts"
        )["verdict"]
        == "invalid"
    )
    broad = targeted_panel_summary([f"G{index}" for index in range(1_500)], assay_kind="rna")
    assert broad["is_targeted_panel"] is False
    declared = targeted_panel_summary(
        [f"G{index}" for index in range(1_500)],
        assay_kind="rna",
        declared_targeted=True,
    )
    assert declared["is_targeted_panel"] is True
    assert declared["detection_basis"] == "declared_technology"


def test_invalid_normalization_is_a_canonicalization_hard_stop(
    tmp_path: Path,
    config: dict,
) -> None:
    source = tmp_path / "ambiguous.h5ad"
    adata = _adata()
    adata.layers["logcounts"] = np.full(adata.shape, 250.5, dtype=np.float32)
    adata.write_h5ad(source)
    model, config_path = _model_and_path(config, tmp_path, source, kind="h5ad", modality="root")
    summary = inspect(model)
    assert any("implausible scale" in blocker for blocker in summary["blockers"])
    with pytest.raises(ContractError, match="fatal representation blockers"):
        canonicalize(model, config_path=config_path, run_root=tmp_path / "blocked-run")


def test_h5ad_inspection_canonicalization_and_qc_are_source_preserving(
    tmp_path: Path,
    config: dict,
) -> None:
    source = tmp_path / "source.h5ad"
    _adata().write_h5ad(source)
    source_hash = sha256(source)
    model, config_path = _model_and_path(config, tmp_path, source, kind="h5ad", modality="root")

    summary = inspect(model)
    assay = summary["assays"]["transcriptome"]
    assert not summary["blockers"]
    assert summary["representations"]["visualization"] == "X_pca_harmony"
    assert summary["representations"]["spatial"] == "spatial"
    assert assay["feature_ids"]["id_convention"] == "gene_symbol"
    assert assay["panel"]["is_targeted_panel"] is True
    assert assay["panel"]["requires_program_coverage_gating"] is True
    assert (
        assay["representations"]["logcounts"]["normalization"]["classification"]
        == "log_normalized"
    )
    assert (
        assay["representations"]["logcounts"]["detection"]["classification"] == "detection_safe"
    )
    covariate = summary["categorical_covariates"]["batch_category"]
    assert covariate["ordered"] is True
    assert covariate["categories"] == ["b2", "b1"]
    assert covariate["colors"] == ["#112233", "#445566"]

    run_root = tmp_path / "run-h5ad"
    manifest = canonicalize(model, config_path=config_path, run_root=run_root)
    assert sha256(source) == source_hash
    assert manifest["source_unchanged"] is True
    assert manifest["qc"]["cell_metrics_sha256"] == sha256(
        Path(manifest["qc"]["cell_metrics_path"])
    )
    canonical = ad.read_h5ad(manifest["canonical_path"])
    assert {
        "cell_curator_original_feature_id",
        "cell_curator_feature_id",
        "cell_curator_gene_symbol",
    }.issubset(canonical.var.columns)
    expected_metrics = {
        "cell_curator_qc_total_counts",
        "cell_curator_qc_n_features_by_counts",
        "cell_curator_qc_pct_counts_mt",
        "cell_curator_qc_pct_counts_ribo",
        "cell_curator_qc_s_score",
        "cell_curator_qc_g2m_score",
        "cell_curator_qc_phase",
    }
    assert expected_metrics.issubset(canonical.obs.columns)
    assert list(canonical.obs["cell_curator_qc_phase"].cat.categories) == [
        "G1",
        "S",
        "G2M",
    ]
    durable = pd.read_csv(manifest["qc"]["cell_metrics_path"], sep="\t")
    assert "transcriptome__pct_counts_mt" in durable
    assert "transcriptome__phase" in durable

    verified = compute_qc(model, run_root=run_root)
    assert verified["canonical_sha256"] == manifest["canonical_sha256"]
    resumed = canonicalize(model, config_path=config_path, run_root=run_root)
    assert resumed["canonical_sha256"] == manifest["canonical_sha256"]
    assert sha256(source) == source_hash


def test_mudata_canonical_qc_preserves_modalities_and_source(
    tmp_path: Path,
    config: dict,
) -> None:
    rna = _adata()
    source = tmp_path / "source.h5mu"
    mdata = md.MuData({"rna": rna})
    for column in rna.obs:
        mdata.obs[column] = rna.obs[column].to_numpy()
    mdata.uns["batch_category_colors"] = rna.uns["batch_category_colors"]
    mdata.write_h5mu(source)
    source_hash = sha256(source)
    model, config_path = _model_and_path(config, tmp_path, source, kind="h5mu", modality="rna")

    run_root = tmp_path / "run-h5mu"
    summary = inspect(model)
    assert summary["input_kind"] == "h5mu"
    manifest = canonicalize(model, config_path=config_path, run_root=run_root)
    assert sha256(source) == source_hash
    canonical = md.read_h5mu(manifest["canonical_path"])
    assert "rna" in canonical.mod
    assert "cell_curator_gene_symbol" in canonical.mod["rna"].var
    assert "cell_curator_qc_total_counts" in canonical.mod["rna"].obs
    assert "cell_curator_qc__transcriptome__total_counts" in canonical.obs
    assert compute_qc(model, run_root=run_root)["source_mutated"] is False
    assert sha256(source) == source_hash

"""Input shapes that real objects have and the synthetic fixtures did not."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse
from test_end_to_end import _config

from cell_curator.config import load_config
from cell_curator.data import inspect


@pytest.mark.parametrize("source", ["X", "raw"])
def test_a_sparse_backed_x_or_raw_can_be_inspected(tmp_path: Path, source: str) -> None:
    """Backed X and raw are h5py sparse datasets, which np.asarray cannot read."""

    rng = np.random.default_rng(7)
    counts = np.zeros((60, 6), dtype=np.float32)
    counts[:30, :2] = rng.poisson(8, size=(30, 2)) + 1
    counts[30:, 2:4] = rng.poisson(8, size=(30, 2)) + 1
    genes = ["A1", "A2", "B1", "B2", "N1", "N2"]
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.log1p(counts)),
        obs=pd.DataFrame(
            {"canonical_cluster": ["0"] * 30 + ["1"] * 30, "lineage": ["immune"] * 60},
            index=[f"c{i}" for i in range(60)],
        ),
        var=pd.DataFrame({"gene_symbol": genes}, index=genes),
    )
    adata.layers["counts"] = sparse.csr_matrix(counts)
    adata.raw = adata.copy()
    adata.obsm["X_pca_harmony"] = rng.normal(size=(60, 3)).astype(np.float32)
    source_path = tmp_path / "input.h5ad"
    adata.write_h5ad(source_path)
    markers = tmp_path / "markers.json"
    markers.write_text(json.dumps({"Type A": {"pos": ["A1"]}}))
    config_path = _config(tmp_path, source_path, markers)
    value = yaml.safe_load(config_path.read_text())
    representation = value["input"]["assays"]["transcriptome"]["representations"]["logcounts"]
    representation["source"] = source
    representation.pop("key", None)
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))

    report = inspect(load_config(config_path))
    assert report["canonical_parent_count"] == 2


def test_numeric_looking_barcodes_survive_to_the_written_object(tmp_path: Path) -> None:
    """Barcodes "0".."59" read back from a TSV as integers unless told otherwise."""

    from test_end_to_end import _confirm, _record_critics

    from cell_curator.pipeline import run_annotation, run_root
    from cell_curator.writeback import write_back

    rng = np.random.default_rng(7)
    counts = np.zeros((60, 6), dtype=np.float32)
    counts[:30, :2] = rng.poisson(8, size=(30, 2)) + 1
    counts[30:, 2:4] = rng.poisson(8, size=(30, 2)) + 1
    genes = ["A1", "A2", "B1", "B2", "N1", "N2"]
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=pd.DataFrame(
            {"canonical_cluster": ["0"] * 30 + ["1"] * 30, "lineage": ["immune"] * 60},
            index=[str(i) for i in range(60)],
        ),
        var=pd.DataFrame({"gene_symbol": genes}, index=genes),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = np.vstack(
        [rng.normal(-3, 0.2, size=(30, 3)), rng.normal(3, 0.2, size=(30, 3))]
    ).astype(np.float32)
    source = tmp_path / "input.h5ad"
    adata.write_h5ad(source)
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            {
                "Type A": {"pos": ["A1", "A2"], "neg": ["B1", "B2"]},
                "Type B": {"pos": ["B1", "B2"], "neg": ["A1", "A2"]},
            }
        )
    )
    config_path = _config(tmp_path, source, markers)
    run_annotation(config_path)
    config = load_config(config_path)
    root = run_root(config)
    _record_critics(config, root, tmp_path)
    output = tmp_path / "out.h5ad"
    write_back(config, root, output_path=output, confirm=_confirm)
    written = ad.read_h5ad(output).obs["cell_curator_type"].astype(str)
    assert list(written.iloc[[0, 59]]) == ["Type A", "Type B"]


def test_10x_feature_names_with_a_few_ensembl_ids_are_gene_symbols() -> None:
    from cell_curator.data import detect_gene_id_convention

    tenx = [f"GENE{i}" for i in range(980)] + [f"ENSG{i:011d}" for i in range(20)]
    assert detect_gene_id_convention(tenx) == "gene_symbol"
    assert detect_gene_id_convention([f"ENSG{i:011d}.1" for i in range(100)]) == "ensembl"
    half = [f"GENE{i}" for i in range(50)] + [f"ENSG{i:011d}" for i in range(50)]
    assert detect_gene_id_convention(half) != "gene_symbol"


def test_a_blank_scope_or_parent_key_is_refused_at_load(tmp_path: Path) -> None:
    from test_end_to_end import _synthetic_input

    source = tmp_path / "input.h5ad"
    _synthetic_input(source)
    markers = tmp_path / "markers.json"
    markers.write_text(json.dumps({"Type A": {"pos": ["A1"]}}))
    config_path = _config(tmp_path, source, markers)
    for key in ("scope", "canonical_parent"):
        value = yaml.safe_load(config_path.read_text())
        value["input"]["keys"][key] = " "
        config_path.write_text(yaml.safe_dump(value, sort_keys=False))
        with pytest.raises(ValueError, match="name obs columns"):
            load_config(config_path)

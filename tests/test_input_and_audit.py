from __future__ import annotations

from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest

from cell_curator import audit_parents
from cell_curator.config import CellCuratorConfig
from cell_curator.data import inspect


def _obs() -> pd.DataFrame:
    barcodes = [f"cell-{index}" for index in range(6)]
    return pd.DataFrame(
        {
            "canonical_cluster": ["0", "0", "0", "1", "1", "1"],
            "lineage": ["immune"] * 6,
            "donor_id": ["d1", "d2", "d1", "d1", "d2", "d2"],
            "capture_id": ["c1", "c2", "c1", "c1", "c2", "c2"],
            "doublet_class": ["singlet"] * 6,
            "stored_1": ["0", "0", "2", "1", "1", "1"],
        },
        index=barcodes,
    )


def _rna(path: Path) -> None:
    obs = _obs()
    adata = ad.AnnData(
        X=np.ones((6, 2), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame({"gene_symbol": ["G1", "G2"]}, index=["G1", "G2"]),
    )
    adata.layers["counts"] = np.ones((6, 2), dtype=np.float32)
    adata.layers["logcounts"] = np.log2(adata.layers["counts"] + 1)
    adata.obsm["X_pca_harmony"] = np.ones((6, 2), dtype=np.float32)
    adata.write(path)


def _multimodal(path: Path) -> None:
    obs = _obs()
    rna = ad.AnnData(
        X=np.ones((6, 2), dtype=np.float32),
        obs=pd.DataFrame(index=obs.index),
        var=pd.DataFrame({"gene_symbol": ["G1", "G2"]}, index=["G1", "G2"]),
    )
    rna.layers["counts"] = np.ones((6, 2), dtype=np.float32)
    rna.layers["logcounts"] = np.ones((6, 2), dtype=np.float32)
    rna.obsm["X_pca_harmony"] = np.ones((6, 2), dtype=np.float32)
    protein = ad.AnnData(
        X=np.ones((6, 2), dtype=np.float32),
        obs=pd.DataFrame(index=obs.index),
        var=pd.DataFrame(index=["P1", "P2"]),
    )
    protein.layers["corrected"] = np.array(
        [[-1.0, 1.0], [0.0, 1.0], [1.0, 2.0]] * 2, dtype=np.float32
    )
    protein.layers["logcounts"] = np.array(
        [[0.0, 1.0], [0.0, 1.0], [1.0, 2.0]] * 2, dtype=np.float32
    )
    mdata = md.MuData({"rna": rna, "protein": protein})
    for column in obs:
        mdata.obs[column] = obs[column].to_numpy()
    mdata.write(path)


def _set_two_parents(config: dict) -> None:
    config["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }


def test_rna_only_h5ad_contract(tmp_path: Path, config: dict) -> None:
    source = tmp_path / "rna.h5ad"
    _rna(source)
    config["input"]["path"] = str(source)
    config["input"]["keys"]["stored_resolutions"] = ["stored_1"]
    _set_two_parents(config)
    manifest = inspect(CellCuratorConfig.model_validate(config))
    assert set(manifest["assays"]) == {"transcriptome"}
    assert manifest["assays"]["transcriptome"]["kind"] == "rna"
    assert not manifest["blockers"]


def test_atac_only_h5ad_is_auditable_but_peak_only_labels_are_blocked(
    tmp_path: Path, config: dict
) -> None:
    source = tmp_path / "atac.h5ad"
    obs = _obs()
    adata = ad.AnnData(
        X=np.ones((6, 2), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["chr1:1-100", "chr1:200-300"]),
    )
    adata.layers["counts"] = np.ones((6, 2), dtype=np.float32)
    adata.write(source)
    config["input"].update({"path": str(source), "kind": "h5ad"})
    config["input"]["assays"] = {
        "chromatin": {
            "modality": "root",
            "kind": "atac",
            "feature_space": "peak",
            "capabilities": ["state"],
            "representations": {
                "accessibility": {
                    "source": "X",
                    "ranking_backend": "project_atac_da",
                    "ontology_compatible": False,
                    "detection": {
                        "source": "layer",
                        "key": "counts",
                        "required": False,
                        "semantics": "fragments > 0",
                    },
                }
            },
        }
    }
    config["evidence"]["lanes"] = [{"assay": "chromatin", "representation": "accessibility"}]
    config["input"]["keys"]["embedding"] = ""
    config["guidance"]["context"]["visualization_key"] = ""
    _set_two_parents(config)
    summary = inspect(CellCuratorConfig.model_validate(config))
    assert summary["assays"]["chromatin"]["kind"] == "atac"
    assert not any("ontology-compatible" in item for item in summary["blockers"])


def _configure_multimodal(config: dict, path: Path) -> None:
    config["input"].update({"path": str(path), "kind": "h5mu"})
    config["input"]["assays"]["transcriptome"]["modality"] = "rna"
    config["input"]["assays"]["protein"] = {
        "modality": "protein",
        "kind": "protein",
        "feature_space": "antibody",
        "capabilities": ["identity", "subtype", "state"],
        "representations": {
            "corrected": {
                "source": "layer",
                "key": "corrected",
                "ranking_backend": "rapids_singlecell",
                "ontology_compatible": False,
                "detection": {
                    "source": "layer",
                    "key": "logcounts",
                    "required": True,
                    "semantics": "protein logcounts > 0",
                },
            }
        },
    }
    config["evidence"]["lanes"].append({"assay": "protein", "representation": "corrected"})
    _set_two_parents(config)


def test_multimodal_mudata_and_signed_effect_representation(
    tmp_path: Path, config: dict
) -> None:
    source = tmp_path / "fixture.h5mu"
    _multimodal(source)
    _configure_multimodal(config, source)
    summary = inspect(CellCuratorConfig.model_validate(config))
    assert set(summary["assays"]) == {"transcriptome", "protein"}
    protein = summary["assays"]["protein"]["representations"]["corrected"]
    assert protein["detection"]["nonnegative"] is True
    assert not summary["blockers"]


def test_signed_required_detection_blocks_acceptance(tmp_path: Path, config: dict) -> None:
    source = tmp_path / "fixture.h5mu"
    _multimodal(source)
    mdata = md.read_h5mu(source)
    mdata.mod["protein"].layers["logcounts"] = np.full((6, 2), -1.0)
    broken = tmp_path / "broken.h5mu"
    mdata.write(broken)
    _configure_multimodal(config, broken)
    summary = inspect(CellCuratorConfig.model_validate(config))
    assert any("requires a nonnegative detection" in item for item in summary["blockers"])


def test_input_contract_fails_on_annotation_without_local_grounding(config: dict) -> None:
    config["run"]["mode"] = "annotation"
    with pytest.raises(ValueError, match="local marker"):
        CellCuratorConfig.model_validate(config)


def test_parent_audit_uses_configured_metadata_and_technical_metrics(
    config: dict,
) -> None:
    config["input"]["guide_columns"] = ["guide_a", "guide_b"]
    config["input"]["technical_metrics"] = {"complexity": "n_features"}
    cells = pd.DataFrame(
        {
            "barcode": ["a", "b", "c", "d"],
            "canonical_cluster": ["0", "0", "1", "1"],
            "lineage": ["immune"] * 4,
            "donor_id": ["d1", "d1", "d1", "d2"],
            "capture_id": ["c1", "c1", "c1", "c2"],
            "doublet_class": ["doublet", "doublet", "singlet", "singlet"],
            "n_features": [50, 50, 100, 100],
            "guide_a": ["T", "B", "B", "B"],
            "guide_b": ["B", "B", "B", "B"],
        }
    )
    result = audit_parents.build_audit(config, cells).set_index("parent_id")
    assert bool(result.loc["0", "doublet_enriched"])
    assert bool(result.loc["0", "donor_dominated"])
    assert bool(result.loc["0", "guide_top_disagreement"])
    assert "median_technical__complexity" in result.columns


def test_parent_audit_allows_missing_replication_axes(config: dict) -> None:
    config["input"]["keys"].update({"donor": "", "capture": "", "doublet_class": ""})
    config["adaptive"]["audit"].update(
        {
            "donor_max_fraction": None,
            "capture_max_fraction": None,
            "doublet_fraction_max": None,
        }
    )
    cells = pd.DataFrame(
        {
            "barcode": ["a", "b"],
            "canonical_cluster": ["0", "0"],
            "lineage": ["immune", "immune"],
        }
    )
    row = audit_parents.build_audit(config, cells).iloc[0]
    assert row["n_donors"] == 0
    assert not bool(row["donor_dominated"])

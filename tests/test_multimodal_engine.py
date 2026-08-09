from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import yaml

from cell_curator.capabilities import build_capability_registry
from cell_curator.config import load_config
from cell_curator.evidence import consolidate_lanes, run_lane
from cell_curator.knowledge import collect_programs, providers_from_config
from cell_curator.models import BackendReceipt, CapabilityLevel


def _modality(
    values: np.ndarray,
    features: list[str],
    barcodes: list[str],
    *,
    layer: str,
) -> ad.AnnData:
    result = ad.AnnData(
        X=values.astype(np.float32),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=features),
    )
    result.layers[layer] = values.astype(np.float32)
    result.layers["counts"] = (values > 0).astype(np.float32)
    return result


def test_executable_multimodal_routes_emit_receipts_and_visible_disagreement(
    tmp_path: Path, config: dict
) -> None:
    barcodes = [f"cell-{index:02d}" for index in range(20)]
    groups = np.array(["0"] * 10 + ["1"] * 10)
    direct = np.column_stack(
        [
            np.r_[np.full(10, 5.0), np.zeros(10)],
            np.r_[np.zeros(10), np.full(10, 5.0)],
            np.r_[np.full(5, 3.0), np.zeros(15)],
        ]
    )
    reversed_panel = direct[:, [1, 0]]
    rna = _modality(direct, ["G1", "G2", "ACT"], barcodes, layer="logcounts")
    rna.var["gene_symbol"] = rna.var_names.astype(str)
    rna.obsm["X_pca_harmony"] = np.column_stack([direct[:, 0], direct[:, 1]])
    protein = _modality(reversed_panel, ["CD3", "CD19"], barcodes, layer="arcsinh")
    atac = _modality(direct[:, :2], ["peak1", "peak2"], barcodes, layer="accessibility")
    activity = _modality(direct[:, :2], ["G1", "G2"], barcodes, layer="activity")
    spatial = _modality(direct[:, :2], ["G1", "G2"], barcodes, layer="logcounts")
    spatial.obsm["spatial"] = np.column_stack([np.arange(20), groups.astype(int) * 100])
    mdata = md.MuData(
        {
            "rna": rna,
            "protein": protein,
            "atac": atac,
            "activity": activity,
            "spatial": spatial,
        }
    )
    mdata.obs["canonical_cluster"] = groups
    mdata.obs["lineage"] = "immune"
    source = tmp_path / "multiome.h5mu"
    mdata.write_h5mu(source)

    mapping = tmp_path / "peak_to_gene.tsv"
    mapping.write_text("peak\tgene\tweight\npeak1\tG1\t1\npeak2\tG2\t1\n")
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            [
                {"label": "Type A", "positive": ["G1"], "negative": ["G2"]},
                {"label": "Type B", "positive": ["G2"], "negative": ["G1"]},
                {
                    "label": "Type A",
                    "feature_space": "antibody",
                    "positive": ["CD3", "CD4"],
                    "negative": ["CD19"],
                },
                {
                    "label": "Type B",
                    "feature_space": "antibody",
                    "positive": ["CD19"],
                    "negative": ["CD3"],
                },
                {
                    "label": "Activated",
                    "axis": "cell_state",
                    "positive": ["ACT"],
                    "feature_space": "gene",
                },
            ]
        )
    )
    config["run"].update(
        {"run_id": "multimodal-v1", "output_root": str(tmp_path / "runs"), "mode": "annotation"}
    )
    config["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": str(markers), "version": "fixture-v1"}
    ]
    config["input"].update({"kind": "h5mu", "path": str(source)})
    config["input"]["keys"].update(
        {"embedding": "X_pca_harmony", "spatial": "spatial", "stored_resolutions": []}
    )
    config["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    config["input"]["assays"] = {
        "rna": {
            "modality": "rna",
            "kind": "rna",
            "feature_space": "gene",
            "feature_id_column": "gene_symbol",
            "representations": {
                "logcounts": {
                    "source": "layer",
                    "key": "logcounts",
                    "ranking_backend": "scipy_mannwhitney",
                    "ontology_compatible": True,
                    "detection": {"source": "layer", "key": "counts", "required": True},
                }
            },
        },
        "protein": {
            "modality": "protein",
            "kind": "protein",
            "feature_space": "antibody",
            "representations": {
                "arcsinh": {
                    "source": "layer",
                    "key": "arcsinh",
                    "ranking_backend": "scipy_mannwhitney",
                    "ontology_compatible": True,
                    "detection": {"source": "layer", "key": "counts"},
                }
            },
        },
        "atac": {
            "modality": "atac",
            "kind": "atac",
            "feature_space": "peak",
            "representations": {
                "accessibility": {
                    "source": "layer",
                    "key": "accessibility",
                    "ranking_backend": "scipy_mannwhitney",
                    "peak_to_gene_path": str(mapping),
                    "detection": {"source": "layer", "key": "counts"},
                }
            },
        },
        "activity": {
            "modality": "activity",
            "kind": "gene_activity",
            "feature_space": "gene_activity",
            "representations": {
                "activity": {
                    "source": "layer",
                    "key": "activity",
                    "feature_space": "gene_activity",
                    "ranking_backend": "scipy_mannwhitney",
                    "detection": {"source": "layer", "key": "counts"},
                }
            },
        },
        "spatial": {
            "modality": "spatial",
            "kind": "spatial",
            "feature_space": "gene",
            "representations": {
                "logcounts": {
                    "source": "layer",
                    "key": "logcounts",
                    "feature_space": "gene",
                    "ranking_backend": "scipy_mannwhitney",
                    "detection": {"source": "layer", "key": "counts"},
                }
            },
        },
    }
    config["evidence"]["lanes"] = [
        {"assay": "rna", "representation": "logcounts", "required_for_acceptance": True},
        {"assay": "protein", "representation": "arcsinh", "required_for_acceptance": False},
        {"assay": "atac", "representation": "accessibility", "required_for_acceptance": False},
        {"assay": "activity", "representation": "activity", "required_for_acceptance": False},
        {"assay": "spatial", "representation": "logcounts", "required_for_acceptance": False},
    ]
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    model = load_config(config_path)
    programs = collect_programs(providers_from_config(model.knowledge.providers))
    root = Path(model.run.output_root) / model.run.run_id

    results = [
        run_lane(
            config=model,
            config_path=config_path,
            run_root=root,
            assay_id=lane.assay,
            representation_name=lane.representation,
            programs=programs,
        )
        for lane in model.evidence.lanes
    ]
    for result in results:
        receipt = BackendReceipt.model_validate_json(result.receipt_path.read_text())
        assert receipt.status == "complete"
        assert receipt.fallback_used is False

    capabilities = {(item.assay, item.level) for item in build_capability_registry(model)}
    assert ("rna", CapabilityLevel.ANNOTATION) in capabilities
    assert ("protein", CapabilityLevel.ANNOTATION) in capabilities
    assert ("atac", CapabilityLevel.CORROBORATION) in capabilities
    assert ("activity", CapabilityLevel.CORROBORATION) in capabilities
    assert ("spatial", CapabilityLevel.CORROBORATION) in capabilities

    protein_evidence = next(
        item.evidence for item in results if item.lane_key == "protein__arcsinh"
    )
    assert "CD4" in set(protein_evidence["missing_positive"])
    spatial_evidence = next(
        item.evidence for item in results if item.lane_key == "spatial__logcounts"
    )
    assert "spatial_neighbor_coherence" in spatial_evidence
    atac_evidence = next(
        item.evidence for item in results if item.lane_key == "atac__accessibility"
    )
    assert set(atac_evidence["feature_space"]) == {"gene"}

    scores, proposals = consolidate_lanes(config=model, run_root=root)
    assert len(proposals) == 2
    assert scores["multimodal_disagreement"].astype(bool).any()
    assert set(proposals["cell_state"]) == {"Activated", "Unknown"}

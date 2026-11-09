"""Synthetic-data checks for the v8 researcher packet and agreement scripts."""

import importlib.util
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


packets = load(ROOT / "benchmarks/scripts/researcher_packets_v8.py")
agree = load(ROOT / "benchmarks/v8-researcher/map_and_agree.py")


def synthetic_adata():
    rng = np.random.default_rng(0)
    genes = ["CD3D", "CD3E", "MS4A1", "CD79A", "LYZ", "GENE1", "GENE2", "GENE3"]
    cluster = np.repeat(["0", "1", "2"], [60, 50, 40])
    x = rng.poisson(0.2, size=(150, len(genes))).astype(np.float32)
    x[cluster == "0", :2] += rng.poisson(4, size=(60, 2))
    x[cluster == "1", 2:4] += rng.poisson(4, size=(50, 2))
    x[cluster == "2", 4] += rng.poisson(4, size=40)
    obs = pd.DataFrame(
        {
            "donor_id": rng.choice(["d1", "d2", "d3"], 150),
            "tissue": np.where(cluster == "2", "spleen", "blood"),
            "assay": "10x 3' v3",
            "bench_cluster": cluster,
            "truth_label": "must not appear",
        },
        index=[f"c{i}" for i in range(150)],
    )
    return ad.AnnData(np.log1p(x), obs=obs, var=pd.DataFrame(index=genes))


def test_packet_finds_cluster_markers_and_hides_extra_columns(tmp_path, monkeypatch):
    h5ad = tmp_path / "input.h5ad"
    synthetic_adata().write_h5ad(h5ad)
    monkeypatch.setattr(sys, "argv", ["x", "toy", str(h5ad), str(tmp_path / "out")])
    monkeypatch.setattr(packets, "N_TOP_GENES", 3)
    packets.main()
    text = (tmp_path / "out/packet.json").read_text()
    assert "must not appear" not in text
    packet = json.loads(text)
    by_cluster = {c["cluster"]: c for c in packet["clusters"]}
    assert {g["gene"] for g in by_cluster["0"]["top_genes"][:2]} == {"CD3D", "CD3E"}
    assert by_cluster["2"]["top_genes"][0]["gene"] == "LYZ"
    assert by_cluster["2"]["tissue_mix"] == {"spleen": 1.0}
    assert by_cluster["0"]["n_cells"] == 60
    assert "CD3D" in by_cluster["0"]["lineage_panel"]
    assert "CD8A" in packet["panel_genes_not_measured"]
    with pytest.raises(FileExistsError):
        packets.main()


def test_fleiss_kappa_known_values():
    perfect = np.array([["a"] * 3, ["b"] * 3, ["a"] * 3, ["b"] * 3])
    assert agree.fleiss_kappa(perfect) == pytest.approx(1.0)
    # Wikipedia-style check: two raters split evenly on every item gives kappa below 0.
    split = np.array([["a", "a", "b"], ["a", "b", "b"], ["a", "a", "b"], ["a", "b", "b"]])
    assert agree.fleiss_kappa(split) < 0.5
    assert np.isnan(agree.fleiss_kappa(np.array([["a"] * 3, ["a"] * 3])))


def test_mapping_and_agreement(tmp_path, monkeypatch):
    obo = tmp_path / "cl.obo"
    obo.write_text(
        "format-version: 1.2\n\n"
        "[Term]\nid: CL:0000000\nname: cell\n\n"
        "[Term]\nid: CL:0000084\nname: T cell\nis_a: CL:0000000\n\n"
        "[Term]\nid: CL:0000624\nname: CD4-positive, alpha-beta T cell\n"
        'synonym: "CD4 T cell" EXACT []\nis_a: CL:0000084\n\n'
        "[Term]\nid: CL:0000625\nname: CD8-positive, alpha-beta T cell\nis_a: CL:0000084\n\n"
        "[Term]\nid: CL:0000236\nname: B cell\nis_a: CL:0000000\n"
        'synonym: "B lymphocyte" EXACT []\n\n'
        "[Term]\nid: CL:0000623\nname: natural killer cell\nis_a: CL:0000000\n"
    )
    monkeypatch.setattr(
        agree,
        "COARSE_LINEAGE",
        {"CL:0000084": "T cell", "CL:0000236": "B cell", "CL:0000623": "natural killer cell"},
    )
    ontology = agree.CellOntology.load(obo)
    assert agree.map_label(ontology, "cd4 t cell.") == ("CL:0000624", "exact")
    assert agree.map_label(ontology, "B lymphocytes")[1] == "unmatched"
    assert agree.map_label(ontology, "B cells") == ("CL:0000236", "exact")
    assert agree.map_label(ontology, "Unresolved") == (agree.UNRESOLVED, "no_call")
    assert agree.map_label(ontology, "T cell (naive)") == (agree.UNRESOLVED, "unmatched")
    assert agree.coarse_id(ontology, "CL:0000625") == "CL:0000084"

    # Three clusters: 0 all CD4; 1 CD4, CD8, T cell; 2 B cell, B cell, unresolved.
    labels = {
        1: ["CD4 T cell", "CD4 T cell", "B cell"],
        2: ["CD4 T cell", "CD8-positive, alpha-beta T cell", "B cell"],
        3: ["CD4 T cell", "T cell", "unresolved"],
    }
    for k, rows in labels.items():
        folder = tmp_path / "resp/toy"
        folder.mkdir(parents=True, exist_ok=True)
        table = pd.DataFrame(
            {"cluster": ["0", "1", "2"], "label": rows, "cl_id": "", "confidence": "high"}
        )
        table.to_csv(folder / f"researcher_{k}.tsv", sep="\t", index=False)
    monkeypatch.setattr(
        sys, "argv", ["x", str(obo), str(tmp_path / "resp"), str(tmp_path / "o")]
    )
    agree.main()
    result = pd.read_csv(tmp_path / "o/agreement.tsv", sep="\t").set_index(["level", "measure"])
    # exact: cluster 0 agrees in 3 of 3 pairs; cluster 1 in 0 (CD4, CD8, T cell); cluster 2 in 1.
    assert result.loc[("exact", "pairwise_agreement"), "value"] == pytest.approx(4 / 9)
    # coarse: cluster 1 is T cell for all three, so it agrees in 3 of 3 pairs.
    assert result.loc[("coarse_lineage", "pairwise_agreement"), "value"] == pytest.approx(7 / 9)
    # ancestor_or_equal: cluster 1 has CD4/T cell and CD8/T cell related, CD4/CD8 not: 2 of 3.
    assert result.loc[("ancestor_or_equal", "pairwise_agreement"), "value"] == pytest.approx(
        (3 + 2 + 1) / 9
    )
    assert result.loc[("exact", "fraction_unresolved"), "value"] == pytest.approx(1 / 9)

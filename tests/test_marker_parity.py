from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from cell_curator.config import CellCuratorConfig
from cell_curator.contracts import sha256
from cell_curator.knowledge import ReferenceProvider
from cell_curator.markers import (
    assemble_marker_programs,
    marker_dotplot_svg,
    profile_marker_panel,
    run_marker_analysis,
    run_marker_analysis_from_config,
    score_aucell,
    score_ucell,
    write_marker_artifacts,
)
from cell_curator.models import MarkerProgram


def _marker_fixture() -> tuple[np.ndarray, pd.Series, pd.Index, list[MarkerProgram]]:
    matrix = np.asarray(
        [
            [8.0, 7.0, 0.0, 0.0],
            [7.0, 8.0, 0.0, 0.0],
            [9.0, 6.0, 0.0, 0.0],
            [6.0, 9.0, 0.0, 0.0],
            [0.0, 0.0, 8.0, 7.0],
            [0.0, 0.0, 7.0, 8.0],
            [0.0, 0.0, 9.0, 6.0],
            [0.0, 0.0, 6.0, 9.0],
        ]
    )
    groups = pd.Series(
        ["c0"] * 4 + ["c1"] * 4,
        index=[f"cell-{index}" for index in range(8)],
        name="cluster",
    )
    features = pd.Index(["A1", "A2", "B1", "B2"])
    programs = [
        MarkerProgram(
            label="Type A",
            parent="immune",
            positive=("A1", "A2", "NOT_ON_PANEL"),
            negative=("B1",),
        ),
        MarkerProgram(
            label="Type B",
            parent="immune",
            positive=("B1", "B2"),
            negative=("A1",),
        ),
    ]
    return matrix, groups, features, programs


def test_assembly_merges_standard_provider_and_signed_model_exports(
    tmp_path: Path,
) -> None:
    standard = tmp_path / "standard.json"
    standard.write_text(
        json.dumps(
            {
                "Type A": {
                    "parent": "immune",
                    "positive": ["A1"],
                    "negative": ["B1"],
                },
                "Unscoped": {"positive": ["HOUSE"]},
            }
        )
    )
    model_export = tmp_path / "model.tsv"
    pd.DataFrame(
        [
            {"cell_type": "Type A", "gene": "A2", "coefficient": 3.0, "parent": "immune"},
            {"cell_type": "Type A", "gene": "B2", "coefficient": -2.0, "parent": "immune"},
            {"cell_type": "Type C", "gene": "C1", "coefficient": 4.0, "parent": "stromal"},
        ]
    ).to_csv(model_export, sep="\t", index=False)

    programs = assemble_marker_programs(
        [standard, model_export],
        parent_scope="immune",
        source_version="fixture-v1",
    )

    assert [program.label for program in programs] == ["Type A"]
    assert programs[0].positive == ("A1", "A2")
    assert programs[0].negative == ("B1", "B2")
    assert "model-export:" in programs[0].source

    with_unscoped = assemble_marker_programs(
        [standard],
        parent_scope="immune",
        include_unscoped=True,
    )
    assert {program.label for program in with_unscoped} == {"Type A", "Unscoped"}


def test_assembly_accepts_json_coefficient_maps_and_labeled_references(
    tmp_path: Path,
) -> None:
    coefficient_map = tmp_path / "coefficients.json"
    coefficient_map.write_text(json.dumps({"Type A": {"A1": 2.0, "B1": -1.0}}))
    exported = assemble_marker_programs([coefficient_map], max_features_per_sign=1)
    assert exported[0].positive == ("A1",)
    assert exported[0].negative == ("B1",)

    reference_path = tmp_path / "reference.h5ad"
    reference = ad.AnnData(
        X=np.asarray([[5.0, 0.0], [4.0, 0.0], [0.0, 5.0], [0.0, 4.0]]),
        obs=pd.DataFrame(
            {"cell_type": pd.Categorical(["A", "A", "B", "B"])},
            index=["r1", "r2", "r3", "r4"],
        ),
        var=pd.DataFrame(index=["A1", "B1"]),
    )
    reference.write_h5ad(reference_path)
    programs = assemble_marker_programs(
        [ReferenceProvider(reference_path, label_key="cell_type", top_n=1)]
    )
    assert {program.label: program.positive for program in programs} == {
        "A": ("A1",),
        "B": ("B1",),
    }


def test_panel_profile_keeps_missing_markers_and_negative_violations_visible() -> None:
    matrix, groups, features, programs = _marker_fixture()
    profile = profile_marker_panel(
        matrix,
        groups,
        features,
        programs,
        detection=matrix,
        detection_threshold=0.25,
        min_positive_fraction=0.5,
    )

    assert profile.mean_matrix.loc["c0", "A1"] == 7.5
    assert profile.detection_matrix.loc["c0", "A1"] == 1.0
    assert np.isnan(profile.mean_matrix.loc["c0", "NOT_ON_PANEL"])
    type_a = profile.panel_capability.set_index("label").loc["Type A"]
    assert bool(type_a["callable"])
    assert type_a["positive_panel_fraction"] == 2 / 3
    assert type_a["missing_positive"] == "NOT_ON_PANEL"

    violations = profile.negative_violations.set_index(["cluster_id", "label", "feature"])
    assert violations.loc[("c0", "Type A", "B1"), "status"] == "not_violated"
    assert violations.loc[("c1", "Type A", "B1"), "status"] == "violated"
    assert violations.loc[("c0", "Type B", "A1"), "violated"] == np.bool_(True)


def test_aucell_ucell_de_overlap_and_consolidation_preserve_route_scores() -> None:
    matrix, groups, features, programs = _marker_fixture()
    result = run_marker_analysis(
        matrix,
        groups,
        features,
        programs,
        detection=matrix,
        detection_threshold=0.25,
        min_positive_fraction=0.5,
        aucell_max_rank_fraction=0.5,
        ucell_max_rank=4,
        de_top_n=4,
        adjusted_pvalue_max=1.0,
    )

    first_cell_auc = result.aucell[result.aucell["cell_id"].eq("cell-0")].set_index("label")
    first_cell_ucell = result.ucell[result.ucell["cell_id"].eq("cell-0")].set_index("label")
    assert (
        first_cell_auc.loc["Type A", "aucell_score"]
        > first_cell_auc.loc["Type B", "aucell_score"]
    )
    assert (
        first_cell_ucell.loc["Type A", "ucell_score"]
        > first_cell_ucell.loc["Type B", "ucell_score"]
    )

    overlaps = result.differential.overlap.set_index(["cluster_id", "label"])
    assert overlaps.loc[("c0", "Type A"), "positive_overlap_n"] == 2
    assert overlaps.loc[("c0", "Type B"), "negative_overlap_n"] == 1
    assert set(result.differential.bottom_up["cluster_id"]) == {"c0", "c1"}

    crosschecks = result.crosschecks.set_index(["cluster_id", "label"])
    assert (
        crosschecks.loc[("c0", "Type A"), "crosscheck_score"]
        > crosschecks.loc[("c0", "Type B"), "crosscheck_score"]
    )
    assert (
        crosschecks.loc[("c1", "Type B"), "crosscheck_score"]
        > crosschecks.loc[("c1", "Type A"), "crosscheck_score"]
    )
    assert crosschecks.loc[("c0", "Type A"), "routes_available"] == (
        "top_down;aucell;ucell;de_overlap"
    )


def test_uncallable_programs_emit_audit_rows_without_scores() -> None:
    matrix, groups, features, _ = _marker_fixture()
    program = MarkerProgram(label="Unavailable", positive=("MISSING",), negative=("A1",))
    profile = profile_marker_panel(matrix, groups, features, [program], detection=matrix)
    assert not bool(profile.panel_capability.loc[0, "callable"])
    assert (
        "insufficient_measured_positive_markers"
        in profile.panel_capability.loc[0, "uncallable_reasons"]
    )
    assert profile.feature_profile["measured"].eq(False).any()  # noqa: E712

    aucell = score_aucell(matrix, features, [program], cell_ids=groups.index)
    ucell = score_ucell(matrix, features, [program], cell_ids=groups.index)
    assert aucell["callable"].eq(False).all()  # noqa: E712
    assert ucell["callable"].eq(False).all()  # noqa: E712
    assert aucell["aucell_score"].isna().all()
    assert ucell["ucell_score"].isna().all()


def test_artifacts_are_durable_and_svg_contains_actual_data(tmp_path: Path) -> None:
    matrix, groups, features, programs = _marker_fixture()
    result = run_marker_analysis(
        matrix,
        groups,
        features,
        programs,
        detection=matrix,
        min_positive_fraction=0.5,
        adjusted_pvalue_max=1.0,
    )
    svg = marker_dotplot_svg(result.profile)
    assert svg.startswith("<svg")
    assert "<circle" in svg
    assert "cluster=c0; label=Type A; feature=A1" in svg
    assert "unmeasured" in svg

    artifacts = write_marker_artifacts(result, tmp_path / "markers", run_id="fixture-v1")
    resumed = write_marker_artifacts(result, tmp_path / "markers", run_id="fixture-v1")
    assert artifacts == resumed
    assert all(path.is_file() for path in artifacts.values())
    manifest = json.loads(artifacts["manifest"].read_text())
    assert manifest["artifact_kind"] == "cell_curator_marker_analysis"
    assert manifest["n_programs"] == 2
    for name, payload in manifest["artifacts"].items():
        assert payload["sha256"] == sha256(artifacts[name])
    assert pd.read_csv(artifacts["crosschecks"], sep="\t").shape[0] == 4


def test_sparse_matrices_use_the_same_marker_routes() -> None:
    matrix, groups, features, programs = _marker_fixture()
    result = run_marker_analysis(
        sparse.csr_matrix(matrix),
        groups,
        features,
        programs,
        detection=sparse.csr_matrix(matrix),
        min_positive_fraction=0.5,
        adjusted_pvalue_max=1.0,
    )
    assert len(result.aucell) == len(groups) * len(programs)
    assert len(result.ucell) == len(groups) * len(programs)
    assert result.profile.detection_matrix.loc["c1", "B1"] == 1.0
    assert result.crosschecks["crosscheck_score"].notna().all()


def test_config_wrapper_scopes_cells_and_preserves_source(
    tmp_path: Path,
    config: dict,
) -> None:
    matrix, groups, features, _ = _marker_fixture()
    source = tmp_path / "source.h5ad"
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            {
                "Type A": {
                    "parent": "immune",
                    "positive": ["A1", "A2"],
                    "negative": ["B1"],
                },
                "Type B": {
                    "parent": "immune",
                    "positive": ["B1", "B2"],
                    "negative": ["A1"],
                },
            }
        )
    )
    obs = pd.DataFrame(
        {
            "canonical_cluster": pd.Categorical(["immune"] * 8),
            "lineage": pd.Categorical(["immune"] * 8),
            "child_cluster": pd.Categorical(groups.to_numpy()),
        },
        index=groups.index,
    )
    adata = ad.AnnData(
        X=matrix.copy(),
        obs=obs,
        var=pd.DataFrame({"gene_symbol": features}, index=features),
    )
    adata.layers["logcounts"] = matrix.copy()
    adata.layers["counts"] = matrix.copy()
    adata.obsm["X_pca_harmony"] = np.c_[np.arange(8), np.arange(8)]
    adata.write_h5ad(source)
    source_before = sha256(source)

    config["input"]["path"] = str(source)
    config["input"]["expected_sha256"] = source_before
    config["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": str(markers), "version": "fixture-v1"}
    ]
    model = CellCuratorConfig.model_validate(config)
    run = run_marker_analysis_from_config(
        model,
        assay="transcriptome",
        representation="logcounts",
        group_key="child_cluster",
        parent_scope="immune",
        output_dir=tmp_path / "marker-output",
        adjusted_pvalue_max=1.0,
    )

    assert len(run.analysis.crosschecks) == 4
    assert run.artifacts["manifest"].is_file()
    assert sha256(source) == source_before

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from conftest import SKILL_ROOT

from cell_curator.config import load_config
from cell_curator.contracts import ContractError, sha256
from cell_curator.guided import (
    approve_assumptions,
    confirm_context,
    construct_plan,
    prepare_scene,
    record_assumption,
    run_directory,
)
from cell_curator.hierarchy import LABEL_COLUMNS, scaffold_labels_table
from cell_curator.pipeline import execute_recursive_hierarchy, initialize


def _write_fixture(
    tmp_path: Path,
    *,
    run_id: str = "recursive-v1",
    interactive: bool = False,
    composition_covariates: list[str] | None = None,
) -> tuple[Path, Path, Path]:
    rng = np.random.default_rng(331)
    genes = ["IMM", "STR", "T", "B", "FIB", "ENDO", "CD4", "CD8", "MB", "PL"]
    counts = rng.poisson(3, size=(24, len(genes))).astype(np.float32)
    obs = pd.DataFrame(
        {
            "canonical_cluster": ["0"] * 12 + ["1"] * 12,
            "lineage": pd.Categorical(["immune"] * 12 + ["stromal"] * 12),
            "numeric_depth": np.linspace(0.125, 9.875, 24),
        },
        index=[f"recursive-{index:03d}" for index in range(24)],
    )
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=obs,
        var=pd.DataFrame({"gene_symbol": genes}, index=genes),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = rng.normal(size=(24, 4)).astype(np.float32)
    source = tmp_path / f"{run_id}.h5ad"
    adata.write_h5ad(source)

    programs = [
        {"label": "Immune", "pos": ["IMM"], "neg": ["STR"]},
        {"label": "Stromal", "pos": ["STR"], "neg": ["IMM"]},
        {"label": "T", "parent": "Immune", "pos": ["T"], "neg": ["B"]},
        {"label": "B", "parent": "Immune", "pos": ["B"], "neg": ["T"]},
        {
            "label": "Fibroblast",
            "parent": "Stromal",
            "pos": ["FIB"],
            "neg": ["ENDO"],
        },
        {
            "label": "Endothelial",
            "parent": "Stromal",
            "pos": ["ENDO"],
            "neg": ["FIB"],
        },
        {"label": "CD4", "parent": "T", "pos": ["CD4"], "neg": ["CD8"]},
        {"label": "CD8", "parent": "T", "pos": ["CD8"], "neg": ["CD4"]},
        {"label": "MemoryB", "parent": "B", "pos": ["MB"], "neg": ["PL"]},
        {"label": "Plasma", "parent": "B", "pos": ["PL"], "neg": ["MB"]},
    ]
    markers = tmp_path / f"{run_id}.markers.json"
    markers.write_text(json.dumps(programs))

    value = yaml.safe_load((SKILL_ROOT / "assets/config.template.yaml").read_text())
    value["run"].update(
        {
            "run_id": run_id,
            "output_root": str(tmp_path / "results"),
            "mode": "annotation",
        }
    )
    value["input"]["path"] = str(source)
    value["canonical_partitions"] = {
        "lineage": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    value["adaptive"]["eligibility"].update(
        {"min_cells": 2, "min_donors": None, "min_captures": None}
    )
    value["adaptive"]["parent_local"]["enabled"] = False
    value["knowledge"].update(
        {
            "providers": [
                {
                    "kind": "local_markers",
                    "path": str(markers),
                    "version": "recursive-fixture-v1",
                }
            ],
            "allow_optional_remote_enhancement": False,
        }
    )
    covariates = composition_covariates or ["lineage"]
    value["guidance"] = {
        "interactive": interactive,
        "preauthorized": not interactive,
        "reviewer": "" if interactive else "recursive-fixture-reviewer",
        "pause_after_each_level": interactive,
        "assumptions": ["The synthetic hierarchy and score tables are authorized."],
        "hierarchy_levels": ["L1", "L2", "L3"],
        "context": {
            "organism": "human",
            "tissue": "synthetic tissue",
            "developmental_stage": "adult",
            "condition": "healthy",
            "technology": "dissociated",
            "experimental_context": "recursive hierarchy integration fixture",
            "composition_covariates": covariates,
            "vocabulary_contract": "local hierarchical signed programs",
            "reference_contract": "local_only",
            "visualization_key": "X_pca_harmony",
            "writeback_prefix": "recursive_label",
        },
    }
    config_path = tmp_path / f"{run_id}.yaml"
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))
    return config_path, source, markers


# The L1 cluster each fixture subcluster refines; a child row must name a cluster its
# parent level actually routed into that scope.
SUBCLUSTER_PARENT = {
    "immune-t": "root-immune",
    "immune-b": "root-immune",
    "immune-cycling": "root-immune",
    "stromal-f": "root-stromal",
    "t-cd4": "immune-t",
}


def _scores(cluster: str, winner: str, runner: str) -> list[dict[str, object]]:
    shared = {"specificity": 0.95, "agreement": 0.95, "coverage": 0.95}
    if cluster in SUBCLUSTER_PARENT:
        shared["parent_cluster_id"] = SUBCLUSTER_PARENT[cluster]
    return [
        {"cluster_id": cluster, "label": winner, "score": 3.0, **shared},
        {"cluster_id": cluster, "label": runner, "score": 0.1, **shared},
    ]


def _score_tables() -> dict[tuple[str, str], pd.DataFrame]:
    return {
        ("L1", "ROOT"): pd.DataFrame(
            [
                *_scores("root-immune", "Immune", "Stromal"),
                *_scores("root-stromal", "Stromal", "Immune"),
            ]
        ),
        ("L2", "Immune"): pd.DataFrame(
            [*_scores("immune-t", "T", "B"), *_scores("immune-b", "B", "T")]
        ),
        ("L2", "Stromal"): pd.DataFrame(_scores("stromal-f", "Fibroblast", "Endothelial")),
        ("L3", "T"): pd.DataFrame(_scores("t-cd4", "CD4", "CD8")),
    }


def _label_row(
    *, level: str, parent: str, cluster: str, label: str, reviewer: str = ""
) -> dict[str, object]:
    return {
        "level": level,
        "parent_scope": parent,
        "cluster_id": cluster,
        "label": label,
        "cell_state": "normal",
        "cl_id": "",
        "confidence": 0.8,
        "margin": 0.7,
        "specificity": 0.8,
        "agreement": 0.8,
        "coverage": 0.8,
        "rationale": f"Automatic rationale for {label}",
        "rejected_alternatives": "alternative",
        "reviewer": reviewer,
    }


def _decisions(
    tables: dict[tuple[str, str], pd.DataFrame],
) -> dict[tuple[str, str], pd.DataFrame]:
    """One label row per cluster: the label with the highest score, as a stand-in consensus."""
    return {
        (level, parent): pd.DataFrame(
            [
                _label_row(level=level, parent=parent, cluster=str(cluster), label=block.loc[block["score"].idxmax(), "label"])
                for cluster, block in table.groupby("cluster_id", sort=True)
            ]
        )
        for (level, parent), table in tables.items()
    }


def test_plan_routes_each_parent_scope_and_rejects_noncategorical_covariates(
    tmp_path: Path,
) -> None:
    config_path, _, _ = _write_fixture(tmp_path)
    prepare_scene(config_path)
    plan = construct_plan(config_path)
    observed = {(route["level"], route["parent_scope"]) for route in plan["routes"]}
    assert observed == {
        ("L1", "ROOT"),
        ("L2", "Immune"),
        ("L2", "Stromal"),
        ("L3", "B"),
        ("L3", "T"),
    }

    invalid_path, _, _ = _write_fixture(
        tmp_path,
        run_id="invalid-covariate-v1",
        interactive=True,
        composition_covariates=["numeric_depth"],
    )
    scene = prepare_scene(invalid_path)
    assert scene["context_complete"] is False
    assert "composition_covariate:numeric_depth" in scene["unresolved"]
    confirmation = tmp_path / "invalid-context-confirmation.json"
    confirmation.write_text(json.dumps(scene["biological_context"]))
    with pytest.raises(ContractError, match="absent/non-categorical"):
        confirm_context(
            invalid_path,
            input_path=confirmation,
            reviewer="context-reviewer",
        )


def test_labels_scaffold_accumulates_scopes_and_preserves_manual_rows(
    tmp_path: Path,
) -> None:
    path = tmp_path / "labels" / "L2.tsv"
    immune = pd.DataFrame(
        [_label_row(level="L2", parent="Immune", cluster="i0", label="T")],
        columns=LABEL_COLUMNS,
    )
    scaffold_labels_table(path, immune)
    manually_edited = pd.read_csv(path, sep="\t", keep_default_na=False)
    manually_edited.loc[0, ["label", "rationale", "reviewer"]] = [
        "Manual T",
        "Reviewer-entered biological rationale.",
        "reviewer-a",
    ]
    path.write_text(manually_edited.to_csv(sep="\t", index=False))

    stromal = pd.DataFrame(
        [_label_row(level="L2", parent="Stromal", cluster="s0", label="Fibroblast")],
        columns=LABEL_COLUMNS,
    )
    accumulated = scaffold_labels_table(path, stromal)
    assert set(accumulated["parent_scope"]) == {"Immune", "Stromal"}

    rebuilt = scaffold_labels_table(
        path,
        pd.DataFrame(
            [
                _label_row(level="L2", parent="Immune", cluster="i0", label="B"),
                _label_row(level="L2", parent="Stromal", cluster="s0", label="Endothelial"),
            ],
            columns=LABEL_COLUMNS,
        ),
    ).set_index(["parent_scope", "cluster_id"])
    assert rebuilt.loc[("Immune", "i0"), "label"] == "Manual T"
    assert (
        rebuilt.loc[("Immune", "i0"), "rationale"] == "Reviewer-entered biological rationale."
    )
    assert rebuilt.loc[("Stromal", "s0"), "label"] == "Endothelial"

    duplicate = pd.concat([stromal, stromal], ignore_index=True)
    with pytest.raises(ContractError, match="duplicate level/parent/cluster"):
        scaffold_labels_table(path, duplicate)


def test_recursive_execution_three_levels_resume_and_scope_regating(
    tmp_path: Path,
) -> None:
    config_path, source, _ = _write_fixture(tmp_path)
    source_digest = sha256(source)
    prepare_scene(config_path)
    construct_plan(config_path)
    score_tables = _score_tables()

    first = execute_recursive_hierarchy(
        config_path, score_tables=score_tables, decisions=_decisions(score_tables)
    )
    assert set(first["executed_scopes"]) == {
        "L1::ROOT",
        "L2::Immune",
        "L2::Stromal",
        "L3::T",
    }
    assert first["hierarchy_validation"]["status"] == "complete"
    assert first["retained_leaves"]["missing_child_score_tables"] == 1
    root = run_directory(load_config(config_path))
    level_two = pd.read_csv(root / "labels/L2.tsv", sep="\t", keep_default_na=False)
    assert set(level_two["parent_scope"]) == {"Immune", "Stromal"}
    retained = pd.read_csv(
        root / "hierarchy/retained_leaves.tsv", sep="\t", keep_default_na=False
    )
    missing_b = retained[retained["reason"].eq("missing_child_score_table")]
    assert missing_b[["requested_level", "requested_parent_scope"]].values.tolist() == [
        ["L3", "B"]
    ]
    assert missing_b.iloc[0]["leaf_label"] == "B"
    assert sha256(source) == source_digest

    checkpoints = {
        identity: Path(value["path"]) for identity, value in first["checkpoints"].items()
    }
    before = {
        identity: (sha256(path), path.stat().st_mtime_ns)
        for identity, path in checkpoints.items()
    }
    second = execute_recursive_hierarchy(
        config_path, score_tables=score_tables, decisions=_decisions(score_tables)
    )
    assert second["executed_scopes"] == [], "valid recursive checkpoints were recomputed"
    assert set(second["resumed_scopes"]) == set(checkpoints)
    assert {
        identity: (sha256(path), path.stat().st_mtime_ns)
        for identity, path in checkpoints.items()
    } == before

    record_assumption(
        config_path,
        scope="L3::T",
        statement="A material CD4/CD8 interpretation change requires renewed L3 review.",
        material=True,
    )
    with pytest.raises(ContractError, match="Gate A blocks biological labels"):
        execute_recursive_hierarchy(
            config_path, score_tables=score_tables, decisions=_decisions(score_tables)
        )
    for identity in ("L1::ROOT", "L2::Immune", "L2::Stromal"):
        assert sha256(checkpoints[identity]) == before[identity][0]

    approve_assumptions(
        config_path,
        reviewer="recursive-fixture-reviewer",
        scope="L3::T",
        preauthorized=True,
    )
    third = execute_recursive_hierarchy(
        config_path, score_tables=score_tables, decisions=_decisions(score_tables)
    )
    assert third["executed_scopes"] == ["L3::T"]
    assert set(third["resumed_scopes"]) == {
        "L1::ROOT",
        "L2::Immune",
        "L2::Stromal",
    }
    for identity in ("L1::ROOT", "L2::Immune", "L2::Stromal"):
        assert sha256(checkpoints[identity]) == before[identity][0]
    assert sha256(checkpoints["L3::T"]) != before["L3::T"][0]


def test_provider_manifest_hashes_local_files_and_disabled_remote_never_networks(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    config_path, _, markers = _write_fixture(tmp_path, run_id="provider-provenance-v1")
    value = yaml.safe_load(config_path.read_text())
    value["knowledge"]["providers"].append(
        {
            "kind": "remote_table",
            "url": "https://invalid.example/markers.tsv",
            "cache_path": str(tmp_path / "absent-remote-cache.tsv"),
            "sha256": "0" * 64,
            "version": "unavailable-fixture-v1",
            "timeout_seconds": 0.1,
        }
    )
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))

    def reject_network(*_args: object, **_kwargs: object) -> object:
        pytest.fail("disabled optional remote enhancement attempted network I/O")

    monkeypatch.setattr("cell_curator.net._open", reject_network)
    initialize(config_path)
    root = run_directory(load_config(config_path))
    manifest = json.loads((root / "knowledge.manifest.json").read_text())
    providers = {item["kind"]: item for item in manifest["providers"]}
    assert providers["local_markers"]["artifact"] == {
        "status": "validated",
        "path": str(markers.resolve()),
        "sha256": sha256(markers),
        "declared_checksum": False,
        "checksum_validated": False,
    }
    assert providers["remote_table"]["artifact"]["status"] == "unavailable"
    assert "disabled" in providers["remote_table"]["artifact"]["reason"]


def test_a_state_program_is_never_assignable_as_an_identity(tmp_path: Path) -> None:
    """The planned vocabulary records every program; assignment takes only identity rows.

    `plan/vocabulary.tsv` is written from all declared marker programs, cell state
    included. The walker reduces it to identity rows and refuses a decision that names
    a cell state program as a cell type.
    """

    config_path, _, markers = _write_fixture(tmp_path, run_id="state-vocabulary-v1")
    programs = json.loads(markers.read_text())
    programs.append(
        {"label": "Cycling", "parent": "Immune", "pos": ["MB"], "axis": "cell_state"}
    )
    markers.write_text(json.dumps(programs))

    prepare_scene(config_path)
    construct_plan(config_path)
    root = run_directory(load_config(config_path))

    planned = pd.read_csv(root / "plan/vocabulary.tsv", sep="\t")
    assert "Cycling" in set(planned["label"]), "the plan must record what was declared"

    tables = _score_tables()
    tables[("L2", "Immune")] = pd.DataFrame(
        [
            *_scores("immune-t", "T", "B"),
            *_scores("immune-b", "B", "T"),
            *_scores("immune-cycling", "Cycling", "T"),
        ]
    )
    with pytest.raises(ContractError, match=r"outside its identity vocabulary: \['Cycling'\]"):
        execute_recursive_hierarchy(
            config_path, score_tables=tables, decisions=_decisions(tables)
        )
    receipt = json.loads((root / "plan/identity_vocabulary_status.json").read_text())
    assert receipt["excluded_cell_state_labels"] == ["Cycling"]


def test_an_operator_vocabulary_label_with_no_approved_assumption_is_refused(
    tmp_path: Path,
) -> None:
    """Gate A fails closed: a label nobody approved cannot be assigned."""

    config_path, _, _ = _write_fixture(tmp_path, run_id="unapproved-vocabulary-v1")
    prepare_scene(config_path)
    construct_plan(config_path)
    root = run_directory(load_config(config_path))
    vocabulary = pd.read_csv(root / "plan/vocabulary.tsv", sep="\t", keep_default_na=False)
    extra = vocabulary[vocabulary["level"].eq("L1")].iloc[[0]].copy()
    extra["label"] = "Neuron"
    operator = tmp_path / "operator-vocabulary.tsv"
    pd.concat([vocabulary, extra]).to_csv(operator, sep="\t", index=False)
    with pytest.raises(ContractError, match="no approved assumption naming them"):
        execute_recursive_hierarchy(
            config_path,
            score_tables=_score_tables(),
            decisions=_decisions(_score_tables()),
            vocabulary_path=operator,
        )


def test_a_changed_decision_reexecutes_and_a_dropped_level_loses_its_table(
    tmp_path: Path,
) -> None:
    config_path, _, _ = _write_fixture(tmp_path, run_id="decision-resume-v1")
    prepare_scene(config_path)
    construct_plan(config_path)
    root = run_directory(load_config(config_path))
    tables = _score_tables()
    decisions = _decisions(tables)
    first = execute_recursive_hierarchy(config_path, score_tables=tables, decisions=decisions)
    assert (root / "labels/L3.tsv").is_file()

    changed = {scope: frame.copy() for scope, frame in decisions.items()}
    changed[("L1", "ROOT")]["rationale"] = "A reviewer reworded the rationale."
    again = execute_recursive_hierarchy(config_path, score_tables=tables, decisions=changed)
    assert again["executed_scopes"] == ["L1::ROOT"], "a changed decision must not resume old labels"

    without_l3 = {scope: frame for scope, frame in tables.items() if scope[0] != "L3"}
    execute_recursive_hierarchy(
        config_path, score_tables=without_l3, decisions=_decisions(without_l3)
    )
    assert not (root / "labels/L3.tsv").exists()
    assert first["executed_scopes"]


def test_a_child_row_for_a_cluster_its_parent_never_routed_is_refused(
    tmp_path: Path,
) -> None:
    config_path, _, _ = _write_fixture(tmp_path, run_id="unrouted-child-v1")
    prepare_scene(config_path)
    construct_plan(config_path)
    tables = _score_tables()
    # root-stromal went to Stromal at L1, so it cannot be refined under Immune.
    tables[("L2", "Immune")] = pd.DataFrame(
        [{**row, "parent_cluster_id": "root-stromal"} for row in _scores("x", "T", "B")]
    )
    with pytest.raises(ContractError, match="did not route there"):
        execute_recursive_hierarchy(
            config_path, score_tables=tables, decisions=_decisions(tables)
        )

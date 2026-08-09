from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from conftest import SKILL_ROOT

import cell_curator.hierarchy as hierarchy_module
from cell_curator import execute_recursive_hierarchy_manifest as public_manifest_execute
from cell_curator.cli import build_parser, dispatch
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


def _scores(cluster: str, winner: str, runner: str) -> list[dict[str, object]]:
    shared = {"specificity": 0.95, "agreement": 0.95, "coverage": 0.95}
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
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    config_path, source, _ = _write_fixture(tmp_path)
    source_digest = sha256(source)
    prepare_scene(config_path)
    construct_plan(config_path)
    score_tables = _score_tables()

    first = execute_recursive_hierarchy(config_path, score_tables=score_tables)
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
    with monkeypatch.context() as resume_patch:
        resume_patch.setattr(
            hierarchy_module,
            "assign_parent_scoped_level",
            lambda *_args, **_kwargs: pytest.fail(
                "valid recursive checkpoints were recomputed"
            ),
        )
        second = execute_recursive_hierarchy(config_path, score_tables=score_tables)
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
    original_assign = hierarchy_module.assign_parent_scoped_level
    recomputed: list[tuple[str, str]] = []

    def spy_assign(*args: object, **kwargs: object) -> pd.DataFrame:
        recomputed.append((str(kwargs["level"]), str(kwargs["parent_scope"])))
        return original_assign(*args, **kwargs)

    monkeypatch.setattr(hierarchy_module, "assign_parent_scoped_level", spy_assign)
    with pytest.raises(ContractError, match="Gate A blocks biological labels"):
        execute_recursive_hierarchy(config_path, score_tables=score_tables)
    assert recomputed == []
    for identity in ("L1::ROOT", "L2::Immune", "L2::Stromal"):
        assert sha256(checkpoints[identity]) == before[identity][0]

    approve_assumptions(
        config_path,
        reviewer="recursive-fixture-reviewer",
        scope="L3::T",
        preauthorized=True,
    )
    third = execute_recursive_hierarchy(config_path, score_tables=score_tables)
    assert recomputed == [("L3", "T")]
    assert set(third["resumed_scopes"]) == {
        "L1::ROOT",
        "L2::Immune",
        "L2::Stromal",
    }
    for identity in ("L1::ROOT", "L2::Immune", "L2::Stromal"):
        assert sha256(checkpoints[identity]) == before[identity][0]
    assert sha256(checkpoints["L3::T"]) != before["L3::T"][0]


def test_recursive_manifest_public_api_and_cli(tmp_path: Path) -> None:
    config_path, source, _ = _write_fixture(tmp_path, run_id="recursive-manifest-v1")
    source_digest = sha256(source)
    prepare_scene(config_path)
    construct_plan(config_path)
    entries: list[dict[str, str]] = []
    for index, ((level, parent), frame) in enumerate(_score_tables().items()):
        score_path = tmp_path / f"scope-{index}.tsv"
        frame.to_csv(score_path, sep="\t", index=False)
        entries.append({"level": level, "parent_scope": parent, "path": score_path.name})
    manifest_path = tmp_path / "recursive-scores.json"
    manifest_path.write_text(
        json.dumps({"schema_version": 3, "score_tables": entries}, indent=2)
    )

    args = build_parser().parse_args(
        [
            "run",
            "execute",
            "--config",
            str(config_path),
            "--score-manifest",
            str(manifest_path),
        ]
    )
    executed = dispatch(args)
    assert set(executed["executed_scopes"]) == {
        "L1::ROOT",
        "L2::Immune",
        "L2::Stromal",
        "L3::T",
    }
    resumed = public_manifest_execute(config_path, manifest_path)
    assert set(resumed["resumed_scopes"]) == set(executed["executed_scopes"])
    assert sha256(source) == source_digest

    invalid_path = tmp_path / "invalid-recursive-scores.json"
    invalid_path.write_text(
        json.dumps(
            {
                "schema_version": 3,
                "score_tables": entries,
                "unsupported": True,
            }
        )
    )
    with pytest.raises(ContractError, match="unsupported fields"):
        public_manifest_execute(config_path, invalid_path)


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

    monkeypatch.setattr("urllib.request.urlopen", reject_network)
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

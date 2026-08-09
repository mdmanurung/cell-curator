from __future__ import annotations

import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from conftest import REPO_ROOT, SKILL_ROOT

import cell_curator.crosschecks as crosschecks_module
import cell_curator.environment as environment_module
from cell_curator.cli import build_parser, dispatch
from cell_curator.config import load_config
from cell_curator.contracts import ContractError
from cell_curator.environment import ComputeProfile, detect_environment
from cell_curator.guided import (
    approve_assumptions,
    confirm_context,
    construct_plan,
    load_assumptions,
    load_state,
    prepare_scene,
    record_assumption,
    require_label_gate,
    run_directory,
)
from cell_curator.pipeline import run_annotation
from cell_curator.routing import (
    CrosscheckMethod,
    HealthContext,
    MethodResource,
    RoutingRequest,
    Technology,
)

REPOSITORY_ROOT = REPO_ROOT
SKILL_PATH = SKILL_ROOT / "SKILL.md"
README_PATH = REPOSITORY_ROOT / "README.md"
DOCUMENTED_COMMAND_PATHS = sorted(
    {
        tuple(part for part in match.groups() if part)
        for match in re.finditer(
            r"^\s*(?:uv run )?cell-curator\s+([A-Za-z][\w-]*)(?:\s+([A-Za-z][\w-]*))?",
            SKILL_PATH.read_text() + "\n" + README_PATH.read_text(),
            re.MULTILINE,
        )
    }
)


def _parse(*arguments: str):
    return build_parser().parse_args(list(arguments))


def _write_interactive_fixture(tmp_path: Path) -> Path:
    rng = np.random.default_rng(91)
    counts = np.zeros((60, 6), dtype=np.float32)
    counts[:30, :2] = rng.poisson(8, size=(30, 2)) + 1
    counts[30:, 2:4] = rng.poisson(8, size=(30, 2)) + 1
    obs = pd.DataFrame(
        {
            "canonical_cluster": ["0"] * 30 + ["1"] * 30,
            "lineage": ["immune"] * 60,
        },
        index=[f"guided-{index:03d}" for index in range(60)],
    )
    genes = ["A1", "A2", "B1", "B2", "N1", "N2"]
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=obs,
        var=pd.DataFrame({"gene_symbol": genes}, index=genes),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = np.vstack(
        [rng.normal(-3, 0.2, size=(30, 3)), rng.normal(3, 0.2, size=(30, 3))]
    ).astype(np.float32)
    source = tmp_path / "interactive-source.h5ad"
    adata.write_h5ad(source)

    markers = tmp_path / "interactive-markers.json"
    markers.write_text(
        json.dumps(
            {
                "Type A": {"pos": ["A1", "A2"], "neg": ["B1", "B2"]},
                "Type B": {"pos": ["B1", "B2"], "neg": ["A1", "A2"]},
            }
        )
    )
    value = yaml.safe_load((SKILL_ROOT / "assets/config.template.yaml").read_text())
    value["run"].update(
        {
            "run_id": "interactive-gate-v1",
            "output_root": str(tmp_path / "results"),
            "mode": "annotation",
        }
    )
    value["input"]["path"] = str(source)
    value["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    value["adaptive"]["eligibility"].update(
        {"min_cells": 5, "min_donors": None, "min_captures": None}
    )
    value["adaptive"]["parent_local"]["enabled"] = False
    value["knowledge"].update(
        {
            "providers": [
                {
                    "kind": "local_markers",
                    "path": str(markers),
                    "version": "interactive-fixture-v1",
                }
            ],
            "allow_optional_remote_enhancement": False,
        }
    )
    value["guidance"] = {
        "interactive": True,
        "preauthorized": False,
        "reviewer": "",
        "pause_after_each_level": True,
        "assumptions": ["The local signed marker programs define the L1 vocabulary."],
        "hierarchy_levels": ["L1"],
        "context": {
            "organism": "human",
            "tissue": "synthetic tissue",
            "developmental_stage": "adult",
            "condition": "healthy",
            "technology": "dissociated",
            "experimental_context": "interactive synthetic gate fixture",
            "composition_covariates": ["lineage"],
            "vocabulary_contract": "local signed marker programs",
            "reference_contract": "local_only",
            "visualization_key": "X_pca_harmony",
            "writeback_prefix": "guided_label",
        },
    }
    config_path = tmp_path / "interactive-config.yaml"
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))
    return config_path


@pytest.mark.parametrize(
    "command_path",
    DOCUMENTED_COMMAND_PATHS,
    ids=lambda value: "-".join(value),
)
def test_every_skill_command_path_resolves_help(
    command_path: tuple[str, ...],
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert len(DOCUMENTED_COMMAND_PATHS) >= 40, "SKILL command extraction unexpectedly shrank"
    with pytest.raises(SystemExit) as stopped:
        build_parser().parse_args([*command_path, "--help"])
    assert stopped.value.code == 0
    rendered = capsys.readouterr().out
    assert "usage: cell-curator" in rendered
    assert all(part in rendered for part in command_path)


def test_lightweight_cli_dispatch_smoke(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    environment = dispatch(_parse("environment", "detect", "--no-cuda-validation"))
    assert environment.detection_only is True
    assert environment.compute_profile in set(ComputeProfile)
    assert environment.usable_cpu_count >= 1
    assert environment.usable_ram_bytes > 0

    request = RoutingRequest(
        level="L1",
        parent_scope="ROOT",
        organism="human",
        tissue="blood",
        developmental_stage="adult",
        health_context=HealthContext.HEALTHY,
        technology=Technology.DISSOCIATED,
        compute_profile=ComputeProfile.LOCAL_CPU,
        usable_ram_bytes=8 * 1024**3,
        existing_predictions=(
            MethodResource(
                name="existing_prediction",
                organisms=("human",),
                tissues=("blood",),
                levels=("L1",),
            ),
        ),
    )
    request_path = tmp_path / "route-request.json"
    request_path.write_text(request.model_dump_json())
    route = dispatch(_parse("route", "select", "--request", str(request_path)))
    assert route.status == "selected"
    assert route.primary_method is CrosscheckMethod.EXISTING_PREDICTION

    confidence = dispatch(
        _parse(
            "cluster",
            "confidence",
            "--margin",
            "0.8",
            "--specificity",
            "0.9",
            "--agreement",
            "0.85",
            "--coverage",
            "0.95",
        )
    )
    assert 0.0 < confidence["advisory_confidence"] < 1.0
    assert confidence["calibrated_probability"] is False

    monkeypatch.setattr(crosschecks_module.importlib.util, "find_spec", lambda _name: None)
    provider_dir = tmp_path / "celltypist-status"
    unavailable = dispatch(
        _parse(
            "crosscheck",
            "celltypist-discover",
            "--output-dir",
            str(provider_dir),
        )
    )
    assert unavailable["status"] == "unavailable"
    assert "celltypist" in unavailable["reason"]
    assert "Install cell-curator[celltypist]" in unavailable["remediation"]
    assert (provider_dir / "provider_status.json").is_file()


def test_interactive_scene_plan_approval_and_material_assumption_revocation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    config_path = _write_interactive_fixture(tmp_path)
    fixed_environment = detect_environment(
        environ={},
        which=lambda _command: None,
        system="Linux",
        machine="x86_64",
        cpu_count=4,
        total_ram_bytes=8 * 1024**3,
        validate_cuda=False,
    )
    monkeypatch.setattr(
        environment_module,
        "detect_environment",
        lambda **_kwargs: fixed_environment,
    )

    scene = prepare_scene(config_path)
    root = run_directory(load_config(config_path))
    scene_state = load_state(root)
    assert scene["context_complete"] is True
    assert scene_state.stage == "scene"
    assert scene_state.labels_permitted is False
    assert not (root / "labels").exists()
    manifest = json.loads((root / "input_contract.manifest.json").read_text())
    canonical = ad.read_h5ad(manifest["canonical_path"])
    assert not any(column.startswith("guided_label") for column in canonical.obs.columns)

    plan = construct_plan(config_path)
    planned_state = load_state(root)
    assert plan["executable"] is True
    assert planned_state.plan_complete is True
    assert planned_state.labels_permitted is False
    with pytest.raises(ContractError, match="Gate A blocks biological labels"):
        require_label_gate(config_path, level="L1", parent="ROOT")

    approval = approve_assumptions(
        config_path,
        reviewer="interactive-reviewer",
        scope="*",
    )
    assert approval["all_material_approved"] is True
    gate_state = require_label_gate(config_path, level="L1", parent="ROOT")
    assert gate_state.labels_permitted is True
    assert gate_state.approved_assumptions_sha256
    assert all(
        row["approval_mode"] == "interactive"
        for row in load_assumptions(root)
        if row["material"]
    )

    changed = record_assumption(
        config_path,
        scope="L1::ROOT",
        statement="A newly observed marker conflict changes the permitted L1 interpretation.",
        material=True,
    )
    revoked_state = load_state(root)
    assert changed["approval_revoked"] is True
    assert revoked_state.assumptions_approved is False
    assert revoked_state.labels_permitted is False
    assert revoked_state.approved_assumptions_sha256 == ""
    with pytest.raises(ContractError, match="Gate A blocks biological labels"):
        require_label_gate(config_path, level="L1", parent="ROOT")


def test_confirmed_context_survives_execute_resume(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    config_path = _write_interactive_fixture(tmp_path)
    value = yaml.safe_load(config_path.read_text())
    value["run"]["run_id"] = "confirmed-context-v1"
    value["guidance"]["context"].update(
        {
            "organism": "",
            "tissue": "",
            "developmental_stage": "",
            "condition": "",
            "technology": "other",
            "experimental_context": "",
            "vocabulary_contract": "",
        }
    )
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))
    fixed_environment = detect_environment(
        environ={},
        which=lambda _command: None,
        system="Linux",
        machine="x86_64",
        cpu_count=4,
        total_ram_bytes=8 * 1024**3,
        validate_cuda=False,
    )
    monkeypatch.setattr(
        environment_module,
        "detect_environment",
        lambda **_kwargs: fixed_environment,
    )

    scene = prepare_scene(config_path)
    assert scene["context_complete"] is False
    confirmation_path = tmp_path / "confirmed-context.json"
    confirmation_path.write_text(
        json.dumps(
            {
                "organism": "human",
                "tissue": "reviewer-corrected tissue",
                "developmental_stage": "adult",
                "condition": "healthy",
                "technology": "dissociated",
                "experimental_context": "reviewer-confirmed fixture",
                "composition_covariates": ["lineage"],
                "vocabulary_contract": "reviewer-approved local markers",
                "reference_contract": "local_only",
                "visualization_key": "X_pca_harmony",
                "writeback_prefix": "guided_label",
            }
        )
    )
    confirmed = confirm_context(
        config_path,
        input_path=confirmation_path,
        reviewer="context-reviewer",
    )
    construct_plan(config_path)
    approve_assumptions(config_path, reviewer="context-reviewer")
    packet = run_annotation(config_path)
    persisted = json.loads((packet.parent / "context.json").read_text())
    assert confirmed["confirmation"] == persisted["confirmation"]
    assert persisted["biological_context"]["tissue"] == "reviewer-corrected tissue"

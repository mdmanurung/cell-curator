from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

import cell_curator
from cell_curator.benchmark import validate_dataset_manifest
from cell_curator.cli import build_parser
from cell_curator.config import CellCuratorConfig, configuration_schema
from cell_curator.contracts import ContractError, sha256


def _leaf_command_paths(
    parser: argparse.ArgumentParser,
    prefix: tuple[str, ...] = (),
) -> list[list[str]]:
    subparsers = [
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    ]
    if not subparsers:
        return [list(prefix)]
    paths: list[list[str]] = []
    for name, child in subparsers[0].choices.items():
        paths.extend(_leaf_command_paths(child, (*prefix, name)))
    return paths


@pytest.mark.parametrize(
    "arguments",
    [["--help"], ["--version"]]
    + [path + ["--help"] for path in _leaf_command_paths(build_parser())],
)
def test_every_cli_command_has_help(arguments: list[str]) -> None:
    with pytest.raises(SystemExit) as stopped:
        build_parser().parse_args(arguments)
    assert stopped.value.code == 0


def test_stable_public_api_is_exposed() -> None:
    expected = {
        "annotation_progress",
        "approve_annotation_assumptions",
        "inspect_input",
        "canonicalize_input",
        "compute_input_qc",
        "detect_compute_environment",
        "execute_recursive_hierarchy",
        "execute_recursive_hierarchy_manifest",
        "initialize_run",
        "prepare_annotation_scene",
        "plan_annotation",
        "route_annotation_crosscheck",
        "run_marker_tools",
        "generate_evidence",
        "propose_annotations",
        "register_assay_representation",
        "register_knowledge_provider",
        "validate_run",
        "build_review_artifacts",
        "check_review_report",
        "preview_writeback",
        "apply_writeback",
        "annotate",
    }
    assert all(callable(getattr(cell_curator, name)) for name in expected)
    assert cell_curator.__version__


def test_all_benchmark_manifests_are_machine_readable() -> None:
    root = Path(__file__).resolve().parents[1]
    manifests = sorted((root / "benchmarks" / "manifests").glob("*.json"))
    assert {path.stem for path in manifests} == {
        "atac",
        "cite-seq",
        "multiome",
        "rna",
        "spatial",
        "synthetic",
    }
    statuses = {validate_dataset_manifest(path)["status"] for path in manifests}
    assert statuses == {"available", "pending-external-resource"}


def test_benchmark_manifest_validation_is_strict(tmp_path: Path) -> None:
    invalid = tmp_path / "invalid.json"
    invalid.write_text(
        json.dumps(
            {
                "schema_version": 2,
                "dataset": "fixture",
                "modality": "rna",
                "source_url": "local://fixture.tsv",
                "sha256": "0" * 64,
                "license": "CC0-1.0",
                "status": "available",
                "local_path": "fixture.tsv",
                "unexpected": True,
            }
        )
    )
    with pytest.raises(ContractError, match="invalid benchmark manifest"):
        validate_dataset_manifest(invalid)


def test_available_benchmark_manifest_is_portable_as_package_data(tmp_path: Path) -> None:
    root = tmp_path / "share" / "cell-curator"
    fixture = root / "benchmarks" / "fixtures" / "predictions.tsv"
    fixture.parent.mkdir(parents=True)
    fixture.write_text("barcode\ttrue_label\tpredicted_label\na\tT\tT\n")
    manifest = root / "benchmarks" / "manifests" / "synthetic.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 3,
                "dataset": "fixture",
                "modality": "synthetic",
                "source_url": "local://benchmarks/fixtures/predictions.tsv",
                "sha256": sha256(fixture),
                "license": "CC0-1.0",
                "status": "available",
                "local_path": "benchmarks/fixtures/predictions.tsv",
            }
        )
    )
    assert validate_dataset_manifest(manifest)["status"] == "available"


def test_shipped_configuration_schema_matches_typed_model() -> None:
    root = Path(__file__).resolve().parents[1]
    shipped = json.loads((root / "schemas/cell-curator-config-v3.schema.json").read_text())
    assert shipped == configuration_schema()


def test_rapids_lane_requires_complete_gpu_contract(config: dict) -> None:
    config["markers"]["allowed_gpu_contracts"] = []
    with pytest.raises(ValueError, match="fail-closed GPU contract"):
        CellCuratorConfig.model_validate(config)

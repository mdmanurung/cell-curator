from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import REPO_ROOT

from cell_curator.annotation import (
    annotation_contract,
    execute_annotation_workflow,
    run_classifier,
)
from cell_curator.config import (
    AnnotationConfig,
    CellCuratorConfig,
    ClassifierConfig,
    anchor_relative_paths,
)
from cell_curator.contracts import ContractError, sha256


@pytest.fixture
def execution(tmp_path: Path, config: dict, monkeypatch: pytest.MonkeyPatch):
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "synthetic configuration; _setup is tested by the integration suite\n"
    )
    model_path = tmp_path / "model.pkl"
    model_path.write_bytes(b"synthetic model identity")
    executable = tmp_path / "isolated-python"
    executable.write_bytes(b"synthetic interpreter identity; subprocess execution is mocked")
    classifier = ClassifierConfig(
        name="immune",
        method="celltypist",
        model=str(model_path),
        python_executable=str(executable),
        required=True,
    )
    model = CellCuratorConfig.model_validate(config).model_copy(
        update={
            "annotation": AnnotationConfig(
                enabled=True,
                assay="transcriptome",
                representation="logcounts",
                classifiers=[classifier],
            )
        }
    )
    root = tmp_path / "run"
    gate_calls = []
    monkeypatch.setattr("cell_curator.annotation._setup", lambda path: (model, root, {}))
    monkeypatch.setattr(
        "cell_curator.guided.require_label_gate",
        lambda path, **scope: gate_calls.append((Path(path), scope)),
    )
    return SimpleNamespace(
        config_path=config_path,
        config=model,
        root=root,
        executable=executable,
        gate_calls=gate_calls,
    )


def worker_receipt(execution, *, status: str = "complete", digest: str | None = None) -> dict:
    directory = execution.root / "annotation" / "classifiers" / "immune"
    directory.mkdir(parents=True)
    cells = directory / "cell_predictions.tsv"
    cells.write_text("barcode\tpredicted_label\ncell-0\tT cell\n")
    value = {
        "status": status,
        "classifier": "immune",
        "required": True,
        "prediction_sha256": sha256(cells) if digest is None else digest,
    }
    (directory / "status.json").write_text(json.dumps(value))
    return value


def forbid_local_input(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "cell_curator.annotation._input",
        lambda *_: pytest.fail("isolated execution must not fall back locally"),
    )


def test_isolated_classifier_uses_the_pinned_interpreter_and_exact_stage_command(
    execution,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    forbid_local_input(monkeypatch)
    commands = []
    expected = worker_receipt(execution)

    def run(command, **kwargs):
        commands.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, "worker complete", "")

    monkeypatch.setattr(subprocess, "run", run)
    assert run_classifier(execution.config_path, "immune") == expected
    assert commands == [
        (
            [
                str(execution.executable),
                "-m",
                "cell_curator",
                "annotation",
                "classifier",
                "--config",
                str(execution.config_path.resolve()),
                "--name",
                "immune",
            ],
            {"text": True, "capture_output": True, "check": False},
        )
    ]
    assert execution.gate_calls == [(execution.config_path, {"level": "L1", "parent": "ROOT"})]
    contract = annotation_contract(execution.config)
    assert contract["classifier_resources"]["immune"][str(execution.executable)] == sha256(
        execution.executable
    )
    execution.executable.write_bytes(b"changed interpreter")
    assert annotation_contract(execution.config) != contract


@pytest.mark.parametrize(
    "problem",
    ["missing_receipt", "malformed_receipt", "missing_cells", "wrong_digest", "unavailable"],
)
def test_isolated_worker_requires_valid_success_evidence(
    execution,
    monkeypatch: pytest.MonkeyPatch,
    problem: str,
) -> None:
    forbid_local_input(monkeypatch)
    if problem != "missing_receipt":
        worker_receipt(
            execution,
            status="unavailable" if problem == "unavailable" else "complete",
            digest="0" * 64 if problem == "wrong_digest" else None,
        )
    if problem == "missing_cells":
        (execution.root / "annotation/classifiers/immune/cell_predictions.tsv").unlink()
    if problem == "malformed_receipt":
        (execution.root / "annotation/classifiers/immune/status.json").write_text(
            "{broken json"
        )
    monkeypatch.setattr(
        subprocess, "run", lambda command, **_: subprocess.CompletedProcess(command, 0, "", "")
    )
    with pytest.raises(ContractError):
        run_classifier(execution.config_path, "immune")


def test_isolated_worker_failure_preserves_the_error_and_never_falls_back(
    execution,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    forbid_local_input(monkeypatch)
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda command, **_: subprocess.CompletedProcess(
            command, 2, "worker stdout", "model load failed"
        ),
    )
    with pytest.raises(ContractError, match="model load failed"):
        run_classifier(execution.config_path, "immune")


def test_current_interpreter_runs_locally_without_recursive_worker_launch(
    execution,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    execution.config.annotation.classifiers[0].python_executable = sys.executable
    query = ad.AnnData(
        np.asarray([[2.0], [3.0]]),
        obs=pd.DataFrame({"canonical_cluster": ["0", "1"]}, index=["c0", "c1"]),
    )
    query.layers["counts"] = query.X.copy()
    loaded = SimpleNamespace(root=query, assays={"transcriptome": query})
    monkeypatch.setattr("cell_curator.annotation._input", lambda *_: loaded)
    closed = []
    monkeypatch.setattr("cell_curator.data.close_input", lambda value: closed.append(value))
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *_args, **_kwargs: pytest.fail("same interpreter must not recurse"),
    )

    def predictions(adata, **_kwargs):
        np.testing.assert_allclose(adata.X, np.log1p(1e4))
        return pd.DataFrame({"barcode": ["c0", "c1"], "predicted_label": ["T", "B"]}), {
            "status": "complete"
        }

    monkeypatch.setattr("cell_curator.crosschecks.celltypist_hypothesis", predictions)
    receipt = run_classifier(execution.config_path, "immune")
    assert receipt["status"] == "complete"
    assert closed == [loaded]
    assert receipt["prediction_sha256"] == sha256(
        execution.root / "annotation/classifiers/immune/cell_predictions.tsv"
    )


def test_distinct_virtualenv_symlink_launches_the_worker(
    execution,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    execution.executable.unlink()
    execution.executable.symlink_to(sys.executable)
    assert execution.executable.resolve() == Path(sys.executable).resolve()
    forbid_local_input(monkeypatch)
    expected = worker_receipt(execution)
    commands = []

    def run(command, **_kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(subprocess, "run", run)
    assert run_classifier(execution.config_path, "immune") == expected
    assert commands[0][0] == str(execution.executable)


def test_relative_virtualenv_python_keeps_its_lexical_environment_path(
    execution,
    tmp_path: Path,
) -> None:
    alias = tmp_path / "venv/bin/python"
    alias.parent.mkdir(parents=True)
    alias.symlink_to(sys.executable)
    execution.config.annotation.classifiers[0].python_executable = "venv/bin/python"
    anchored = anchor_relative_paths(execution.config, config_dir=tmp_path)
    assert anchored.annotation.classifiers[0].python_executable == str(alias.absolute())
    assert (
        Path(anchored.annotation.classifiers[0].python_executable).resolve()
        == Path(sys.executable).resolve()
    )


def workflow_dependency(monkeypatch: pytest.MonkeyPatch, *, available: bool) -> None:
    original = importlib.util.find_spec
    monkeypatch.setattr(
        importlib.util,
        "find_spec",
        lambda name: (
            object()
            if name == "snakemake" and available
            else None
            if name == "snakemake"
            else original(name)
        ),
    )


def test_workflow_launcher_uses_snakemake_and_checks_its_final_receipt(
    execution,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workflow_dependency(monkeypatch, available=True)
    for name in ("collect_level", "resolve_level", "finalize_annotation"):
        monkeypatch.setattr(
            f"cell_curator.annotation.{name}",
            lambda *_: pytest.fail("launcher must use the DAG"),
        )
    commands = []

    def run(command, **kwargs):
        commands.append((command, kwargs))
        manifest = execution.root / "annotation/annotation.manifest.json"
        manifest.parent.mkdir(parents=True)
        output = manifest.parent / "cell_labels.tsv"
        output.write_text("barcode\tL1\ncell-0\tT\n")
        manifest.write_text(
            json.dumps(
                {
                    "inputs": {str(execution.config_path): sha256(execution.config_path)},
                    "artifacts": {str(output): sha256(output)},
                }
            )
        )
        return subprocess.CompletedProcess(command, 0, "done", "")

    monkeypatch.setattr(subprocess, "run", run)
    manifest = execute_annotation_workflow(execution.config_path, cores=2)
    assert manifest == execution.root / "annotation/annotation.manifest.json"
    assert commands == [
        (
            [
                sys.executable,
                "-m",
                "snakemake",
                "--snakefile",
                str(REPO_ROOT / "workflow/annotation.smk"),
                "--configfile",
                str(execution.config_path.resolve()),
                "--cores",
                "2",
                "--directory",
                str(execution.root),
            ],
            {"text": True, "capture_output": True, "check": False},
        )
    ]


@pytest.mark.parametrize("cores", [0, -1, True, 1.5])
def test_workflow_launcher_rejects_invalid_core_counts(
    execution,
    monkeypatch: pytest.MonkeyPatch,
    cores,
) -> None:
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *_args, **_kwargs: pytest.fail("invalid cores must not launch"),
    )
    with pytest.raises(ContractError, match="positive integer"):
        execute_annotation_workflow(execution.config_path, cores=cores)


def test_workflow_launcher_reports_missing_dependency_without_a_local_bypass(
    execution,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workflow_dependency(monkeypatch, available=False)
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *_args, **_kwargs: pytest.fail("missing workflow dependency must not launch"),
    )
    with pytest.raises(ContractError, match="workflow extra.*Snakemake"):
        execute_annotation_workflow(execution.config_path)


@pytest.mark.parametrize("problem", ["failed_process", "missing_manifest", "stale_artifact"])
def test_workflow_launcher_refuses_incomplete_or_changed_results(
    execution,
    monkeypatch: pytest.MonkeyPatch,
    problem: str,
) -> None:
    workflow_dependency(monkeypatch, available=True)

    def run(command, **_kwargs):
        if problem == "stale_artifact":
            directory = execution.root / "annotation"
            directory.mkdir(parents=True)
            output = directory / "cell_labels.tsv"
            output.write_text("changed result")
            (directory / "annotation.manifest.json").write_text(
                json.dumps({"inputs": {}, "artifacts": {str(output): "0" * 64}})
            )
        return subprocess.CompletedProcess(
            command, 2 if problem == "failed_process" else 0, "", "DAG failed"
        )

    monkeypatch.setattr(subprocess, "run", run)
    with pytest.raises(ContractError):
        execute_annotation_workflow(execution.config_path)


@pytest.mark.parametrize("stage", ["classifier", "workflow"])
def test_execution_entrypoints_check_gate_a_before_launching(
    execution,
    monkeypatch: pytest.MonkeyPatch,
    stage: str,
) -> None:
    def refuse(*_args, **_kwargs):
        raise ContractError("assumptions are not approved")

    monkeypatch.setattr("cell_curator.guided.require_label_gate", refuse)
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *_args, **_kwargs: pytest.fail("unapproved execution must not launch"),
    )
    with pytest.raises(ContractError, match="not approved"):
        if stage == "classifier":
            run_classifier(execution.config_path, "immune")
        else:
            execute_annotation_workflow(execution.config_path)

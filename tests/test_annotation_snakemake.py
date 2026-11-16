"""The annotation DAG must schedule every level and track source/model changes."""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest
import yaml
from conftest import REPO_ROOT

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="Snakemake is an optional workflow extra",
)


def _fixture(tmp_path: Path, config: dict) -> tuple[Path, Path, Path]:
    # Dry-run checks the dependency graph; these are declared artifacts, not data.
    source = tmp_path / "source.h5ad"
    corpus = tmp_path / "corpus.json"
    model = tmp_path / "celltypist.pkl"
    consensus = tmp_path / "consensus"
    (consensus / "work").mkdir(parents=True)
    for path in (
        source,
        corpus,
        model,
        *(
            consensus / name
            for name in ("consensus.tsv", "run_manifest.json", "work/cells.tsv")
        ),
    ):
        path.write_text("{}")
    config["input"]["path"] = str(source)
    config["knowledge"]["providers"] = [
        {
            "kind": "corpus",
            "path": str(corpus),
            "version": "fixture-v1",
            "sha256": "0" * 64,
        }
    ]
    config["annotation"] = {
        "enabled": True,
        "assay": "transcriptome",
        "representation": "logcounts",
        "classifiers": [{"name": "matched", "method": "celltypist", "model": str(model)}],
        "consensus_directory": str(consensus),
    }
    config["guidance"]["hierarchy_levels"] = ["L1"]
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    root = Path(config["run"]["output_root"]) / config["run"]["run_id"]
    for name in [
        "context.json",
        "plan.json",
        "assumptions.json",
        "input_contract.manifest.json",
        "canonical/input.h5ad",
    ]:
        output = root / name
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("{}")
    return path, root, model


def _run(path: Path, cwd: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    env["XDG_CACHE_HOME"] = str(cwd / ".cache")
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            str(REPO_ROOT / "workflow/annotation.smk"),
            "--configfile",
            str(path),
            "--cores",
            "1",
            *extra,
        ],
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        timeout=120,
    )


def test_annotation_dag_schedules_classifiers_and_the_first_level(tmp_path, config):
    path, root, _ = _fixture(tmp_path, config)
    result = _run(path, tmp_path, "annotation_review", "--dry-run")
    assert result.returncode == 0, result.stdout + result.stderr
    output = result.stdout + result.stderr
    assert "annotation_classifier" in output
    assert "annotation_finalize" in output
    assert "annotation_build_review" in output
    for level in ["L1"]:
        assert str(root / "annotation" / "levels" / level / "decision.manifest.json") in output
        assert str(root / "annotation" / "levels" / level / "marker_evidence.svg") in output
    assert not (root / "annotation").exists(), "dry-run must not execute annotation"


def test_annotation_dag_tracks_model_staleness_without_approving_or_writing_back(
    tmp_path, config
):
    path, root, model = _fixture(tmp_path, config)
    outputs = [
        "annotation/classifiers/matched/status.json",
        "annotation/classifiers/matched/cell_predictions.tsv",
        "annotation/annotation.manifest.json",
        "annotation/cell_labels.tsv",
        "annotation/report.html",
        "hierarchy/execution.manifest.json",
        "mapping/parent_child_mapping.tsv",
        "mapping/mapping.validation.json",
    ]
    for level in ["L1"]:
        outputs.extend(
            f"annotation/levels/{level}/{name}"
            for name in [
                "evidence.manifest.json",
                "markers.tsv",
                "scores.tsv",
                "states.tsv",
                "marker_evidence.svg",
                "labels.tsv",
                "decision.manifest.json",
            ]
        )
    for name in outputs:
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("{}")
    baseline = time.time() - 100
    for item in tmp_path.rglob("*"):
        if item.is_file():
            os.utime(item, (baseline, baseline))
    for name in outputs:
        step = 1 if "/classifiers/" in name else 5
        for index, level in enumerate(["L1"], 2):
            if f"/levels/{level}/" in name:
                step = index
        os.utime(root / name, (baseline + step, baseline + step))
    result = _run(path, tmp_path, "--dry-run", "--rerun-triggers", "mtime")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Nothing to be done" in result.stdout + result.stderr
    future = baseline + 10
    os.utime(model, (future, future))
    result = _run(path, tmp_path, "--dry-run", "--rerun-triggers", "mtime")
    assert result.returncode == 0, result.stdout + result.stderr
    output = result.stdout + result.stderr
    assert "annotation_classifier" in output and "annotation_finalize" in output
    assert "annotation_level_evidence" in output and "annotation_level_decision" in output
    assert not (root / "gates").exists()
    assert not (root / "writeback").exists()


def test_annotation_dag_rejects_config_overlays(tmp_path, config):
    path, _, _ = _fixture(tmp_path, config)
    result = _run(path, tmp_path, "--dry-run", "--config", "annotation={enabled: false}")
    assert result.returncode != 0
    assert "overlays are prohibited" in result.stdout + result.stderr

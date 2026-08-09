from __future__ import annotations

import json
import os
import shutil
import site
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml


def _run(command: list[str], *, cwd: Path, env: dict[str, str] | None = None) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        check=False,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    return completed.stdout


def _distributable_only(directory: str, entries: list[str]) -> set[str]:
    """Copy only what a fresh clone would carry into a build.

    Naming transient directories one at a time is unsafe: every new dot-directory
    (a second virtualenv, a tool cache, a scratch area) silently joins the copy
    and can make this test fail for reasons unrelated to packaging, or race a
    live writer mid-copy. Nothing under a dot-directory is a build input, so
    exclude the whole class instead of enumerating members.
    """

    del directory
    excluded = {"build", "dist", "__pycache__"}
    return {
        entry
        for entry in entries
        if entry.startswith(".") or entry in excluded or entry.endswith(".egg-info")
    }


def test_clean_wheel_install_import_cli_and_package_data(tmp_path: Path) -> None:
    repository = Path(__file__).resolve().parents[1]
    source_tree = tmp_path / "source-tree"
    shutil.copytree(repository, source_tree, ignore=_distributable_only)
    distribution = tmp_path / "dist"
    distribution.mkdir()
    _run(
        [
            sys.executable,
            "-m",
            "build",
            "--wheel",
            "--no-isolation",
            "--outdir",
            str(distribution),
            str(source_tree),
        ],
        cwd=tmp_path,
    )
    wheels = list(distribution.glob("cell_curator-*.whl"))
    assert len(wheels) == 1

    environment = tmp_path / "clean-venv"
    uv = shutil.which("uv")
    assert uv is not None, "the locked-environment test requires the project's uv runner"
    clean_env = os.environ.copy()
    clean_env.pop("PYTHONPATH", None)
    clean_env["UV_CACHE_DIR"] = str(tmp_path / "uv-cache")
    _run(
        [
            uv,
            "venv",
            "--python",
            sys.executable,
            str(environment),
        ],
        cwd=tmp_path,
        env=clean_env,
    )
    python = environment / "bin" / "python"
    cli = environment / "bin" / "cell-curator"
    clean_env["PATH"] = f"{environment / 'bin'}{os.pathsep}{clean_env.get('PATH', '')}"
    _run(
        [uv, "pip", "install", "--python", str(python), "--no-deps", str(wheels[0])],
        cwd=tmp_path,
        env=clean_env,
    )
    new_site = Path(
        _run(
            [str(python), "-c", "import site; print(site.getsitepackages()[0])"],
            cwd=tmp_path,
            env=clean_env,
        ).strip()
    )
    (new_site / "locked-dependencies.pth").write_text(site.getsitepackages()[0] + "\n")
    probe = _run(
        [
            str(python),
            "-c",
            (
                "import importlib.metadata as m,json,cell_curator;"
                "d=m.distribution('cell-curator');f=list(d.files or []);"
                "s=lambda suffix:next(str(d.locate_file(x)) for x in f if str(x).endswith(suffix));"
                "print(json.dumps({'version':cell_curator.__version__,"
                "'module':cell_curator.__file__,"
                "'files':[str(x) for x in f],"
                "'snakefile':s('share/cell-curator/workflow/Snakefile'),"
                "'example':s('share/cell-curator/workflow/config/preauthorized.example.yaml'),"
                "'asset':s('share/cell-curator/assets/config.template.yaml')}))"
            ),
        ],
        cwd=tmp_path,
        env=clean_env,
    )
    installed = json.loads(probe)
    assert installed["version"] == "0.2.0"
    assert str(environment) in installed["module"]
    files = installed["files"]
    expected_families = (
        "share/cell-curator/SKILL.md",
        "share/cell-curator/agents/",
        "share/cell-curator/assets/",
        "share/cell-curator/references/",
        "share/cell-curator/workflow/Snakefile",
        "share/cell-curator/workflow/config/",
        "share/cell-curator/workflow/profiles/slurm/",
        "share/cell-curator/workflow/scripts/",
        "share/cell-curator/benchmarks/",
        "share/cell-curator/benchmarks/manifests/",
        "share/cell-curator/benchmarks/fixtures/",
        "share/cell-curator/benchmarks/results/",
        "share/cell-curator/schemas/",
        "cell_curator/py.typed",
    )
    assert all(any(family in item for item in files) for family in expected_families)
    assert shutil.which("cell-curator", path=clean_env["PATH"]) == str(cli)
    assert "usage: cell-curator" in _run([str(cli), "--help"], cwd=tmp_path, env=clean_env)

    counts = np.array([[5, 0], [4, 0], [0, 5], [0, 4]], dtype=np.float32)
    fixture = ad.AnnData(
        X=np.log1p(counts),
        obs=pd.DataFrame(
            {
                "canonical_cluster": ["0", "0", "1", "1"],
                "lineage": ["immune"] * 4,
            },
            index=[f"wheel-{index}" for index in range(4)],
        ),
        var=pd.DataFrame({"gene_symbol": ["A", "B"]}, index=["A", "B"]),
    )
    fixture.layers["counts"] = counts
    fixture.layers["logcounts"] = np.log1p(counts)
    fixture.obsm["X_pca_harmony"] = np.array(
        [[-2, 0], [-1, 0], [1, 0], [2, 0]], dtype=np.float32
    )
    input_path = tmp_path / "wheel-fixture.h5ad"
    fixture.write_h5ad(input_path)
    markers = tmp_path / "wheel-markers.json"
    markers.write_text(json.dumps({"A type": {"positive": ["A"]}}))
    workflow_config = yaml.safe_load(Path(installed["asset"]).read_text())
    workflow_config["run"].update(
        {
            "run_id": "wheel-zero-candidate-run",
            "output_root": str(tmp_path / "wheel-results"),
            "mode": "evidence-only",
        }
    )
    workflow_config["input"]["path"] = str(input_path)
    workflow_config["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    workflow_config["adaptive"]["eligibility"].update(
        {"min_cells": 2, "min_donors": None, "min_captures": None}
    )
    workflow_config["adaptive"]["parent_local"]["enabled"] = False
    workflow_config["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": str(markers), "version": "wheel-fixture"}
    ]
    workflow_config["guidance"] = {
        "interactive": False,
        "preauthorized": True,
        "reviewer": "wheel-fixture-reviewer",
        "pause_after_each_level": False,
        "assumptions": ["Wheel fixture is pre-authorized."],
        "hierarchy_levels": ["L1"],
        "context": {
            "organism": "human",
            "tissue": "synthetic",
            "developmental_stage": "adult",
            "condition": "healthy",
            "technology": "dissociated",
            "experimental_context": "installed-wheel zero-candidate execution",
            "composition_covariates": ["lineage"],
            "vocabulary_contract": "local fixture",
            "reference_contract": "local_only",
            "visualization_key": "X_pca_harmony",
            "writeback_prefix": "wheel_label",
        },
    }
    config_path = tmp_path / "wheel-config.yaml"
    config_path.write_text(yaml.safe_dump(workflow_config, sort_keys=False))
    _run(
        [
            str(python),
            "-m",
            "snakemake",
            "--snakefile",
            installed["snakefile"],
            "all_evidence",
            "--configfile",
            str(config_path),
            "--cores",
            "1",
        ],
        cwd=tmp_path,
        env=clean_env,
    )
    workflow_root = tmp_path / "wheel-results" / "wheel-zero-candidate-run"
    assert (workflow_root / "comparators" / "groups").is_dir()
    comparator_manifest = json.loads(
        (workflow_root / "comparators" / "comparator_manifest.json").read_text()
    )
    assert comparator_manifest["zero_applicable_lane_run"] is True
    assert (workflow_root / "evidence" / "summary" / "evidence.manifest.json").is_file()

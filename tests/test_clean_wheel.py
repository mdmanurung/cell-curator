from __future__ import annotations

import json
import os
import shutil
import site
import subprocess
import sys
from pathlib import Path


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


def test_clean_wheel_install_import_cli_and_package_data(tmp_path: Path) -> None:
    repository = Path(__file__).resolve().parents[1]
    source_tree = tmp_path / "source-tree"
    shutil.copytree(
        repository,
        source_tree,
        ignore=shutil.ignore_patterns(
            ".git", ".venv", "build", "*.egg-info", "__pycache__", ".pytest_cache"
        ),
    )
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
                "f=m.files('cell-curator') or [];"
                "print(json.dumps({'version':cell_curator.__version__,"
                "'module':cell_curator.__file__,"
                "'skill':any(str(x).endswith('share/cell-curator/SKILL.md') for x in f),"
                "'schema':any('share/cell-curator/schemas/' in str(x) for x in f),"
                "'typed':any(str(x).endswith('cell_curator/py.typed') for x in f)}))"
            ),
        ],
        cwd=tmp_path,
        env=clean_env,
    )
    installed = json.loads(probe)
    assert installed == {
        "version": "0.2.0",
        "module": installed["module"],
        "skill": True,
        "schema": True,
        "typed": True,
    }
    assert str(environment) in installed["module"]
    assert "usage: cell-curator" in _run([str(cli), "--help"], cwd=tmp_path, env=clean_env)

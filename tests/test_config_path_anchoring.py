"""Relative configuration paths must not depend on the process working directory.

A notebook's cwd is whatever the kernel was started in, and it moves: an
`os.chdir` in one cell, a kernel restart from a different launcher, a script
invoked from elsewhere. Resolving `run.output_root` and `input.path` against the
cwd at the moment each artifact is opened made a run silently depend on that,
and the failure surfaced as "input object does not exist" rather than as a path
problem.

Paths are therefore anchored to the configuration file's own directory. These
tests pin that behaviour, the hash-safety property that makes it safe, and the
one case the ambiguity guard cannot detect.
"""

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


def _write_relative_fixture(home: Path) -> Path:
    """Write a config whose input, output, and provider paths are all relative."""

    rng = np.random.default_rng(11)
    counts = np.zeros((40, 6), dtype=np.float32)
    counts[:20, :2] = rng.poisson(9, size=(20, 2)) + 1
    counts[20:, 2:4] = rng.poisson(9, size=(20, 2)) + 1
    genes = ["A1", "A2", "B1", "B2", "N1", "N2"]
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=pd.DataFrame(
            {
                "canonical_cluster": ["0"] * 20 + ["1"] * 20,
                "lineage": pd.Categorical(["immune"] * 40),
            },
            index=[f"cell-{index:03d}" for index in range(40)],
        ),
        var=pd.DataFrame({"gene_symbol": genes}, index=genes),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = np.vstack(
        [rng.normal(-3, 0.2, size=(20, 3)), rng.normal(3, 0.2, size=(20, 3))]
    ).astype(np.float32)
    adata.write_h5ad(home / "source.h5ad")

    (home / "markers.json").write_text(
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
            "run_id": "anchoring-v1",
            "output_root": "results",  # relative on purpose
            "mode": "evidence-only",
        }
    )
    value["input"]["path"] = "source.h5ad"  # relative on purpose
    value["input"]["keys"]["scope"] = "lineage"
    value["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    value["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": "markers.json", "version": "anchoring-v1"}
    ]
    config_path = home / "config.yaml"
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))
    return config_path


def test_relative_paths_anchor_to_the_config_not_the_working_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    home = tmp_path / "project"
    home.mkdir()
    config_path = _write_relative_fixture(home)

    monkeypatch.chdir(home)
    first = load_config(config_path)

    # The notebook failure mode: cwd moves between calls.
    elsewhere = tmp_path / "somewhere-else"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    second = load_config(config_path)

    assert first.run.output_root == second.run.output_root
    assert first.input.path == second.input.path
    assert Path(second.input.path) == (home / "source.h5ad").resolve()
    assert Path(second.run.output_root) == (home / "results").resolve()
    assert Path(second.knowledge.providers[0].path) == (home / "markers.json").resolve()
    assert Path(second.input.path).is_file(), "the anchored input must actually resolve"


def test_a_real_phase_still_runs_after_the_working_directory_moves(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End to end, not just string comparison: a phase must execute post-chdir."""

    from cell_curator.pipeline import initialize, run_root

    home = tmp_path / "project"
    home.mkdir()
    config_path = _write_relative_fixture(home)

    monkeypatch.chdir(home)
    initialize(config_path)

    elsewhere = tmp_path / "somewhere-else"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    from cell_curator.api import canonicalize_input

    canonicalize_input(config_path)
    root = run_root(load_config(config_path))
    assert root == (home / "results" / "anchoring-v1").resolve()
    assert (root / "prepared" / "cells.tsv").is_file()


def test_anchoring_never_changes_the_recorded_configuration_hash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Provenance hashes the file bytes, so in-memory anchoring must not shift them."""

    home = tmp_path / "project"
    home.mkdir()
    config_path = _write_relative_fixture(home)
    before = sha256(config_path)

    monkeypatch.chdir(home)
    load_config(config_path)
    monkeypatch.chdir(tmp_path)
    load_config(config_path)

    assert sha256(config_path) == before


def test_an_ambiguous_legacy_run_is_refused_rather_than_silently_restarted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    home = tmp_path / "project"
    home.mkdir()
    config_path = _write_relative_fixture(home)

    # A run recorded where the old cwd-relative behaviour would have put it.
    legacy = tmp_path / "launched-from-here"
    (legacy / "results" / "anchoring-v1").mkdir(parents=True)
    (legacy / "results" / "anchoring-v1" / "run_state.json").write_text("{}")

    monkeypatch.chdir(legacy)
    with pytest.raises(ContractError, match="relative and ambiguous"):
        load_config(config_path)


def test_the_ambiguity_guard_is_best_effort_and_stays_silent_once_cwd_moves_on(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Document the limit: the guard cannot find a legacy run it cannot see.

    It compares the anchored location against what *the current* cwd would have
    produced. Once the caller is no longer in the directory the legacy run was
    created from, there is nothing to compare against and a fresh run starts. A
    reader must not mistake this guard for an exhaustive search.
    """

    home = tmp_path / "project"
    home.mkdir()
    config_path = _write_relative_fixture(home)

    legacy = tmp_path / "launched-from-here"
    (legacy / "results" / "anchoring-v1").mkdir(parents=True)
    (legacy / "results" / "anchoring-v1" / "run_state.json").write_text("{}")

    unrelated = tmp_path / "third-directory"
    unrelated.mkdir()
    monkeypatch.chdir(unrelated)

    config = load_config(config_path)  # no raise
    assert Path(config.run.output_root) == (home / "results").resolve()


def test_absolute_paths_are_left_exactly_as_written(tmp_path: Path) -> None:
    home = tmp_path / "project"
    home.mkdir()
    config_path = _write_relative_fixture(home)
    value = yaml.safe_load(config_path.read_text())
    absolute_root = (tmp_path / "explicit-results").resolve()
    value["run"]["output_root"] = str(absolute_root)
    value["input"]["path"] = str((home / "source.h5ad").resolve())
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))

    config = load_config(config_path)
    assert Path(config.run.output_root) == absolute_root


def test_the_session_object_resolves_its_output_root_immediately(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CellCurator caches config_path, so it cannot wait for load_config to anchor."""

    from cell_curator import CellCurator

    home = tmp_path / "project"
    home.mkdir()
    monkeypatch.chdir(home)

    rng = np.random.default_rng(3)
    counts = rng.poisson(5, size=(20, 4)).astype(np.float32)
    genes = ["A1", "A2", "B1", "B2"]
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=pd.DataFrame(
            {"leiden": ["0"] * 10 + ["1"] * 10},
            index=[f"c-{index:02d}" for index in range(20)],
        ),
        var=pd.DataFrame({"gene_symbol": genes}, index=genes),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = rng.normal(size=(20, 3)).astype(np.float32)

    curator = CellCurator(
        adata=adata,
        organism="human",
        tissue="synthetic tissue",
        run_id="session-anchor-v1",
        cluster_key="leiden",
        output_root="relative-results",
    )
    assert curator.output_root.is_absolute()
    assert curator.config_path.is_absolute()

    monkeypatch.chdir(tmp_path)
    assert curator.config_path.is_file(), "the cached config path must survive a chdir"
    assert curator.validate().run.run_id == "session-anchor-v1"

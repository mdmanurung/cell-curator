from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from test_annotation_pipeline import annotation_fixture

from cell_curator.annotation import _ranked_markers, _ranking_positions, collect_level
from cell_curator.config import MarkerConfig

CONTEXT = {"input_sha256": "abc", "detection": None}


def _data(sizes: dict[str, int]) -> tuple[np.ndarray, pd.Series, pd.Index]:
    rng = np.random.default_rng(1)
    names = [f"c{i:04d}" for i in range(sum(sizes.values()))]
    groups = pd.Series(np.repeat(list(sizes), list(sizes.values())), index=names)
    matrix = rng.poisson(2.0, (len(groups), 6)).astype(float)
    matrix[groups.to_numpy() == "a", 0] += 3
    return matrix, groups, pd.Index([f"g{i}" for i in range(6)])


def _rank(matrix, groups, features, tmp_path, *, cap=None, cache=True):
    config = SimpleNamespace(markers=MarkerConfig(
        ranking_max_cells_per_cluster=cap, ranking_cache_dir=str(tmp_path / "cache") if cache else ""))
    return _ranked_markers(matrix, groups, features, matrix, config, CONTEXT, "parent_remainder")


def test_cap_none_keeps_every_cell_and_matches_the_plain_computation(tmp_path: Path) -> None:
    from cell_curator.evidence import compute_bottom_up_markers

    matrix, groups, features = _data({"a": 30, "b": 40})
    assert _ranking_positions(groups, None, 0).tolist() == list(range(70))
    plain = compute_bottom_up_markers(matrix, groups, features, top_n=None, methods=("wilcoxon", "logreg"),
                                      detection=matrix, marker_config=MarkerConfig().model_dump(mode="json"))
    got = _rank(matrix, groups, features, tmp_path, cache=False)
    pd.testing.assert_frame_equal(got, plain)
    assert got.attrs["runtime"]["ranking_cache"] == "off"


def test_cap_subsamples_reproducibly_per_cluster() -> None:
    _, groups, _ = _data({"a": 30, "b": 8})
    first = _ranking_positions(groups, 10, 0)
    assert first.tolist() == _ranking_positions(groups, 10, 0).tolist()
    assert groups.iloc[first].value_counts().to_dict() == {"a": 10, "b": 8}
    assert first.tolist() != _ranking_positions(groups, 10, 1).tolist()
    # A cluster's draw ignores the other clusters and the row order.
    shuffled = groups.iloc[::-1]
    assert set(shuffled.iloc[_ranking_positions(shuffled, 10, 0)].index) == set(groups.iloc[first].index)
    only_a = groups[groups.eq("a")]
    assert set(only_a.iloc[_ranking_positions(only_a, 10, 0)].index) == set(groups.iloc[first].index) - set(groups[groups.eq("b")].index)


def test_cache_miss_then_hit_gives_identical_frames_and_key_changes_miss(tmp_path: Path) -> None:
    matrix, groups, features = _data({"a": 30, "b": 40})
    miss = _rank(matrix, groups, features, tmp_path, cap=15)
    hit = _rank(matrix, groups, features, tmp_path, cap=15)
    assert (miss.attrs["runtime"]["ranking_cache"], hit.attrs["runtime"]["ranking_cache"]) == ("miss", "hit")
    pd.testing.assert_frame_equal(miss, hit, check_exact=True)
    assert hit.attrs["runtime"]["ranked_cells_per_cluster"] == {"a": 15, "b": 15}
    assert _rank(matrix, groups, features, tmp_path, cap=20).attrs["runtime"]["ranking_cache"] == "miss"
    moved = groups.copy()
    moved.iloc[0] = "b"
    assert _rank(matrix, moved, features, tmp_path, cap=15).attrs["runtime"]["ranking_cache"] == "miss"
    assert len(list((tmp_path / "cache").glob("*.pkl"))) == 3
    assert not list((tmp_path / "cache").glob(".*partial*"))


def test_pipeline_records_cap_seed_counts_and_cache_status(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config_path, root, _ = annotation_fixture(tmp_path, monkeypatch, markers={
        "ranking_max_cells_per_cluster": 5, "ranking_seed": 3, "ranking_cache_dir": str(tmp_path / "rc")})
    manifest = collect_level(config_path, "L1")
    runtime = manifest["marker_runtime"]
    assert manifest["ranking"]["max_cells_per_cluster"] == 5
    assert all(max(item["ranked_cells_per_cluster"].values()) <= 5 for item in runtime)
    assert "miss" in {item["ranking_cache"] for item in runtime}  # a mirrored sibling pair can hit within a run
    scores = pd.read_csv(root / "annotation" / "levels" / "L1" / "scores.tsv", sep="\t")
    assert len(scores) > 0

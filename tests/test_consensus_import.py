from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest
from consensus_fixture import write_consensus_run
from test_annotation_pipeline import annotation_fixture

from cell_curator.annotation import _check_receipt, collect_level
from cell_curator.consensus_import import import_consensus_level
from cell_curator.contracts import ContractError, sha256
from cell_curator.hierarchy import LABEL_COLUMNS

UNKNOWN_CLUSTER = "3"


@pytest.fixture
def run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    config_path, root, expected = annotation_fixture(tmp_path, monkeypatch)
    collect_level(config_path, "L1")
    consensus = write_consensus_run(
        config_path,
        root,
        {cluster: values[0] for cluster, values in expected.items()},
        unknown={UNKNOWN_CLUSTER},
    )
    return config_path, root, consensus


def test_consensus_labels_become_the_level_decisions(run) -> None:
    config_path, root, consensus = run
    receipt = import_consensus_level(config_path, "L1")
    target = root / "annotation" / "levels" / "L1"
    labels = pd.read_csv(target / "labels.tsv", sep="\t", keep_default_na=False, dtype=str)
    assert list(labels.columns) == list(LABEL_COLUMNS)
    assert set(labels["parent_scope"]) == {"ROOT"} and set(labels["level"]) == {"L1"}
    unknown = labels[labels["cluster_id"] == UNKNOWN_CLUSTER].iloc[0]
    assert (unknown["label"], unknown["cl_id"]) == ("Unknown", "")
    assert unknown["rationale"].startswith("consensus abstain: ")
    assert set(labels["confidence"].astype(float)) == {0.4, 0.9}
    assert receipt["decided_by"] == "annotation-consensus"
    assert _check_receipt(target / "decision.manifest.json") == receipt
    assert sha256(target / "consensus" / "consensus.tsv") == sha256(consensus / "consensus.tsv")


def test_replaying_the_import_rewrites_nothing(run) -> None:
    config_path, root, consensus = run
    first = import_consensus_level(config_path, "L1")
    labels = root / "annotation" / "levels" / "L1" / "labels.tsv"
    digest = sha256(labels)
    assert import_consensus_level(config_path, "L1") == first
    assert sha256(labels) == digest


def test_only_the_first_level_can_be_imported(run) -> None:
    config_path, _, consensus = run
    with pytest.raises(ContractError, match="first level 'L1' only"):
        import_consensus_level(config_path, "L2")


def _drop_a_cell(out: Path) -> None:
    cells = pd.read_csv(out / "work" / "cells.tsv", sep="\t", dtype=str)
    cells.iloc[1:].to_csv(out / "work" / "cells.tsv", sep="\t", index=False)


def _move_a_cell(out: Path) -> None:
    cells = pd.read_csv(out / "work" / "cells.tsv", sep="\t", dtype=str)
    cells.loc[0, "cluster_id"] = "7" if cells.loc[0, "cluster_id"] != "7" else "6"
    cells.to_csv(out / "work" / "cells.tsv", sep="\t", index=False)


def _edit_consensus(**changes: str) -> Callable[[Path], None]:
    def edit(out: Path) -> None:
        table = pd.read_csv(out / "consensus.tsv", sep="\t", keep_default_na=False, dtype=str)
        for column, value in changes.items():
            table.loc[0, column] = value
        table.to_csv(out / "consensus.tsv", sep="\t", index=False)

    return edit


def _drop_a_cluster(out: Path) -> None:
    table = pd.read_csv(out / "consensus.tsv", sep="\t", keep_default_na=False, dtype=str)
    table.iloc[1:].to_csv(out / "consensus.tsv", sep="\t", index=False)


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        (_drop_a_cell, "different cells"),
        (_move_a_cell, "different clusters than the input"),
        (_drop_a_cluster, "differ from the clusters"),
        (_edit_consensus(label="not a vocabulary label"), "outside the approved L1 vocabulary"),
        (_edit_consensus(cl_id="CL:9999999"), "cl_id differs from the vocabulary"),
        (_edit_consensus(cluster_cell_share="1.5"), "between 0 and 1"),
    ],
    ids=[
        "missing-cell",
        "moved-cell",
        "missing-cluster",
        "bad-label",
        "bad-cl-id",
        "bad-share",
    ],
)
def test_a_consensus_run_that_does_not_fit_the_run_is_refused(run, edit, message) -> None:
    config_path, root, consensus = run
    edit(consensus)
    with pytest.raises(ContractError, match=message):
        import_consensus_level(config_path, "L1")
    assert not (root / "annotation" / "levels" / "L1" / "decision.manifest.json").exists()


def test_a_config_without_a_consensus_directory_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config_path, _, _ = annotation_fixture(tmp_path, monkeypatch, consensus=False)
    collect_level(config_path, "L1")
    with pytest.raises(ContractError, match="names no annotation-consensus run"):
        import_consensus_level(config_path, "L1")


def test_a_corpus_config_with_a_second_level_is_refused_at_load(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import yaml

    from cell_curator.config import load_config

    config_path, _, _ = annotation_fixture(tmp_path, monkeypatch)
    value = yaml.safe_load(config_path.read_text())
    value["guidance"]["hierarchy_levels"] = ["L1", "L2"]
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))
    with pytest.raises(Exception, match="L1 only"):
        load_config(config_path)

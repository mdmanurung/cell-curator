from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest
from test_annotation_pipeline import annotation_fixture

from cell_curator.annotation import _check_receipt, _input, collect_level
from cell_curator.config import load_config
from cell_curator.consensus_import import import_consensus_level
from cell_curator.contracts import ContractError, sha256
from cell_curator.data import close_input
from cell_curator.hierarchy import LABEL_COLUMNS

CONSENSUS_COLUMNS = [
    "cluster_id", "label", "cl_id", "decision_rule", "rationale",
    "rejected_alternatives", "cluster_cell_share", "cluster_cell_lead",
]  # fmt: skip
UNKNOWN_CLUSTER = "3"


def consensus_run(
    tmp_path: Path, config_path: Path, root: Path, expected: dict[str, list[str]]
) -> Path:
    """A stand-in for an annotation-consensus output: its tables, not its computation."""
    plan = json.loads((root / "plan.json").read_text())
    vocabulary = pd.read_csv(plan["vocabulary_path"], sep="\t", keep_default_na=False)
    cl_id = dict(zip(vocabulary["label"], vocabulary["cl_id"], strict=True))
    rows = [
        {
            "cluster_id": cluster,
            "label": "Unknown" if cluster == UNKNOWN_CLUSTER else values[0],
            "cl_id": "" if cluster == UNKNOWN_CLUSTER else cl_id[values[0]],
            "decision_rule": "abstain" if cluster == UNKNOWN_CLUSTER else "all_agree",
            "rationale": f"synthetic call for cluster {cluster}",
            "rejected_alternatives": "[]",
            "cluster_cell_share": 0.4 if cluster == UNKNOWN_CLUSTER else 0.9,
            "cluster_cell_lead": 0.1 if cluster == UNKNOWN_CLUSTER else 0.5,
        }
        for cluster, values in sorted(expected.items())
    ]
    out = tmp_path / "consensus"
    (out / "work").mkdir(parents=True)
    pd.DataFrame(rows, columns=CONSENSUS_COLUMNS).to_csv(
        out / "consensus.tsv", sep="\t", index=False
    )
    loaded = _input(load_config(config_path), root)
    try:
        obs = loaded.root.obs
        cells = pd.DataFrame(
            {
                "barcode": loaded.root.obs_names.astype(str),
                "cluster_id": obs[load_config(config_path).input.keys.canonical_parent]
                .astype(str)
                .to_numpy(),
            }
        )
    finally:
        close_input(loaded)
    cells.to_csv(out / "work" / "cells.tsv", sep="\t", index=False)
    (out / "run_manifest.json").write_text(json.dumps({"note": "synthetic"}) + "\n")
    return out


@pytest.fixture
def run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    config_path, root, expected = annotation_fixture(tmp_path, monkeypatch)
    collect_level(config_path, "L1")
    return config_path, root, consensus_run(tmp_path, config_path, root, expected)


def test_consensus_labels_become_the_level_decisions(run) -> None:
    config_path, root, consensus = run
    receipt = import_consensus_level(config_path, "L1", consensus)
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
    first = import_consensus_level(config_path, "L1", consensus)
    labels = root / "annotation" / "levels" / "L1" / "labels.tsv"
    digest = sha256(labels)
    assert import_consensus_level(config_path, "L1", consensus) == first
    assert sha256(labels) == digest


def test_only_the_first_level_can_be_imported(run) -> None:
    config_path, _, consensus = run
    with pytest.raises(ContractError, match="first level 'L1' only"):
        import_consensus_level(config_path, "L2", consensus)


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
        import_consensus_level(config_path, "L1", consensus)
    assert not (root / "annotation" / "levels" / "L1" / "decision.manifest.json").exists()

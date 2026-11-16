"""Import the labels of an ``annotation-consensus`` run as the decisions of one hierarchy level.

The consensus skill (``skills/annotation-consensus``) chooses the label of each cluster. This
module only checks that run against the cell-curator run and writes the same ``labels.tsv`` and
``decision.manifest.json`` that write-back, export and review already read. It assigns nothing.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import pandas as pd

from .annotation import _check_receipt, _frame, _input, _receipt, _setup
from .contracts import ContractError, sha256
from .hierarchy import LABEL_COLUMNS

CONSENSUS_COLUMNS = (
    "cluster_id",
    "label",
    "cl_id",
    "decision_rule",
    "rationale",
    "rejected_alternatives",
    "cluster_cell_share",
    "cluster_cell_lead",
)
UNKNOWN = "Unknown"
ROOT = "ROOT"


def _read(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise ContractError(f"consensus run lacks {path.name}: {path}")
    return pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str)


def _share(table: pd.DataFrame, column: str) -> pd.Series:
    values = pd.to_numeric(table[column], errors="coerce")
    if values.isna().any() or ((values < 0) | (values > 1)).any():
        raise ContractError(f"consensus column {column} must hold numbers between 0 and 1")
    return values


def _check_cells(consensus_dir: Path, config: Any, root: Path) -> None:
    """The consensus run must describe the same cells and clusters as this run's input."""

    from .data import close_input

    theirs = _read(consensus_dir / "work" / "cells.tsv")
    loaded = _input(config, root)
    try:
        ours = pd.DataFrame(
            {
                "barcode": loaded.root.obs_names.astype(str),
                "cluster_id": loaded.root.obs[config.input.keys.canonical_parent]
                .astype(str)
                .to_numpy(),
            }
        )
    finally:
        close_input(loaded)
    if theirs["barcode"].duplicated().any() or ours["barcode"].duplicated().any():
        raise ContractError(
            "cells must be unique by barcode in the consensus run and the input"
        )
    joined = ours.merge(
        theirs,
        on="barcode",
        how="outer",
        suffixes=("", "_consensus"),
        indicator=True,
        validate="1:1",
    )
    if not joined["_merge"].eq("both").all():
        raise ContractError(
            "the consensus run and the input hold different cells: "
            f"{int(joined['_merge'].eq('left_only').sum())} only in the input, "
            f"{int(joined['_merge'].eq('right_only').sum())} only in the consensus run"
        )
    if (joined["cluster_id"] != joined["cluster_id_consensus"]).any():
        raise ContractError(
            "the consensus run assigns cells to different clusters than the input"
        )


def import_consensus_level(
    config_path: str | Path, level: str, consensus_dir: str | Path
) -> dict[str, Any]:
    """Write ``levels/<level>/labels.tsv`` and its decision receipt from a consensus run."""

    from .guided import require_label_gate

    config, root, plan = _setup(config_path)
    levels = list(config.guidance.hierarchy_levels)
    if level != levels[0]:
        raise ContractError(
            f"consensus import covers the first level {levels[0]!r} only, not {level!r}"
        )
    consensus_dir = Path(consensus_dir)
    target = root / "annotation" / "levels" / level
    _check_receipt(target / "evidence.manifest.json")
    require_label_gate(config_path, level=level, parent=ROOT)

    consensus = _read(consensus_dir / "consensus.tsv")
    missing = set(CONSENSUS_COLUMNS).difference(consensus.columns)
    if missing:
        raise ContractError(f"consensus.tsv misses columns: {sorted(missing)}")
    if consensus["cluster_id"].duplicated().any():
        raise ContractError("consensus.tsv holds one row per cluster; a cluster is repeated")
    scores = _read(target / "scores.tsv")
    if set(consensus["cluster_id"]) != set(scores["cluster_id"]):
        raise ContractError(
            "the consensus clusters differ from the clusters of this run's evidence"
        )
    _check_cells(consensus_dir, config, root)

    vocabulary = pd.read_csv(
        plan["vocabulary_path"], sep="\t", keep_default_na=False, dtype=str
    )
    allowed = vocabulary[
        vocabulary["axis"].eq("cell_type")
        & vocabulary["level"].eq(level)
        & vocabulary["parent_scope"].eq(ROOT)
    ]
    outside = sorted(set(consensus["label"]).difference({UNKNOWN, *allowed["label"]}))
    if outside:
        raise ContractError(
            f"consensus labels are outside the approved {level} vocabulary: {outside[:5]}"
        )
    known = dict(zip(allowed["label"], allowed["cl_id"], strict=True))
    wrong_id = consensus[
        consensus["label"].ne(UNKNOWN) & consensus["cl_id"].ne(consensus["label"].map(known))
    ]
    if not wrong_id.empty:
        raise ContractError(
            f"consensus cl_id differs from the vocabulary for clusters {sorted(wrong_id['cluster_id'])[:5]}"
        )

    # `confidence` is the share of the cluster's cells that vote for the label; it is not a calibrated probability.
    labels = pd.DataFrame(
        {
            "level": level,
            "parent_scope": ROOT,
            "cluster_id": consensus["cluster_id"],
            "label": consensus["label"],
            "cell_state": "",
            "cl_id": consensus["cl_id"],
            "confidence": _share(consensus, "cluster_cell_share"),
            "margin": _share(consensus, "cluster_cell_lead"),
            "specificity": "",
            "agreement": "",
            "coverage": "",
            "rationale": "consensus "
            + consensus["decision_rule"]
            + ": "
            + consensus["rationale"],
            "rejected_alternatives": consensus["rejected_alternatives"],
            "reviewer": "",
        }
    ).loc[:, list(LABEL_COLUMNS)]

    copied = target / "consensus"
    copied.mkdir(parents=True, exist_ok=True)
    for name in ("consensus.tsv", "run_manifest.json"):
        if not (copied / name).is_file() or sha256(copied / name) != sha256(
            consensus_dir / name
        ):
            shutil.copyfile(consensus_dir / name, copied / name)
    _frame(target / "labels.tsv", labels)
    return _receipt(
        target / "decision.manifest.json",
        artifacts=[target / "labels.tsv"],
        inputs=[
            target / "evidence.manifest.json",
            target / "scores.tsv",
            copied / "consensus.tsv",
            copied / "run_manifest.json",
        ],
        artifact_kind="cell_curator_corpus_level_decisions",
        level=level,
        decided_by="annotation-consensus",
    )

"""A stand-in annotation-consensus output for tests: its tables, not its computation."""

from __future__ import annotations

import json
from collections.abc import Mapping, Set
from pathlib import Path

import pandas as pd

from cell_curator.annotation import _input
from cell_curator.config import load_config
from cell_curator.data import close_input

COLUMNS = [
    "cluster_id", "label", "cl_id", "decision_rule", "rationale",
    "rejected_alternatives", "cluster_cell_share", "cluster_cell_lead",
]  # fmt: skip


def write_consensus_run(
    config_path: Path,
    root: Path,
    labels: Mapping[str, str],
    *,
    unknown: Set[str] = frozenset(),
) -> Path:
    """Write consensus.tsv, work/cells.tsv and run_manifest.json to the configured directory.

    `labels` maps each cluster id to its L1 label; a cluster in `unknown` is abstained.
    """
    config = load_config(config_path)
    out = Path(config.annotation.consensus_directory)
    plan = json.loads((root / "plan.json").read_text())
    vocabulary = pd.read_csv(plan["vocabulary_path"], sep="\t", keep_default_na=False)
    cl_id = dict(zip(vocabulary["label"], vocabulary["cl_id"], strict=True))
    rows = [
        {
            "cluster_id": cluster,
            "label": "Unknown" if cluster in unknown else label,
            "cl_id": "" if cluster in unknown else cl_id[label],
            "decision_rule": "abstain" if cluster in unknown else "all_agree",
            "rationale": f"synthetic call for cluster {cluster}",
            "rejected_alternatives": "[]",
            "cluster_cell_share": 0.4 if cluster in unknown else 0.9,
            "cluster_cell_lead": 0.1 if cluster in unknown else 0.5,
        }
        for cluster, label in sorted(labels.items())
    ]
    (out / "work").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=COLUMNS).to_csv(out / "consensus.tsv", sep="\t", index=False)
    loaded = _input(config, root)
    try:
        cells = pd.DataFrame(
            {
                "barcode": loaded.root.obs_names.astype(str),
                "cluster_id": loaded.root.obs[config.input.keys.canonical_parent]
                .astype(str)
                .to_numpy(),
            }
        )
    finally:
        close_input(loaded)
    cells.to_csv(out / "work" / "cells.tsv", sep="\t", index=False)
    (out / "run_manifest.json").write_text(json.dumps({"note": "synthetic"}) + "\n")
    return out

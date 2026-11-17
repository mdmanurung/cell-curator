# Purpose: write the two PBMC3K consensus configs into RUN_DIR: run 1 with the frozen demo
#          values (rule h4, marker co-voter) and run 2 with plain majority (no marker vote).
# Run:     .venv/bin/python 03_make-configs.py RUN_DIR
# Inputs:  RUN_DIR/pbmc3k.h5ad                         [01_prep-pbmc3k.py]
#          RUN_DIR/tools/{celltypist,popv,singler}.tsv [02_*]
#          vocabulary.tsv                              34 L1 corpus terms, copied from
#                                                      runs/v8-dev-aida-split-a-003/plan/ (sha256 in vocabulary.sha256)
# Outputs: RUN_DIR/config_h4.yaml, RUN_DIR/config_majority.yaml

import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ONTOLOGY = HERE.parents[3] / "benchmarks/data/aida-phase1-v2-pbmc-hvg2000/cl.obo"
run_dir = Path(sys.argv[1])

# Values of run 1 are the frozen demo values (config.example.yaml); run 2 changes only
# the rule and the marker weight. The user chose both runs on 2026-11-13.
BASE = {
    "schema_version": 1,
    "inputs": {
        "cells": str(run_dir / "pbmc3k.h5ad"),
        "cluster_key": "cluster",
        "vocabulary": str(HERE / "vocabulary.tsv"),
        "ontology": str(ONTOLOGY),
        "tools": {
            t: str(run_dir / "tools" / f"{t}.tsv") for t in ("celltypist", "popv", "singler")
        },
    },
    "votes": {
        "mode": "compute",
        "compute": {
            "matrix_source": "X",
            "matrix_key": "",
            "feature_id_column": "",
            "gene_map": "",
            "level": "L1",
            "parent_scope": "ROOT",
        },
    },
    "smoothing": {"k": 20, "embedding_key": "X_pca"},
    "consensus": {
        "rule": "h4",
        "marker_weight": 1.0,
        "gate": 0.53,
        "tie_break_delta": 0.0,
        "snap_to_vocabulary": True,
        "marker_scope": "cell",
        "tool_order": ["celltypist", "popv", "singler"],
        "values_status": "demo values frozen on AIDA split A; approved for the PBMC3K showcase by the user, 2026-11-13; rule h4 (tie-break dropped) chosen by the user, 2026-11-17",
    },
    "cluster_gate": {"share": 0.5, "lead": 0.2, "fallback_margin": 0.40},
    "determinism": {"threads": 1, "float_digits": 6},
}

for name, rule, weight in (("h4", "h4", 1.0), ("majority", "majority", 0.0)):
    cfg = yaml.safe_load(yaml.safe_dump(BASE))
    cfg["run"] = {"name": f"pbmc3k-{name}"}
    cfg["consensus"].update(rule=rule, marker_weight=weight)
    if name == "majority":
        cfg["consensus"]["values_status"] = (
            "plain majority baseline (no marker vote); approved for the PBMC3K showcase by the user, 2026-11-13"
        )
    (run_dir / f"config_{name}.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))

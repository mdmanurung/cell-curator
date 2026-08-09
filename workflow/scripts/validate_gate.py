"""Validate the real guided assumptions gate and publish its durable receipt."""

import traceback
from pathlib import Path
from typing import Any

from cell_curator.contracts import sha256
from cell_curator.guided import load_assumptions, require_label_gate
from cell_curator.provenance import atomic_publish, publish_json

workflow_context: Any = globals()["snakemake"]


def validate_gate() -> None:
    config_path = Path(str(workflow_context.input.config))
    assumptions_path = Path(str(workflow_context.input.assumptions))
    plan_path = Path(str(workflow_context.input.plan))
    level = str(workflow_context.params.level)
    parent = str(workflow_context.params.parent)

    state = require_label_gate(config_path, level=level, parent=parent)
    rows = load_assumptions(assumptions_path.parent)
    applicable_scopes = {"global", f"{level}::{parent}"}
    approved = [
        row
        for row in rows
        if bool(row.get("material", True)) and str(row.get("scope")) in applicable_scopes
    ]
    publish_json(
        Path(str(workflow_context.output.receipt)),
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_validated_assumption_gate",
            "status": "complete",
            "run_id": state.run_id,
            "level": level,
            "parent": parent,
            "labels_permitted": state.labels_permitted,
            "interactive": state.interactive,
            "preauthorized": state.preauthorized,
            "approval_modes": sorted({str(row.get("approval_mode", "")) for row in approved}),
            "reviewers": sorted({str(row.get("reviewer", "")) for row in approved}),
            "approved_material_assumptions": len(approved),
            "input_sha256": state.input_sha256,
            "config_sha256": state.config_sha256,
            "plan_sha256": sha256(plan_path),
            "assumptions_sha256": sha256(assumptions_path),
        },
    )


try:
    validate_gate()
except Exception:
    atomic_publish(Path(str(workflow_context.log[0])), traceback.format_exc().encode())
    raise
else:
    atomic_publish(
        Path(str(workflow_context.log[0])),
        b"validated cell-curator assumptions gate\n",
    )

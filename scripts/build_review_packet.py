#!/usr/bin/env python3
"""Assemble a versioned immutable review packet after all evidence gates complete."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from contracts import ContractError, read_config, require_blank_labels, sha256


REQUIRED_FILES = (
    "context.md",
    "assumptions.md",
    "plan.md",
    "run_state.json",
    "input_contract.manifest.json",
    "adaptive_parent_audit.tsv",
    "discovery/candidate_registry.tsv",
    "discovery/parent_triage.tsv",
    "discovery/discovery.manifest.json",
    "contrasts/comparator_manifest.tsv",
    "gpu_markers/gpu_receipts.validated.json",
    "parent_decision_summary.tsv",
    "child_evidence_summary.tsv",
    "parent_child_mapping.tsv",
    "leaf_identity_totals.tsv",
    "annotation_review_template.tsv",
    "critic_findings.tsv",
    "critic_reconciliation.tsv",
    "final_provenance.manifest.json",
)
OPTIONAL_DIRECTORIES = (
    "evidence_cards",
    "parent_cards",
    "figures",
    "summaries",
    "resolution_persistence",
)


def _labels_present(frame: pd.DataFrame) -> bool:
    columns = [
        column
        for column in frame.columns
        if str(column).startswith("biological_label_")
        or column in {"cell_type", "cl_id"}
    ]
    return any(
        frame[column].fillna("").astype(str).str.strip().ne("").any()
        for column in columns
    )


def validate_review_gate(config: dict, run_root: Path) -> dict:
    missing = [relative for relative in REQUIRED_FILES if not (run_root / relative).is_file()]
    if missing:
        raise ContractError(f"review packet inputs are missing: {missing}")
    state = json.loads((run_root / "run_state.json").read_text())
    mode = str(config["run"]["mode"])
    if state.get("mode") != mode:
        raise ContractError("run state mode drifted from configuration")
    tables = [
        pd.read_csv(run_root / relative, sep="\t", dtype=str, keep_default_na=False)
        for relative in (
            "discovery/candidate_registry.tsv",
            "parent_decision_summary.tsv",
            "child_evidence_summary.tsv",
            "annotation_review_template.tsv",
        )
    ]
    labels_present = any(_labels_present(frame) for frame in tables)
    if mode == "evidence-only":
        for frame in tables:
            require_blank_labels(frame)
        if labels_present:
            raise ContractError("evidence-only mode contains prohibited labels or CL IDs")
    elif labels_present:
        gates = (
            "gate_a_approved",
            "knowledgebase_available",
            "critic_reconciled",
            "human_approval",
        )
        if not all(state.get(gate) is True for gate in gates):
            raise ContractError("biological labels exist before all review gates passed")
    findings = pd.read_csv(
        run_root / "critic_findings.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    reconciliation = pd.read_csv(
        run_root / "critic_reconciliation.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    if not findings.empty:
        if "finding_id" not in findings or "finding_id" not in reconciliation:
            raise ContractError("critic artifacts lack finding IDs")
        unresolved = set(findings["finding_id"].astype(str)).difference(
            reconciliation.loc[
                reconciliation.get("status", pd.Series(dtype=str)).astype(str).eq("RESOLVED"),
                "finding_id",
            ].astype(str)
        )
        if unresolved:
            raise ContractError(f"critic findings remain unresolved: {sorted(unresolved)}")
    mapping_validation = run_root / "mapping.validation.json"
    if mapping_validation.is_file():
        value = json.loads(mapping_validation.read_text())
        if value.get("status") != "complete":
            raise ContractError("mapping validation is incomplete")
    return {
        "mode": mode,
        "labels_present": labels_present,
        "critic_findings_reconciled": True,
        "writeback_permitted": bool(
            mode == "annotation" and state.get("human_approval") is True
        ),
    }


def build_packet(config_path: Path, run_root: Path) -> Path:
    config = read_config(config_path)
    gate = validate_review_gate(config, run_root)
    packet = run_root / "review_packet"
    if packet.exists():
        raise ContractError(f"refusing to overwrite completed packet: {packet}")
    stage = packet.with_name(f".{packet.name}.partial.{os.getpid()}")
    if stage.exists():
        raise ContractError(f"stale packet stage exists: {stage}")
    stage.mkdir(parents=True)
    try:
        copied: list[Path] = []
        for relative in REQUIRED_FILES:
            source = run_root / relative
            target = stage / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            copied.append(target)
        if (run_root / "mapping.validation.json").is_file():
            shutil.copy2(
                run_root / "mapping.validation.json",
                stage / "mapping.validation.json",
            )
            copied.append(stage / "mapping.validation.json")
        for name in OPTIONAL_DIRECTORIES:
            source = run_root / name
            if source.is_dir():
                shutil.copytree(source, stage / name)
                copied.extend(path for path in (stage / name).rglob("*") if path.is_file())
        manifest = {
            "schema_version": 1,
            "artifact_kind": "adaptive_multimodal_review_packet",
            "status": "complete",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "run_id": config["run"]["run_id"],
            "mode": gate["mode"],
            "labels_present": gate["labels_present"],
            "writeback_permitted": gate["writeback_permitted"],
            "critic_findings_reconciled": True,
            "immutable": True,
            "files": {
                str(path.relative_to(stage)): {
                    "sha256": sha256(path),
                    "bytes": path.stat().st_size,
                }
                for path in sorted(copied)
            },
        }
        manifest_path = stage / "review_packet.manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        os.replace(stage, packet)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return packet


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        build_packet(args.config, args.run_root)
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()

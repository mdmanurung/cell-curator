"""Derive refinement gates, parent outcomes, and immutable cell mappings from evidence."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from .capabilities import build_capability_registry
from .config import CellCuratorConfig
from .contracts import ContractError, as_bool, sha256, values_sha256
from .evidence import lane_key, validate_backend_receipt
from .models import BackendReceipt, CapabilityLevel, MarkerProgram, ParentOutcome
from .provenance import atomic_publish, publish_json

COMPARATOR_COLUMNS = [
    "schema_version",
    "run_id",
    "lane_id",
    "scope",
    "parent_id",
    "candidate_key",
    "candidate_child_uid",
    "comparator_role",
    "assay",
    "representation",
    "lane_key",
    "backend",
    "required_for_acceptance",
    "groups_path",
    "target_n_cells",
    "comparator_n_cells",
    "target_barcode_sha256",
    "comparator_barcode_sha256",
]


def _row_positions(index: pd.Index) -> dict[str, int]:
    return {str(value): position for position, value in enumerate(index)}


def _subset_matrix(matrix: Any, positions: np.ndarray) -> Any:
    return matrix[positions]


def _binary_separation(matrix: Any, target: np.ndarray, reference: np.ndarray) -> float:
    target_values = matrix[target]
    reference_values = matrix[reference]
    target_mean = np.asarray(target_values.mean(axis=0)).ravel()
    reference_mean = np.asarray(reference_values.mean(axis=0)).ravel()
    difference = target_mean - reference_mean
    if sparse.issparse(matrix):
        scale = np.sqrt(np.asarray(matrix.power(2).mean(axis=0)).ravel() + 1e-8)
    else:
        scale = np.sqrt(np.mean(np.asarray(matrix) ** 2, axis=0) + 1e-8)
    standardized = difference / scale
    return float(np.sqrt(np.mean(standardized**2)))


def build_comparator_manifest(
    *,
    config: CellCuratorConfig,
    run_root: Path,
) -> pd.DataFrame:
    """Freeze every candidate x comparator x assay lane from computed membership."""

    cells = pd.read_csv(run_root / "prepared" / "cells.tsv", sep="\t", keep_default_na=False)
    registry = pd.read_csv(
        run_root / "discovery" / "candidate_registry.tsv", sep="\t", keep_default_na=False
    )
    membership = pd.read_csv(
        run_root / "discovery" / "candidate_membership.tsv", sep="\t", keep_default_na=False
    )
    # Snakemake declares this directory as a checkpoint output. It must exist
    # even when every parent is retained and there are zero candidate lanes.
    (run_root / "comparators" / "groups").mkdir(parents=True, exist_ok=True)
    member_sets = {
        str(candidate): set(
            block.loc[
                block["candidate_membership"].astype(str).str.lower().isin({"true", "1"}),
                "barcode",
            ].astype(str)
        )
        for candidate, block in membership.groupby("candidate_key", observed=True)
    }
    rows: list[dict[str, Any]] = []
    for candidate in registry.to_dict("records"):
        candidate_key = str(candidate["candidate_key"])
        block = membership[membership["candidate_key"].astype(str).eq(candidate_key)]
        target = member_sets.get(candidate_key, set())
        parent_cells = set(block["barcode"].astype(str))
        comparator_sets = {
            "parent_remainder": parent_cells - target,
            "closest_sibling": member_sets.get(str(candidate["closest_sibling"]), set()),
        }
        for role in config.comparators.roles:
            comparator = comparator_sets[role]
            if not target or not comparator or target.intersection(comparator):
                raise ContractError(
                    f"candidate {candidate_key} has an invalid {role} comparator"
                )
            selected = cells[cells["barcode"].astype(str).isin(target | comparator)].copy()
            selected["analysis_group"] = (
                selected["barcode"]
                .astype(str)
                .map(
                    lambda barcode, target_set=target: (
                        "candidate" if barcode in target_set else "reference"
                    )
                )
            )
            selected["analysis_role"] = selected["analysis_group"].map(
                {"candidate": "candidate_child", "reference": role}
            )
            selected = selected.sort_values("barcode", ignore_index=True)
            groups_path = run_root / "comparators" / "groups" / candidate_key / f"{role}.tsv"
            atomic_publish(
                groups_path,
                selected.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
            )
            for lane in config.evidence.lanes:
                representation = config.input.assays[lane.assay].representations[
                    lane.representation
                ]
                configured_key = lane_key(lane.assay, lane.representation, lane.lane_key)
                lane_id = f"{candidate_key}__{role}__{configured_key}"
                rows.append(
                    {
                        "schema_version": 3,
                        "run_id": config.run.run_id,
                        "lane_id": lane_id,
                        "scope": str(candidate["scope"]),
                        "parent_id": str(candidate["parent_id"]),
                        "candidate_key": candidate_key,
                        "candidate_child_uid": str(candidate["candidate_child_uid"]),
                        "comparator_role": role,
                        "assay": lane.assay,
                        "representation": lane.representation,
                        "lane_key": configured_key,
                        "backend": representation.ranking_backend,
                        "required_for_acceptance": lane.required_for_acceptance,
                        "groups_path": str(groups_path),
                        "target_n_cells": len(target),
                        "comparator_n_cells": len(comparator),
                        "target_barcode_sha256": values_sha256(sorted(target)),
                        "comparator_barcode_sha256": values_sha256(sorted(comparator)),
                    }
                )
    result = pd.DataFrame(rows, columns=COMPARATOR_COLUMNS)
    expected = len(registry) * len(config.comparators.roles) * len(config.evidence.lanes)
    if len(result) != expected or result["lane_id"].duplicated().any():
        raise ContractError("computed comparator manifest is incomplete or duplicated")
    path = run_root / "comparators" / "comparator_manifest.tsv"
    atomic_publish(path, result.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    publish_json(
        run_root / "comparators" / "comparator_manifest.json",
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_dynamic_comparator_manifest",
            "run_id": config.run.run_id,
            "status": "complete",
            "n_candidates": len(registry),
            "n_lanes": len(result),
            "zero_applicable_lane_run": len(result) == 0,
            "manifest_sha256": sha256(path),
        },
    )
    return result


def derive_candidate_evidence(
    *,
    config: CellCuratorConfig,
    run_root: Path,
    programs: list[MarkerProgram],
) -> pd.DataFrame:
    registry_path = run_root / "discovery" / "candidate_registry.tsv"
    membership_path = run_root / "discovery" / "candidate_membership.tsv"
    comparator_path = run_root / "comparators" / "comparator_manifest.tsv"
    if (
        not registry_path.is_file()
        or not membership_path.is_file()
        or not comparator_path.is_file()
    ):
        raise ContractError(
            "adaptive candidate registry, membership, or comparator manifest is missing"
        )
    registry = pd.read_csv(registry_path, sep="\t", keep_default_na=False)
    if registry.empty:
        result = pd.DataFrame(
            columns=[
                "schema_version",
                "scope",
                "parent_id",
                "candidate_key",
                "candidate_child_uid",
                "accepted_child_evidence_gate",
                "failed_gates",
                "cell_type",
                "cl_id",
            ]
        )
        path = run_root / "decisions" / "child_evidence_summary.tsv"
        atomic_publish(path, result.to_csv(sep="\t", index=False, lineterminator="\n").encode())
        publish_json(
            run_root / "decisions" / "candidate_evidence.manifest.json",
            {
                "schema_version": 3,
                "artifact_kind": "cell_curator_candidate_evidence",
                "run_id": config.run.run_id,
                "status": "complete",
                "derived_from_raw_artifacts": True,
                "n_candidates": 0,
                "n_accepted": 0,
                "zero_candidate_run": True,
                "comparator_manifest_sha256": sha256(comparator_path),
                "output_sha256": sha256(path),
            },
        )
        return result
    comparator_manifest = pd.read_csv(comparator_path, sep="\t", keep_default_na=False)
    capability_index = {
        (item.assay, item.representation): item for item in build_capability_registry(config)
    }
    input_contract = json.loads((run_root / "input_contract.manifest.json").read_text())
    rows: list[dict[str, Any]] = []
    for candidate in registry.to_dict("records"):
        candidate_key = str(candidate["candidate_key"])
        failed: list[str] = []
        lane_rows: list[dict[str, Any]] = []
        winners: list[str] = []
        state_winners: list[str] = []
        role_results: dict[str, list[bool]] = {role: [] for role in config.comparators.roles}
        for lane in config.evidence.lanes:
            capability = capability_index[(lane.assay, lane.representation)]
            configured_key = lane_key(lane.assay, lane.representation, lane.lane_key)
            lane_passes: list[bool] = []
            for role in config.comparators.roles:
                selected = comparator_manifest[
                    comparator_manifest["candidate_key"].astype(str).eq(candidate_key)
                    & comparator_manifest["comparator_role"].astype(str).eq(role)
                    & comparator_manifest["lane_key"].astype(str).eq(configured_key)
                ]
                if len(selected) != 1:
                    raise ContractError(
                        f"comparator manifest lacks {candidate_key}/{role}/{configured_key}"
                    )
                lane_id = str(selected.iloc[0]["lane_id"])
                lane_dir = run_root / "evidence" / "candidate_lanes" / lane_id
                evidence_path = lane_dir / "program_evidence.tsv"
                markers_path = lane_dir / "bottom_up_markers.tsv"
                comparison_path = lane_dir / "comparison.json"
                receipt_path = lane_dir / "receipt.json"
                missing = [
                    path
                    for path in (evidence_path, markers_path, comparison_path, receipt_path)
                    if not path.is_file()
                ]
                if missing:
                    raise ContractError(f"candidate evidence lane is incomplete: {missing}")
                receipt = BackendReceipt.model_validate_json(receipt_path.read_text())
                if (
                    receipt.status != "complete"
                    or receipt.fallback_used
                    or receipt.input_sha256 != input_contract["input_sha256"]
                    or receipt.config_sha256 != input_contract["config_sha256"]
                    or receipt.output_sha256 != sha256(evidence_path)
                ):
                    raise ContractError(
                        f"candidate evidence receipt failed validation: {lane_id}"
                    )
                validate_backend_receipt(receipt, capability)
                scored = pd.read_csv(evidence_path, sep="\t", keep_default_na=False)
                markers = pd.read_csv(markers_path, sep="\t", keep_default_na=False)
                comparison = json.loads(comparison_path.read_text())
                if comparison.get("groups_sha256") != sha256(
                    Path(str(selected.iloc[0]["groups_path"]))
                ):
                    raise ContractError(f"candidate comparison group hash drifted: {lane_id}")
                parent_name = config.evidence.parent_vocabulary.get(
                    str(candidate["parent_id"]), str(candidate["parent_id"])
                )
                target = scored[
                    scored["cluster_id"].astype(str).eq("candidate")
                    & scored.get("axis", pd.Series("cell_type", index=scored.index))
                    .astype(str)
                    .eq("cell_type")
                    & scored["parent_label"]
                    .astype(str)
                    .isin({"", str(candidate["parent_id"]), parent_name})
                ].sort_values(["score", "label"], ascending=[False, True])
                states = scored[
                    scored["cluster_id"].astype(str).eq("candidate")
                    & scored.get("axis", pd.Series("cell_type", index=scored.index))
                    .astype(str)
                    .eq("cell_state")
                ].sort_values(["score", "label"], ascending=[False, True])
                marker_target = markers[markers["cluster_id"].astype(str).eq("candidate")]
                if marker_target.empty:
                    marker_significant = False
                elif capability.requires_gpu:
                    if "method" not in marker_target:
                        raise ContractError(f"RAPIDS marker table lacks method rows: {lane_id}")
                    wilcoxon = marker_target[
                        marker_target["method"].astype(str).eq("wilcoxon")
                        & (pd.to_numeric(marker_target["effect"], errors="coerce") > 0)
                        & (
                            pd.to_numeric(
                                marker_target["adjusted_pvalue"], errors="coerce"
                            )
                            <= config.evidence.adjusted_pvalue_max
                        )
                    ]
                    logreg = marker_target[
                        marker_target["method"].astype(str).eq("logreg")
                        & (pd.to_numeric(marker_target["effect"], errors="coerce") > 0)
                    ]
                    marker_significant = bool(
                        set(wilcoxon["feature"].astype(str)).intersection(
                            logreg["feature"].astype(str)
                        )
                    )
                else:
                    marker_significant = bool(
                        (
                            (pd.to_numeric(marker_target["effect"], errors="coerce") > 0)
                            & (
                                pd.to_numeric(
                                    marker_target["adjusted_pvalue"], errors="coerce"
                                )
                                <= config.evidence.adjusted_pvalue_max
                            )
                        ).any()
                    )
                separation = float(comparison["separation"])
                if target.empty:
                    winner = cl_id = ""
                    coverage = violation = 0.0
                    passed = False
                else:
                    best = target.iloc[0]
                    winner = str(best["label"])
                    cl_id = str(best["cl_id"])
                    coverage = float(best["positive_coverage"])
                    violation = float(best["negative_violation"])
                    program_coverage = float(best["program_coverage"])
                    detection_valid = as_bool(best["detection_valid"])
                    passed = bool(
                        marker_significant
                        and separation >= 0.1
                        and detection_valid
                        and program_coverage >= config.evidence.min_positive_coverage
                        and coverage >= config.evidence.min_positive_coverage
                        and violation <= config.evidence.max_negative_violation
                    )
                if passed and capability.level is CapabilityLevel.ANNOTATION:
                    winners.append(winner)
                    if not states.empty:
                        state_winners.append(str(states.iloc[0]["label"]))
                lane_passes.append(passed)
                role_results[role].append(passed)
                lane_rows.append(
                    {
                        "lane_id": lane_id,
                        "lane_key": configured_key,
                        "comparator_role": role,
                        "capability_level": capability.level.value,
                        "pass": passed,
                        "winner": winner,
                        "cl_id": cl_id,
                        "separation": separation,
                        "positive_coverage": coverage,
                        "program_coverage": (
                            program_coverage if not target.empty else 0.0
                        ),
                        "detection_valid": detection_valid if not target.empty else False,
                        "negative_violation": violation,
                        "marker_significant": marker_significant,
                        "receipt_sha256": sha256(receipt_path),
                    }
                )
            if lane.required_for_acceptance and not all(lane_passes):
                failed.append(f"lane:{configured_key}")
        stability_pass = bool(
            float(candidate["stability_score"])
            > config.adaptive.parent_local.stability_threshold_strictly_greater_than
        )
        limit = config.adaptive.audit.technical_effect_size_max
        technical_pass = bool(
            limit is None or float(candidate["technical_effect_size"]) < limit
        )
        gate_values = {
            "stability_pass": stability_pass,
            "candidate_vs_parent_remainder_pass": bool(
                role_results.get("parent_remainder")
            )
            and all(role_results["parent_remainder"]),
            "candidate_vs_closest_sibling_pass": bool(
                role_results.get("closest_sibling")
            )
            and all(role_results["closest_sibling"]),
            "positive_program_pass": any(
                item["detection_valid"]
                and item["program_coverage"] >= config.evidence.min_positive_coverage
                and item["positive_coverage"] >= config.evidence.min_positive_coverage
                for item in lane_rows
            ),
            "expected_negative_program_pass": all(
                item["detection_valid"]
                and item["negative_violation"] <= config.evidence.max_negative_violation
                for item in lane_rows
            )
            and bool(lane_rows),
            "cell_level_separation_pass": all(item["separation"] >= 0.1 for item in lane_rows),
            "source_corroboration_pass": bool(winners),
            "technical_confounding_pass": technical_pass,
        }
        for gate in config.evidence.required_candidate_gates:
            if gate not in gate_values:
                failed.append(f"unsupported_required_gate:{gate}")
            elif not gate_values[gate]:
                failed.append(gate)
        if not stability_pass:
            failed.append("stability")
        if not technical_pass:
            failed.append("technical_confounding")
        winner_counts = (
            pd.Series(winners).value_counts().rename_axis("label").reset_index(name="n")
            if winners
            else pd.DataFrame(columns=["label", "n"])
        )
        if not winner_counts.empty:
            winner_counts = winner_counts.sort_values(["n", "label"], ascending=[False, True])
            consensus_label = str(winner_counts.iloc[0]["label"])
            if int(winner_counts.iloc[0]["n"]) / len(winners) < 0.5:
                failed.append("multimodal_label_disagreement")
        else:
            consensus_label = ""
            failed.append("no_annotation_capable_identity")
        if config.run.mode.value != "annotation":
            failed.append("run_mode_evidence_only")
        cl_ids = sorted(
            {
                str(item["cl_id"])
                for item in lane_rows
                if item["winner"] == consensus_label and str(item["cl_id"])
            }
        )
        if len(cl_ids) > 1:
            failed.append("conflicting_ontology_accessions")
        state_counts = pd.Series(state_winners).value_counts() if state_winners else None
        consensus_state = (
            str(state_counts.index[0]) if state_counts is not None else "not_assessed"
        )
        accepted = not failed
        rows.append(
            {
                "schema_version": 3,
                "scope": str(candidate["scope"]),
                "parent_id": str(candidate["parent_id"]),
                "candidate_key": candidate_key,
                "candidate_child_uid": str(candidate["candidate_child_uid"]),
                "n_cells": int(candidate["n_cells_in_parent"]),
                "stability_score": float(candidate["stability_score"]),
                "separation_score": float(candidate["separation_score"]),
                "technical_effect_size": float(candidate["technical_effect_size"]),
                **gate_values,
                "accepted_child_evidence_gate": accepted,
                "failed_gates": ";".join(sorted(set(failed))),
                "cell_type": consensus_label if accepted else "",
                "cell_state": consensus_state,
                "cl_id": cl_ids[0] if accepted and len(cl_ids) == 1 else "",
                "lane_evidence_json": json.dumps(lane_rows, sort_keys=True),
                "registry_sha256": sha256(registry_path),
                "membership_sha256": sha256(membership_path),
            }
        )
    result = pd.DataFrame(rows)
    decisions_dir = run_root / "decisions"
    decisions_dir.mkdir(parents=True, exist_ok=True)
    path = decisions_dir / "child_evidence_summary.tsv"
    atomic_publish(path, result.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    publish_json(
        decisions_dir / "candidate_evidence.manifest.json",
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_candidate_evidence",
            "run_id": config.run.run_id,
            "status": "complete",
            "derived_from_raw_artifacts": True,
            "n_candidates": len(result),
            "n_accepted": int(result["accepted_child_evidence_gate"].sum())
            if len(result)
            else 0,
            "comparator_manifest_sha256": sha256(comparator_path),
            "output_sha256": sha256(path),
        },
    )
    return result


def decide_parents_from_artifacts(
    *,
    config: CellCuratorConfig,
    run_root: Path,
) -> pd.DataFrame:
    cells = pd.read_csv(run_root / "prepared" / "cells.tsv", sep="\t", keep_default_na=False)
    audit_path = run_root / "adaptive_parent_audit.tsv"
    audit = (
        pd.read_csv(audit_path, sep="\t", keep_default_na=False)
        if audit_path.is_file()
        else pd.DataFrame()
    )
    children_path = run_root / "decisions" / "child_evidence_summary.tsv"
    children = pd.read_csv(children_path, sep="\t", keep_default_na=False)
    rows = []
    for (scope, parent), _ in cells.groupby(["scope", "parent_id"], observed=True):
        block = children[
            children["scope"].astype(str).eq(str(scope))
            & children["parent_id"].astype(str).eq(str(parent))
        ]
        accepted = block[
            block["accepted_child_evidence_gate"].astype(str).str.lower().isin({"true", "1"})
        ]
        audit_row = (
            audit[
                audit["scope"].astype(str).eq(str(scope))
                & audit["parent_id"].astype(str).eq(str(parent))
            ]
            if not audit.empty
            else pd.DataFrame()
        )
        if not accepted.empty:
            outcome = ParentOutcome.ACCEPT_SPLIT
            reason = "one or more nested children passed two-comparator, modality, stability, and technical gates"
        elif not audit_row.empty and bool(
            str(audit_row.iloc[0].get("doublet_enriched", "false")).lower() in {"true", "1"}
            or str(audit_row.iloc[0].get("incompatible_program_enriched", "false")).lower()
            in {"true", "1"}
        ):
            outcome = ParentOutcome.DOUBLET_MIXED
            reason = "computed doublet or incompatible-program audit evidence"
        elif not audit_row.empty and any(
            str(audit_row.iloc[0].get(key, "false")).lower() in {"true", "1"}
            for key in ("donor_dominated", "capture_dominated", "qc_outlier")
        ):
            outcome = ParentOutcome.TECHNICAL_UNRESOLVED
            reason = "computed donor, capture, or QC structure dominated"
        else:
            outcome = ParentOutcome.RETAIN_PARENT
            reason = "no candidate passed the complete provenance-derived evidence gate"
        rows.append(
            {
                "schema_version": 3,
                "scope": str(scope),
                "parent_id": str(parent),
                "parent_uid": f"{scope}::{parent}",
                "recommended_outcome": outcome.value,
                "decision_outcome": "",
                "accepted_child_ids_recommended": ";".join(
                    accepted["candidate_child_uid"].astype(str)
                ),
                "accepted_child_ids": "",
                "decision_reasons": reason,
                "manual_review_status": "PENDING_REVIEW",
            }
        )
    result = pd.DataFrame(rows)
    if result.duplicated(["scope", "parent_id"]).any():
        raise ContractError("every canonical parent must have exactly one decision")
    path = run_root / "decisions" / "parent_decision_summary.tsv"
    atomic_publish(path, result.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    return result


def build_mapping(
    *,
    config: CellCuratorConfig,
    run_root: Path,
    use_recommendations: bool = True,
) -> pd.DataFrame:
    cells = pd.read_csv(run_root / "prepared" / "cells.tsv", sep="\t", keep_default_na=False)
    proposals = pd.read_csv(
        run_root / "evidence" / "summary" / "annotation_proposals.tsv",
        sep="\t",
        keep_default_na=False,
    )
    decisions = pd.read_csv(
        run_root / "decisions" / "parent_decision_summary.tsv",
        sep="\t",
        keep_default_na=False,
    )
    child_summary = pd.read_csv(
        run_root / "decisions" / "child_evidence_summary.tsv",
        sep="\t",
        keep_default_na=False,
    )
    membership = pd.read_csv(
        run_root / "discovery" / "candidate_membership.tsv",
        sep="\t",
        keep_default_na=False,
    )
    proposal_index = proposals.set_index(proposals["cluster_id"].astype(str), drop=False)
    child_by_uid = (
        child_summary.set_index("candidate_child_uid", drop=False)
        if not child_summary.empty
        else None
    )
    rows: list[dict[str, Any]] = []
    for decision in decisions.to_dict("records"):
        scope = str(decision["scope"])
        parent = str(decision["parent_id"])
        parent_uid = str(decision["parent_uid"])
        outcome = str(
            decision["recommended_outcome"]
            if use_recommendations
            else decision["decision_outcome"]
        )
        accepted_value = str(
            decision["accepted_child_ids_recommended"]
            if use_recommendations
            else decision["accepted_child_ids"]
        )
        accepted = [item for item in accepted_value.split(";") if item]
        parent_cells = cells[
            cells["scope"].astype(str).eq(scope) & cells["parent_id"].astype(str).eq(parent)
        ]
        assignments: dict[str, tuple[str, str, str]] = {}
        for child_uid in accepted:
            if child_by_uid is None or child_uid not in child_by_uid.index:
                raise ContractError(
                    f"accepted child is absent from evidence summary: {child_uid}"
                )
            child = child_by_uid.loc[child_uid]
            candidate_key = str(child["candidate_key"])
            selected = membership[
                membership["candidate_key"].astype(str).eq(candidate_key)
                & membership["candidate_membership"].astype(str).str.lower().isin({"true", "1"})
            ]
            for barcode in selected["barcode"].astype(str):
                if barcode in assignments:
                    raise ContractError("accepted children overlap")
                assignments[barcode] = (
                    child_uid,
                    str(child["cell_type"]),
                    str(child.get("cell_state", "normal")),
                )
        parent_proposal = proposal_index.loc[parent] if parent in proposal_index.index else None
        parent_label = (
            str(parent_proposal["cell_type"]) if parent_proposal is not None else "Unknown"
        )
        parent_state = (
            str(parent_proposal["cell_state"]) if parent_proposal is not None else "normal"
        )
        for barcode in parent_cells["barcode"].astype(str):
            if barcode in assignments:
                leaf, cell_type, cell_state = assignments[barcode]
                role = "accepted_child"
            elif accepted:
                leaf = f"{parent_uid}::remainder"
                cell_type, cell_state = parent_label, parent_state
                role = "parent_remainder"
            else:
                leaf = parent_uid
                cell_type, cell_state = parent_label, parent_state
                role = "retained_parent"
            rows.append(
                {
                    "schema_version": 3,
                    "run_id": config.run.run_id,
                    "barcode": barcode,
                    "scope": scope,
                    "canonical_parent": parent,
                    "parent_uid": parent_uid,
                    "leaf_id": leaf,
                    "cell_type": cell_type,
                    "cell_state": cell_state,
                    "membership_role": role,
                    "decision_outcome": outcome,
                }
            )
    mapping = pd.DataFrame(rows).sort_values("barcode", ignore_index=True)
    validate_mapping(cells, mapping)
    path = run_root / "mapping" / "parent_child_mapping.tsv"
    atomic_publish(path, mapping.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    totals = (
        mapping.groupby("scope", observed=True)
        .agg(
            n_cells=("barcode", "size"),
            n_canonical_parents=("parent_uid", "nunique"),
            n_leaves=("leaf_id", "nunique"),
        )
        .reset_index()
    )
    atomic_publish(
        run_root / "mapping" / "leaf_identity_totals.tsv",
        totals.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    publish_json(
        run_root / "mapping" / "mapping.validation.json",
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_final_mapping_validation",
            "status": "complete",
            "n_cells": len(mapping),
            "n_leaves": int(mapping["leaf_id"].nunique()),
            "mapping_sha256": sha256(path),
            "canonical_membership_sha256": values_sha256(
                [
                    f"{row.barcode}\t{row.scope}\t{row.parent_id}"
                    for row in cells.itertuples(index=False)
                ]
            ),
            "invariants": {
                "one_row_per_input_cell": True,
                "canonical_parent_membership_unchanged": True,
                "parent_counts_conserved": True,
                "unsplit_membership_unchanged": True,
            },
        },
    )
    return mapping


def validate_mapping(cells: pd.DataFrame, mapping: pd.DataFrame) -> None:
    if mapping["barcode"].duplicated().any() or set(mapping["barcode"].astype(str)) != set(
        cells["barcode"].astype(str)
    ):
        raise ContractError("mapping must contain exactly one row per input cell")
    expected = cells.set_index(cells["barcode"].astype(str))
    observed = mapping.set_index(mapping["barcode"].astype(str)).loc[expected.index]
    if not observed["scope"].astype(str).equals(expected["scope"].astype(str)):
        raise ContractError("mapping changed canonical scope")
    if not observed["canonical_parent"].astype(str).equals(expected["parent_id"].astype(str)):
        raise ContractError("mapping changed canonical parent membership")
    expected_counts = (
        expected.assign(
            scope=expected["scope"].astype(str),
            parent_id=expected["parent_id"].astype(str),
        )
        .groupby(["scope", "parent_id"], observed=True)
        .size()
    )
    observed_counts = (
        observed.assign(
            scope=observed["scope"].astype(str),
            canonical_parent=observed["canonical_parent"].astype(str),
        )
        .groupby(["scope", "canonical_parent"], observed=True)
        .size()
    )
    observed_counts.index.names = expected_counts.index.names
    if not expected_counts.equals(observed_counts):
        raise ContractError("parent counts are not conserved")
    unsplit = observed[observed["membership_role"].eq("retained_parent")]
    expected_leaf = (
        unsplit["scope"].astype(str) + "::" + unsplit["canonical_parent"].astype(str)
    )
    if not unsplit["leaf_id"].astype(str).equals(expected_leaf):
        raise ContractError("unsplit parent membership changed")

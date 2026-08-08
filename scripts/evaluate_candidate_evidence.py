#!/usr/bin/env python3
"""Apply multimodal, technical, stability, and doublet gates to candidate evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from contracts import (
    OUTCOMES,
    ContractError,
    as_bool,
    atomic_write_tsv,
    evidence_lanes,
    lane_gate_column,
    read_config,
    require_blank_labels,
    require_columns,
)


def evaluate_candidates(config: dict, evidence: pd.DataFrame) -> pd.DataFrame:
    configured_gates = list(config["evidence"].get("required_candidate_gates", []))
    lane_gates = [
        lane_gate_column(lane["lane_key"])
        for lane in evidence_lanes(config)
        if lane["required_for_acceptance"]
    ]
    boolean_gates = tuple(dict.fromkeys([*configured_gates, *lane_gates]))
    required = {
        "scope",
        "candidate_key",
        "parent_id",
        "parent_uid",
        "candidate_child_uid",
        "n_cells",
        *boolean_gates,
    }
    eligibility = config["adaptive"]["eligibility"]
    audit = config["adaptive"]["audit"]
    conditional = {
        "n_donors": eligibility.get("min_donors"),
        "n_captures": eligibility.get("min_captures"),
        "max_donor_fraction": audit.get("donor_max_fraction"),
        "max_capture_fraction": audit.get("capture_max_fraction"),
        "technical_effect_size": audit.get("technical_effect_size_max"),
        "doublet_fraction": audit.get("doublet_fraction_max"),
        "incompatible_coexpression_fraction": audit.get(
            "incompatible_program_fraction_min"
        ),
    }
    required.update(name for name, threshold in conditional.items() if threshold is not None)
    require_columns(evidence, required, artifact="candidate evidence")
    require_blank_labels(evidence)
    rows = []
    for row in evidence.to_dict("records"):
        gates = {name: as_bool(row[name]) for name in boolean_gates}
        gates["minimum_cells"] = int(row["n_cells"]) >= int(
            eligibility["min_cells"]
        )
        optional_checks = (
            ("minimum_donors", "n_donors", eligibility.get("min_donors"), ">="),
            ("minimum_captures", "n_captures", eligibility.get("min_captures"), ">="),
            ("donor_nondominance", "max_donor_fraction", audit.get("donor_max_fraction"), "<"),
            ("capture_nondominance", "max_capture_fraction", audit.get("capture_max_fraction"), "<"),
            ("technical_confounding", "technical_effect_size", audit.get("technical_effect_size_max"), "<"),
            ("doublet", "doublet_fraction", audit.get("doublet_fraction_max"), "<"),
            ("incompatible_coexpression", "incompatible_coexpression_fraction", audit.get("incompatible_program_fraction_min"), "<"),
        )
        for gate_name, column, threshold, operator in optional_checks:
            if threshold is None:
                continue
            value = float(row[column])
            gates[gate_name] = value >= float(threshold) if operator == ">=" else value < float(threshold)
        failed = [name for name, passed in gates.items() if not passed]
        accepted = not failed
        row.update(
            {
                "accepted_child_evidence_gate": accepted,
                "failed_gates": ";".join(failed),
                "recommended_candidate_action": "ACCEPT CHILD" if accepted else "REJECT CHILD",
                "confidence": "high" if accepted else "low",
                "biological_label_L1": "",
                "biological_label_L2": "",
                "biological_label_L3": "",
                "cl_id": "",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def decide_parents(audit: pd.DataFrame, children: pd.DataFrame, config: dict) -> pd.DataFrame:
    require_columns(
        audit,
        {
            "scope",
            "parent_id",
            "parent_uid",
            "audit_flagged",
            "doublet_fraction",
            "incompatible_program_fraction",
            "donor_dominated",
            "capture_dominated",
            "qc_outlier",
        },
        artifact="parent audit",
    )
    audit_thresholds = config["adaptive"]["audit"]
    rows = []
    for parent in audit.to_dict("records"):
        block = children[
            children["scope"].astype(str).eq(str(parent["scope"]))
            & children["parent_id"].astype(str).eq(str(parent["parent_id"]))
        ]
        accepted = block[block["accepted_child_evidence_gate"].map(as_bool)]
        accepted_ids = ";".join(accepted["candidate_child_uid"].astype(str))
        if not accepted.empty:
            outcome = "ACCEPT SPLIT"
            reasons = "one or more candidate children passed every evidence gate"
        elif (
            float(parent["doublet_fraction"]) >= float(audit_thresholds["doublet_fraction_max"])
            or float(parent["incompatible_program_fraction"])
            >= float(audit_thresholds.get("incompatible_program_fraction_min", 0.15))
        ):
            outcome = "DOUBLET/MIXED"
            reasons = "doublet enrichment or incompatible programs co-occur within parent cells"
        elif any(
            as_bool(parent[field])
            for field in ("donor_dominated", "capture_dominated", "qc_outlier")
        ):
            outcome = "TECHNICAL/UNRESOLVED"
            reasons = "donor, capture, QC, or complexity structure dominates"
        else:
            outcome = "RETAIN PARENT"
            reasons = "no candidate child passed the complete evidence gate"
        if outcome not in OUTCOMES:
            raise ContractError("parent decision is outside the frozen outcome set")
        rows.append(
            {
                "scope": str(parent["scope"]),
                "parent_id": str(parent["parent_id"]),
                "parent_uid": str(parent["parent_uid"]),
                "recommended_outcome": outcome,
                "decision_outcome": "",
                "accepted_child_ids_recommended": accepted_ids,
                "accepted_child_ids": "",
                "decision_reasons": reasons,
                "confidence": "high" if outcome in {"ACCEPT SPLIT", "DOUBLET/MIXED"} else "moderate",
                "manual_review_status": "PENDING_REVIEW",
                "biological_label_L1": "",
                "biological_label_L2": "",
                "biological_label_L3": "",
                "cl_id": "",
            }
        )
    result = pd.DataFrame(rows)
    if result.duplicated(["scope", "parent_id"]).any() or len(result) != len(audit):
        raise ContractError("every parent must have exactly one decision row")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--candidate-evidence", required=True, type=Path)
    parser.add_argument("--child-summary", required=True, type=Path)
    parser.add_argument("--parent-decisions", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        config = read_config(args.config)
        evidence = pd.read_csv(args.candidate_evidence, sep="\t", dtype=str, keep_default_na=False)
        audit = pd.read_csv(args.audit, sep="\t", dtype=str, keep_default_na=False)
        children = evaluate_candidates(config, evidence)
        decisions = decide_parents(audit, children, config)
        atomic_write_tsv(args.child_summary, children)
        atomic_write_tsv(args.parent_decisions, decisions)
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()

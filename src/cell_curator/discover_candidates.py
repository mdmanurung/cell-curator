#!/usr/bin/env python3
"""Discover persistent stored or bounded parent-local candidate children."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .contracts import (
    ContractError,
    as_bool,
    atomic_write_json,
    atomic_write_tsv,
    natural_key,
    read_config,
    require_columns,
    values_sha256,
)


@dataclass(frozen=True)
class Child:
    scope: str
    parent_id: str
    resolution: float
    child_id: str
    members: frozenset[str]
    parent_purity: float
    parent_fraction: float
    n_donors: int
    n_captures: int


def jaccard(left: frozenset[str], right: frozenset[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 0.0


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[b] = a


def _stored_children(config: dict, parents: pd.DataFrame, stored: pd.DataFrame) -> list[Child]:
    require_columns(
        stored,
        {"barcode", "scope", "resolution", "child_id"},
        artifact="stored memberships",
    )
    thresholds = config["adaptive"]["stored_candidates"]
    eligibility = config["adaptive"]["eligibility"]
    keys = config["input"]["keys"]
    donor_key = str(keys.get("donor", "")).strip()
    capture_key = str(keys.get("capture", "")).strip()
    indexed = parents.set_index("barcode", verify_integrity=True)
    rows: list[Child] = []
    for (scope, resolution, child_id), frame in stored.groupby(
        ["scope", "resolution", "child_id"], observed=True, sort=False
    ):
        members = pd.Index(frame["barcode"].astype(str))
        if set(members).difference(indexed.index):
            raise ContractError("stored memberships contain unknown barcodes")
        metadata = indexed.loc[members]
        parent_counts = metadata["parent_id"].astype(str).value_counts()
        parent_id = str(parent_counts.index[0])
        inside = metadata[metadata["parent_id"].astype(str).eq(parent_id)]
        parent_total = int(
            (
                (parents["scope"].astype(str) == str(scope))
                & (parents["parent_id"].astype(str) == parent_id)
            ).sum()
        )
        purity = len(inside) / len(metadata)
        fraction = len(inside) / parent_total
        donors = inside[donor_key].astype(str).nunique() if donor_key else 0
        captures = inside[capture_key].astype(str).nunique() if capture_key else 0
        donor_ok = eligibility.get("min_donors") is None or donors >= int(
            eligibility["min_donors"]
        )
        capture_ok = eligibility.get("min_captures") is None or captures >= int(
            eligibility["min_captures"]
        )
        if not (
            len(inside) >= int(eligibility["min_cells"])
            and donor_ok
            and capture_ok
            and purity >= float(thresholds["min_parent_purity"])
            and fraction >= float(thresholds["min_parent_fraction"])
            and fraction <= float(thresholds.get("max_parent_fraction", 0.95))
        ):
            continue
        rows.append(
            Child(
                scope=str(scope),
                parent_id=parent_id,
                resolution=float(resolution),
                child_id=str(child_id),
                members=frozenset(inside.index.astype(str)),
                parent_purity=purity,
                parent_fraction=fraction,
                n_donors=donors,
                n_captures=captures,
            )
        )
    return rows


def collapse_stored_families(config: dict, children: list[Child]) -> list[list[Child]]:
    threshold = float(config["adaptive"]["stored_candidates"]["adjacent_jaccard_min"])
    finder = UnionFind(len(children))
    by_scope_parent: dict[tuple[str, str], list[int]] = {}
    for index, child in enumerate(children):
        by_scope_parent.setdefault((child.scope, child.parent_id), []).append(index)
    for indices in by_scope_parent.values():
        resolutions = sorted({children[index].resolution for index in indices})
        for left_resolution, right_resolution in zip(
            resolutions[:-1], resolutions[1:], strict=True
        ):
            left = [index for index in indices if children[index].resolution == left_resolution]
            right = [
                index for index in indices if children[index].resolution == right_resolution
            ]
            for a in left:
                scores = [(jaccard(children[a].members, children[b].members), b) for b in right]
                if not scores:
                    continue
                score, b = max(scores)
                if score >= threshold:
                    finder.union(a, b)
    groups: dict[int, list[Child]] = {}
    for index, child in enumerate(children):
        groups.setdefault(finder.find(index), []).append(child)
    min_matches = int(config["adaptive"]["stored_candidates"]["min_adjacent_matches"])
    return [
        sorted(group, key=lambda item: item.resolution)
        for group in groups.values()
        if len({item.resolution for item in group}) >= min_matches + 1
    ]


def _candidate_rows_from_stored(
    config: dict, parents: pd.DataFrame, families: list[list[Child]]
) -> tuple[list[dict], list[dict], list[dict]]:
    registry: list[dict] = []
    membership: list[dict] = []
    ledger: list[dict] = []
    counter: dict[tuple[str, str], int] = {}
    for family_number, family in enumerate(
        sorted(
            families,
            key=lambda group: (
                group[0].scope,
                natural_key(group[0].parent_id),
                group[0].resolution,
            ),
        ),
        1,
    ):
        representative = family[0]
        key = (representative.scope, representative.parent_id)
        counter[key] = counter.get(key, 0) + 1
        candidate_key = f"{representative.scope}__p{representative.parent_id}__c{counter[key]}"
        parent_frame = parents[
            parents["scope"].astype(str).eq(representative.scope)
            & parents["parent_id"].astype(str).eq(representative.parent_id)
        ]
        closest = sorted(parent_frame["closest_sibling"].astype(str).unique())
        if len(closest) != 1 or not closest[0] or closest[0] == representative.parent_id:
            raise ContractError("closest sibling must be one different parent per candidate")
        registry.append(
            {
                "scope": representative.scope,
                "candidate_key": candidate_key,
                "parent_id": representative.parent_id,
                "parent_uid": f"{representative.scope}::{representative.parent_id}",
                "candidate_child_uid": f"{representative.scope}::{representative.parent_id}::candidate_{counter[key]}",
                "source_kind": "stored_resolution",
                "source_resolution": representative.resolution,
                "source_child_id": representative.child_id,
                "n_cells_in_parent": len(representative.members),
                "n_parent_remainder": len(parent_frame) - len(representative.members),
                "n_donors": representative.n_donors,
                "n_captures": representative.n_captures,
                "parent_purity": representative.parent_purity,
                "parent_fraction": representative.parent_fraction,
                "stability_pass": True,
                "eligibility_counts_pass": True,
                "closest_sibling": closest[0],
                "membership_sha256": values_sha256(sorted(representative.members)),
                "triage_action": "FULL_TESTING",
                "biological_label_L1": "",
                "biological_label_L2": "",
                "biological_label_L3": "",
                "cl_id": "",
            }
        )
        for barcode in parent_frame["barcode"].astype(str):
            membership.append(
                {
                    "scope": representative.scope,
                    "candidate_key": candidate_key,
                    "parent_id": representative.parent_id,
                    "barcode": barcode,
                    "candidate_membership": barcode in representative.members,
                }
            )
        ledger.append(
            {
                "family_id": f"stored_family_{family_number}",
                "scope": representative.scope,
                "parent_id": representative.parent_id,
                "n_source_resolutions": len({item.resolution for item in family}),
                "resolutions": ";".join(f"{item.resolution:g}" for item in family),
                "minimum_adjacent_jaccard": min(
                    jaccard(left.members, right.members)
                    for left, right in zip(family[:-1], family[1:], strict=True)
                ),
                "family_action": "FULL_TESTING",
                "candidate_key": candidate_key,
            }
        )
    return registry, membership, ledger


def _local_candidates(
    config: dict,
    parents: pd.DataFrame,
    local: pd.DataFrame,
    parents_with_stored: set[tuple[str, str]],
    counters: dict[tuple[str, str], int],
) -> tuple[list[dict], list[dict], dict[tuple[str, str], str]]:
    if local.empty:
        return [], [], {}
    require_columns(
        local,
        {
            "barcode",
            "scope",
            "parent_id",
            "resolution",
            "local_child_id",
            "stable",
            "biologically_resolving",
            "technical_dominated",
            "eligible",
        },
        artifact="parent-local candidates",
    )
    keys = config["input"]["keys"]
    donor_key = str(keys.get("donor", "")).strip()
    capture_key = str(keys.get("capture", "")).strip()
    registry: list[dict] = []
    membership: list[dict] = []
    evaluations: dict[tuple[str, str], str] = {}
    for (scope, parent_id), frame in local.groupby(
        ["scope", "parent_id"], observed=True, sort=False
    ):
        key = (str(scope), str(parent_id))
        if key in parents_with_stored:
            evaluations[key] = "stored candidate family already routed"
            continue
        chosen: tuple[float, pd.DataFrame] | None = None
        messages: list[str] = []
        for resolution in sorted(frame["resolution"].astype(float).unique()):
            partition = frame[frame["resolution"].astype(float).eq(resolution)].copy()
            clusters = partition["local_child_id"].astype(str).nunique()
            stable = partition["stable"].map(as_bool).all()
            resolving = partition["biologically_resolving"].map(as_bool).all()
            technical = partition["technical_dominated"].map(as_bool).any()
            eligible_clusters = (
                partition.assign(_eligible=partition["eligible"].map(as_bool))
                .groupby("local_child_id", observed=True)["_eligible"]
                .all()
                .sum()
            )
            valid = (
                stable
                and resolving
                and not technical
                and clusters >= 2
                and eligible_clusters >= 2
            )
            messages.append(
                f"{resolution:g}:n_clusters={clusters},eligible={eligible_clusters},valid={str(valid).lower()}"
            )
            if valid:
                chosen = (resolution, partition)
                break
        evaluations[key] = ";".join(messages)
        if chosen is None:
            continue
        resolution, partition = chosen
        groups = {
            str(child): set(block["barcode"].astype(str))
            for child, block in partition.groupby("local_child_id", observed=True)
        }
        remainder_id = max(groups, key=lambda child: (len(groups[child]), natural_key(child)))
        parent_frame = parents[
            parents["scope"].astype(str).eq(str(scope))
            & parents["parent_id"].astype(str).eq(str(parent_id))
        ]
        if set(parent_frame["barcode"].astype(str)) != set(partition["barcode"].astype(str)):
            raise ContractError("local partition does not cover its parent exactly")
        closest = sorted(parent_frame["closest_sibling"].astype(str).unique())
        if len(closest) != 1 or not closest[0] or closest[0] == str(parent_id):
            raise ContractError("closest sibling must be one different parent per candidate")
        for child_id in sorted(groups, key=natural_key):
            if child_id == remainder_id:
                continue
            counters[key] = counters.get(key, 0) + 1
            candidate_key = f"{scope}__p{parent_id}__c{counters[key]}"
            members = groups[child_id]
            selected = parent_frame[parent_frame["barcode"].astype(str).isin(members)]
            registry.append(
                {
                    "scope": str(scope),
                    "candidate_key": candidate_key,
                    "parent_id": str(parent_id),
                    "parent_uid": f"{scope}::{parent_id}",
                    "candidate_child_uid": f"{scope}::{parent_id}::candidate_{counters[key]}",
                    "source_kind": "parent_local",
                    "source_resolution": resolution,
                    "source_child_id": child_id,
                    "n_cells_in_parent": len(members),
                    "n_parent_remainder": len(parent_frame) - len(members),
                    "n_donors": selected[donor_key].astype(str).nunique() if donor_key else 0,
                    "n_captures": selected[capture_key].astype(str).nunique()
                    if capture_key
                    else 0,
                    "parent_purity": 1.0,
                    "parent_fraction": len(members) / len(parent_frame),
                    "stability_pass": True,
                    "eligibility_counts_pass": True,
                    "closest_sibling": closest[0],
                    "membership_sha256": values_sha256(sorted(members)),
                    "triage_action": "FULL_TESTING",
                    "biological_label_L1": "",
                    "biological_label_L2": "",
                    "biological_label_L3": "",
                    "cl_id": "",
                }
            )
            for barcode in parent_frame["barcode"].astype(str):
                membership.append(
                    {
                        "scope": str(scope),
                        "candidate_key": candidate_key,
                        "parent_id": str(parent_id),
                        "barcode": barcode,
                        "candidate_membership": barcode in members,
                    }
                )
    return registry, membership, evaluations


def run_discovery(
    *,
    config_path: Path,
    parents_path: Path,
    stored_path: Path,
    local_path: Path,
    audit_path: Path,
    output_dir: Path,
) -> dict:
    config = read_config(config_path)
    parents = pd.read_csv(parents_path, sep="\t", dtype=str, keep_default_na=False)
    keys = config["input"]["keys"]
    required_parent_columns = {"barcode", "scope", "parent_id", "closest_sibling"}
    required_parent_columns.update(
        str(keys[name]) for name in ("donor", "capture") if str(keys.get(name, "")).strip()
    )
    require_columns(
        parents,
        required_parent_columns,
        artifact="parent membership",
    )
    if parents["barcode"].duplicated().any():
        raise ContractError("parent membership must have one row per input cell")
    stored = pd.read_csv(stored_path, sep="\t", dtype=str, keep_default_na=False)
    stored["resolution"] = stored["resolution"].astype(float)
    local = (
        pd.read_csv(local_path, sep="\t", dtype=str, keep_default_na=False)
        if local_path.is_file() and local_path.stat().st_size
        else pd.DataFrame()
    )
    audit = pd.read_csv(audit_path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(audit, {"scope", "parent_id", "audit_flagged"}, artifact="parent audit")
    children = _stored_children(config, parents, stored)
    families = collapse_stored_families(config, children)
    registry, membership, ledger = _candidate_rows_from_stored(config, parents, families)
    parents_with_stored = {(row["scope"], row["parent_id"]) for row in registry}
    counters: dict[tuple[str, str], int] = {}
    for row in registry:
        key = (row["scope"], row["parent_id"])
        counters[key] = counters.get(key, 0) + 1
    local_registry, local_membership, evaluations = _local_candidates(
        config, parents, local, parents_with_stored, counters
    )
    registry.extend(local_registry)
    membership.extend(local_membership)
    registry_frame = pd.DataFrame(registry)
    membership_frame = pd.DataFrame(membership)
    ledger_frame = pd.DataFrame(ledger)
    candidate_parents = {(str(row["scope"]), str(row["parent_id"])) for row in registry}
    triage_rows = []
    for row in audit.itertuples(index=False):
        key = (str(row.scope), str(row.parent_id))
        flagged = as_bool(row.audit_flagged)
        if key in candidate_parents:
            action = "FULL_TESTING"
        elif flagged:
            action = "REJECT_NO_DISCRETE_CANDIDATE_EVIDENCE"
        else:
            action = "NO_ACTION_NOT_FLAGGED"
        triage_rows.append(
            {
                "scope": key[0],
                "parent_id": key[1],
                "audit_flagged": flagged,
                "triage_action": action,
                "local_resolution_evaluations": evaluations.get(key, ""),
            }
        )
    triage = pd.DataFrame(triage_rows)
    if not registry_frame.empty:
        registry_frame = registry_frame.sort_values(
            ["scope", "parent_id", "candidate_key"], ignore_index=True
        )
        membership_frame = membership_frame.sort_values(
            ["scope", "candidate_key", "barcode"], ignore_index=True
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_tsv(output_dir / "candidate_registry.tsv", registry_frame)
    atomic_write_tsv(output_dir / "candidate_membership.tsv", membership_frame)
    atomic_write_tsv(output_dir / "parent_triage.tsv", triage)
    atomic_write_tsv(output_dir / "stored_family_ledger.tsv", ledger_frame)
    manifest = {
        "schema_version": 3,
        "artifact_kind": "adaptive_candidate_discovery",
        "status": "complete",
        "created_utc": datetime.now(UTC).isoformat(),
        "run_id": config["run"]["run_id"],
        "n_candidates": len(registry_frame),
        "n_stored_families": len(ledger_frame),
        "n_local_candidates": len(local_registry),
        "invariants": {
            "every_flagged_parent_routed": not triage[triage["audit_flagged"]]["triage_action"]
            .eq("")
            .any(),
            "biological_labels_blank": True,
            "stable_one_cluster_rejected": True,
            "local_candidates_promoted_to_full_registry": True,
        },
    }
    atomic_write_json(output_dir / "discovery.manifest.json", manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--parents", required=True, type=Path)
    parser.add_argument("--stored-memberships", required=True, type=Path)
    parser.add_argument("--local-candidates", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        run_discovery(
            config_path=args.config,
            parents_path=args.parents,
            stored_path=args.stored_memberships,
            local_path=args.local_candidates,
            audit_path=args.audit,
            output_dir=args.output_dir,
        )
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()

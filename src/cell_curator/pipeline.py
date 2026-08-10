"""Public phase API and end-to-end independent annotation execution."""

from __future__ import annotations

import json
import os
import shutil
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .audit_parents import build_audit
from .benchmark import write_pending_status
from .config import CellCuratorConfig, load_config
from .contracts import (
    ContractError,
    require_safe_path_component,
    safe_path_component,
    sha256,
)
from .decisions import (
    build_comparator_manifest,
    build_mapping,
    decide_parents_from_artifacts,
    derive_candidate_evidence,
)
from .evidence import consolidate_lanes, lane_key, run_candidate_lane, run_lane
from .knowledge import collect_programs, programs_digest, providers_from_config
from .models import RunState
from .ontology import ontology_from_config, validate_program_ontology
from .provenance import atomic_publish, atomic_replace, publish_json, replace_json
from .refinement import discover_refinement_candidates
from .review import (
    build_critic_packet,
    build_evidence_cards,
    build_html_report,
    validate_reconciliation,
)


def run_root(config: CellCuratorConfig) -> Path:
    return Path(config.run.output_root) / config.run.run_id


def _state_path(root: Path) -> Path:
    return root / "run_state.json"


def _write_control_json(path: Path, value: dict[str, Any]) -> None:
    """Atomically replace the intentionally mutable run-state control document."""

    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.name}.partial.{os.getpid()}")
    try:
        partial.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def _load_state(root: Path) -> RunState:
    path = _state_path(root)
    if not path.is_file():
        raise ContractError("run state is missing; initialize the run first")
    return RunState.model_validate_json(path.read_text())


def _advance(root: Path, phase: str, *, blockers: list[str] | None = None) -> RunState:
    state = _load_state(root)
    completed = tuple(dict.fromkeys((*state.completed_phases, phase)))
    updated = state.model_copy(
        update={
            "phase": phase,
            "completed_phases": completed,
            "blockers": tuple(blockers if blockers is not None else state.blockers),
            "updated_utc": datetime.now(UTC).isoformat(),
        }
    )
    _write_control_json(_state_path(root), updated.model_dump(mode="json"))
    return updated


def _programs(config: CellCuratorConfig):
    programs = collect_programs(
        providers_from_config(
            config.knowledge.providers,
            allow_remote_enhancement=config.knowledge.allow_optional_remote_enhancement,
        )
    )
    ontology = ontology_from_config(config)
    if ontology is not None:
        validate_program_ontology(programs, ontology)
        return programs
    # Labels remain usable without an ontology, but unvalidated accessions must not
    # leak into proposals or write-back.
    return [
        program.model_copy(update={"cl_id": None, "parent_cl_id": None}) for program in programs
    ]


def _knowledge_provider_records(config: CellCuratorConfig) -> list[dict[str, Any]]:
    """Record observed local provider artifacts and verify declared checksums."""

    records: list[dict[str, Any]] = []
    for provider in config.knowledge.providers:
        record = provider.model_dump(mode="json")
        artifact_path = (
            Path(provider.cache_path)
            if provider.kind in {"remote_table", "remote_ontology"}
            else Path(provider.path)
        )
        ontology_status_path = run_root(config) / "knowledge" / "ontology_provider_status.json"
        ontology_status = (
            json.loads(ontology_status_path.read_text())
            if provider.kind == "remote_ontology" and ontology_status_path.is_file()
            else {}
        )
        if (
            provider.kind == "remote_ontology"
            and ontology_status.get("status") == "unavailable"
        ):
            record["artifact"] = {
                "status": "unavailable",
                "path": str(artifact_path),
                "reason": ontology_status.get("reason", "remote ontology validation failed"),
                "observed_sha256": sha256(artifact_path) if artifact_path.is_file() else "",
                "checksum_validated": False,
                "status_receipt": str(ontology_status_path),
            }
        elif artifact_path.is_file():
            observed = sha256(artifact_path)
            if provider.sha256 and observed != provider.sha256:
                raise ContractError(
                    f"{provider.kind} provider checksum mismatch for {artifact_path}: "
                    f"expected {provider.sha256}, observed {observed}"
                )
            record["artifact"] = {
                "status": "validated",
                "path": str(artifact_path.resolve()),
                "sha256": observed,
                "declared_checksum": bool(provider.sha256),
                "checksum_validated": bool(provider.sha256),
            }
        elif provider.kind == "remote_ontology":
            record["artifact"] = {
                "status": "unavailable",
                "path": str(artifact_path),
                "reason": ontology_status.get(
                    "reason", "optional remote ontology cache is not locally available"
                ),
                "checksum_validated": False,
                "status_receipt": (
                    str(ontology_status_path) if ontology_status_path.is_file() else ""
                ),
            }
        elif provider.kind == "remote_table" and not (
            config.knowledge.allow_optional_remote_enhancement
        ):
            record["artifact"] = {
                "status": "unavailable",
                "path": str(artifact_path),
                "reason": (
                    "optional remote enhancement is disabled and the checksum-pinned "
                    "cache is not locally available"
                ),
                "checksum_validated": False,
            }
        else:
            raise ContractError(
                f"configured {provider.kind} provider artifact does not exist: {artifact_path}"
            )
        records.append(record)
    return records


def _scope_filename(parent: str) -> str:
    return safe_path_component(parent)


def _replace_if_changed(path: Path, payload: bytes) -> str:
    """Replace a living artifact only when its bytes actually changed."""

    if path.is_file() and path.read_bytes() == payload:
        return sha256(path)
    return atomic_replace(path, payload)


def _hierarchy_level_number(level: str) -> int:
    try:
        number = int(str(level).removeprefix("L"))
    except ValueError as exc:
        raise ContractError(f"hierarchy level must use L<number>: {level}") from exc
    if number < 1:
        raise ContractError(f"hierarchy level must be positive: {level}")
    return number


def _validate_execution_vocabulary(
    vocabulary: pd.DataFrame,
) -> tuple[list[str], dict[str, str | None]]:
    from .hierarchy import validate_parent_vocabulary

    validate_parent_vocabulary(vocabulary)
    values = vocabulary.loc[:, ["level", "parent_scope", "label"]].astype(str)
    if values.duplicated().any():
        row = values.loc[values.duplicated(keep=False)].iloc[0]
        raise ContractError(
            "hierarchy vocabulary duplicates level/parent/label key "
            f"{tuple(row[column] for column in values.columns)}"
        )
    label_levels = values.groupby("label", observed=True)["level"].nunique()
    repeated = sorted(label_levels[label_levels > 1].index.astype(str))
    if repeated:
        raise ContractError(f"hierarchy labels cannot be reused across levels: {repeated}")
    levels = sorted(values["level"].unique(), key=_hierarchy_level_number)
    if not levels:
        raise ContractError("recursive hierarchy execution requires a non-empty vocabulary")
    first_parents = set(values.loc[values["level"].eq(levels[0]), "parent_scope"])
    if first_parents != {"ROOT"}:
        raise ContractError(
            f"the first hierarchy level must have exactly parent ROOT, observed {sorted(first_parents)}"
        )
    next_level: dict[str, str | None] = {}
    for index, level in enumerate(levels):
        following = levels[index + 1] if index + 1 < len(levels) else None
        next_level[level] = following
        if following is None:
            continue
        parents = set(values.loc[values["level"].eq(following), "parent_scope"])
        permitted = set(values.loc[values["level"].eq(level), "label"])
        invalid = sorted(parents.difference(permitted))
        if invalid:
            raise ContractError(
                f"strict hierarchy nesting failed: parents at {following} are not labels "
                f"from {level}: {invalid}"
            )
    return levels, next_level


def _normalize_score_table(
    value: pd.DataFrame | str | Path,
    *,
    level: str,
    parent: str,
) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        frame = value.copy(deep=True)
    else:
        source = Path(value)
        if not source.is_file():
            raise ContractError(f"hierarchy score table does not exist: {source}")
        separator = "\t" if source.suffix.casefold() in {".tsv", ".txt"} else ","
        frame = pd.read_csv(source, sep=separator, keep_default_na=False)
    required = {"cluster_id", "label", "score"}
    missing = required.difference(frame.columns)
    if missing:
        raise ContractError(
            f"hierarchy score table for {level}::{parent} misses columns: {sorted(missing)}"
        )
    frame["cluster_id"] = frame["cluster_id"].astype(str)
    frame["label"] = frame["label"].astype(str)
    if frame["cluster_id"].str.strip().eq("").any() or frame["label"].str.strip().eq("").any():
        raise ContractError(f"hierarchy score table for {level}::{parent} has blank keys")
    numeric_score = pd.to_numeric(frame["score"], errors="coerce")
    if not np.isfinite(numeric_score.to_numpy(dtype=float)).all():
        raise ContractError(
            f"hierarchy score table for {level}::{parent} has non-finite scores"
        )
    frame["score"] = numeric_score
    if frame.duplicated(["cluster_id", "label"]).any():
        raise ContractError(
            f"hierarchy score table for {level}::{parent} duplicates cluster/label rows"
        )
    return frame.sort_values(["cluster_id", "label"], kind="stable").reset_index(drop=True)


def _normalize_state_table(
    value: Mapping[str, Mapping[str, Any]] | str | Path,
    *,
    level: str,
    parent: str,
) -> dict[str, dict[str, Any]]:
    if isinstance(value, (str, Path)):
        source = Path(value)
        if not source.is_file():
            raise ContractError(f"hierarchy state table does not exist: {source}")
        loaded = json.loads(source.read_text())
    else:
        loaded = dict(value)
    if not isinstance(loaded, dict) or any(
        not isinstance(item, dict) for item in loaded.values()
    ):
        raise ContractError(
            f"hierarchy state table for {level}::{parent} must map clusters to objects"
        )
    return {str(cluster): dict(state) for cluster, state in sorted(loaded.items())}


def _valid_scope_checkpoint(
    path: Path,
    *,
    level: str,
    parent: str,
    state: RunState,
    assumption_digest: str,
    artifacts: Mapping[str, str],
) -> bool:
    if not path.is_file():
        return False
    try:
        checkpoint = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    if any(
        (
            checkpoint.get("status") != "complete",
            checkpoint.get("level") != level,
            checkpoint.get("parent") != parent,
            checkpoint.get("input_sha256") != state.input_sha256,
            checkpoint.get("config_sha256") != state.config_sha256,
            checkpoint.get("assumptions_sha256") != assumption_digest,
            checkpoint.get("artifacts") != dict(artifacts),
        )
    ):
        return False
    checksums = checkpoint.get("artifact_sha256", {})
    return all(
        Path(artifact).is_file() and checksums.get(name) == sha256(Path(artifact))
        for name, artifact in artifacts.items()
    )


def execute_recursive_hierarchy(
    config_path: str | Path,
    *,
    score_tables: Mapping[tuple[str, str], pd.DataFrame | str | Path],
    vocabulary_path: str | Path | None = None,
    state_tables: Mapping[tuple[str, str], Mapping[str, Mapping[str, Any]] | str | Path]
    | None = None,
    minimum_confidence: float | None = None,
    minimum_margin: float | None = None,
) -> dict[str, Any]:
    """Execute and resume arbitrary-depth, strictly parent-scoped annotation.

    Score tables are addressed only by their exact ``(level, parent_scope)`` key.
    Each executable scope gets independent input/label shards and a content-hashed
    checkpoint.  A reachable child without a score table is retained as its
    already-supported parent leaf rather than borrowing evidence from a sibling.
    """

    from .guided import (
        assumption_scope_sha256,
        checkpoint_level,
        load_state,
        require_label_gate,
        save_state,
    )
    from .hierarchy import (
        LABEL_COLUMNS,
        assign_parent_scoped_level,
        scaffold_labels_table,
        validate_hierarchy_nesting,
    )

    source = Path(config_path)
    config = load_config(source)
    root = run_root(config)
    state = load_state(root)
    if state.config_sha256 != sha256(source):
        raise ContractError(
            "recursive hierarchy configuration changed after scene setting; use a new run_id"
        )
    plan_path = root / "plan.json"
    if not plan_path.is_file():
        raise ContractError("recursive hierarchy execution requires Stage 2 plan.json")
    plan = json.loads(plan_path.read_text())
    selected_vocabulary = Path(vocabulary_path or plan.get("vocabulary_path", ""))
    if not selected_vocabulary.is_file():
        raise ContractError(
            f"recursive hierarchy vocabulary does not exist: {selected_vocabulary}"
        )
    vocabulary = pd.read_csv(selected_vocabulary, sep="\t", keep_default_na=False)
    levels, next_level = _validate_execution_vocabulary(vocabulary)
    available_scopes = {
        (str(row.level), str(row.parent_scope))
        for row in vocabulary[["level", "parent_scope"]]
        .drop_duplicates()
        .itertuples(index=False)
    }
    normalized_scores: dict[tuple[str, str], pd.DataFrame | str | Path] = {}
    for raw_scope, value in score_tables.items():
        if not isinstance(raw_scope, tuple) or len(raw_scope) != 2:
            raise ContractError("score table keys must be (level, parent_scope) tuples")
        scope = (str(raw_scope[0]), str(raw_scope[1]))
        if scope not in available_scopes:
            raise ContractError(f"score table scope is absent from the vocabulary: {scope}")
        if scope in normalized_scores:
            raise ContractError(f"score table scope is duplicated after normalization: {scope}")
        normalized_scores[scope] = value
    root_scope = (levels[0], "ROOT")
    if root_scope not in normalized_scores:
        raise ContractError(f"root hierarchy score table is required: {root_scope}")
    normalized_states: dict[tuple[str, str], Mapping[str, Mapping[str, Any]] | str | Path] = {}
    for raw_scope, value in (state_tables or {}).items():
        scope = (str(raw_scope[0]), str(raw_scope[1]))
        if scope not in normalized_scores:
            raise ContractError(f"state table has no matching score scope: {scope}")
        normalized_states[scope] = value

    pending: list[tuple[str, str]] = [root_scope]
    origins: dict[tuple[str, str], list[dict[str, str]]] = {}
    visited: set[tuple[str, str]] = set()
    outputs: dict[tuple[str, str], pd.DataFrame] = {}
    executed: list[str] = []
    resumed: list[str] = []
    retained: list[dict[str, str]] = []
    checkpoint_paths: dict[str, Path] = {}
    while pending:
        level, parent = pending.pop(0)
        scope = (level, parent)
        if scope in visited:
            continue
        visited.add(scope)
        state = require_label_gate(source, level=level, parent=parent)
        assumption_digest = assumption_scope_sha256(root, level=level, parent=parent)
        if scope not in normalized_scores:
            for origin in origins.get(scope, []):
                retained.append(
                    {
                        **origin,
                        "requested_level": level,
                        "requested_parent_scope": parent,
                        "leaf_label": parent,
                        "reason": "missing_child_score_table",
                    }
                )
            continue

        scores = _normalize_score_table(normalized_scores[scope], level=level, parent=parent)
        safe_parent = _scope_filename(parent)
        score_path = root / "hierarchy" / "scores" / level / f"{safe_parent}.tsv"
        _replace_if_changed(
            score_path,
            scores.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
        )
        artifacts: dict[str, str] = {
            "scores": str(score_path),
            "vocabulary": str(selected_vocabulary),
        }
        state_by_cluster: dict[str, dict[str, Any]] = {}
        if scope in normalized_states:
            state_by_cluster = _normalize_state_table(
                normalized_states[scope], level=level, parent=parent
            )
            state_path = root / "hierarchy" / "states" / level / f"{safe_parent}.json"
            state_payload = (
                json.dumps(state_by_cluster, indent=2, sort_keys=True, default=str) + "\n"
            ).encode()
            _replace_if_changed(state_path, state_payload)
            artifacts["states"] = str(state_path)
        labels_path = root / "labels" / "scopes" / level / f"{safe_parent}.tsv"
        artifacts["labels"] = str(labels_path)
        checkpoint_path = root / "checkpoints" / level / f"{safe_parent}.json"
        identity = f"{level}::{parent}"
        if _valid_scope_checkpoint(
            checkpoint_path,
            level=level,
            parent=parent,
            state=state,
            assumption_digest=assumption_digest,
            artifacts=artifacts,
        ):
            assigned = pd.read_csv(labels_path, sep="\t", keep_default_na=False)
            if set(LABEL_COLUMNS).difference(assigned.columns):
                raise ContractError(
                    f"resumed labels shard has an invalid schema: {labels_path}"
                )
            if set(assigned["level"].astype(str)) != {level} or set(
                assigned["parent_scope"].astype(str)
            ) != {parent}:
                raise ContractError(f"resumed labels shard violates its scope: {labels_path}")
            resumed.append(identity)
        else:
            assigned = assign_parent_scoped_level(
                scores,
                vocabulary,
                level=level,
                parent_scope=parent,
                minimum_confidence=(
                    config.evidence.min_confidence
                    if minimum_confidence is None
                    else minimum_confidence
                ),
                minimum_margin=(
                    config.evidence.min_margin if minimum_margin is None else minimum_margin
                ),
                state_by_cluster=state_by_cluster,
            )
            _replace_if_changed(
                labels_path,
                assigned.to_csv(
                    sep="\t", index=False, columns=list(LABEL_COLUMNS), lineterminator="\n"
                ).encode(),
            )
            checkpoint_level(
                source,
                level=level,
                parent=parent,
                artifacts=artifacts,
                next_action="continue to the next approved parent-scoped hierarchy scope",
            )
            executed.append(identity)
        outputs[scope] = assigned.loc[:, list(LABEL_COLUMNS)].copy()
        checkpoint_paths[identity] = checkpoint_path

        following = next_level[level]
        available_parents = (
            set(
                vocabulary.loc[
                    vocabulary["level"].astype(str).eq(following), "parent_scope"
                ].astype(str)
            )
            if following is not None
            else set()
        )
        child_scopes: set[tuple[str, str]] = set()
        for row in assigned.to_dict("records"):
            label = str(row["label"])
            origin = {
                "source_level": level,
                "source_parent_scope": parent,
                "source_cluster_id": str(row["cluster_id"]),
            }
            if label == "Unknown":
                retained.append(
                    {
                        **origin,
                        "requested_level": level,
                        "requested_parent_scope": parent,
                        "leaf_label": "Unknown",
                        "reason": "unknown_assignment",
                    }
                )
            elif following is not None and label in available_parents:
                child_scope = (following, label)
                origins.setdefault(child_scope, []).append(origin)
                child_scopes.add(child_scope)
            else:
                retained.append(
                    {
                        **origin,
                        "requested_level": level,
                        "requested_parent_scope": parent,
                        "leaf_label": label,
                        "reason": "terminal_vocabulary_leaf",
                    }
                )
        pending.extend(sorted(child_scopes, key=lambda item: item[1]))

    level_tables: dict[str, dict[str, Any]] = {}
    consolidated: list[pd.DataFrame] = []
    for level in levels:
        frames = [frame for (item_level, _), frame in outputs.items() if item_level == level]
        if not frames:
            continue
        proposals = pd.concat(frames, ignore_index=True)
        labels_path = root / "labels" / f"{level}.tsv"
        durable = scaffold_labels_table(labels_path, proposals, retain_unmentioned=False)
        consolidated.append(durable)
        level_tables[level] = {
            "path": str(labels_path),
            "sha256": sha256(labels_path),
            "parent_scopes": sorted(set(durable["parent_scope"].astype(str))),
            "n_rows": len(durable),
        }
    hierarchy_validation = validate_hierarchy_nesting(consolidated)
    retained_columns = [
        "source_level",
        "source_parent_scope",
        "source_cluster_id",
        "requested_level",
        "requested_parent_scope",
        "leaf_label",
        "reason",
    ]
    retained_frame = pd.DataFrame(retained, columns=retained_columns).drop_duplicates()
    if not retained_frame.empty:
        retained_frame = retained_frame.sort_values(retained_columns, kind="stable")
    retained_path = root / "hierarchy" / "retained_leaves.tsv"
    _replace_if_changed(
        retained_path,
        retained_frame.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    manifest = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_recursive_hierarchy_execution",
        "status": "complete",
        "run_id": config.run.run_id,
        "input_sha256": state.input_sha256,
        "config_sha256": state.config_sha256,
        "vocabulary_path": str(selected_vocabulary),
        "vocabulary_sha256": sha256(selected_vocabulary),
        "executed_scopes": executed,
        "resumed_scopes": resumed,
        "reachable_scopes": sorted(f"{level}::{parent}" for level, parent in visited),
        "unused_score_scopes": sorted(
            f"{level}::{parent}" for level, parent in set(normalized_scores).difference(visited)
        ),
        "checkpoints": {
            identity: {"path": str(path), "sha256": sha256(path)}
            for identity, path in sorted(checkpoint_paths.items())
        },
        "level_tables": level_tables,
        "retained_leaves": {
            "path": str(retained_path),
            "sha256": sha256(retained_path),
            "n_rows": len(retained_frame),
            "missing_child_score_tables": int(
                retained_frame["reason"].eq("missing_child_score_table").sum()
            ),
        },
        "hierarchy_validation": hierarchy_validation,
    }
    manifest_path = root / "hierarchy" / "execution.manifest.json"
    replace_json(manifest_path, manifest)
    latest_state = load_state(root)
    progress = dict(latest_state.progress)
    progress.update(
        {
            "completed_levels": list(level_tables),
            "recursive_execution_manifest": str(manifest_path),
            "retained_leaf_count": len(retained_frame),
        }
    )
    save_state(
        root,
        latest_state.model_copy(
            update={
                "stage": "execute",
                "phase": "HIERARCHY_EXECUTED",
                "current_level": levels[-1],
                "current_parent": "",
                "next_action": "review, reconcile, and approve write-back",
                "progress": progress,
            }
        ),
    )
    return manifest


def execute_recursive_hierarchy_manifest(
    config_path: str | Path,
    manifest_path: str | Path,
) -> dict[str, Any]:
    """Execute a recursive hierarchy from a strict, portable JSON manifest.

    Relative artifact paths are resolved from the manifest directory.  The
    manifest is only an invocation contract; the executor copies and hashes the
    actual score, state, vocabulary, label, and checkpoint artifacts in the run
    directory.
    """

    source = Path(manifest_path)
    if not source.is_file():
        raise ContractError(f"recursive hierarchy manifest does not exist: {source}")
    try:
        value = json.loads(source.read_text())
    except json.JSONDecodeError as exc:
        raise ContractError(f"recursive hierarchy manifest is invalid JSON: {source}") from exc
    if not isinstance(value, dict):
        raise ContractError("recursive hierarchy manifest must be a JSON object")
    allowed = {
        "schema_version",
        "vocabulary_path",
        "score_tables",
        "minimum_confidence",
        "minimum_margin",
    }
    extra = set(value).difference(allowed)
    if extra:
        raise ContractError(
            f"recursive hierarchy manifest has unsupported fields: {sorted(extra)}"
        )
    if value.get("schema_version") != 3:
        raise ContractError("recursive hierarchy manifest requires schema_version 3")
    entries = value.get("score_tables")
    if not isinstance(entries, list) or not entries:
        raise ContractError("recursive hierarchy manifest requires non-empty score_tables")

    def artifact_path(raw: object, *, field: str) -> Path:
        if not isinstance(raw, str) or not raw.strip():
            raise ContractError(f"recursive hierarchy manifest requires nonblank {field}")
        path = Path(raw)
        return path if path.is_absolute() else source.parent / path

    score_tables: dict[tuple[str, str], Path] = {}
    state_tables: dict[tuple[str, str], Path] = {}
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ContractError(f"score_tables[{index}] must be a JSON object")
        entry_allowed = {"level", "parent_scope", "path", "state_path"}
        entry_extra = set(entry).difference(entry_allowed)
        if entry_extra:
            raise ContractError(
                f"score_tables[{index}] has unsupported fields: {sorted(entry_extra)}"
            )
        level = entry.get("level")
        parent = entry.get("parent_scope")
        if not isinstance(level, str) or not level.strip():
            raise ContractError(f"score_tables[{index}] requires nonblank level")
        if not isinstance(parent, str) or not parent.strip():
            raise ContractError(f"score_tables[{index}] requires nonblank parent_scope")
        scope = (level, parent)
        if scope in score_tables:
            raise ContractError(f"recursive hierarchy manifest duplicates scope: {scope}")
        score_tables[scope] = artifact_path(entry.get("path"), field="score table path")
        if "state_path" in entry:
            state_tables[scope] = artifact_path(
                entry.get("state_path"), field="state table path"
            )

    confidence = value.get("minimum_confidence")
    margin = value.get("minimum_margin")
    if confidence is not None and (
        isinstance(confidence, bool)
        or not isinstance(confidence, (int, float))
        or not 0 <= float(confidence) <= 1
    ):
        raise ContractError("minimum_confidence must be a number between 0 and 1")
    if margin is not None and (
        isinstance(margin, bool) or not isinstance(margin, (int, float)) or float(margin) < 0
    ):
        raise ContractError("minimum_margin must be a nonnegative number")
    vocabulary = value.get("vocabulary_path")
    return execute_recursive_hierarchy(
        config_path,
        score_tables=score_tables,
        vocabulary_path=(
            artifact_path(vocabulary, field="vocabulary_path")
            if vocabulary is not None
            else None
        ),
        state_tables=state_tables,
        minimum_confidence=float(confidence) if confidence is not None else None,
        minimum_margin=float(margin) if margin is not None else None,
    )


def initialize(config_path: Path) -> dict[str, Any]:
    """Execute Stage 1 and persist local knowledge without opening the label gate."""

    from .guided import load_state, prepare_scene

    config = load_config(config_path)
    root = run_root(config)
    state_path = root / "run_state.json"
    contract_path = root / "input_contract.manifest.json"
    if state_path.is_file() and contract_path.is_file():
        state = load_state(root)
        if state.config_sha256 != sha256(config_path):
            raise ContractError(
                "existing run state belongs to a different configuration; choose a new run_id"
            )
        if state.input_sha256 != sha256(Path(config.input.path)):
            raise ContractError(
                "existing run state belongs to a changed input; choose a new run_id"
            )
    else:
        prepare_scene(config_path)
    manifest = json.loads(contract_path.read_text())
    programs = _programs(config)
    ontology = ontology_from_config(config)
    publish_json(
        root / "knowledge.manifest.json",
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_local_knowledge",
            "run_id": config.run.run_id,
            "n_programs": len(programs),
            "program_digest": programs_digest(programs),
            "providers": _knowledge_provider_records(config),
            "external_service_required": False,
            "ontology": (
                {
                    "status": "validated",
                    "source": ontology.source,
                    "n_terms": len(ontology.terms),
                    "sha256": sha256(Path(ontology.source)),
                }
                if ontology is not None
                else {
                    "status": "unavailable",
                    "limitation": "labels are permitted without CL IDs",
                }
            ),
        },
    )
    return manifest


def audit(config_path: Path) -> pd.DataFrame:
    config = load_config(config_path)
    root = run_root(config)
    cells = pd.read_csv(root / "prepared" / "cells.tsv", sep="\t", keep_default_na=False)
    frame = build_audit(config.model_dump(mode="json"), cells)
    path = root / "adaptive_parent_audit.tsv"
    atomic_publish(path, frame.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    technical_columns = [
        column
        for column in frame.columns
        if column
        in {
            "schema_version",
            "scope",
            "parent_id",
            "parent_uid",
            "n_cells",
            "n_donors",
            "n_captures",
            "max_donor_fraction",
            "max_capture_fraction",
            "max_sex_fraction",
            "donor_dominated",
            "capture_dominated",
            "qc_outlier",
        }
        or column.startswith("median_technical__")
    ]
    doublet_columns = [
        "schema_version",
        "scope",
        "parent_id",
        "parent_uid",
        "n_cells",
        "doublet_fraction",
        "doublet_enriched",
        "incompatible_program_fraction",
        "incompatible_program_enriched",
    ]
    adaptive_dir = root / "adaptive"
    atomic_publish(
        adaptive_dir / "technical_summary.tsv",
        frame[technical_columns].to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    atomic_publish(
        adaptive_dir / "doublet_summary.tsv",
        frame[doublet_columns].to_csv(sep="\t", index=False, lineterminator="\n").encode(),
    )
    _advance(root, "PARENTS_AUDITED")
    return frame


def refine(config_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    config = load_config(config_path)
    result = discover_refinement_candidates(config=config, run_root=run_root(config))
    build_comparator_manifest(config=config, run_root=run_root(config))
    from .hierarchy import profile_refined_clusters

    profile_refined_clusters(config_path)
    _advance(run_root(config), "CANDIDATES_DISCOVERED")
    return result


def assign_hierarchy_level(
    config_path: str | Path,
    *,
    scores_path: str | Path,
    level: str,
    parent: str,
    vocabulary_path: str | Path | None = None,
    state_path: str | Path | None = None,
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    """Assign one parent-scoped hierarchy level and checkpoint the labels table.

    Shared by the CLI and the Python API. The level is validated as a safe path
    component here rather than at the CLI boundary, because it is interpolated
    into the output filename — a caller reaching this directly must not be able
    to skip that check.
    """

    from .guided import checkpoint_level, require_label_gate
    from .hierarchy import assign_parent_scoped_level, scaffold_labels_table

    level = require_safe_path_component(level, field="hierarchy level")
    config = load_config(config_path)
    root = run_root(config)
    require_label_gate(config_path, level=level, parent=parent)
    vocabulary = pd.read_csv(
        Path(vocabulary_path)
        if vocabulary_path is not None
        else root / "plan" / "vocabulary.tsv",
        sep="\t",
        keep_default_na=False,
    )
    scores = pd.read_csv(Path(scores_path), sep="\t", keep_default_na=False)
    state_by_cluster = (
        json.loads(Path(state_path).read_text()) if state_path is not None else None
    )
    assigned = assign_parent_scoped_level(
        scores,
        vocabulary,
        level=level,
        parent_scope=parent,
        minimum_confidence=config.evidence.min_confidence,
        minimum_margin=config.evidence.min_margin,
        state_by_cluster=state_by_cluster,
    )
    output = (
        Path(output_path) if output_path is not None else root / "labels" / f"{level}.tsv"
    )
    labels = scaffold_labels_table(output, assigned)
    checkpoint_level(
        config_path,
        level=level,
        parent=parent,
        artifacts={"labels": str(output), "scores": str(Path(scores_path))},
        next_action="review and reconcile this hierarchy level",
    )
    return {"assigned": assigned, "labels": labels, "output": output}


def evidence_lane(
    config_path: Path,
    *,
    assay: str,
    representation: str,
) -> Path:
    config = load_config(config_path)
    programs = _programs(config)
    result = run_lane(
        config=config,
        config_path=config_path,
        run_root=run_root(config),
        assay_id=assay,
        representation_name=representation,
        programs=programs,
    )
    return result.receipt_path


def candidate_evidence_lane(config_path: Path, *, comparison_lane_id: str) -> Path:
    config = load_config(config_path)
    result = run_candidate_lane(
        config=config,
        config_path=config_path,
        run_root=run_root(config),
        comparison_lane_id=comparison_lane_id,
        programs=_programs(config),
    )
    return result.receipt_path


def evidence_all(config_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    config = load_config(config_path)
    root = run_root(config)
    receipts = []
    for lane in config.evidence.lanes:
        receipt = evidence_lane(
            config_path,
            assay=lane.assay,
            representation=lane.representation,
        )
        receipts.append(
            {
                "lane_key": lane_key(lane.assay, lane.representation, lane.lane_key),
                "path": str(receipt),
                "sha256": sha256(receipt),
            }
        )
    publish_json(
        root / "evidence" / "receipts.validated.json",
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_validated_backend_receipts",
            "run_id": config.run.run_id,
            "status": "complete",
            "n_applicable_lanes": len(receipts),
            "zero_lane_run": len(receipts) == 0,
            "receipts": receipts,
            "all_gpu_fallback_flags_false": True,
        },
    )
    from .guided import require_assumption_gate

    level = config.guidance.hierarchy_levels[0]
    require_assumption_gate(config_path, level=level, parent="ROOT")
    result = consolidate_lanes(config=config, run_root=root)
    from .hierarchy import labels_from_existing_proposals, scaffold_labels_table

    labels = labels_from_existing_proposals(result[1], level=level, parent_scope="ROOT")
    scaffold_labels_table(root / "labels" / f"{level}.tsv", labels)
    _advance(root, "EVIDENCE_COMPLETE")
    return result


def decide(config_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    config = load_config(config_path)
    root = run_root(config)
    from .guided import require_assumption_gate

    require_assumption_gate(
        config_path, level=config.guidance.hierarchy_levels[0], parent="ROOT"
    )
    comparator_manifest = pd.read_csv(
        root / "comparators" / "comparator_manifest.tsv", sep="\t", keep_default_na=False
    )
    for comparison_lane_id in comparator_manifest["lane_id"].astype(str):
        candidate_evidence_lane(config_path, comparison_lane_id=comparison_lane_id)
    children = derive_candidate_evidence(
        config=config,
        run_root=root,
        programs=_programs(config),
    )
    parents = decide_parents_from_artifacts(config=config, run_root=root)
    _advance(root, "DECISIONS_PROPOSED")
    return children, parents


def map_cells(config_path: Path) -> pd.DataFrame:
    config = load_config(config_path)
    root = run_root(config)
    from .guided import require_assumption_gate

    require_assumption_gate(
        config_path, level=config.guidance.hierarchy_levels[0], parent="ROOT"
    )
    mapping = build_mapping(config=config, run_root=root)
    _advance(root, "MAPPING_VALIDATED")
    return mapping


def review(config_path: Path, *, strict: bool = True) -> Path:
    config = load_config(config_path)
    root = run_root(config)
    from .guided import checkpoint_level, require_assumption_gate

    level = config.guidance.hierarchy_levels[0]
    require_assumption_gate(config_path, level=level, parent="ROOT")
    build_evidence_cards(config, root)
    build_critic_packet(config, root)
    level_report = build_html_report(config, root, strict=strict, level=level)
    report = build_html_report(config, root, strict=strict)
    _advance(root, "REVIEW_PACKET_READY")
    checkpoint_level(
        config_path,
        level=level,
        parent="ROOT",
        artifacts={
            "labels": str(root / "labels" / f"{level}.tsv"),
            "level_report": str(level_report),
            "full_report": str(report),
            "critic_packet": str(root / "review" / "critic_packet.json"),
        },
        next_action=(
            "pause for level approval before descending"
            if config.guidance.interactive and config.guidance.pause_after_each_level
            else "continue to the next planned hierarchy scope or final validation"
        ),
    )
    return report


def validate(config_path: Path) -> dict[str, Any]:
    config = load_config(config_path)
    root = run_root(config)
    required = [
        "context.json",
        "environment.json",
        "assumptions.json",
        "plan.json",
        "input_contract.manifest.json",
        "qc/qc.manifest.json",
        "qc/cell_metrics.tsv",
        "adaptive_parent_audit.tsv",
        "adaptive/technical_summary.tsv",
        "adaptive/doublet_summary.tsv",
        "discovery/discovery.manifest.json",
        "comparators/comparator_manifest.json",
        "evidence/receipts.validated.json",
        "evidence/summary/evidence.manifest.json",
        "evidence/summary/program_cooccurrence.tsv",
        "evidence/summary/source_corroboration.tsv",
        "decisions/candidate_evidence.manifest.json",
        "decisions/parent_decision_summary.tsv",
        "mapping/mapping.validation.json",
        "review/critic_packet.json",
        "review/colors.json",
        "review/report.manifest.json",
        "review/report.html",
    ]
    missing = [item for item in required if not (root / item).is_file()]
    if missing:
        raise ContractError(f"run validation misses artifacts: {missing}")
    mapping_validation = json.loads((root / "mapping" / "mapping.validation.json").read_text())
    if mapping_validation.get("status") != "complete":
        raise ContractError("mapping validation is incomplete")
    from .data import compute_qc
    from .guided import assumption_scope_sha256, load_assumptions, load_state
    from .hierarchy import validate_hierarchy_nesting
    from .review import check_report_completeness

    compute_qc(config, run_root=root)
    state = load_state(root)
    assumptions = load_assumptions(root)
    pending_assumptions = [
        item["assumption_id"]
        for item in assumptions
        if bool(item.get("material", True)) and item.get("status") != "APPROVED"
    ]
    if (
        pending_assumptions
        or not state.assumptions_approved
        or state.approved_assumptions_sha256 != sha256(root / "assumptions.json")
    ):
        raise ContractError(
            "run validation requires the current material assumptions gate; "
            f"pending={pending_assumptions}"
        )
    labels_paths = sorted((root / "labels").glob("*.tsv"))
    if not labels_paths:
        raise ContractError("run validation requires durable hierarchy labels tables")
    label_frames = [pd.read_csv(path, sep="\t", keep_default_na=False) for path in labels_paths]
    mapping = pd.read_csv(
        root / "mapping" / "parent_child_mapping.tsv", sep="\t", keep_default_na=False
    )
    hierarchy_validation = validate_hierarchy_nesting(
        label_frames,
        cell_membership=mapping[["barcode", "leaf_id"]],
    )
    check_report_completeness(config, root, require_report=True)
    first_level = config.guidance.hierarchy_levels[0]
    first_checkpoint = root / "checkpoints" / first_level / "ROOT.json"
    if not first_checkpoint.is_file():
        raise ContractError(
            f"run validation requires a hierarchy checkpoint: {first_checkpoint}"
        )
    checkpoint = json.loads(first_checkpoint.read_text())
    expected_checkpoint_assumptions = assumption_scope_sha256(
        root, level=first_level, parent="ROOT"
    )
    if (
        checkpoint.get("status") != "complete"
        or checkpoint.get("assumptions_sha256") != expected_checkpoint_assumptions
    ):
        raise ContractError("hierarchy checkpoint is incomplete or uses stale assumptions")
    for name, artifact_path in checkpoint.get("artifacts", {}).items():
        artifact = Path(str(artifact_path))
        if not artifact.is_file() or sha256(artifact) != checkpoint.get(
            "artifact_sha256", {}
        ).get(name):
            raise ContractError(f"hierarchy checkpoint artifact is absent or stale: {name}")
    plan = json.loads((root / "plan.json").read_text())
    for route in plan.get("routes", []):
        if not route.get("assessments") or any(
            not assessment.get("reasons") for assessment in route["assessments"]
        ):
            raise ContractError("method routing omitted unavailable/unsuitable reasons")
    reconciliation = validate_reconciliation(root)
    ontology = ontology_from_config(config)
    proposals = pd.read_csv(
        root / "evidence" / "summary" / "annotation_proposals.tsv",
        sep="\t",
        keep_default_na=False,
    )
    supplied_cl_ids = proposals[
        proposals.get("cl_id", pd.Series(index=proposals.index, dtype=str))
        .astype(str)
        .str.strip()
        .ne("")
    ]
    if not supplied_cl_ids.empty and ontology is None:
        raise ContractError("CL IDs are present but no local ontology was configured")
    if ontology is not None:
        for proposal in supplied_cl_ids.to_dict("records"):
            ontology.validate_term(
                str(proposal["cl_id"]), expected_name=str(proposal["cell_type"])
            )
    source_sha = sha256(Path(config.input.path))
    input_manifest = json.loads((root / "input_contract.manifest.json").read_text())
    if source_sha != input_manifest["input_sha256"]:
        raise ContractError("source input changed after the run was frozen")
    manifest = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_run_validation",
        "run_id": config.run.run_id,
        "status": "complete",
        "promotion_status": "computationally_complete_pending_human_review",
        "source_unchanged": True,
        "critics_reconciled": reconciliation["all_reconciled"],
        "ontology_status": "validated" if ontology is not None else "unavailable-cl-free",
        "mapping_invariants": {
            **mapping_validation["invariants"],
            **hierarchy_validation["invariants"],
        },
        "artifacts": {
            item: {"sha256": sha256(root / item), "bytes": (root / item).stat().st_size}
            for item in required
        },
        "writeback_permitted": False,
    }
    publish_json(root / "final_provenance.manifest.json", manifest)
    _advance(root, "COMPUTATIONALLY_COMPLETE")
    return manifest


def build_review_packet(config_path: Path) -> Path:
    config = load_config(config_path)
    root = run_root(config)
    packet = root / "review_packet"
    if packet.exists():
        manifest_path = packet / "review_packet.manifest.json"
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text())
            files = manifest.get("files")
            if (
                manifest.get("status") != "complete"
                or manifest.get("run_id") != config.run.run_id
                or manifest.get("immutable") is not True
                or not isinstance(files, dict)
            ):
                raise ContractError("existing review packet manifest is invalid")
            observed = {
                str(path.relative_to(packet))
                for path in packet.rglob("*")
                if path.is_file() and path != manifest_path
            }
            if observed != set(files):
                raise ContractError("existing review packet file inventory drifted")
            for relative, record in files.items():
                frozen = packet / relative
                source = root / relative
                if (
                    not source.is_file()
                    or sha256(frozen) != record.get("sha256")
                    or frozen.stat().st_size != record.get("bytes")
                    or sha256(source) != record.get("sha256")
                ):
                    raise ContractError(
                        f"existing review packet or source artifact drifted: {relative}"
                    )
            return packet
        try:
            # Snakemake creates the parent of a declared file output before the
            # rule starts. Accept only that empty directory; never replace a
            # partially populated packet.
            packet.rmdir()
        except OSError as exc:
            raise ContractError("incomplete review packet directory already exists") from exc
    validate(config_path)
    stage = packet.with_name(f".{packet.name}.partial.{os.getpid()}")
    stage.mkdir(parents=True)
    include = [
        "context.json",
        "context.md",
        "environment.json",
        "assumptions.json",
        "assumptions.md",
        "plan.json",
        "plan.md",
        "run_state.json",
        "input_contract.manifest.json",
        "knowledge.manifest.json",
        "qc",
        "adaptive_parent_audit.tsv",
        "discovery",
        "evidence",
        "decisions",
        "labels",
        "checkpoints",
        "mapping",
        "review",
        "final_provenance.manifest.json",
    ]
    try:
        for relative in include:
            source = root / relative
            target = stage / relative
            if source.is_dir():
                shutil.copytree(source, target)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
        files = sorted(path for path in stage.rglob("*") if path.is_file())
        manifest = {
            "schema_version": 3,
            "artifact_kind": "cell_curator_immutable_review_packet",
            "run_id": config.run.run_id,
            "status": "complete",
            "immutable": True,
            "files": {
                str(path.relative_to(stage)): {
                    "sha256": sha256(path),
                    "bytes": path.stat().st_size,
                }
                for path in files
            },
        }
        (stage / "review_packet.manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )
        os.replace(stage, packet)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return packet


def run_annotation(config_path: str | Path) -> Path:
    """Execute all computational phases and return the immutable review packet."""

    path = Path(config_path)
    from .guided import construct_plan, load_state, require_assumption_gate

    config = load_config(path)
    if len(config.guidance.hierarchy_levels) != 1:
        raise ContractError(
            "the standard pipeline supports exactly one hierarchy level; "
            "use `cell-curator run execute --score-manifest MANIFEST.json` "
            "for recursive multi-level annotation"
        )
    existing_packet = run_root(config) / "review_packet" / "review_packet.manifest.json"
    if existing_packet.is_file():
        return build_review_packet(path)
    initialize(path)
    state = load_state(run_root(config))
    if not state.plan_complete:
        if config.guidance.interactive:
            raise ContractError(
                "interactive run is paused after scene setting; run `cell-curator run plan`, "
                "review/approve assumptions, then retry execute"
            )
        construct_plan(path)
    require_assumption_gate(path, level=config.guidance.hierarchy_levels[0], parent="ROOT")
    audit(path)
    refine(path)
    evidence_all(path)
    decide(path)
    map_cells(path)
    review(path, strict=True)
    validate(path)
    benchmarks = run_root(config) / "benchmarks"
    for modality in ("rna", "spatial", "cite-seq", "atac", "multiome"):
        write_pending_status(
            benchmarks / f"public_{modality}.pending.json",
            dataset=f"public-{modality}-benchmark",
            modality=modality,
            reason="Public dataset and external comparator exports were not supplied to this local run.",
            reproduce=[
                "cell-curator benchmark evaluate --predictions <neutral.tsv> --system <name> --dataset <dataset> --output <result.json>"
            ],
        )
    return build_review_packet(path)

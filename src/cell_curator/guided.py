"""Three-stage guided workflow state, assumptions, planning, and human gates."""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from .capabilities import annotation_capable, build_capability_registry
from .config import BiologicalContextConfig, CellCuratorConfig, load_config
from .contracts import (
    ContractError,
    require_safe_path_component,
    safe_path_component,
    sha256,
)
from .data import canonicalize, inspect
from .knowledge import collect_programs, programs_digest, providers_from_config
from .models import RunMode, RunState
from .provenance import atomic_publish, atomic_replace, publish_json, replace_json


def run_directory(config: CellCuratorConfig) -> Path:
    return Path(config.run.output_root) / config.run.run_id


def _state_path(root: Path) -> Path:
    return root / "run_state.json"


def load_state(root: Path) -> RunState:
    path = _state_path(root)
    if not path.is_file():
        raise ContractError(
            "run state is missing; call prepare_annotation_scene(config_path) "
            "(CLI: `cell-curator run scene`) first"
        )
    return RunState.model_validate_json(path.read_text())


def save_state(root: Path, state: RunState) -> RunState:
    updated = state.model_copy(update={"updated_utc": datetime.now(UTC).isoformat()})
    replace_json(_state_path(root), updated.model_dump(mode="json"))
    return updated


def _context_unresolved(config: CellCuratorConfig, inspection: dict[str, Any]) -> list[str]:
    unresolved = [f"context.{item}" for item in config.guidance.context.unresolved()]
    if not config.input.keys.embedding:
        unresolved.append("input.keys.embedding")
    if not (
        config.guidance.context.visualization_key.strip() or config.input.keys.embedding.strip()
    ):
        unresolved.append("guidance.context.visualization_key")
    categorical_value = inspection.get("categorical_covariates", {})
    categorical = (
        set(map(str, categorical_value))
        if isinstance(categorical_value, dict)
        else {
            str(item.get("column", item)) if isinstance(item, dict) else str(item)
            for item in categorical_value
        }
    )
    for column in config.guidance.context.composition_covariates:
        if column not in categorical:
            unresolved.append(f"composition_covariate:{column}")
    if inspection.get("blockers"):
        unresolved.extend(f"input:{item}" for item in inspection["blockers"])
    return sorted(set(unresolved))


def _render_context(value: dict[str, Any]) -> str:
    context = value["biological_context"]
    lines = [
        "# cell-curator scene context",
        "",
        f"- Run: `{value['run_id']}`",
        f"- Input: `{value['input']['path']}`",
        f"- Input SHA-256: `{value['input']['sha256']}`",
        f"- Cluster key: `{value['contracts']['cluster_key']}`",
        f"- Anchor embedding: `{value['contracts']['embedding_key']}`",
        f"- Visualization: `{value['contracts']['visualization_key']}`",
        f"- Spatial coordinates: `{value['contracts']['spatial_key'] or 'not configured'}`",
        f"- Organism: {context['organism'] or 'UNRESOLVED'}",
        f"- Tissue: {context['tissue'] or 'UNRESOLVED'}",
        f"- Developmental stage: {context['developmental_stage'] or 'UNRESOLVED'}",
        f"- Condition: {context['condition'] or 'UNRESOLVED'}",
        f"- Technology: {context['technology']}",
        f"- Experimental context: {context['experimental_context'] or 'UNRESOLVED'}",
        f"- Composition covariates: {', '.join(context['composition_covariates']) or 'none selected'}",
        f"- Vocabulary contract: {context['vocabulary_contract'] or 'UNRESOLVED'}",
        f"- Reference contract: {context['reference_contract'] or 'UNRESOLVED'}",
        f"- Write-back prefix: `{context['writeback_prefix']}`",
        f"- Compute profile: `{value['compute']['compute_profile']}`",
        f"- Context complete: `{value['context_complete']}`",
        f"- Unresolved: {'; '.join(value['unresolved']) if value['unresolved'] else 'none'}",
        "",
        "Automated annotations are hypotheses. Missing panel features are unmeasured. "
        "Cell type and state remain separate, and Unknown is a successful outcome.",
    ]
    return "\n".join(lines) + "\n"


def _assumption_id(scope: str, statement: str) -> str:
    return hashlib.sha256(f"{scope}\0{statement}".encode()).hexdigest()[:16]


def _ledger_path(root: Path) -> Path:
    return root / "assumptions.json"


def load_assumptions(root: Path) -> list[dict[str, Any]]:
    path = _ledger_path(root)
    if not path.is_file():
        return []
    value = json.loads(path.read_text())
    if not isinstance(value, dict) or not isinstance(value.get("assumptions"), list):
        raise ContractError("assumptions ledger does not satisfy its structured schema")
    return list(value["assumptions"])


def _write_assumptions(root: Path, rows: list[dict[str, Any]]) -> str:
    ordered = sorted(rows, key=lambda item: (item["scope"], item["assumption_id"]))
    payload = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_assumptions_ledger",
        "assumptions": ordered,
    }
    replace_json(_ledger_path(root), payload)
    digest = sha256(_ledger_path(root))
    lines = ["# Annotation assumptions ledger", "", f"Digest: `{digest}`", ""]
    if not ordered:
        lines.append("No assumptions have been recorded.")
    for row in ordered:
        lines.extend(
            [
                f"## {row['assumption_id']} · {row['scope']}",
                "",
                f"- Status: `{row['status']}`",
                f"- Material: `{row['material']}`",
                f"- Statement: {row['statement']}",
                f"- Reviewer: {row.get('reviewer') or 'pending'}",
                f"- Approval mode: {row.get('approval_mode') or 'pending'}",
                "",
            ]
        )
    atomic_replace(root / "assumptions.md", ("\n".join(lines) + "\n").encode())
    return digest


def add_assumptions(
    root: Path,
    values: list[tuple[str, str, bool]],
) -> tuple[list[dict[str, Any]], str]:
    existing = {item["assumption_id"]: item for item in load_assumptions(root)}
    for scope, statement, material in values:
        statement = statement.strip()
        if not statement:
            continue
        identifier = _assumption_id(scope, statement)
        if identifier not in existing:
            existing[identifier] = {
                "assumption_id": identifier,
                "scope": scope,
                "statement": statement,
                "material": bool(material),
                "status": "PENDING",
                "reviewer": "",
                "approval_mode": "",
                "approved_utc": "",
            }
    rows = list(existing.values())
    return rows, _write_assumptions(root, rows)


def list_assumptions(config_path: str | Path) -> dict[str, Any]:
    """Return the living assumptions ledger for a configured run."""

    return {"assumptions": load_assumptions(run_directory(load_config(config_path)))}


def record_assumption(
    config_path: str | Path,
    *,
    scope: str,
    statement: str,
    material: bool = True,
) -> dict[str, Any]:
    """Add a living assumption and revoke approval when it is materially new."""

    config = load_config(config_path)
    root = run_directory(config)
    state = load_state(root)
    before = {_item["assumption_id"] for _item in load_assumptions(root)}
    rows, digest = add_assumptions(root, [(scope, statement, material)])
    after = {_item["assumption_id"] for _item in rows}
    added = sorted(after.difference(before))
    if added and material:
        state = state.model_copy(
            update={
                "assumptions_sha256": digest,
                "assumptions_approved": False,
                "approved_assumptions_sha256": "",
                "labels_permitted": False,
                "next_action": "approve the materially changed assumptions ledger",
            }
        )
    else:
        state = state.model_copy(
            update={
                "assumptions_sha256": digest,
                "approved_assumptions_sha256": (
                    digest if state.assumptions_approved else state.approved_assumptions_sha256
                ),
            }
        )
    save_state(root, state)
    return {
        "schema_version": 3,
        "artifact_kind": "cell_curator_assumption_update",
        "status": "complete",
        "added": added,
        "material": material,
        "approval_revoked": bool(added and material),
        "assumptions_sha256": digest,
    }


def prepare_scene(config_path: str | Path) -> dict[str, Any]:
    """Inspect, detect, freeze, and persist Stage 1 without producing labels."""

    source = Path(config_path)
    config = load_config(source)
    root = run_directory(config)
    root.mkdir(parents=True, exist_ok=True)
    inspection = inspect(config)
    from .environment import detect_environment

    compute_report = detect_environment(
        configured_cpu_resources=config.hpc.cpu_resources,
        configured_gpu_resources=config.hpc.gpu_resources,
    )
    compute = compute_report.model_dump(mode="json")
    manifest = canonicalize(config, config_path=source, run_root=root)
    unresolved = _context_unresolved(config, inspection)
    context = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_scene_context",
        "run_id": config.run.run_id,
        "input": {
            "path": manifest["input_path"],
            "sha256": manifest["input_sha256"],
            "kind": config.input.kind,
            "n_cells": manifest["n_cells"],
            "assays": inspection["assays"],
        },
        "contracts": {
            "cluster_key": config.input.keys.canonical_parent,
            "scope_key": config.input.keys.scope,
            "embedding_key": config.input.keys.embedding,
            "visualization_key": config.guidance.context.visualization_key
            or config.input.keys.embedding,
            "spatial_key": config.input.keys.spatial,
            "stored_resolutions": config.input.keys.stored_resolutions,
        },
        "biological_context": config.guidance.context.model_dump(mode="json"),
        "compute": compute,
        "interactive": config.guidance.interactive,
        "preauthorized": config.guidance.preauthorized,
        "unresolved": unresolved,
        "context_complete": not unresolved,
    }
    replace_json(root / "context.json", context)
    atomic_replace(root / "context.md", _render_context(context).encode())
    assumptions = [("global", item, True) for item in config.guidance.assumptions] + [
        (
            "global",
            "Automated predictions are hypotheses and disagreement remains visible.",
            True,
        ),
        (
            "global",
            "Panel-absent features are unmeasured rather than negative evidence.",
            True,
        ),
        ("global", "Unknown is a permitted successful outcome.", True),
        ("global", "Cell type and cell state are assigned on separate axes.", True),
    ]
    _, ledger_digest = add_assumptions(root, assumptions)
    capabilities = build_capability_registry(config)
    blockers = list(unresolved)
    fatal_blockers = list(inspection.get("blockers", []))
    blockers.extend(fatal_blockers)
    if not annotation_capable(capabilities):
        blockers.append(
            "no annotation-capable evidence lane is executable; biological calls must remain Unknown"
        )
    effective_mode = (
        RunMode.EVIDENCE_ONLY
        if fatal_blockers or config.run.mode is RunMode.EVIDENCE_ONLY
        else RunMode.ANNOTATION
    )
    if _state_path(root).is_file():
        prior = load_state(root)
        if prior.input_sha256 != manifest["input_sha256"] or prior.config_sha256 != sha256(
            source
        ):
            raise ContractError(
                "existing run state belongs to a different input/configuration; choose a new run_id"
            )
        state = prior.model_copy(
            update={
                "mode": effective_mode,
                "stage": "scene" if not prior.plan_complete else prior.stage,
                "phase": "SCENE_SET",
                "blockers": tuple(sorted(set(blockers))),
                "context_complete": not unresolved,
                "compute_profile": compute_report.compute_profile.value,
                "compute": compute,
                "interactive": config.guidance.interactive,
                "preauthorized": config.guidance.preauthorized,
                "assumptions_sha256": ledger_digest,
                "labels_permitted": (
                    prior.assumptions_approved
                    and not blockers
                    and effective_mode is RunMode.ANNOTATION
                ),
                "next_action": (
                    "resolve scene context blockers"
                    if blockers
                    else "construct and review the annotation plan"
                ),
            }
        )
    else:
        state = RunState(
            run_id=config.run.run_id,
            mode=effective_mode,
            phase="SCENE_SET",
            input_path=str(Path(config.input.path).resolve()),
            input_sha256=manifest["input_sha256"],
            config_sha256=sha256(source),
            canonical_membership_sha256=manifest["canonical_membership_sha256"],
            blockers=tuple(sorted(set(blockers))),
            completed_phases=("INPUT_FROZEN", "SCENE_SET"),
            labels_permitted=False,
            stage="scene",
            next_action=(
                "resolve scene context blockers"
                if blockers
                else "construct and review the annotation plan"
            ),
            compute_profile=compute_report.compute_profile.value,
            compute=compute,
            context_complete=not unresolved,
            interactive=config.guidance.interactive,
            preauthorized=config.guidance.preauthorized,
            assumptions_sha256=ledger_digest,
            progress={"completed_levels": [], "checkpoints": []},
        )
    save_state(root, state)
    replace_json(root / "environment.json", compute)
    return context


def confirm_context(
    config_path: str | Path,
    *,
    input_path: str | Path,
    reviewer: str,
) -> dict[str, Any]:
    """Persist an interactive confirm-or-correct response for inferred scene facts."""

    if not reviewer.strip():
        raise ContractError("context confirmation requires a named reviewer")
    config = load_config(config_path)
    root = run_directory(config)
    context_path = root / "context.json"
    if not context_path.is_file():
        prepare_scene(config_path)
    context = json.loads(context_path.read_text())
    supplied = json.loads(Path(input_path).read_text())
    if not isinstance(supplied, dict):
        raise ContractError("context confirmation must be a JSON object")
    biological_value = supplied.get("biological_context", supplied)
    try:
        biological = BiologicalContextConfig.model_validate(biological_value)
    except ValueError as exc:
        raise ContractError(f"invalid biological context confirmation: {exc}") from exc
    unresolved = [f"context.{item}" for item in biological.unresolved()]
    if unresolved:
        raise ContractError("context confirmation remains incomplete: " + ", ".join(unresolved))
    inspection = inspect(config)
    categorical_value = inspection.get("categorical_covariates", {})
    categorical = (
        set(map(str, categorical_value))
        if isinstance(categorical_value, dict)
        else {
            str(item.get("column", item)) if isinstance(item, dict) else str(item)
            for item in categorical_value
        }
    )
    invalid_covariates = sorted(set(biological.composition_covariates).difference(categorical))
    if invalid_covariates:
        raise ContractError(
            f"confirmed composition covariates are absent/non-categorical: {invalid_covariates}"
        )
    other_unresolved = [
        item
        for item in context.get("unresolved", [])
        if not str(item).startswith(("context.", "composition_covariate:"))
    ]
    context["biological_context"] = biological.model_dump(mode="json")
    context["unresolved"] = other_unresolved
    context["context_complete"] = not other_unresolved
    context["confirmation"] = {
        "reviewer": reviewer.strip(),
        "confirmed_utc": datetime.now(UTC).isoformat(),
        "mode": "interactive",
    }
    replace_json(context_path, context)
    atomic_replace(root / "context.md", _render_context(context).encode())
    state = load_state(root)
    decisions = (*state.decisions, context["confirmation"] | {"kind": "context_confirmation"})
    remaining_blockers = tuple(
        item
        for item in state.blockers
        if not str(item).startswith(("context.", "composition_covariate:"))
    )
    state = state.model_copy(
        update={
            "context_complete": not other_unresolved,
            "blockers": remaining_blockers,
            "decisions": decisions,
            "next_action": (
                "construct and review the annotation plan"
                if not other_unresolved
                else "resolve remaining input-contract blockers"
            ),
        }
    )
    save_state(root, state)
    return context


def _vocabulary(programs: list[Any]) -> tuple[pd.DataFrame, list[str]]:
    by_label: dict[str, set[str]] = {}
    for program in programs:
        by_label.setdefault(program.label, set()).add(program.parent or "ROOT")
    duplicated = {label: parents for label, parents in by_label.items() if len(parents) > 1}
    if duplicated:
        raise ContractError(
            "hierarchy labels must have exactly one parent scope: "
            + "; ".join(
                f"{label} -> {sorted(parents)}" for label, parents in duplicated.items()
            )
        )
    labels = set(by_label)
    unresolved = sorted(
        {parent for parents in by_label.values() for parent in parents if parent != "ROOT"}
        - labels
    )
    depth: dict[str, int] = {}

    def resolve_depth(label: str, trail: tuple[str, ...] = ()) -> int:
        if label in depth:
            return depth[label]
        if label in trail:
            raise ContractError(
                f"marker-program hierarchy contains a cycle: {trail + (label,)}"
            )
        parent = next(iter(by_label[label]))
        value = (
            1
            if parent == "ROOT" or parent not in by_label
            else resolve_depth(parent, trail + (label,)) + 1
        )
        depth[label] = value
        return value

    rows = []
    for program in programs:
        rows.append(
            {
                "level": f"L{resolve_depth(program.label)}",
                "parent_scope": program.parent or "ROOT",
                "label": program.label,
                "axis": program.axis,
                "cl_id": program.cl_id or "",
                "parent_cl_id": program.parent_cl_id or "",
                "positive_markers": "|".join(program.positive),
                "negative_markers": "|".join(program.negative),
                "source": program.source,
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["level", "parent_scope", "axis", "label"], kind="stable"
    ), unresolved


def construct_plan(config_path: str | Path) -> dict[str, Any]:
    """Build Stage 2 vocabulary, routes, assumptions, and unresolved decisions."""

    source = Path(config_path)
    config = load_config(source)
    root = run_directory(config)
    if not (root / "context.json").is_file():
        prepare_scene(source)
    state = load_state(root)
    if not state.context_complete:
        raise ContractError(
            "cannot plan with unresolved scene context; inspect "
            "annotation_progress(config_path) (CLI: `cell-curator run status`)"
        )
    programs = collect_programs(
        providers_from_config(
            config.knowledge.providers,
            allow_remote_enhancement=config.knowledge.allow_optional_remote_enhancement,
        )
    )
    vocabulary, missing_parents = (
        _vocabulary(programs)
        if programs
        else (
            pd.DataFrame(
                columns=[
                    "level",
                    "parent_scope",
                    "label",
                    "axis",
                    "cl_id",
                    "parent_cl_id",
                    "positive_markers",
                    "negative_markers",
                    "source",
                ]
            ),
            [],
        )
    )
    vocab_path = root / "plan" / "vocabulary.tsv"
    atomic_publish(vocab_path, vocabulary.to_csv(sep="\t", index=False).encode())
    scene = json.loads((root / "context.json").read_text())
    biological_context = scene["biological_context"]
    from .data import close_input, load_input
    from .environment import EnvironmentReport
    from .routing import (
        HealthContext,
        MethodResource,
        RoutingRequest,
        Technology,
        route_crosscheck,
    )

    loaded = load_input(config, backed=True)
    try:
        obs_columns = list(map(str, loaded.root.obs.columns))
    finally:
        close_input(loaded)
    existing = sorted(
        column
        for column in obs_columns
        if column.endswith("_pred") or "predicted" in column.casefold()
    )
    references = [item for item in config.knowledge.providers if item.kind == "reference"]

    existing_resources = tuple(
        MethodResource(name=column, version="input-column") for column in existing
    )
    reference_artifacts: list[dict[str, Any]] = []
    reference_resources_list: list[MethodResource] = []
    for item in references:
        reference_path = Path(item.path)
        if not reference_path.is_file():
            raise ContractError(
                f"configured labeled reference does not exist: {reference_path}"
            )
        observed_reference_sha256 = sha256(reference_path)
        if item.sha256 and observed_reference_sha256 != item.sha256:
            raise ContractError(
                f"reference provider checksum mismatch for {reference_path}: "
                f"expected {item.sha256}, observed {observed_reference_sha256}"
            )
        reference_artifacts.append(
            {
                "path": str(reference_path.resolve()),
                "version": item.version,
                "sha256": observed_reference_sha256,
                "declared_checksum": bool(item.sha256),
                "checksum_validated": bool(item.sha256),
            }
        )
        reference_resources_list.append(
            MethodResource(
                name=str(reference_path.resolve()),
                version=item.version,
                organisms=(str(biological_context["organism"]),),
                checksum_validated=bool(item.sha256),
            )
        )
    reference_resources = tuple(reference_resources_list)
    # Package installation alone does not prove that a compatible model or
    # bounded reference artifact exists.  These routes remain unavailable until
    # an actual discovered resource is supplied to the router/API.
    celltypist_resources: tuple[MethodResource, ...] = ()
    model_resources: tuple[MethodResource, ...] = ()
    scimilarity_resources: tuple[MethodResource, ...] = ()
    census_resources: tuple[MethodResource, ...] = ()
    condition = str(biological_context["condition"]).casefold()
    if any(token in condition for token in ("tumor", "cancer", "malignan")):
        health_context = HealthContext.TUMOR
    elif any(token in condition for token in ("healthy", "normal", "control")):
        health_context = HealthContext.HEALTHY
    elif condition:
        health_context = HealthContext.DISEASE
    else:
        health_context = HealthContext.UNKNOWN
    environment = EnvironmentReport.model_validate(state.compute)
    technology = Technology(str(biological_context["technology"]))

    def level_key(value: object) -> tuple[int, str]:
        rendered = str(value)
        try:
            return int(rendered.removeprefix("L")), rendered
        except ValueError:
            return 10**9, rendered

    levels = (
        sorted(map(str, vocabulary["level"].unique()), key=level_key)
        if not vocabulary.empty
        else list(config.guidance.hierarchy_levels)
    )
    route_scopes = (
        sorted(
            {
                (str(row.level), str(row.parent_scope))
                for row in vocabulary[["level", "parent_scope"]]
                .drop_duplicates()
                .itertuples(index=False)
            },
            key=lambda item: (*level_key(item[0]), item[1]),
        )
        if not vocabulary.empty
        else [(str(level), "ROOT") for level in levels]
    )
    routes = []
    for level, parent_scope in route_scopes:
        request = RoutingRequest.from_environment(
            environment=environment,
            level=str(level),
            parent_scope=parent_scope,
            organism=str(biological_context["organism"]),
            tissue=str(biological_context["tissue"]),
            developmental_stage=str(biological_context["developmental_stage"]),
            health_context=health_context,
            technology=technology,
            allow_network=config.knowledge.allow_optional_remote_enhancement,
            existing_predictions=existing_resources,
            reference_models=model_resources,
            labeled_references=reference_resources,
            celltypist_models=celltypist_resources,
            azimuth_references=(
                reference_resources
                if str(biological_context["organism"]).casefold() in {"human", "homo sapiens"}
                else ()
            ),
            scimilarity_models=scimilarity_resources,
            census_references=census_resources,
        )
        routes.append(route_crosscheck(request).model_dump(mode="json"))
    unresolved = [f"vocabulary parent is undefined: {item}" for item in missing_parents]
    if config.run.mode is RunMode.ANNOTATION and vocabulary.empty:
        unresolved.append("no local marker/reference vocabulary is available")
    assumption_values: list[tuple[str, str, bool]] = []
    for row in vocabulary.to_dict("records"):
        scope = f"{row['level']}::{row['parent_scope']}"
        assumption_values.append(
            (
                scope,
                f"{row['label']} is a permitted child of {row['parent_scope']} with positive markers "
                f"{row['positive_markers'] or 'none'} and negative markers {row['negative_markers'] or 'none'}.",
                True,
            )
        )
    _, digest = add_assumptions(root, assumption_values)
    plan = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_annotation_plan",
        "run_id": config.run.run_id,
        "program_digest": programs_digest(programs),
        "vocabulary_path": str(vocab_path),
        "levels": levels,
        "routes": routes,
        "reference_artifacts": reference_artifacts,
        "existing_prediction_columns": existing,
        "unresolved": unresolved,
        "executable": not unresolved,
        "assumptions_sha256": digest,
    }
    publish_json(root / "plan.json", plan)
    lines = [
        "# Executable annotation plan",
        "",
        f"- Levels: {', '.join(levels) or 'none'}",
        f"- Vocabulary rows: {len(vocabulary)}",
        f"- Assumptions digest: `{digest}`",
        f"- Unresolved: {'; '.join(unresolved) if unresolved else 'none'}",
        "",
        "## Primary cross-check routes",
        "",
    ]
    for route in routes:
        primary = route.get("primary_method") or "unavailable"
        lines.append(
            f"- {route['level']}::{route['parent_scope']}: `{primary}` — {route['next_action']}"
        )
    lines.extend(
        [
            "",
            "## Execution",
            "",
            "For each parent scope: score, cross-check, assign or Unknown, review, reconcile, "
            "checkpoint, then descend only into children of that parent.",
        ]
    )
    atomic_publish(root / "plan.md", ("\n".join(lines) + "\n").encode())
    state = state.model_copy(
        update={
            "stage": "plan",
            "phase": "PLAN_READY" if not unresolved else "PLAN_BLOCKED",
            "plan_complete": not unresolved,
            "current_level": levels[0] if levels else "",
            "current_parent": "ROOT",
            "assumptions_sha256": digest,
            "blockers": tuple(sorted(set((*state.blockers, *unresolved)))),
            "next_action": (
                "resolve plan assumptions"
                if unresolved
                else "approve the living assumptions ledger"
            ),
        }
    )
    save_state(root, state)
    if not config.guidance.interactive and config.guidance.preauthorized and not unresolved:
        approve_assumptions(
            source,
            reviewer=config.guidance.reviewer,
            scope="*",
            preauthorized=True,
        )
        plan["gate"] = "auto-approved-recorded"
    return plan


def approve_assumptions(
    config_path: str | Path,
    *,
    reviewer: str,
    scope: str = "*",
    preauthorized: bool = False,
) -> dict[str, Any]:
    """Record a human or pre-authorized approval without removing the gate."""

    if not reviewer.strip():
        raise ContractError("assumption approval requires a named reviewer")
    config = load_config(config_path)
    root = run_directory(config)
    state = load_state(root)
    if not state.plan_complete:
        raise ContractError("assumptions cannot be approved before an executable plan exists")
    if preauthorized and (config.guidance.interactive or not config.guidance.preauthorized):
        raise ContractError("pre-authorization was not declared in the non-interactive config")
    rows = load_assumptions(root)
    matched = 0
    now = datetime.now(UTC).isoformat()
    for row in rows:
        if scope == "*" or row["scope"] in {"global", scope}:
            row.update(
                {
                    "status": "APPROVED",
                    "reviewer": reviewer.strip(),
                    "approval_mode": (
                        "preauthorized_non_interactive" if preauthorized else "interactive"
                    ),
                    "approved_utc": now,
                }
            )
            matched += 1
    if not matched:
        raise ContractError(f"no assumptions match approval scope {scope!r}")
    digest = _write_assumptions(root, rows)
    all_material = all(
        row["status"] == "APPROVED" for row in rows if bool(row.get("material", True))
    )
    scopes = tuple(
        sorted(
            set(
                (*state.signed_off_scopes, scope),
            )
        )
    )
    state = state.model_copy(
        update={
            "assumptions_approved": all_material,
            "approved_assumptions_sha256": digest if all_material else "",
            "assumptions_sha256": digest,
            "signed_off_scopes": scopes,
            "labels_permitted": all_material and state.mode is RunMode.ANNOTATION,
            "next_action": (
                "execute the first hierarchy level"
                if all_material
                else "approve remaining assumptions"
            ),
        }
    )
    save_state(root, state)
    return {
        "schema_version": 3,
        "artifact_kind": "cell_curator_assumption_gate",
        "status": "complete" if all_material else "partial",
        "scope": scope,
        "reviewer": reviewer.strip(),
        "approval_mode": "preauthorized_non_interactive" if preauthorized else "interactive",
        "approved_count": matched,
        "all_material_approved": all_material,
        "assumptions_sha256": digest,
    }


def require_assumption_gate(
    config_path: str | Path,
    *,
    level: str = "",
    parent: str = "ROOT",
) -> RunState:
    """Require complete scene/plan state and approved applicable assumptions."""

    config = load_config(config_path)
    root = run_directory(config)
    state = load_state(root)
    if not state.context_complete or not state.plan_complete:
        raise ContractError("execution requires complete scene and plan stages")
    ledger_path = _ledger_path(root)
    if not ledger_path.is_file() or sha256(ledger_path) != state.assumptions_sha256:
        raise ContractError(
            "assumptions ledger and run state disagree; rebuild or approve the current ledger"
        )
    rows = load_assumptions(root)
    scopes = {"global"}
    if level:
        scopes.add(f"{level}::{parent}")
    pending = [
        row["assumption_id"]
        for row in rows
        if bool(row.get("material", True))
        and row["scope"] in scopes
        and row["status"] != "APPROVED"
    ]
    if pending:
        raise ContractError(
            f"Gate A blocks biological labels until assumptions are approved; pending={pending}"
        )
    if not config.guidance.interactive and not (
        config.guidance.preauthorized and state.preauthorized
    ):
        raise ContractError("non-interactive execution lacks recorded pre-authorization")
    return state


def require_label_gate(
    config_path: str | Path,
    *,
    level: str = "",
    parent: str = "ROOT",
) -> RunState:
    """Require the assumptions gate and annotation mode before emitting labels."""

    state = require_assumption_gate(config_path, level=level, parent=parent)
    if state.mode is not RunMode.ANNOTATION:
        raise ContractError("Gate A blocks biological labels while the run is evidence-only")
    return state


def assumption_scope_sha256(
    root: str | Path,
    *,
    level: str = "",
    parent: str = "ROOT",
    require_approved: bool = True,
) -> str:
    """Hash only the material assumptions applicable to one hierarchy scope."""

    scopes = {"global"}
    if level:
        scopes.add(f"{level}::{parent}")
    applicable = [
        {
            "assumption_id": str(row["assumption_id"]),
            "scope": str(row["scope"]),
            "statement": str(row["statement"]),
            "material": bool(row.get("material", True)),
            "status": str(row.get("status", "PENDING")),
        }
        for row in load_assumptions(Path(root))
        if bool(row.get("material", True)) and str(row.get("scope")) in scopes
    ]
    applicable.sort(key=lambda row: (row["scope"], row["assumption_id"]))
    pending = [row["assumption_id"] for row in applicable if row["status"] != "APPROVED"]
    if require_approved and pending:
        raise ContractError(
            f"Gate A scope digest cannot be issued with pending assumptions: {pending}"
        )
    payload = json.dumps(
        {"schema_version": 3, "scopes": sorted(scopes), "assumptions": applicable},
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def checkpoint_level(
    config_path: str | Path,
    *,
    level: str,
    parent: str,
    artifacts: dict[str, str],
    next_action: str,
) -> dict[str, Any]:
    level = require_safe_path_component(level, field="hierarchy level")
    config = load_config(config_path)
    root = run_directory(config)
    state = require_assumption_gate(config_path, level=level, parent=parent)
    missing = sorted(name for name, path in artifacts.items() if not Path(path).is_file())
    if missing:
        raise ContractError(
            f"cannot checkpoint {level}::{parent}; declared artifacts are absent: {missing}"
        )
    checksums = {name: sha256(Path(path)) for name, path in artifacts.items()}
    scope_digest = assumption_scope_sha256(root, level=level, parent=parent)
    checkpoint = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_hierarchy_checkpoint",
        "level": level,
        "parent": parent,
        "input_sha256": state.input_sha256,
        "config_sha256": state.config_sha256,
        "assumptions_sha256": scope_digest,
        "assumption_scopes": ["global", f"{level}::{parent}"],
        "artifacts": artifacts,
        "artifact_sha256": checksums,
        "status": "complete",
    }
    safe_parent = safe_path_component(parent)
    path = root / "checkpoints" / level / f"{safe_parent}.json"
    replace_json(path, checkpoint)
    progress = dict(state.progress)
    completed = list(progress.get("checkpoints", []))
    identity = f"{level}::{parent}"
    if identity not in completed:
        completed.append(identity)
    progress["checkpoints"] = completed
    state = state.model_copy(
        update={
            "stage": "execute",
            "phase": "LEVEL_CHECKPOINTED",
            "current_level": level,
            "current_parent": parent,
            "next_action": next_action,
            "progress": progress,
        }
    )
    save_state(root, state)
    return checkpoint


def status(config_path: str | Path) -> dict[str, Any]:
    config = load_config(config_path)
    root = run_directory(config)
    state = load_state(root)
    return {
        **state.model_dump(mode="json"),
        "current_position": {
            "stage": state.stage,
            "level": state.current_level or None,
            "parent": state.current_parent or None,
        },
        "blockers": list(state.blockers),
        "next_action": state.next_action,
        "artifacts": {
            "context": str(root / "context.json"),
            "assumptions": str(root / "assumptions.json"),
            "plan": str(root / "plan.json"),
        },
    }

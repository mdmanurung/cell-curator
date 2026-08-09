"""Evidence cards, critic packets, reconciliation, and self-contained HTML reports."""

from __future__ import annotations

import html
import json
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import CellCuratorConfig
from .contracts import (
    ContractError,
    require_safe_path_component,
    safe_path_component,
    sha256,
)
from .data import close_input, dense, load_input, resolve_obsm
from .guided import load_assumptions, load_state
from .provenance import atomic_publish, atomic_replace, replace_json

CRITIC_COLUMNS = (
    "finding_id",
    "target_kind",
    "target_id",
    "critic_pass",
    "severity",
    "claim",
    "evidence_path",
    "proposed_resolution",
)
RECONCILIATION_COLUMNS = (
    "finding_id",
    "status",
    "resolution",
    "reviewer",
)
TERMINAL_RECONCILIATION_STATES = frozenset({"RESOLVED", "DISMISSED"})
RECONCILIATION_STATES = TERMINAL_RECONCILIATION_STATES | {"OPEN"}
CRITIC_SEVERITIES = frozenset({"BLOCKER", "ERROR", "WARNING", "INFO"})


def build_evidence_cards(config: CellCuratorConfig, run_root: Path) -> list[Path]:
    labels_paths = sorted((run_root / "labels").glob("*.tsv"))
    if labels_paths:
        proposals = pd.concat(
            [pd.read_csv(path, sep="\t", keep_default_na=False) for path in labels_paths],
            ignore_index=True,
        ).rename(columns={"label": "cell_type"})
    else:
        proposals = pd.read_csv(
            run_root / "evidence" / "summary" / "annotation_proposals.tsv",
            sep="\t",
            keep_default_na=False,
        )
        proposals["level"] = config.guidance.hierarchy_levels[0]
        proposals["parent_scope"] = "ROOT"
    scores = pd.read_csv(
        run_root / "evidence" / "summary" / "multimodal_label_scores.tsv",
        sep="\t",
        keep_default_na=False,
    )
    cards = run_root / "review" / "evidence_cards"
    cards.mkdir(parents=True, exist_ok=True)
    outputs = []
    occupied_paths: set[Path] = set()
    for proposal in proposals.to_dict("records"):
        cluster = str(proposal["cluster_id"])
        level = str(proposal.get("level", config.guidance.hierarchy_levels[0]))
        parent_scope = str(proposal.get("parent_scope", "ROOT"))
        block = scores[scores["cluster_id"].astype(str).eq(cluster)]
        if "level" in block:
            block = block[block["level"].astype(str).eq(level)]
        if "parent_scope" in block:
            block = block[block["parent_scope"].astype(str).eq(parent_scope)]
        block = block.head(10)
        lines = [
            f"# {level} · {parent_scope} · cluster {cluster}",
            "",
            f"- Parent scope: `{parent_scope}`",
            f"- Proposed cell type: `{proposal['cell_type']}`",
            f"- Cell state: `{proposal['cell_state']}`",
            f"- Confidence: `{float(proposal['confidence']):.3f}`",
            f"- Margin: `{float(proposal['margin']):.3f}`",
            f"- Normalized entropy: `{float(proposal.get('normalized_entropy', 1.0 - float(proposal.get('agreement', 0.0)))):.3f}`",
            f"- Abstained: `{proposal.get('abstained', proposal['cell_type'] == 'Unknown')}`",
            f"- Rationale: {proposal['rationale']}",
            f"- Rejected alternatives: {proposal.get('rejected_alternatives', 'none measured')}",
            "",
            "## Actual multimodal evidence",
            "",
        ]
        if block.empty:
            lines.append("No measurable signed-program evidence.")
        else:
            columns = [
                "label",
                "score",
                "probability",
                "positive_coverage",
                "negative_violation",
                "lane_winners",
                "multimodal_disagreement",
            ]
            columns = [column for column in columns if column in block]
            lines.extend(_markdown_table(block[columns]))
        path = _label_card_path(run_root, level, parent_scope, cluster)
        if path in occupied_paths:
            raise ContractError(
                "evidence-card path collision after identifier sanitization for "
                f"{level}::{parent_scope}::{cluster}"
            )
        occupied_paths.add(path)
        # Cards are derived views of the durable labels table. They are intentionally
        # replaceable so a manual rationale edit is preserved on the next scaffold.
        atomic_replace(path, ("\n".join(lines) + "\n").encode())
        outputs.append(path)
    child_path = run_root / "decisions" / "child_evidence_summary.tsv"
    if child_path.is_file():
        children = pd.read_csv(child_path, sep="\t", keep_default_na=False)
        for child in children.to_dict("records"):
            key = str(child["candidate_key"])
            lines = [
                f"# Candidate child {key}",
                "",
                f"- Parent: `{child['scope']}::{child['parent_id']}`",
                f"- Accepted evidence gate: `{child['accepted_child_evidence_gate']}`",
                f"- Failed gates: `{child['failed_gates']}`",
                f"- Proposed type: `{child['cell_type']}`",
                f"- Stability: `{float(child['stability_score']):.3f}`",
                f"- Separation: `{float(child['separation_score']):.3f}`",
                f"- Technical effect: `{float(child['technical_effect_size']):.3f}`",
                "",
                "## Comparator evidence",
                "",
                "```json",
                json.dumps(
                    json.loads(str(child["lane_evidence_json"])), indent=2, sort_keys=True
                ),
                "```",
            ]
            path = cards / f"candidate_{_safe(key)}.md"
            atomic_replace(path, ("\n".join(lines) + "\n").encode())
            outputs.append(path)
    return outputs


def _safe(value: str) -> str:
    return safe_path_component(value)


def _label_card_path(run_root: Path, level: str, parent_scope: str, cluster: str) -> Path:
    return (
        run_root
        / "review"
        / "evidence_cards"
        / f"{_safe(level)}__{_safe(parent_scope)}__cluster_{_safe(cluster)}.md"
    )


def hierarchy_level_sort_key(value: str) -> tuple[int, str]:
    text = str(value)
    suffix = text.removeprefix("L")
    return (int(suffix), text) if suffix.isdigit() else (10**9, text)


def resolve_durable_mapping_values(
    mapping: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    max_level: str | None = None,
) -> dict[str, np.ndarray]:
    """Resolve cells to the deepest reviewed label along each strict parent path."""

    required_mapping = {"barcode", "leaf_id", "canonical_parent"}
    required_labels = {
        "level",
        "parent_scope",
        "cluster_id",
        "label",
        "cell_state",
    }
    if missing := required_mapping.difference(mapping.columns):
        raise ContractError(f"durable mapping resolution misses columns: {sorted(missing)}")
    if missing := required_labels.difference(labels.columns):
        raise ContractError(f"durable label resolution misses columns: {sorted(missing)}")
    levels = sorted(labels["level"].astype(str).unique(), key=hierarchy_level_sort_key)
    if max_level is not None:
        maximum = hierarchy_level_sort_key(max_level)
        levels = [level for level in levels if hierarchy_level_sort_key(level) <= maximum]
    if not levels:
        raise ContractError("durable mapping resolution has no applicable hierarchy levels")
    by_scope: dict[tuple[str, str], pd.DataFrame] = {
        (str(level), str(parent)): block
        for (level, parent), block in labels.groupby(
            ["level", "parent_scope"], observed=True, sort=False
        )
    }
    cell_types: list[str] = []
    cell_states: list[str] = []
    for row in mapping.to_dict("records"):
        base_identifiers = {
            str(row[column]).strip()
            for column in ("leaf_id", "canonical_parent", "parent_id", "parent_uid")
            if column in row and str(row[column]).strip()
        }
        identifiers = {
            identifier
            for value in base_identifiers
            for identifier in (value, *value.split("::"))
            if identifier
        }
        parent = "ROOT"
        selected_label = ""
        selected_state = ""
        for level in levels:
            block = by_scope.get((level, parent))
            if block is None:
                continue
            matched = block[block["cluster_id"].astype(str).isin(identifiers)]
            if len(matched) > 1:
                raise ContractError(
                    f"cell {row.get('barcode')} matches multiple durable labels at "
                    f"{level}::{parent}"
                )
            if matched.empty:
                continue
            selected = matched.iloc[0]
            selected_label = str(selected["label"]).strip()
            selected_state = str(selected["cell_state"]).strip()
            parent = selected_label
            if selected_label == "Unknown":
                break
        if not selected_label or not selected_state:
            raise ContractError(
                f"durable labels do not resolve barcode {row.get('barcode')} within its parent path"
            )
        cell_types.append(selected_label)
        cell_states.append(selected_state)
    return {
        "cell_type": np.asarray(cell_types, dtype=object),
        "cell_state": np.asarray(cell_states, dtype=object),
        "leaf_id": mapping["leaf_id"].astype(str).to_numpy(),
    }


def _markdown_table(frame: pd.DataFrame) -> list[str]:
    columns = list(map(str, frame.columns))
    rows = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in frame.astype(str).itertuples(index=False, name=None):
        rows.append("| " + " | ".join(value.replace("|", "\\|") for value in row) + " |")
    return rows


def build_critic_packet(config: CellCuratorConfig, run_root: Path) -> Path:
    proposals_path = run_root / "evidence" / "summary" / "annotation_proposals.tsv"
    scores_path = run_root / "evidence" / "summary" / "multimodal_label_scores.tsv"
    decisions_path = run_root / "decisions" / "parent_decision_summary.tsv"
    candidates_path = run_root / "decisions" / "child_evidence_summary.tsv"
    label_paths = sorted((run_root / "labels").glob("*.tsv"))

    def entry(path: Path) -> dict[str, str]:
        return {"path": str(path.relative_to(run_root)), "sha256": sha256(path)}

    packet = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_critic_packet",
        "run_id": config.run.run_id,
        "instructions": [
            "Try to refute each proposed identity from raw evidence, not from its card alone.",
            "Check positive and expected-negative programs, alternative identities, hierarchy, mixed/doublet explanations, and technical confounding.",
            "Record one structured finding per issue; do not directly mutate annotations.",
        ],
        "critic_passes": config.review.critic_passes,
        "artifacts": {
            "proposals": entry(proposals_path),
            "scores": entry(scores_path),
            "parent_decisions": entry(decisions_path),
            "candidate_evidence": entry(candidates_path),
            "labels": [entry(path) for path in label_paths],
        },
    }
    path = run_root / "review" / "critic_packet.json"
    # The packet is a derived index over immutable evidence and must be refreshable
    # when a durable labels table is reviewed or extended.
    replace_json(path, packet)
    findings_path = run_root / "review" / "critic_findings.tsv"
    reconciliation_path = run_root / "review" / "critic_reconciliation.tsv"
    if not findings_path.exists():
        findings = pd.DataFrame(columns=CRITIC_COLUMNS)
        atomic_publish(
            findings_path,
            findings.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
        )
    if not reconciliation_path.exists():
        reconciliation = pd.DataFrame(columns=RECONCILIATION_COLUMNS)
        atomic_publish(
            reconciliation_path,
            reconciliation.to_csv(sep="\t", index=False, lineterminator="\n").encode(),
        )
    return path


def _replace_review_table(path: Path, frame: pd.DataFrame) -> Path:
    """Atomically replace mutable human review state without touching evidence artifacts."""

    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.name}.partial.{os.getpid()}")
    try:
        partial.write_text(frame.to_csv(sep="\t", index=False, lineterminator="\n"))
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)
    return path


def _validate_findings_frame(run_root: Path, frame: pd.DataFrame) -> pd.DataFrame:
    if set(frame.columns) != set(CRITIC_COLUMNS):
        raise ContractError("critic findings do not satisfy the exact structured schema")
    normalized = frame.loc[:, list(CRITIC_COLUMNS)].copy()
    if normalized["finding_id"].astype(str).str.strip().eq("").any():
        raise ContractError("every critic finding requires a non-empty finding_id")
    if normalized["finding_id"].duplicated().any():
        raise ContractError("critic finding IDs must be unique")
    severities = set(normalized["severity"].astype(str).str.upper())
    invalid_severities = severities.difference(CRITIC_SEVERITIES)
    if invalid_severities:
        raise ContractError(f"invalid critic severities: {sorted(invalid_severities)}")
    required_text = (
        "target_kind",
        "target_id",
        "critic_pass",
        "claim",
        "evidence_path",
        "proposed_resolution",
    )
    for column in required_text:
        if normalized[column].astype(str).str.strip().eq("").any():
            raise ContractError(f"critic findings require non-empty {column}")
    root = run_root.resolve()
    for value in normalized["evidence_path"].astype(str):
        supplied = Path(value)
        evidence_path = (supplied if supplied.is_absolute() else run_root / supplied).resolve()
        if not evidence_path.is_relative_to(root):
            raise ContractError(
                f"critic evidence path escapes the run directory: {evidence_path}"
            )
        if not evidence_path.is_file():
            raise ContractError(f"critic evidence path is absent: {evidence_path}")
    normalized["severity"] = normalized["severity"].astype(str).str.upper()
    return normalized


def import_critic_findings(run_root: Path, input_path: Path) -> dict[str, Any]:
    """Validate and atomically import platform-neutral structured critic findings."""

    if not input_path.is_file():
        raise ContractError(f"critic findings input is absent: {input_path}")
    frame = pd.read_csv(input_path, sep="\t", keep_default_na=False)
    missing = set(CRITIC_COLUMNS).difference(frame.columns)
    extra = set(frame.columns).difference(CRITIC_COLUMNS)
    if missing or extra:
        raise ContractError(
            f"critic findings columns must match the schema; missing={sorted(missing)}, "
            f"extra={sorted(extra)}"
        )
    frame = _validate_findings_frame(run_root, frame)
    destination = _replace_review_table(run_root / "review" / "critic_findings.tsv", frame)
    return {
        "schema_version": 3,
        "artifact_kind": "cell_curator_critic_findings_import",
        "status": "complete",
        "n_findings": len(frame),
        "path": str(destination),
        "sha256": sha256(destination),
    }


def import_reconciliation(run_root: Path, input_path: Path) -> dict[str, Any]:
    """Validate and atomically import reviewer dispositions, then enforce completeness."""

    findings_path = run_root / "review" / "critic_findings.tsv"
    if not findings_path.is_file():
        raise ContractError("critic findings must be imported before reconciliation")
    if not input_path.is_file():
        raise ContractError(f"critic reconciliation input is absent: {input_path}")
    frame = pd.read_csv(input_path, sep="\t", keep_default_na=False)
    missing = set(RECONCILIATION_COLUMNS).difference(frame.columns)
    extra = set(frame.columns).difference(RECONCILIATION_COLUMNS)
    if missing or extra:
        raise ContractError(
            f"critic reconciliation columns must match the schema; missing={sorted(missing)}, "
            f"extra={sorted(extra)}"
        )
    frame = frame.loc[:, list(RECONCILIATION_COLUMNS)].copy()
    if frame["finding_id"].duplicated().any():
        raise ContractError("critic reconciliation IDs must be unique")
    frame["status"] = frame["status"].astype(str).str.upper()
    invalid_states = set(frame["status"]).difference(RECONCILIATION_STATES)
    if invalid_states:
        raise ContractError(f"invalid critic reconciliation states: {sorted(invalid_states)}")
    terminal = frame["status"].isin(TERMINAL_RECONCILIATION_STATES)
    if not terminal.all():
        raise ContractError("critic findings remain unresolved")
    for column in ("resolution", "reviewer"):
        if frame.loc[terminal, column].astype(str).str.strip().eq("").any():
            raise ContractError(f"terminal critic dispositions require non-empty {column}")
    findings = _validate_findings_frame(
        run_root,
        pd.read_csv(findings_path, sep="\t", keep_default_na=False),
    )
    expected = set(findings["finding_id"].astype(str))
    observed = set(frame["finding_id"].astype(str))
    if observed != expected:
        raise ContractError(
            "critic reconciliation IDs must exactly match findings; "
            f"missing={sorted(expected - observed)}, extra={sorted(observed - expected)}"
        )
    destination = _replace_review_table(
        run_root / "review" / "critic_reconciliation.tsv", frame
    )
    result = validate_reconciliation(run_root)
    return {
        "schema_version": 3,
        "artifact_kind": "cell_curator_critic_reconciliation_import",
        **result,
        "path": str(destination),
        "sha256": sha256(destination),
    }


def validate_reconciliation(run_root: Path) -> dict[str, Any]:
    findings_path = run_root / "review" / "critic_findings.tsv"
    reconciliation_path = run_root / "review" / "critic_reconciliation.tsv"
    findings = pd.read_csv(findings_path, sep="\t", keep_default_na=False)
    reconciliation = pd.read_csv(reconciliation_path, sep="\t", keep_default_na=False)
    findings = _validate_findings_frame(run_root, findings)
    if set(reconciliation.columns) != set(RECONCILIATION_COLUMNS):
        raise ContractError(
            "critic reconciliation does not satisfy the exact structured schema"
        )
    if reconciliation["finding_id"].duplicated().any():
        raise ContractError("critic finding IDs must be unique")
    indexed = reconciliation.set_index("finding_id")
    expected = set(findings["finding_id"].astype(str))
    observed = set(indexed.index.astype(str))
    if observed != expected:
        raise ContractError(
            "critic reconciliation IDs must exactly match findings; "
            f"missing={sorted(expected - observed)}, extra={sorted(observed - expected)}"
        )
    if findings.empty:
        return {"status": "complete", "n_findings": 0, "all_reconciled": True}
    ordered = indexed.loc[findings["finding_id"].astype(str)]
    terminal = ordered["status"].astype(str).str.upper().isin(TERMINAL_RECONCILIATION_STATES)
    if not terminal.all():
        raise ContractError("critic findings remain unresolved")
    for column in ("resolution", "reviewer"):
        if ordered[column].astype(str).str.strip().eq("").any():
            raise ContractError(f"reconciled critic findings require a non-empty {column}")
    return {"status": "complete", "n_findings": len(findings), "all_reconciled": True}


def validate_critic_packet(run_root: Path) -> dict[str, Any]:
    """Verify that a structured critic packet still names the evidence it hashed."""

    packet_path = run_root / "review" / "critic_packet.json"
    if not packet_path.is_file():
        raise ContractError("strict review requires a structured critic packet")
    try:
        packet = json.loads(packet_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"critic packet is not valid JSON: {exc}") from exc
    if (
        packet.get("schema_version") != 3
        or packet.get("artifact_kind") != "cell_curator_critic_packet"
    ):
        raise ContractError("critic packet does not satisfy its versioned schema")
    artifacts = packet.get("artifacts")
    if not isinstance(artifacts, dict):
        raise ContractError("critic packet artifacts must be a mapping")
    required_artifacts = {
        "proposals",
        "scores",
        "parent_decisions",
        "candidate_evidence",
        "labels",
    }
    missing_artifacts = required_artifacts.difference(artifacts)
    if missing_artifacts or not isinstance(artifacts.get("labels"), list):
        raise ContractError(
            f"critic packet misses required evidence artifacts: {sorted(missing_artifacts)}"
        )
    entries: list[dict[str, Any]] = []
    for value in artifacts.values():
        if isinstance(value, dict):
            entries.append(value)
        elif isinstance(value, list) and all(isinstance(item, dict) for item in value):
            entries.extend(value)
        else:
            raise ContractError("critic packet contains an invalid artifact entry")
    root = run_root.resolve()
    verified = 0
    for entry in entries:
        path_value = str(entry.get("path", "")).strip()
        expected = str(entry.get("sha256", "")).strip()
        if not path_value or not expected:
            raise ContractError("critic packet artifacts require path and SHA-256")
        supplied = Path(path_value)
        path = (supplied if supplied.is_absolute() else run_root / supplied).resolve()
        if not path.is_relative_to(root):
            raise ContractError(f"critic packet artifact escapes the run directory: {path}")
        if not path.is_file() or sha256(path) != expected:
            raise ContractError(f"critic packet artifact is absent or hash-stale: {path}")
        verified += 1
    return {
        "status": "complete",
        "n_artifacts": verified,
        "packet_sha256": sha256(packet_path),
    }


def _report_dependency_hashes(
    run_root: Path,
    labels_paths: list[Path],
    labels: pd.DataFrame,
    *,
    level: str | None,
) -> dict[str, str | None]:
    """Hash every durable input or derived review artifact represented by a report."""

    paths = [
        *labels_paths,
        run_root / "context.json",
        run_root / "assumptions.json",
        run_root / "evidence" / "summary" / "multimodal_label_scores.tsv",
        run_root / "evidence" / "summary" / "cell_state_evidence.tsv",
        run_root / "evidence" / "summary" / "program_cooccurrence.tsv",
        run_root / "evidence" / "summary" / "source_corroboration.tsv",
        run_root / "decisions" / "parent_decision_summary.tsv",
        run_root / "decisions" / "child_evidence_summary.tsv",
        run_root / "mapping" / "parent_child_mapping.tsv",
        run_root / "mapping" / "mapping.validation.json",
        run_root / "adaptive_parent_audit.tsv",
        run_root / "adaptive" / "technical_summary.tsv",
        run_root / "adaptive" / "doublet_summary.tsv",
        run_root / "discovery" / "cluster_profiles.tsv",
        run_root / "review" / (f"colors_{level}.json" if level else "colors.json"),
        run_root / "review" / "critic_packet.json",
        run_root / "review" / "critic_findings.tsv",
        run_root / "review" / "critic_reconciliation.tsv",
    ]
    card_paths = [
        _label_card_path(
            run_root,
            str(row.level),
            str(row.parent_scope),
            str(row.cluster_id),
        )
        for row in labels.itertuples(index=False)
    ]
    if len(set(card_paths)) != len(card_paths):
        raise ContractError("durable labels map to colliding evidence-card paths")
    paths.extend(card_paths)
    return {
        str(path.relative_to(run_root)): sha256(path) if path.is_file() else None
        for path in paths
    }


_PALETTE = (
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
    "#000000",
    "#44AA99",
    "#DDCC77",
)


def _categorical_colors(adata: Any, column: str, values: np.ndarray) -> dict[str, str]:
    labels = list(dict.fromkeys(map(str, values)))
    if column in adata.obs:
        series = adata.obs[column]
        categories = (
            list(map(str, series.cat.categories))
            if isinstance(series.dtype, pd.CategoricalDtype)
            else list(dict.fromkeys(series.astype(str)))
        )
        colors = adata.uns.get(f"{column}_colors")
        if colors is not None and len(colors) == len(categories):
            inherited = dict(zip(categories, map(str, colors), strict=True))
            return {
                label: inherited.get(label, _PALETTE[index % len(_PALETTE)])
                for index, label in enumerate(labels)
            }
    return {label: _PALETTE[index % len(_PALETTE)] for index, label in enumerate(labels)}


def _coordinate_svg(
    coordinates: np.ndarray,
    groups: np.ndarray,
    colors: dict[str, str],
    *,
    label: str,
) -> str:
    if not len(coordinates):
        return ""
    coordinates = np.asarray(coordinates)[:, :2]
    minimum = coordinates.min(axis=0)
    span = np.maximum(coordinates.max(axis=0) - minimum, 1e-12)
    scaled = 20 + 560 * (coordinates - minimum) / span
    points = "".join(
        f'<circle cx="{x:.1f}" cy="{600 - y:.1f}" r="2.1" fill="{html.escape(colors[str(group)])}" opacity="0.65"><title>{html.escape(str(group))}</title></circle>'
        for (x, y), group in zip(scaled, groups, strict=True)
    )
    legend = "".join(
        f'<span style="white-space:nowrap;margin-right:.8rem"><i style="display:inline-block;width:.7rem;height:.7rem;background:{html.escape(color)};margin-right:.25rem"></i>{html.escape(group)}</span>'
        for group, color in colors.items()
    )
    return (
        f'<svg viewBox="0 0 600 600" role="img" aria-label="{html.escape(label)}">'
        f"{points}</svg><div>{legend}</div>"
    )


def _embedding_svg(
    config: CellCuratorConfig,
    *,
    label_by_cluster: dict[str, str] | None = None,
    label_by_barcode: dict[str, str] | None = None,
    preferred_color_column: str = "",
    coordinate_key: str = "",
) -> tuple[str, dict[str, str], list[str]]:
    loaded = load_input(config, backed=False)
    try:
        key = coordinate_key or config.input.keys.embedding
        if key:
            coordinates, barcodes, _ = resolve_obsm(loaded, key)
            coordinates = coordinates[:, :2]
        else:
            adata = next(iter(loaded.assays.values()))
            matrix = dense(adata.X)
            centered = matrix - matrix.mean(axis=0)
            u, singular, _ = np.linalg.svd(centered, full_matrices=False)
            coordinates = u[:, :2] * singular[:2]
            barcodes = adata.obs_names.astype(str)
        root_obs = loaded.root.obs.copy()
        root_obs.index = root_obs.index.astype(str)
        canonical = (
            root_obs.loc[barcodes, config.input.keys.canonical_parent].astype(str).to_numpy()
        )
        if label_by_barcode is not None:
            groups = np.asarray(
                [label_by_barcode.get(str(barcode), "Unknown") for barcode in barcodes]
            )
        elif label_by_cluster:
            groups = np.asarray([label_by_cluster.get(item, "Unknown") for item in canonical])
        else:
            groups = canonical
        color_column = (
            config.input.keys.canonical_parent
            if label_by_cluster is None and label_by_barcode is None
            else preferred_color_column
        )
        color_source = loaded.root
        if color_column not in color_source.obs:
            color_source = next(
                (
                    assay
                    for assay in loaded.assays.values()
                    if color_column in assay.obs
                ),
                loaded.root,
            )
        colors = _categorical_colors(color_source, color_column, groups)
        observed = list(dict.fromkeys(map(str, groups)))
        if color_column in color_source.obs and isinstance(
            color_source.obs[color_column].dtype, pd.CategoricalDtype
        ):
            existing = list(map(str, color_source.obs[color_column].cat.categories))
            categories = [item for item in existing if item in observed]
            categories.extend(item for item in observed if item not in categories)
        else:
            categories = observed
        colors = {item: colors[item] for item in categories}
    finally:
        close_input(loaded)
    return (
        _coordinate_svg(
            coordinates, groups, colors, label="Frozen embedding with reviewed labels"
        ),
        colors,
        categories,
    )


def _html_table(frame: pd.DataFrame, *, limit: int = 500) -> str:
    if frame.empty:
        return "<p>No rows.</p>"
    return frame.head(limit).to_html(index=False, escape=True, border=0, classes="data")


def _uncertainty_svg(labels: pd.DataFrame) -> str:
    confidence = pd.to_numeric(
        labels.get("confidence", pd.Series(dtype=float)), errors="coerce"
    )
    confidence = confidence.dropna().clip(0, 1)
    if confidence.empty:
        return "<p>Confidence unavailable.</p>"
    counts, edges = np.histogram(confidence, bins=np.linspace(0, 1, 11))
    maximum = max(int(counts.max()), 1)
    bars = []
    for index, count in enumerate(counts):
        x = 35 + index * 52
        height = 220 * int(count) / maximum
        bars.append(
            f'<rect x="{x}" y="{250 - height:.1f}" width="42" height="{height:.1f}" fill="#4477AA"><title>{edges[index]:.1f}–{edges[index + 1]:.1f}: {count}</title></rect>'
        )
    return (
        '<svg viewBox="0 0 600 280" role="img" aria-label="Advisory confidence distribution">'
        + "".join(bars)
        + '<line x1="30" x2="570" y1="250" y2="250" stroke="#333"/></svg>'
    )


def _score_dotplot_svg(scores: pd.DataFrame, *, limit_labels: int = 18) -> str:
    required = {"cluster_id", "label", "score"}
    if scores.empty or required.difference(scores.columns):
        return "<p>Marker dotplot unavailable.</p>"
    ranked = scores.sort_values("score", ascending=False)
    labels = list(dict.fromkeys(ranked["label"].astype(str)))[:limit_labels]
    clusters = sorted(ranked["cluster_id"].astype(str).unique())
    block = ranked[
        ranked["label"].astype(str).isin(labels)
        & ranked["cluster_id"].astype(str).isin(clusters)
    ]
    lookup = {
        (str(row.cluster_id), str(row.label)): row for row in block.itertuples(index=False)
    }
    width = max(500, 100 + len(labels) * 35)
    height = max(260, 80 + len(clusters) * 38)
    score_values = pd.to_numeric(block["score"], errors="coerce")
    low, high = float(score_values.min()), float(score_values.max())
    span = max(high - low, 1e-12)
    marks = []
    for row_index, cluster in enumerate(clusters):
        for column_index, label in enumerate(labels):
            row = lookup.get((cluster, label))
            if row is None:
                continue
            score = float(row.score)
            fraction = float(getattr(row, "positive_coverage", 0.0))
            red = int(230 * (score - low) / span)
            radius = 3 + 8 * max(0.0, min(1.0, fraction))
            marks.append(
                f'<circle cx="{90 + 35 * column_index}" cy="{45 + 38 * row_index}" r="{radius:.1f}" fill="rgb({red},80,{230 - red})"><title>{html.escape(cluster)} · {html.escape(label)} · score={score:.3f} · coverage={fraction:.3f}</title></circle>'
            )
    xlabels = "".join(
        f'<text x="{90 + 35 * index}" y="{height - 8}" transform="rotate(-55 {90 + 35 * index},{height - 8})" font-size="10">{html.escape(label)}</text>'
        for index, label in enumerate(labels)
    )
    ylabels = "".join(
        f'<text x="5" y="{49 + 38 * index}" font-size="11">{html.escape(cluster)}</text>'
        for index, cluster in enumerate(clusters)
    )
    return (
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Marker score dotplot">'
        + "".join(marks)
        + xlabels
        + ylabels
        + "</svg>"
    )


def _hierarchy_html(labels: pd.DataFrame) -> str:
    if labels.empty:
        return "<p>No hierarchy labels.</p>"
    parts = ["<ul class=tree>"]
    for (level, parent), block in labels.groupby(
        ["level", "parent_scope"], observed=True, sort=True
    ):
        children = ", ".join(sorted(block["label"].astype(str).unique()))
        parts.append(
            f"<li><strong>{html.escape(str(level))}</strong> · {html.escape(str(parent))} → {html.escape(children)}</li>"
        )
    parts.append("</ul>")
    return "".join(parts)


def _spatial_svg(
    config: CellCuratorConfig,
    *,
    label_by_cluster: dict[str, str] | None,
    label_by_barcode: dict[str, str] | None,
    colors: dict[str, str],
) -> str:
    if not config.input.keys.spatial:
        return "<p>Spatial coordinates not configured.</p>"
    loaded = load_input(config, backed=False)
    try:
        key = config.input.keys.spatial
        coordinates, barcodes, _ = resolve_obsm(loaded, key)
        root_obs = loaded.root.obs.copy()
        root_obs.index = root_obs.index.astype(str)
        canonical = root_obs.loc[barcodes, config.input.keys.canonical_parent].astype(str)
        if label_by_barcode is not None:
            groups = np.asarray(
                [label_by_barcode.get(str(barcode), "Unknown") for barcode in barcodes]
            )
        else:
            cluster_lookup = label_by_cluster or {}
            groups = np.asarray([cluster_lookup.get(item, "Unknown") for item in canonical])
        local_colors = {
            label: colors.get(label, _PALETTE[index % len(_PALETTE)])
            for index, label in enumerate(dict.fromkeys(groups))
        }
        return _coordinate_svg(
            coordinates, groups, local_colors, label="Spatial label context"
        )
    finally:
        close_input(loaded)


def _composition_panels(
    config: CellCuratorConfig,
    *,
    label_by_cluster: dict[str, str] | None,
    label_by_barcode: dict[str, str] | None,
    covariates: list[str],
) -> str:
    if not covariates:
        return "<p>No composition covariates selected.</p>"
    loaded = load_input(config, backed=False)
    panels = []
    try:
        canonical = loaded.root.obs[config.input.keys.canonical_parent].astype(str)
        if label_by_barcode is not None:
            values = [
                label_by_barcode.get(str(barcode), "Unknown")
                for barcode in loaded.root.obs_names
            ]
        else:
            cluster_lookup = label_by_cluster or {}
            values = [cluster_lookup.get(item, "Unknown") for item in canonical]
        labels = pd.Series(values, index=loaded.root.obs_names, name="cell_type")
        for column in covariates:
            if column not in loaded.root.obs:
                panels.append(f"<h3>{html.escape(column)}</h3><p>Column absent.</p>")
                continue
            table = pd.crosstab(
                labels,
                loaded.root.obs[column].astype(str),
                normalize="index",
            ).reset_index()
            panels.append(f"<h3>{html.escape(column)}</h3>{_html_table(table)}")
    finally:
        close_input(loaded)
    return "".join(panels)


def check_report_completeness(
    config: CellCuratorConfig,
    run_root: Path,
    *,
    level: str | None = None,
    require_report: bool = True,
) -> dict[str, Any]:
    if level is not None:
        level = require_safe_path_component(level, field="report level")
    all_labels_paths = sorted((run_root / "labels").glob("*.tsv"))
    labels_paths = [run_root / "labels" / f"{level}.tsv"] if level else all_labels_paths
    missing_files = [str(path) for path in labels_paths if not path.is_file()]
    if not labels_paths or missing_files:
        raise ContractError(f"strict report requires durable labels tables: {missing_files}")
    labels = pd.concat(
        [pd.read_csv(path, sep="\t", keep_default_na=False) for path in labels_paths],
        ignore_index=True,
    )
    dependency_labels_paths = (
        [
            path
            for path in all_labels_paths
            if hierarchy_level_sort_key(path.stem) <= hierarchy_level_sort_key(level)
        ]
        if level
        else all_labels_paths
    )
    dependency_labels = pd.concat(
        [
            pd.read_csv(path, sep="\t", keep_default_na=False)
            for path in dependency_labels_paths
        ],
        ignore_index=True,
    )
    required = {
        "level",
        "parent_scope",
        "cluster_id",
        "label",
        "confidence",
        "rationale",
        "rejected_alternatives",
    }
    missing_columns = required.difference(dependency_labels.columns)
    if missing_columns:
        raise ContractError(f"strict report labels miss columns: {sorted(missing_columns)}")
    durable_keys = ["level", "parent_scope", "cluster_id"]
    if dependency_labels.duplicated(durable_keys).any():
        raise ContractError("strict report labels duplicate a level/parent/cluster key")
    if level is not None and set(labels["level"].astype(str)) != {level}:
        raise ContractError(f"strict level report contains labels outside {level!r}")
    text_columns = ["parent_scope", "cluster_id", "label", "rationale", "rejected_alternatives"]
    blank = {
        column: int(
            dependency_labels[column]
            .astype(str)
            .str.strip()
            .isin({"", "TODO", "PENDING"})
            .sum()
        )
        for column in text_columns
    }
    blank = {column: count for column, count in blank.items() if count}
    confidence = pd.to_numeric(dependency_labels["confidence"], errors="coerce")
    if blank or confidence.isna().any() or ((confidence < 0) | (confidence > 1)).any():
        raise ContractError(
            f"strict report has incomplete labels; blank={blank}, invalid_confidence={int(confidence.isna().sum())}"
        )
    state = load_state(run_root)
    assumptions = load_assumptions(run_root)
    pending = [
        item["assumption_id"]
        for item in assumptions
        if bool(item.get("material", True)) and item.get("status") != "APPROVED"
    ]
    if pending or not state.assumptions_approved:
        raise ContractError(f"strict report requires approved assumptions; pending={pending}")
    assumptions_path = run_root / "assumptions.json"
    if (
        not state.approved_assumptions_sha256
        or not assumptions_path.is_file()
        or sha256(assumptions_path) != state.approved_assumptions_sha256
    ):
        raise ContractError(
            "strict report assumptions ledger changed after its recorded approval"
        )
    report_path = run_root / "review" / (f"report_{level}.html" if level else "report.html")
    if require_report and not report_path.is_file():
        raise ContractError(
            f"strict check-only validation requires existing report: {report_path}"
        )
    if require_report:
        manifest_path = (
            run_root
            / "review"
            / (f"report_{level}.manifest.json" if level else "report.manifest.json")
        )
        if not manifest_path.is_file():
            raise ContractError(f"strict report manifest is absent: {manifest_path}")
        manifest = json.loads(manifest_path.read_text())
        expected_panels = {
            "embedding",
            "marker_dotplot",
            "hierarchy",
            "uncertainty",
            "rejected_alternatives",
            "qc",
            "spatial",
            "composition",
            "critics",
        }
        missing_panels = expected_panels.difference(manifest.get("panels", []))
        expected_dependencies = _report_dependency_hashes(
            run_root, dependency_labels_paths, dependency_labels, level=level
        )
        recorded_dependencies = manifest.get("dependency_sha256")
        if (
            manifest.get("status") != "complete"
            or manifest.get("report_sha256") != sha256(report_path)
            or manifest.get("source_sha256") != sha256(Path(config.input.path))
            or missing_panels
            or recorded_dependencies != expected_dependencies
        ):
            raise ContractError(
                "strict report manifest is incomplete or stale; "
                f"missing_panels={sorted(missing_panels)}, "
                f"dependencies_match={recorded_dependencies == expected_dependencies}"
            )
        document = report_path.read_text()
        required_headings = (
            "Reviewed embedding",
            "Marker dotplot",
            "Hierarchy",
            "Uncertainty and confidence",
            "Rejected alternatives",
            "Technical audit",
            "Refined child profiles",
            "Spatial context",
            "Composition panels",
            "Critic findings",
            "Critic reconciliation",
            "Provenance",
        )
        absent_headings = [heading for heading in required_headings if heading not in document]
        external_resource = re.search(
            r"<(?:script|link|img)\b[^>]*(?:src|href)\s*=\s*['\"](?!data:|#)[^'\"]+",
            document,
            flags=re.IGNORECASE,
        )
        if absent_headings or external_resource:
            raise ContractError(
                "strict self-contained report is incomplete; "
                f"missing_sections={absent_headings}"
            )
        colors_path = run_root / "review" / (f"colors_{level}.json" if level else "colors.json")
        critic_path = run_root / "review" / "critic_packet.json"
        expected_cards = [
            _label_card_path(
                run_root,
                str(row.level),
                str(row.parent_scope),
                str(row.cluster_id),
            )
            for row in dependency_labels.itertuples(index=False)
        ]
        if (
            not colors_path.is_file()
            or not critic_path.is_file()
            or any(not path.is_file() for path in expected_cards)
        ):
            raise ContractError(
                "strict report requires colors, critic packet, and an evidence card per label"
            )
        validate_critic_packet(run_root)
        validate_reconciliation(run_root)
        colors = json.loads(colors_path.read_text())
        for contract in colors.get("columns", {}).values():
            categories = list(map(str, contract.get("categories", [])))
            values = list(map(str, contract.get("colors", [])))
            if len(categories) != len(values) or len(set(categories)) != len(categories):
                raise ContractError("strict report color categories and values are misaligned")
    return {
        "schema_version": 3,
        "artifact_kind": "cell_curator_report_completeness",
        "status": "complete",
        "level": level or "full-run",
        "n_labels": len(labels),
        "assumptions_approved": True,
        "report": str(report_path),
    }


def build_html_report(
    config: CellCuratorConfig,
    run_root: Path,
    *,
    strict: bool = False,
    level: str | None = None,
    check_only: bool = False,
) -> Path:
    if level is not None:
        level = require_safe_path_component(level, field="report level")
    if check_only:
        if not strict:
            raise ContractError("report --check-only requires --strict")
        result = check_report_completeness(config, run_root, level=level, require_report=True)
        return Path(result["report"])
    all_labels_paths = sorted((run_root / "labels").glob("*.tsv"))
    labels_paths = [run_root / "labels" / f"{level}.tsv"] if level else all_labels_paths
    if labels_paths and all(path.is_file() for path in labels_paths):
        labels = pd.concat(
            [pd.read_csv(path, sep="\t", keep_default_na=False) for path in labels_paths],
            ignore_index=True,
        )
    else:
        proposals = pd.read_csv(
            run_root / "evidence" / "summary" / "annotation_proposals.tsv",
            sep="\t",
            keep_default_na=False,
        )
        proposals["level"] = config.guidance.hierarchy_levels[0]
        proposals["parent_scope"] = "ROOT"
        labels = proposals.rename(columns={"cell_type": "label"})
    if all_labels_paths:
        all_labels = pd.concat(
            [pd.read_csv(path, sep="\t", keep_default_na=False) for path in all_labels_paths],
            ignore_index=True,
        )
    else:
        all_labels = labels.copy()
    dependency_labels_paths = (
        [
            path
            for path in all_labels_paths
            if hierarchy_level_sort_key(path.stem) <= hierarchy_level_sort_key(level)
        ]
        if level
        else all_labels_paths
    )
    dependency_labels = (
        all_labels[
            all_labels["level"]
            .astype(str)
            .map(
                lambda value: hierarchy_level_sort_key(value) <= hierarchy_level_sort_key(level)
            )
        ]
        if level
        else all_labels
    )
    scores = pd.read_csv(
        run_root / "evidence" / "summary" / "multimodal_label_scores.tsv",
        sep="\t",
        keep_default_na=False,
    )
    report_scores = scores
    if level is not None and "level" in scores:
        report_scores = scores[scores["level"].astype(str).eq(level)]
    state_path = run_root / "evidence" / "summary" / "cell_state_evidence.tsv"
    states = (
        pd.read_csv(state_path, sep="\t", keep_default_na=False)
        if state_path.is_file()
        else pd.DataFrame()
    )
    cooccurrence = pd.read_csv(
        run_root / "evidence" / "summary" / "program_cooccurrence.tsv",
        sep="\t",
        keep_default_na=False,
    )
    corroboration = pd.read_csv(
        run_root / "evidence" / "summary" / "source_corroboration.tsv",
        sep="\t",
        keep_default_na=False,
    )
    parents = pd.read_csv(
        run_root / "decisions" / "parent_decision_summary.tsv", sep="\t", keep_default_na=False
    )
    candidates = pd.read_csv(
        run_root / "decisions" / "child_evidence_summary.tsv", sep="\t", keep_default_na=False
    )
    mapping_path = run_root / "mapping" / "parent_child_mapping.tsv"
    mapping = (
        pd.read_csv(mapping_path, sep="\t", keep_default_na=False)
        if mapping_path.is_file()
        else pd.DataFrame()
    )
    audit_path = run_root / "adaptive_parent_audit.tsv"
    audit = (
        pd.read_csv(audit_path, sep="\t", keep_default_na=False)
        if audit_path.is_file()
        else pd.DataFrame()
    )
    technical = pd.read_csv(
        run_root / "adaptive" / "technical_summary.tsv",
        sep="\t",
        keep_default_na=False,
    )
    doublets = pd.read_csv(
        run_root / "adaptive" / "doublet_summary.tsv",
        sep="\t",
        keep_default_na=False,
    )
    critic_findings_path = run_root / "review" / "critic_findings.tsv"
    critic_reconciliation_path = run_root / "review" / "critic_reconciliation.tsv"
    critic_findings = (
        pd.read_csv(critic_findings_path, sep="\t", keep_default_na=False)
        if critic_findings_path.is_file()
        else pd.DataFrame(columns=CRITIC_COLUMNS)
    )
    critic_reconciliation = (
        pd.read_csv(critic_reconciliation_path, sep="\t", keep_default_na=False)
        if critic_reconciliation_path.is_file()
        else pd.DataFrame(columns=RECONCILIATION_COLUMNS)
    )
    cluster_profiles_path = run_root / "discovery" / "cluster_profiles.tsv"
    cluster_profiles = (
        pd.read_csv(cluster_profiles_path, sep="\t", keep_default_na=False)
        if cluster_profiles_path.is_file()
        else pd.DataFrame()
    )
    if strict:
        check_report_completeness(config, run_root, level=level, require_report=False)
        validation = json.loads((run_root / "mapping" / "mapping.validation.json").read_text())
        mapping_path = run_root / "mapping" / "parent_child_mapping.tsv"
        if validation.get("status") != "complete" or validation.get("mapping_sha256") != sha256(
            mapping_path
        ):
            raise ContractError("strict report requires a validated mapping")
        validate_critic_packet(run_root)
        validate_reconciliation(run_root)
    style = """
body{font-family:system-ui,sans-serif;max-width:1400px;margin:2rem auto;padding:0 1rem;color:#17212b}
h1,h2{color:#123a52}.grid{display:grid;grid-template-columns:minmax(300px,1fr) minmax(500px,2fr);gap:1.5rem}
svg{width:100%;height:auto;background:#f7fafb;border:1px solid #d7e0e5}.data{border-collapse:collapse;width:100%;font-size:.82rem}
.data th,.data td{border:1px solid #d7e0e5;padding:.35rem;text-align:left}.data th{background:#eaf1f4;position:sticky;top:0}
.scroll{overflow:auto;max-height:34rem}.flag{color:#9b2c2c;font-weight:700}.ok{color:#276749;font-weight:700}
code{background:#eef2f4;padding:.1rem .25rem}section{margin:2rem 0}.tree{line-height:1.7}
"""
    unknown = int(labels["label"].astype(str).eq("Unknown").sum())
    disagreement = int(
        scores.get("multimodal_disagreement", pd.Series(dtype=bool))
        .astype(str)
        .str.lower()
        .eq("true")
        .sum()
    )
    context_path = run_root / "context.json"
    covariates = list(config.guidance.context.composition_covariates)
    prefix = config.guidance.context.writeback_prefix
    visualization_key = config.guidance.context.visualization_key or config.input.keys.embedding
    if context_path.is_file():
        context = json.loads(context_path.read_text())
        biological = context.get("biological_context", {})
        covariates = list(biological.get("composition_covariates", covariates))
        prefix = str(biological.get("writeback_prefix", prefix))
        visualization_key = str(
            context.get("contracts", {}).get("visualization_key", visualization_key)
        )
    display_level = level or config.guidance.hierarchy_levels[0]
    first = labels[labels["level"].astype(str).eq(display_level)]
    if first.empty:
        raise ContractError(f"report contains no durable labels for level {display_level!r}")
    label_by_cluster = dict(
        zip(first["cluster_id"].astype(str), first["label"].astype(str), strict=True)
    )
    label_by_barcode: dict[str, str] | None = None
    if not mapping.empty:
        resolved = resolve_durable_mapping_values(
            mapping,
            all_labels,
            max_level=level,
        )
        label_by_barcode = dict(
            zip(
                mapping["barcode"].astype(str),
                map(str, resolved["cell_type"]),
                strict=True,
            )
        )
    embedding, colors, categories = _embedding_svg(
        config,
        label_by_cluster=label_by_cluster,
        label_by_barcode=label_by_barcode,
        preferred_color_column=f"{prefix}_type",
        coordinate_key=visualization_key,
    )
    color_manifest = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_report_colors",
        "columns": {
            f"{prefix}_type": {
                "categories": categories,
                "colors": [colors[item] for item in categories],
            }
        },
    }
    colors_path = run_root / "review" / (f"colors_{level}.json" if level else "colors.json")
    replace_json(colors_path, color_manifest)
    spatial = _spatial_svg(
        config,
        label_by_cluster=label_by_cluster,
        label_by_barcode=label_by_barcode,
        colors=colors,
    )
    composition = _composition_panels(
        config,
        label_by_cluster=label_by_cluster,
        label_by_barcode=label_by_barcode,
        covariates=covariates,
    )
    report_scope = level or "full run"
    mapping_summary = (
        f"{len(mapping)} cells, {mapping['leaf_id'].nunique()} leaves."
        if not mapping.empty and "leaf_id" in mapping
        else "Mapping not yet constructed."
    )
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>cell-curator {html.escape(config.run.run_id)}</title><style>{style}</style></head><body>
<h1>cell-curator review: {html.escape(config.run.run_id)}</h1>
<p><strong>Scope:</strong> {html.escape(report_scope)} · <strong>Mode:</strong> {html.escape(config.run.mode.value)} · <strong>Unknown:</strong> {unknown}/{len(labels)} · <strong>Multimodal disagreements:</strong> {disagreement}</p>
<div class="grid"><section><h2>Reviewed embedding</h2>{embedding}</section><section><h2>Labels source of truth</h2><div class="scroll">{_html_table(labels)}</div></section></div>
<section><h2>Hierarchy</h2>{_hierarchy_html(labels)}</section>
<section><h2>Marker dotplot</h2>{_score_dotplot_svg(report_scores)}</section>
<section><h2>Uncertainty and confidence</h2>{_uncertainty_svg(labels)}</section>
<section><h2>Rejected alternatives</h2><div class="scroll">{_html_table(labels[[column for column in ["level", "parent_scope", "cluster_id", "label", "rejected_alternatives", "rationale"] if column in labels]])}</div></section>
<section><h2>Multimodal label evidence</h2><div class="scroll">{_html_table(scores)}</div></section>
<section><h2>Independent cell-state evidence</h2><div class="scroll">{_html_table(states)}</div></section>
<section><h2>Program co-occurrence</h2><p>Cluster-aggregate co-support; this does not imply single-cell co-detection.</p><div class="scroll">{_html_table(cooccurrence)}</div></section>
<section><h2>Independent-source corroboration</h2><div class="scroll">{_html_table(corroboration)}</div></section>
<section><h2>Adaptive candidates</h2><div class="scroll">{_html_table(candidates)}</div></section>
<section><h2>Parent decisions</h2><div class="scroll">{_html_table(parents)}</div></section>
<section><h2>Technical audit</h2><div class="scroll">{_html_table(audit)}</div></section>
<section><h2>Refined child profiles</h2><div class="scroll">{_html_table(cluster_profiles)}</div></section>
<section><h2>Technical composition summary</h2><div class="scroll">{_html_table(technical)}</div></section>
<section><h2>Doublet and incompatible-program summary</h2><div class="scroll">{_html_table(doublets)}</div></section>
<section><h2>Spatial context</h2>{spatial}</section>
<section><h2>Composition panels</h2>{composition}</section>
<section><h2>Critic findings</h2><div class="scroll">{_html_table(critic_findings)}</div></section>
<section><h2>Critic reconciliation</h2><div class="scroll">{_html_table(critic_reconciliation)}</div></section>
<section><h2>Immutable mapping</h2><p>{mapping_summary}</p><div class="scroll">{_html_table(mapping, limit=200)}</div></section>
<section><h2>Provenance</h2><ul><li>Input: <code>{html.escape(str(config.input.path))}</code></li><li>Input SHA-256: <code>{sha256(Path(config.input.path))}</code></li><li>Mapping SHA-256: <code>{sha256(mapping_path) if mapping_path.is_file() else "unavailable"}</code></li></ul></section>
</body></html>"""
    path = run_root / "review" / (f"report_{level}.html" if level else "report.html")
    atomic_replace(path, document.encode())
    manifest_path = (
        run_root
        / "review"
        / (f"report_{level}.manifest.json" if level else "report.manifest.json")
    )
    dependency_hashes = _report_dependency_hashes(
        run_root,
        dependency_labels_paths or labels_paths,
        dependency_labels,
        level=level,
    )
    replace_json(
        manifest_path,
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_review_report",
            "status": "complete",
            "run_id": config.run.run_id,
            "strict": strict,
            "level": level or "full-run",
            "report_sha256": sha256(path),
            "source_sha256": sha256(Path(config.input.path)),
            "dependency_sha256": dependency_hashes,
            "n_unknown": unknown,
            "n_multimodal_disagreements": disagreement,
            "panels": [
                "embedding",
                "marker_dotplot",
                "hierarchy",
                "uncertainty",
                "rejected_alternatives",
                "qc",
                "spatial",
                "composition",
                "critics",
            ],
            "colors": str(colors_path),
        },
    )
    if strict:
        check_report_completeness(config, run_root, level=level, require_report=True)
    return path

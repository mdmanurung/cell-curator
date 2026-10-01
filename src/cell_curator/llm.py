"""A gated LLM evidence lane.

An LLM never labels a cluster on its own say-so. It may only break a near-tie
that the measured evidence already left open, and only for a candidate that
passes its own marker gates. The flow:

1. A completed run abstains (or backs off) on some clusters. :func:`build_packet`
   describes those clusters: their near-tied candidates, top markers, cluster
   statistics, and, for clusters of at least ``evidence.llm.min_cells`` cells,
   pseudobulk means of the top markers. Nothing per cell, no barcodes, no donors.
2. A model answers, through ``cell-curator llm propose`` (an OpenAI-compatible
   endpoint reached through net.py) or a host agent. Either way the answer is a
   :class:`ProposalSet` carrying the packet's sha256.
3. A new run names the proposals file in ``evidence.llm.proposals``. Plan binds
   its hash. Aggregation rebuilds the packet from the new run's own evidence and
   refuses proposals made against a different one, then applies
   :func:`apply_proposals`.

Precedence: a validated proposal that every model agrees on, then ontology
backoff, then Unknown.
"""

from __future__ import annotations

import hashlib
import json
import os
import urllib.request
from pathlib import Path
from typing import Any, Literal

import pandas as pd
from pydantic import Field

from .config import CellCuratorConfig
from .contracts import ContractError, as_bool, measured_float, sha256
from .models import AnnotationDecision, StrictModel
from .net import urlopen

TEMPLATE_VERSION = "cell-curator-llm-v1"
TOP_MARKERS = 20
PREFIX = "LLM tie-break"


class Proposal(StrictModel):
    cluster_id: str
    label: str
    cl_id: str = ""
    cited_genes: tuple[str, ...] = ()
    rationale: str = ""
    model: str


class ProposalSet(StrictModel):
    schema_version: Literal[3] = 3
    payload_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    template: str = TEMPLATE_VERSION
    proposals: list[Proposal]


def packet_sha256(packet: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical(packet)).hexdigest()


def _canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def _round(value: object) -> float | None:
    number = measured_float(value)
    return None if number is None else round(number, 4)


def _blocked(decision: AnnotationDecision) -> bool:
    return decision.abstained or decision.rationale.startswith("Backed off")


def build_packet(
    config: CellCuratorConfig,
    decisions: list[AnnotationDecision],
    aggregate: pd.DataFrame,
    bottom_up: pd.DataFrame,
    cluster_sizes: dict[str, int],
) -> dict[str, Any]:
    """What a model may see: blocked clusters, described by cluster aggregates only."""

    margin = config.evidence.min_margin
    clusters = []
    for decision in sorted(decisions, key=lambda item: item.cluster_id):
        if not _blocked(decision):
            continue
        rows = aggregate[aggregate["cluster_id"].astype(str).eq(decision.cluster_id)]
        rows = rows.assign(_p=pd.to_numeric(rows["probability"], errors="coerce"))
        if rows.empty or rows["_p"].isna().all():
            continue
        near = rows[rows["_p"] >= rows["_p"].max() - margin].sort_values(
            ["_p", "label"], ascending=[False, True]
        )
        if len(near) < 2:
            continue
        markers = bottom_up[
            bottom_up["cluster_id"].astype(str).eq(decision.cluster_id)
            & (pd.to_numeric(bottom_up["effect"], errors="coerce") > 0)
        ]
        markers = markers.assign(_r=pd.to_numeric(markers["rank"], errors="coerce"))
        markers = markers.sort_values(["_r", "feature"]).drop_duplicates("feature")
        markers = markers.head(TOP_MARKERS)
        n_cells = int(cluster_sizes.get(decision.cluster_id, 0))
        entry: dict[str, Any] = {
            "cluster_id": decision.cluster_id,
            "n_cells": n_cells,
            "candidates": [
                {
                    "label": str(row["label"]),
                    "cl_id": str(row.get("cl_id", "")),
                    "vote_share": _round(row["_p"]),
                    "positive_coverage": _round(row["positive_coverage"]),
                    "negative_violation": _round(row["negative_violation"]),
                }
                for row in near.to_dict("records")
            ],
            "top_markers": [str(item) for item in markers["feature"]],
        }
        if n_cells >= config.evidence.llm.min_cells:
            entry["pseudobulk"] = {
                str(row["feature"]): {
                    "cluster_mean": _round(row["target_mean"]),
                    "rest_mean": _round(row["reference_mean"]),
                    "cluster_fraction_detected": _round(row["target_fraction_detected"]),
                }
                for row in markers.to_dict("records")
            }
        clusters.append(entry)
    context = config.guidance.context
    return {
        "schema_version": 3,
        "template": TEMPLATE_VERSION,
        "context": {
            "organism": context.organism,
            "tissue": context.tissue,
            "developmental_stage": context.developmental_stage,
            "condition": context.condition,
            "technology": context.technology,
        },
        "clusters": clusters,
    }


def apply_proposals(
    config: CellCuratorConfig,
    packet: dict[str, Any],
    proposal_set: ProposalSet,
    decisions: list[AnnotationDecision],
    aggregate: pd.DataFrame,
) -> tuple[list[AnnotationDecision], list[dict[str, Any]]]:
    """Let agreeing, gate-passing proposals break near-ties. Returns decisions and outcomes."""

    if proposal_set.payload_sha256 != packet_sha256(packet):
        raise ContractError(
            "LLM proposals were made against a different packet than this run measured "
            f"({proposal_set.payload_sha256[:12]} != {packet_sha256(packet)[:12]})"
        )
    by_cluster: dict[str, list[Proposal]] = {}
    for proposal in proposal_set.proposals:
        by_cluster.setdefault(proposal.cluster_id, []).append(proposal)
    entries = {entry["cluster_id"]: entry for entry in packet["clusters"]}
    outcomes: list[dict[str, Any]] = []
    updated = {decision.cluster_id: decision for decision in decisions}
    for cluster_id, proposals in sorted(by_cluster.items()):
        models = sorted({item.model for item in proposals})
        entry = entries.get(cluster_id)
        label, reason = _judge(config, entry, proposals, aggregate)
        outcomes.append(
            {
                "cluster_id": cluster_id,
                "label": label,
                "models": models,
                "status": "refused" if reason else "accepted",
                "reason": reason,
            }
        )
        if reason or entry is None:
            if cluster_id in updated:
                decision = updated[cluster_id]
                updated[cluster_id] = decision.model_copy(
                    update={
                        "rationale": f"{decision.rationale}; LLM proposal "
                        f"{label or 'set'} from {', '.join(models)} refused: {reason}"
                    }
                )
            continue
        candidate = next(item for item in entry["candidates"] if item["label"] == label)
        cited = {gene for item in proposals for gene in item.cited_genes}
        decision = updated[cluster_id]
        updated[cluster_id] = decision.model_copy(
            update={
                "cell_type": label,
                "cl_id": candidate["cl_id"] or None,
                "confidence": float(candidate["vote_share"] or 0.0),
                "abstained": False,
                "rationale": (
                    f"{PREFIX} by {', '.join(models)} among near-tied candidates "
                    f"{', '.join(item['label'] for item in entry['candidates'])}; cites "
                    f"{', '.join(sorted(cited))}. Measured vote share "
                    f"{candidate['vote_share']}; the markers admitted it, the model chose it. "
                    f"Before: {decision.rationale}"
                ),
            }
        )
    return [updated[decision.cluster_id] for decision in decisions], outcomes


def _judge(
    config: CellCuratorConfig,
    entry: dict[str, Any] | None,
    proposals: list[Proposal],
    aggregate: pd.DataFrame,
) -> tuple[str, str]:
    """The proposed label and why it is refused, or an empty reason if it stands."""

    labels = {item.label for item in proposals}
    label = next(iter(labels)) if len(labels) == 1 else ""
    if entry is None:
        return label, "cluster is not blocked by the gates, so there is no tie to break"
    if len(labels) != 1:
        return label, f"models disagree: {sorted(labels)}"
    candidate = next((item for item in entry["candidates"] if item["label"] == label), None)
    if candidate is None:
        return label, "label is not one of the near-tied candidates"
    if any(item.cl_id and item.cl_id != candidate["cl_id"] for item in proposals):
        return label, "CL id does not match the candidate's"
    row = aggregate[
        aggregate["cluster_id"].astype(str).eq(entry["cluster_id"])
        & aggregate["label"].astype(str).eq(label)
    ].iloc[0]
    thresholds = config.evidence
    coverage = measured_float(row["positive_coverage"])
    program = measured_float(row["program_coverage"])
    violation = measured_float(row["negative_violation"])
    if not (
        coverage is not None
        and program is not None
        and violation is not None
        and coverage >= thresholds.min_positive_coverage
        and program >= thresholds.min_positive_coverage
        and violation <= thresholds.max_negative_violation
        and as_bool(row["required_program_evidence_complete"])
        and int(float(row["annotation_lanes"])) >= 1
    ):
        return label, "the candidate fails its own coverage or negative-marker gates"
    cited = {gene for item in proposals for gene in item.cited_genes}
    uncited = sorted(cited - set(entry["top_markers"]))
    if not cited or uncited:
        return (
            label,
            f"cited genes are not positive markers of this cluster: {uncited or 'none cited'}",
        )
    return label, ""


def _schema(entry: dict[str, Any]) -> dict[str, Any]:
    labels = [item["label"] for item in entry["candidates"]] + ["Unknown"]
    return {
        "type": "object",
        "additionalProperties": False,
        "required": ["label", "cited_genes", "rationale"],
        "properties": {
            "label": {"type": "string", "enum": labels},
            "cited_genes": {"type": "array", "items": {"type": "string"}},
            "rationale": {"type": "string"},
        },
    }


def _prompt(context: dict[str, Any], entry: dict[str, Any]) -> str:
    return (
        "You are breaking a tie in single-cell cell type annotation. The measured "
        "evidence admits the candidates below and cannot separate them. Choose one, or "
        "Unknown if the markers do not support any. Cite only genes from top_markers "
        "that support your choice.\n\n"
        f"Study: {json.dumps(context, sort_keys=True)}\n"
        f"Cluster: {json.dumps(entry, sort_keys=True)}"
    )


def propose(
    config: CellCuratorConfig,
    packet: dict[str, Any],
    *,
    base_url: str,
    model: str,
    out_dir: Path,
    api_key_env: str = "",
    timeout: float = 120.0,
) -> ProposalSet:
    """Ask one OpenAI-compatible endpoint (Ollama, vLLM, llama.cpp, OpenAI, ...).

    The request goes through net.py under ``knowledge.network``, follows no
    redirect, and is sent at temperature 0 with a JSON schema. The exact request
    body and raw response are stored under ``out_dir``; the API key is read from
    the environment and appears only in the Authorization header.
    """

    url = base_url.rstrip("/") + "/chat/completions"
    headers = {"Content-Type": "application/json", "User-Agent": "cell-curator"}
    if api_key_env:
        key = os.environ.get(api_key_env, "")
        if not key:
            raise ContractError(f"environment variable {api_key_env} is empty")
        headers["Authorization"] = f"Bearer {key}"
    out_dir.mkdir(parents=True, exist_ok=True)
    proposals = []
    for entry in packet["clusters"]:
        body = _canonical(
            {
                "model": model,
                "temperature": 0,
                "messages": [{"role": "user", "content": _prompt(packet["context"], entry)}],
                "response_format": {
                    "type": "json_schema",
                    "json_schema": {
                        "name": "tie_break",
                        "strict": True,
                        "schema": _schema(entry),
                    },
                },
            }
        )
        request = urllib.request.Request(url, data=body, headers=headers, method="POST")
        with urlopen(
            request, timeout=timeout, policy=config.knowledge.policy, redirects=False
        ) as response:
            raw = response.read()
        stem = out_dir / f"cluster-{entry['cluster_id']}"
        stem.with_suffix(".request.json").write_bytes(body)
        stem.with_suffix(".response.json").write_bytes(raw)
        try:
            answer = json.loads(json.loads(raw)["choices"][0]["message"]["content"])
        except (KeyError, IndexError, TypeError, ValueError) as exc:
            raise ContractError(
                f"cluster {entry['cluster_id']}: the model's answer is not the schema: {exc}"
            ) from exc
        if answer.get("label") == "Unknown":
            continue
        proposals.append(
            Proposal(
                cluster_id=entry["cluster_id"],
                label=str(answer.get("label", "")),
                cl_id=next(
                    (
                        c["cl_id"]
                        for c in entry["candidates"]
                        if c["label"] == answer.get("label")
                    ),
                    "",
                ),
                cited_genes=tuple(str(gene) for gene in answer.get("cited_genes", [])),
                rationale=str(answer.get("rationale", "")),
                model=model,
            )
        )
    return ProposalSet(payload_sha256=packet_sha256(packet), proposals=proposals)


def merge(sets: list[ProposalSet]) -> ProposalSet:
    """Combine answers from several models into one set, for the agreement rule."""

    digests = {item.payload_sha256 for item in sets}
    if len(digests) != 1:
        raise ContractError("proposal sets answer different packets")
    return ProposalSet(
        payload_sha256=digests.pop(),
        proposals=[proposal for item in sets for proposal in item.proposals],
    )


def lane_bottom_up(config: CellCuratorConfig, run_root: Path) -> pd.DataFrame:
    frames = []
    for lane in config.evidence.lanes:
        from .evidence import lane_key

        path = (
            run_root
            / "evidence"
            / "lanes"
            / lane_key(lane.assay, lane.representation, lane.lane_key)
            / "bottom_up_markers.tsv"
        )
        if path.is_file():
            frames.append(pd.read_csv(path, sep="\t", keep_default_na=False))
    return (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(columns=["cluster_id", "feature", "rank", "effect"])
    )


def cluster_sizes(config: CellCuratorConfig) -> dict[str, int]:
    from .data import close_input, load_input

    loaded = load_input(config, backed=True)
    try:
        counts = loaded.root.obs[config.input.keys.canonical_parent].astype(str).value_counts()
    finally:
        close_input(loaded)
    return {str(key): int(value) for key, value in counts.items()}


def packet_from_run(config: CellCuratorConfig, run_root: Path) -> dict[str, Any]:
    """Rebuild the packet from a completed run's published tables."""

    summary = run_root / "evidence" / "summary"
    proposals = pd.read_csv(
        summary / "annotation_proposals.tsv", sep="\t", keep_default_na=False, dtype=str
    )
    decisions = [
        AnnotationDecision.model_validate(
            {
                **row,
                "cl_id": row["cl_id"] or None,
                "parent_label": row.get("parent_label") or None,
                "abstained": as_bool(row["abstained"]),
                "required_capabilities_complete": as_bool(
                    row["required_capabilities_complete"]
                ),
                "rejected_alternatives": (),
                "schema_version": 3,
            }
        )
        for row in proposals.to_dict("records")
    ]
    aggregate = pd.read_csv(
        summary / "multimodal_label_scores.tsv", sep="\t", keep_default_na=False
    )
    return build_packet(
        config, decisions, aggregate, lane_bottom_up(config, run_root), cluster_sizes(config)
    )


def consume(
    config: CellCuratorConfig,
    run_root: Path,
    decisions: list[AnnotationDecision],
    aggregate: pd.DataFrame,
) -> list[AnnotationDecision]:
    """Apply the plan-bound proposals file, if the config names one."""

    path = config.evidence.llm.proposals
    if not path:
        return decisions
    plan = json.loads((run_root / "plan.json").read_text())
    if plan.get("llm_proposals_sha256") != sha256(Path(path)):
        raise ContractError(
            "the LLM proposals file changed after plan bound it; start a new run.run_id"
        )
    packet = build_packet(
        config, decisions, aggregate, lane_bottom_up(config, run_root), cluster_sizes(config)
    )
    proposal_set = ProposalSet.model_validate_json(Path(path).read_text())
    updated, outcomes = apply_proposals(config, packet, proposal_set, decisions, aggregate)
    from .provenance import replace_json

    replace_json(run_root / "llm" / "packet.json", packet)
    replace_json(
        run_root / "llm" / "receipt.json",
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_llm_receipt",
            "payload_sha256": packet_sha256(packet),
            "proposals_sha256": sha256(Path(path)),
            "template": TEMPLATE_VERSION,
            "outcomes": outcomes,
        },
    )
    return updated

"""Bounded imports of saved, normalized source records; imports never approve labels.

Input is a JSON list or ``{"labels": [...], "context": {...}}``. CAP and Allen
adapters accept the documented field aliases below, not raw GraphQL or HTML.
The exact input is retained beside the snapshot, including unrecognized fields.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Literal

from pydantic import ValidationError

from .contracts import ContractError
from .corpus import (
    CorpusAssertion,
    CorpusExclusion,
    CorpusLabel,
    CorpusSnapshot,
    CorpusSource,
    load_corpus,
)
from .provenance import atomic_publish, publish_json


def _strings(value: Any) -> tuple[str, ...]:
    if value is None or value == "":
        return ()
    if isinstance(value, str):
        return tuple(dict.fromkeys(item.strip() for item in re.split(r"[,;|]", value) if item.strip()))
    if isinstance(value, list) and all(isinstance(item, str) and item.strip() for item in value):
        return tuple(dict.fromkeys(item.strip() for item in value))
    raise ValueError("marker and scope fields must be strings or lists of strings")


def _first(row: Mapping[str, Any], *keys: str, default: Any = "") -> Any:
    for key in keys:
        if key in row:
            return row[key]
    return default


def _record(
    row: Mapping[str, Any], index: int, *, source_id: str, source_kind: str,
    context: Mapping[str, Any],
) -> tuple[CorpusLabel, list[str]]:
    name = _first(row, "name", "cell_label", "label")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("source record requires name, cell_label, or label")
    local_id = str(_first(row, "source_record_id", "label_id", "id", default=index))
    label_id = local_id if local_id.startswith(source_id + ":") else f"{source_id}:{local_id}"
    description = _first(row, "description", "rationale")
    omissions = []
    if not description:
        description = f"Imported source label {name}; biological definition was not supplied."
        omissions.append("description")
    positive = _strings(_first(row, "positive", "canonical_marker_genes", "markers", default=[]))
    negative = _strings(row.get("negative", []))
    if not positive:
        omissions.append("canonical_positive_markers")
    if "negative" not in row:
        omissions.append("negative_markers_unreported")
    source_level = str(_first(row, "source_level", "level", default=""))
    if not source_level:
        omissions.append("source_level")
    canonical_level = row.get("canonical_level", "L1")
    if "canonical_level" not in row:
        omissions.append("canonical_level_unreviewed")
    scope = {}
    for key, field in (
        ("species", "species"), ("tissue", "tissues"), ("assay", "assays"),
        ("stage", "stages"), ("condition", "conditions"),
    ):
        scope[field] = _strings(row.get(field, row.get(key, context.get(key, []))))
    if not scope["species"]:
        omissions.append("species")
    measurement = row.get("measurement", "gene_expression")
    assertions = []
    reported_markers = _strings(row.get("marker_gene_evidence", []))
    if reported_markers:
        assertions.append(CorpusAssertion(
            assertion_id=f"{label_id}:reported-markers", source_id=source_id,
            kind="reported_marker", features=reported_markers,
            statement="Source-reported marker evidence; not recomputed by this import.",
            measurement=measurement,
        ))
    cl_id = _first(row, "cl_id", "cell_ontology_term_id", "ontology_term_id", default=None)
    relation = row.get("ontology_relation", "unresolved")
    if cl_id and "ontology_relation" not in row:
        omissions.append("ontology_relation_unreviewed")
    limitations = [*_strings(row.get("limitations", [])),
        f"Imported {source_kind} record; source labels and levels require review before use."]
    label = CorpusLabel(
        label_id=label_id, name=name, full_name=_first(row, "full_name", "cell_fullname"),
        description=description, source_id=source_id, source_label=row.get("source_label", name),
        source_record_id=local_id, source_parent_id=str(row.get("parent_id") or ""),
        source_level=source_level, level=canonical_level,
        axis=row.get("axis", "cell_type"), positive=positive, negative=negative,
        synonyms=_strings(row.get("synonyms", [])), cl_id=cl_id, ontology_relation=relation,
        measurement=measurement, feature_space=row.get("feature_space", "gene"),
        status="imported", limitations=tuple(limitations), assertions=tuple(assertions),
        exclusions=tuple(CorpusExclusion.model_validate(item) for item in row.get("exclusions", [])),
        species=scope["species"], tissues=scope["tissues"], assays=scope["assays"],
        stages=scope["stages"], conditions=scope["conditions"],
    )
    return label, omissions


def import_corpus_records(
    base_path: str | Path,
    input_path: str | Path,
    output_path: str | Path,
    *,
    source_kind: Literal["cap", "allen", "generic", "scgate"],
    source_id: str,
    name: str,
    url: str,
    version: str,
    retrieved_at: str,
    license: str = "",
    max_records: int = 10000,
    max_bytes: int = 10 * 1024 * 1024,
    scgate: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Import a saved source into a new immutable snapshot with an omission receipt."""
    if source_kind not in {"cap", "allen", "generic", "scgate"}:
        raise ContractError(f"unsupported corpus source kind: {source_kind}")
    if max_bytes < 1 or max_records < 1:
        raise ContractError("corpus import limits must be positive")
    base_path, input_path, output_path = Path(base_path), Path(input_path), Path(output_path)
    source_output = output_path.with_suffix(output_path.suffix + ".source.json")
    receipt_output = output_path.with_suffix(output_path.suffix + ".import-receipt.json")
    if any(path.resolve() in {base_path.resolve(), input_path.resolve()}
           for path in (output_path, source_output, receipt_output)):
        raise ContractError("corpus import cannot overwrite its base or source")
    base = load_corpus(base_path)
    if source_kind == "scgate":
        # `input_path` is a local scGate_models clone; the converted records are the retained source.
        from .corpus_scgate import scgate_payload

        options = dict(scgate or {})
        raw = scgate_payload(input_path, list(options.get("collections", [])),
                             species=str(options.get("species", "")), tissue=str(options.get("tissue", "")),
                             tag=str(options.get("tag", "")))
    else:
        if input_path.stat().st_size > max_bytes:
            raise ContractError("saved corpus source exceeds max_bytes")
        with input_path.open("rb") as handle:
            raw = handle.read(max_bytes + 1)
    if len(raw) > max_bytes:
        raise ContractError("saved corpus source exceeds max_bytes")
    try:
        payload = json.loads(raw)
    except (ValueError, UnicodeDecodeError) as exc:
        raise ContractError("saved corpus source must be UTF-8 JSON") from exc
    context = payload.get("context", {}) if isinstance(payload, dict) else {}
    rows = payload.get("labels") if isinstance(payload, dict) else payload
    if not isinstance(rows, list) or not isinstance(context, dict):
        raise ContractError("saved corpus source requires a label list and optional context mapping")
    if not rows or len(rows) > max_records:
        raise ContractError("saved corpus source has zero records or exceeds max_records")
    source = CorpusSource(
        source_id=source_id, name=name, url=url, version=version,
        sha256=hashlib.sha256(raw).hexdigest(), retrieved_at=retrieved_at, license=license,
        limitations=(
            ("Converted from a local scGate_models clone; the repository has no licence file, so nothing "
             "from it is shipped with this package.",)
            if source_kind == "scgate"
            else ("Saved normalized excerpt; not a complete atlas or live API response.",)
        ),
    )
    existing = {item.source_id: item for item in base.sources}
    if source_id in existing and existing[source_id] != source:
        raise ContractError("source_id already identifies different content; use a new source_id")
    labels: list[CorpusLabel] = []
    omissions: list[dict[str, Any]] = []
    failed: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        try:
            if not isinstance(row, dict):
                raise ValueError("source label must be an object")
            label, missing = _record(row, index, source_id=source_id, source_kind=source_kind, context=context)
            labels.append(label)
            if missing:
                omissions.append({"record_index": index, "label_id": label.label_id, "fields": missing})
        except (ValueError, TypeError) as exc:
            failed.append({"record_index": index, "reason": str(exc)})
    if not labels:
        raise ContractError(f"saved corpus source contains no valid labels: {failed}")
    all_labels = {item.label_id: item for item in base.labels}
    ids = [item.label_id for item in labels]
    if len(set(ids)) != len(ids) or set(ids).intersection(all_labels):
        raise ContractError("imported corpus label IDs repeat or already exist")
    imported = {item.label_id: item for item in labels}
    # Parentage is preserved separately until the source hierarchy is explicitly mapped.
    for label in labels:
        if label.source_parent_id:
            omissions.append({"label_id": label.label_id, "fields": ["canonical_parent_unreviewed"],
                              "source_parent_id": label.source_parent_id})
    try:
        revised = CorpusSnapshot(
            corpus_id=base.corpus_id,
            revision=f"{base.revision}+{source_id}-{source.sha256[:12]}",
            parent_revision=base.revision,
            sources=tuple(existing[key] for key in sorted(existing) if key != source_id) + (source,),
            labels=tuple(all_labels[key] for key in sorted(all_labels)) + tuple(imported[key] for key in sorted(imported)),
            mappings=base.mappings,
        )
    except ValidationError as exc:
        raise ContractError(f"invalid corpus import: {exc}") from exc
    atomic_publish(source_output, raw)
    output_digest = publish_json(output_path, revised.model_dump(mode="json"))
    receipt = {
        "schema_version": 3, "source_kind": source_kind, "source_id": source_id,
        "source_sha256": source.sha256, "base_sha256": base.digest,
        "output_sha256": output_digest, "corpus_digest": revised.digest,
        "status": "partial" if omissions or failed else "imported",
        "input_count": len(rows), "imported_count": len(labels), "failed_count": len(failed),
        "omissions": omissions, "failed_records": failed, "network_attempted": False,
        "source_copy": str(source_output), "review_required": True,
    }
    publish_json(receipt_output, receipt)
    return receipt

"""Versioned, source-backed label definitions and explicitly reviewed revisions."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Literal

from pydantic import ConfigDict, Field, ValidationError, field_validator, model_validator

from .contracts import ContractError, sha256
from .models import MarkerProgram, StrictModel
from .ontology import CellOntology
from .provenance import publish_json

ReviewStatus = Literal["imported", "proposed", "reviewed", "disputed", "superseded", "retired"]
Relation = Literal["exact", "broader", "narrower", "related", "conflicting", "unresolved"]
Measurement = Literal[
    "gene_expression", "protein_expression", "gene_activity", "accessibility", "transcript_isoform"
]


class _CorpusModel(StrictModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    @field_validator("*", mode="before")
    @classmethod
    def strip_strings(cls, value: Any) -> Any:
        if isinstance(value, str):
            return value.strip()
        if isinstance(value, (list, tuple)) and any(isinstance(item, str) for item in value):
            if any(isinstance(item, str) and not item.strip() for item in value):
                raise ValueError("corpus lists cannot contain blank strings")
            return tuple(item.strip() if isinstance(item, str) else item for item in value)
        return value


class CorpusSource(_CorpusModel):
    source_id: str = Field(min_length=1)
    name: str = Field(min_length=1)
    url: str = Field(min_length=1)
    version: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    retrieved_at: str = Field(min_length=1)
    license: str = ""
    depends_on: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()


class CorpusAssertion(_CorpusModel):
    """A source assertion; observed evidence stays distinguishable from prior knowledge."""

    assertion_id: str = Field(min_length=1)
    source_id: str = Field(min_length=1)
    statement: str = Field(min_length=1)
    kind: Literal[
        "definition", "canonical_marker", "observed_marker", "reported_marker", "exclusion", "threshold", "limitation"
    ]
    features: tuple[str, ...] = ()
    measurement: Measurement = "gene_expression"
    threshold: float | None = Field(default=None, allow_inf_nan=False)
    units: str = ""
    run_id: str = ""

    @model_validator(mode="after")
    def evidence_contract(self) -> CorpusAssertion:
        if self.kind == "observed_marker" and not self.run_id:
            raise ValueError("observed marker assertions require a run_id")
        if self.threshold is not None and not self.units:
            raise ValueError("a corpus threshold requires units")
        return self


class CorpusExclusion(_CorpusModel):
    """A population this label must not match, tested through its own marker genes."""

    name: str = Field(min_length=1)
    positive: tuple[str, ...] = Field(min_length=1)
    negative: tuple[str, ...] = ()


class CorpusObservation(_CorpusModel):
    run_id: str = Field(min_length=1)
    dataset_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    evidence_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    conclusion: Literal["supports", "contradicts", "unresolved"]
    validation: Literal["descriptive", "reviewed", "independent"] = "descriptive"
    sample_count: int | None = Field(default=None, ge=1)
    donor_count: int | None = Field(default=None, ge=1)
    rationale: str = Field(min_length=1)


class CorpusLabel(_CorpusModel):
    label_id: str = Field(min_length=1)
    name: str = Field(min_length=1)
    full_name: str = ""
    description: str = Field(min_length=1)
    source_id: str = Field(min_length=1)
    source_label: str = Field(min_length=1)
    source_record_id: str = ""
    source_parent_id: str = ""
    source_level: str = ""
    level: str = Field(default="L1", pattern=r"^L[1-9][0-9]*$")
    parent_id: str | None = None
    axis: Literal["cell_type", "cell_state"] = "cell_type"
    species: tuple[str, ...] = ()
    tissues: tuple[str, ...] = ()
    assays: tuple[str, ...] = ()
    stages: tuple[str, ...] = ()
    conditions: tuple[str, ...] = ()
    positive: tuple[str, ...] = ()
    negative: tuple[str, ...] = ()
    synonyms: tuple[str, ...] = ()
    cl_id: str | None = Field(default=None, pattern=r"^CL:[0-9]{7}$")
    ontology_relation: Relation = "unresolved"
    feature_space: str = "gene"
    measurement: Measurement = "gene_expression"
    status: ReviewStatus = "imported"
    reviewer: str = ""
    rationale: str = ""
    study_id: str = ""
    superseded_by: str = ""
    exclusions: tuple[CorpusExclusion, ...] = ()
    alternatives: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()
    assertions: tuple[CorpusAssertion, ...] = ()
    experience: tuple[CorpusObservation, ...] = ()

    @model_validator(mode="after")
    def label_contract(self) -> CorpusLabel:
        if self.status == "reviewed" and not (self.reviewer and self.rationale):
            raise ValueError("reviewed corpus labels require reviewer and rationale")
        if self.status == "reviewed" and not self.species:
            raise ValueError("reviewed corpus labels require an explicit species scope")
        if self.status == "superseded" and not self.superseded_by:
            raise ValueError("a superseded corpus label must name the label that replaced it")
        if self.status == "reviewed" and not self.tissues:
            raise ValueError("reviewed corpus labels require an explicit tissue scope")
        overlap = {value.casefold() for value in self.positive}.intersection(
            value.casefold() for value in self.negative
        )
        if overlap:
            raise ValueError(f"contradictory marker signs in corpus label {self.label_id}")
        if self.ontology_relation == "exact" and not self.cl_id:
            raise ValueError("an exact ontology mapping requires cl_id")
        if self.study_id and not self.parent_id:
            raise ValueError("a study-specific label requires an accepted parent")
        expected_spaces = {
            "gene_expression": "gene",
            "gene_activity": "gene",
            "protein_expression": "protein",
            "accessibility": "peak",
            "transcript_isoform": "transcript",
        }
        if self.feature_space != expected_spaces[self.measurement]:
            raise ValueError("measurement and feature_space disagree")
        assertion_ids = [item.assertion_id for item in self.assertions]
        if len(assertion_ids) != len(set(assertion_ids)):
            raise ValueError("duplicate corpus assertion IDs")
        return self


class CorpusMapping(_CorpusModel):
    source_label_id: str = Field(min_length=1)
    target_label_id: str = Field(min_length=1)
    relation: Relation
    rationale: str = Field(min_length=1)
    reviewer: str = ""


class CorpusReview(_CorpusModel):
    decision: Literal["accept", "reject", "dispute"]
    reviewer: str = Field(min_length=1)
    rationale: str = Field(min_length=1)
    proposal_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    base_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


def _check_graph(parents: Mapping[str, Sequence[str]], *, kind: str) -> None:
    complete: set[str] = set()
    for start in parents:
        pending: list[tuple[str, bool]] = [(start, False)]
        active: set[str] = set()
        while pending:
            node, expanded = pending.pop()
            if node not in parents:
                raise ValueError(f"unknown {kind} reference: {node}")
            if expanded:
                active.remove(node)
                complete.add(node)
                continue
            if node in active:
                raise ValueError(f"cycle in corpus {kind}: {node}")
            if node in complete:
                continue
            active.add(node)
            pending.append((node, True))
            pending.extend((parent, False) for parent in parents[node])


class CorpusSnapshot(_CorpusModel):
    schema_version: Literal[3] = 3
    corpus_id: str = Field(min_length=1)
    revision: str = Field(min_length=1)
    parent_revision: str | None = None
    sources: tuple[CorpusSource, ...]
    labels: tuple[CorpusLabel, ...]
    mappings: tuple[CorpusMapping, ...] = ()
    review: CorpusReview | None = None

    @model_validator(mode="after")
    def validate_records(self) -> CorpusSnapshot:
        sources = {item.source_id: item for item in self.sources}
        labels = {item.label_id: item for item in self.labels}
        if len(sources) != len(self.sources) or len(labels) != len(self.labels):
            raise ValueError("duplicate corpus source or label IDs")
        _check_graph({key: value.depends_on for key, value in sources.items()}, kind="source")
        _check_graph(
            {key: (value.parent_id,) if value.parent_id else () for key, value in labels.items()},
            kind="parent",
        )
        depths: dict[str, int] = {}
        for identifier in labels:
            chain = []
            cursor: str | None = identifier
            while cursor is not None and cursor not in depths:
                chain.append(cursor)
                cursor = labels[cursor].parent_id
            depth = depths[cursor] if cursor is not None else 0
            for node in reversed(chain):
                depth += 1
                depths[node] = depth
        for label in self.labels:
            for source_id in (label.source_id, *(item.source_id for item in label.assertions)):
                if source_id not in sources:
                    raise ValueError(f"unknown corpus source: {source_id}")
            if label.parent_id:
                parent = labels[label.parent_id]
                if parent.axis != "cell_type":
                    raise ValueError("corpus parents must be cell identities")
                if int(label.level[1:]) <= int(parent.level[1:]):
                    raise ValueError("corpus child level must follow its parent")
                if parent.study_id and parent.study_id != label.study_id:
                    raise ValueError("corpus label crosses study-specific parent scope")
            if label.status == "reviewed" and label.axis == "cell_type" and int(label.level[1:]) != depths[label.label_id]:
                raise ValueError("reviewed corpus identity level must equal its hierarchy depth")
        pairs: set[tuple[str, str]] = set()
        exact: dict[str, tuple[str, ...]] = {}
        for mapping in self.mappings:
            pair = (mapping.source_label_id, mapping.target_label_id)
            if pair in pairs or pair[0] == pair[1]:
                raise ValueError("duplicate or self-referential corpus mapping")
            pairs.add(pair)
            if any(key not in labels for key in pair):
                raise ValueError("corpus mapping refers to an unknown label")
            if mapping.relation == "exact":
                if not mapping.reviewer:
                    raise ValueError("exact corpus mappings require a reviewer")
                if pair[0] in exact:
                    raise ValueError("ambiguous exact corpus mapping")
                source, target = (labels[key] for key in pair)
                if source.axis != target.axis or source.measurement != target.measurement:
                    raise ValueError("exact corpus mapping crosses identity or measurement axes")
                for key, field in (
                    ("species", "species"), ("tissue", "tissues"),
                    ("assay", "assays"), ("stage", "stages"), ("condition", "conditions"),
                ):
                    left, right = getattr(source, field), getattr(target, field)
                    if left and right and not {_context_value(key, item) for item in left}.intersection(
                        _context_value(key, item) for item in right
                    ):
                        raise ValueError("exact corpus mapping has incompatible biological scope")
                if source.cl_id and target.cl_id and source.cl_id != target.cl_id:
                    raise ValueError("exact corpus mapping has contradictory ontology IDs")
                exact[pair[0]] = (pair[1],)
        _check_graph({key: exact.get(key, ()) for key in labels}, kind="exact mapping")
        return self

    @property
    def digest(self) -> str:
        """Hash the normalized snapshot, independent of JSON whitespace and record order."""
        value = self.model_dump(mode="json")
        value["sources"] = sorted(value["sources"], key=lambda item: item["source_id"])
        value["labels"] = sorted(value["labels"], key=lambda item: item["label_id"])
        value["mappings"] = sorted(
            value["mappings"], key=lambda item: (item["source_label_id"], item["target_label_id"])
        )
        return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()

    def select(
        self, context: Mapping[str, str] | None = None, *, callable_only: bool = True
    ) -> list[CorpusLabel]:
        """Select eligible definitions; unknown context cannot satisfy a limited scope."""
        if context is not None and not isinstance(context, Mapping):
            raise ContractError("corpus context must be a mapping")
        context = context or {}
        labels = {label.label_id: label for label in self.labels}

        def matches(label: CorpusLabel) -> bool:
            for key, values in (
                ("species", label.species), ("tissue", label.tissues),
                ("assay", label.assays), ("stage", label.stages),
                ("condition", label.conditions),
            ):
                if values and _context_value(key, str(context.get(key, ""))) not in {
                    _context_value(key, value) for value in values
                }:
                    return False
            return not label.study_id or context.get("study_id") == label.study_id

        selected = []
        for label in sorted(self.labels, key=lambda item: (int(item.level[1:]), item.label_id)):
            if not matches(label):
                continue
            if callable_only:
                ancestors = [label]
                parent = label.parent_id
                while parent:
                    ancestors.append(labels[parent])
                    parent = labels[parent].parent_id
                if any(
                    item.status != "reviewed" or not item.positive or not matches(item)
                    for item in ancestors
                ):
                    continue
                if not label.positive:
                    continue
            selected.append(label)
        return selected

    def programs(
        self,
        context: Mapping[str, str] | None = None,
        *,
        callable_only: bool = True,
        ontology: CellOntology | None = None,
        measurements: Sequence[str] | None = None,
    ) -> list[MarkerProgram]:
        """Export reviewed callable definitions to the existing signed-program engine."""
        if not callable_only:
            raise ContractError("unreviewed corpus records cannot be exported as callable programs")
        selected = {item.label_id: item for item in self.select(context)}
        all_labels = {item.label_id: item for item in self.labels}
        exact = {
            item.source_label_id: item.target_label_id
            for item in self.mappings if item.relation == "exact"
        }
        merged: dict[str, list[CorpusLabel]] = {}
        for label in selected.values():
            target = label.label_id
            while target in exact:
                target = exact[target]
            if target not in selected:
                raise ContractError("an exact corpus target is not callable in this context")
            merged.setdefault(target, []).append(label)
        programs = []
        names: set[tuple[str, str, str]] = set()
        aliases: dict[tuple[str, str, str], str] = {}
        for target, records in sorted(merged.items()):
            label = selected[target]
            if measurements is not None and label.measurement not in measurements:
                raise ContractError(f"required corpus measurement is unavailable: {label.measurement}")
            if label.measurement == "transcript_isoform":
                raise ContractError("transcript-isoform corpus rules have no executable marker backend")
            parent = all_labels.get(label.parent_id or "")
            key = (label.axis, parent.name if parent else "", label.name)
            if key in names:
                raise ContractError("ambiguous corpus labels require an explicit exact mapping")
            names.add(key)
            for name in (label.name, *label.synonyms):
                alias_key = (label.axis, parent.name if parent else "", name.casefold())
                if alias_key in aliases and aliases[alias_key] != target:
                    raise ContractError("ambiguous corpus alias requires an explicit exact mapping")
                aliases[alias_key] = target
            accession = label.cl_id if label.ontology_relation == "exact" and ontology else None
            if accession and ontology:
                ontology.validate_term(accession, expected_name=label.full_name or label.name)
            parent_accession = (
                parent.cl_id if accession and parent and parent.ontology_relation == "exact" else None
            )
            if parent_accession and parent and ontology and accession:
                ontology.validate_term(parent_accession, expected_name=parent.full_name or parent.name)
                if not ontology.is_descendant(accession, parent_accession):
                    raise ContractError("corpus ontology mapping contradicts its parent")
            # Decision row 9: sources mapped to one label share only the positives they all
            # list. Each record keeps its own list, so the source's view stays inspectable.
            shared = set.intersection(*({gene.casefold() for gene in item.positive} for item in records))
            positive = {gene for item in records for gene in item.positive if gene.casefold() in shared}
            if not positive:
                raise ContractError(f"exact corpus mappings share no positive marker: {target}")
            negative = {gene for item in records for gene in item.negative}
            if {gene.casefold() for gene in positive}.intersection(gene.casefold() for gene in negative):
                raise ContractError(f"conflicting reviewed corpus definitions: {target}")
            try:
                programs.append(MarkerProgram(
                    label=label.name, axis=label.axis,
                    positive=tuple(sorted(positive)),
                    negative=tuple(sorted(negative)),
                    parent=parent.name if parent else None,
                    parent_cl_id=parent_accession, cl_id=accession,
                    feature_space=label.feature_space,
                    source=f"corpus:{self.corpus_id}:{self.digest}:{target}",
                    source_version=self.revision,
                ))
            except ValidationError as exc:
                raise ContractError(f"conflicting reviewed corpus definitions: {target}: {exc}") from exc
        return programs


def _context_value(key: str, value: str) -> str:
    normalized = value.strip().casefold()
    if key == "species":
        return {
            "human": "homo sapiens", "9606": "homo sapiens", "ncbitaxon:9606": "homo sapiens",
            "mouse": "mus musculus", "10090": "mus musculus", "ncbitaxon:10090": "mus musculus",
        }.get(normalized, normalized)
    return normalized


def load_corpus(path: str | Path, expected_sha256: str = "") -> CorpusSnapshot:
    path = Path(path)
    if not path.is_file():
        raise ContractError(f"corpus source does not exist: {path}")
    if expected_sha256 and sha256(path) != expected_sha256:
        raise ContractError(f"corpus checksum mismatch: {path}")
    try:
        return CorpusSnapshot.model_validate_json(path.read_text())
    except (ValueError, OSError) as exc:
        raise ContractError(f"invalid corpus snapshot: {exc}") from exc


class CorpusProposal(_CorpusModel):
    schema_version: Literal[3] = 3
    corpus_id: str
    base_revision: str
    base_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    proposer: str = Field(min_length=1)
    reason: str = Field(min_length=1)
    labels: tuple[CorpusLabel, ...] = ()
    sources: tuple[CorpusSource, ...] = ()
    mappings: tuple[CorpusMapping, ...] = ()

    @model_validator(mode="after")
    def unique_changes(self) -> CorpusProposal:
        for values in (
            [item.label_id for item in self.labels],
            [item.source_id for item in self.sources],
            [(item.source_label_id, item.target_label_id) for item in self.mappings],
        ):
            if len(values) != len(set(values)):
                raise ValueError("duplicate changes in corpus proposal")
        return self


def propose_corpus(
    base_path: str | Path,
    proposal_path: str | Path,
    *,
    labels: Sequence[CorpusLabel | Mapping[str, Any]] = (),
    sources: Sequence[CorpusSource | Mapping[str, Any]] = (),
    mappings: Sequence[CorpusMapping | Mapping[str, Any]] = (),
    reason: str,
    proposer: str,
) -> CorpusProposal:
    if Path(base_path).resolve() == Path(proposal_path).resolve():
        raise ContractError("corpus proposal cannot overwrite its base")
    base = load_corpus(base_path)
    proposed_labels = tuple(
        CorpusLabel.model_validate({
            **(item.model_dump() if isinstance(item, CorpusLabel) else dict(item)),
            "status": "proposed", "reviewer": "", "rationale": "",
        }) for item in labels
    )
    proposal = CorpusProposal(
        corpus_id=base.corpus_id, base_revision=base.revision, base_sha256=base.digest,
        proposer=proposer, reason=reason, labels=proposed_labels,
        sources=tuple(CorpusSource.model_validate(item) for item in sources),
        mappings=tuple(CorpusMapping.model_validate(item) for item in mappings),
    )
    if not (proposal.labels or proposal.sources or proposal.mappings):
        raise ContractError("an empty corpus proposal changes nothing")
    _merge_proposal(base, proposal, revision="proposal-validation")
    publish_json(Path(proposal_path), proposal.model_dump(mode="json"))
    return proposal


def _merge_proposal(
    base: CorpusSnapshot, proposal: CorpusProposal, *, revision: str,
    review: CorpusReview | None = None,
) -> CorpusSnapshot:
    def merge(old: Sequence[Any], new: Sequence[Any], key: str) -> list[Any]:
        rows = {getattr(item, key): item for item in old}
        added = [getattr(item, key) for item in new]
        if len(added) != len(set(added)):
            raise ContractError(f"duplicate {key} in corpus proposal")
        rows.update({getattr(item, key): item for item in new})
        return [rows[item] for item in sorted(rows)]

    def retained(old: Sequence[CorpusLabel], new: Sequence[CorpusLabel]) -> list[CorpusLabel]:
        # Decision row 20: a replaced record stays, marked superseded and linked to its
        # successor, under an id that names the revision it came from.
        replaced = {item.label_id for item in new}
        return [
            item.model_copy(update={
                "label_id": f"{item.label_id}~{base.revision}", "status": "superseded",
                "superseded_by": item.label_id,
            })
            for item in old if item.label_id in replaced and item.status != "superseded"
        ]

    labels = proposal.labels
    if review:
        labels = tuple(CorpusLabel.model_validate({
            **item.model_dump(),
            "status": "reviewed" if review.decision == "accept" else "disputed",
            "reviewer": review.reviewer, "rationale": review.rationale,
        }) for item in labels)
    old_mappings = {(item.source_label_id, item.target_label_id): item for item in base.mappings}
    for item in proposal.mappings:
        old_mappings[(item.source_label_id, item.target_label_id)] = item
    try:
        return CorpusSnapshot(
            corpus_id=base.corpus_id, revision=revision, parent_revision=base.revision,
            sources=tuple(merge(base.sources, proposal.sources, "source_id")),
            labels=tuple(merge([*base.labels, *retained(base.labels, labels)], labels, "label_id")),
            mappings=tuple(old_mappings[key] for key in sorted(old_mappings)), review=review,
        )
    except ValueError as exc:
        raise ContractError(f"invalid proposed corpus revision: {exc}") from exc


def review_corpus(
    base_path: str | Path,
    proposal_path: str | Path,
    output_path: str | Path,
    *,
    reviewer: str,
    rationale: str,
    decision: Literal["accept", "reject", "dispute"],
    revision: str,
) -> dict[str, Any]:
    """Apply an explicit reviewer decision into a new immutable revision and receipt."""
    base = load_corpus(base_path)
    proposal_path, output_path = Path(proposal_path), Path(output_path)
    receipt_path = output_path.with_suffix(output_path.suffix + ".receipt.json")
    if any(path.resolve() in {Path(base_path).resolve(), proposal_path.resolve()}
           for path in (output_path, receipt_path)):
        raise ContractError("corpus review cannot overwrite its base or proposal")
    if not revision.strip() or revision == base.revision:
        raise ContractError("corpus review requires a new nonblank revision")
    try:
        proposal = CorpusProposal.model_validate_json(proposal_path.read_text())
        review = CorpusReview(
            decision=decision, reviewer=reviewer, rationale=rationale,
            proposal_sha256=sha256(proposal_path), base_sha256=base.digest,
        )
    except (ValueError, OSError) as exc:
        raise ContractError(f"invalid corpus review: {exc}") from exc
    if (proposal.corpus_id, proposal.base_revision, proposal.base_sha256) != (
        base.corpus_id, base.revision, base.digest,
    ):
        raise ContractError("corpus proposal belongs to a different or changed base")
    if decision == "reject":
        revised = CorpusSnapshot.model_validate({
            **base.model_dump(), "revision": revision,
            "parent_revision": base.revision, "review": review.model_dump(),
        })
    else:
        revised = _merge_proposal(base, proposal, revision=revision, review=review)
    output_digest = publish_json(output_path, revised.model_dump(mode="json"))
    receipt = {
        "schema_version": 3, "corpus_id": base.corpus_id,
        "base_revision": base.revision, "revision": revision,
        **review.model_dump(mode="json"), "output_sha256": output_digest,
        "corpus_digest": revised.digest,
        "label_changes": [],
    }
    original = {item.label_id: item for item in base.labels}
    for label in revised.labels:
        previous = original.get(label.label_id)
        if previous == label:
            continue
        receipt["label_changes"].append({
            "label_id": label.label_id, "change": "changed" if previous else "added",
            "fields": sorted(
                key for key, value in label.model_dump().items()
                if previous is None or previous.model_dump()[key] != value
            ),
        })
    publish_json(receipt_path, receipt)
    return receipt

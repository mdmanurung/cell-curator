"""Local Cell Ontology parsing and semantic hierarchy validation."""

from __future__ import annotations

import csv
import json
import re
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from .contracts import ContractError, sha256

if TYPE_CHECKING:
    from .config import CellCuratorConfig
    from .models import MarkerProgram


@dataclass(frozen=True)
class OntologyTerm:
    accession: str
    name: str
    parents: tuple[str, ...] = ()
    synonyms: tuple[str, ...] = ()
    obsolete: bool = False
    replaced_by: str | None = None


class CellOntology:
    def __init__(self, terms: dict[str, OntologyTerm], *, source: str) -> None:
        self.terms = terms
        self.source = source
        self._names: dict[str, set[str]] = {}
        for accession, term in terms.items():
            for name in (term.name, *term.synonyms):
                self._names.setdefault(name.casefold(), set()).add(accession)

    @classmethod
    def load(cls, path: str | Path) -> CellOntology:
        source = Path(path)
        if not source.is_file():
            raise ContractError(f"ontology source does not exist: {source}")
        suffix = source.suffix.lower()
        if suffix == ".obo":
            terms = cls._parse_obo(source)
        elif suffix in {".owl", ".rdf", ".xml"}:
            terms = cls._parse_owl(source)
        elif suffix in {".tsv", ".csv"}:
            terms = cls._parse_table(source)
        else:
            raise ContractError("ontology source must be OBO, OWL/RDF/XML, TSV, or CSV")
        if not terms:
            raise ContractError("ontology source contains no terms")
        return cls(terms, source=str(source.resolve()))

    @staticmethod
    def _parse_obo(path: Path) -> dict[str, OntologyTerm]:
        terms: dict[str, OntologyTerm] = {}
        current: dict[str, Any] | None = None

        def commit() -> None:
            nonlocal current
            if current and current.get("id") and current.get("name"):
                accession = str(current["id"])
                terms[accession] = OntologyTerm(
                    accession=accession,
                    name=str(current["name"]),
                    parents=tuple(map(str, current.get("parents", []))),
                    synonyms=tuple(map(str, current.get("synonyms", []))),
                    obsolete=bool(current.get("obsolete", False)),
                    replaced_by=(
                        str(current["replaced_by"]) if current.get("replaced_by") else None
                    ),
                )
            current = None

        for raw in path.read_text().splitlines():
            line = raw.strip()
            if line == "[Term]":
                commit()
                current = {"parents": [], "synonyms": []}
                continue
            if line.startswith("["):
                commit()
                continue
            if current is None or ": " not in line:
                continue
            key, value = line.split(": ", 1)
            if key in {"id", "name", "replaced_by"}:
                current[key] = value.split(" ! ", 1)[0]
            elif key == "is_a":
                current["parents"].append(value.split(" ! ", 1)[0])
            elif key == "synonym":
                match = re.match(r'"([^"]+)"', value)
                if match:
                    current["synonyms"].append(match.group(1))
            elif key == "is_obsolete":
                current["obsolete"] = value.lower() == "true"
        commit()
        return terms

    @staticmethod
    def _parse_table(path: Path) -> dict[str, OntologyTerm]:
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter=delimiter))
        terms = {}
        for row in rows:
            accession = str(row.get("id") or row.get("accession") or "").strip()
            name = str(row.get("name") or "").strip()
            if not accession or not name:
                raise ContractError("ontology table rows require id/accession and name")

            def split(value: str | None) -> tuple[str, ...]:
                return tuple(
                    item.strip() for item in re.split(r"[;|]", value or "") if item.strip()
                )

            terms[accession] = OntologyTerm(
                accession=accession,
                name=name,
                parents=split(row.get("parents", "")),
                synonyms=split(row.get("synonyms", "")),
                obsolete=str(row.get("obsolete", "")).lower() in {"true", "1", "yes"},
                replaced_by=str(row.get("replaced_by") or "").strip() or None,
            )
        return terms

    @staticmethod
    def _parse_owl(path: Path) -> dict[str, OntologyTerm]:
        root = ET.parse(path).getroot()
        rdf_about = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"
        rdf_resource = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"
        terms: dict[str, OntologyTerm] = {}
        for element in root.iter():
            about = element.attrib.get(rdf_about, "")
            if not about or "CL_" not in about:
                continue
            accession = about.rsplit("/", 1)[-1].replace("CL_", "CL:")
            name = ""
            parents: list[str] = []
            synonyms: list[str] = []
            obsolete = False
            replaced_by = None
            for child in element:
                local = child.tag.rsplit("}", 1)[-1]
                if local in {"label", "prefLabel"} and child.text:
                    name = child.text.strip()
                elif local == "subClassOf":
                    resource = child.attrib.get(rdf_resource, "")
                    if "CL_" in resource:
                        parents.append(resource.rsplit("/", 1)[-1].replace("CL_", "CL:"))
                elif "synonym" in local.lower() and child.text:
                    synonyms.append(child.text.strip())
                elif local == "deprecated" and child.text:
                    obsolete = child.text.strip().lower() == "true"
                elif local == "replacedBy":
                    resource = child.attrib.get(rdf_resource, "")
                    if "CL_" in resource:
                        replaced_by = resource.rsplit("/", 1)[-1].replace("CL_", "CL:")
            if name:
                terms[accession] = OntologyTerm(
                    accession, name, tuple(parents), tuple(synonyms), obsolete, replaced_by
                )
        return terms

    def resolve(self, value: str) -> OntologyTerm:
        value = value.strip()
        if value in self.terms:
            return self.terms[value]
        matches = self._names.get(value.casefold(), set())
        if not matches:
            raise ContractError(f"ontology term is unknown: {value}")
        if len(matches) > 1:
            raise ContractError(f"ontology name is ambiguous: {value} -> {sorted(matches)}")
        return self.terms[next(iter(matches))]

    def ancestors(self, accession: str) -> set[str]:
        term = self.resolve(accession)
        seen: set[str] = set()
        stack = list(term.parents)
        while stack:
            parent = stack.pop()
            if parent in seen:
                continue
            seen.add(parent)
            if parent in self.terms:
                stack.extend(self.terms[parent].parents)
        return seen

    def validate_term(
        self, accession: str, *, expected_name: str | None = None
    ) -> OntologyTerm:
        term = self.resolve(accession)
        if term.obsolete:
            replacement = f"; replaced by {term.replaced_by}" if term.replaced_by else ""
            raise ContractError(f"ontology term is obsolete: {accession}{replacement}")
        if expected_name:
            normalized = {term.name.casefold(), *(item.casefold() for item in term.synonyms)}
            if expected_name.casefold() not in normalized:
                raise ContractError(
                    f"label {expected_name!r} is inconsistent with {accession} ({term.name!r})"
                )
        return term

    def is_descendant(self, child: str, parent: str) -> bool:
        child_term = self.validate_term(child)
        parent_term = self.validate_term(parent)
        return (
            child_term.accession == parent_term.accession
            or parent_term.accession in self.ancestors(child_term.accession)
        )


def ontology_from_config(config: CellCuratorConfig) -> CellOntology | None:
    providers = [
        item
        for item in config.knowledge.providers
        if item.kind in {"ontology", "remote_ontology"}
    ]
    if len(providers) > 1:
        raise ContractError("configure at most one authoritative ontology provider")
    if not providers:
        return None
    provider = providers[0]
    if provider.kind == "ontology":
        return CellOntology.load(provider.path)
    from .knowledge import RemoteTableProvider
    from .provenance import replace_json

    status_path = (
        Path(config.run.output_root)
        / config.run.run_id
        / "knowledge"
        / "ontology_provider_status.json"
    )
    cache_path = Path(provider.cache_path)
    if not config.knowledge.allow_optional_remote_enhancement and not cache_path.is_file():
        replace_json(
            status_path,
            {
                "schema_version": 3,
                "artifact_kind": "cell_curator_optional_ontology_status",
                "status": "unavailable",
                "reason": (
                    "optional remote enhancement is disabled and the checksum-pinned "
                    f"ontology cache is absent: {cache_path}"
                ),
                "remediation": (
                    "Provide the checksum-pinned ontology cache locally or authorize remote "
                    "enhancement; CL-free annotation remains available."
                ),
                "network_attempted": False,
            },
        )
        return None
    try:
        cache = RemoteTableProvider(
            provider.url,
            provider.cache_path,
            expected_sha256=provider.sha256,
            version=provider.version,
            timeout_seconds=provider.timeout_seconds,
            allow_network=config.knowledge.allow_optional_remote_enhancement,
        ).provision()
        ontology = CellOntology.load(cache)
    except (ContractError, OSError) as exc:
        status = {
            "schema_version": 3,
            "artifact_kind": "cell_curator_optional_ontology_status",
            "status": "unavailable",
            "reason": str(exc),
            "remediation": (
                "Provide the checksum-pinned ontology cache locally or restore authorized "
                "network access; CL-free annotation remains available."
            ),
        }
        replace_json(status_path, status)
        return None
    replace_json(
        status_path,
        {
            "schema_version": 3,
            "artifact_kind": "cell_curator_optional_ontology_status",
            "status": "complete",
            "cache_path": str(cache),
            "sha256": sha256(cache),
            "version": provider.version,
            "n_terms": len(ontology.terms),
        },
    )
    return ontology


def validate_program_ontology(programs: list[MarkerProgram], ontology: CellOntology) -> None:
    for program in programs:
        if not program.cl_id:
            if program.parent_cl_id:
                raise ContractError(
                    f"marker program {program.label!r} has parent_cl_id without cl_id"
                )
            continue
        ontology.validate_term(program.cl_id, expected_name=program.label)
        if program.parent_cl_id:
            ontology.validate_term(program.parent_cl_id)
            if not ontology.is_descendant(program.cl_id, program.parent_cl_id):
                raise ContractError(
                    f"marker program {program.label!r} is not a descendant of "
                    f"{program.parent_cl_id}"
                )


def validate_durable_label_ontology(
    labels: pd.DataFrame,
    ontology: CellOntology | None,
) -> dict[str, Any]:
    """Semantically validate every supplied CL mapping in durable labels tables.

    Blank CL identifiers remain valid CL-free annotation. Once an identifier is
    supplied, its term, name/synonym, obsolescence, and parent compatibility are
    all enforced against the same local ontology snapshot.
    """

    required = {"level", "parent_scope", "cluster_id", "label"}
    missing = required.difference(labels.columns)
    if missing:
        raise ContractError(f"durable labels miss ontology fields: {sorted(missing)}")
    frame = labels.copy()
    if "cl_id" not in frame:
        frame["cl_id"] = ""
    if "parent_cl_id" not in frame:
        frame["parent_cl_id"] = ""
    for column in ("label", "parent_scope", "cl_id", "parent_cl_id"):
        frame[column] = frame[column].astype(str).str.strip()
    for column in ("cl_id", "parent_cl_id"):
        frame.loc[
            frame[column].str.casefold().isin({"none", "nan", "na", "<na>"}),
            column,
        ] = ""
    supplied = frame["cl_id"].ne("")
    supplied_parent = frame["parent_cl_id"].ne("")
    if supplied_parent.any() and (~supplied & supplied_parent).any():
        raise ContractError("durable labels cannot supply parent_cl_id without cl_id")
    if not supplied.any() and not supplied_parent.any():
        return {
            "status": "complete",
            "ontology_status": "unavailable-cl-free"
            if ontology is None
            else "available-cl-free",
            "n_supplied_cl_ids": 0,
            "n_validated_parent_links": 0,
        }
    if ontology is None:
        raise ContractError("CL IDs cannot be written without a validated local ontology")
    unknown_with_id = frame.loc[supplied & frame["label"].eq("Unknown")]
    if not unknown_with_id.empty:
        raise ContractError("Unknown durable labels cannot carry a biological CL identifier")

    identifiers_by_label = (
        frame.loc[supplied]
        .groupby("label", observed=True)["cl_id"]
        .agg(lambda values: sorted(set(values)))
    )
    inconsistent = {
        str(label): values for label, values in identifiers_by_label.items() if len(values) != 1
    }
    if inconsistent:
        raise ContractError(
            f"durable labels map one label name to multiple CL identifiers: {inconsistent}"
        )
    label_to_accession = {
        str(label): str(values[0]) for label, values in identifiers_by_label.items()
    }
    for label, accession in label_to_accession.items():
        ontology.validate_term(accession, expected_name=label)

    validated_parent_links = 0
    for row in frame.loc[supplied].itertuples(index=False):
        child = str(row.cl_id)
        parent_scope = str(row.parent_scope)
        explicit_parent = str(row.parent_cl_id)
        inferred_parent = label_to_accession.get(parent_scope, "")
        if explicit_parent:
            ontology.validate_term(explicit_parent)
            if inferred_parent and explicit_parent != inferred_parent:
                raise ContractError(
                    f"durable label {row.label!r} supplies parent {explicit_parent}, but "
                    f"parent scope {parent_scope!r} maps to {inferred_parent}"
                )
        parent = explicit_parent or inferred_parent
        if parent_scope.casefold() == "root":
            if explicit_parent:
                raise ContractError("ROOT-scoped durable labels cannot supply parent_cl_id")
            continue
        if parent:
            if not ontology.is_descendant(child, parent):
                raise ContractError(
                    f"durable label {row.label!r} ({child}) is not compatible with "
                    f"parent scope {parent_scope!r} ({parent})"
                )
            validated_parent_links += 1
    return {
        "status": "complete",
        "ontology_status": "validated",
        "ontology_source": ontology.source,
        "n_supplied_cl_ids": int(supplied.sum()),
        "n_unique_cl_ids": len(set(frame.loc[supplied, "cl_id"])),
        "n_validated_parent_links": validated_parent_links,
    }


def validate_terms(
    ontology_path: str | Path,
    terms: Sequence[str],
    *,
    expected_names: Sequence[str | None] | None = None,
) -> dict[str, Any]:
    """Validate ontology accessions against a local ontology file.

    Shared by the CLI and the Python API so both apply the same semantic checks
    rather than each assembling their own.
    """

    ontology = CellOntology.load(ontology_path)
    expected = list(expected_names) if expected_names is not None else [None] * len(terms)
    if len(expected) != len(terms):
        raise ContractError(
            "expected_names must have exactly one entry per term when provided"
        )
    return {
        "status": "complete",
        "source": ontology.source,
        "n_terms": len(ontology.terms),
        "terms": [
            ontology.validate_term(term, expected_name=name)
            for term, name in zip(terms, expected, strict=True)
        ],
    }


def provision_from_config(config: CellCuratorConfig) -> dict[str, Any]:
    """Report the configured ontology provider's provisioning status.

    A recorded status file always wins: provisioning is an immutable run
    artifact, so a later read must not re-derive a different answer.
    """

    ontology = ontology_from_config(config)
    root = Path(config.run.output_root) / config.run.run_id
    status_path = root / "knowledge" / "ontology_provider_status.json"
    if status_path.is_file():
        value = json.loads(status_path.read_text())
        return dict(value) if isinstance(value, dict) else {"status": "unknown"}
    return {
        "schema_version": 3,
        "status": "unavailable" if ontology is None else "complete",
        "reason": "no ontology provider is configured" if ontology is None else "",
        "remediation": "Configure a local ontology or remote_ontology provider.",
        "n_terms": len(ontology.terms) if ontology is not None else 0,
    }

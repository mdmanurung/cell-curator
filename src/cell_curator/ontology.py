"""Local Cell Ontology parsing and semantic hierarchy validation."""

from __future__ import annotations

import csv
import json
import re
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from .contracts import ContractError, safe_path_component, sha256
from .net import urlopen

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
    alt_ids: tuple[str, ...] = ()


def normalize_accession(value: str) -> str:
    """Return the ``PREFIX:LOCAL`` CURIE for a CURIE, underscore form, or OBO IRI."""

    text = str(value).strip().rsplit("/", 1)[-1].rsplit("#", 1)[-1]
    match = re.fullmatch(r"([A-Za-z]+)[:_](\w+)", text)
    return f"{match.group(1).upper()}:{match.group(2)}" if match else str(value).strip()


class CellOntology:
    def __init__(self, terms: dict[str, OntologyTerm], *, source: str) -> None:
        self.terms = terms
        self.source = source
        self._names: dict[str, set[str]] = {}
        self._aliases: dict[str, str] = {
            alias: accession for accession, term in terms.items() for alias in term.alt_ids
        }
        for accession, term in terms.items():
            for name in (term.name, *term.synonyms):
                self._names.setdefault(name.casefold(), set()).add(accession)
        self._roots: frozenset[str] = frozenset(
            accession
            for accession, term in terms.items()
            if not any(parent in terms for parent in term.parents)
        )
        self._depths: dict[str, int] = {}

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
                    alt_ids=tuple(map(str, current.get("alt_ids", []))),
                )
            current = None

        def identifier(value: str) -> str:
            # `is_a: CL:0000084 {is_inferred="true"} ! T cell`: the trailing
            # qualifier block and comment are not part of the identifier, and
            # keeping them disconnected 1,320 edges of the real Cell Ontology.
            return value.split(" ! ", 1)[0].split(" {", 1)[0].strip()

        for raw in path.read_text().splitlines():
            line = raw.strip()
            if line == "[Term]":
                commit()
                current = {"parents": [], "synonyms": [], "alt_ids": []}
                continue
            if line.startswith("["):
                commit()
                continue
            if current is None or ": " not in line:
                continue
            key, value = line.split(": ", 1)
            if key in {"id", "replaced_by"}:
                current[key] = identifier(value)
            elif key == "name":
                current[key] = value.strip()
            elif key == "is_a":
                current["parents"].append(identifier(value))
            elif key == "alt_id":
                current["alt_ids"].append(identifier(value))
            elif key == "synonym":
                # Only EXACT synonyms name the same term; BROAD, NARROW, and RELATED
                # ones name other terms and made common names ambiguous.
                match = re.match(r'"([^"]+)"\s+(EXACT|BROAD|NARROW|RELATED)\b', value)
                if match and match.group(2) == "EXACT":
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
                elif local == "hasExactSynonym" and child.text:
                    synonyms.append(child.text.strip())
                elif local == "deprecated" and child.text:
                    obsolete = child.text.strip().lower() == "true"
                elif local in {"replacedBy", "IAO_0100001"}:  # IAO_0100001: term replaced by
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
        accession = normalize_accession(value)
        accession = self._aliases.get(accession, accession)
        if accession in self.terms:
            return self.terms[accession]
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

    def roots(self) -> frozenset[str]:
        """Report the terms that have no parent inside the supplied graph.

        "Root" is a property of the snapshot, not of the Cell Ontology. A term
        is a root here exactly when none of its `is_a` parents were loaded, which
        is what an OBO, OWL, or table export of a subtree looks like: the top of
        the file is a cut, not a biological top. Consumers that reason about
        ancestry need this set because the whole graph shares its roots by
        construction, so a shared root is the one piece of ancestry that
        distinguishes nothing.
        """

        return self._roots

    def depth(self, accession: str) -> int:
        """Return the longest `is_a` chain from a term to a root of this graph.

        Depth is the number of edges, so a root of the supplied graph is 0. The
        longest chain rather than a particular one is what makes "deeper"
        well defined for a DAG: a term with two parents of different depth sits
        below the deeper one, and any statement about how far below an
        agreement point something lies has to use the same measure on both
        sides or it is not a distance.

        Results are memoized because an evaluation over a prediction table asks
        for the same handful of terms thousands of times. The walk is iterative
        rather than recursive so that a malformed snapshot with a parent cycle
        cannot exhaust the interpreter stack, and a cycle resolves to depth 0:
        a cycle is a broken hierarchy, and inventing a depth for it would let a
        corrupt file produce a plausible-looking number.
        """

        start = self.resolve(accession).accession
        if start in self._depths:
            return self._depths[start]
        pending: list[tuple[str, bool]] = [(start, False)]
        active: set[str] = set()
        while pending:
            current, expanded = pending.pop()
            if current in self._depths:
                continue
            parents = [p for p in self.terms[current].parents if p in self.terms]
            if not parents:
                self._depths[current] = 0
                continue
            if expanded:
                self._depths[current] = 1 + max(self._depths[p] for p in parents)
                continue
            if current in active:
                self._depths[current] = 0
                continue
            active.add(current)
            pending.append((current, True))
            pending.extend((p, False) for p in parents if p not in self._depths)
        return self._depths[start]

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
            policy=config.knowledge.policy,
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


_LIVE_STATUS_KIND = "cell_curator_live_ontology_status"


def live_ontology_provider(config: CellCuratorConfig) -> Any:
    """Return the single configured live ontology provider, or ``None``."""

    providers = [item for item in config.knowledge.providers if item.kind == "live_ontology"]
    if len(providers) > 1:
        raise ContractError("configure at most one live ontology provider")
    return providers[0] if providers else None


def _live_term_from_payload(payload: object, accession: str) -> OntologyTerm | None:
    """Read one OLS-shaped term document without assuming a single response shape."""

    if not isinstance(payload, dict):
        return None
    candidates: list[dict[str, Any]] = []
    embedded = payload.get("_embedded")
    if isinstance(embedded, dict) and isinstance(embedded.get("terms"), list):
        candidates = [item for item in embedded["terms"] if isinstance(item, dict)]
    elif isinstance(payload.get("terms"), list):
        candidates = [item for item in payload["terms"] if isinstance(item, dict)]
    elif payload.get("obo_id") or payload.get("short_form"):
        candidates = [payload]
    for item in candidates:
        obo = str(item.get("obo_id") or item.get("short_form") or "").replace("_", ":")
        if obo and obo.casefold() != accession.casefold():
            continue
        replaced = item.get("term_replaced_by")
        synonyms = item.get("synonyms") or item.get("obo_synonym") or []
        return OntologyTerm(
            accession=obo or accession,
            name=str(item.get("label") or ""),
            synonyms=tuple(
                str(value.get("name") if isinstance(value, dict) else value)
                for value in synonyms
            ),
            obsolete=bool(item.get("is_obsolete", False)),
            replaced_by=str(replaced) if replaced else None,
        )
    return None


def _live_findings(
    live: OntologyTerm,
    local: OntologyTerm | None,
    expected_name: str | None,
) -> list[str]:
    findings: list[str] = []
    if live.obsolete:
        replacement = f"; replaced by {live.replaced_by}" if live.replaced_by else ""
        findings.append(f"live service reports the term obsolete{replacement}")
    if local is not None:
        if local.obsolete != live.obsolete:
            findings.append(
                "local obsolescence disagrees with the live service "
                f"(local={local.obsolete}, live={live.obsolete})"
            )
        if live.name and local.name and live.name.casefold() != local.name.casefold():
            findings.append(
                f"local name {local.name!r} disagrees with live name {live.name!r}"
            )
    if expected_name:
        known = {live.name.casefold(), *(value.casefold() for value in live.synonyms)}
        if expected_name.casefold() not in known:
            findings.append(
                f"label {expected_name!r} is not the live name or a live synonym of "
                f"{live.accession}"
            )
    return findings


def corroborate_terms_live(
    config: CellCuratorConfig,
    accessions: Sequence[str],
    *,
    expected_names: Sequence[str | None] | None = None,
    ontology: CellOntology | None = None,
    run_root: Path | None = None,
) -> dict[str, Any]:
    """Corroborate local ontology terms against a live service, never depending on it.

    A live service has no stable digest, so reproducibility comes from the answer
    this run stored rather than from the endpoint: each response is written under
    ``knowledge/ontology_live/`` in canonical JSON and hashed, and the receipt
    records the endpoint, the query, that digest, and the retrieval time. A later
    reader reproduces the run from its own artifact.

    Nothing here may ground a label. ``require_local_grounding_for_labels`` is
    pinned true, so this reports agreement and disagreement and leaves the local
    answer standing — including when the network is absent, which is an expected
    state for an optional provider and not a failure.
    """

    from .provenance import replace_json

    names: list[str | None] = list(expected_names or [None] * len(accessions))
    if len(names) != len(accessions):
        raise ContractError("expected_names must be supplied once per accession")
    root = Path(run_root or Path(config.run.output_root) / config.run.run_id)
    status_path = root / "knowledge" / "ontology_live_status.json"

    def unavailable(reason: str, remediation: str, *, attempted: bool) -> dict[str, Any]:
        value = {
            "schema_version": 3,
            "artifact_kind": _LIVE_STATUS_KIND,
            "status": "unavailable",
            "reason": reason,
            "remediation": remediation,
            "network_attempted": attempted,
            "corroborated": [],
            "automated_output_role": "hypothesis_not_verdict",
        }
        replace_json(status_path, value)
        return value

    provider = live_ontology_provider(config)
    if provider is None:
        return unavailable(
            "no live ontology provider is configured",
            "Add a live_ontology provider to knowledge.providers to corroborate CL terms.",
            attempted=False,
        )
    if not config.knowledge.allow_optional_remote_enhancement:
        return unavailable(
            "optional remote enhancement is disabled",
            "Set knowledge.network to institution or online to corroborate CL terms; "
            "local validation and CL-free annotation remain available.",
            attempted=False,
        )

    rows: list[dict[str, Any]] = []
    for accession, expected_name in zip(accessions, names, strict=True):
        query = f"{provider.url}?{urllib.parse.urlencode({'obo_id': accession})}"
        local: OntologyTerm | None = None
        if ontology is not None:
            local = ontology.terms.get(accession)
        row: dict[str, Any] = {
            "accession": accession,
            "expected_name": expected_name or "",
            "url": query,
            "local_name": local.name if local else "",
            "local_obsolete": bool(local.obsolete) if local else None,
        }
        try:
            request = urllib.request.Request(
                query, headers={"User-Agent": "cell-curator/0.2", "Accept": "application/json"}
            )
            with urlopen(
                request, timeout=provider.timeout_seconds, policy=config.knowledge.policy
            ) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except (OSError, ValueError) as exc:
            rows.append({**row, "status": "unavailable", "findings": [str(exc)]})
            continue
        artifact = (
            root / "knowledge" / "ontology_live" / f"{safe_path_component(accession)}.json"
        )
        digest = replace_json(artifact, payload)
        row.update(
            {
                "artifact": str(artifact),
                "sha256": digest,
                "retrieved_at": datetime.now(UTC).isoformat(timespec="seconds"),
            }
        )
        live = _live_term_from_payload(payload, accession)
        if live is None:
            rows.append(
                {
                    **row,
                    "status": "absent",
                    "findings": [f"live service returned no term for {accession}"],
                }
            )
            continue
        findings = _live_findings(live, local, expected_name)
        rows.append(
            {
                **row,
                "status": "disagreement" if findings else "agreement",
                "live_name": live.name,
                "live_obsolete": live.obsolete,
                "live_replaced_by": live.replaced_by or "",
                "findings": findings,
            }
        )

    value = {
        "schema_version": 3,
        "artifact_kind": _LIVE_STATUS_KIND,
        "status": "complete",
        "reason": "",
        "remediation": "",
        "network_attempted": True,
        "endpoint": provider.url,
        "version": provider.version,
        "n_agreement": sum(1 for row in rows if row["status"] == "agreement"),
        "n_disagreement": sum(1 for row in rows if row["status"] == "disagreement"),
        "n_unavailable": sum(1 for row in rows if row["status"] in {"unavailable", "absent"}),
        "corroborated": rows,
        "automated_output_role": "hypothesis_not_verdict",
    }
    replace_json(status_path, value)
    return value


# Ancestors too broad to be a call: every PBMC label shares them, so backing off
# to one says nothing a reader did not know.
BACKOFF_FLOOR = frozenset(
    {
        "CL:0000000",  # cell
        "CL:0000548",  # animal cell
        "CL:0000255",  # eukaryotic cell
        "CL:0002371",  # somatic cell
        "CL:0000226",  # single nucleate cell
        "CL:0000081",  # blood cell
        "CL:0000988",  # hematopoietic cell
        "CL:0000738",  # leukocyte
        "CL:0000842",  # mononuclear cell
    }
)
BACKOFF_ASSUMPTION = (
    "When near-tied cell type candidates share one deepest Cell Ontology ancestor "
    "that is not a generic term (cell, leukocyte, mononuclear cell, and similar), "
    "and their combined vote share reaches the confidence threshold, the cluster "
    "is labelled with that ancestor instead of Unknown."
)


def backoff_term(
    ontology: CellOntology,
    candidates: list[tuple[str, float]],
    *,
    min_margin: float,
    min_confidence: float,
) -> tuple[OntologyTerm, float, list[str]] | None:
    """Deepest ancestor shared by the near-tied candidates, or None.

    ``candidates`` are (CL id, vote share) pairs, best first. The near-tie set is
    every candidate within ``min_margin`` of the best. The ancestor must be the
    unique deepest common one, outside :data:`BACKOFF_FLOOR`, and the candidates
    under it must together hold ``min_confidence`` of the votes.
    """

    if not candidates:
        return None
    top = candidates[0][1]
    tied = [(cl_id, share) for cl_id, share in candidates if share >= top - min_margin]
    if len(tied) < 2 or any(not cl_id for cl_id, _ in tied):
        return None
    try:
        closure = {
            cl_id: ontology.ancestors(cl_id) | {ontology.resolve(cl_id).accession}
            for cl_id, _ in candidates
            if cl_id
        }
    except ContractError:
        return None
    common = {
        item
        for item in set.intersection(*(closure[cl_id] for cl_id, _ in tied))
        # CL only: snapshots import UBERON and BFO, and a tie across lineages
        # would otherwise back off to "anatomical structure".
        if item.startswith("CL:") and item in ontology.terms and item not in BACKOFF_FLOOR
    }
    if not common:
        return None
    deepest = max(ontology.depth(item) for item in common)
    lowest = [item for item in common if ontology.depth(item) == deepest]
    if len(lowest) != 1:
        return None
    mass = sum(share for cl_id, share in candidates if cl_id and lowest[0] in closure[cl_id])
    if mass < min_confidence:
        return None
    return ontology.terms[lowest[0]], min(mass, 1.0), [cl_id for cl_id, _ in tied]

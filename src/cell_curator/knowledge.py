"""Independent local-first knowledge providers for marker and reference evidence."""

from __future__ import annotations

import hashlib
import json
import os
import re
import urllib.request
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Literal, Protocol, cast

import numpy as np
import pandas as pd
from scipy import sparse

from .config import CellCuratorConfig, ProviderConfig
from .contracts import ContractError, sha256
from .models import MarkerProgram
from .net import OFFLINE, NetworkPolicy, urlopen


def _features(value: object) -> tuple[str, ...]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ()
    if isinstance(value, (list, tuple)):
        return tuple(dict.fromkeys(str(item).strip() for item in value if str(item).strip()))
    return tuple(
        dict.fromkeys(
            token.strip() for token in re.split(r"[;,|]", str(value)) if token.strip()
        )
    )


class MarkerProvider(Protocol):
    """Provider contract; implementations return immutable typed programs."""

    @property
    def identity(self) -> str: ...

    def programs(self) -> list[MarkerProgram]: ...


class LocalMarkerProvider:
    def __init__(self, path: str | Path, *, version: str = "unversioned") -> None:
        self.path = Path(path)
        self.version = version

    @property
    def identity(self) -> str:
        return f"local-markers:{self.path.resolve()}:{self.version}"

    def programs(self) -> list[MarkerProgram]:
        if not self.path.is_file():
            raise ContractError(f"marker source does not exist: {self.path}")
        suffix = self.path.suffix.lower()
        if suffix == ".json":
            return self._json_programs()
        if suffix in {".tsv", ".csv"}:
            separator = "\t" if suffix == ".tsv" else ","
            return self._table_programs(
                pd.read_csv(self.path, sep=separator, keep_default_na=False)
            )
        raise ContractError("local marker source must be JSON, TSV, or CSV")

    def _json_programs(self) -> list[MarkerProgram]:
        value = json.loads(self.path.read_text())
        if isinstance(value, dict):
            rows = []
            for label, payload in value.items():
                if not isinstance(payload, dict):
                    raise ContractError(f"marker program {label!r} must be a mapping")
                rows.append({"label": label, **payload})
        elif isinstance(value, list):
            rows = value
        else:
            raise ContractError("marker JSON must be a mapping or list")
        return [self._program(row) for row in rows]

    def _table_programs(self, frame: pd.DataFrame) -> list[MarkerProgram]:
        if "label" not in frame:
            raise ContractError("marker table requires a label column")
        return [self._program(row) for row in frame.to_dict("records")]

    def _program(self, row: dict) -> MarkerProgram:
        axis_value = str(row.get("axis", "cell_type"))
        if axis_value not in {"cell_type", "cell_state"}:
            raise ContractError(f"marker program axis is invalid: {axis_value!r}")
        axis = cast(Literal["cell_type", "cell_state"], axis_value)
        return MarkerProgram(
            label=str(row["label"]),
            axis=axis,
            positive=_features(row.get("positive", row.get("pos", ()))),
            negative=_features(row.get("negative", row.get("neg", ()))),
            parent=str(row["parent"]).strip() if row.get("parent") else None,
            parent_cl_id=(
                str(row["parent_cl_id"]).strip() if row.get("parent_cl_id") else None
            ),
            cl_id=str(row["cl_id"]).strip() if row.get("cl_id") else None,
            feature_space=str(row.get("feature_space", "gene")),
            source=self.identity,
            source_version=self.version,
        )


class CorpusProvider:
    """Expose only context-matched, reviewed records from a pinned snapshot."""

    def __init__(self, config: ProviderConfig, *, context: dict[str, str] | None = None,
                 measurements: tuple[str, ...] = (), ontology: Any = None) -> None:
        from .corpus import load_corpus

        self.snapshot = load_corpus(config.path, expected_sha256=config.sha256)
        if self.snapshot.revision != config.version:
            raise ContractError("corpus revision differs from the configured version")
        self.context = merge_corpus_context(context or {}, config.context)
        self.measurements = measurements
        self.ontology = ontology

    @property
    def identity(self) -> str:
        return f"corpus:{self.snapshot.corpus_id}:{self.snapshot.revision}:{self.snapshot.digest}"

    def programs(self) -> list[MarkerProgram]:
        return self.snapshot.programs(context=self.context, ontology=self.ontology,
                                      measurements=self.measurements or None)


def corpus_context(config: CellCuratorConfig) -> dict[str, str]:
    context = config.guidance.context
    assay = config.input.assays.get(config.annotation.assay)
    if assay is None and len(config.input.assays) == 1:
        assay = next(iter(config.input.assays.values()))
    return {"species": context.organism, "tissue": context.tissue,
            "stage": context.developmental_stage, "condition": context.condition,
            "study_id": config.annotation.study_id,
            "assay": config.annotation.corpus_assay or (assay.kind if assay else context.technology)}


def merge_corpus_context(context: dict[str, str], overrides: dict[str, str]) -> dict[str, str]:
    """Provider scopes cannot replace or infer biological facts about a dataset."""
    aliases = {"human": "homo sapiens", "9606": "homo sapiens", "ncbitaxon:9606": "homo sapiens",
               "mouse": "mus musculus", "10090": "mus musculus", "ncbitaxon:10090": "mus musculus"}
    for key, value in overrides.items():
        known = context.get(key, "").strip().casefold()
        proposed = value.strip().casefold()
        if key == "species":
            known, proposed = aliases.get(known, known), aliases.get(proposed, proposed)
        if key in {"species", "tissue", "stage", "condition"} and known in {"", "unknown"} and proposed not in {"", "unknown"}:
            raise ContractError(f"corpus provider cannot infer unknown dataset {key}")
        if known and known != "unknown" and proposed != known:
            raise ContractError(f"corpus provider context contradicts dataset {key}")
    return {**context, **overrides}


def run_providers(config: CellCuratorConfig) -> list[MarkerProvider]:
    from .ontology import ontology_from_config

    measurements = {"rna": "gene_expression", "protein": "protein_expression",
                    "gene_activity": "gene_activity", "atac": "accessibility",
                    "spatial": "gene_expression"}
    assays = (config.input.assays[config.annotation.assay],) if config.annotation.enabled else tuple(config.input.assays.values())
    return providers_from_config(
        config.knowledge.providers, policy=config.knowledge.policy,
        context=corpus_context(config),
        measurements=tuple(sorted({measurements[assay.kind] for assay in assays
                                   if assay.kind in measurements})),
        ontology=ontology_from_config(config) if any(item.kind == "corpus" for item in config.knowledge.providers) else None,
    )


class CachedTableProvider(LocalMarkerProvider):
    """Versioned database export using the same signed table contract."""

    @property
    def identity(self) -> str:
        return f"cached-table:{self.path.resolve()}:{self.version}"


class RemoteTableProvider(CachedTableProvider):
    """Optional checksum-pinned remote table with a durable local cache."""

    def __init__(
        self,
        url: str,
        cache_path: str | Path,
        *,
        expected_sha256: str,
        version: str,
        timeout_seconds: float = 30.0,
        policy: NetworkPolicy = OFFLINE,
    ) -> None:
        super().__init__(cache_path, version=version)
        self.url = url
        self.expected_sha256 = expected_sha256
        self.timeout_seconds = timeout_seconds
        self.policy = policy

    @property
    def identity(self) -> str:
        return f"remote-table:{self.url}:{self.version}"

    def provision(self) -> Path:
        if self.path.is_file():
            self._validate_checksum()
            return self.path
        if self.policy.mode == "offline":
            raise ContractError(
                "optional remote enhancement is disabled and the checksum-pinned "
                f"cache is absent: {self.path}"
            )
        if not self.url or not self.expected_sha256:
            raise ContractError("remote providers require a URL and SHA-256")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        partial = self.path.with_name(f".{self.path.name}.partial.{os.getpid()}")
        request = urllib.request.Request(
            self.url,
            headers={"User-Agent": "cell-curator/0.2"},
        )
        try:
            with urlopen(request, timeout=self.timeout_seconds, policy=self.policy) as response:
                partial.write_bytes(response.read())
            observed = sha256(partial)
            if observed != self.expected_sha256:
                raise ContractError(
                    f"remote table checksum mismatch: expected {self.expected_sha256}, got {observed}"
                )
            os.replace(partial, self.path)
        finally:
            partial.unlink(missing_ok=True)
        return self.path

    def _validate_checksum(self) -> None:
        if self.expected_sha256 and sha256(self.path) != self.expected_sha256:
            raise ContractError(f"cached remote table checksum mismatch: {self.path}")

    def programs(self) -> list[MarkerProgram]:
        self.provision()
        return super().programs()


class ReferenceProvider:
    """Derive reproducible positive programs from a labeled local reference object."""

    def __init__(
        self,
        path: str | Path,
        *,
        label_key: str,
        version: str = "unversioned",
        modality: str = "",
        top_n: int = 30,
    ) -> None:
        self.path = Path(path)
        self.label_key = label_key
        self.version = version
        self.modality = modality
        self.top_n = top_n

    @property
    def identity(self) -> str:
        return f"reference:{self.path.resolve()}:{self.label_key}:{self.version}"

    def programs(self) -> list[MarkerProgram]:
        if not self.path.is_file():
            raise ContractError(f"reference does not exist: {self.path}")
        if self.path.suffix.lower() == ".h5mu":
            import mudata as md

            mdata = md.read_h5mu(self.path)
            if not self.modality or self.modality not in mdata.mod:
                raise ContractError("reference MuData provider requires a valid modality")
            adata = mdata.mod[self.modality]
            labels = (
                adata.obs[self.label_key]
                if self.label_key in adata.obs
                else mdata.obs[self.label_key]
            )
        else:
            import anndata as ad

            adata = ad.read_h5ad(self.path)
            if self.label_key not in adata.obs:
                raise ContractError(f"reference label key is absent: {self.label_key}")
            labels = adata.obs[self.label_key]
        matrix = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
        features = np.asarray(adata.var_names.astype(str))
        programs: list[MarkerProgram] = []
        for label in sorted(labels.dropna().astype(str).unique()):
            selected = labels.astype(str).eq(label).to_numpy()
            if not selected.any() or selected.all():
                continue
            effect = matrix[selected].mean(axis=0) - matrix[~selected].mean(axis=0)
            order = np.lexsort((features, -effect))
            positive = tuple(features[order[: self.top_n]])
            programs.append(
                MarkerProgram(
                    label=label,
                    positive=positive,
                    feature_space="gene",
                    source=self.identity,
                    source_version=self.version,
                )
            )
        return programs


def providers_from_config(
    configs: Iterable[ProviderConfig],
    *,
    policy: NetworkPolicy = OFFLINE,
    context: dict[str, str] | None = None,
    measurements: tuple[str, ...] = (),
    ontology: Any = None,
) -> list[MarkerProvider]:
    providers: list[MarkerProvider] = []
    for item in configs:
        if item.kind == "local_markers":
            providers.append(LocalMarkerProvider(item.path, version=item.version))
        elif item.kind == "corpus":
            providers.append(CorpusProvider(item, context=context, measurements=measurements,
                                            ontology=ontology))
        elif item.kind == "cached_table":
            providers.append(CachedTableProvider(item.path, version=item.version))
        elif item.kind == "remote_table":
            if policy.mode == "offline" and not Path(item.cache_path).is_file():
                # The run-level knowledge manifest records this configured provider
                # as unavailable.  It is deliberately omitted from execution so a
                # local vocabulary can proceed without hidden network I/O.
                continue
            providers.append(
                RemoteTableProvider(
                    item.url,
                    item.cache_path,
                    expected_sha256=item.sha256,
                    version=item.version,
                    timeout_seconds=item.timeout_seconds,
                    policy=policy,
                )
            )
        elif item.kind == "reference":
            providers.append(
                ReferenceProvider(
                    item.path,
                    label_key=item.label_key,
                    version=item.version,
                )
            )
    return providers


def collect_programs(providers: Iterable[MarkerProvider]) -> list[MarkerProgram]:
    """Merge equivalent programs while keeping positive/negative provenance additive."""

    merged: dict[tuple[str, str, str, str | None, str | None], MarkerProgram] = {}
    corpus_keys: set[tuple[str, str, str, str | None, str | None]] = set()
    for provider in providers:
        for program in provider.programs():
            key = (
                program.axis,
                program.label,
                program.feature_space,
                program.parent,
                program.parent_cl_id,
            )
            if key in merged and (isinstance(provider, CorpusProvider) or key in corpus_keys):
                raise ContractError("colliding corpus definitions require an explicit exact mapping in one snapshot")
            if isinstance(provider, CorpusProvider):
                corpus_keys.add(key)
            if key not in merged:
                merged[key] = program
                continue
            current = merged[key]
            if current.cl_id and program.cl_id and current.cl_id != program.cl_id:
                raise ContractError(f"conflicting CL IDs for marker program {program.label!r}")
            # Validated, not model_copy: two sources can disagree on a gene's sign.
            merged[key] = MarkerProgram.model_validate(
                {
                    **current.model_dump(),
                    "positive": tuple(dict.fromkeys((*current.positive, *program.positive))),
                    "negative": tuple(dict.fromkeys((*current.negative, *program.negative))),
                    "cl_id": current.cl_id or program.cl_id,
                    "source": f"{current.source};{program.source}",
                }
            )
    return sorted(merged.values(), key=lambda item: (item.axis, item.parent or "", item.label))


def programs_digest(programs: Iterable[MarkerProgram]) -> str:
    payload = "\n".join(
        program.model_dump_json(exclude_none=False)
        for program in sorted(
            programs, key=lambda item: (item.axis, item.parent or "", item.label)
        )
    )
    return hashlib.sha256(payload.encode()).hexdigest()

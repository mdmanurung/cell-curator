"""Versioned, strict configuration models."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import Field, ValidationError, field_validator, model_validator

from .models import SCHEMA_VERSION, RunMode, StrictModel


class DetectionConfig(StrictModel):
    source: Literal["X", "x", "raw", "layer"]
    key: str = ""
    required: bool = False
    semantics: str = "> 0"

    @model_validator(mode="after")
    def require_layer_key(self) -> DetectionConfig:
        if self.source.lower() == "layer" and not self.key.strip():
            raise ValueError("detection layer key must not be blank")
        return self


class RepresentationConfig(StrictModel):
    source: Literal["X", "x", "raw", "layer"]
    key: str = ""
    feature_space: str | None = None
    ranking_backend: str
    ontology_compatible: bool = False
    detection: DetectionConfig | None = None
    peak_to_gene_path: str = ""
    motif_mapping_path: str = ""
    weight: float = Field(default=1.0, gt=0)

    @model_validator(mode="after")
    def require_layer_key(self) -> RepresentationConfig:
        if self.source.lower() == "layer" and not self.key.strip():
            raise ValueError("representation layer key must not be blank")
        if self.peak_to_gene_path and self.motif_mapping_path:
            raise ValueError(
                "a representation must choose peak-to-gene or motif mapping, not both"
            )
        if not self.ranking_backend.strip():
            raise ValueError("ranking_backend must not be blank")
        return self


class GraphRepresentationConfig(StrictModel):
    source: Literal["obsm", "X", "x", "layer"]
    key: str = ""


class AssayConfig(StrictModel):
    modality: str = "root"
    kind: Literal["rna", "protein", "atac", "gene_activity", "spatial", "other"]
    feature_space: str
    feature_id_column: str = ""
    gene_mapping_path: str = ""
    capabilities: list[str] = Field(default_factory=list)
    graph_representations: list[GraphRepresentationConfig] = Field(default_factory=list)
    representations: dict[str, RepresentationConfig]


class InputKeysConfig(StrictModel):
    canonical_parent: str
    scope: str
    stored_resolutions: list[str] = Field(default_factory=list)
    embedding: str = ""
    spatial: str = ""
    donor: str = ""
    capture: str = ""
    doublet_class: str = ""
    doublet_score: str = ""
    sex: str = ""
    batch: str = ""


class InputConfig(StrictModel):
    kind: Literal["h5ad", "h5mu"]
    path: str
    expected_sha256: str = ""
    assays: dict[str, AssayConfig]
    keys: InputKeysConfig
    required_obs_keys: list[str] = Field(default_factory=list)
    guide_columns: list[str] = Field(default_factory=list)
    technical_metrics: dict[str, str] = Field(default_factory=dict)

    @model_validator(mode="after")
    def valid_source(self) -> InputConfig:
        if not self.path.strip():
            raise ValueError("input.path must not be blank")
        if self.expected_sha256 and not re.fullmatch(r"[0-9a-f]{64}", self.expected_sha256):
            raise ValueError("input.expected_sha256 must be a lowercase SHA-256 hex digest")
        if not self.assays:
            raise ValueError("input.assays must declare at least one assay")
        return self


class RunConfig(StrictModel):
    run_id: str
    output_root: str
    mode: RunMode = RunMode.EVIDENCE_ONLY
    immutable: bool = True

    @model_validator(mode="after")
    def safe_run_path(self) -> RunConfig:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", self.run_id):
            raise ValueError("run.run_id must be a safe 1-128 character path component")
        if not self.output_root.strip():
            raise ValueError("run.output_root must not be blank")
        if not self.immutable:
            raise ValueError("run.immutable must remain true")
        return self


class ProviderConfig(StrictModel):
    kind: Literal[
        "local_markers",
        "ontology",
        "remote_ontology",
        "reference",
        "cached_table",
        "remote_table",
    ]
    path: str = ""
    url: str = ""
    sha256: str = ""
    version: str = "unversioned"
    label_key: str = ""
    parent_key: str = ""
    cache_path: str = ""
    timeout_seconds: float = Field(default=30.0, gt=0)

    @model_validator(mode="after")
    def provider_contract(self) -> ProviderConfig:
        if (
            self.kind in {"local_markers", "ontology", "reference", "cached_table"}
            and not self.path
        ):
            raise ValueError(f"{self.kind} provider requires path")
        if self.kind == "reference" and not self.label_key:
            raise ValueError("reference provider requires label_key")
        if self.kind in {"remote_table", "remote_ontology"}:
            missing = [
                name
                for name, value in {
                    "url": self.url,
                    "cache_path": self.cache_path,
                    "sha256": self.sha256,
                    "version": self.version,
                }.items()
                if not value or value == "unversioned"
            ]
            if missing:
                raise ValueError(f"{self.kind} provider misses: {', '.join(missing)}")
            if not re.fullmatch(r"[0-9a-f]{64}", self.sha256):
                raise ValueError(f"{self.kind} sha256 must be a lowercase SHA-256 digest")
        return self


class KnowledgeConfig(StrictModel):
    providers: list[ProviderConfig] = Field(default_factory=list)
    require_local_grounding_for_labels: bool = True
    allow_optional_remote_enhancement: bool = True

    @field_validator("require_local_grounding_for_labels")
    @classmethod
    def preserve_local_grounding(cls, value: bool) -> bool:
        if not value:
            raise ValueError("knowledge.require_local_grounding_for_labels must remain true")
        return value


class CanonicalPartitionConfig(StrictModel):
    parent_key: str
    expected_parent_count: int = Field(ge=1)


class EligibilityConfig(StrictModel):
    min_cells: int = Field(default=50, ge=2)
    min_donors: int | None = Field(default=None, ge=1)
    min_captures: int | None = Field(default=None, ge=1)
    rare_exceptions_require_approval: bool = True

    @field_validator("rare_exceptions_require_approval")
    @classmethod
    def preserve_rare_exception_gate(cls, value: bool) -> bool:
        if not value:
            raise ValueError(
                "adaptive.eligibility.rare_exceptions_require_approval must remain true"
            )
        return value


class AuditConfig(StrictModel):
    guide_entropy_min: float = Field(default=0.6, ge=0, le=1)
    meaningful_minority_fraction: float = Field(default=0.1, ge=0, le=1)
    doublet_fraction_max: float | None = Field(default=None, ge=0, le=1)
    donor_max_fraction: float | None = Field(default=None, ge=0, le=1)
    capture_max_fraction: float | None = Field(default=None, ge=0, le=1)
    technical_effect_size_max: float | None = Field(default=0.15, ge=0)
    technical_metric_robust_z_max: float = Field(default=3.0, ge=0)
    incompatible_program_fraction_min: float | None = Field(default=0.15, ge=0, le=1)
    incompatible_program_columns: list[str] = Field(default_factory=list)


class StoredCandidateConfig(StrictModel):
    min_parent_purity: float = Field(default=0.9, ge=0, le=1)
    min_parent_fraction: float = Field(default=0.03, gt=0, lt=1)
    max_parent_fraction: float = Field(default=0.95, gt=0, le=1)
    adjacent_jaccard_min: float = Field(default=0.8, ge=0, le=1)
    min_adjacent_matches: int = Field(default=1, ge=1)


class ParentLocalConfig(StrictModel):
    enabled: bool = True
    max_passes_per_parent: int = Field(default=1, ge=0, le=1)
    resolutions: list[float] = Field(default_factory=lambda: [0.2, 0.4, 0.8])
    stability_metric: Literal["ari", "jaccard"] = "ari"
    stability_threshold_strictly_greater_than: float = Field(default=0.8, ge=0, le=1)
    reference_seeds: list[int] = Field(default_factory=lambda: [0, 1, 2, 3, 4])
    reject_one_cluster: bool = True
    select: Literal["lowest_stable_biologically_resolving"] = (
        "lowest_stable_biologically_resolving"
    )
    knn_k: int = Field(default=15, ge=2)
    bootstrap_fraction: float = Field(default=0.8, gt=0.5, le=1)

    @field_validator("reject_one_cluster")
    @classmethod
    def preserve_one_cluster_rejection(cls, value: bool) -> bool:
        if not value:
            raise ValueError("adaptive.parent_local.reject_one_cluster must remain true")
        return value


class AdaptiveConfig(StrictModel):
    outcomes: list[str]
    eligibility: EligibilityConfig
    audit: AuditConfig
    stored_candidates: StoredCandidateConfig
    parent_local: ParentLocalConfig


class ComparatorConfig(StrictModel):
    roles: list[Literal["parent_remainder", "closest_sibling"]]
    closest_sibling_strategy: Literal["embedding_centroid"] = "embedding_centroid"


class GpuContractConfig(StrictModel):
    name_contains: str
    compute_capability: str


class MarkerConfig(StrictModel):
    rapids_singlecell_version: str = ""
    rapids_distribution: str = ""
    cuda_major: int | None = None
    cpu_fallback: bool = False
    allowed_gpu_contracts: list[GpuContractConfig] = Field(default_factory=list)
    wilcoxon: dict[str, Any] = Field(default_factory=dict)
    logreg: dict[str, Any] = Field(default_factory=dict)

    @field_validator("cpu_fallback")
    @classmethod
    def prohibit_gpu_fallback(cls, value: bool) -> bool:
        if value:
            raise ValueError("CPU fallback is prohibited for RAPIDS lanes")
        return value


class EvidenceLaneConfig(StrictModel):
    assay: str
    representation: str
    lane_key: str | None = None
    required_for_acceptance: bool = True


class EvidenceConfig(StrictModel):
    lanes: list[EvidenceLaneConfig] = Field(default_factory=list)
    adjusted_pvalue_max: float = Field(default=0.05, gt=0, le=1)
    require_both_marker_methods: bool = True
    require_both_comparators: bool = True
    required_candidate_gates: list[str] = Field(default_factory=list)
    parent_vocabulary: dict[str, str] = Field(default_factory=dict)
    min_positive_coverage: float = Field(default=0.35, ge=0, le=1)
    max_negative_violation: float = Field(default=0.25, ge=0, le=1)
    min_confidence: float = Field(default=0.55, ge=0, le=1)
    min_margin: float = Field(default=0.10, ge=0)
    max_normalized_entropy: float = Field(default=0.85, ge=0, le=1)

    @model_validator(mode="after")
    def preserve_independent_crosschecks(self) -> EvidenceConfig:
        if not self.require_both_marker_methods:
            raise ValueError("evidence.require_both_marker_methods must remain true")
        if not self.require_both_comparators:
            raise ValueError("evidence.require_both_comparators must remain true")
        return self


class ReviewConfig(StrictModel):
    critic_passes: list[str] = Field(default_factory=list)
    require_reconciliation: bool = True
    require_explicit_human_approval: bool = True
    type_and_state_separate: bool = True

    @model_validator(mode="after")
    def preserve_scientific_safeguards(self) -> ReviewConfig:
        if not self.require_reconciliation:
            raise ValueError("review.require_reconciliation must remain true")
        if not self.require_explicit_human_approval:
            raise ValueError("review.require_explicit_human_approval must remain true")
        if not self.type_and_state_separate:
            raise ValueError("review.type_and_state_separate must remain true")
        return self


class HpcConfig(StrictModel):
    cpu_resources: dict[str, Any] = Field(default_factory=dict)
    gpu_resources: dict[str, Any] = Field(default_factory=dict)


class BiologicalContextConfig(StrictModel):
    """Biological and reporting facts that cannot be inferred safely from the object."""

    organism: str = ""
    tissue: str = ""
    developmental_stage: str = ""
    condition: str = ""
    technology: Literal[
        "dissociated", "whole_transcriptome_spatial", "targeted_panel", "other"
    ] = "other"
    experimental_context: str = ""
    composition_covariates: list[str] = Field(default_factory=list)
    vocabulary_contract: str = ""
    reference_contract: str = "local_only"
    visualization_key: str = ""
    writeback_prefix: str = "cell_curator"

    def unresolved(self) -> list[str]:
        required = {
            "organism": self.organism,
            "tissue": self.tissue,
            "developmental_stage": self.developmental_stage,
            "condition": self.condition,
            "experimental_context": self.experimental_context,
            "vocabulary_contract": self.vocabulary_contract,
            "reference_contract": self.reference_contract,
            "writeback_prefix": self.writeback_prefix,
        }
        unresolved = [name for name, value in required.items() if not value.strip()]
        if self.technology == "other":
            unresolved.append("technology")
        return unresolved


class GuidanceConfig(StrictModel):
    """Interactive/non-interactive gate contract for the three-stage workflow."""

    interactive: bool = True
    preauthorized: bool = False
    reviewer: str = ""
    pause_after_each_level: bool = True
    assumptions: list[str] = Field(default_factory=list)
    hierarchy_levels: list[str] = Field(default_factory=lambda: ["L1"])
    context: BiologicalContextConfig = Field(default_factory=BiologicalContextConfig)

    @model_validator(mode="after")
    def validate_preauthorization(self) -> GuidanceConfig:
        if not self.hierarchy_levels or any(not item.strip() for item in self.hierarchy_levels):
            raise ValueError("guidance.hierarchy_levels must contain non-blank level names")
        if len(set(self.hierarchy_levels)) != len(self.hierarchy_levels):
            raise ValueError("guidance.hierarchy_levels must be unique and ordered")
        if self.preauthorized and self.interactive:
            raise ValueError("preauthorized runs must set guidance.interactive=false")
        if not self.interactive:
            missing = self.context.unresolved()
            if not self.preauthorized:
                missing.append("guidance.preauthorized=true")
            if not self.reviewer.strip():
                missing.append("guidance.reviewer")
            if not self.assumptions:
                missing.append("guidance.assumptions")
            if missing:
                raise ValueError(
                    "non-interactive operation requires complete recorded context and "
                    "pre-authorization; missing: " + ", ".join(missing)
                )
        return self


class CellCuratorConfig(StrictModel):
    schema_version: int
    run: RunConfig
    knowledge: KnowledgeConfig = Field(default_factory=KnowledgeConfig)
    input: InputConfig
    canonical_partitions: dict[str, CanonicalPartitionConfig]
    adaptive: AdaptiveConfig
    comparators: ComparatorConfig
    markers: MarkerConfig
    evidence: EvidenceConfig
    review: ReviewConfig
    hpc: HpcConfig = Field(default_factory=HpcConfig)
    guidance: GuidanceConfig = Field(default_factory=GuidanceConfig)

    @field_validator("schema_version")
    @classmethod
    def current_schema(cls, value: int) -> int:
        if value != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}")
        return value

    @model_validator(mode="after")
    def references_exist(self) -> CellCuratorConfig:
        if (
            self.run.mode is RunMode.ANNOTATION
            and self.knowledge.require_local_grounding_for_labels
        ):
            grounding_kinds = {"local_markers", "reference", "cached_table"}
            if not any(item.kind in grounding_kinds for item in self.knowledge.providers):
                raise ValueError(
                    "annotation mode requires a local marker, labeled reference, or cached-table "
                    "provider for biological grounding"
                )
        gpu_lane = False
        lane_keys: set[str] = set()
        for lane in self.evidence.lanes:
            if lane.assay not in self.input.assays:
                raise ValueError(f"evidence lane refers to unknown assay {lane.assay!r}")
            assay = self.input.assays[lane.assay]
            if lane.representation not in assay.representations:
                raise ValueError(
                    f"assay {lane.assay!r} lacks representation {lane.representation!r}"
                )
            backend = assay.representations[lane.representation].ranking_backend.lower()
            gpu_lane = gpu_lane or backend in {"rapids_singlecell", "rapids-singlecell"}
            key = lane.lane_key or f"{lane.assay}__{lane.representation}"
            if key in lane_keys:
                raise ValueError(f"evidence lane key is duplicated: {key}")
            lane_keys.add(key)
        if gpu_lane:
            missing = []
            if not self.markers.rapids_singlecell_version:
                missing.append("markers.rapids_singlecell_version")
            if not self.markers.rapids_distribution:
                missing.append("markers.rapids_distribution")
            if self.markers.cuda_major is None:
                missing.append("markers.cuda_major")
            if not self.markers.allowed_gpu_contracts:
                missing.append("markers.allowed_gpu_contracts")
            if missing:
                raise ValueError(
                    "RAPIDS lanes require a complete fail-closed GPU contract: "
                    + ", ".join(missing)
                )
        return self


def load_config(path: str | Path) -> CellCuratorConfig:
    source = Path(path)
    try:
        value = yaml.safe_load(source.read_text())
    except (OSError, yaml.YAMLError) as exc:
        raise ValueError(f"could not read configuration {source}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError("configuration must be a YAML mapping")
    try:
        return CellCuratorConfig.model_validate(value)
    except ValidationError as exc:
        raise ValueError(f"invalid cell-curator configuration: {exc}") from exc


def configuration_schema() -> dict[str, Any]:
    """Return the canonical machine-readable schema for configuration version 3."""

    schema = CellCuratorConfig.model_json_schema()
    schema["$id"] = "urn:cell-curator:schema:configuration:3"
    schema["x-cell-curator-schema-version"] = SCHEMA_VERSION
    return schema

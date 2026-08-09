"""Typed public models shared by the annotation engine and artifact readers."""

from __future__ import annotations

from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

SCHEMA_VERSION: Literal[3] = 3


class StrictModel(BaseModel):
    """Base model that rejects misspelled or undeclared fields."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class CapabilityLevel(StrEnum):
    ANNOTATION = "annotation-capable"
    CORROBORATION = "corroboration-only"
    QC = "qc-only"
    UNAVAILABLE = "unavailable"


class AssayKind(StrEnum):
    RNA = "rna"
    PROTEIN = "protein"
    ATAC = "atac"
    GENE_ACTIVITY = "gene_activity"
    SPATIAL = "spatial"
    OTHER = "other"


class RunMode(StrEnum):
    ANNOTATION = "annotation"
    EVIDENCE_ONLY = "evidence-only"


class ParentOutcome(StrEnum):
    ACCEPT_SPLIT = "ACCEPT SPLIT"
    RETAIN_PARENT = "RETAIN PARENT"
    DOUBLET_MIXED = "DOUBLET/MIXED"
    TECHNICAL_UNRESOLVED = "TECHNICAL/UNRESOLVED"


class MarkerProgram(StrictModel):
    schema_version: Literal[3] = SCHEMA_VERSION
    label: str
    axis: Literal["cell_type", "cell_state"] = "cell_type"
    positive: tuple[str, ...] = ()
    negative: tuple[str, ...] = ()
    parent: str | None = None
    parent_cl_id: str | None = None
    cl_id: str | None = None
    feature_space: str = "gene"
    source: str = "local"
    source_version: str = "unversioned"

    @field_validator("label")
    @classmethod
    def nonblank_label(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("marker program label must not be blank")
        return value.strip()


class Capability(StrictModel):
    schema_version: Literal[3] = SCHEMA_VERSION
    assay: str
    representation: str
    assay_kind: AssayKind
    feature_space: str
    level: CapabilityLevel
    claims: tuple[str, ...]
    reasons: tuple[str, ...] = ()
    backend: str
    requires_gpu: bool = False


class EvidenceRecord(StrictModel):
    schema_version: Literal[3] = SCHEMA_VERSION
    run_id: str
    scope: str
    cluster_id: str
    label: str
    axis: Literal["cell_type", "cell_state"] = "cell_type"
    assay: str
    representation: str
    capability_level: CapabilityLevel
    score: float
    positive_coverage: float
    negative_violation: float
    measured_positive: int
    measured_negative: int
    missing_positive: tuple[str, ...] = ()
    missing_negative: tuple[str, ...] = ()
    detection_valid: bool
    source: str


class AnnotationDecision(StrictModel):
    schema_version: Literal[3] = SCHEMA_VERSION
    run_id: str
    scope: str
    cluster_id: str
    parent_label: str | None = None
    cell_type: str
    cell_state: str = "normal"
    cl_id: str | None = None
    confidence: float = Field(ge=0.0, le=1.0)
    margin: float
    normalized_entropy: float = Field(ge=0.0, le=1.0)
    abstained: bool
    rationale: str
    rejected_alternatives: tuple[str, ...] = ()
    required_capabilities_complete: bool
    approval_status: Literal["PENDING", "APPROVED", "REJECTED"] = "PENDING"


class MappingRecord(StrictModel):
    schema_version: Literal[3] = SCHEMA_VERSION
    run_id: str
    barcode: str
    scope: str
    canonical_parent: str
    leaf_id: str
    cell_type: str
    cell_state: str
    membership_role: Literal["retained_parent", "accepted_child", "parent_remainder"]


class BackendReceipt(StrictModel):
    schema_version: Literal[3] = SCHEMA_VERSION
    artifact_kind: str = "cell_curator_backend_receipt"
    run_id: str
    lane_key: str
    backend: str
    capability_level: CapabilityLevel
    status: Literal["complete", "unavailable", "failed"]
    input_sha256: str
    config_sha256: str
    output_sha256: str | None = None
    artifact_sha256: dict[str, str] = Field(default_factory=dict)
    feature_order_sha256: str | None = None
    group_membership_sha256: str | None = None
    fallback_used: bool = False
    runtime: dict[str, Any] = Field(default_factory=dict)
    limitations: tuple[str, ...] = ()
    created_utc: str = Field(default_factory=lambda: datetime.now(UTC).isoformat())


class RunState(StrictModel):
    schema_version: Literal[3] = SCHEMA_VERSION
    run_id: str
    mode: RunMode
    phase: str
    input_path: str
    input_sha256: str
    config_sha256: str
    canonical_membership_sha256: str
    blockers: tuple[str, ...] = ()
    completed_phases: tuple[str, ...] = ()
    assumptions_approved: bool = False
    critics_reconciled: bool = False
    human_approval: bool = False
    labels_permitted: bool = False
    writeback_permitted: bool = False
    stage: Literal["scene", "plan", "execute", "review", "writeback", "complete"] = "scene"
    current_level: str = ""
    current_parent: str = ""
    next_action: str = "complete scene setting"
    compute_profile: str = ""
    compute: dict[str, Any] = Field(default_factory=dict)
    context_complete: bool = False
    plan_complete: bool = False
    interactive: bool = True
    preauthorized: bool = False
    assumptions_sha256: str = ""
    approved_assumptions_sha256: str = ""
    signed_off_scopes: tuple[str, ...] = ()
    decisions: tuple[dict[str, Any], ...] = ()
    progress: dict[str, Any] = Field(default_factory=dict)
    updated_utc: str = Field(default_factory=lambda: datetime.now(UTC).isoformat())


class BenchmarkResult(StrictModel):
    schema_version: Literal[3] = SCHEMA_VERSION
    system: str
    dataset: str
    n: int
    accuracy: float
    macro_f1: float
    micro_f1: float
    coverage: float
    selective_accuracy: float
    coverage_risk: float
    hierarchical_accuracy: float | None = None
    unknown_precision: float | None = None
    unknown_recall: float | None = None
    unknown_f1: float | None = None
    brier_score: float | None = None
    expected_calibration_error: float | None = None
    hierarchy_error: float | None = None
    ontology_distance_error: float | None = None
    mixed_detection_f1: float | None = None
    split_decision_accuracy: float | None = None
    repeat_reproducibility: float | None = None
    input_invariance: float | None = None
    unsplit_parent_invariance: float | None = None
    runtime_seconds: float | None = None
    peak_memory_mb: float | None = None
    artifact_completeness: float | None = None
    artifact_path: Path | None = None


class BenchmarkArtifact(StrictModel):
    """One checksum-pinned acquisition artifact in a benchmark recipe."""

    name: str
    url: str
    sha256: str

    @model_validator(mode="after")
    def validate_artifact(self) -> BenchmarkArtifact:
        if not self.name.strip():
            raise ValueError("benchmark artifact name must not be blank")
        if not self.url.startswith("https://"):
            raise ValueError("benchmark artifact URLs must use HTTPS")
        if self.sha256 != "pending-after-download" and not _is_sha256(self.sha256):
            raise ValueError(
                "benchmark artifact checksum must be SHA-256 or explicitly pending"
            )
        return self


class BenchmarkDatasetManifest(StrictModel):
    """Versioned, neutral public or local benchmark acquisition contract."""

    schema_version: Literal[3] = SCHEMA_VERSION
    dataset: str
    modality: str
    source_url: str
    artifacts: tuple[BenchmarkArtifact, ...] = ()
    sha256: str
    license: str
    status: Literal["available", "pending-external-resource"]
    blocker: str | None = None
    reproduce: tuple[str, ...] = ()
    local_path: str | None = None
    scenarios: tuple[str, ...] = ()

    @model_validator(mode="after")
    def validate_status_contract(self) -> BenchmarkDatasetManifest:
        for field_name in ("dataset", "modality", "source_url", "license"):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(f"benchmark {field_name} must not be blank")
        if not (
            self.source_url.startswith("https://") or self.source_url.startswith("local://")
        ):
            raise ValueError("benchmark source_url must use HTTPS or local://")
        if self.status == "available":
            if not self.local_path:
                raise ValueError("available benchmarks require local_path")
            if not _is_sha256(self.sha256):
                raise ValueError("available benchmarks require an exact SHA-256 checksum")
        else:
            if not self.blocker or not self.blocker.strip():
                raise ValueError("pending benchmarks require an explicit blocker")
            if not self.reproduce:
                raise ValueError("pending benchmarks require reproduction commands")
            if self.sha256 != "pending-after-artifact-selection":
                raise ValueError("pending benchmark checksum must be explicitly pending")
        if any(not command.strip() for command in self.reproduce):
            raise ValueError("benchmark reproduction commands must not be blank")
        return self


def _is_sha256(value: str) -> bool:
    return len(value) == 64 and all(
        character in "0123456789abcdef" for character in value.lower()
    )

"""Independent, evidence-driven multimodal cell annotation."""

from importlib.metadata import PackageNotFoundError, version

from .api import (
    annotate,
    annotation_progress,
    apply_writeback,
    approve_annotation_assumptions,
    build_review_artifacts,
    canonicalize_input,
    check_review_report,
    compute_input_qc,
    detect_compute_environment,
    execute_recursive_hierarchy,
    execute_recursive_hierarchy_manifest,
    generate_evidence,
    initialize_run,
    inspect_input,
    plan_annotation,
    prepare_annotation_scene,
    preview_writeback,
    propose_annotations,
    register_assay_representation,
    register_knowledge_provider,
    route_annotation_crosscheck,
    run_marker_tools,
    validate_run,
)
from .config import CellCuratorConfig, load_config
from .models import (
    AnnotationDecision,
    BackendReceipt,
    BenchmarkArtifact,
    BenchmarkDatasetManifest,
    BenchmarkResult,
    Capability,
    CapabilityLevel,
    EvidenceRecord,
    MappingRecord,
    MarkerProgram,
    RunState,
)
from .session import CellCurator

try:
    __version__ = version("cell-curator")
except PackageNotFoundError:  # Source checkout before installation.
    __version__ = "0+unknown"

__all__ = [
    "AnnotationDecision",
    "BackendReceipt",
    "BenchmarkArtifact",
    "BenchmarkDatasetManifest",
    "BenchmarkResult",
    "Capability",
    "CapabilityLevel",
    "CellCurator",
    "CellCuratorConfig",
    "EvidenceRecord",
    "MappingRecord",
    "MarkerProgram",
    "RunState",
    "annotation_progress",
    "annotate",
    "apply_writeback",
    "approve_annotation_assumptions",
    "build_review_artifacts",
    "canonicalize_input",
    "check_review_report",
    "compute_input_qc",
    "detect_compute_environment",
    "execute_recursive_hierarchy",
    "execute_recursive_hierarchy_manifest",
    "generate_evidence",
    "initialize_run",
    "inspect_input",
    "load_config",
    "plan_annotation",
    "prepare_annotation_scene",
    "preview_writeback",
    "propose_annotations",
    "register_assay_representation",
    "register_knowledge_provider",
    "route_annotation_crosscheck",
    "run_marker_tools",
    "validate_run",
]

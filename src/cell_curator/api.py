"""Stable Python API for independent cell-curator workflows.

The functions import their implementation lazily so applications can inspect the
package version and configuration models without importing the scientific stack.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .config import CellCuratorConfig, ProviderConfig, RepresentationConfig
from .contracts import ContractError


def inspect_input(config_path: str | Path) -> dict[str, Any]:
    """Inspect an AnnData or MuData source without changing it."""

    from .config import load_config
    from .data import inspect

    return inspect(load_config(config_path))


def canonicalize_input(
    config_path: str | Path, run_directory: str | Path | None = None
) -> dict[str, Any]:
    """Create the hash-frozen working copy and canonical membership artifacts."""

    from .config import load_config
    from .data import canonicalize
    from .pipeline import run_root

    config = load_config(config_path)
    return canonicalize(
        config,
        config_path=Path(config_path),
        run_root=Path(run_directory) if run_directory is not None else run_root(config),
    )


def compute_input_qc(
    config_path: str | Path, run_directory: str | Path | None = None
) -> dict[str, Any]:
    """Compute or validate source-preserving QC on the canonical working copy."""

    from .config import load_config
    from .data import canonicalize, compute_qc
    from .pipeline import run_root

    config = load_config(config_path)
    root = Path(run_directory) if run_directory is not None else run_root(config)
    if not (root / "qc" / "qc.manifest.json").is_file():
        return canonicalize(config, config_path=Path(config_path), run_root=root)["qc"]
    return compute_qc(config, run_root=root)


def detect_compute_environment(config_path: str | Path | None = None) -> Any:
    """Detect CPU, RAM, Apple Silicon, CUDA, and scheduler state without mutation."""

    from .environment import detect_environment

    cpu: dict[str, Any] = {}
    gpu: dict[str, Any] = {}
    if config_path is not None:
        from .config import load_config

        config = load_config(config_path)
        cpu, gpu = config.hpc.cpu_resources, config.hpc.gpu_resources
    return detect_environment(
        configured_cpu_resources=cpu,
        configured_gpu_resources=gpu,
    )


def route_annotation_crosscheck(request: Any) -> Any:
    """Select one deterministic automated cross-check for a hierarchy scope."""

    from .routing import RoutingRequest, route_crosscheck

    model = (
        request
        if isinstance(request, RoutingRequest)
        else RoutingRequest.model_validate(request)
    )
    return route_crosscheck(model)


def prepare_annotation_scene(config_path: str | Path) -> dict[str, Any]:
    """Persist Stage 1 inspection, context, compute, assumptions, and progress."""

    from .guided import prepare_scene

    return prepare_scene(config_path)


def plan_annotation(config_path: str | Path) -> dict[str, Any]:
    """Persist Stage 2 vocabulary, evidence routes, and unresolved assumptions."""

    from .guided import construct_plan

    return construct_plan(config_path)


def approve_annotation_assumptions(
    config_path: str | Path,
    *,
    reviewer: str,
    scope: str = "*",
    preauthorized: bool = False,
) -> dict[str, Any]:
    """Record interactive or declared non-interactive assumption approval."""

    from .guided import approve_assumptions

    return approve_assumptions(
        config_path,
        reviewer=reviewer,
        scope=scope,
        preauthorized=preauthorized,
    )


def annotation_progress(config_path: str | Path) -> dict[str, Any]:
    """Return current stage, scope, blockers, progress, and next action."""

    from .guided import status

    return status(config_path)


def run_marker_tools(
    config_path: str | Path,
    *,
    assay: str,
    representation: str,
    group_key: str,
    parent_scope: str | None = None,
    output_directory: str | Path | None = None,
) -> Any:
    """Run the complete signed marker profile and cross-check surface."""

    from .markers import run_marker_analysis_from_config

    return run_marker_analysis_from_config(
        config_path,
        assay=assay,
        representation=representation,
        group_key=group_key,
        parent_scope=parent_scope,
        output_dir=output_directory,
    )


def execute_recursive_hierarchy(
    config_path: str | Path,
    *,
    score_tables: Mapping[tuple[str, str], Any],
    vocabulary_path: str | Path | None = None,
    state_tables: Mapping[tuple[str, str], Any] | None = None,
    minimum_confidence: float | None = None,
    minimum_margin: float | None = None,
) -> dict[str, Any]:
    """Execute and resume arbitrary-depth parent-scoped annotation."""

    from .pipeline import execute_recursive_hierarchy as execute

    return execute(
        config_path,
        score_tables=score_tables,
        vocabulary_path=vocabulary_path,
        state_tables=state_tables,
        minimum_confidence=minimum_confidence,
        minimum_margin=minimum_margin,
    )


def execute_recursive_hierarchy_manifest(
    config_path: str | Path,
    manifest_path: str | Path,
) -> dict[str, Any]:
    """Execute recursive annotation from a strict JSON artifact manifest."""

    from .pipeline import execute_recursive_hierarchy_manifest as execute

    return execute(config_path, manifest_path)


def check_review_report(config_path: str | Path, *, level: str | None = None) -> dict[str, Any]:
    """Run strict read-only report validation after approval."""

    from .config import load_config
    from .pipeline import run_root
    from .review import check_report_completeness

    config = load_config(config_path)
    return check_report_completeness(config, run_root(config), level=level)


def initialize_run(config_path: str | Path) -> dict[str, Any]:
    """Initialize or resume an immutable annotation run."""

    from .pipeline import initialize

    return initialize(Path(config_path))


def register_assay_representation(
    config: CellCuratorConfig,
    *,
    assay: str,
    name: str,
    representation: RepresentationConfig | dict[str, Any],
) -> CellCuratorConfig:
    """Return a validated config copy with one newly registered assay representation."""

    if assay not in config.input.assays:
        raise ContractError(f"cannot register a representation for unknown assay: {assay}")
    if not name.strip():
        raise ContractError("representation name must not be blank")
    declaration = (
        representation
        if isinstance(representation, RepresentationConfig)
        else RepresentationConfig.model_validate(representation)
    )
    input_value = config.input.model_dump(mode="python")
    representations = input_value["assays"][assay]["representations"]
    if name in representations:
        raise ContractError(f"assay representation is already registered: {assay}/{name}")
    representations[name] = declaration.model_dump(mode="python")
    return config.model_copy(
        update={"input": config.input.__class__.model_validate(input_value)}
    )


def register_knowledge_provider(
    config: CellCuratorConfig,
    provider: ProviderConfig | dict[str, Any],
) -> CellCuratorConfig:
    """Return a validated config copy with one local or optional cached provider."""

    declaration = (
        provider
        if isinstance(provider, ProviderConfig)
        else ProviderConfig.model_validate(provider)
    )
    providers = list(config.knowledge.providers)
    identity = (declaration.kind, declaration.path, declaration.version, declaration.url)
    if identity in {(item.kind, item.path, item.version, item.url) for item in providers}:
        raise ContractError("knowledge provider is already registered")
    providers.append(declaration)
    knowledge = config.knowledge.model_copy(update={"providers": providers})
    return config.model_copy(update={"knowledge": knowledge})


def generate_evidence(config_path: str | Path) -> tuple[Any, Any]:
    """Run and consolidate all configured, executable assay evidence lanes."""

    from .pipeline import evidence_all

    return evidence_all(Path(config_path))


def propose_annotations(config_path: str | Path) -> tuple[Any, Any]:
    """Derive child and parent annotation decisions from immutable evidence."""

    from .pipeline import decide

    return decide(Path(config_path))


def validate_run(config_path: str | Path) -> dict[str, Any]:
    """Validate completeness, hashes, mapping conservation, and review state."""

    from .pipeline import validate

    return validate(Path(config_path))


def build_review_artifacts(config_path: str | Path, *, strict: bool = True) -> Path:
    """Build evidence cards, critic packets, and the self-contained HTML report."""

    from .pipeline import review

    return review(Path(config_path), strict=strict)


def preview_writeback(
    config_path: str | Path, ontology_path: str | Path | None = None
) -> dict[str, Any]:
    """Return a source-preserving write-back diff after computational checks."""

    from .config import load_config
    from .ontology import CellOntology, ontology_from_config
    from .pipeline import run_root
    from .writeback import validate_writeback

    config = load_config(config_path)
    ontology = (
        CellOntology.load(Path(ontology_path))
        if ontology_path is not None
        else ontology_from_config(config)
    )
    return validate_writeback(config, run_root(config), ontology=ontology)


def apply_writeback(
    config_path: str | Path,
    *,
    output_path: str | Path,
    approval_path: str | Path,
    ontology_path: str | Path | None = None,
    in_place: bool = False,
) -> Path:
    """Apply an explicitly approved, atomic write-back to a new object by default."""

    from .config import load_config
    from .ontology import CellOntology, ontology_from_config
    from .pipeline import run_root
    from .writeback import write_back

    config = load_config(config_path)
    ontology = (
        CellOntology.load(Path(ontology_path))
        if ontology_path is not None
        else ontology_from_config(config)
    )
    return write_back(
        config,
        run_root(config),
        output_path=Path(output_path),
        approval_path=Path(approval_path),
        ontology=ontology,
        in_place=in_place,
    )


def annotate(config_path: str | Path) -> Path:
    """Execute every computational phase and freeze a review packet."""

    from .pipeline import run_annotation

    return run_annotation(Path(config_path))

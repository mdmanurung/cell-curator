"""Stable Python API for independent cell-curator workflows.

The functions import their implementation lazily so applications can inspect the
package version and configuration models without importing the scientific stack.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
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


def detect_compute_environment(
    config_path: str | Path | None = None, *, validate_cuda: bool = True
) -> Any:
    """Detect CPU, RAM, Apple Silicon, CUDA, and scheduler state without mutation.

    Set ``validate_cuda=False`` to skip the CUDA probe, matching the CLI's
    ``--no-cuda-validation``.
    """

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
        validate_cuda=validate_cuda,
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
    config_path: str | Path,
    ontology_path: str | Path | None = None,
    *,
    output_path: str | Path | None = None,
    in_place: bool = False,
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
    return validate_writeback(
        config,
        run_root(config),
        ontology=ontology,
        output_path=Path(output_path) if output_path is not None else None,
        in_place=in_place,
    )


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


# --- configuration, capabilities, and environment -----------------------------


def configuration_schema() -> dict[str, Any]:
    """Return the machine-readable schema for configuration version 3."""

    from .config import configuration_schema as schema

    return schema()


def list_capabilities(config_path: str | Path) -> list[Any]:
    """Return the assay-aware capability registry a configuration resolves to."""

    from .capabilities import build_capability_registry
    from .config import load_config

    return build_capability_registry(load_config(config_path))


# --- knowledge and ontology ---------------------------------------------------


def validate_knowledge_providers(config_path: str | Path) -> dict[str, Any]:
    """Validate configured local knowledge providers and any configured ontology."""

    from .config import load_config
    from .knowledge import collect_programs, providers_from_config
    from .ontology import ontology_from_config, validate_program_ontology

    config = load_config(config_path)
    programs = collect_programs(
        providers_from_config(
            config.knowledge.providers,
            allow_remote_enhancement=config.knowledge.allow_optional_remote_enhancement,
        )
    )
    ontology = ontology_from_config(config)
    if ontology is not None:
        validate_program_ontology(programs, ontology)
    return {
        "status": "complete",
        "external_service_required": False,
        "n_programs": len(programs),
        "labels": [item.label for item in programs],
        "ontology_status": "validated" if ontology is not None else "unavailable",
        "ontology_terms": len(ontology.terms) if ontology is not None else 0,
    }


def validate_ontology_terms(
    ontology_path: str | Path,
    terms: Sequence[str],
    *,
    expected_names: Sequence[str | None] | None = None,
) -> dict[str, Any]:
    """Validate ontology accessions against a local ontology file."""

    from .ontology import validate_terms

    return validate_terms(ontology_path, terms, expected_names=expected_names)


def provision_ontology(config_path: str | Path) -> dict[str, Any]:
    """Report the configured ontology provider's provisioning status."""

    from .config import load_config
    from .ontology import provision_from_config

    return provision_from_config(load_config(config_path))


# --- guided workflow ----------------------------------------------------------


def confirm_annotation_context(
    config_path: str | Path,
    *,
    input_path: str | Path,
    reviewer: str,
) -> dict[str, Any]:
    """Confirm or correct the scene context recorded in Stage 1."""

    from .guided import confirm_context

    return confirm_context(config_path, input_path=Path(input_path), reviewer=reviewer)


def list_annotation_assumptions(config_path: str | Path) -> dict[str, Any]:
    """Return the living assumptions ledger for a configured run."""

    from .guided import list_assumptions

    return list_assumptions(config_path)


def record_annotation_assumption(
    config_path: str | Path,
    *,
    scope: str,
    statement: str,
    material: bool = True,
) -> dict[str, Any]:
    """Append one assumption to the living ledger."""

    from .guided import record_assumption

    return record_assumption(
        config_path, scope=scope, statement=statement, material=material
    )


# --- pipeline phases ----------------------------------------------------------


def audit_parent_clusters(config_path: str | Path) -> Any:
    """Flag frozen parent clusters whose composition warrants a closer look."""

    from .pipeline import audit

    return audit(Path(config_path))


def discover_subcluster_candidates(config_path: str | Path) -> tuple[Any, Any, Any]:
    """Discover bounded subcluster candidates and build the comparator manifest."""

    from .pipeline import refine

    return refine(Path(config_path))


def run_evidence_lane(
    config_path: str | Path, *, assay: str, representation: str
) -> Path:
    """Execute one immutable assay/representation evidence lane."""

    from .pipeline import evidence_lane

    return evidence_lane(Path(config_path), assay=assay, representation=representation)


def run_candidate_evidence_lane(config_path: str | Path, *, lane_id: str) -> Path:
    """Execute one manifest-selected candidate comparison lane."""

    from .pipeline import candidate_evidence_lane

    return candidate_evidence_lane(Path(config_path), comparison_lane_id=lane_id)


def map_annotated_cells(config_path: str | Path) -> Any:
    """Derive the durable cell-to-label mapping from frozen decisions."""

    from .pipeline import map_cells

    return map_cells(Path(config_path))


def freeze_review_packet(config_path: str | Path) -> Path:
    """Freeze the immutable review packet directory for human sign-off."""

    from .pipeline import build_review_packet

    return build_review_packet(Path(config_path))


def assign_hierarchy_level(
    config_path: str | Path,
    *,
    scores_path: str | Path,
    level: str,
    parent: str,
    vocabulary_path: str | Path | None = None,
    state_path: str | Path | None = None,
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    """Assign one parent-scoped hierarchy level and checkpoint the labels table."""

    from .pipeline import assign_hierarchy_level as assign

    return assign(
        config_path,
        scores_path=scores_path,
        level=level,
        parent=parent,
        vocabulary_path=vocabulary_path,
        state_path=state_path,
        output_path=output_path,
    )


# --- markers and clusters -----------------------------------------------------


def build_marker_vocabulary(
    config_path: str | Path,
    *,
    source: Iterable[str | Path] = (),
    parent_scope: str | None = None,
    parent_cl_id: str | None = None,
    include_unscoped: bool = False,
    axis: str | None = None,
    feature_space: str | None = None,
    output_path: str | Path | None = None,
) -> tuple[list[Any], Path]:
    """Assemble and publish the signed marker vocabulary for a configuration."""

    from .markers import assemble_marker_programs_from_config

    return assemble_marker_programs_from_config(
        config_path,
        source=source,
        parent_scope=parent_scope,
        parent_cl_id=parent_cl_id,
        include_unscoped=include_unscoped,
        axis=axis,
        feature_space=feature_space,
        output_path=output_path,
    )


def profile_refined_clusters(config_path: str | Path) -> dict[str, Any]:
    """Summarize refined-child QC, composition, and spatial context."""

    from .hierarchy import profile_refined_clusters as profile

    return profile(config_path)


def advisory_confidence(
    *, margin: float, specificity: float, agreement: float, coverage: float
) -> float:
    """Combine four evidence dimensions into an explicitly uncalibrated score."""

    from .hierarchy import advisory_confidence as confidence

    return confidence(
        margin=margin,
        specificity=specificity,
        agreement=agreement,
        coverage=coverage,
    )


def validate_hierarchy_nesting(
    level_assignments: Iterable[Any], *, cell_membership: Any | None = None
) -> dict[str, Any]:
    """Enforce one parent per label and optional one-cell-one-leaf membership."""

    from .hierarchy import validate_hierarchy_nesting as validate

    return validate(level_assignments, cell_membership=cell_membership)


# --- automated cross-checks ---------------------------------------------------


def existing_prediction(
    *,
    input_path: str | Path,
    prediction_column: str,
    cluster_column: str,
    modality: str = "",
    output_dir: str | Path | None = None,
) -> tuple[Any, Any, dict[str, Any]]:
    """Turn an existing per-cell prediction column into one auditable hypothesis."""

    from .crosschecks import existing_prediction_file

    return existing_prediction_file(
        input_path=input_path,
        prediction_column=prediction_column,
        cluster_column=cluster_column,
        modality=modality,
        output_dir=output_dir,
    )


def discover_celltypist_models(*, output_dir: str | Path | None = None) -> dict[str, Any]:
    """List locally cached CellTypist models without downloading anything."""

    from .crosschecks import discover_celltypist_models as discover

    return discover(output_dir=output_dir)


def predict_celltypist(
    *,
    input_path: str | Path,
    model: str,
    cluster_column: str = "",
    modality: str = "",
    output_dir: str | Path | None = None,
) -> tuple[Any, dict[str, Any]]:
    """Run a locally cached CellTypist model as one hypothesis."""

    from .crosschecks import predict_celltypist_file

    return predict_celltypist_file(
        input_path=input_path,
        model=model,
        cluster_column=cluster_column,
        modality=modality,
        output_dir=output_dir,
    )


def map_reference(
    *,
    query_path: str | Path,
    reference_path: str | Path,
    reference_label_column: str,
    query_cluster_column: str = "",
    query_modality: str = "",
    reference_modality: str = "",
    output_dir: str | Path | None = None,
    **kwargs: Any,
) -> tuple[Any, Any, dict[str, Any]]:
    """Transfer labels from a labeled reference by nearest neighbours."""

    from .crosschecks import map_reference_files

    return map_reference_files(
        query_path=query_path,
        reference_path=reference_path,
        reference_label_column=reference_label_column,
        query_cluster_column=query_cluster_column,
        query_modality=query_modality,
        reference_modality=reference_modality,
        output_dir=output_dir,
        **kwargs,
    )


def map_curated_human_reference(
    *,
    query_path: str | Path,
    reference_path: str | Path,
    reference_label_column: str,
    query_cluster_column: str = "",
    query_modality: str = "",
    reference_modality: str = "",
    output_dir: str | Path | None = None,
    **kwargs: Any,
) -> tuple[Any, Any, dict[str, Any]]:
    """Map a query onto a curated human reference (Azimuth-style)."""

    from .crosschecks import map_curated_human_reference_files

    return map_curated_human_reference_files(
        query_path=query_path,
        reference_path=reference_path,
        reference_label_column=reference_label_column,
        query_cluster_column=query_cluster_column,
        query_modality=query_modality,
        reference_modality=reference_modality,
        output_dir=output_dir,
        **kwargs,
    )


def map_scimilarity(
    *,
    input_path: str | Path,
    model_path: str | Path,
    modality: str = "",
    output_dir: str | Path | None = None,
) -> tuple[Any, dict[str, Any]]:
    """Run a local SCimilarity model as one hypothesis."""

    from .crosschecks import map_scimilarity_file

    return map_scimilarity_file(
        input_path=input_path,
        model_path=model_path,
        modality=modality,
        output_dir=output_dir,
    )


def fetch_census_reference(
    *,
    organism: str,
    tissue: str,
    output_path: str | Path,
    census_version: str,
    max_cells: int = 50_000,
    label_columns: Iterable[str] = ("cell_type",),
    expected_sha256: str = "",
    timeout_seconds: float = 300.0,
) -> dict[str, Any]:
    """Acquire a bounded CELLxGENE Census slice with a versioned cache receipt."""

    from .crosschecks import fetch_census_reference as fetch

    return fetch(
        organism=organism,
        tissue=tissue,
        output_path=Path(output_path),
        census_version=census_version,
        max_cells=max_cells,
        label_columns=label_columns,
        expected_sha256=expected_sha256,
        timeout_seconds=timeout_seconds,
    )


def train_scanvi_reference_model(
    *,
    reference_path: str | Path,
    labels_key: str,
    output_dir: str | Path,
    modality: str = "",
    counts_layer: str = "counts",
    max_epochs: int = 200,
    seed: int = 0,
) -> dict[str, Any]:
    """Train a scANVI reference model from a labeled local reference object."""

    from .crosschecks import train_scanvi_reference_file

    return train_scanvi_reference_file(
        reference_path=reference_path,
        labels_key=labels_key,
        output_dir=output_dir,
        modality=modality,
        counts_layer=counts_layer,
        max_epochs=max_epochs,
        seed=seed,
    )


def map_scarches(
    *,
    query_path: str | Path,
    model_dir: str | Path,
    output_dir: str | Path,
    modality: str = "",
    max_epochs: int = 100,
) -> tuple[Any, dict[str, Any]]:
    """Map a query onto a trained reference model with scArches surgery."""

    from .crosschecks import map_scarches_file

    return map_scarches_file(
        query_path=query_path,
        model_dir=model_dir,
        output_dir=output_dir,
        modality=modality,
        max_epochs=max_epochs,
    )


# --- review, critics, and reporting ------------------------------------------


def build_evidence_cards(config_path: str | Path) -> list[Path]:
    """Scaffold one evidence card per proposed label."""

    from .config import load_config
    from .pipeline import run_root
    from .review import build_evidence_cards as build

    config = load_config(config_path)
    return build(config, run_root(config))


def build_html_report(
    config_path: str | Path, *, strict: bool = False, level: str | None = None
) -> Path:
    """Build the self-contained HTML review report."""

    from .config import load_config
    from .pipeline import run_root
    from .review import build_html_report as build

    config = load_config(config_path)
    return build(config, run_root(config), strict=strict, level=level)


def build_critic_packet(config_path: str | Path) -> Path:
    """Build the platform-neutral critic packet for adversarial review."""

    from .config import load_config
    from .pipeline import run_root
    from .review import build_critic_packet as build

    config = load_config(config_path)
    return build(config, run_root(config))


def import_critic_findings(config_path: str | Path, input_path: str | Path) -> dict[str, Any]:
    """Validate and atomically import structured critic findings."""

    from .config import load_config
    from .pipeline import run_root
    from .review import import_critic_findings as import_findings

    config = load_config(config_path)
    return import_findings(run_root(config), Path(input_path))


def import_critic_reconciliation(
    config_path: str | Path, input_path: str | Path
) -> dict[str, Any]:
    """Validate and atomically import reviewer dispositions for critic findings."""

    from .config import load_config
    from .pipeline import run_root
    from .review import import_reconciliation

    config = load_config(config_path)
    return import_reconciliation(run_root(config), Path(input_path))


def validate_critic_reconciliation(config_path: str | Path) -> dict[str, Any]:
    """Enforce that every critic finding carries a reviewer disposition."""

    from .config import load_config
    from .pipeline import run_root
    from .review import validate_reconciliation

    config = load_config(config_path)
    return validate_reconciliation(run_root(config))


# --- benchmarks ---------------------------------------------------------------


def evaluate_benchmark_predictions(
    prediction_path: str | Path,
    *,
    system: str,
    dataset: str,
    output_path: str | Path,
) -> Any:
    """Score a neutral-format prediction table against its recorded truth."""

    from .benchmark import evaluate_file

    return evaluate_file(
        Path(prediction_path),
        system=system,
        dataset=dataset,
        output_path=Path(output_path),
    )


def validate_dataset_manifest(path: str | Path) -> dict[str, Any]:
    """Validate one public benchmark dataset manifest."""

    from .benchmark import validate_dataset_manifest as validate

    return validate(Path(path))

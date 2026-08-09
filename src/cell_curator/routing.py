"""Deterministic, evidence-preserving automated cross-check routing.

The router only chooses among resources the caller has already discovered and
validated.  It performs no download, installation, model execution, or scheduler
submission.  Every known method receives an assessment so unavailable and
biologically unsuitable routes remain visible in the audit trail.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from .environment import ComputeProfile, EnvironmentReport
from .models import StrictModel


class HealthContext(StrEnum):
    HEALTHY = "healthy"
    DISEASE = "disease"
    TUMOR = "tumor"
    UNKNOWN = "unknown"


class Technology(StrEnum):
    DISSOCIATED = "dissociated"
    WHOLE_TRANSCRIPTOME_SPATIAL = "whole_transcriptome_spatial"
    TARGETED_PANEL = "targeted_panel"


class CrosscheckMethod(StrEnum):
    EXISTING_PREDICTION = "existing_prediction"
    REFERENCE_MODEL_SCARCHES = "reference_model_scarches"
    REFERENCE_KNN_GPU = "reference_knn_gpu"
    REFERENCE_KNN_CPU = "reference_knn_cpu"
    CELLTYPIST = "celltypist"
    AZIMUTH_COMPATIBLE = "azimuth_compatible"
    SCIMILARITY = "scimilarity"
    CENSUS_REFERENCE_KNN = "census_reference_knn"


ALL_METHODS: tuple[CrosscheckMethod, ...] = tuple(CrosscheckMethod)


class MethodResource(StrictModel):
    """A discovered prediction column, reference, model, or bounded provider."""

    name: str
    available: bool = True
    version: str = "unversioned"
    organisms: tuple[str, ...] = ()
    tissues: tuple[str, ...] = ()
    developmental_stages: tuple[str, ...] = ()
    health_contexts: tuple[HealthContext, ...] = ()
    technologies: tuple[Technology, ...] = ()
    levels: tuple[str, ...] = ()
    parent_scopes: tuple[str, ...] = ()
    gene_overlap_count: int | None = Field(default=None, ge=0)
    minimum_gene_overlap: int = Field(default=50, ge=1)
    estimated_memory_bytes: int | None = Field(default=None, gt=0)
    requires_gpu: bool = False
    network_required: bool = False
    checksum_validated: bool = False
    remediation: str = ""

    @model_validator(mode="after")
    def valid_identity(self) -> MethodResource:
        if not self.name.strip():
            raise ValueError("method resource name must not be blank")
        return self


class RoutingRequest(StrictModel):
    """Complete context needed to choose one cross-check for one parent level."""

    level: str
    parent_scope: str = "root"
    organism: str
    tissue: str = ""
    developmental_stage: str = "unknown"
    health_context: HealthContext = HealthContext.UNKNOWN
    technology: Technology = Technology.DISSOCIATED
    compute_profile: ComputeProfile
    usable_ram_bytes: int = Field(gt=0)
    gpu_validated: bool = False
    scheduler_ambiguity: bool = False
    allow_network: bool = False
    prefer_gpu_knn: bool = True
    existing_predictions: tuple[MethodResource, ...] = ()
    reference_models: tuple[MethodResource, ...] = ()
    labeled_references: tuple[MethodResource, ...] = ()
    celltypist_models: tuple[MethodResource, ...] = ()
    azimuth_references: tuple[MethodResource, ...] = ()
    scimilarity_models: tuple[MethodResource, ...] = ()
    census_references: tuple[MethodResource, ...] = ()

    @model_validator(mode="after")
    def complete_context(self) -> RoutingRequest:
        if not self.level.strip():
            raise ValueError("routing level must not be blank")
        if not self.parent_scope.strip():
            raise ValueError("routing parent_scope must not be blank")
        if not self.organism.strip():
            raise ValueError("routing organism must not be blank")
        if self.compute_profile is ComputeProfile.LOCAL_GPU and not self.gpu_validated:
            raise ValueError(
                "local_gpu routing requires validated real-GPU evidence; CPU substitution is disabled"
            )
        return self

    @classmethod
    def from_environment(
        cls,
        *,
        environment: EnvironmentReport,
        level: str,
        organism: str,
        parent_scope: str = "root",
        **values: object,
    ) -> RoutingRequest:
        """Create a request without losing environment ambiguity or GPU proof."""

        payload = dict(values)
        payload.update(
            level=level,
            parent_scope=parent_scope,
            organism=organism,
            compute_profile=environment.compute_profile,
            usable_ram_bytes=environment.usable_ram_bytes,
            gpu_validated=environment.gpu_validated,
            scheduler_ambiguity=environment.requires_confirmation,
        )
        return cls.model_validate(payload)


class MethodAssessment(StrictModel):
    method: CrosscheckMethod
    executable: bool
    suitable: bool
    selected: bool = False
    resource: str | None = None
    resource_version: str | None = None
    reasons: tuple[str, ...]
    remediation: str | None = None
    resources_considered: tuple[str, ...] = ()


class RoutingDecision(StrictModel):
    level: str
    parent_scope: str
    status: Literal["selected", "unavailable", "blocked"]
    primary_method: CrosscheckMethod | None
    primary_resource: str | None
    assessments: tuple[MethodAssessment, ...]
    blockers: tuple[str, ...] = ()
    next_action: str
    deterministic: Literal[True] = True
    detection_only: Literal[True] = True


_RESOURCE_FIELDS: dict[CrosscheckMethod, str] = {
    CrosscheckMethod.EXISTING_PREDICTION: "existing_predictions",
    CrosscheckMethod.REFERENCE_MODEL_SCARCHES: "reference_models",
    CrosscheckMethod.REFERENCE_KNN_GPU: "labeled_references",
    CrosscheckMethod.REFERENCE_KNN_CPU: "labeled_references",
    CrosscheckMethod.CELLTYPIST: "celltypist_models",
    CrosscheckMethod.AZIMUTH_COMPATIBLE: "azimuth_references",
    CrosscheckMethod.SCIMILARITY: "scimilarity_models",
    CrosscheckMethod.CENSUS_REFERENCE_KNN: "census_references",
}

_REMEDIATION: dict[CrosscheckMethod, str] = {
    CrosscheckMethod.EXISTING_PREDICTION: (
        "Provide a level-compatible prediction column already present in the input."
    ),
    CrosscheckMethod.REFERENCE_MODEL_SCARCHES: (
        "Provide a context-matched scANVI reference model and run on a validated real CUDA GPU."
    ),
    CrosscheckMethod.REFERENCE_KNN_GPU: (
        "Provide a labeled local reference with harmonized genes and validate a real CUDA GPU."
    ),
    CrosscheckMethod.REFERENCE_KNN_CPU: (
        "Provide a labeled local reference with harmonized genes and declared label metadata."
    ),
    CrosscheckMethod.CELLTYPIST: (
        "Install the optional CellTypist extra and make a context-matched model available."
    ),
    CrosscheckMethod.AZIMUTH_COMPATIBLE: (
        "Configure an optional curated human reference mapper with a versioned compatible reference."
    ),
    CrosscheckMethod.SCIMILARITY: (
        "Provision the optional checksum-verified SCimilarity-compatible human model bundle."
    ),
    CrosscheckMethod.CENSUS_REFERENCE_KNN: (
        "Authorize network access and configure a bounded, versioned CELLxGENE Census subset request."
    ),
}


def _normal(value: str) -> str:
    normalized = value.strip().casefold().replace("_", " ").replace("-", " ")
    aliases = {
        "homo sapiens": "human",
        "mus musculus": "mouse",
        "adult human": "adult",
        "adult mouse": "adult",
    }
    return aliases.get(normalized, normalized)


def _matches(value: str, allowed: tuple[str, ...]) -> bool:
    if not allowed:
        return True
    normalized = {_normal(item) for item in allowed}
    return "any" in normalized or _normal(value) in normalized


def _specificity(resource: MethodResource, request: RoutingRequest) -> int:
    return sum(
        bool(value)
        for value in (
            resource.organisms,
            resource.tissues,
            resource.developmental_stages,
            resource.health_contexts,
            resource.technologies,
            resource.levels,
            resource.parent_scopes,
        )
    )


def _resource_reasons(
    method: CrosscheckMethod,
    resource: MethodResource,
    request: RoutingRequest,
) -> tuple[bool, bool, list[str]]:
    """Return executable, suitable, and reasons for one declared resource."""

    executable = resource.available
    suitable = True
    reasons: list[str] = []
    if not resource.available:
        reasons.append(f"resource {resource.name!r} is unavailable")
        return False, False, reasons
    if not _matches(request.organism, resource.organisms):
        suitable = False
        reasons.append(
            f"resource {resource.name!r} does not support organism {request.organism!r}"
        )
    if request.tissue and not _matches(request.tissue, resource.tissues):
        suitable = False
        reasons.append(f"resource {resource.name!r} does not match tissue {request.tissue!r}")
    if resource.developmental_stages and not _matches(
        request.developmental_stage, resource.developmental_stages
    ):
        suitable = False
        reasons.append(
            f"resource {resource.name!r} does not match developmental stage "
            f"{request.developmental_stage!r}"
        )
    if resource.health_contexts and request.health_context not in resource.health_contexts:
        suitable = False
        reasons.append(
            f"resource {resource.name!r} does not match {request.health_context.value} context"
        )
    if resource.technologies and request.technology not in resource.technologies:
        suitable = False
        reasons.append(
            f"resource {resource.name!r} is unsuitable for {request.technology.value} technology"
        )
    if resource.levels and request.level not in resource.levels:
        suitable = False
        reasons.append(
            f"resource {resource.name!r} is not calibrated for level {request.level!r}"
        )
    if resource.parent_scopes and request.parent_scope not in resource.parent_scopes:
        suitable = False
        reasons.append(
            f"resource {resource.name!r} is outside parent scope {request.parent_scope!r}"
        )
    if (
        resource.estimated_memory_bytes
        and resource.estimated_memory_bytes > request.usable_ram_bytes
    ):
        executable = False
        reasons.append(
            f"resource {resource.name!r} needs {resource.estimated_memory_bytes} bytes but only "
            f"{request.usable_ram_bytes} bytes are usable"
        )
    if resource.network_required and not request.allow_network:
        executable = False
        reasons.append(
            f"resource {resource.name!r} requires network access that was not authorized"
        )

    gpu_method = method in {
        CrosscheckMethod.REFERENCE_MODEL_SCARCHES,
        CrosscheckMethod.REFERENCE_KNN_GPU,
    }
    if (gpu_method or resource.requires_gpu) and not request.gpu_validated:
        executable = False
        reasons.append(
            "a validated real CUDA GPU is required; CPU substitution is disabled for this method"
        )
    if gpu_method and request.compute_profile not in {
        ComputeProfile.LOCAL_GPU,
        ComputeProfile.SCHEDULER,
    }:
        executable = False
        reasons.append(f"compute profile {request.compute_profile.value!r} is not GPU-capable")

    if method is CrosscheckMethod.AZIMUTH_COMPATIBLE and _normal(request.organism) != "human":
        suitable = False
        reasons.append("the configured Azimuth-compatible route is curated for human data only")
    if method is CrosscheckMethod.SCIMILARITY:
        if _normal(request.organism) != "human":
            suitable = False
            reasons.append("SCimilarity-compatible models are human-only")
        if request.technology is Technology.TARGETED_PANEL:
            suitable = False
            reasons.append("SCimilarity is unsuitable for targeted-panel expression")

    needs_overlap = method not in {
        CrosscheckMethod.EXISTING_PREDICTION,
        CrosscheckMethod.SCIMILARITY,
    }
    if request.technology is Technology.TARGETED_PANEL and needs_overlap:
        if resource.gene_overlap_count is None:
            suitable = False
            reasons.append(
                f"targeted-panel overlap is undeclared for resource {resource.name!r}"
            )
        elif resource.gene_overlap_count < resource.minimum_gene_overlap:
            suitable = False
            reasons.append(
                f"targeted-panel overlap {resource.gene_overlap_count} is below the declared "
                f"minimum {resource.minimum_gene_overlap}"
            )

    if method is CrosscheckMethod.CENSUS_REFERENCE_KNN and not resource.checksum_validated:
        executable = False
        reasons.append("the bounded Census subset lacks a validated version/checksum receipt")
    if not reasons:
        reasons.append(f"resource {resource.name!r} matches the declared biological context")
    return executable, suitable, reasons


def _assess(method: CrosscheckMethod, request: RoutingRequest) -> MethodAssessment:
    resources = tuple(getattr(request, _RESOURCE_FIELDS[method]))
    if not resources:
        return MethodAssessment(
            method=method,
            executable=False,
            suitable=False,
            reasons=("no resource was declared for this method",),
            remediation=_REMEDIATION[method],
        )
    evaluated = [
        (resource, *_resource_reasons(method, resource, request)) for resource in resources
    ]
    suitable = [item for item in evaluated if item[1] is True and item[2] is True]
    considered = tuple(resource.name for resource in resources)
    if suitable:
        resource, executable, is_suitable, reasons = sorted(
            suitable,
            key=lambda item: (-_specificity(item[0], request), item[0].name.casefold()),
        )[0]
        preserved_reasons = list(reasons)
        for other, _executable, _is_suitable, other_reasons in evaluated:
            if other is resource:
                continue
            preserved_reasons.extend(
                f"resource {other.name!r} skipped: {reason}" for reason in other_reasons
            )
        return MethodAssessment(
            method=method,
            executable=executable,
            suitable=is_suitable,
            resource=resource.name,
            resource_version=resource.version,
            reasons=tuple(preserved_reasons),
            resources_considered=considered,
        )
    combined = tuple(
        reason for resource, _executable, _suitable, reasons in evaluated for reason in reasons
    )
    remediation = next(
        (resource.remediation for resource in resources if resource.remediation.strip()),
        _REMEDIATION[method],
    )
    return MethodAssessment(
        method=method,
        executable=any(item[1] for item in evaluated),
        suitable=False,
        reasons=combined or ("all declared resources were unsuitable",),
        remediation=remediation,
        resources_considered=considered,
    )


def _priority(request: RoutingRequest) -> tuple[CrosscheckMethod, ...]:
    gpu, cpu = (
        (CrosscheckMethod.REFERENCE_KNN_GPU, CrosscheckMethod.REFERENCE_KNN_CPU)
        if request.prefer_gpu_knn
        else (CrosscheckMethod.REFERENCE_KNN_CPU, CrosscheckMethod.REFERENCE_KNN_GPU)
    )
    if request.technology is Technology.WHOLE_TRANSCRIPTOME_SPATIAL:
        return (
            CrosscheckMethod.EXISTING_PREDICTION,
            CrosscheckMethod.REFERENCE_MODEL_SCARCHES,
            CrosscheckMethod.AZIMUTH_COMPATIBLE,
            gpu,
            cpu,
            CrosscheckMethod.CELLTYPIST,
            CrosscheckMethod.SCIMILARITY,
            CrosscheckMethod.CENSUS_REFERENCE_KNN,
        )
    if request.technology is Technology.TARGETED_PANEL:
        return (
            CrosscheckMethod.EXISTING_PREDICTION,
            gpu,
            cpu,
            CrosscheckMethod.AZIMUTH_COMPATIBLE,
            CrosscheckMethod.CELLTYPIST,
            CrosscheckMethod.REFERENCE_MODEL_SCARCHES,
            CrosscheckMethod.CENSUS_REFERENCE_KNN,
            CrosscheckMethod.SCIMILARITY,
        )
    return (
        CrosscheckMethod.EXISTING_PREDICTION,
        CrosscheckMethod.REFERENCE_MODEL_SCARCHES,
        gpu,
        cpu,
        CrosscheckMethod.CELLTYPIST,
        CrosscheckMethod.AZIMUTH_COMPATIBLE,
        CrosscheckMethod.SCIMILARITY,
        CrosscheckMethod.CENSUS_REFERENCE_KNN,
    )


def route_crosscheck(request: RoutingRequest) -> RoutingDecision:
    """Choose exactly one executable primary, or fail honestly with all skip reasons."""

    assessments_by_method = {method: _assess(method, request) for method in ALL_METHODS}
    if request.scheduler_ambiguity:
        blocker = (
            "scheduler/login-node ambiguity is unresolved; confirm scheduler or local_cpu before "
            "selecting a compute route"
        )
        return RoutingDecision(
            level=request.level,
            parent_scope=request.parent_scope,
            status="blocked",
            primary_method=None,
            primary_resource=None,
            assessments=tuple(assessments_by_method[method] for method in ALL_METHODS),
            blockers=(blocker,),
            next_action="Resolve the compute-profile ambiguity and rerun routing.",
        )

    selected_method = next(
        (
            method
            for method in _priority(request)
            if assessments_by_method[method].executable
            and assessments_by_method[method].suitable
        ),
        None,
    )
    if selected_method is None:
        return RoutingDecision(
            level=request.level,
            parent_scope=request.parent_scope,
            status="unavailable",
            primary_method=None,
            primary_resource=None,
            assessments=tuple(assessments_by_method[method] for method in ALL_METHODS),
            blockers=("no executable context-compatible automated cross-check is available",),
            next_action=(
                "Provide an existing prediction, compatible local labeled reference/model, or "
                "explicitly authorize and provision one optional provider; marker evidence may "
                "still yield Unknown without an automated cross-check."
            ),
        )

    primary = assessments_by_method[selected_method]
    assessments: list[MethodAssessment] = []
    for method in ALL_METHODS:
        assessment = assessments_by_method[method]
        if method is selected_method:
            assessment = assessment.model_copy(update={"selected": True})
        elif assessment.executable and assessment.suitable:
            assessment = assessment.model_copy(
                update={
                    "reasons": (
                        *assessment.reasons,
                        f"not selected because {selected_method.value} has higher deterministic "
                        f"priority for {request.technology.value}",
                    )
                }
            )
        assessments.append(assessment)
    return RoutingDecision(
        level=request.level,
        parent_scope=request.parent_scope,
        status="selected",
        primary_method=selected_method,
        primary_resource=primary.resource,
        assessments=tuple(assessments),
        next_action=(
            f"Run {selected_method.value} as the sole primary automated hypothesis for "
            f"{request.level} within {request.parent_scope}; preserve all other evidence as "
            "independent disagreement."
        ),
    )


def route_hierarchy(requests: tuple[RoutingRequest, ...]) -> tuple[RoutingDecision, ...]:
    """Route multiple level/parent scopes without allowing duplicate state keys."""

    keys = [(request.level, request.parent_scope) for request in requests]
    if len(set(keys)) != len(keys):
        raise ValueError("hierarchy routing requests duplicate a level/parent scope")
    return tuple(route_crosscheck(request) for request in requests)


__all__ = [
    "ALL_METHODS",
    "CrosscheckMethod",
    "HealthContext",
    "MethodAssessment",
    "MethodResource",
    "RoutingDecision",
    "RoutingRequest",
    "Technology",
    "route_crosscheck",
    "route_hierarchy",
]

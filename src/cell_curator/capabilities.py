"""Assay-aware capability routing with scientifically explicit claim scopes."""

from __future__ import annotations

from .config import CellCuratorConfig
from .models import AssayKind, Capability, CapabilityLevel

GPU_BACKENDS = {"rapids_singlecell", "rapids-singlecell"}


def _route_level(
    *,
    kind: AssayKind,
    feature_space: str,
    ontology_compatible: bool,
    has_mapping: bool,
) -> tuple[CapabilityLevel, tuple[str, ...], tuple[str, ...]]:
    feature_space = feature_space.lower()
    if kind is AssayKind.RNA and feature_space in {"gene", "transcript"}:
        return (
            CapabilityLevel.ANNOTATION,
            ("cell_identity", "subtype", "state"),
            ("direct gene-expression identity evidence",),
        )
    if kind is AssayKind.GENE_ACTIVITY:
        return (
            CapabilityLevel.CORROBORATION,
            ("cell_identity", "regulatory_state"),
            ("gene activity is indirect regulatory evidence",),
        )
    if kind is AssayKind.PROTEIN:
        return (
            CapabilityLevel.ANNOTATION
            if ontology_compatible
            else CapabilityLevel.CORROBORATION,
            ("cell_identity", "subtype", "state"),
            (
                "panel-aware protein identity evidence"
                if ontology_compatible
                else "protein evidence lacks an ontology-compatible signed panel",
            ),
        )
    if kind is AssayKind.ATAC:
        if feature_space in {"gene", "gene_activity"} or has_mapping:
            return (
                CapabilityLevel.CORROBORATION,
                ("cell_identity", "regulatory_state"),
                ("accessibility is indirect even after peak-to-gene or motif mapping",),
            )
        return (
            CapabilityLevel.CORROBORATION,
            ("regulatory_state", "subpopulation_separation"),
            ("peak accessibility cannot directly establish a gene-marker label",),
        )
    if kind is AssayKind.SPATIAL:
        if feature_space in {"gene", "transcript"} and ontology_compatible:
            return (
                CapabilityLevel.ANNOTATION,
                ("cell_identity", "state", "spatial_context", "neighborhood_coherence"),
                ("direct spatial expression plus coordinate-context evidence",),
            )
        return (
            CapabilityLevel.CORROBORATION,
            ("spatial_context", "neighborhood_coherence"),
            ("coordinates corroborate but do not define identity",),
        )
    return (
        CapabilityLevel.UNAVAILABLE,
        (),
        ("no executable identity or corroboration route is registered",),
    )


def build_capability_registry(config: CellCuratorConfig) -> list[Capability]:
    capabilities: list[Capability] = []
    for lane in config.evidence.lanes:
        assay = config.input.assays[lane.assay]
        representation = assay.representations[lane.representation]
        feature_space = representation.feature_space or assay.feature_space
        if representation.peak_to_gene_path:
            feature_space = "gene"
        elif representation.motif_mapping_path:
            feature_space = "motif"
        has_mapping = bool(
            representation.peak_to_gene_path or representation.motif_mapping_path
        )
        level, claims, reasons = _route_level(
            kind=AssayKind(assay.kind),
            feature_space=feature_space,
            ontology_compatible=representation.ontology_compatible,
            has_mapping=has_mapping,
        )
        capabilities.append(
            Capability(
                assay=lane.assay,
                representation=lane.representation,
                assay_kind=AssayKind(assay.kind),
                feature_space=feature_space,
                level=level,
                claims=claims,
                reasons=reasons,
                backend=representation.ranking_backend,
                requires_gpu=representation.ranking_backend.lower() in GPU_BACKENDS,
            )
        )
    return capabilities


def annotation_capable(capabilities: list[Capability]) -> bool:
    return any(item.level is CapabilityLevel.ANNOTATION for item in capabilities)

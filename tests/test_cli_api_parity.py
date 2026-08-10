"""Every CLI capability must also be reachable from Python.

The package grew CLI-first, and for a while most of its capabilities could only
be run by shelling out — which makes it unusable from a script or a notebook
without `subprocess`. The map below is the contract that keeps that from
reopening: adding a CLI action without a Python entry point fails here with a
missing-key diff rather than quietly reintroducing the gap.

The relation is deliberately not one-to-one. Several actions share one
implementation (`pipeline refine` and `cluster subcluster`; the three `markers`
scoring actions), and `config validate` is served by `load_config`, which lives
in `config.py` rather than `api.py`. So the mapping is written out rather than
derived.
"""

from __future__ import annotations

import argparse

import cell_curator
from cell_curator.cli import build_parser

# CLI action -> the public name in `cell_curator` that serves it.
CLI_TO_PUBLIC_NAME: dict[tuple[str, ...], str] = {
    ("assumptions", "add"): "record_annotation_assumption",
    ("assumptions", "approve"): "approve_annotation_assumptions",
    ("assumptions", "list"): "list_annotation_assumptions",
    ("benchmark", "evaluate"): "evaluate_benchmark_predictions",
    ("benchmark", "validate-manifest"): "validate_dataset_manifest",
    ("capabilities",): "list_capabilities",
    ("cluster", "confidence"): "advisory_confidence",
    ("cluster", "profile"): "profile_refined_clusters",
    ("cluster", "subcluster"): "discover_subcluster_candidates",
    ("config", "schema"): "configuration_schema",
    ("config", "validate"): "load_config",
    ("critic", "import-findings"): "import_critic_findings",
    ("critic", "packet"): "build_critic_packet",
    ("critic", "reconcile"): "import_critic_reconciliation",
    ("critic", "validate"): "validate_critic_reconciliation",
    ("crosscheck", "azimuth"): "map_curated_human_reference",
    ("crosscheck", "celltypist-discover"): "discover_celltypist_models",
    ("crosscheck", "celltypist-predict"): "predict_celltypist",
    ("crosscheck", "obs-pred"): "existing_prediction",
    ("crosscheck", "reference-knn"): "map_reference",
    ("crosscheck", "scimilarity"): "map_scimilarity",
    ("data", "canonicalize"): "canonicalize_input",
    ("data", "inspect"): "inspect_input",
    ("data", "qc"): "compute_input_qc",
    ("environment", "detect"): "detect_compute_environment",
    ("evidence-card", "scaffold"): "build_evidence_cards",
    ("hierarchy", "assign"): "assign_hierarchy_level",
    ("hierarchy", "validate"): "validate_hierarchy_nesting",
    ("knowledge", "validate"): "validate_knowledge_providers",
    ("markers", "assemble"): "build_marker_vocabulary",
    ("markers", "compute"): "run_marker_tools",
    ("markers", "profile"): "run_marker_tools",
    ("markers", "score"): "run_marker_tools",
    ("ontology", "provision"): "provision_ontology",
    ("ontology", "validate"): "validate_ontology_terms",
    ("pipeline", "audit"): "audit_parent_clusters",
    ("pipeline", "candidate-evidence-lane"): "run_candidate_evidence_lane",
    ("pipeline", "decide"): "propose_annotations",
    ("pipeline", "evidence"): "generate_evidence",
    ("pipeline", "evidence-lane"): "run_evidence_lane",
    ("pipeline", "map"): "map_annotated_cells",
    ("pipeline", "packet"): "freeze_review_packet",
    ("pipeline", "refine"): "discover_subcluster_candidates",
    ("pipeline", "review"): "build_review_artifacts",
    ("progress", "status"): "annotation_progress",
    ("reference", "acquire-census"): "fetch_census_reference",
    ("reference", "map-scarches"): "map_scarches",
    ("reference", "train-scanvi"): "train_scanvi_reference_model",
    ("report", "build"): "build_html_report",
    ("report", "check"): "check_review_report",
    ("route", "select"): "route_annotation_crosscheck",
    ("run", "confirm-context"): "confirm_annotation_context",
    ("run", "execute"): "annotate",
    ("run", "init"): "initialize_run",
    ("run", "plan"): "plan_annotation",
    ("run", "scene"): "prepare_annotation_scene",
    ("run", "status"): "annotation_progress",
    ("run", "validate"): "validate_run",
    ("write-back", "apply"): "apply_writeback",
    ("write-back", "preview"): "preview_writeback",
}


def _leaf_actions(
    parser: argparse.ArgumentParser, prefix: tuple[str, ...] = ()
) -> list[tuple[str, ...]]:
    subparsers = [
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    ]
    if not subparsers:
        return [prefix]
    leaves: list[tuple[str, ...]] = []
    for action in subparsers:
        for name, sub in action.choices.items():
            leaves.extend(_leaf_actions(sub, (*prefix, name)))
    return leaves


def test_every_cli_action_has_a_documented_python_entry_point() -> None:
    actions = set(_leaf_actions(build_parser()))
    mapped = set(CLI_TO_PUBLIC_NAME)
    assert actions == mapped, (
        "the CLI and the Python parity map have diverged. "
        f"CLI actions with no Python entry point: {sorted(actions - mapped)}. "
        f"Map entries for actions that no longer exist: {sorted(mapped - actions)}."
    )


def test_every_mapped_name_is_exported_and_callable() -> None:
    for action, name in sorted(CLI_TO_PUBLIC_NAME.items()):
        assert name in cell_curator.__all__, f"{name} (for {action}) is not in __all__"
        assert callable(getattr(cell_curator, name)), f"{name} is not callable"


def test_the_parity_map_covers_the_documented_action_count() -> None:
    """A floor, so a mass deletion of CLI actions cannot make parity trivial."""

    assert len(CLI_TO_PUBLIC_NAME) >= 60

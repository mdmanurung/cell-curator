"""Unified, thin command-line interface for cell-curator."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__
from .contracts import ContractError

LOG = logging.getLogger("cell_curator")


def _json(value: Any) -> None:
    if hasattr(value, "model_dump"):
        value = value.model_dump(mode="json")
    print(json.dumps(value, indent=2, sort_keys=True, default=str))


def _config_argument(parser: argparse.ArgumentParser, *, required: bool = True) -> None:
    parser.add_argument("--config", required=required, type=Path)


def _object_arguments(parser: argparse.ArgumentParser, *, name: str = "input") -> None:
    parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument(f"--{name}-modality", default="")


def _marker_arguments(parser: argparse.ArgumentParser) -> None:
    _config_argument(parser)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--representation", required=True)
    parser.add_argument("--group-key", required=True)
    parser.add_argument("--parent")
    parser.add_argument("--parent-membership-key")
    parser.add_argument("--parent-cl-id")
    parser.add_argument("--include-unscoped", action="store_true")
    parser.add_argument("--source", action="append", type=Path, default=[])
    parser.add_argument("--output-dir", type=Path)


def _reference_mapping_arguments(parser: argparse.ArgumentParser) -> None:
    _object_arguments(parser, name="query")
    _object_arguments(parser, name="reference")
    parser.add_argument("--reference-label", required=True)
    parser.add_argument("--query-cluster", default="")
    parser.add_argument("--query-layer", default="")
    parser.add_argument("--reference-layer", default="")
    parser.add_argument("--query-symbol-column", default="")
    parser.add_argument("--reference-symbol-column", default="")
    parser.add_argument("--query-symbol-map", type=Path)
    parser.add_argument("--reference-symbol-map", type=Path)
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument("--components", type=int, default=50)
    parser.add_argument("--max-memory-mb", type=int, default=8192)
    parser.add_argument("--output-dir", required=True, type=Path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cell-curator",
        description="Independent evidence-driven multimodal cell-type annotation",
    )
    parser.add_argument("--version", action="version", version=f"cell-curator {__version__}")
    parser.add_argument(
        "--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR")
    )
    commands = parser.add_subparsers(dest="command", required=True)

    config_command = commands.add_parser("config", help="validate strict configuration")
    config_sub = config_command.add_subparsers(dest="action", required=True)
    validate_config = config_sub.add_parser("validate", help="validate a configuration")
    _config_argument(validate_config)
    schema_config = config_sub.add_parser("schema", help="emit the configuration schema")
    schema_config.add_argument("--output", type=Path)

    capability_command = commands.add_parser(
        "capabilities", help="show executable assay capability routing"
    )
    _config_argument(capability_command)

    environment = commands.add_parser(
        "environment", help="detect compute resources without changing the environment"
    )
    environment_sub = environment.add_subparsers(dest="action", required=True)
    environment_detect = environment_sub.add_parser(
        "detect", help="detect CPU, RAM, GPU, and scheduler"
    )
    _config_argument(environment_detect, required=False)
    environment_detect.add_argument("--no-cuda-validation", action="store_true")

    route = commands.add_parser("route", help="select one deterministic cross-check route")
    route_sub = route.add_subparsers(dest="action", required=True)
    route_select = route_sub.add_parser(
        "select", help="route one hierarchy level from a JSON request"
    )
    route_select.add_argument("--request", required=True, type=Path)

    data_command = commands.add_parser(
        "data", help="inspect, canonicalize, or compute source-preserving QC"
    )
    data_sub = data_command.add_subparsers(dest="action", required=True)
    for action, help_text in (
        ("inspect", "inspect AnnData or MuData contracts"),
        ("canonicalize", "create or validate a canonical working copy"),
        ("qc", "compute or validate read-only QC artifacts"),
    ):
        action_parser = data_sub.add_parser(action, help=help_text)
        _config_argument(action_parser)
        if action != "inspect":
            action_parser.add_argument("--run-root", type=Path)

    knowledge_command = commands.add_parser(
        "knowledge", help="validate configured local knowledge providers"
    )
    knowledge_sub = knowledge_command.add_subparsers(dest="action", required=True)
    knowledge_validate = knowledge_sub.add_parser("validate")
    _config_argument(knowledge_validate)

    ontology_command = commands.add_parser(
        "ontology", help="validate local Cell Ontology names, synonyms, and accessions"
    )
    ontology_sub = ontology_command.add_subparsers(dest="action", required=True)
    ontology_validate = ontology_sub.add_parser("validate")
    ontology_validate.add_argument("--ontology", required=True, type=Path)
    ontology_validate.add_argument("--term", action="append", default=[])
    ontology_validate.add_argument("--expected-name", action="append", default=[])
    ontology_provision = ontology_sub.add_parser(
        "provision", help="provision a configured checksum-pinned optional ontology cache"
    )
    _config_argument(ontology_provision)

    run_command = commands.add_parser("run", help="operate the three-stage guided workflow")
    run_sub = run_command.add_subparsers(dest="action", required=True)
    for action in ("scene", "init", "plan", "status", "validate"):
        action_parser = run_sub.add_parser(action)
        _config_argument(action_parser)
    execute = run_sub.add_parser(
        "execute",
        help="execute the standard pipeline or a recursive score manifest",
    )
    _config_argument(execute)
    execute.add_argument(
        "--score-manifest",
        type=Path,
        help="strict JSON manifest for arbitrary-depth parent-scoped execution",
    )
    confirm = run_sub.add_parser("confirm-context", help="confirm or correct scene context")
    _config_argument(confirm)
    confirm.add_argument("--input", required=True, type=Path)
    confirm.add_argument("--reviewer", required=True)

    assumptions = commands.add_parser("assumptions", help="manage the living assumptions gate")
    assumptions_sub = assumptions.add_subparsers(dest="action", required=True)
    assumptions_list = assumptions_sub.add_parser("list")
    _config_argument(assumptions_list)
    assumptions_add = assumptions_sub.add_parser("add")
    _config_argument(assumptions_add)
    assumptions_add.add_argument("--scope", default="global")
    assumptions_add.add_argument("--statement", required=True)
    assumptions_add.add_argument("--non-material", action="store_true")
    assumptions_approve = assumptions_sub.add_parser("approve")
    _config_argument(assumptions_approve)
    assumptions_approve.add_argument("--reviewer", required=True)
    assumptions_approve.add_argument("--scope", default="*")
    assumptions_approve.add_argument("--preauthorized", action="store_true")

    pipeline_command = commands.add_parser(
        "pipeline", help="execute one immutable workflow phase"
    )
    pipeline_sub = pipeline_command.add_subparsers(dest="action", required=True)
    for action in ("audit", "refine", "evidence", "decide", "map", "review", "packet"):
        action_parser = pipeline_sub.add_parser(action)
        _config_argument(action_parser)
    lane_parser = pipeline_sub.add_parser("evidence-lane")
    _config_argument(lane_parser)
    lane_parser.add_argument("--assay", required=True)
    lane_parser.add_argument("--representation", required=True)
    candidate_lane_parser = pipeline_sub.add_parser("candidate-evidence-lane")
    _config_argument(candidate_lane_parser)
    candidate_lane_parser.add_argument("--lane-id", required=True)

    markers = commands.add_parser(
        "markers", help="assemble, profile, and score signed marker programs"
    )
    markers_sub = markers.add_subparsers(dest="action", required=True)
    marker_assemble = markers_sub.add_parser(
        "assemble", help="assemble local/reference/model-export programs"
    )
    _config_argument(marker_assemble)
    marker_assemble.add_argument("--source", action="append", type=Path, default=[])
    marker_assemble.add_argument("--parent")
    marker_assemble.add_argument("--parent-cl-id")
    marker_assemble.add_argument("--include-unscoped", action="store_true")
    marker_assemble.add_argument("--axis", choices=("cell_type", "cell_state"))
    marker_assemble.add_argument("--feature-space")
    marker_assemble.add_argument("--output", type=Path)
    for action, help_text in (
        ("compute", "compute bottom-up differential markers and specificity"),
        ("profile", "profile top-down means, detection, and panel capability"),
        ("score", "compute AUCell, UCell, DE overlap, and consolidated scores"),
    ):
        marker_parser = markers_sub.add_parser(action, help=help_text)
        _marker_arguments(marker_parser)

    crosscheck = commands.add_parser("crosscheck", help="run independent automated hypotheses")
    crosscheck_sub = crosscheck.add_subparsers(dest="action", required=True)
    obs_pred = crosscheck_sub.add_parser(
        "obs-pred", help="summarize an existing prediction column"
    )
    _object_arguments(obs_pred)
    obs_pred.add_argument("--prediction-column", required=True)
    obs_pred.add_argument("--cluster-column", required=True)
    obs_pred.add_argument("--output-dir", required=True, type=Path)
    typist_discover = crosscheck_sub.add_parser(
        "celltypist-discover", help="discover CellTypist models"
    )
    typist_discover.add_argument("--output-dir", type=Path)
    typist_predict = crosscheck_sub.add_parser(
        "celltypist-predict", help="run CellTypist as a hypothesis"
    )
    _object_arguments(typist_predict)
    typist_predict.add_argument("--model", required=True)
    typist_predict.add_argument("--cluster-column", default="")
    typist_predict.add_argument("--output-dir", required=True, type=Path)
    reference_knn = crosscheck_sub.add_parser(
        "reference-knn", help="map a labeled reference by CPU or validated GPU kNN"
    )
    _reference_mapping_arguments(reference_knn)
    reference_knn.add_argument("--backend", choices=("cpu", "gpu"), default="cpu")
    azimuth = crosscheck_sub.add_parser(
        "azimuth", help="map an optional curated human reference equivalent"
    )
    _reference_mapping_arguments(azimuth)
    scimilarity = crosscheck_sub.add_parser(
        "scimilarity", help="run an optional SCimilarity-compatible model"
    )
    _object_arguments(scimilarity)
    scimilarity.add_argument("--model-path", required=True, type=Path)
    scimilarity.add_argument("--output-dir", required=True, type=Path)

    reference = commands.add_parser(
        "reference", help="acquire, train, or map reference resources"
    )
    reference_sub = reference.add_subparsers(dest="action", required=True)
    census = reference_sub.add_parser(
        "acquire-census", help="acquire a bounded versioned Census subset"
    )
    census.add_argument("--organism", required=True)
    census.add_argument("--tissue", required=True)
    census.add_argument("--output", required=True, type=Path)
    census.add_argument("--version", required=True)
    census.add_argument("--max-cells", type=int, default=50_000)
    census.add_argument("--label-column", action="append", default=[])
    census.add_argument("--expected-sha256", default="")
    census.add_argument("--timeout-seconds", type=float, default=300.0)
    train = reference_sub.add_parser(
        "train-scanvi", help="train scVI/scANVI under the real-GPU contract"
    )
    _object_arguments(train, name="reference")
    train.add_argument("--labels-key", required=True)
    train.add_argument("--counts-layer", default="counts")
    train.add_argument("--max-epochs", type=int, default=200)
    train.add_argument("--seed", type=int, default=0)
    train.add_argument("--output-dir", required=True, type=Path)
    scarches = reference_sub.add_parser(
        "map-scarches", help="map a query with scArches uncertainty"
    )
    _object_arguments(scarches, name="query")
    scarches.add_argument("--model-dir", required=True, type=Path)
    scarches.add_argument("--max-epochs", type=int, default=100)
    scarches.add_argument("--output-dir", required=True, type=Path)

    cluster = commands.add_parser(
        "cluster", help="refine, profile, and score parent-local clusters"
    )
    cluster_sub = cluster.add_subparsers(dest="action", required=True)
    subcluster = cluster_sub.add_parser(
        "subcluster", help="run stored and parent-local refinement"
    )
    _config_argument(subcluster)
    profile = cluster_sub.add_parser(
        "profile", help="summarize refined-child QC/composition/spatial context"
    )
    _config_argument(profile)
    confidence = cluster_sub.add_parser("confidence", help="calculate advisory confidence")
    confidence.add_argument("--margin", required=True, type=float)
    confidence.add_argument("--specificity", required=True, type=float)
    confidence.add_argument("--agreement", required=True, type=float)
    confidence.add_argument("--coverage", required=True, type=float)

    hierarchy = commands.add_parser(
        "hierarchy", help="assign and validate strict parent-scoped labels"
    )
    hierarchy_sub = hierarchy.add_subparsers(dest="action", required=True)
    hierarchy_assign = hierarchy_sub.add_parser(
        "assign", help="assign one level inside one approved parent scope"
    )
    _config_argument(hierarchy_assign)
    hierarchy_assign.add_argument("--scores", required=True, type=Path)
    hierarchy_assign.add_argument("--vocabulary", type=Path)
    hierarchy_assign.add_argument("--level", required=True)
    hierarchy_assign.add_argument("--parent", required=True)
    hierarchy_assign.add_argument("--state-json", type=Path)
    hierarchy_assign.add_argument("--output", type=Path)
    hierarchy_validate = hierarchy_sub.add_parser(
        "validate", help="validate strict nesting and one-cell-one-leaf"
    )
    hierarchy_validate.add_argument("--labels", action="append", required=True, type=Path)
    hierarchy_validate.add_argument("--cell-membership", type=Path)

    cards = commands.add_parser(
        "evidence-card", help="build evidence cards from actual durable values"
    )
    cards_sub = cards.add_subparsers(dest="action", required=True)
    cards_scaffold = cards_sub.add_parser("scaffold")
    _config_argument(cards_scaffold)

    report_command = commands.add_parser(
        "report", help="build or strictly check self-contained HTML reports"
    )
    report_sub = report_command.add_subparsers(dest="action", required=True)
    report_build = report_sub.add_parser("build")
    _config_argument(report_build)
    report_build.add_argument("--strict", action="store_true")
    report_build.add_argument("--level")
    report_check = report_sub.add_parser("check")
    _config_argument(report_check)
    report_check.add_argument("--level")
    report_check.add_argument(
        "--strict", action="store_true", help="required completeness mode (always enforced)"
    )

    critic_command = commands.add_parser(
        "critic", help="exchange and reconcile platform-neutral critic packets"
    )
    critic_sub = critic_command.add_subparsers(dest="action", required=True)
    critic_packet = critic_sub.add_parser("packet")
    _config_argument(critic_packet)
    critic_findings = critic_sub.add_parser("import-findings")
    _config_argument(critic_findings)
    critic_findings.add_argument("--input", required=True, type=Path)
    critic_reconcile = critic_sub.add_parser("reconcile")
    _config_argument(critic_reconcile)
    critic_reconcile.add_argument("--input", required=True, type=Path)
    critic_validate = critic_sub.add_parser("validate")
    _config_argument(critic_validate)

    progress = commands.add_parser(
        "progress", help="show current stage, scope, blockers, and next action"
    )
    progress_sub = progress.add_subparsers(dest="action", required=True)
    progress_status = progress_sub.add_parser("status")
    _config_argument(progress_status)

    write_command = commands.add_parser(
        "write-back", help="preview or apply approved atomic write-back"
    )
    write_sub = write_command.add_subparsers(dest="action", required=True)
    preview = write_sub.add_parser("preview")
    _config_argument(preview)
    preview.add_argument("--ontology", type=Path)
    preview.add_argument("--output", type=Path)
    preview.add_argument("--in-place", action="store_true")
    apply = write_sub.add_parser("apply")
    _config_argument(apply)
    apply.add_argument("--output", type=Path)
    apply.add_argument("--approval", required=True, type=Path)
    apply.add_argument("--ontology", type=Path)
    apply.add_argument("--in-place", action="store_true")

    benchmark_command = commands.add_parser(
        "benchmark", help="evaluate neutral exported predictions"
    )
    benchmark_sub = benchmark_command.add_subparsers(dest="action", required=True)
    evaluate = benchmark_sub.add_parser("evaluate")
    evaluate.add_argument("--predictions", required=True, type=Path)
    evaluate.add_argument("--system", required=True)
    evaluate.add_argument("--dataset", required=True)
    evaluate.add_argument("--output", required=True, type=Path)
    manifest = benchmark_sub.add_parser("validate-manifest")
    manifest.add_argument("--manifest", required=True, type=Path)
    return parser


def _marker_run(args: argparse.Namespace) -> Any:
    from .markers import run_marker_analysis_from_config

    result = run_marker_analysis_from_config(
        args.config,
        assay=args.assay,
        representation=args.representation,
        group_key=args.group_key,
        parent_scope=args.parent,
        parent_membership_key=args.parent_membership_key,
        parent_cl_id=args.parent_cl_id,
        include_unscoped_programs=args.include_unscoped,
        additional_sources=args.source,
        output_dir=args.output_dir,
    )
    selected = {
        "compute": ("bottom_up", "de_overlap"),
        "profile": (
            "feature_profile",
            "mean_matrix",
            "detection_matrix",
            "panel_capability",
            "program_summary",
            "negative_violations",
            "dotplot",
        ),
        "score": ("aucell", "ucell", "de_overlap", "crosschecks"),
    }[args.action]
    return {
        "status": "complete",
        "action": args.action,
        "hypothesis_not_verdict": True,
        "artifacts": {
            key: str(result.artifacts[key]) for key in selected if key in result.artifacts
        },
        "manifest": str(result.artifacts["manifest"]),
    }


def _reference_kwargs(args: argparse.Namespace) -> dict[str, Any]:
    from .crosschecks import load_symbol_map

    return {
        "query_path": args.query,
        "reference_path": args.reference,
        "reference_label_column": args.reference_label,
        "query_cluster_column": args.query_cluster,
        "query_modality": args.query_modality,
        "reference_modality": args.reference_modality,
        "output_dir": args.output_dir,
        "query_layer": args.query_layer,
        "reference_layer": args.reference_layer,
        "query_symbol_column": args.query_symbol_column,
        "reference_symbol_column": args.reference_symbol_column,
        "query_symbol_map": load_symbol_map(args.query_symbol_map),
        "reference_symbol_map": load_symbol_map(args.reference_symbol_map),
        "k": args.k,
        "components": args.components,
        "max_memory_mb": args.max_memory_mb,
    }


def dispatch(args: argparse.Namespace) -> Any:
    if args.command == "config":
        from .config import configuration_schema, load_config
        from .provenance import publish_json

        if args.action == "schema":
            schema = configuration_schema()
            if args.output is not None:
                publish_json(args.output, schema)
                return {"schema_version": 3, "path": str(args.output)}
            return schema
        return load_config(args.config).model_dump(mode="json")
    if args.command == "capabilities":
        from .capabilities import build_capability_registry
        from .config import load_config

        return [
            item.model_dump(mode="json")
            for item in build_capability_registry(load_config(args.config))
        ]
    if args.command == "environment":
        from .environment import detect_environment

        cpu: dict[str, Any] = {}
        gpu: dict[str, Any] = {}
        if args.config:
            from .config import load_config

            config = load_config(args.config)
            cpu, gpu = config.hpc.cpu_resources, config.hpc.gpu_resources
        return detect_environment(
            configured_cpu_resources=cpu,
            configured_gpu_resources=gpu,
            validate_cuda=not args.no_cuda_validation,
        )
    if args.command == "route":
        from .routing import RoutingRequest, route_crosscheck

        request = RoutingRequest.model_validate(json.loads(args.request.read_text()))
        return route_crosscheck(request)
    if args.command == "data":
        from .config import load_config
        from .data import canonicalize, compute_qc, inspect
        from .pipeline import run_root

        config = load_config(args.config)
        if args.action == "inspect":
            return inspect(config)
        root = args.run_root or run_root(config)
        if args.action == "canonicalize":
            return canonicalize(config, config_path=args.config, run_root=root)
        manifest = root / "qc" / "qc.manifest.json"
        if manifest.is_file():
            return compute_qc(config, run_root=root)
        result = canonicalize(config, config_path=args.config, run_root=root)
        return result["qc"]
    if args.command == "knowledge":
        from .api import validate_knowledge_providers

        return validate_knowledge_providers(args.config)
    if args.command == "ontology":
        if args.action == "provision":
            from .api import provision_ontology

            return provision_ontology(args.config)

        from .api import validate_ontology_terms

        # A flag-worded arity check belongs to the CLI presentation; the shared
        # implementation carries its own argument-worded equivalent.
        if args.expected_name and len(args.expected_name) != len(args.term):
            raise ContractError("--expected-name must be repeated once per --term")
        result = validate_ontology_terms(
            args.ontology,
            args.term,
            expected_names=args.expected_name or [None] * len(args.term),
        )
        return {**result, "terms": [term.__dict__ for term in result["terms"]]}
    if args.command == "run":
        from .config import load_config
        from .guided import confirm_context, construct_plan, prepare_scene
        from .guided import status as guided_status
        from .pipeline import initialize, run_annotation, validate

        if args.action == "scene":
            return prepare_scene(args.config)
        if args.action == "init":
            return initialize(args.config)
        if args.action == "confirm-context":
            return confirm_context(args.config, input_path=args.input, reviewer=args.reviewer)
        if args.action == "plan":
            return construct_plan(args.config)
        if args.action == "execute":
            if args.score_manifest is not None:
                from .pipeline import execute_recursive_hierarchy_manifest

                return execute_recursive_hierarchy_manifest(args.config, args.score_manifest)
            return {"review_packet": str(run_annotation(args.config))}
        if args.action == "status":
            return guided_status(args.config)
        load_config(args.config)
        return validate(args.config)
    if args.command == "assumptions":
        from .guided import approve_assumptions

        if args.action == "list":
            from .api import list_annotation_assumptions

            return list_annotation_assumptions(args.config)
        if args.action == "add":
            from .api import record_annotation_assumption

            return record_annotation_assumption(
                args.config,
                scope=args.scope,
                statement=args.statement,
                material=not args.non_material,
            )
        return approve_assumptions(
            args.config,
            reviewer=args.reviewer,
            scope=args.scope,
            preauthorized=args.preauthorized,
        )
    if args.command == "pipeline":
        from .pipeline import (
            audit,
            build_review_packet,
            candidate_evidence_lane,
            decide,
            evidence_all,
            evidence_lane,
            map_cells,
            refine,
            review,
        )

        if args.action == "audit":
            return {"n_parents": len(audit(args.config))}
        if args.action == "refine":
            registry, _, evaluations = refine(args.config)
            return {"n_candidates": len(registry), "n_evaluations": len(evaluations)}
        if args.action == "evidence-lane":
            return {
                "receipt": str(
                    evidence_lane(
                        args.config, assay=args.assay, representation=args.representation
                    )
                )
            }
        if args.action == "candidate-evidence-lane":
            return {
                "receipt": str(
                    candidate_evidence_lane(args.config, comparison_lane_id=args.lane_id)
                )
            }
        if args.action == "evidence":
            scores, proposals = evidence_all(args.config)
            return {"n_scores": len(scores), "n_proposals": len(proposals)}
        if args.action == "decide":
            children, parents = decide(args.config)
            return {"n_children": len(children), "n_parents": len(parents)}
        if args.action == "map":
            return {"n_cells": len(map_cells(args.config))}
        if args.action == "review":
            return {"report": str(review(args.config))}
        return {"review_packet": str(build_review_packet(args.config))}
    if args.command == "markers":
        if args.action != "assemble":
            return _marker_run(args)
        from .api import build_marker_vocabulary

        programs, output = build_marker_vocabulary(
            args.config,
            source=args.source,
            parent_scope=args.parent,
            parent_cl_id=args.parent_cl_id,
            include_unscoped=args.include_unscoped,
            axis=args.axis,
            feature_space=args.feature_space,
            output_path=args.output,
        )
        return {"status": "complete", "n_programs": len(programs), "output": str(output)}
    if args.command == "crosscheck":
        from .crosschecks import (
            discover_celltypist_models,
            existing_prediction_file,
            map_curated_human_reference_files,
            map_reference_files,
            map_scimilarity_file,
            predict_celltypist_file,
        )

        if args.action == "obs-pred":
            cells, clusters, status = existing_prediction_file(
                input_path=args.input,
                modality=args.input_modality,
                prediction_column=args.prediction_column,
                cluster_column=args.cluster_column,
                output_dir=args.output_dir,
            )
            return {**status, "n_cells": len(cells), "n_clusters": len(clusters)}
        if args.action == "celltypist-discover":
            return discover_celltypist_models(output_dir=args.output_dir)
        if args.action == "celltypist-predict":
            cells, status = predict_celltypist_file(
                input_path=args.input,
                modality=args.input_modality,
                model=args.model,
                cluster_column=args.cluster_column,
                output_dir=args.output_dir,
            )
            return {**status, "n_cells": len(cells)}
        if args.action == "reference-knn":
            cells, clusters, status = map_reference_files(
                **_reference_kwargs(args), backend=args.backend
            )
            return {**status, "n_cells": len(cells), "n_clusters": len(clusters)}
        if args.action == "azimuth":
            cells, clusters, status = map_curated_human_reference_files(
                **_reference_kwargs(args), backend="cpu"
            )
            return {**status, "n_cells": len(cells), "n_clusters": len(clusters)}
        cells, status = map_scimilarity_file(
            input_path=args.input,
            modality=args.input_modality,
            model_path=args.model_path,
            output_dir=args.output_dir,
        )
        return {**status, "n_cells": len(cells)}
    if args.command == "reference":
        from .crosschecks import (
            fetch_census_reference,
            map_scarches_file,
            train_scanvi_reference_file,
        )

        if args.action == "acquire-census":
            return fetch_census_reference(
                organism=args.organism,
                tissue=args.tissue,
                output_path=args.output,
                census_version=args.version,
                max_cells=args.max_cells,
                label_columns=args.label_column or ("cell_type",),
                expected_sha256=args.expected_sha256,
                timeout_seconds=args.timeout_seconds,
            )
        if args.action == "train-scanvi":
            return train_scanvi_reference_file(
                reference_path=args.reference,
                modality=args.reference_modality,
                labels_key=args.labels_key,
                output_dir=args.output_dir,
                counts_layer=args.counts_layer,
                max_epochs=args.max_epochs,
                seed=args.seed,
            )
        cells, status = map_scarches_file(
            query_path=args.query,
            modality=args.query_modality,
            model_dir=args.model_dir,
            output_dir=args.output_dir,
            max_epochs=args.max_epochs,
        )
        return {**status, "n_cells": len(cells)}
    if args.command == "cluster":
        if args.action == "confidence":
            from .hierarchy import advisory_confidence

            return {
                "advisory_confidence": advisory_confidence(
                    margin=args.margin,
                    specificity=args.specificity,
                    agreement=args.agreement,
                    coverage=args.coverage,
                ),
                "calibrated_probability": False,
            }
        if args.action == "subcluster":
            from .pipeline import refine

            registry, membership, evaluations = refine(args.config)
            return {
                "status": "complete",
                "n_candidates": len(registry),
                "n_membership_rows": len(membership),
                "n_evaluations": len(evaluations),
            }
        from .hierarchy import profile_refined_clusters

        return profile_refined_clusters(args.config)
    if args.command == "hierarchy":
        from .hierarchy import (
            validate_hierarchy_nesting,
        )

        if args.action == "validate":
            frames = [
                pd.read_csv(path, sep="\t", keep_default_na=False) for path in args.labels
            ]
            membership = (
                pd.read_csv(args.cell_membership, sep="\t", keep_default_na=False)
                if args.cell_membership
                else None
            )
            return validate_hierarchy_nesting(frames, cell_membership=membership)
        from .api import assign_hierarchy_level

        result = assign_hierarchy_level(
            args.config,
            scores_path=args.scores,
            level=args.level,
            parent=args.parent,
            vocabulary_path=args.vocabulary,
            state_path=args.state_json,
            output_path=args.output,
        )
        return {
            "status": "complete",
            "n_labels": len(result["assigned"]),
            "output": str(result["output"]),
        }
    if args.command == "evidence-card":
        from .api import build_evidence_cards

        outputs = build_evidence_cards(args.config)
        return {"status": "complete", "n_cards": len(outputs), "cards": list(map(str, outputs))}
    if args.command == "report":
        if args.action == "check":
            from .api import check_review_report

            # The completeness dict is the real result; the previous CLI shape
            # discarded everything but the report path.
            return check_review_report(args.config, level=args.level)
        from .api import build_html_report

        path = build_html_report(args.config, strict=args.strict, level=args.level)
        return {"status": "complete", "report": str(path)}
    if args.command == "critic":
        from .config import load_config
        from .pipeline import run_root
        from .review import (
            build_critic_packet,
            import_critic_findings,
            import_reconciliation,
            validate_reconciliation,
        )

        config = load_config(args.config)
        root = run_root(config)
        if args.action == "packet":
            return {"status": "complete", "packet": str(build_critic_packet(config, root))}
        if args.action == "import-findings":
            return import_critic_findings(root, args.input)
        if args.action == "reconcile":
            return import_reconciliation(root, args.input)
        return validate_reconciliation(root)
    if args.command == "progress":
        from .guided import status as progress_status

        return progress_status(args.config)
    if args.command == "write-back":
        from .config import load_config
        from .ontology import CellOntology, ontology_from_config
        from .pipeline import run_root
        from .writeback import validate_writeback, write_back

        config = load_config(args.config)
        root = run_root(config)
        ontology = (
            CellOntology.load(args.ontology) if args.ontology else ontology_from_config(config)
        )
        if args.action == "preview":
            return validate_writeback(
                config,
                root,
                ontology=ontology,
                output_path=args.output,
                in_place=args.in_place,
            )
        output = args.output
        if output is None:
            output = (
                Path(config.input.path)
                if args.in_place
                else root / "writeback" / f"annotated.{config.input.kind}"
            )
        return {
            "output": str(
                write_back(
                    config,
                    root,
                    output_path=output,
                    approval_path=args.approval,
                    ontology=ontology,
                    in_place=args.in_place,
                )
            )
        }
    if args.command == "benchmark":
        from .benchmark import evaluate_file, validate_dataset_manifest

        if args.action == "evaluate":
            return evaluate_file(
                args.predictions,
                system=args.system,
                dataset=args.dataset,
                output_path=args.output,
            ).model_dump(mode="json")
        return validate_dataset_manifest(args.manifest)
    raise ContractError(f"unsupported command: {args.command}")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    try:
        result = dispatch(args)
    except (ContractError, ValueError, OSError, RuntimeError) as exc:
        LOG.error("%s", exc)
        raise SystemExit(2) from exc
    if result is not None:
        _json(result)


if __name__ == "__main__":
    main()

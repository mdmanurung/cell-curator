"""Corpus-backed, parent-scoped annotation after human-directed scene and plan.

Run with --configfile CONFIG.yaml. This DAG validates existing approvals through
the Python API. It never approves assumptions or performs write-back.
"""

import re
from pathlib import Path

import yaml
from snakemake.exceptions import WorkflowError

from cell_curator.config import load_config


if len(workflow.configfiles) != 1:
    raise WorkflowError("pass exactly one strict run configuration with --configfile CONFIG.yaml")

CONFIG = Path(workflow.configfiles[0]).resolve()
if yaml.safe_load(CONFIG.read_text()) != config:
    raise WorkflowError("Snakemake config overlays are prohibited; use one materialized configuration")
ANCHORED = load_config(CONFIG)
if not ANCHORED.annotation.enabled:
    raise WorkflowError("workflow/annotation.smk requires annotation.enabled=true")
RUN = Path(ANCHORED.run.output_root) / ANCHORED.run.run_id
LEVELS = ANCHORED.guidance.hierarchy_levels
CLASSIFIERS = {item.name: item for item in ANCHORED.annotation.classifiers}
CLASSIFIER_STATUSES = [RUN / "annotation" / "classifiers" / name / "status.json" for name in CLASSIFIERS]
CLASSIFIER_CELLS = [RUN / "annotation" / "classifiers" / name / "cell_predictions.tsv" for name in CLASSIFIERS]
DECISIONS = [RUN / "annotation" / "levels" / level / "decision.manifest.json" for level in LEVELS]
LEVEL_LABELS = [RUN / "annotation" / "levels" / level / "labels.tsv" for level in LEVELS]
CONSENSUS_DIRECTORY = Path(ANCHORED.annotation.consensus_directory) if ANCHORED.annotation.consensus_directory else None
if CONSENSUS_DIRECTORY is None:
    raise WorkflowError("set annotation.consensus_directory to the output of an annotation-consensus run")
CONSENSUS = [CONSENSUS_DIRECTORY / name for name in ("consensus.tsv", "run_manifest.json", "work/cells.tsv")]
CANONICAL = RUN / "canonical" / f"input.{ANCHORED.input.kind}"
EXTERNAL = [RUN / "context.json", RUN / "plan.json", RUN / "assumptions.json", RUN / "input_contract.manifest.json"]
KNOWLEDGE = sorted({Path(provider.path) for provider in ANCHORED.knowledge.providers if provider.path})
MAPPINGS = sorted({Path(assay.gene_mapping_path) for assay in ANCHORED.input.assays.values() if assay.gene_mapping_path})


def classifier_inputs(wildcards):
    item = CLASSIFIERS[str(wildcards.classifier)]
    model = Path(item.model)
    models = sorted(path for path in model.rglob("*") if path.is_file()) if model.is_dir() else []
    return [*(models or [model]), *([Path(item.reference)] if item.reference else []),
            *([Path(item.python_executable)] if item.python_executable else [])]


def previous_decision(wildcards):
    index = LEVELS.index(str(wildcards.level))
    return [DECISIONS[index - 1], LEVEL_LABELS[index - 1]] if index else []


rule annotation:
    input:
        RUN / "annotation" / "annotation.manifest.json",
        RUN / "annotation" / "cell_labels.tsv",
        RUN / "annotation" / "report.html",
        RUN / "hierarchy" / "execution.manifest.json",
        RUN / "mapping" / "parent_child_mapping.tsv",
        RUN / "mapping" / "mapping.validation.json",


rule annotation_review:
    input:
        RUN / "review" / "report.html",
        RUN / "review" / "report.manifest.json",
        RUN / "review" / "colors.json",
        RUN / "review" / "critic_packet.json",
        RUN / "review" / "evidence_cards",


rule annotation_classifier:
    input:
        config=CONFIG,
        source=ANCHORED.input.path,
        canonical=CANONICAL,
        gates=EXTERNAL,
        knowledge=KNOWLEDGE,
        mappings=MAPPINGS,
        models=classifier_inputs,
    output:
        status=RUN / "annotation" / "classifiers" / "{classifier}" / "status.json",
        cells=RUN / "annotation" / "classifiers" / "{classifier}" / "cell_predictions.tsv",
    wildcard_constraints:
        classifier="|".join(re.escape(name) for name in CLASSIFIERS) if CLASSIFIERS else "a^",
    run:
        from cell_curator.annotation import run_classifier
        run_classifier(Path(str(input.config)), str(wildcards.classifier))


rule annotation_level_evidence:
    input:
        config=CONFIG,
        source=ANCHORED.input.path,
        canonical=CANONICAL,
        gates=EXTERNAL,
        knowledge=KNOWLEDGE,
        mappings=MAPPINGS,
        statuses=CLASSIFIER_STATUSES,
        predictions=CLASSIFIER_CELLS,
        parent=previous_decision,
    output:
        manifest=RUN / "annotation" / "levels" / "{level}" / "evidence.manifest.json",
        markers=RUN / "annotation" / "levels" / "{level}" / "markers.tsv",
        scores=RUN / "annotation" / "levels" / "{level}" / "scores.tsv",
        states=RUN / "annotation" / "levels" / "{level}" / "states.tsv",
        plot=RUN / "annotation" / "levels" / "{level}" / "marker_evidence.svg",
    wildcard_constraints:
        level="|".join(re.escape(level) for level in LEVELS),
    run:
        from cell_curator.annotation import collect_level
        collect_level(Path(str(input.config)), str(wildcards.level))


rule annotation_level_decision:
    input:
        config=CONFIG,
        manifest=RUN / "annotation" / "levels" / "{level}" / "evidence.manifest.json",
        markers=RUN / "annotation" / "levels" / "{level}" / "markers.tsv",
        scores=RUN / "annotation" / "levels" / "{level}" / "scores.tsv",
        states=RUN / "annotation" / "levels" / "{level}" / "states.tsv",
        plot=RUN / "annotation" / "levels" / "{level}" / "marker_evidence.svg",
        gates=EXTERNAL,
        knowledge=KNOWLEDGE,
        predictions=CLASSIFIER_CELLS,
        statuses=CLASSIFIER_STATUSES,
        parent=previous_decision,
        consensus=CONSENSUS,
    output:
        labels=RUN / "annotation" / "levels" / "{level}" / "labels.tsv",
        manifest=RUN / "annotation" / "levels" / "{level}" / "decision.manifest.json",
    wildcard_constraints:
        level="|".join(re.escape(level) for level in LEVELS),
    run:
        from cell_curator.consensus_import import import_consensus_level
        import_consensus_level(Path(str(input.config)), str(wildcards.level))


rule annotation_finalize:
    input:
        config=CONFIG,
        source=ANCHORED.input.path,
        canonical=CANONICAL,
        gates=EXTERNAL,
        knowledge=KNOWLEDGE,
        decisions=DECISIONS,
        labels=LEVEL_LABELS,
        predictions=CLASSIFIER_CELLS,
        statuses=CLASSIFIER_STATUSES,
    output:
        manifest=RUN / "annotation" / "annotation.manifest.json",
        cells=RUN / "annotation" / "cell_labels.tsv",
        report=RUN / "annotation" / "report.html",
        hierarchy=RUN / "hierarchy" / "execution.manifest.json",
        mapping=RUN / "mapping" / "parent_child_mapping.tsv",
        mapping_validation=RUN / "mapping" / "mapping.validation.json",
    run:
        from cell_curator.annotation import finalize_annotation
        finalize_annotation(Path(str(input.config)))


rule annotation_build_review:
    input:
        config=CONFIG,
        gates=EXTERNAL,
        source=ANCHORED.input.path,
        knowledge=KNOWLEDGE,
        final=rules.annotation_finalize.output.manifest,
        cells=rules.annotation_finalize.output.cells,
        mapping=rules.annotation_finalize.output.mapping,
        mapping_validation=rules.annotation_finalize.output.mapping_validation,
        hierarchy=rules.annotation_finalize.output.hierarchy,
        decisions=DECISIONS,
        labels=LEVEL_LABELS,
        critic_ledgers=lambda wildcards: [path for path in (
            RUN / "review" / "critic_findings.tsv", RUN / "review" / "critic_reconciliation.tsv"
        ) if path.is_file()],
    output:
        report=RUN / "review" / "report.html",
        manifest=RUN / "review" / "report.manifest.json",
        colors=RUN / "review" / "colors.json",
        packet=RUN / "review" / "critic_packet.json",
        cards=directory(RUN / "review" / "evidence_cards"),
    run:
        from cell_curator.review import build_evidence_cards, build_critic_packet, build_html_report
        build_evidence_cards(ANCHORED, RUN)
        build_critic_packet(ANCHORED, RUN)
        build_html_report(ANCHORED, RUN, strict=True)

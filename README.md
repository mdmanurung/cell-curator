# cell-curator

`cell-curator` is an independent, local-first Python package and Codex skill for
evidence-driven annotation of pre-clustered AnnData and MuData objects. It combines
signed cell-type and cell-state programs, bottom-up markers, ontology validation,
assay-aware multimodal evidence, deterministic uncertainty and `Unknown` calls,
adaptive parent-local refinement, immutable provenance, human review, and guarded
write-back.

The package has no dependency on another annotation package, CLI, agent platform,
MCP, atlas, or web service. Local marker, ontology, cached-table, and labeled-reference
providers are sufficient for complete offline runs. Optional remote providers are
checksum-pinned and cached.

## Install and inspect

```bash
uv sync --frozen --extra dev --extra workflow
uv run cell-curator --version
uv run cell-curator --help
```

Copy [`assets/config.template.yaml`](assets/config.template.yaml), declare the frozen
input and local knowledge sources, and run the three guided stages:

```bash
uv run cell-curator config validate --config config.yaml
uv run cell-curator environment detect --config config.yaml
uv run cell-curator data inspect --config config.yaml
uv run cell-curator data canonicalize --config config.yaml
uv run cell-curator data qc --config config.yaml
uv run cell-curator run scene --config config.yaml
uv run cell-curator run confirm-context --config config.yaml \
  --input context-confirmation.json --reviewer REVIEWER
uv run cell-curator run plan --config config.yaml
uv run cell-curator assumptions approve --config config.yaml --reviewer REVIEWER
uv run cell-curator run execute --config config.yaml
```

Interactive runs stop before the first biological label until the context and
assumptions are approved. A non-interactive run must declare complete context and
pre-authorization in configuration; it records the same gate and decision history.
`run execute` runs the configured evidence/refinement pipeline and builds review
artifacts. Supply `--score-manifest hierarchy-scores.json` to execute arbitrary-depth
exact `(level, parent_scope)` score tables with scope gates, retained-parent leaves,
and content-hashed resume. Neither route writes labels back to the source object.
Reconcile critics and run strict read-only validation before the separate write-back
approval:

```bash
uv run cell-curator write-back preview --config config.yaml
uv run cell-curator write-back apply --config config.yaml \
  --approval approval.json --output annotated.h5ad
```

Bind `approval.json` to the run ID and the source, mapping, durable-label, and
preview-diff hashes emitted by `write-back preview`.

## Stable Python API

```python
from cell_curator import (
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
```

These functions and the CLI call the same package implementation. The Snakemake DAG
in [`workflow/Snakefile`](workflow/Snakefile) exposes separate `all_evidence`,
`review_packet`, `finalization`, and `writeback` targets; write-back is not the default.

See [`SKILL.md`](SKILL.md) for the agent workflow and the one-level
[`references/`](references/) documents for input, assay, refinement, evidence,
decision, artifact, compute, and validation contracts.

## Scientific boundaries

- Cell type and cell state are independent axes.
- Missing panel features are unmeasured, not negative.
- ATAC peaks are accessibility evidence; gene activity is indirect evidence.
- A configured RAPIDS lane requires the declared real GPU/software contract and fails
  closed. It never substitutes a CPU result.
- `Unknown` is a valid result when coverage, coherence, margin, or required capability
  is inadequate.
- Labels may remain CL-free when no ontology is available; supplied CL identifiers
  always require local semantic validation.
- No superiority claim is made without measured neutral benchmark results. Public
  benchmark manifests record pending external inputs and exact reproduction commands.

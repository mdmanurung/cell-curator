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

## Installation

`cell-curator` requires Python 3.11 or newer. The current release is installed from
source; no PyPI release is assumed.

### With uv (recommended)

```bash
git clone https://github.com/mdmanurung/cell-curator.git
cd cell-curator
uv sync --frozen --extra workflow
uv run cell-curator --version
uv run cell-curator --help
```

Add the development environment when contributing:

```bash
uv sync --frozen --extra dev --extra workflow
```

### With pip

```bash
git clone https://github.com/mdmanurung/cell-curator.git
cd cell-curator
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ".[workflow]"
cell-curator --version
```

The core install omits optional knowledge providers and workflow execution. Install
only the capabilities needed for the run, for example `.[celltypist]`,
`.[reference-models]`, `.[census]`, `.[scimilarity]`, or combinations such as
`.[workflow,celltypist]`. The `gpu` extra is not a portability mode: a configured
RAPIDS lane also requires a compatible real NVIDIA GPU, driver, CUDA runtime, and
declared RAPIDS contract. It fails closed when that contract is not met and never
falls back to CPU.

## Quick start

Copy the strict template, then edit the input, assay representations, context,
knowledge providers, and output run ID:

```bash
cp assets/config.template.yaml config.yaml
uv run cell-curator config validate --config config.yaml
```

Every input must be pre-clustered and immutable for the duration of a run. Set
`input.expected_sha256` when the input is supplied externally, use a new `run.run_id`
for a changed input or configuration, and keep `run.immutable: true`.

### Guided CLI run

Interactive mode is the safest default. It deliberately stops before the first
biological label until a reviewer has confirmed the context and assumptions:

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
uv run cell-curator write-back preview --config config.yaml --output annotated.h5ad
uv run cell-curator write-back apply --config config.yaml \
  --approval approval.json --output annotated.h5ad
```

Bind `approval.json` to every field emitted for approval by `write-back preview`:
the run, source, mapping, durable labels, assumptions, report, critic artifacts,
diff, exact output path, and in-place intent.

### Snakemake execution

The DAG is authoritative for dependencies and resumability. Its default aggregate
target stops at evidence; review, finalization, and write-back remain explicit:

```bash
uv run snakemake --snakefile workflow/Snakefile all_evidence \
  --configfile config.yaml --cores 4
uv run snakemake --snakefile workflow/Snakefile review_packet \
  --configfile config.yaml --cores 4
uv run snakemake --snakefile workflow/Snakefile finalization \
  --configfile config.yaml --cores 4
```

Use `--dry-run` before scheduling a new configuration. For HPC execution, adapt
[`workflow/profiles/slurm/config.yaml`](workflow/profiles/slurm/config.yaml) to the
site rather than hard-coding cluster resources in the workflow. The DAG uses the
locked environment that invoked Snakemake, which keeps an installed wheel independent
of a source-tree `-e .` path. RAPIDS lanes request one `gpu` resource and validate the
real device at runtime. The separate `writeback` target remains blocked until its
run-bound approval exists.

## Use-case guide

Choose the closest starting point below and encode only evidence the assay actually
measures. In every case, `input.assays` names the modalities and representations,
while `evidence.lanes` selects the representations allowed to support a decision.

| Use case | Configuration starting point | Important boundary |
| --- | --- | --- |
| Clustered scRNA-seq (`.h5ad`) | `input.kind: h5ad`; one root `rna` assay; declare log-normalized values for ranking and counts for detection | Counts and log values are different evidence; missing genes are unmeasured, not negative |
| CITE-seq or other multimodal data (`.h5mu`) | `input.kind: h5mu`; declare RNA and protein modalities as separate assays and evidence lanes | Preserve cross-modal disagreement for review instead of averaging it away |
| scATAC-seq or multiome | Declare `atac` peak accessibility separately from any `gene_activity` representation | Peaks are direct accessibility evidence; gene activity is indirect and cannot silently become RNA expression |
| Spatial transcriptomics | Declare the transcriptome assay and set `input.keys.spatial` to the coordinate representation | Spatial proximity is contextual evidence, not proof of cell identity |
| Evidence-only audit | Set `run.mode: evidence-only` and provide frozen cluster assignments | Produces auditable evidence and uncertainty without forcing biological labels |
| Hierarchical annotation and selective refinement | Declare `guidance.hierarchy_levels`; start from stored partitions and enable parent-local refinement only where gates pass | Every split is scoped to its parent and may validly retain the parent or remain unresolved |

### 1. RNA-only cluster annotation

Use this for a pre-clustered AnnData object with gene-level RNA measurements. Point
the template's root assay at the actual expression and detection representations:

```yaml
input:
  kind: h5ad
  path: data/clustered.h5ad
  assays:
    transcriptome:
      modality: root
      kind: rna
      feature_space: gene
      representations:
        logcounts:
          source: layer
          key: logcounts
          feature_space: gene
          ranking_backend: scipy_mannwhitney
          ontology_compatible: true
          detection:
            source: layer
            key: counts
            required: true
            semantics: counts > 0
```

Set `input.keys.canonical_parent` to the frozen cluster column and register local
marker programs or ontology tables under `knowledge.providers`. Then run the guided
CLI or the `all_evidence` DAG target.

### 2. CITE-seq and multimodal annotation

Use `h5mu` and declare each modality independently. RNA may support ontology markers;
protein may provide orthogonal identity or state evidence. Each evidence lane must
refer to a declared `(assay, representation)` pair. Do not copy RNA detection rules
onto antibody-derived counts, and do not discard an RNA/protein conflict: it belongs
in the critic and review artifacts.

### 3. ATAC, multiome, and spatial evidence

For ATAC or multiome data, declare peak accessibility and derived gene activity as
different assays or representations with explicit feature spaces. Only a compatible
gene-level representation may enter a gene-marker lane. For spatial data, set the
spatial key and preserve coordinates in the canonical object, but keep neighborhood
or location evidence separate from the cell-type and cell-state calls.

### 4. Evidence-only review

Set `run.mode: evidence-only` when the goal is to inspect marker coherence,
alternatives, technical confounding, or split stability without assigning durable
labels. This is useful for evaluating an existing clustering or preparing a review
packet. `Unknown`, `RETAIN PARENT`, and technical/unresolved outcomes are successful,
explicit results when the acceptance gates are not met.

### 5. Hierarchical refinement

List the hierarchy levels in `guidance.hierarchy_levels`. The normal execution path
first audits compatible stored resolutions, then performs parent-local reclustering
only for eligible scopes. The standard pipeline intentionally accepts one declared
level so it cannot finalize a partial hierarchy. For two or more levels, pass
`--score-manifest hierarchy-scores.json` to `run execute`; each entry must identify
its exact `(level, parent_scope)` and pass the same scope, evidence, and review gates.

## Outputs and approval gates

Artifacts are written beneath `run.output_root/run.run_id` with configuration,
input, rule, and output hashes. The main milestones are evidence receipts, hierarchy
decisions, critic reconciliation, the immutable review packet, final provenance, and
the write-back preview. Annotation execution never mutates the source object. Only a
separate approval bound to the preview and mapping hashes can authorize a new
annotated `.h5ad` or `.h5mu` output.

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

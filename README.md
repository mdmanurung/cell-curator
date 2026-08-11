# cell-curator

Evidence-driven cell-type annotation for **pre-clustered** AnnData and MuData — as a
Claude Code plugin, a Codex skill, a Python API, and a CLI.

It does not hand you labels. It computes the evidence a reviewer needs to assign
them: signed marker programs, bottom-up markers, ontology checks, multimodal
cross-checks, and calibrated uncertainty — then holds the run at explicit gates until
a person signs off. `Unknown` is a valid answer.

- **Local-first.** Local markers, ontology, and labeled references are enough for a
  complete offline run. Nothing contacts a third-party inference service.
- **Nothing mutates your data.** Labels reach your object only through a separate,
  hash-bound write-back approval.
- **Auditable.** Every run records config, input, rule, and output hashes.
- **No lock-in.** No dependency on another annotation package, agent platform, MCP,
  atlas, or web service.

Requires Python 3.11+ and clusters you already trust — cell-curator annotates a
clustering, it does not produce one.

## Install

**1. The Python runtime** — required in every case, including for the plugin, which
shells out to the `cell-curator` CLI for all deterministic work.

```bash
git clone https://github.com/mdmanurung/cell-curator.git
cd cell-curator
uv sync --frozen --extra workflow
uv run cell-curator --version
```

Or install the built wheel without cloning (not on PyPI — see
[releases](https://github.com/mdmanurung/cell-curator/releases/latest)):

```bash
uv pip install https://github.com/mdmanurung/cell-curator/releases/download/v0.2.0/cell_curator-0.2.0-py3-none-any.whl
```

A wheel install still ships the config template, skill, and workflow — the commands
below say `skills/…` and `workflow/…` because they assume a clone, but from a wheel
the same files live under `<prefix>/share/cell-curator/` (`SKILL.md`, `assets/`,
`references/`, `agents/`, `workflow/`, `schemas/`, `benchmarks/`).

**2. The agent skill.** For Claude Code:

```
/plugin marketplace add https://github.com/mdmanurung/cell-curator
/plugin install cell-curator@cell-curator
```

For Codex, link the skill directory into the skills path:

```bash
ln -s "$PWD/skills/cell-curator" ~/.codex/skills/cell-curator
```

<details>
<summary>pip instead of uv, and optional extras</summary>

```bash
python3.11 -m venv .venv && source .venv/bin/activate
python -m pip install ".[workflow]"
```

The core install omits optional knowledge providers and workflow execution. Add only
what a run needs: `celltypist`, `reference-models`, `census`, `scimilarity`,
`workflow`, or combinations such as `.[workflow,celltypist]`.

`gpu` is not a portability mode. A configured RAPIDS lane also requires a real
NVIDIA GPU, driver, CUDA runtime, and declared RAPIDS contract. It fails closed when
that contract is unmet and never falls back to CPU.
</details>

## Quick start — drive it with Claude

Install the plugin and the runtime, then describe the job in plain language. The
agent picks routes, explains each gate, and asks you to confirm before it labels
anything:

> Annotate `data/pancreas_clustered.h5ad` with cell-curator. Human pancreas, adult,
> healthy donors, 10x 3' dissociated. Clusters are in `obs["leiden_1"]`, log-normalized
> values in `layers["logcounts"]` and raw counts in `layers["counts"]`. Use my marker
> programs in `markers/pancreas.json`. Start evidence-only — I want to see marker
> coherence and uncertainty before we commit to labels.

Useful follow-ups once it is running:

> Which clusters came back `Unknown`, and what was missing for each?

> Cluster 7 looks impure. Audit it for donor dominance and doublets before we split it.

> Build the review packet and show me the diff you would write back.

The agent reads [`skills/cell-curator/SKILL.md`](skills/cell-curator/SKILL.md) and the
contracts in [`references/`](skills/cell-curator/references/). It cannot skip a gate
by running non-interactively — that mode requires the answers declared up front, and
records the same ledger.

## Quick start — drive it from Python

Same engine, no agent. `CellCurator` builds the strict configuration from keyword
arguments and calls the same functions the CLI and Snakemake DAG call:

```python
from cell_curator import CellCurator

cur = CellCurator(
    adata=adata,                      # or path="pancreas.h5ad"
    organism="human",
    tissue="pancreas",
    cluster_key="leiden_1",
    markers="markers/pancreas.json",  # or a {label: {positive, negative}} mapping
    run_id="pancreas-v1",
    output_root="/abs/path/results",
    mode="annotation",                # or "evidence-only"
    # Context the package refuses to guess. Omit any and prepare() stops.
    developmental_stage="adult",
    condition="healthy",
    experimental_context="pancreas atlas v1",
    # Gates are declared, not skipped. Drop these two to confirm interactively.
    preauthorized=True,
    reviewer="your-name",
    assumptions=["The signed marker programs define the L1 vocabulary."],
)

cur.inspect()                          # structure and keys; changes nothing
cur.prepare()                          # freeze + hash the input, then scene and plan
evidence, bottom_up = cur.evidence()   # the tables you reason over
packet = cur.annotate()                # remaining phases; freezes a review packet
cur.validate_run()
```

`annotate()` is idempotent — calling it after the granular steps resumes rather than
repeats — and writes nothing into your object. That is `write_back()`, and it requires
an approval file.

A runnable version on synthetic data is in
[`examples/quickstart_synthetic.ipynb`](examples/quickstart_synthetic.ipynb). The test
suite executes it, so it cannot drift.

## Configure for your data

Every route starts from one strict config. Copy the template and validate it before
loading anything:

```bash
cp skills/cell-curator/assets/config.template.yaml config.yaml
uv run cell-curator config validate --config config.yaml
```

`input.assays` names the modalities and representations present; `evidence.lanes`
selects which of them may support a decision. Encode only evidence the assay actually
measures:

| Your data | Start from | Boundary that matters |
| --- | --- | --- |
| Clustered scRNA-seq (`.h5ad`) | `input.kind: h5ad`, one root `rna` assay; log-normalized values for ranking, counts for detection | Counts and log values are different evidence; a missing gene is unmeasured, not negative |
| CITE-seq / multimodal (`.h5mu`) | `input.kind: h5mu`, RNA and protein as separate assays and lanes | Preserve cross-modal disagreement for review instead of averaging it away |
| scATAC-seq or multiome | `atac` peak accessibility declared separately from any `gene_activity` representation | Peaks are direct accessibility evidence; gene activity is indirect and never silently becomes RNA |
| Spatial transcriptomics | The transcriptome assay plus `input.keys.spatial` | Proximity is context, not identity |
| Any assay, no labels yet | `run.mode: evidence-only` with frozen clusters | Auditable evidence and uncertainty without forcing a biological label |
| Hierarchical annotation | `guidance.hierarchy_levels`; enable parent-local refinement only where gates pass | Every split is scoped to its parent and may validly retain it or stay unresolved |

Minimal RNA assay block — point it at your actual layers:

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

Set `input.keys.canonical_parent` to the frozen cluster column, register marker or
ontology tables under `knowledge.providers`, and keep `run.immutable: true`. Use a new
`run.run_id` whenever the input or configuration changes, and set
`input.expected_sha256` when the input arrives from elsewhere.

Per-assay rules — which evidence each modality may support, and what stays blocked —
are in [assay-routing.md](skills/cell-curator/references/assay-routing.md).

**Two or more hierarchy levels.** The standard pipeline accepts one declared level so
it cannot finalize a partial hierarchy. For deeper vocabularies, pass
`--score-manifest hierarchy-scores.json` to `run execute`; each entry names its exact
`(level, parent_scope)` and clears the same scope, evidence, and review gates.

## What a run produces

Artifacts land under `run.output_root/run.run_id`, hashed at each step: evidence
receipts, hierarchy decisions, critic reconciliation, an immutable review packet,
final provenance, and a write-back preview.

Annotation never touches the source object. A new annotated `.h5ad`/`.h5mu` requires a
separate approval bound to the preview and mapping hashes:

```bash
uv run cell-curator write-back preview --config config.yaml --output annotated.h5ad
uv run cell-curator write-back apply --config config.yaml \
  --approval approval.json --output annotated.h5ad
```

`approval.json` must bind every field `write-back preview` emits — run, source,
mapping, durable labels, assumptions, report, critics, diff, exact output path, and
in-place intent. Schemas are in
[artifact-schemas.md](skills/cell-curator/references/artifact-schemas.md).

## Command line and Snakemake

The guided CLI is the safest unattended route. It stops before the first biological
label until context and assumptions are confirmed:

```bash
uv run cell-curator config validate      --config config.yaml
uv run cell-curator environment detect   --config config.yaml
uv run cell-curator data inspect         --config config.yaml
uv run cell-curator data canonicalize    --config config.yaml
uv run cell-curator data qc              --config config.yaml
uv run cell-curator run scene            --config config.yaml
uv run cell-curator run confirm-context  --config config.yaml \
  --input context-confirmation.json --reviewer REVIEWER
uv run cell-curator run plan             --config config.yaml
uv run cell-curator assumptions approve  --config config.yaml --reviewer REVIEWER
uv run cell-curator run execute          --config config.yaml
```

`cell-curator --help` lists all 21 command groups.

For a resumable DAG, Snakemake owns dependencies and staleness. Its default aggregate
target stops at evidence; the later stages stay explicit:

```bash
uv run snakemake --snakefile workflow/Snakefile all_evidence \
  --configfile config.yaml --cores 4
uv run snakemake --snakefile workflow/Snakefile review_packet \
  --configfile config.yaml --cores 4
uv run snakemake --snakefile workflow/Snakefile finalization \
  --configfile config.yaml --cores 4
```

`--dry-run` first on a new configuration. The `writeback` target stays blocked until
its run-bound approval exists. For HPC, adapt
[`workflow/profiles/slurm/config.yaml`](workflow/profiles/slurm/config.yaml) rather than
hard-coding cluster resources into the workflow. RAPIDS lanes request one `gpu`
resource and validate the real device at runtime.

## Python API

Beyond the quick start, the object exposes each phase. Refinement is the part worth
knowing — cluster impurity is investigated, not assumed:

```python
cur.audit_parents()        # per-parent impurity: donor/capture dominance,
                           # doublet and incompatible-program fractions, QC outliers
cur.propose_subclusters()  # bounded candidates: reuse a matching stored resolution,
                           # else one parent-local pass
cur.decide_subclusters()   # ACCEPT SPLIT | RETAIN PARENT | DOUBLET/MIXED |
                           # TECHNICAL/UNRESOLVED, from frozen evidence
cur.refine_clusters()      # all four, in order
```

Canonical parents stay frozen throughout. A parent whose impurity turns out to be
technical or doublet-driven is reported as such rather than carved into subtypes.

Three things worth knowing in a notebook: an in-memory `AnnData` is written to the run
directory and hashed before anything reads it; relative paths anchor to the config
file's directory rather than the process working directory, so a `chdir` or kernel
restart will not break a run; and anything without a keyword argument can be set
through `overrides=`, deep-merged onto the template last.

Long phases are silent by default:

```python
import logging
logging.getLogger("cell_curator").setLevel(logging.INFO)
logging.getLogger("cell_curator").addHandler(logging.StreamHandler())
```

<details>
<summary>Stable function API — the same implementation, without the object</summary>

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
</details>

## Scientific boundaries

- Cell type and cell state are independent axes.
- Missing panel features are unmeasured, not negative.
- ATAC peaks are accessibility evidence; gene activity is indirect evidence.
- A configured RAPIDS lane requires its declared real GPU/software contract and fails
  closed. It never substitutes a CPU result.
- `Unknown` is a valid result when coverage, coherence, margin, or required capability
  is inadequate.
- Labels may remain CL-free when no ontology is available; supplied CL identifiers
  always require local semantic validation.
- No superiority claim is made without measured neutral benchmark results. Public
  benchmark manifests record pending external inputs and exact reproduction commands.

## Documentation

| Topic | Document |
| --- | --- |
| Agent workflow, stage by stage | [SKILL.md](skills/cell-curator/SKILL.md) |
| Input, IDs, matrices, QC | [input-contract.md](skills/cell-curator/references/input-contract.md) |
| Assay and provider routing | [assay-routing.md](skills/cell-curator/references/assay-routing.md) |
| Parent-local refinement | [adaptive-subclustering.md](skills/cell-curator/references/adaptive-subclustering.md) |
| Signed and multimodal evidence | [multimodal-evidence.md](skills/cell-curator/references/multimodal-evidence.md) |
| Assumptions, decisions, critics | [decision-gates.md](skills/cell-curator/references/decision-gates.md) |
| Durable artifact schemas | [artifact-schemas.md](skills/cell-curator/references/artifact-schemas.md) |
| Compute, GPU contract, resume | [compute-and-snakemake.md](skills/cell-curator/references/compute-and-snakemake.md) |
| Strict completion criteria | [validation.md](skills/cell-curator/references/validation.md) |
| Neutral capability mapping | [quadbio-parity.md](skills/cell-curator/references/quadbio-parity.md) |

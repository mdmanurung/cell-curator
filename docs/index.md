```{toctree}
:hidden:
:caption: Getting started

Install and quick starts <_generated/readme>
Agent workflow <_generated/skill>
```

```{toctree}
:hidden:
:caption: Python API

_generated/api
```

```{toctree}
:hidden:
:caption: Contracts

_generated/reference/input-contract
_generated/reference/assay-routing
_generated/reference/adaptive-subclustering
_generated/reference/multimodal-evidence
_generated/reference/decision-gates
_generated/reference/artifact-schemas
_generated/reference/compute-and-snakemake
_generated/reference/validation
_generated/reference/quadbio-parity
```

# cell-curator

Evidence-driven cell-type annotation for **pre-clustered** AnnData and MuData — as a
Claude Code plugin, a Codex skill, a Python API, and a CLI.

Start with clusters you already trust. cell-curator assembles the evidence for each
label, preserves inadequate evidence as `Unknown`, and pauses for review before
anything is written back.

## Choose how to start

::::{grid}
:gutter: 2

:::{grid-item-card} Start with the plugin
:columns: 12 12 7 7
:link: _generated/readme
:link-type: doc

Install the runtime and agent skill, then adapt the example prompt with your input
path, cluster key, assay layers, biological context, and marker programs.
:::

:::{grid-item-card} Run from Python
:columns: 12 12 5 5
:link: _generated/api
:link-type: doc

Use `CellCurator` in a notebook or script, or call the stable functions shared by
the CLI and Snakemake workflow.
:::

:::{grid-item-card} Configure your assay
:columns: 12 12 5 5
:link: _generated/reference/assay-routing
:link-type: doc

Declare RNA, protein, ATAC or gene activity, and spatial evidence without crossing
modality or representation boundaries.
:::

:::{grid-item-card} Understand the gates
:columns: 12 12 7 7
:link: _generated/reference/decision-gates
:link-type: doc

See when a run must stop, how `Unknown` is retained, and what must be approved before
labels can reach the source object.
:::

::::

## The rules it will not bend

- Cell type and cell state are independent axes.
- Missing panel features are unmeasured, not negative.
- A configured RAPIDS lane fails closed rather than substituting a CPU result.
- `Unknown` is a valid result when coverage, coherence, margin, or capability is
  inadequate.
- Nothing is written back to your object without a separate, hash-bound approval.

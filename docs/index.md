```{toctree}
:hidden:
:caption: Getting started

Overview <_generated/readme>
```

```{toctree}
:hidden:
:caption: Python API

_generated/api
```

```{toctree}
:hidden:
:caption: Contracts

_generated/skill
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

It does not hand you labels. It computes the evidence a reviewer needs to assign
them, then holds the run at explicit gates until a person signs off. `Unknown` is a
valid answer.

## Where to start

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} Overview
:link: _generated/readme
:link-type: doc

Install, the two quick starts, configuration for your assay, and what a run
produces. Start here.
:::

:::{grid-item-card} Python API
:link: _generated/api
:link-type: doc

`CellCurator` and the stable function API, generated from the source.
:::

:::{grid-item-card} Agent workflow
:link: _generated/skill
:link-type: doc

The stage-by-stage instructions Claude and Codex follow when driving the skill.
:::

:::{grid-item-card} Contracts
:link: _generated/reference/input-contract
:link-type: doc

Input, assay routing, refinement, evidence, gates, artifact schemas, compute, and
strict completion criteria.
:::

::::

## The rules it will not bend

- Cell type and cell state are independent axes.
- Missing panel features are unmeasured, not negative.
- A configured RAPIDS lane fails closed rather than substituting a CPU result.
- `Unknown` is a valid result when coverage, coherence, margin, or capability is
  inadequate.
- Nothing is written back to your object without a separate, hash-bound approval.

```{toctree}
:hidden:
:caption: Getting started

Install and examples <_generated/readme>
Agent workflow <_generated/skill>
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Python API

_generated/api
```

```{toctree}
:hidden:
:caption: Methods and validation

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

This package annotates pre-clustered AnnData and MuData from marker, ontology,
reference, and multimodal evidence. It records uncertainty and requires review before
write-back.

## Start by task

::::{grid}
:gutter: 2

:::{grid-item-card} Use an agent
:columns: 12 12 7 7
:link: _generated/readme
:link-type: doc

Install the runtime and skill, then provide the input, assay, cluster key, context,
and marker programs.
:::

:::{grid-item-card} Use Python
:columns: 12 12 5 5
:link: _generated/api
:link-type: doc

Run `CellCurator` from a notebook or script, or call the underlying functions.
:::

:::{grid-item-card} Declare assays
:columns: 12 12 5 5
:link: _generated/reference/assay-routing
:link-type: doc

Define RNA, protein, ATAC, gene-activity, and spatial evidence without mixing their
measurement semantics.
:::

:::{grid-item-card} Review decisions
:columns: 12 12 7 7
:link: _generated/reference/decision-gates
:link-type: doc

Review uncertainty, retained parent labels, and the approvals required for
write-back.
:::

::::

## Analytical constraints

- Cell type and cell state are separate axes.
- Missing features are unmeasured, not negative evidence.
- A required RAPIDS result is produced only under its declared GPU contract.
- `Unknown` is retained when the evidence is inadequate.
- Write-back requires a person to type the token of the exact reviewed state.

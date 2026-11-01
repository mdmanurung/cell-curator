# Installation

cell-curator needs Python 3.11 or newer and a clustering you already made.

## Install

```bash
pip install cell-curator            # from PyPI, from 0.3.0 on
uv tool install "cell-curator[workflow]"   # the CLI as a standalone tool
```

To work from a checkout:

```bash
git clone https://github.com/mdmanurung/cell-curator.git
cd cell-curator
uv sync --frozen --extra workflow
uv run cell-curator --version
```

## Extras

| Extra | Adds | Also needs |
| --- | --- | --- |
| `workflow` | Snakemake and its Slurm executor | a Slurm partition and account in `workflow/profiles/slurm/config.yaml` |
| `celltypist` | the CellTypist cross-check | a model already in CellTypist's cache; cell-curator never downloads one |
| `reference-models` | scVI, scANVI, and scArches | a CUDA GPU; there is no CPU path |
| `census` | CELLxGENE Census reference slices | `knowledge.network: online` |
| `scimilarity` | the SCimilarity cross-check | a local model directory |
| `gpu` | RAPIDS marker ranking | an NVIDIA GPU that meets the RAPIDS contract |

None is needed for the core marker and ontology workflow.

## What the wheel installs besides the package

The wheel puts the agent skill, config template, marker panels, workflow, schemas,
and benchmark manifests under `<prefix>/share/cell-curator/`. Find the prefix with:

```bash
python -c "import sysconfig; print(sysconfig.get_path('data'))"
```

The shipped human PBMC marker panel is at
`<prefix>/share/cell-curator/assets/markers/human-pbmc.json`.

## Check the install

```bash
cell-curator --version
cell-curator config validate --config <prefix>/share/cell-curator/assets/config.template.yaml
```

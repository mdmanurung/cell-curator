# Compute, Snakemake, and immutability

## Workflow boundaries

`workflow/Snakefile` defines the manifest-driven DAG. Its phases are:

1. freeze the input and audit parent clusters on CPU;
2. discover candidates and assess local graph stability on CPU;
3. prepare binary comparisons on CPU;
4. run backend-specific marker methods, including validated RAPIDS routes;
5. compute enrichment, absolute evidence, technical summaries, and critic inputs;
6. build the review packet;
7. validate the final mapping;
8. apply an approved write-back.

Review, finalization, and write-back are explicit targets. Supply the configuration
with `--configfile` and use `workflow/profiles/slurm` for scheduler execution.
Snakemake manages dependencies, stale artifacts, and resume.

## Immutable artifacts

Each rule writes to a partial file on the output filesystem, validates it, and then
publishes it atomically. An existing output may be reused only when its manifest
matches the requested run. Configuration changes require a new run ID.

Manifests record input, configuration, rule, script, registry, barcode, and output
hashes; package versions; commands; parameters; runtime; and fallbacks. RAPIDS
receipts also record the device, compute capability, visible-device count, driver,
CUDA runtime, RAPIDS distribution, native API, resolved options, and per-method group
coverage.

## GPU contract

Every RAPIDS rule validates its GPU contract before opening the input. A failed
check creates no output or partial evidence and leaves the downstream target
incomplete. Other marker methods use separate rules; they do not replace a required
RAPIDS result. Declare GPU models and scheduler resources in the configuration.

Validate the runtime inside the allocated job. Pin CPU and GPU environments with
lockfiles or immutable definitions, and record whether each environment was solved
or reused. Python 3.11 or newer is required.

## Validated RAPIDS environment

The GPU route was validated on an NVIDIA RTX PRO 6000 Blackwell Server Edition MIG
`1g.24gb` slice with compute capability 12.0, driver 580.159.04, and CUDA 13.2. The
validation covered contract rejection, Snakemake GPU resource requests, and native
Wilcoxon and logistic-regression execution through
`rapids_singlecell.tl.rank_genes_groups`.

### Release selection

`tl.rank_genes_groups` first appeared in rapids-singlecell 0.14.

| Release | Distribution | Python | Constraint |
| --- | --- | --- | --- |
| 0.13.x or earlier | wheel | 3.11 or newer | No native Wilcoxon; unsuitable for this route. |
| 0.14.0 to 0.14.1 | wheel | 3.12 or newer | Native Wilcoxon; no `use_continuity` or `multi_gpu`. |
| 0.15 or later | source | 3.12 or newer | Full signature; build requires a complete CUDA toolkit. |

A wheel-only environment can be installed without `nvcc`:

```bash
uv venv --python 3.12 .venv-gpu
uv pip install --python .venv-gpu \
  --extra-index-url https://pypi.nvidia.com --index-strategy unsafe-best-match \
  "rapids-singlecell[rapids13]==0.14.1" \
  "cudf-cu13==25.12.*" "cuml-cu13==25.12.*" "cuvs-cu13==25.12.*" "cugraph-cu13==25.12.*"
```

Pin RAPIDS to 25.12. RAPIDS 26.x removed
`cuml.thirdparty_adapters.check_array`, which rapids-singlecell 0.13 and 0.14 use.

### Options and group coverage

Optional ranking parameters are checked against the installed function signature
before the input is opened. A requested option is applied only when supported; the
resolution is recorded under `runtime.ranking_options`.

Wilcoxon reports every comparison group. Binary logistic regression reports one
coefficient vector, so only one group is emitted. Receipts record this distinction
under `runtime.ranking_group_coverage`; no result is inferred for an unreported
group. On the validated two-cluster run, coverage was `wilcoxon: ["0", "1"]` and
`logreg: ["0"]`.

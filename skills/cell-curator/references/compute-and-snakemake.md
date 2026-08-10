# Compute, Snakemake, and immutability

## Phase boundaries

Use explicit aggregate targets:

1. input freeze and parent audits (CPU);
2. candidate discovery and local graphs/stability (CPU, with local resources);
3. binary comparison preparation (CPU);
4. backend-specific marker lanes (validated GPU for RAPIDS lanes);
5. applicable enrichment, absolute evidence, technical summaries, and critic inputs;
6. explicit review-packet build;
7. explicit final mapping validation;
8. separate, human-approved write-back.

Review packets and finalizers are never default targets. Snakemake owns
dependencies and resume; do not rerun valid immutable artifacts manually.

## Immutable writer contract

Each rule writes to a versioned path through a same-filesystem partial file,
validates postconditions, then atomically publishes. Existing outputs cause a
failure unless their manifest proves an exact resumable match. Never delete a
valid artifact to satisfy a rerun. A new configuration produces a new run ID.

## Manifests

Record input/config/rule/script hashes, package versions, command/parameters,
candidate registry hash, ordered and set barcode hashes, output hashes, runtime,
and fallbacks. RAPIDS receipts additionally record device model, compute
capability, visible-device count, driver, CUDA runtime, RAPIDS distribution,
and native API.

## GPU fail-closed behavior

GPU validation is the first action in every RAPIDS marker rule. If it fails:

- do not read the input matrix;
- do not create an output directory;
- do not create partial evidence;
- do not schedule a CPU replacement;
- leave the downstream completion target unsatisfied.

Non-RAPIDS methods use separate rules and receipts; they never serve as a CPU
fallback for a required RAPIDS lane. GPU models and scheduler resources are
configuration inputs. Do not hard-code a partition, model, account, or wall time.

## Environment handling

Pin environments with a lockfile or immutable environment definition. Validate
the runtime inside the allocated job, not only on the submission host. Record
whether the environment was solved or reused. Keep CPU and GPU phases separate
when their dependency stacks differ. The package requires Python 3.11 or newer
plus its locked scientific dependencies; do not run it with a
legacy login-shell Python by accident.

## Workflow

`workflow/Snakefile` is the executable manifest-driven DAG. Supply the strict
configuration with `--configfile`, retain its evidence/review/finalization/write-back
boundaries, and use `workflow/profiles/slurm` for configured scheduler execution.

## Real-hardware GPU validation status

Validated on an NVIDIA RTX PRO 6000 Blackwell Server Edition MIG `1g.24gb` slice
(compute capability 12.0, driver 580.159.04, CUDA runtime 13.2).

Confirmed working on that hardware:

- the declared/observed contract check in `rank_markers_gpu.inspect_runtime`
  validates version, distribution, CUDA major, single-device visibility, and the
  device/compute-capability allowlist;
- a deliberately wrong `allowed_gpu_contracts` entry is rejected on the live
  device rather than merely parsed;
- the Snakemake DAG requests a GPU for a RAPIDS lane (`evidence_gpu` returns 1);
- **native GPU Wilcoxon and logreg both execute** through
  `rapids_singlecell.tl.rank_genes_groups`.

### Choosing a rapids-singlecell release

`tl.rank_genes_groups` — the single entry point for both `method="wilcoxon"` and
`method="logreg"` — arrived in **0.14**. Releases before it expose only
`rank_genes_groups_logreg` and contain no Wilcoxon at all, so the lane cannot
satisfy its receipt contract on them; `native_ranking_callable` fails closed with
the installed version named.

Version selection is a real constraint, not a formality:

| Release | Ships | Python | Notes |
| --- | --- | --- | --- |
| ≤ 0.13.x | wheel | ≥ 3.11 | no `rank_genes_groups`, no Wilcoxon — unusable for this lane |
| 0.14.0–0.14.1 | wheel | ≥ 3.12 | has Wilcoxon; **no** `use_continuity` or `multi_gpu` |
| ≥ 0.15 | sdist only | ≥ 3.12 | full signature, but building it needs a complete CUDA toolkit (`nvcc`, `cuda_runtime.h`, `crt/host_config.h`); the fragmented pip CUDA wheels do not satisfy CMake's toolkit detection |

Working wheel-only environment, no `nvcc` required:

```bash
uv venv --python 3.12 .venv-gpu
uv pip install --python .venv-gpu \
  --extra-index-url https://pypi.nvidia.com --index-strategy unsafe-best-match \
  "rapids-singlecell[rapids13]==0.14.1" \
  "cudf-cu13==25.12.*" "cuml-cu13==25.12.*" "cuvs-cu13==25.12.*" "cugraph-cu13==25.12.*"
```

Pin RAPIDS to 25.12: `rapids-singlecell` 0.13/0.14 are incompatible with RAPIDS
26.x, which removed `cuml.thirdparty_adapters.check_array`.

### Optional statistical switches are resolved against the real signature

`rank_genes_groups` gained parameters over time and also accepts `**kwds`, which
would swallow an unknown name without applying it. So each optional switch is
checked against the installed signature: applied when declared, accepted as a
no-op only when it was requested *off* (absent means not applied, which is the
same result), and refused outright when requested *on* against a release that
cannot honour it. The resolution is recorded in the receipt under
`runtime.ranking_options`, and it runs before the input is opened so an
unhonourable option fails fast.

### Still open

Result parsing does not yet match 0.14.1's `uns` layout. For a two-group
one-vs-rest comparison, `uns[key]["names"]` came back as a recarray with a single
field `'0'` instead of one field per group, so reading group `'1'` raises
`ValueError: no field of name 1`. The GPU-gated test carries this as a strict
xfail so it flips to a failure once fixed. Diagnosing it needs a GPU debug run
over the returned structure — and the ≥0.15 layout may differ again.

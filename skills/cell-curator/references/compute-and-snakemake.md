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

Validated on 2026-08-09 on an NVIDIA RTX PRO 6000 Blackwell Server Edition
MIG `1g.24gb` slice (compute capability 12.0, driver 580.159.04, CUDA runtime
13.2) with `rapids-singlecell` 0.13.4 and RAPIDS `cu13` 25.12.

Confirmed working on that hardware:

- the declared/observed contract check in `rank_markers_gpu.inspect_runtime`
  validates version, distribution, CUDA major, single-device visibility, and
  the device/compute-capability allowlist;
- a deliberately wrong `allowed_gpu_contracts` entry is rejected on the live
  device rather than merely parsed;
- the Snakemake DAG requests a GPU for a RAPIDS lane (`evidence_gpu` returns 1).

Confirmed blocked, and the reason the lane still cannot execute:

- `rapids-singlecell` exposes no `tl.rank_genes_groups` and contains no
  Wilcoxon implementation at all — `tl.rank_genes_groups_logreg` is the only
  ranking entry point. `rank_markers_gpu.inspect_runtime` and
  `evidence.compute_gpu_bottom_up_markers` bind `rsc.tl.rank_genes_groups`
  with a scanpy-style signature, so both raise `AttributeError` on any real
  GPU. The `BackendReceipt` contract additionally requires both native
  `wilcoxon` and `logreg` marker methods, which that library cannot supply.

Resolving this is a scientific decision, not a rename: either supply a GPU
Wilcoxon (for example a CuPy implementation carrying its own tie and
continuity correction) or restate the RAPIDS receipt contract around the
methods the library actually provides. Until then the lane fails closed, which
is the intended safety behaviour but is not the same as a working GPU path.

Environment note for reproduction: `rapids-singlecell` 0.13.4 is the newest
release supporting Python 3.11, and it is incompatible with RAPIDS 26.x
(`cuml.thirdparty_adapters.check_array` was removed), so pin RAPIDS to 25.12.
Newer `rapids-singlecell` needs Python 3.12+, where `scikit-misc` has no wheel
and its CMake build fails on this cluster.

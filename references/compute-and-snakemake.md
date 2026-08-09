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

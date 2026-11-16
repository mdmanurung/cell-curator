# Container for annotation-consensus

An Apptainer image that holds the environment only. The skill code stays outside and is bound read-only, so one image runs any commit of the skill. Pinned: the base image digest, 78 hashed packages and one `cell_curator` wheel built from a git commit.

## Build and run

```bash
cd skills/annotation-consensus/container
bash build_wheel.sh                                  # wheelhouse/cell_curator-*.whl, cell-curator.txt
sbatch -o _build/sbatch_%j.log build_sif.sbatch      # _build/consensus.sif and consensus.record.txt (about 6 min)
bash run_container.sh _build/consensus.sif OUT        # the demo; OUT must not exist
bash run_container.sh _build/consensus.sif OUT verify --demo
```

`run_container.sh` uses `--cleanenv --containall`, binds the repo read-only, and passes the SIF sha256 in `CONSENSUS_SIF_SHA256`. The skill writes it to `run_manifest.json` (`container_sif_sha256`), which is not hashed.

## What is pinned

| Item | Value |
| --- | --- |
| Base image | `python@sha256:dbf1de47…9d0c` (python:3.11.12-slim-bookworm, resolved 2026-11-12) |
| Packages | `lock/requirements.lock`: hashes for every line, installed with `--require-hashes --only-binary=:all: --no-deps`. Versions equal the host `.venv` (`lock/host_freeze.txt`) |
| `cell_curator` | wheel from `build_wheel.sh`; the hash goes to `cell-curator.txt`. `SOURCE_DATE_EPOCH` is the time of the last commit that changed `src`, `pyproject.toml`, `README.md`, `LICENSE`, `schemas` or `workflow`, so the hash does not move for other commits |
| `connection-pool 0.0.3` | no wheel on PyPI. Wheel built from the checked sdist (`wheelhouse/SDIST_…src`) with a fixed epoch; the lock line holds the wheel hash. The setuptools backend for that one build is not pinned |
| BLAS kernel | `OPENBLAS_CORETYPE=SANDYBRIDGE`, the lowest common kernel of the cluster (one node type has AVX but no AVX2). Threads are 1 |

Rebuild the wheel and the SIF after any change under `src/`, because the image carries the package.

## Measured on 2026-11-16 (demo, 654 cells)

SIF `2deafd07f26c…`, built from `fd819fb` in 5 min 33 s. The table compares each run's 18 hashed files with the committed `expected_hashes.json`.

| Run | Where | Files differing |
| --- | --- | --- |
| login | `res-hpc-exe102` | 0 of 18 |
| login, no network (`unshare -rn`) | `res-hpc-exe102` | 0 of 18 |
| short partition | `res-hpc-path02`, Xeon Gold 6234 | 0 of 18 |
| highmem partition | `res-hpc-exe032`, Xeon E5-4650 v3 | 0 of 18 |

`verify --demo` inside the container reports "18 stored files and the rerun are byte-identical". The container passed the version check with no `--allow-version-drift`.

`table fallback: the evidence is sha256 equality, not a plotted quantity`.

## Limits

- Two builds of the same definition give different SIF hashes (the squashfs stores build times), with the same file contents. Record the SIF hash per build; `run_manifest.json` does.
- The demo is small. Equal hashes across these CPUs do not show the same for AIFI-scale input.
- Without `OPENBLAS_CORETYPE` the demo hashes were not measured. numpy's own SIMD dispatch is not pinned.
- pixi was not needed and was not built.
- `--score` and the AIDA and AIFI regressions were not run in the container.

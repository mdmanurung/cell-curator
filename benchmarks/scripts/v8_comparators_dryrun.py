"""Dry-run check for the v8 comparator scripts: paths and argument parsing only.

Purpose: catch a wrong path, an existing output dir or a bad argument before a job is
         queued. It loads no data, trains nothing and submits nothing.
Run:     python benchmarks/scripts/v8_comparators_dryrun.py \
             SAMPLED_TSV TRAIN_DIR QUERY_H5AD OUT_ROOT
         TRAIN_DIR may be missing before v8_prepare_training has run; the comparator
         checks are then reported as pending, not failed.
Inputs:
    SAMPLED_TSV, TRAIN_DIR, QUERY_H5AD   as the sbatch files take them
    OUT_ROOT   parent of the output paths; the names used below must not exist in it
Outputs:
    printed    one line per check: ok, FAIL or pending; exit status = number of FAILs
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENVS = Path("/exports/para-lipg-hpc/mdmanurung/envs")
PY = ENVS / "cc-celltypist/bin/python"
R_BIN = Path("/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript")
R_LIB = Path("/exports/para-lipg-hpc/mdmanurung/R/4.5/SingleR")
ATLAS = Path(
    "/exports/para-lipg-hpc/mdmanurung/cc-v8-work/acquire-2026-11-09/downloads/"
    "e522d2cd-7927-4e59-a4ed-064009569279.h5ad"
)
SBATCH = {  # file: runs on the GPU partition
    "v8_prepare_training.sbatch": False,
    "v8_celltypist.sbatch": False,
    "v8_singler.sbatch": False,
    "v8_scanvi.sbatch": True,
    "v8_popv.sbatch": True,
}


def report(ok, text):
    print(f"{'ok' if ok else 'FAIL':8s} {text}")
    return 0 if ok else 1


def main():
    sampled, train, query, out_root = (Path(a) for a in sys.argv[1:5])
    built = (train / "train.h5ad").exists()
    # label, interpreter, script, arguments, needs a built TRAIN_DIR
    runs = [
        (
            "prepare_training",
            PY,
            "v8_prepare_training.py",
            [sampled, out_root / "train"],
            False,
        ),
        ("prepare_query", PY, "v8_prepare_query.py", [query, out_root / "query.h5ad"], False),
        ("export_for_r", PY, "v8_export_for_r.py", [query, out_root / "export"], False),
        ("celltypist", PY, "v8_celltypist.py", [train, query, out_root / "ct"], True),
        # the R script takes export folders; the job makes them, so the check borrows dirs
        ("singler", R_BIN, "v8_singler.R", [train, train, out_root / "sr.tsv"], True),
    ]
    failures = 0
    for label, interpreter, script, args, needs_train in runs:
        if needs_train and not built:
            print(f"pending  {label} --check: {train}/train.h5ad is not built yet")
            continue
        command = [str(interpreter), str(HERE / script), *map(str, args), "--check"]
        result = subprocess.run(command, capture_output=True, text=True)
        failures += report(
            result.returncode == 0,
            f"{label} --check: {(result.stdout + result.stderr).strip()[-160:]}",
        )
    for name, gpu in SBATCH.items():
        syntax = subprocess.run(["bash", "-n", str(HERE / name)], capture_output=True)
        text = (HERE / name).read_text()
        on_gpu = "--partition=gpu" in text and "--gres=gpu:1g.24gb:1" in text
        failures += report(
            syntax.returncode == 0 and on_gpu == gpu, f"{name}: bash -n, gpu={gpu}"
        )
    for path in (ATLAS, ENVS / "cc-gpu/bin/cell-curator", ENVS / "cc-popv/bin/python"):
        failures += report(path.exists(), f"exists {path}")
    for path in (R_BIN, R_LIB, HERE / "run_popv.py", HERE / "comparator_to_neutral.py"):
        failures += report(path.exists(), f"exists {path}")
    sys.exit(failures)


if __name__ == "__main__":
    main()

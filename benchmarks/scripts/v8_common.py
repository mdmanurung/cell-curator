"""Helpers shared by the protocol v8 comparator scripts.

Purpose: one argument and path check for every v8 comparator script, and the one
         conversion (log1p CP10k back to counts) that scANVI and popV need on a query.
Run:     imported by v8_*.py; not run directly.
Inputs:  none.
Outputs: none.
"""

import argparse
from pathlib import Path


def parse_args(description, inputs, outputs):
    """Positional paths: `inputs` must exist, `outputs` must not. `--check` stops there."""
    parser = argparse.ArgumentParser(description=description)
    for name in (*inputs, *outputs):
        parser.add_argument(name, type=Path)
    parser.add_argument("--check", action="store_true", help="validate paths, then exit")
    args = parser.parse_args()
    missing = [n for n in inputs if not getattr(args, n).exists()]
    existing = [n for n in outputs if getattr(args, n).exists()]
    if missing or existing:
        raise SystemExit(f"missing inputs {missing}; outputs that already exist {existing}")
    if args.check:
        print("check ok:", {n: str(getattr(args, n)) for n in (*inputs, *outputs)})
        raise SystemExit(0)
    return args


def recover_counts(lognorm):
    """Counts from a log1p CP10k matrix.

    Each cell's smallest non-zero value is one count scaled by 1e4 / total, so dividing
    the cell by it gives integers. The check fails if that does not hold.
    """
    import numpy as np
    from scipy import sparse

    expm1 = sparse.csr_matrix(lognorm, dtype=np.float64)
    expm1.data = np.expm1(expm1.data)
    filled = np.diff(expm1.indptr) > 0
    smallest = np.minimum.reduceat(expm1.data, expm1.indptr[:-1][filled])
    scale = np.ones(expm1.shape[0])
    scale[filled] = smallest
    counts = sparse.diags(1.0 / scale) @ expm1
    counts = sparse.csr_matrix(counts)
    if np.abs(counts.data - np.round(counts.data)).max() > 0.1:
        raise ValueError("X is not log1p CP10k: recovered counts are not integers")
    counts.data = np.round(counts.data).astype(np.float32)
    return counts

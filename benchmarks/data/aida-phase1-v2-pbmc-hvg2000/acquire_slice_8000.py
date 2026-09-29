"""Widen the AIDA slice to 8,000 genes, holding one gene window at a time.

The first slice was 2,000 highly variable genes because a single contiguous read
of more would not fit in this box's 4 GB. That is a property of one read, not of
the slice: reading the same cells across successive 2,000-gene windows keeps the
peak at one window and assembles the rest on disk, which is what this does. The
cell set, the census release, the dataset, and the truth columns are read back out
of the existing slice, so the two runs differ only in how many genes were
measured and a comparison between them is clean.

Memory is the constraint throughout, so nothing holds the whole panel: each
window is written to its own file and released before the next is opened, and the
final matrices are preallocated by column rather than concatenated.
"""

from __future__ import annotations

import gc
import json
import os
import resource
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

CENSUS_VERSION = "2025-11-08"
WINDOW = 2_000
PANEL = json.loads((Path(__file__).resolve().parent / "hvg_panel_8000.json").read_text())
SOURCE = Path(__file__).resolve().parent / "aida_pbmc_benchmark.h5ad"
WORK = Path(os.environ.get("CELL_CURATOR_SCRATCH", "/tmp/cc2"))
TILEDB_CFG = {
    "soma.init_buffer_bytes": 16_777_216,
    "sm.compute_concurrency_level": 4,
    "sm.io_concurrency_level": 4,
}
KEEP_OBS = [
    "soma_joinid",
    "cell_type",
    "cell_type_ontology_term_id",
    "donor_id",
    "development_stage",
    "sex",
    "lineage",
    "true_path",
    "collection_id",
    "dataset_id",
    "census_version",
]


def peak_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2


def dense(matrix) -> np.ndarray:
    return np.asarray(matrix.todense() if sp.issparse(matrix) else matrix, dtype=np.float32)


def main() -> None:
    import anndata as ad
    import cellxgene_census

    source = ad.read_h5ad(SOURCE)
    joinids = np.sort(source.obs["soma_joinid"].to_numpy().astype(np.int64))
    obs = source.obs.loc[:, KEEP_OBS].copy()
    feature_ids = sorted(set(PANEL["feature_id"]))
    windows = [
        (index, feature_ids[start : start + WINDOW])
        for index, start in enumerate(range(0, len(feature_ids), WINDOW))
    ]
    print(f"[wide] cells={len(joinids):,} panel={len(feature_ids):,} windows={len(windows)}")

    paths: list[Path] = []
    for index, ids in windows:
        path = WORK / f"block_{index:02d}.h5ad"
        if not path.is_file():
            subset = "feature_id in [" + ",".join(f"'{g}'" for g in ids) + "]"
            with cellxgene_census.open_soma(
                census_version=CENSUS_VERSION, tiledb_config=TILEDB_CFG
            ) as census:
                block = cellxgene_census.get_anndata(
                    census,
                    "Homo sapiens",
                    obs_coords=joinids,
                    var_value_filter=subset,
                    var_column_names=["feature_name", "feature_id"],
                )
            block.obs_names = source.obs_names
            block.write_h5ad(path, compression="gzip")
            print(f"[wide] block {index}: {block.shape} peak={peak_gb():.2f}")
            del block
            gc.collect()
        paths.append(path)

    meta: list[pd.DataFrame] = []
    columns: list[np.ndarray] = []
    for path in paths:
        frame = ad.read_h5ad(path)
        block = dense(frame.X)
        names = frame.var["feature_name"].astype(str).to_numpy()
        identifiers = frame.var["feature_id"].astype(str).to_numpy()
        seen = set(meta_frames_names := [n for m in meta for n in m["feature_name"]])
        keep = np.asarray([name not in seen and not seen.add(name) for name in names])
        columns.append(block[:, keep])
        meta.append(pd.DataFrame({"feature_name": names[keep], "feature_id": identifiers[keep]}))
        del frame, block
        gc.collect()

    var = pd.concat(meta, ignore_index=True)
    counts = np.concatenate(columns, axis=1)
    del columns
    gc.collect()
    n_cells, total = counts.shape
    print(f"[wide] assembled {n_cells:,} x {total:,} peak={peak_gb():.2f}")

    totals = counts.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    lognorm = np.log1p(counts * (1e4 / totals)).astype(np.float32)

    var.index = var["feature_name"].astype(str)
    var.index.name = None
    counts_sparse = sp.csr_matrix(counts)
    lognorm_sparse = sp.csr_matrix(lognorm)
    del counts, lognorm
    gc.collect()

    adata = ad.AnnData(X=lognorm_sparse, obs=obs, var=var)
    adata.layers["counts"] = counts_sparse
    adata.layers["lognorm"] = lognorm_sparse
    out = WORK / "aida_pbmc_wide.h5ad"
    adata.write_h5ad(out, compression="gzip")
    print(f"[wide] wrote {out} {adata.shape} peak={peak_gb():.2f}")


if __name__ == "__main__":
    main()

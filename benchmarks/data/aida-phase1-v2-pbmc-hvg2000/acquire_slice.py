"""
Write a bounded, Cell-Ontology-labelled healthy-human-PBMC benchmark slice for cell-curator.

Dataset : AIDA Phase 1 Data Freeze v2 -- Kock et al., Cell 188(8):2288-2306 (2025)
          doi:10.1016/j.cell.2025.02.017
  collection_id : ced320a1-29f3-47c1-a735-513c7084d508
  dataset_id    : c838aec3-03ef-4398-b882-0e3912abfff0
  census_version: 2025-11-08  (stable; census.cellxgene.cziscience.com/cellxgene-census/v1/release.json)

Prerequisite: run pbmc_hvg_panel.py first (writes /tmp/ccbench/hvg_panel.json).

Memory strategy, and why it looks odd
-------------------------------------
tiledbsoma's peak RSS for an X read is driven by the number of Census tiles a read
touches, not by the returned nnz. Measured on this box (AIDA, 1,265,624 cells):

    20,000 cells, full 61,497 genes, one call ......... 5.51 GB   (over budget)
    20,000 cells, 2,000 genes via explicit var_coords .. 9.75 GB   (over budget)
    20,000 cells, 2,000 genes, one call, scattered .... 9.20 GB   (over budget)
    20,000 cells, 2,000 genes, contiguous slice ....... 2.61 GB   (fits)
    7 contiguous windows x 2,400 cells + small
      scattered reads for the rare cell types .......... 2.2 GB  (fits, used here)

Census obs is NOT blocked by cell type (844,254 interleaved runs for 32 types), so a
single contiguous slice covers only ~6 CL types. The fix is to straddle the joinid
axis with several contiguous windows -- each read is cheap, and their union spans
111 donors and 28+ cell types. The rarest CL types (hematopoietic stem cell, pre-cDC,
CD141+ DC, double-negative Treg, ILC) are too sparse to appear in any window, so they
are pulled in with small scattered reads capped at RARE_CHUNK cells.
"""

import gc
import hashlib
import json
import os
import resource

import anndata as ad
import cellxgene_census
import numpy as np
import scipy.sparse as sp

CENSUS_VERSION = "2025-11-08"
COLLECTION_ID = "ced320a1-29f3-47c1-a735-513c7084d508"
DATASET_ID = "c838aec3-03ef-4398-b882-0e3912abfff0"

N_WINDOWS, PER_WINDOW = 7, 2_400
RARE_N = 3_000          # cell types smaller than this are back-filled explicitly
RARE_KEEP = 100          # cap per rare type, so the back-fill does not distort the prior
RARE_CHUNK = 500        # scattered-read cap; bigger chunks blow the memory budget
SEED = 0
PANEL = "/tmp/ccbench/hvg_panel.json"
OUT = "/tmp/ccbench/aida_pbmc_benchmark.h5ad"
TILEDB_CFG = {
    "soma.init_buffer_bytes": 16_777_216,
    "sm.compute_concurrency_level": 4,
    "sm.io_concurrency_level": 4,
}

# AIDA's own broad-lineage -> sub-population scheme. AIDA reports 8 major types:
# B, CD34+ HSPC, myeloid, NK, plasma, T (platelets are QC-separated from the rest).
LINEAGE = {
    "CL:0001054": "monocyte", "CL:0002396": "monocyte", "CL:0000763": "myeloid",
    "CL:0002399": "myeloid", "CL:0000784": "myeloid", "CL:0002394": "myeloid",
    "CL:0002010": "myeloid",
    "CL:0000895": "T cell", "CL:0000624": "T cell", "CL:0000625": "T cell",
    "CL:0000900": "T cell", "CL:0000904": "T cell", "CL:0000905": "T cell",
    "CL:0000909": "T cell", "CL:0000794": "T cell", "CL:0000934": "T cell",
    "CL:0000084": "T cell", "CL:0000798": "T cell", "CL:0000940": "T cell",
    "CL:0000815": "T cell", "CL:0011024": "T cell",
    "CL:0000939": "NK cell", "CL:0000938": "NK cell", "CL:0000623": "NK cell",
    "CL:0000788": "B cell", "CL:0000787": "B cell", "CL:0000236": "B cell",
    "CL:0000785": "B cell", "CL:0000786": "plasma cell",
    "CL:0001065": "innate lymphoid cell", "CL:0000233": "platelet",
    "CL:0000037": "hematopoietic stem cell",
}
COARSE = {"monocyte": "myeloid"}


def coarse(cl: str) -> str:
    base = LINEAGE.get(cl)
    if base is None:
        raise SystemExit(f"extend LINEAGE for CL id {cl}")
    return COARSE.get(base, base)


def peak() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def rss() -> float:
    for line in open("/proc/self/status"):
        if line.startswith("VmRSS"):
            return int(line.split()[1]) / 1e6
    return 0.0


def main() -> None:
    panel = json.load(open(PANEL))
    hvg_ids = sorted(set(panel["feature_id"]))
    var_filter = "feature_id in [" + ",".join(f"'{g}'" for g in hvg_ids) + "]"
    print(f"[1] HVG panel: {len(hvg_ids):,} Ensembl gene ids (census {CENSUS_VERSION})   "
          f"peak={peak():.2f} rss={rss():.2f}")

    rng = np.random.default_rng(SEED)
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION,
                                    tiledb_config=TILEDB_CFG) as census:
        obs = cellxgene_census.get_obs(
            census, "Homo sapiens",
            value_filter=(f"dataset_id == '{DATASET_ID}' and is_primary_data == True "
                          "and disease == 'normal' and suspension_type == 'cell' "
                          "and tissue == 'blood'"),
            column_names=["soma_joinid", "cell_type", "cell_type_ontology_term_id",
                          "donor_id", "development_stage", "sex"])
        obs["cell_type"] = obs["cell_type"].astype("category").cat.remove_unused_categories()
        obs["cell_type_ontology_term_id"] = \
            obs["cell_type_ontology_term_id"].astype("category").cat.remove_unused_categories()
        obs = obs.sort_values("soma_joinid").reset_index(drop=True)
        cl = obs["cell_type_ontology_term_id"].astype(str)
        print(f"[2] healthy primary blood cells={len(obs):,}  CL cell types={cl.nunique()}   "
              f"peak={peak():.2f} rss={rss():.2f}")
        unmapped = set(cl.unique()) - set(LINEAGE)
        if unmapped:
            raise SystemExit(f"extend LINEAGE for: {sorted(unmapped)}")
        sizes = cl.value_counts()
        print(f"    CL cell types below {RARE_N:,} cells: "
              f"{sorted(sizes[sizes < RARE_N].index)}")

        cache_dir = "/tmp/ccbench/_blocks"
        os.makedirs(cache_dir, exist_ok=True)
        _seq = [0]

        def read(ids: np.ndarray) -> ad.AnnData:
            """One Census read: contiguous joinids, sliced to the HVG panel (cached)."""
            cached = os.path.join(cache_dir, f"block{_seq[0]:02d}.h5ad")
            _seq[0] += 1
            if os.path.isfile(cached):
                return ad.read_h5ad(cached)
            block = cellxgene_census.get_anndata(
                census, "Homo sapiens", obs_coords=np.sort(ids),
                var_value_filter=var_filter,
                var_column_names=["feature_name", "feature_id"],
                obs_column_names=["soma_joinid", "cell_type",
                                  "cell_type_ontology_term_id", "donor_id",
                                  "development_stage", "sex"])
            cols = np.flatnonzero(np.isin(block.var["feature_id"].astype(str).to_numpy(),
                                          hvg_ids))
            out = block[:, cols].copy()
            out.write_h5ad(cached, compression="gzip")
            return out

        blocks = []
        # (a) straddle the joinid axis with contiguous windows
        for k, start in enumerate(np.linspace(0, len(obs) - PER_WINDOW, N_WINDOWS).astype(int)):
            win = obs.iloc[start:start + PER_WINDOW]
            blocks.append(read(win["soma_joinid"].to_numpy()))
            print(f"    window {k + 1}/{N_WINDOWS} @{start:,}: {PER_WINDOW:,} cells, "
                  f"{win['donor_id'].nunique()} donors, "
                  f"{win['cell_type_ontology_term_id'].nunique()} CL types   "
                  f"peak={peak():.2f} rss={rss():.2f}")
            del win
            gc.collect()

        # (b) back-fill the rare CL types with small scattered reads
        rare_rows = obs.index[cl.isin(sizes[sizes < RARE_N].index).to_numpy()].to_numpy()
        rare_sub = obs.loc[rare_rows]
        rare_ids = np.sort(np.concatenate([
            rng.choice(g["soma_joinid"].to_numpy(),
                       min(RARE_KEEP, len(g)), replace=False)
            for _, g in rare_sub.groupby(rare_sub["cell_type_ontology_term_id"].astype(str),
                                         observed=True)]))
        del rare_rows, rare_sub
        gc.collect()
        for k in range(0, rare_ids.size, RARE_CHUNK):
            blocks.append(read(rare_ids[k:k + RARE_CHUNK]))
        print(f"    rare-cell-type back-fill: {rare_ids.size:,} cells "
              f"(<= {RARE_KEEP} per rare CL type) in "
              f"{-(-rare_ids.size // RARE_CHUNK)} scattered reads   "
              f"peak={peak():.2f} rss={rss():.2f}")
        del obs, cl, sizes, rare_ids
        gc.collect()


    # ad.concat(join="outer") does not carry var identity, so keep it aside
    var_index = blocks[0].var.loc[:, [c for c in ("feature_name", "feature_id")
                                       if c in blocks[0].var.columns]].copy()
    assert all(list(b.var.index) == list(var_index.index) for b in blocks), \
        "HVG panel differs between reads"
    adata = ad.concat(blocks, join="outer", index_unique=None)
    del blocks
    gc.collect()
    adata.var = var_index
    # a rare cell type back-filled cell can also sit inside a window -> keep one copy
    n_before = adata.n_obs
    adata = adata[~adata.obs["soma_joinid"].astype(np.int64).duplicated()].copy()
    print(f"    dropped {n_before - adata.n_obs:,} cells duplicated across reads")

    # derived layers + the 2-level truth path the benchmark harness consumes
    adata.var_names = adata.var["feature_name"].astype(str).to_numpy()
    adata.var_names_make_unique()
    counts = adata.X.tocsr().astype(np.float32)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    adata.layers["counts"] = counts
    adata.layers["lognorm"] = sp.csr_matrix(
        np.log1p(counts.multiply((1e4 / totals)[:, None]).tocsr()))
    adata.obs["cell_type_ontology_term_id"] = \
        adata.obs["cell_type_ontology_term_id"].astype(str)
    adata.obs["lineage"] = [coarse(c) for c in adata.obs["cell_type_ontology_term_id"]]
    adata.obs["true_path"] = adata.obs["lineage"] + "/" + adata.obs["cell_type"].astype(str)
    adata.obs["collection_id"] = COLLECTION_ID
    adata.obs["dataset_id"] = DATASET_ID
    adata.obs["census_version"] = CENSUS_VERSION
    adata.obs_names = [f"{d}_{i}" for d, i in zip(adata.obs["donor_id"].astype(str),
                                                  adata.obs_names)]
    adata.obs_names_make_unique()
    for col in adata.obs.columns:
        if adata.obs[col].dtype.name == "category":
            adata.obs[col] = adata.obs[col].astype("category").cat.remove_unused_categories()

    print(f"[3] AnnData {adata.shape} nnz={adata.X.nnz:,} "
          f"density={adata.X.nnz / (adata.shape[0] * adata.shape[1]):.2%} "
          f"layers={list(adata.layers)}")
    print(f"    genes[:12]={list(adata.var_names[:12])}")
    print(f"    donors={adata.obs['donor_id'].nunique()}  "
          f"CL types={adata.obs['cell_type_ontology_term_id'].nunique()}  "
          f"lineages={sorted(adata.obs['lineage'].unique())}")
    counts_by_type = (adata.obs.groupby(["cell_type", "cell_type_ontology_term_id"],
                                        observed=True).size().sort_values(ascending=False))
    for (label, cl_id), n in counts_by_type.items():
        print(f"      {n:>6,}  {cl_id:<26} {label}")

    adata.write_h5ad(OUT, compression="gzip")
    digest = hashlib.sha256(open(OUT, "rb").read()).hexdigest()
    print(f"[4] wrote {OUT}  {os.path.getsize(OUT) / 1e6:.1f} MB")
    print(f"    sha256 = {digest}")
    print(f"[5] FINAL PEAK RSS = {peak():.2f} GB / 4.00 GB budget")


if __name__ == "__main__":
    main()

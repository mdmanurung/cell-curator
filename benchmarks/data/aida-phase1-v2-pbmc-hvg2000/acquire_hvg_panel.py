"""Pass 1 (own process): rank HVGs on a small stratified healthy-PBMC seed.

Isolated so its peak heap is returned to the OS before the benchmark slice is
fetched. Writes the chosen var soma_joinids + feature names to JSON.
"""
import gc, json, resource
import cellxgene_census, numpy as np

CENSUS_VERSION = "2025-11-08"
DATASET_ID    = "c838aec3-03ef-4398-b882-0e3912abfff0"
N_HVG, N_SEED_CELLS, FLOOR, SEED = 2_000, 1_000, 10, 0
OUT = "/tmp/ccbench/hvg_panel.json"
TILEDB_CFG = {"soma.init_buffer_bytes": 16_777_216,
              "sm.compute_concurrency_level": 4,
              "sm.io_concurrency_level": 4}

def peak(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
def rss():
    for line in open("/proc/self/status"):
        if line.startswith("VmRSS"):
            return int(line.split()[1]) / 1e6

rng = np.random.default_rng(SEED)
with cellxgene_census.open_soma(census_version=CENSUS_VERSION,
                                tiledb_config=TILEDB_CFG) as census:
    # (a) gene panel, addressed explicitly afterwards so pass 2 can reuse it
    var = census["census_data"]["homo_sapiens"].ms["RNA"]["var"].read(
        coords=(slice(None),)).concat().to_pandas()
    var = var[(var["feature_type"] == "protein_coding") & (var["n_measured_obs"] > 1000)]
    panel_ids = np.sort(var["soma_joinid"].to_numpy())
    panel = var.sort_values("soma_joinid")
    panel_names = panel["feature_name"].astype(str).to_numpy()
    panel_fids = panel["feature_id"].astype(str).to_numpy()
    del var, panel; gc.collect()
    print(f"[panel] gene panel = {panel_ids.size:,} protein-coding genes   peak={peak():.2f} rss={rss():.2f}")

    # (b) stratified seed cells
    obs = cellxgene_census.get_obs(
        census, "Homo sapiens",
        value_filter=(f"dataset_id == '{DATASET_ID}' and is_primary_data == True "
                      "and disease == 'normal' and suspension_type == 'cell' and tissue == 'blood'"),
        column_names=["soma_joinid", "cell_type_ontology_term_id"])
    groups = [g.index.to_numpy() for _, g in
              obs.groupby("cell_type_ontology_term_id", observed=True)]
    take = np.concatenate([rng.choice(g, min(FLOOR, len(g)), replace=False) for g in groups])
    rest = np.setdiff1d(obs.index.to_numpy(), take)
    rows = np.sort(np.concatenate([take, rng.choice(rest, max(0, N_SEED_CELLS - take.size), replace=False)]))
    seed_ids = np.sort(obs.loc[rows, "soma_joinid"].to_numpy())
    del obs, groups, take, rest, rows; gc.collect()
    print(f"[panel] seed cells = {seed_ids.size:,} over {len(set())or ''}CL types   peak={peak():.2f} rss={rss():.2f}")

    # (c) seed counts, then drop everything else before computing moments
    seed = cellxgene_census.get_anndata(census, "Homo sapiens", obs_coords=seed_ids,
                                        var_coords=panel_ids, obs_column_names=[])
    X = seed.X.tocsc(); del seed; gc.collect()
    print(f"[panel] seed matrix {X.shape} nnz={X.nnz:,}   peak={peak():.2f} rss={rss():.2f}")

    # (d) per-gene mean / variance of the raw counts, without densifying
    data, indices, indptr = X.data.astype(np.float64), X.indices, X.indptr
    n_cells = X.shape[0]
    sum_g = np.bincount(np.repeat(np.arange(X.shape[1]), np.diff(indptr)),
                        weights=data, minlength=X.shape[1])
    sq_g = np.bincount(np.repeat(np.arange(X.shape[1]), np.diff(indptr)),
                       weights=data * data, minlength=X.shape[1])
    detected = np.bincount(np.repeat(np.arange(X.shape[1]), np.diff(indptr)),
                           minlength=X.shape[1])
    del data, indices, indptr, X; gc.collect()
    mean = sum_g / n_cells
    var_g = np.maximum(sq_g / n_cells - mean**2, 0.0)
    disp = var_g / np.maximum(mean, 1e-8)
    keep = (detected >= 5) & (mean > 0.05)
    picked = np.where(keep)[0][np.argsort(-disp[keep])[:N_HVG]]

    json.dump({"var_soma_joinid": [int(x) for x in panel_ids[picked]],
               "feature_name": [str(x) for x in panel_names[picked]],
               "feature_id": [str(x) for x in panel_fids[picked]],
               "n_hvg": int(len(picked)), "census_version": CENSUS_VERSION,
               "dataset_id": DATASET_ID, "seed": SEED,
               "panel_size": int(len(panel_ids))},
              open(OUT, "w"), indent=1)
    print(f"[panel] {len(picked)} HVGs selected from {int(keep.sum()):,} variable genes -> {OUT}")
    print(f"[panel] first 15: {', '.join(panel_names[picked][:15])}")
    print(f"[panel] PEAK = {peak():.2f} GB   rss={rss():.2f} GB")

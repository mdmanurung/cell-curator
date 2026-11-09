"""Add a counts layer to a v8 query so scANVI and popV can read it.

Purpose: the v8 query input.h5ad holds log1p CP10k only; scANVI and popV model counts.
Run:     python benchmarks/scripts/v8_prepare_query.py QUERY_H5AD OUT_H5AD [--check]
Inputs:
    QUERY_H5AD   input.h5ad [prepare_novel_v8.py]: log1p CP10k in X, symbols as
                 var_names, obs donor_id and bench_cluster
Outputs:
    OUT_H5AD     same cells and genes; layers["counts"] as given, or recovered from X
                 (see v8_common.recover_counts)
"""

from v8_common import parse_args, recover_counts


def main():
    args = parse_args("counts layer for a v8 query", ["query_h5ad"], ["out_h5ad"])
    import anndata as ad

    query = ad.read_h5ad(args.query_h5ad)
    if not query.var_names.is_unique:
        raise ValueError("query gene symbols are not unique")
    # v8 inputs carry raw counts (prepare_novel_v8.py); recover only for older inputs
    if "counts" not in query.layers:
        query.layers["counts"] = recover_counts(query.X)
    query.write_h5ad(args.out_h5ad)
    print(f"{query.n_obs} cells, {query.n_vars} genes")


if __name__ == "__main__":
    main()

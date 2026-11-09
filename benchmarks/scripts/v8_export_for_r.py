"""Write an h5ad's log-normalized matrix in a form R reads without Matrix Market text.

Purpose: SingleR trains on about 200k cells; a text .mtx of that is several GB.
Run:     python benchmarks/scripts/v8_export_for_r.py H5AD OUT_DIR [--check]
Inputs:
    H5AD    train.h5ad or a query input.h5ad (log1p CP10k in X)
Outputs (OUT_DIR, must not exist):
    p.bin, i.bin, x.bin   cells x genes CSR arrays = genes x cells CSC (int32, int32, float32)
    dim.txt               n_genes, n_cells
    genes.tsv, barcodes.tsv
    labels.tsv            barcode, cell_type (only when obs has cell_type)
"""

from v8_common import parse_args


def main():
    args = parse_args("R export of a log-normalized h5ad", ["h5ad"], ["out_dir"])
    import anndata as ad
    import numpy as np
    from scipy import sparse

    data = ad.read_h5ad(args.h5ad)
    matrix = sparse.csr_matrix(data.X)
    args.out_dir.mkdir(parents=True)
    matrix.indptr.astype(np.int32).tofile(args.out_dir / "p.bin")
    matrix.indices.astype(np.int32).tofile(args.out_dir / "i.bin")
    matrix.data.astype(np.float32).tofile(args.out_dir / "x.bin")
    (args.out_dir / "dim.txt").write_text(f"{data.n_vars}\n{data.n_obs}\n")
    (args.out_dir / "genes.tsv").write_text("\n".join(data.var_names) + "\n")
    (args.out_dir / "barcodes.tsv").write_text("\n".join(data.obs_names) + "\n")
    if "cell_type" in data.obs:
        data.obs[["cell_type"]].rename_axis("barcode").reset_index().to_csv(
            args.out_dir / "labels.tsv", sep="\t", index=False
        )


if __name__ == "__main__":
    main()

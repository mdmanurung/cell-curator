"""Write one comparator's predictions in the neutral benchmark format.

Purpose: v8 truth has other column names than comparator_to_neutral.py expects, and
         the Cell Ontology ids of predicted labels come from the training label map.
         This renames the truth columns, then calls comparator_to_neutral.py.
Run:     .venv/bin/python benchmarks/scripts/v8_to_neutral.py \
             PREDICTIONS QUERY_H5AD TRUTH_TSV LABEL_MAP OUT_TSV [--check]
         Scoring step: it reads the truth file, so run it only when scoring is due.
Inputs:
    PREDICTIONS  TSV: barcode, predicted_label, uncertainty
    QUERY_H5AD   input.h5ad; its cells are every cell the comparator was asked to label
    TRUTH_TSV    prepare_novel_v8.py scoring/truth.tsv: barcode, cell_type, cell_type_ontology_term_id
    LABEL_MAP    TRAIN_DIR/label_map.tsv: label, cl_id
Outputs:
    OUT_TSV      neutral TSV (see comparator_to_neutral.py); OUT_TSV.truth.tsv and
                 OUT_TSV.barcodes.txt are its inputs
"""

import subprocess
import sys
from pathlib import Path

from v8_common import parse_args

# Parameters -----------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
ONTOLOGY = (  # cl.obo release 2026-06-08, pinned by the protocol
    HERE.parents[0] / "data/aida-phase1-v2-pbmc-hvg2000/cl.obo"
)


def main():
    args = parse_args(
        "neutral table for a v8 comparator",
        ["predictions", "query_h5ad", "truth_tsv", "label_map"],
        ["out_tsv"],
    )
    import anndata as ad
    import pandas as pd

    truth = pd.read_csv(args.truth_tsv, sep="\t", dtype=str, keep_default_na=False)
    truth = truth.rename(
        columns={"cell_type": "true_label", "cell_type_ontology_term_id": "true_cl_id"}
    )
    truth_path = Path(f"{args.out_tsv}.truth.tsv")
    truth.to_csv(truth_path, sep="\t", index=False)
    barcodes = ad.read_h5ad(args.query_h5ad, backed="r").obs_names
    barcodes_path = Path(f"{args.out_tsv}.barcodes.txt")
    barcodes_path.write_text("\n".join(barcodes) + "\n")
    subprocess.run(
        [
            sys.executable, str(HERE / "comparator_to_neutral.py"),
            str(args.predictions), str(truth_path), str(ONTOLOGY), str(args.out_tsv),
            "--label-map", str(args.label_map), "--barcodes", str(barcodes_path),
        ],
        check=True,
    )  # fmt: skip


if __name__ == "__main__":
    main()

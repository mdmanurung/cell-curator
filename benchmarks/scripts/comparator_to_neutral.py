"""Turn one comparator's per-cell predictions into a neutral benchmark table.

Purpose: give every external system (CellTypist, SingleR, scANVI, popV, and the LLM
         annotators) the same path into `cell-curator benchmark evaluate`, with one
         label-to-Cell-Ontology rule declared before any output is seen.
Run:     uv run python benchmarks/scripts/comparator_to_neutral.py \
             PREDICTIONS TRUTH ONTOLOGY OUT [--label-map MAP_TSV]
Inputs:
    PREDICTIONS   TSV: barcode, predicted_label, and optional uncertainty in [0, 1]
                  [crosscheck celltypist-predict, reference map-scarches, or a
                  benchmarks/scripts runner]
    TRUTH         TSV: barcode, true_label, true_cl_id [ground_truth.tsv]
    ONTOLOGY      Cell Ontology OBO snapshot [cl.obo, release pinned in the protocol]
    MAP_TSV       optional TSV: label, cl_id, written before the comparator ran
    BARCODES      optional, every cell the system was asked to label [splits/B/barcodes.tsv]
Outputs:
    OUT           neutral TSV: barcode, true_label, predicted_label, true_cl_id,
                  predicted_cl_id, and confidence when the comparator gave uncertainty
    printed       rows, called fraction, and the fraction of called labels with no CL id

Label rule: a declared map wins; otherwise the label must resolve in the ontology as
an exact name or EXACT synonym, case-insensitively, after trying the singular of a
label ending in " cells". A label that does not resolve keeps
its text and gets no CL id, so it scores as wrong on exact metrics and is skipped by
the ontology metrics. The unmapped fraction is printed so it is reported, not hidden.
"""

import argparse
from pathlib import Path

import pandas as pd

from cell_curator.contracts import ContractError
from cell_curator.ontology import CellOntology

parser = argparse.ArgumentParser()
parser.add_argument("predictions", type=Path)
parser.add_argument("truth", type=Path)
parser.add_argument("ontology", type=Path)
parser.add_argument("out", type=Path)
parser.add_argument("--label-map", type=Path)
parser.add_argument(
    "--barcodes",
    type=Path,
    help="every cell the system was asked to label, one per line; a missing one is Unknown",
)
args = parser.parse_args()

if args.out.exists():
    raise ValueError(f"refusing to overwrite {args.out}")

# Join predictions to truth --------------------------------------------------------

predictions = pd.read_csv(args.predictions, sep="\t", keep_default_na=False, dtype=str)
if args.barcodes:
    # A cell the system declined to label is an abstention, not a cell to drop.
    asked = pd.Index(args.barcodes.read_text().split(), name="barcode")
    if not predictions["barcode"].isin(asked).all():
        raise ValueError("predictions contain barcodes the system was not asked to label")
    missing = asked.difference(predictions["barcode"])
    filler = {"barcode": missing, "predicted_label": "Unknown"}
    if "uncertainty" in predictions:
        filler["uncertainty"] = "1.0"
    predictions = pd.concat([predictions, pd.DataFrame(filler)], ignore_index=True)
    print(f"{len(missing)} cells the system did not label are Unknown")
truth = pd.read_csv(args.truth, sep="\t", keep_default_na=False, dtype=str)
table = truth.merge(predictions, on="barcode", how="inner", validate="1:1")
if len(table) != len(predictions):
    raise ValueError(f"{len(predictions) - len(table)} predicted barcodes are not in the truth")

# Map labels to the Cell Ontology ---------------------------------------------------

ontology = CellOntology.load(args.ontology)
declared = (
    pd.read_csv(args.label_map, sep="\t", dtype=str).set_index("label")["cl_id"].to_dict()
    if args.label_map
    else {}
)


def cl_id(label: str) -> str:
    if label == "Unknown":
        return ""
    if label in declared:
        return ontology.resolve(declared[label]).accession
    # Free-text annotators answer in plurals ("NK cells"); the singular is tried too.
    candidates = [label, label[:-1]] if label.endswith(" cells") else [label]
    for candidate in candidates:
        try:
            return ontology.resolve(candidate).accession
        except ContractError:
            continue
    return ""


ids = {label: cl_id(label) for label in table["predicted_label"].unique()}
table["predicted_cl_id"] = table["predicted_label"].map(ids)

# Write ------------------------------------------------------------------------------

columns = ["barcode", "true_label", "predicted_label", "true_cl_id", "predicted_cl_id"]
if "uncertainty" in table and table["uncertainty"].ne("").all():
    table["confidence"] = 1.0 - table["uncertainty"].astype(float)
    columns.append("confidence")
args.out.parent.mkdir(parents=True, exist_ok=True)
table[columns].to_csv(args.out, sep="\t", index=False)

called = table["predicted_label"].ne("Unknown")
unmapped = called & table["predicted_cl_id"].eq("")
print(f"rows: {len(table)}")
print(f"called: {called.mean():.4f}")
print(f"called with no CL id: {unmapped.sum() / max(called.sum(), 1):.4f}")
print("unmapped labels:", sorted(table.loc[unmapped, "predicted_label"].unique()))

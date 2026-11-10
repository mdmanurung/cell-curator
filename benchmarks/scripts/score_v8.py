# Purpose: score one method on one protocol v8 dataset with the package's own v6 metrics
#          (benchmark.py), plus report-only strata by the truth label's relation to the
#          38 corpus terms (amendment 3 item 4, amendment 2 item 2). No metric is
#          reimplemented here.
# Run:     .venv/bin/python benchmarks/scripts/score_v8.py \
#              PREDICTIONS TRUTH OUT_DIR --method NAME --dataset NAME
# Inputs:
#   PREDICTIONS  TSV, one row per cell (a cluster call is repeated on each of its
#                cells, as the v6 export does): barcode, predicted_cl_id, and
#                optional predicted_label, confidence, abstain [method runner]
#   TRUTH        TSV: barcode, cell_type, cell_type_ontology_term_id [dataset prep;
#                novel-set truth is read only at scoring time]
#   CORPUS_IDS   benchmarks/data/aifi-hiha-dev-split/corpus_cl_ids.tsv (id, name)
#   ONTOLOGY     benchmarks/data/aida-phase1-v2-pbmc-hvg2000/cl.obo, the pinned
#                Cell Ontology snapshot of v5 and v6
# Outputs:
#   OUT_DIR/{dataset}.{method}.json  v6 metrics, risk-coverage curves (strict,
#                descendant, relaxed), strata, and input hashes. Never overwritten.

import argparse
from pathlib import Path

import pandas as pd

from cell_curator.benchmark import evaluate_predictions, risk_coverage_sweep
from cell_curator.contracts import sha256
from cell_curator.ontology import CellOntology
from cell_curator.provenance import publish_json

# Parameters -----------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
CORPUS_IDS = ROOT / "benchmarks/data/aifi-hiha-dev-split/corpus_cl_ids.tsv"
ONTOLOGY = ROOT / "benchmarks/data/aida-phase1-v2-pbmc-hvg2000/cl.obo"  # pinned in v5/v6
UNKNOWN = "Unknown"  # the label benchmark.py reads as an abstention
CURVES = ("strict", "descendant", "relaxed")  # v6 published strict and descendant
ABSTAIN_ON_CURVE = True  # v6 sweeps used --abstain-unknown: an abstention is not a call
STRATUM_FIELDS = (  # v6 fields repeated per stratum
    "n",
    "accuracy",
    "coverage",
    "selective_accuracy",
    "descendant_accuracy",
    "relaxed_accuracy",
)
TRUE_VALUES = {"1", "true", "yes", "t"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("predictions", type=Path)
    parser.add_argument("truth", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--method", required=True)
    parser.add_argument("--dataset", required=True)
    args = parser.parse_args()

    ontology = CellOntology.load(ONTOLOGY)

    # Join predictions to truth ------------------------------------------------------
    truth = pd.read_csv(args.truth, sep="\t", keep_default_na=False, dtype=str)
    preds = pd.read_csv(args.predictions, sep="\t", keep_default_na=False, dtype=str)
    if not preds["barcode"].isin(truth["barcode"]).all():
        raise ValueError("predictions contain barcodes that are not in the truth")
    if "predicted_label" not in preds:
        preds["predicted_label"] = ""
    if "confidence" not in preds:
        preds["confidence"] = ""
    if "abstain" not in preds:
        preds["abstain"] = ""
    # A truth cell with no prediction row is an abstention, as in comparator_to_neutral.
    frame = truth.merge(preds, on="barcode", how="left", validate="1:1")
    n_unpredicted = int(frame["predicted_cl_id"].isna().sum())
    frame = frame.fillna({"predicted_cl_id": "", "predicted_label": "", "confidence": ""})

    # An abstention has no term and the label Unknown (benchmark.py contract) --------
    abstain = (
        frame["abstain"].str.lower().isin(TRUE_VALUES)
        | frame["predicted_cl_id"].str.strip().eq("")
        | frame["predicted_label"].eq(UNKNOWN)
    )
    names = {
        accession: ontology.resolve(accession).name
        for accession in frame.loc[~abstain, "predicted_cl_id"].unique()
    }
    frame["predicted_label"] = frame["predicted_label"].where(
        frame["predicted_label"].ne(""), frame["predicted_cl_id"].map(names)
    )
    frame.loc[abstain, ["predicted_label", "predicted_cl_id"]] = [UNKNOWN, ""]
    # One unmeasured confidence makes the whole column absent, so no curve is faked.
    has_confidence = frame["confidence"].ne("").all()
    frame["confidence"] = pd.to_numeric(frame["confidence"], errors="coerce")  # "" -> NaN
    frame = frame.rename(
        columns={"cell_type": "true_label", "cell_type_ontology_term_id": "true_cl_id"}
    )
    columns = ["barcode", "true_label", "predicted_label", "true_cl_id", "predicted_cl_id"]
    if has_confidence:
        columns.append("confidence")
    frame = frame[columns]

    # v6 metrics, computed by the package ----------------------------------------------
    result = evaluate_predictions(
        frame, system=args.method, dataset=args.dataset, ontology=ontology
    )
    curves = {}
    if has_confidence:
        for agreement in CURVES:
            points = risk_coverage_sweep(
                frame,
                agreement=agreement,
                ontology=ontology,
                abstain_unknown=ABSTAIN_ON_CURVE,
            )
            curves[agreement] = [point.model_dump(mode="json") for point in points]

    # Strata by the truth label's relation to the corpus terms -------------------------
    corpus = pd.read_csv(CORPUS_IDS, sep="\t", dtype=str)
    corpus_ids = {ontology.resolve(term).accession for term in corpus["id"]}
    corpus_ancestors = set().union(*(ontology.ancestors(term) for term in corpus_ids))
    true_ids = {term: ontology.resolve(term).accession for term in frame["true_cl_id"].unique()}
    relation = {}
    for term, accession in true_ids.items():
        if accession in corpus_ids:
            relation[term] = "corpus_term"
        elif accession in corpus_ancestors:
            relation[term] = "coarser_parent"
        elif corpus_ids & ontology.ancestors(accession):
            relation[term] = "finer_child"
        else:
            relation[term] = "outside"
    frame["stratum"] = frame["true_cl_id"].map(relation)

    strata = {}
    for name in ("corpus_term", "coarser_parent", "finer_child", "outside"):
        part = frame[frame["stratum"].eq(name)].reset_index(drop=True)
        if part.empty:
            strata[name] = {"n": 0}
        elif name == "outside":
            # Abstaining on an unrelated label is the right answer here.
            abstained = int(part["predicted_label"].eq(UNKNOWN).sum())
            strata[name] = {
                "n": len(part),
                "n_abstained": abstained,
                "abstention_correct": abstained / len(part),
            }
        else:
            scored = evaluate_predictions(
                part, system=args.method, dataset=args.dataset, ontology=ontology
            ).model_dump(mode="json")
            strata[name] = {field: scored[field] for field in STRATUM_FIELDS}
        strata[name]["true_cl_ids"] = sorted(
            term for term, group in relation.items() if group == name
        )

    # Write one file per method and dataset --------------------------------------------
    output = {
        **result.model_dump(mode="json"),
        "risk_coverage": curves,
        "risk_coverage_abstain_unknown": ABSTAIN_ON_CURVE,
        "strata": strata,
        "n_truth_cells_without_prediction": n_unpredicted,
        "prediction_path": args.predictions.name,
        "prediction_sha256": sha256(args.predictions),
        "truth_path": args.truth.name,
        "truth_sha256": sha256(args.truth),
        "corpus_ids_sha256": sha256(CORPUS_IDS),
        "ontology_source": ONTOLOGY.name,
        "ontology_sha256": sha256(ONTOLOGY),
    }
    publish_json(args.out_dir / f"{args.dataset}.{args.method}.json", output)


if __name__ == "__main__":
    main()

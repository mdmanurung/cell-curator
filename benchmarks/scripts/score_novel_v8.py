# Purpose: the one-shot scoring of a protocol v8 novel set: consensus (strict and best guess) and the
#          three tools, on the blood cells of a pooled run, with both pass bars.
#          Decisions (plans/v8-agent-decisions.md, 2026-11-17): every non-blood cell is set to Unknown
#          for every method; the bars are computed on blood cells only; a truth term the ontology
#          cannot resolve is reported, not fatal; failures are reported, never tuned.
# Run:     .venv/bin/python benchmarks/scripts/score_novel_v8.py DATASET INPUT_H5AD TRUTH CONSENSUS_DIR \
#              CELLTYPIST.tsv POPV.tsv SINGLER.tsv OUT_DIR [--tissue-column tissue --blood-value blood]
#          This is the only step that reads the truth. Run it once, after the consensus runs verify.
# Inputs:  INPUT_H5AD     prepared input; obs holds the tissue column (read for barcodes and tissue only)
#          TRUTH          scoring/truth.tsv: barcode, cell_type, cell_type_ontology_term_id
#          CONSENSUS_DIR  annotation-consensus output: calls.tsv, consensus.tsv, cell_flags.tsv
#          tool tables    barcode, predicted_label, predicted_cl_id, confidence [v8_tools_for_consensus.py]
# Outputs (OUT_DIR, must not exist):
#   predictions/*.tsv         the masked predictions that were scored (all cells)
#   blood_barcodes.txt        the scored population
#   scores/{DATASET}.{method}.json         score_v8.py output on blood cells
#   scores_allmasked/{DATASET}_allmasked.{method}.json   report only: all cells, non-blood as Unknown
#   compare/                  matched_coverage.tsv and risk_coverage.png (compare_v8.py, blood cells)
#   bars.json, bars.tsv       both pass bars and the recorded-reason check

import argparse
import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

# Parameters -----------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
UNKNOWN = "Unknown"
TOOLS = ("celltypist", "popv", "singler")
PROTOCOL_MIN_DESCENDANT = 0.60  # README amendment 1 item 3 (the v6 "Not comparable" line)
CONSENSUS_COVERAGE = 0.95  # plans/v8-consensus-mode.md row 1 of the round-1 decisions
CONSENSUS_MARGIN = 0.01  # "within 1 point of the best single tool"
CURVE = "descendant"  # the agreement level the v6 stop rules used


def curve(score: dict) -> pd.DataFrame:
    points = pd.DataFrame(score["risk_coverage"][CURVE])
    points = points[points["risk"].notna()]
    return points.assign(accuracy=1 - points["risk"])[["coverage", "accuracy"]]


def at_coverage(score: dict, target: float) -> float:
    points = curve(score)
    reachable = points[points["coverage"] >= target - 1e-9]
    return float(reachable["accuracy"].max()) if len(reachable) else float("nan")


def dominates(rival: dict, mine: dict) -> bool:
    """The rival's curve lies at or above mine at every point of my range, and above at one."""
    gaps = [at_coverage(rival, c) - a for c, a in curve(mine)[["coverage", "accuracy"]].values]
    if not gaps or any(pd.isna(g) or g < -1e-12 for g in gaps):
        return False
    return any(g > 1e-12 for g in gaps)


def prediction_table(frame: pd.DataFrame, blood: set[str]) -> pd.DataFrame:
    """barcode, predicted_cl_id, predicted_label, confidence; non-blood cells are Unknown."""
    out = frame[["barcode", "predicted_cl_id", "predicted_label", "confidence"]].copy()
    masked = ~out["barcode"].isin(blood)
    out.loc[masked, ["predicted_cl_id", "predicted_label"]] = ["", UNKNOWN]
    return out


def bars(scores: dict[str, dict], consensus_dir: Path, blood: set[str]) -> dict:
    mine = scores["consensus_strict"]
    rivals = {name: scores[name] for name in TOOLS}
    # Protocol bar ---------------------------------------------------------------------
    not_comparable = bool(mine["descendant_accuracy"] < PROTOCOL_MIN_DESCENDANT)
    dominators = sorted(name for name, score in rivals.items() if dominates(score, mine))
    protocol = {
        "descendant_accuracy": mine["descendant_accuracy"],
        "not_comparable": not_comparable,
        "dominated_by": dominators,
        "pass": (not not_comparable) and not dominators,
    }
    # Consensus bar: within 1 point of the best tool at coverage >= 0.95, and a reason per call --
    best = {name: at_coverage(score, CONSENSUS_COVERAGE) for name, score in rivals.items()}
    best_tool = max(best, key=lambda name: -1 if pd.isna(best[name]) else best[name])
    mine_at = at_coverage(mine, CONSENSUS_COVERAGE)
    reached = not pd.isna(mine_at) and not pd.isna(best[best_tool])
    gap = float(best[best_tool] - mine_at) if reached else None
    guess_at = at_coverage(scores["consensus_best_guess"], CONSENSUS_COVERAGE)  # report only
    consensus_table = pd.read_csv(
        consensus_dir / "consensus.tsv", sep="\t", keep_default_na=False, dtype=str
    )
    flags = pd.read_csv(
        consensus_dir / "cell_flags.tsv", sep="\t", keep_default_na=False, dtype=str
    )
    clusters_without_reason = int(
        (
            consensus_table["rationale"].str.strip().eq("")
            | consensus_table["decision_rule"].str.strip().eq("")
        ).sum()
    )
    flagged_without_reason = int(flags["reason"].str.strip().eq("").sum())
    consensus = {
        "coverage_target": CONSENSUS_COVERAGE,
        "consensus_accuracy_at_target": None if pd.isna(mine_at) else mine_at,
        "best_tool": best_tool,
        "best_guess_accuracy_at_target_report_only": None if pd.isna(guess_at) else guess_at,
        "strict_coverage": mine["coverage"],
        "best_tool_accuracy_at_target": None if pd.isna(best[best_tool]) else best[best_tool],
        "gap": gap,
        "target_reached_by_all": reached,
        "clusters_without_reason": clusters_without_reason,
        "flagged_cells_without_reason": flagged_without_reason,
        # None (NA) when strict calls never reach the target coverage: the bar is not met or failed
        "pass": None
        if not reached
        else bool(
            gap <= CONSENSUS_MARGIN + 1e-12
            and clusters_without_reason == 0
            and flagged_without_reason == 0
        ),
    }
    return {"n_blood_cells": len(blood), "protocol_bar": protocol, "consensus_bar": consensus}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset")
    parser.add_argument("input_h5ad", type=Path)
    parser.add_argument("truth", type=Path)
    parser.add_argument("consensus_dir", type=Path)
    parser.add_argument("celltypist", type=Path)
    parser.add_argument("popv", type=Path)
    parser.add_argument("singler", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--tissue-column", default="tissue")
    parser.add_argument("--blood-value", default="blood")
    args = parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"{args.out_dir} exists; this step runs once and never overwrites")

    obs = ad.read_h5ad(args.input_h5ad, backed="r").obs
    blood = set(obs.index[obs[args.tissue_column].astype(str).eq(args.blood_value)])
    if not blood:
        raise SystemExit(f"no cell has {args.tissue_column} == {args.blood_value}")
    for name in ("predictions", "scores", "scores_allmasked", "compare"):
        (args.out_dir / name).mkdir(parents=True)
    (args.out_dir / "blood_barcodes.txt").write_text("\n".join(sorted(blood)) + "\n")

    # Predictions of the five methods, non-blood set to Unknown -----------------------------
    calls = pd.read_csv(
        args.consensus_dir / "calls.tsv", sep="\t", keep_default_na=False, dtype=str
    )
    strict = calls.rename(
        columns={
            "label_strict": "predicted_label",
            "cl_id_strict": "predicted_cl_id",
            "winner_share": "confidence",
        }
    )
    best_guess = calls.rename(
        columns={
            "label_best_guess": "predicted_label",
            "cl_id_best_guess": "predicted_cl_id",
            "winner_share": "confidence",
        }
    )
    tables = {"consensus_strict": strict, "consensus_best_guess": best_guess}
    for name, path in zip(TOOLS, (args.celltypist, args.popv, args.singler), strict=True):
        tables[name] = pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str)
    scorer = HERE / "score_v8.py"
    for name, frame in tables.items():
        path = args.out_dir / "predictions" / f"{name}.tsv"
        prediction_table(frame, blood).to_csv(path, sep="\t", index=False)
        for folder, dataset, extra in (
            ("scores", args.dataset, ["--restrict", str(args.out_dir / "blood_barcodes.txt")]),
            ("scores_allmasked", f"{args.dataset}_allmasked", []),
        ):
            subprocess.run(
                [
                    sys.executable,
                    str(scorer),
                    str(path),
                    str(args.truth),
                    str(args.out_dir / folder),
                    "--method",
                    name,
                    "--dataset",
                    dataset,
                    *extra,
                ],
                check=True,
            )

    # Matched-coverage comparison and both bars, on blood cells ------------------------------
    files = sorted((args.out_dir / "scores").glob("*.json"))
    subprocess.run(
        [
            sys.executable,
            str(HERE / "compare_v8.py"),
            str(args.out_dir / "compare"),
            *map(str, files),
        ],
        check=True,
    )
    scores = {path.stem.split(".", 1)[1]: json.loads(path.read_text()) for path in files}
    result = bars(scores, args.consensus_dir, blood)
    (args.out_dir / "bars.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    rows = [
        {
            "bar": bar,
            **{
                k: json.dumps(v) if isinstance(v, (list, dict)) else v
                for k, v in values.items()
            },
        }
        for bar, values in result.items()
        if isinstance(values, dict)
    ]
    pd.DataFrame(rows).to_csv(args.out_dir / "bars.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

# Purpose: compare methods on one v8 dataset at matched coverage, from score_v8.py output.
#          Each method is read at the operating point that covers at least the target
#          coverage with the highest accuracy among calls; a method that never reaches
#          the target gets NA, never an extrapolated value.
# Run:     .venv/bin/python benchmarks/scripts/compare_v8.py OUT_DIR SCORE.json [SCORE.json ...]
# Inputs:  score_v8.py JSON files for one dataset, one per method [score_v8.py]
# Outputs:
#   OUT_DIR/matched_coverage.tsv   method, target coverage, coverage, descendant and strict accuracy
#   OUT_DIR/risk_coverage.png      descendant and strict curves, one line per method

import argparse
import json
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Parameters -----------------------------------------------------------------------
AGREEMENTS = ("descendant", "strict")  # the two agreement levels v6 published
TARGETS = (
    None  # None: use the smallest full-run coverage among the methods (the most selective one)
)


def curve(score: dict, agreement: str) -> pd.DataFrame:
    points = pd.DataFrame(score["risk_coverage"][agreement])
    points = points[points["risk"].notna()]
    return points.assign(accuracy=1 - points["risk"])[["coverage", "accuracy"]]


def at_coverage(score: dict, agreement: str, target: float) -> float:
    points = curve(score, agreement)
    reachable = points[points["coverage"] >= target - 1e-9]
    return float(reachable["accuracy"].max()) if len(reachable) else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("scores", type=Path, nargs="+")
    args = parser.parse_args()

    scores = {Path(path).stem: json.loads(Path(path).read_text()) for path in args.scores}
    datasets = {score["dataset"] for score in scores.values()}
    if len(datasets) != 1:
        raise ValueError(f"scores cover more than one dataset: {sorted(datasets)}")
    target = TARGETS or min(score["coverage"] for score in scores.values())

    # Table at the matched coverage ----------------------------------------------------
    rows = [
        {
            "method": score["system"],
            "target_coverage": round(target, 3),
            "coverage": round(score["coverage"], 3),
            **{f"{a}_accuracy_at_target": at_coverage(score, a, target) for a in AGREEMENTS},
            "strict_accuracy_full": score["accuracy"],
        }
        for score in scores.values()
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows)
    table.to_csv(args.out_dir / "matched_coverage.tsv", sep="\t", index=False)

    # Risk-coverage curves --------------------------------------------------------------
    fig, axes = plt.subplots(1, len(AGREEMENTS), figsize=(8, 3.4), sharey=True)
    for axis, agreement in zip(axes, AGREEMENTS, strict=True):
        for score in scores.values():
            points = curve(score, agreement).sort_values("coverage")
            axis.plot(
                points["coverage"], points["accuracy"], marker="o", ms=3, label=score["system"]
            )
        axis.axvline(target, color="grey", ls=":", lw=0.8)
        axis.set(title=f"{agreement} agreement", xlabel="coverage", xlim=(0, 1.05))
    axes[0].set_ylabel("accuracy among calls")
    axes[0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(args.out_dir / "risk_coverage.png", dpi=160)


if __name__ == "__main__":
    main()

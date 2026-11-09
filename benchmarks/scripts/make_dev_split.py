# Purpose: hold out 20% of the Human Immune Health Atlas donors as the development split
#          for the protocol v8 tusoskill iterations (decisions Q21, amendment 2).
# Run:     .venv/bin/python benchmarks/scripts/make_dev_split.py OBS_SUMMARY_HIHA OUT_DIR
# Inputs:  obs_summary_hiha.tsv   from acquire_v8.py (donor_id, sex, development_stage, disease)
# Outputs: OUT_DIR/donor_split.tsv   donor_id, stratum, split ("dev" or "corpus")
#          OUT_DIR/split_summary.tsv donors per stratum and split

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Parameters -----------------------------------------------------------------------
DEV_FRACTION = 0.2  # protocol v8, Q21
SEED = 1  # same seed as the cell cap, protocol v8, Q11
MIN_STRATUM = 5  # a stratum with fewer donors stays entirely in the corpus
AGE_BANDS = [0, 20, 45, 120]  # the atlas has donors at 11-12, 25-36 and 55-65 years
AGE_LABELS = ["child", "young_adult", "older_adult"]


def main():
    obs_path, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=False)  # a split is never overwritten

    # One row per donor ---------------------------------------------------------------
    donors = (
        pd.read_csv(obs_path, sep="\t")
        .groupby("donor_id")
        .agg(
            sex=("sex", "first"),
            stage=("development_stage", "first"),
            disease=("disease", "first"),
        )
        .reset_index()
    )
    donors["age"] = donors["stage"].str.extract(r"(\d+)-year", expand=False).astype(float)
    donors["stratum"] = (
        pd.cut(donors["age"], AGE_BANDS, labels=AGE_LABELS, right=False).astype(str)
        + "|"
        + donors["sex"]
        + "|"
        + donors["disease"].map({"normal": "CMV-", "cytomegalovirus infection": "CMV+"})
    )

    # Seeded draw inside each stratum -------------------------------------------------
    rng = np.random.default_rng(SEED)
    donors["split"] = "corpus"
    for _, rows in donors.groupby("stratum"):
        n_dev = round(DEV_FRACTION * len(rows)) if len(rows) >= MIN_STRATUM else 0
        chosen = rng.choice(sorted(rows["donor_id"]), size=n_dev, replace=False)
        donors.loc[donors["donor_id"].isin(chosen), "split"] = "dev"

    donors[["donor_id", "stratum", "split"]].to_csv(
        out_dir / "donor_split.tsv", sep="\t", index=False
    )
    (donors.groupby(["stratum", "split"]).size().rename("donors").reset_index()).to_csv(
        out_dir / "split_summary.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()

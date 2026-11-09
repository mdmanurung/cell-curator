# Purpose: fetch the protocol v8 data and write the pre-run tables (donor overlap,
#          donors per label, Tabula Sapiens subset counts). No corpus is learned here.
# Run:     sbatch benchmarks/scripts/acquire_v8.sbatch   (writes to a new RUN_DIR)
# Inputs:  CELLxGENE curation API, collections e9360edf-... (AIFI) and e5f58829-... (TS)
# Outputs: RUN_DIR/downloads/*.h5ad            the 8 AIFI files and the TS Immune file
#          RUN_DIR/downloads.tsv               dataset_id, file, bytes, sha256
#          RUN_DIR/aifi_donors_per_label.tsv   donors and cells per CL term
#          RUN_DIR/aifi_donor_attributes.tsv   one row per donor: age group, sex, CMV
#          RUN_DIR/ts_subset_counts.tsv        cells and donors per tissue and assay
#          RUN_DIR/donor_overlap.tsv           AIFI vs TS donor IDs shared
#          RUN_DIR/versions.txt

import hashlib
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

import anndata as ad
import pandas as pd

# Parameters -----------------------------------------------------------------------
AIFI_COLLECTION = "e9360edf-b0b7-4e01-bce8-e596814f13e7"  # protocol v8, Q2
TS_COLLECTION = "e5f58829-1a66-40b5-a624-9046778e74f5"
TS_DATASET = "c5d88abe-f23a-45fa-a534-788985e93dad"  # TS Immune, protocol v8, Q3
TS_TISSUES = ["blood", "spleen", "lymph node", "inguinal lymph node", "mesenteric lymph node"]
TS_ASSAYS = ["10x 3' v3", "10x 5' v2"]  # 10x only, protocol v8, Q3
API = "https://api.cellxgene.cziscience.com/curation/v1/collections/"


# Helpers (called 3+ times) --------------------------------------------------------
def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    run_dir = Path(sys.argv[1])
    run_dir.mkdir(parents=True, exist_ok=False)  # a run dir is never overwritten
    (run_dir / "downloads").mkdir()

    # Asset URLs --------------------------------------------------------------------
    aifi = json.load(urllib.request.urlopen(API + AIFI_COLLECTION, timeout=60))["datasets"]
    ts = json.load(urllib.request.urlopen(API + TS_COLLECTION, timeout=60))["datasets"]
    wanted = aifi + [d for d in ts if d["dataset_id"] == TS_DATASET]
    if len(wanted) != 9:
        raise ValueError(f"expected 8 AIFI + 1 TS datasets, found {len(wanted)}")

    # Download (resumable) and hash -------------------------------------------------
    rows = []
    for dataset in wanted:
        asset = next(a for a in dataset["assets"] if a["filetype"] == "H5AD")
        path = run_dir / "downloads" / f"{dataset['dataset_id']}.h5ad"
        subprocess.run(["curl", "-fL", "-C", "-", "-o", str(path), asset["url"]], check=True)
        if path.stat().st_size != asset["filesize"]:
            raise ValueError(f"{path.name}: size {path.stat().st_size} != {asset['filesize']}")
        rows.append((dataset["dataset_id"], path.name, path.stat().st_size, sha256(path)))
    pd.DataFrame(rows, columns=["dataset_id", "file", "bytes", "sha256"]).to_csv(
        run_dir / "downloads.tsv", sep="\t", index=False
    )

    # Read obs only (backed) --------------------------------------------------------
    aifi_obs = pd.concat(
        [
            ad.read_h5ad(
                run_dir / "downloads" / f"{d['dataset_id']}.h5ad", backed="r"
            ).obs.assign(dataset_id=d["dataset_id"])
            for d in aifi
        ]
    )
    ts_obs = ad.read_h5ad(run_dir / "downloads" / f"{TS_DATASET}.h5ad", backed="r").obs

    # AIFI donors per CL term (decides whether Q7 and Q11 leave enough donors) ------
    (
        aifi_obs.groupby("cell_type", observed=True)
        .agg(donors=("donor_id", "nunique"), cells=("donor_id", "size"))
        .reset_index()
        .to_csv(run_dir / "aifi_donors_per_label.tsv", sep="\t", index=False)
    )

    # One row per AIFI donor: the strata for the 20% dev split (protocol v8, Q21) ----
    (
        aifi_obs.groupby("donor_id", observed=True)
        .agg(
            dataset_id=("dataset_id", "first"),
            sex=("sex", "first"),
            development_stage=("development_stage", "first"),
            disease=("disease", "first"),
            cells=("donor_id", "size"),
        )
        .reset_index()
        .to_csv(run_dir / "aifi_donor_attributes.tsv", sep="\t", index=False)
    )

    # TS subset: what survives the tissue and assay filter ----------------------------
    ts_kept = ts_obs[ts_obs["tissue"].isin(TS_TISSUES) & ts_obs["assay"].isin(TS_ASSAYS)]
    (
        ts_kept.groupby(["tissue", "assay"], observed=True)
        .agg(cells=("donor_id", "size"), donors=("donor_id", "nunique"))
        .reset_index()
        .to_csv(run_dir / "ts_subset_counts.tsv", sep="\t", index=False)
    )

    # Donor overlap: IDs shared between the two studies (warns only, protocol v8, Q5) --
    shared = sorted(set(aifi_obs["donor_id"].astype(str)) & set(ts_obs["donor_id"].astype(str)))
    pd.DataFrame({"shared_donor_id": shared}).to_csv(
        run_dir / "donor_overlap.tsv", sep="\t", index=False
    )

    # Versions ----------------------------------------------------------------------
    (run_dir / "versions.txt").write_text(f"python {sys.version}\nanndata {ad.__version__}\n")


if __name__ == "__main__":
    main()

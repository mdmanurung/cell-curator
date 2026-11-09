# Purpose: fetch the protocol v8 data (amendment 1) and write the pre-run tables:
#          one long obs summary per source and the donor overlap between sources.
#          No corpus is learned here.
# Run:     sbatch benchmarks/scripts/acquire_v8.sbatch RUN_DIR   (RUN_DIR must not exist)
# Inputs:  CELLxGENE curation API, collections named in SOURCES below
# Outputs: RUN_DIR/downloads/<dataset_id>.h5ad   every dataset in SOURCES
#          RUN_DIR/downloads.tsv                 source, dataset_id, file, bytes, sha256
#          RUN_DIR/obs_summary_<source>.tsv      cells per donor, cell type, tissue, assay,
#                                                sex, development stage, disease
#          RUN_DIR/donor_overlap.tsv             donor IDs shared by each pair of sources
#          RUN_DIR/versions.txt

import hashlib
import itertools
import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

# Parameters -----------------------------------------------------------------------
API = "https://api.cellxgene.cziscience.com/curation/v1/collections/"
# source -> (collection id, which of its datasets), protocol v8 and amendment 1
SOURCES = {
    "sound_life": ("e9360edf-b0b7-4e01-bce8-e596814f13e7", lambda d: True),
    "hiha": (
        "77f9d7e9-5675-49c3-abed-ce02f39eef1b",
        lambda d: d["title"] == "Human Immune Health Atlas",
    ),
    "ts": (
        "e5f58829-1a66-40b5-a624-9046778e74f5",
        lambda d: d["dataset_id"] == "c5d88abe-f23a-45fa-a534-788985e93dad",
    ),
    "indonesia": (
        "d1e0e64d-6d2a-4a3e-b7f4-43ed909a9d9c",
        lambda d: d["title"].startswith("All cells"),
    ),
    "cross_tissue": (
        "62ef75e4-cbea-454e-a0ce-998ec40223d3",
        lambda d: d["title"] == "Global",
    ),
}
EXPECTED_DATASETS = {"sound_life": 8, "hiha": 1, "ts": 1, "indonesia": 1, "cross_tissue": 1}
OBS_COLUMNS = [
    "donor_id",
    "cell_type",
    "tissue",
    "assay",
    "sex",
    "development_stage",
    "disease",
]


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

    # Pick the datasets and check how many each source gives -------------------------
    chosen = {}
    for source, (collection, keep) in SOURCES.items():
        # curl, because the uv Python on compute nodes has no CA bundle for urllib
        reply = subprocess.run(
            ["curl", "-fsSL", "--max-time", "60", API + collection],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        datasets = json.loads(reply)["datasets"]
        chosen[source] = [d for d in datasets if keep(d)]
        if len(chosen[source]) != EXPECTED_DATASETS[source]:
            raise ValueError(f"{source}: expected {EXPECTED_DATASETS[source]} datasets")

    # Download (resumable), check size, hash ------------------------------------------
    rows = []
    for source, datasets in chosen.items():
        for dataset in datasets:
            asset = next(a for a in dataset["assets"] if a["filetype"] == "H5AD")
            path = run_dir / "downloads" / f"{dataset['dataset_id']}.h5ad"
            subprocess.run(
                ["curl", "-fL", "-C", "-", "-o", str(path), asset["url"]], check=True
            )
            if path.stat().st_size != asset["filesize"]:
                raise ValueError(
                    f"{path.name}: size {path.stat().st_size} != {asset['filesize']}"
                )
            rows.append(
                (source, dataset["dataset_id"], path.name, path.stat().st_size, sha256(path))
            )
    pd.DataFrame(rows, columns=["source", "dataset_id", "file", "bytes", "sha256"]).to_csv(
        run_dir / "downloads.tsv", sep="\t", index=False
    )

    # One long obs summary per source, read from obs only (backed) --------------------
    donors = {}
    for source, datasets in chosen.items():
        # sex and age labels differ between studies, so keep them as stored
        summary = pd.concat(
            [
                ad.read_h5ad(run_dir / "downloads" / f"{d['dataset_id']}.h5ad", backed="r")
                .obs[OBS_COLUMNS]
                .astype(str)
                .value_counts()
                .rename("cells")
                .reset_index()
                .assign(dataset_id=d["dataset_id"])
                for d in datasets
            ]
        )
        summary.to_csv(run_dir / f"obs_summary_{source}.tsv", sep="\t", index=False)
        donors[source] = set(summary["donor_id"])

    # Donor overlap between every pair of sources (warns only for novel sets, Q5) -------
    pd.DataFrame(
        [
            (a, b, len(donors[a] & donors[b]), ";".join(sorted(donors[a] & donors[b])[:20]))
            for a, b in itertools.combinations(donors, 2)
        ],
        columns=["source_a", "source_b", "shared_donor_ids", "first_20_ids"],
    ).to_csv(run_dir / "donor_overlap.tsv", sep="\t", index=False)

    # Versions -----------------------------------------------------------------------
    (run_dir / "versions.txt").write_text(f"python {sys.version}\nanndata {ad.__version__}\n")


if __name__ == "__main__":
    main()

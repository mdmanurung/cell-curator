# Purpose: learn the protocol v8 corpus proposal from the Human Immune Health Atlas:
#          cell-level Wilcoxon on a donor-capped sample, top positives and sibling-based
#          negatives per Cell Ontology term. The output is a proposal; no record is
#          reviewed or callable until the user reviews it in chat (AGENTS.md, gate 3).
# Run:     sbatch benchmarks/scripts/learn_corpus_v8.sbatch RUN_DIR   (RUN_DIR must not exist)
# Inputs:  the atlas h5ad            from acquire_v8.py (dataset e522d2cd-...)
#          donor_split.tsv           benchmarks/data/aifi-hiha-dev-split/
# Outputs: RUN_DIR/sampled_cells.tsv   one row per sampled cell: barcode, donor, label
#          RUN_DIR/ranked.tsv.gz       the full Wilcoxon table, every label and gene
#          RUN_DIR/markers.tsv         one row per label and gene: role, score, log fold
#                                      change, FDR, detection inside and outside the label
#          RUN_DIR/siblings.tsv        each label's AIFI_L1 group and its two nearest siblings
#          RUN_DIR/corpus_records.json records for `cell-curator corpus import --source-kind generic`
#          RUN_DIR/corpus_donors.txt   the corpus donors, for derived_from
#          RUN_DIR/versions.txt

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]

# Parameters -----------------------------------------------------------------------
H5AD = Path(
    "/exports/para-lipg-hpc/mdmanurung/cc-v8-work/acquire-2026-11-09/downloads/"
    "e522d2cd-7927-4e59-a4ed-064009569279.h5ad"
)
DONOR_SPLIT = ROOT / "benchmarks/data/aifi-hiha-dev-split/donor_split.tsv"
CELLS_PER_DONOR_LABEL = 100  # protocol v8, Q11
SEED = 1  # protocol v8, Q11
N_POSITIVE = 10  # protocol v8, Q12
N_NEGATIVE = 5  # protocol v8, Q12
FDR = 0.05  # BH over genes within each label, protocol v8, Q12 and the FDR-family row
N_SIBLINGS = 2  # v5 rule: the two nearest siblings in the same broad group
BARRED_GENES = r"^(RPL|RPS|MT-)"  # ribosomal and mitochondrial genes, amendment 5
# housekeeping genes barred as markers (user, 2026-11-09, amendment 6)
BARRED_LIST = [
    "ACTB", "ACTG1", "B2M", "CFL1", "EEF1A1", "EEF1B2", "EEF1D", "EEF1G", "EEF2", "FTH1",
    "FTL", "GAPDH", "H2AZ2", "H3-3A", "HMGB1", "HMGB2", "HMGN2", "HNRNPA1", "HSP90AB1",
    "MYL6", "NPM1", "OAZ1", "PFN1", "PPIA", "RACK1", "TAGLN2", "TMSB10", "TMSB4X", "TPT1",
    "UBA52", "UBB", "UBC",
]  # fmt: skip
# atlas CL terms that hold one AIFI L3 state; mapping approved by the user (amendment 6)
STATE_RECORDS = {
    "T cell": ("cycling T cell", "CL:4033069", "exact"),  # AIFI L3 Proliferating T cell
    "CD8-positive, alpha-beta memory T cell": (
        "ISG+ memory CD8 T cell",
        "CL:0000909",
        "narrower",
    ),  # AIFI L3 ISG+ memory CD8 T cell
}
MAX_OWN_DETECTION = 0.2  # a negative is detected in under 20% of the label, amendment 5
CONTEXT = {
    "species": "human",
    "tissue": "blood",
    "assay": "10x 3' v3",
}  # stage, condition: review


def main():
    run_dir = Path(sys.argv[1])
    run_dir.mkdir(parents=True, exist_ok=False)  # a run dir is never overwritten

    # Corpus donors only: the dev donors are never read --------------------------------
    split = pd.read_csv(DONOR_SPLIT, sep="\t")
    corpus_donors = set(split.loc[split["split"] == "corpus", "donor_id"])
    (run_dir / "corpus_donors.txt").write_text("\n".join(sorted(corpus_donors)) + "\n")
    atlas = ad.read_h5ad(H5AD, backed="r")
    # the integer index is each cell's row in the file, so sorting it keeps the file order
    obs = atlas.obs[["donor_id", "cell_type", "cell_type_ontology_term_id", "AIFI_L1"]]
    obs = obs.astype(str).rename_axis("barcode").reset_index()
    obs = obs[obs["donor_id"].isin(corpus_donors)]

    # At most CELLS_PER_DONOR_LABEL cells per donor and label ---------------------------
    rng = np.random.default_rng(SEED)
    keep = []
    for _, rows in obs.groupby(["donor_id", "cell_type"], sort=True):
        keep.append(
            rng.choice(
                rows.index.to_numpy(), min(len(rows), CELLS_PER_DONOR_LABEL), replace=False
            )
        )
    sampled = obs.loc[np.sort(np.concatenate(keep))]
    sampled.to_csv(run_dir / "sampled_cells.tsv", sep="\t", index=False)

    # Raw counts of the sample, then log1p of counts per 10,000 over all genes ----------
    data = ad.AnnData(
        X=atlas.raw.X[sampled.index.to_numpy()],
        obs=sampled.set_index("barcode"),
        var=pd.DataFrame(
            index=atlas.var_names, data={"symbol": atlas.var["feature_name"].astype(str).values}
        ),
    )
    # a symbol shared by several genes stays with the gene that has the most counts
    totals = pd.Series(np.asarray(data.X.sum(0)).ravel(), index=data.var_names)
    order = totals.sort_values(ascending=False).index
    rank = pd.Series(range(len(order)), index=order)[data.var_names]
    first = (
        data.var.assign(rank=rank.values).groupby("symbol")["rank"].transform("min")
        == rank.values
    )
    data.var["gene"] = np.where(
        first, data.var["symbol"], data.var["symbol"] + "--" + data.var_names
    )
    data.var_names = data.var["gene"].values
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)

    # Wilcoxon one-versus-rest per label (BH inside each label) -----------------------
    sc.tl.rank_genes_groups(data, "cell_type", method="wilcoxon", pts=True, n_genes=data.n_vars)
    labels = sorted(data.obs["cell_type"].unique())
    ranked = pd.concat(
        [
            sc.get.rank_genes_groups_df(data, group=label).assign(cell_type=label)
            for label in labels
        ]
    ).rename(columns={"names": "gene"})
    pts_in = data.uns["rank_genes_groups"]["pts"]
    pts_out = data.uns["rank_genes_groups"]["pts_rest"]
    ranked["pct_in"] = [
        pts_in.loc[g, c] for g, c in zip(ranked["gene"], ranked["cell_type"], strict=True)
    ]
    ranked["pct_out"] = [
        pts_out.loc[g, c] for g, c in zip(ranked["gene"], ranked["cell_type"], strict=True)
    ]

    ranked.to_csv(run_dir / "ranked.tsv.gz", sep="\t", index=False)

    # Positives: top scores that are enriched and significant, barred genes excluded -------
    positives = {
        label: ranked[
            (ranked["cell_type"] == label)
            & (ranked["pvals_adj"] < FDR)
            & (ranked["logfoldchanges"] > 0)
            & ~ranked["gene"].str.match(BARRED_GENES)
            & ~ranked["gene"].isin(BARRED_LIST)
        ]
        .sort_values("scores", ascending=False)["gene"]
        .head(N_POSITIVE)
        .tolist()
        for label in labels
    }

    # Siblings: same AIFI_L1 group, the two with the most similar mean profile --------------
    union = sorted({g for genes in positives.values() for g in genes})
    means = pd.DataFrame(
        {
            label: np.asarray(data[data.obs["cell_type"] == label, union].X.mean(0)).ravel()
            for label in labels
        },
        index=union,
    )
    similarity = means.corr()
    group = sampled.groupby("cell_type")["AIFI_L1"].agg(lambda s: s.value_counts().index[0])
    sibling_rows, negatives = [], {}
    for label in labels:
        peers = [p for p in labels if p != label and group[p] == group[label]]
        # a label alone in its group takes the nearest labels anywhere (amendment 5)
        peers = peers or [p for p in labels if p != label]
        nearest = (
            similarity.loc[label, peers]
            .sort_values(ascending=False)
            .head(N_SIBLINGS)
            .index.tolist()
        )
        # sibling positives, once each, absent from the label, rare inside the label
        candidates = dict.fromkeys(g for sib in nearest for g in positives[sib])
        negatives[label] = [
            g
            for g in candidates
            if g not in positives[label] and pts_in.loc[g, label] < MAX_OWN_DETECTION
        ][:N_NEGATIVE]
        sibling_rows.append((label, group[label], ";".join(nearest)))
    pd.DataFrame(sibling_rows, columns=["cell_type", "aifi_l1", "nearest_siblings"]).to_csv(
        run_dir / "siblings.tsv", sep="\t", index=False
    )

    # Tables and records --------------------------------------------------------------
    roles = pd.concat(
        [
            ranked[(ranked["cell_type"] == label) & ranked["gene"].isin(genes)].assign(
                role=role
            )
            for label in labels
            for role, genes in (("positive", positives[label]), ("negative", negatives[label]))
        ]
    )
    # per-donor support: share of donors whose label cells agree with the marker's role
    support = []
    for label, gene, role, out in roles[["cell_type", "gene", "role", "pct_out"]].itertuples(
        index=False
    ):
        cells = data[data.obs["cell_type"] == label, gene]
        per_donor = (
            pd.Series(np.asarray(cells.X.toarray()).ravel() > 0, index=cells.obs["donor_id"])
            .groupby(level=0)
            .mean()
        )
        agrees = per_donor > out if role == "positive" else per_donor < MAX_OWN_DETECTION
        support.append(agrees.mean())
    roles["donor_support"] = support
    roles.to_csv(run_dir / "markers.tsv", sep="\t", index=False)
    ids = sampled.groupby("cell_type")["cell_type_ontology_term_id"].first()
    counts = sampled.groupby("cell_type").agg(
        cells=("donor_id", "size"), donors=("donor_id", "nunique")
    )
    records = [
        {
            "source_record_id": ids[name].replace(":", "_"),
            "name": STATE_RECORDS.get(name, (name,))[0],
            "source_label": name,
            "description": (
                f"Human Immune Health Atlas cell type {name} ({ids[name]}). Positives: top {N_POSITIVE} "
                f"genes by Wilcoxon score at FDR {FDR} from {counts.loc[name, 'cells']} cells of "
                f"{counts.loc[name, 'donors']} corpus donors. Negatives: positives of the nearest "
                f"siblings in {group[name]} that are not this label's positives."
            ),
            "positive": positives[name],
            "negative": negatives[name],
            "cl_id": STATE_RECORDS[name][1] if name in STATE_RECORDS else ids[name],
            # the ontology check runs at export
            "ontology_relation": STATE_RECORDS[name][2] if name in STATE_RECORDS else "exact",
            "canonical_level": "L1",
            "source_level": "CL term",
        }
        for name in labels
    ]
    (run_dir / "corpus_records.json").write_text(
        json.dumps({"context": CONTEXT, "labels": records}, indent=2, sort_keys=True)
    )
    (run_dir / "versions.txt").write_text(
        f"scanpy {sc.__version__}\nanndata {ad.__version__}\n"
    )


if __name__ == "__main__":
    main()

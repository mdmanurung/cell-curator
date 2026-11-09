# Purpose: learn the protocol v8 corpus proposal from the Human Immune Health Atlas:
#          cell-level Wilcoxon on a donor-capped sample, top positives and sibling-based
#          negatives per Cell Ontology term. The output is a proposal; no record is
#          reviewed or callable until the user reviews it in chat (AGENTS.md, gate 3).
# Run:     sbatch benchmarks/scripts/learn_corpus_v8.sbatch RUN_DIR [--candidate CAND.json]
#              [--reuse-ranked PREVIOUS_RUN_DIR]          (RUN_DIR must not exist)
#          With neither option the output is the frozen learn-2026-11-09c selection.
#          --reuse-ranked skips the Wilcoxon and reads PREVIOUS_RUN_DIR/ranked.tsv.gz; it needs
#          the same sample (checked against PREVIOUS_RUN_DIR/sampled_cells.tsv).
#          --candidate CAND.json sets selection parameters for a development candidate; keys in
#          SELECTION_KEYS below. A candidate corpus is never the frozen corpus.
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
#          RUN_DIR/params.json         the selection parameters this run used
#          RUN_DIR/ranked_vs_siblings.tsv.gz   only when ranking_contrast is vs_siblings

import argparse
import json
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
# Selection parameters a candidate JSON may set. The defaults are the frozen v8 values above.
SELECTION_KEYS = {
    "n_positive": N_POSITIVE,
    "n_negative": N_NEGATIVE,
    "fdr": FDR,
    "max_own_detection": MAX_OWN_DETECTION,
    "n_siblings": N_SIBLINGS,
    "barred_add": [],  # extra genes barred from positives and negatives
    "barred_pattern_add": "",  # extra regex, joined to BARRED_GENES with |
    "min_logfc": 0.0,  # positives need a log fold change above this (0 is the frozen rule)
    "min_pct_in": 0.0,  # positives are detected in at least this share of the label's cells
    "ranking_contrast": "one_vs_rest",  # or "vs_siblings": positives ranked against siblings
    "negative_source": "sibling_positives",  # or "sibling_contrast" (needs vs_siblings)
}
CONTEXT = {
    "species": "human",
    "tissue": "blood",
    "assay": "10x 3' v3",
}  # stage, condition: review


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--reuse-ranked", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir
    run_dir.mkdir(parents=True, exist_ok=False)  # a run dir is never overwritten

    # Selection parameters: the frozen defaults, then the candidate's changes ---------------
    candidate = json.loads(args.candidate.read_text()) if args.candidate else {}
    if set(candidate) - set(SELECTION_KEYS) - {"name", "note", "engine"}:
        raise ValueError(f"unknown candidate keys: {set(candidate) - set(SELECTION_KEYS)}")
    params = {**SELECTION_KEYS, **candidate}
    if params["ranking_contrast"] not in {"one_vs_rest", "vs_siblings"}:
        raise ValueError("ranking_contrast is one_vs_rest or vs_siblings")
    if params["negative_source"] not in {"sibling_positives", "sibling_contrast"}:
        raise ValueError("negative_source is sibling_positives or sibling_contrast")
    if (
        params["negative_source"] == "sibling_contrast"
        and params["ranking_contrast"] != "vs_siblings"
    ):
        raise ValueError("sibling_contrast negatives need ranking_contrast vs_siblings")
    n_positive, n_negative, fdr = params["n_positive"], params["n_negative"], params["fdr"]
    min_logfc, min_pct_in = params["min_logfc"], params["min_pct_in"]
    max_own, n_siblings = params["max_own_detection"], params["n_siblings"]
    barred_pattern = "|".join(
        [BARRED_GENES]
        + ([f"(?:{params['barred_pattern_add']})"] if params["barred_pattern_add"] else [])
    )
    barred_list = BARRED_LIST + list(params["barred_add"])
    (run_dir / "params.json").write_text(json.dumps(params, indent=2, sort_keys=True))

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

    # Wilcoxon one-versus-rest per label (BH inside each label), or the saved table -------
    labels = sorted(data.obs["cell_type"].unique())
    if args.reuse_ranked:
        previous = pd.read_csv(args.reuse_ranked / "sampled_cells.tsv", sep="\t", dtype=str)
        if previous["barcode"].tolist() != sampled["barcode"].tolist():
            raise ValueError("the sample differs from the one behind --reuse-ranked")
        ranked = pd.read_csv(
            args.reuse_ranked / "ranked.tsv.gz",
            sep="\t",
            keep_default_na=False,
            na_values=[""],
            float_precision="round_trip",
        )
        (run_dir / "ranked_from.txt").write_text(f"{args.reuse_ranked}\n")
    else:
        sc.tl.rank_genes_groups(
            data, "cell_type", method="wilcoxon", pts=True, n_genes=data.n_vars
        )
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
    pts_in = ranked.pivot(index="gene", columns="cell_type", values="pct_in")

    # Positives (one-versus-rest): top scores that are enriched and significant ------------
    def pick_positives(table, label):
        return (
            table[
                (table["pvals_adj"] < fdr)
                & (table["logfoldchanges"] > min_logfc)
                & (pts_in.loc[table["gene"], label].to_numpy() >= min_pct_in)
                & ~table["gene"].str.match(barred_pattern)
                & ~table["gene"].isin(barred_list)
            ]
            .sort_values("scores", ascending=False)["gene"]
            .head(n_positive)
            .tolist()
        )

    by_label = {label: rows for label, rows in ranked.groupby("cell_type")}
    positives = {label: pick_positives(by_label[label], label) for label in labels}

    # Siblings: same AIFI_L1 group, the nearest by mean profile over the positives ----------
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
    nearest = {}
    for label in labels:
        peers = [p for p in labels if p != label and group[p] == group[label]]
        # a label alone in its group takes the nearest labels anywhere (amendment 5)
        peers = peers or [p for p in labels if p != label]
        nearest[label] = (
            similarity.loc[label, peers]
            .sort_values(ascending=False)
            .head(n_siblings)
            .index.tolist()
        )
    pd.DataFrame(
        [(label, group[label], ";".join(nearest[label])) for label in labels],
        columns=["cell_type", "aifi_l1", "nearest_siblings"],
    ).to_csv(run_dir / "siblings.tsv", sep="\t", index=False)

    # Candidate contrast: positives ranked against the nearest siblings only ----------------
    contrast = {}
    if params["ranking_contrast"] == "vs_siblings":
        for label in labels:
            sub = data[data.obs["cell_type"].isin([label, *nearest[label]])].copy()
            sub.obs["arm"] = np.where(sub.obs["cell_type"] == label, "label", "siblings")
            sc.tl.rank_genes_groups(
                sub,
                "arm",
                groups=["label"],
                reference="siblings",
                method="wilcoxon",
                n_genes=sub.n_vars,
            )
            contrast[label] = sc.get.rank_genes_groups_df(sub, group="label").rename(
                columns={"names": "gene"}
            )
        pd.concat([t.assign(cell_type=label) for label, t in contrast.items()]).to_csv(
            run_dir / "ranked_vs_siblings.tsv.gz", sep="\t", index=False
        )
        positives = {label: pick_positives(contrast[label], label) for label in labels}

    # Negatives: absent from the label, rare inside it, at most n_negative -----------------
    negatives = {}
    for label in labels:
        if params["negative_source"] == "sibling_contrast":
            # genes the siblings express above the label, strongest first
            table = contrast[label]
            candidates = table[
                (table["pvals_adj"] < fdr)
                & (table["logfoldchanges"] < 0)
                & ~table["gene"].str.match(barred_pattern)
                & ~table["gene"].isin(barred_list)
            ].sort_values("scores")["gene"]
        else:
            # sibling positives, once each
            candidates = pd.Series(
                list(dict.fromkeys(g for sib in nearest[label] for g in positives[sib]))
            )
        negatives[label] = [
            g
            for g in candidates
            if g not in positives[label] and pts_in.loc[g, label] < max_own
        ][:n_negative]

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
        agrees = per_donor > out if role == "positive" else per_donor < max_own
        support.append(agrees.mean())
    roles["donor_support"] = support
    roles.to_csv(run_dir / "markers.tsv", sep="\t", index=False)
    ids = sampled.groupby("cell_type")["cell_type_ontology_term_id"].first()
    counts = sampled.groupby("cell_type").agg(
        cells=("donor_id", "size"), donors=("donor_id", "nunique")
    )
    contrast_text = " against the nearest siblings" if contrast else ""
    negative_text = (
        "genes the nearest siblings express above this label"
        if params["negative_source"] == "sibling_contrast"
        else "positives of the nearest siblings"
    )
    records = [
        {
            "source_record_id": ids[name].replace(":", "_"),
            "name": STATE_RECORDS.get(name, (name,))[0],
            "source_label": name,
            "description": (
                f"Human Immune Health Atlas cell type {name} ({ids[name]}). Positives: top {n_positive} "
                f"genes by Wilcoxon score{contrast_text} at FDR {fdr} from {counts.loc[name, 'cells']} "
                f"cells of {counts.loc[name, 'donors']} corpus donors. Negatives: {negative_text} "
                f"in {group[name]} that are not this label's positives."
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

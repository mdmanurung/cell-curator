"""Map researcher free-text labels to Cell Ontology ids and measure agreement.

Purpose: protocol v8 consistency measure for the three independent researcher agents.
         The rule is in mapping_rule.md and is mechanical: no judgement, no fuzzy match.
Run:     .venv/bin/python benchmarks/v8-researcher/map_and_agree.py CL_OBO RESPONSES_DIR OUT_DIR
         OUT_DIR must not exist
Inputs:  CL_OBO         pinned Cell Ontology [benchmarks/data/aida-phase1-v2-pbmc-hvg2000/cl.obo]
         RESPONSES_DIR/<dataset>/researcher_<1..3>.tsv
                        columns: cluster, label, cl_id, confidence [the researchers, prompt.md]
Outputs: OUT_DIR/mapped.tsv     one row per dataset, researcher, cluster: label, stated id,
                                mapped id, status, coarse-lineage id
         OUT_DIR/agreement.tsv  dataset, level, measure, value, n_clusters
         OUT_DIR/unmapped.tsv   labels that did not map, with counts (report these)
"""

import itertools
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from cell_curator.contracts import ContractError
from cell_curator.ontology import CellOntology

# Parameters -----------------------------------------------------------------------
RESEARCHERS = (1, 2, 3)
UNRESOLVED = "UNRESOLVED"
# Labels that mean "no call". Compared after the normalisation below.
NO_CALL = {"unresolved", "unknown", "unassigned", "na", "n/a", ""}
# Coarse lineage terms for the parent-level measure (proposed, see mapping_rule.md).
# A CL id maps to its deepest ancestor-or-self in this list; if none, it stays itself.
COARSE_LINEAGE = {
    "CL:0000084": "T cell",
    "CL:0000236": "B cell",
    "CL:0000623": "natural killer cell",
    "CL:0001065": "innate lymphoid cell",
    "CL:0000786": "plasma cell",
    "CL:0000576": "monocyte",
    "CL:0000235": "macrophage",
    "CL:0000451": "dendritic cell",
    "CL:0000094": "granulocyte",
    "CL:0000233": "platelet",
    "CL:0000232": "erythrocyte",
    "CL:0000037": "hematopoietic stem cell",
    "CL:0000066": "epithelial cell",
    "CL:0000115": "endothelial cell",
    "CL:0000057": "fibroblast",
}


# Mapping rule -----------------------------------------------------------------------


def normalise(text):
    """Trim, collapse spaces, drop a trailing period; case is handled by the ontology."""
    return re.sub(r"\s+", " ", str(text)).strip().rstrip(".")


def map_stated(ontology, stated, text):
    """User rule 2026-11-09: a valid, current stated CL id wins; else map the free text."""
    if stated.startswith("CL:"):
        try:
            term = ontology.resolve(stated)
        except ContractError:
            term = None
        if term is not None and term.accession == stated and not term.obsolete:
            return stated, "stated_id"
    return map_label(ontology, text)


def map_label(ontology, text):
    """Return (cl_id, status). Status is one of exact, no_call, unmatched, ambiguous, obsolete."""
    label = normalise(text)
    if label.casefold() in NO_CALL:
        return UNRESOLVED, "no_call"
    # Plurals ("NK cells") get the singular tried second, as comparator_to_neutral.py does.
    candidates = [label, label[:-1]] if label.casefold().endswith(" cells") else [label]
    for candidate in candidates:
        try:
            term = ontology.resolve(candidate)
        except ContractError as err:
            if "ambiguous" in str(err):
                return UNRESOLVED, "ambiguous"
            continue
        return (UNRESOLVED, "obsolete") if term.obsolete else (term.accession, "exact")
    return UNRESOLVED, "unmatched"


def coarse_id(ontology, cl_id):
    if cl_id == UNRESOLVED:
        return UNRESOLVED
    lineage = [a for a in ontology.ancestors(cl_id) | {cl_id} if a in COARSE_LINEAGE]
    return max(lineage, key=ontology.depth) if lineage else cl_id


# Agreement ------------------------------------------------------------------------------


def fleiss_kappa(ratings):
    """ratings: clusters x raters array of category labels. NaN when chance agreement is 1."""
    categories = np.unique(ratings)
    counts = np.stack([(ratings == c).sum(axis=1) for c in categories], axis=1)
    n_raters = ratings.shape[1]
    p_item = ((counts**2).sum(axis=1) - n_raters) / (n_raters * (n_raters - 1))
    p_cat = counts.sum(axis=0) / counts.sum()
    p_chance = (p_cat**2).sum()
    return np.nan if p_chance == 1 else (p_item.mean() - p_chance) / (1 - p_chance)


def pairwise_agreement(ratings, same):
    pairs = itertools.combinations(range(ratings.shape[1]), 2)
    return float(np.mean([same(ratings[:, i], ratings[:, j]).mean() for i, j in pairs]))


def agreement_rows(ontology, wide, dataset):
    """wide: one row per cluster, columns r1..r3 holding CL ids or UNRESOLVED."""
    rows = []
    exact = wide[[f"r{k}" for k in RESEARCHERS]].to_numpy()
    coarse = wide[[f"c{k}" for k in RESEARCHERS]].to_numpy()
    all_resolved = (exact != UNRESOLVED).all(axis=1)

    def related(a, b):  # one id is the other, or its ancestor
        return np.array(
            [
                x == y
                or (
                    UNRESOLVED not in (x, y)
                    and (x in ontology.ancestors(y) or y in ontology.ancestors(x))
                )
                for x, y in zip(a, b, strict=True)
            ]
        )

    def add(level, measure, value, n):
        rows.append((dataset, level, measure, value, n))

    same = lambda a, b: a == b  # noqa: E731
    add("exact", "pairwise_agreement", pairwise_agreement(exact, same), len(wide))
    add("exact", "fleiss_kappa", fleiss_kappa(exact), len(wide))
    add("coarse_lineage", "pairwise_agreement", pairwise_agreement(coarse, same), len(wide))
    add("coarse_lineage", "fleiss_kappa", fleiss_kappa(coarse), len(wide))
    add(
        "ancestor_or_equal", "pairwise_agreement", pairwise_agreement(exact, related), len(wide)
    )
    # Clusters all three researchers resolved: agreement among real calls only.
    if all_resolved.any():
        add(
            "exact_all_resolved",
            "pairwise_agreement",
            pairwise_agreement(exact[all_resolved], same),
            int(all_resolved.sum()),
        )
        add(
            "exact_all_resolved",
            "fleiss_kappa",
            fleiss_kappa(exact[all_resolved]),
            int(all_resolved.sum()),
        )
        add(
            "coarse_all_resolved",
            "pairwise_agreement",
            pairwise_agreement(coarse[all_resolved], same),
            int(all_resolved.sum()),
        )
        add(
            "coarse_all_resolved",
            "fleiss_kappa",
            fleiss_kappa(coarse[all_resolved]),
            int(all_resolved.sum()),
        )
    add("exact", "fraction_unresolved", float((exact == UNRESOLVED).mean()), len(wide))
    return rows


# Run ------------------------------------------------------------------------------------


def main():
    obo, responses, out_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    out_dir.mkdir(parents=True)

    ontology = CellOntology.load(obo)
    for term in COARSE_LINEAGE:
        ontology.resolve(term)  # fails if a coarse term is not in the pinned ontology

    tables = []
    for path in sorted(responses.glob("*/researcher_*.tsv")):
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        table.insert(0, "researcher", int(path.stem.split("_")[1]))
        table.insert(0, "dataset", path.parent.name)
        tables.append(table)
    mapped = pd.concat(tables, ignore_index=True)
    mapped[["mapped_cl_id", "status"]] = [
        map_stated(ontology, s, t)
        for s, t in zip(mapped["cl_id"], mapped["label"], strict=True)
    ]
    mapped["coarse_cl_id"] = [coarse_id(ontology, c) for c in mapped["mapped_cl_id"]]
    # A valid stated id is used; the mismatch with the text mapping is still reported.
    mapped = mapped.rename(columns={"cl_id": "stated_cl_id"})

    rows = []
    for dataset, part in mapped.groupby("dataset"):
        # validate: each researcher must label each cluster exactly once
        wide = part.pivot(
            index="cluster", columns="researcher", values=["mapped_cl_id", "coarse_cl_id"]
        )
        if wide.isna().any().any() or wide.shape[1] != 2 * len(RESEARCHERS):
            raise ValueError(f"{dataset}: every researcher must label every cluster once")
        wide.columns = [
            f"{'r' if kind == 'mapped_cl_id' else 'c'}{k}" for kind, k in wide.columns
        ]
        rows += agreement_rows(ontology, wide, dataset)
    # Stated id against mapped id, over calls that gave an id: reported, never corrected.
    stated = mapped[mapped["stated_cl_id"].ne("")]
    for dataset, part in stated.groupby("dataset"):
        mismatch = (part["stated_cl_id"] != part["mapped_cl_id"]).mean()
        rows.append(
            (dataset, "stated_id", "fraction_mismatch_with_mapped", mismatch, len(part))
        )
    agreement = pd.DataFrame(
        rows, columns=["dataset", "level", "measure", "value", "n_clusters"]
    )

    mapped.to_csv(out_dir / "mapped.tsv", sep="\t", index=False)
    agreement.to_csv(out_dir / "agreement.tsv", sep="\t", index=False)
    unmapped = mapped[mapped["status"].isin(["unmatched", "ambiguous", "obsolete"])]
    unmapped.groupby(["dataset", "status", "label"]).size().rename("n").reset_index().to_csv(
        out_dir / "unmapped.tsv", sep="\t", index=False
    )
    print(agreement.to_string(index=False))


if __name__ == "__main__":
    main()

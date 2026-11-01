"""Build the protocol v5 PBMC vocabulary from the published Azimuth level-2 table.

Purpose: replace the self-authored v3 vocabulary with one whose positives, labels, and
         Cell Ontology terms all come from Azimuth (Hao et al., Cell 2021), and whose
         negatives come from one mechanical rule fixed before any v5 run.
Run:     uv run python benchmarks/data/aida-pbmc-protocol-v5/build_vocabulary.py
Inputs:
    azimuth_references_2026-11-02.html   https://azimuth.hubmapconsortium.org/references/
                                         saved 2026-11-02 (sha256 in VOCABULARY_PROVENANCE.md)
Outputs:
    azimuth_pbmc_l2.tsv     label, expanded_label, cl_id, cl_name, markers, as published
    pbmc_markers.json       the vocabulary in the shipped marker-file format

Negative rule: each type's L1 group is Azimuth's level-1 group (fixed below). Its two
nearest siblings in that group are the two with the highest Jaccard overlap of positive
markers, ties broken by table order. Negatives are the siblings' positives that are not
the type's own positives, in published order, first sibling first, at most 6.
"""

import csv
import html
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAGE = HERE / "azimuth_references_2026-11-02.html"
MAX_NEGATIVES = 6
N_SIBLINGS = 2
# Azimuth PBMC level-1 groups (celltype.l1), as defined by the reference.
L1_GROUP = {
    "B intermediate": "B", "B memory": "B", "B naive": "B", "Plasmablast": "B",
    "CD4 CTL": "CD4 T", "CD4 Naive": "CD4 T", "CD4 Proliferating": "CD4 T",
    "CD4 TCM": "CD4 T", "CD4 TEM": "CD4 T", "Treg": "CD4 T",
    "CD8 Naive": "CD8 T", "CD8 Proliferating": "CD8 T", "CD8 TCM": "CD8 T",
    "CD8 TEM": "CD8 T",
    "dnT": "other T", "gdT": "other T", "MAIT": "other T",
    "NK": "NK", "NK Proliferating": "NK", "NK_CD56bright": "NK",
    "CD14 Mono": "Mono", "CD16 Mono": "Mono",
    "ASDC": "DC", "cDC1": "DC", "cDC2": "DC", "pDC": "DC",
    "Eryth": "other", "HSPC": "other", "ILC": "other", "Platelet": "other",
}  # fmt: skip

if (HERE / "pbmc_markers.json").exists():
    raise ValueError("pbmc_markers.json exists; the v5 vocabulary is frozen")

# Extract the published level-2 table --------------------------------------------------

page = PAGE.read_text()
pbmc = page[
    [m.start() for m in re.finditer("Human - PBMC", page)][1] : [
        m.start() for m in re.finditer("Human - Motor Cortex", page)
    ][1]
]
section = next(part for part in pbmc.split("<summary>") if part.startswith("celltype.l2"))
table = re.search(r"<table>(.*?)</table>", section).group(1)


def text(cell: str) -> str:
    return html.unescape(re.sub("<[^>]+>", "", cell)).strip()


rows = []
for row in re.findall(r"<tr>(.*?)</tr>", table)[1:]:
    cells = re.findall(r"<td[^>]*>(.*?)</td>", row)
    rows.append(
        {
            "label": text(cells[0]),
            "expanded_label": text(cells[1]),
            "cl_id": re.search(r"obo/(CL_\d+)", cells[2]).group(1).replace("_", ":"),
            "cl_name": text(cells[2]),
            "markers": text(cells[3]),
        }
    )
if len({row["cl_name"] for row in rows}) != len(rows):
    raise ValueError("two level-2 labels share a Cell Ontology name")
if {row["label"] for row in rows} != set(L1_GROUP):
    raise ValueError("the published level-2 labels differ from the declared L1 groups")

with open(HERE / "azimuth_pbmc_l2.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

# Apply the negative rule ----------------------------------------------------------------

positives = {row["label"]: [g.strip() for g in row["markers"].split(",")] for row in rows}
order = [row["label"] for row in rows]


def overlap(left: str, right: str) -> float:
    a, b = set(positives[left]), set(positives[right])
    return len(a & b) / len(a | b)


vocabulary = {}
for row in rows:
    label = row["label"]
    siblings = [s for s in order if s != label and L1_GROUP[s] == L1_GROUP[label]]
    nearest = sorted(siblings, key=lambda s: (-overlap(label, s), order.index(s)))
    negatives = []
    for sibling in nearest[:N_SIBLINGS]:
        for gene in positives[sibling]:
            if gene not in positives[label] and gene not in negatives:
                negatives.append(gene)
    # Programs are named by their Cell Ontology term, the form AIDA's truth labels use,
    # so exact-label metrics compare like with like. The Azimuth label is kept beside it.
    vocabulary[row["cl_name"]] = {
        "azimuth_label": label,
        "cl_id": row["cl_id"],
        "pos": positives[label],
        "negative": negatives[:MAX_NEGATIVES],
    }

(HERE / "pbmc_markers.json").write_text(json.dumps(vocabulary, indent=2) + "\n")
print(f"{len(vocabulary)} programs written")

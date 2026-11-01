# Provenance of the protocol v5 PBMC vocabulary

**Source.** The Human - PBMC `celltype.l2` table on
<https://azimuth.hubmapconsortium.org/references/>, the reference of Hao et al.,
*Cell* 2021. Saved on 2026-11-02 as `azimuth_references_2026-11-02.html`
(sha256 `3e52b81b169bb613ed47f5396d57622966d501cdddca80b1e54dabac9d1c02b8`).

**What was taken as published, with no edit.** All 30 level-2 labels, their ten
markers each in published order (the positives), and the Cell Ontology term
Azimuth links for each label. `azimuth_pbmc_l2.tsv` (sha256 `2b574e7fcee05e8f12067c2951f219c292a91253b916c8da74f28737811cfa4b`) is that
table, extracted by `build_vocabulary.py`.

**What the evaluator decided, all before any v5 run.**

1. Programs are named by their Cell Ontology term name, so exact-label metrics
   compare with AIDA's labels in the same form. The Azimuth label is kept as
   `azimuth_label`.
2. Negatives come from one rule, applied by code: the positives of the two
   level-2 types in the same Azimuth level-1 group with the highest Jaccard
   overlap of positives, minus the type's own positives, in published order,
   at most 6.
3. The level-1 grouping is written out in `build_vocabulary.py`.

`pbmc_markers.json` (sha256 `a549c09e899170667a773500a7285ddbac58d650352e623bb285d301c0f4d2ca`) is the result. No gene was added, removed,
or reordered by hand.

**Known consequences, recorded before running.** All come from the vocabulary and
the 2,000-gene panel, not from AIDA labels or clusters.

- Eryth, CD4 Proliferating, CD8 Proliferating, and NK Proliferating have at most
  3 of 10 positives in the panel, and ILC has 4. All five are below
  `min_positive_coverage` 0.5 and are expected to abstain.
- Platelet's negatives are erythroid genes, because Eryth is its nearest
  sibling. No NK gene is a platelet negative, so the vocabulary does not target
  the v12 error that called an NK cluster platelet.
- NK's negatives are proliferation genes, because NK Proliferating is its
  nearest sibling.
- Azimuth maps NK Proliferating to the generic `natural killer cell`
  (CL:0000623), which is also an AIDA truth term. A cluster called NK
  Proliferating therefore matches AIDA's generic NK label exactly.

**Crosswalk to the 28 AIDA truth types** (cells, from `ground_truth.tsv` and
the pinned `cl.obo`; no annotation was run to compute it):

| Best available Azimuth term for the truth term | Cells | Fraction |
| --- | --- | --- |
| exact | 8,878 | 0.510 |
| more specific only | 5,483 | 0.315 |
| coarser only | 2,001 | 0.115 |
| none (CD8 cytotoxic, double-negative Treg, plasma cell) | 1,032 | 0.059 |

Descendant accuracy credits exact and more specific calls, so the vocabulary
caps it at 0.825.

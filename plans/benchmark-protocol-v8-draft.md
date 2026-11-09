# Protocol v8 (draft): AIFI corpus, tested on Tabula Sapiens immune cells

**Status: DECLARED 2026-11-09** in `benchmarks/README.md` ("Protocol v8"); this file is the working record. Written 2026-11-09. Nothing here is declared until the
user reviews it and it moves into `benchmarks/README.md` under a dated "Declared"
line, before any v8 run. No data is downloaded and no run has happened.
Separate from v7 (user, 2026-11-09): v7 stays as drafted.

## Claim and showcase goals

The package learns a reviewed, ontology-checked corpus from one atlas and annotates
another with evidence-gated labels. Showcase goals (user, 2026-11-09): accuracy
against comparators, ontology and evidence quality, and a runnable tutorial.
"Effort saved" was dropped (Q18). The v6 stop rules apply (Q14): v6 failed
*Not comparable* and *Dominated*, and v8 does not claim a result before it runs.

## Data

- **Corpus source.** Allen Institute Sound Life atlas, CELLxGENE collection
  `e9360edf-b0b7-4e01-bce8-e596814f13e7`: 8 datasets, 13,789,548 cells, 46 GB H5AD,
  blood, 38 Cell Ontology terms (Q1, Q2). Paper: doi 10.1038/s41586-025-09686-5.
- **Query.** Tabula Sapiens - Immune, CELLxGENE dataset
  `c5d88abe-f23a-45fa-a534-788985e93dad` (592,317 cells, 24 donors). Keep blood,
  spleen and lymph node, 10x assays only (Q3). Truth is `cell_ontology_class`,
  scored as in v6 (Q4).
- Both counts come from the CELLxGENE curation API on 2026-11-09, not from the files.
- Work directory `/exports/para-lipg-hpc/mdmanurung/cc-v8-work`, outside the repo.

## Decisions (user, 2026-11-09)

| Q | Decision | Source |
| --- | --- | --- |
| Q1 | Corpus labels are the 38 CL terms as shipped | user |
| Q2 | Use all 8 CELLxGENE files | user |
| Q3 | TS blood + spleen + lymph node, 10x only | user |
| Q4 | Truth `cell_ontology_class`, v6 ontology scoring | user |
| Q5, Q10 | Package `check_test_isolation` stays a hard failure; my study-level overlap check warns | user |
| Q6 | Comparators (CellTypist, SingleR, scANVI, popV) retrained on the AIFI corpus | user |
| Q7 | Corpus markers by cell-level Wilcoxon on a capped sample | user (not my pick) |
| Q8, Q15 | One Leiden run on the TS subset, resolution 1.0 | user |
| Q9 | Truth labels outside the 38 terms form a separate stratum; abstention scores correct | user |
| Q11 | Cap 100 cells per donor per label, seed 1 | user |
| Q12 | Top 10 positives, 5 negatives per label, FDR 0.05 | user |
| Q13 | No batch correction before clustering | user |
| Q14 | v6 endpoint (descendant accuracy, risk-coverage) and four stop rules | user |
| Q16 | Preprocessing: log1p of counts per 10,000 over all genes, then PCA, as v6 | `inherited: benchmarks/README.md` v6, user confirmed |
| Q17 | 2000 HVG, 50 PCs, 15 neighbours, Leiden seed 0 | user |
| Q19 | `knowledge.network: offline`, LLM evidence off | user |
| Q20 | tusoskill iterations score on AIFI held-out donors; TS scored once at the end | user |
| Q21 | Dev split: 20% of AIFI donors, whole donors, stratified by age group, sex, CMV | user |
| Q22 | 50 tusoskill iterations, no early stop | user |

Q11, Q12, Q15 and Q17 numbers are scanpy or Seurat conventions with no source
checked. They are the user's choices and are frozen as given.

## Open before freezing

1. **Marker cap with few donors.** Q7 and Q11 together: a label with fewer than about
   5 donors is learned from very few donors. Report donors per label before any run.
2. **Comparator abstention (Q9).** SingleR, scANVI, popV and CellTypist may not
   abstain. If not, the outside-vocabulary rows are not comparable and are reported
   for cell-curator alone.
3. **Pretrained-model leakage.** Q6 retrains comparators on AIFI. Whether the shipped
   CellTypist and popV models include Tabula Sapiens is unverified (model cards not
   opened). It matters only if pretrained models are also reported.
4. **Corpus review.** 38 records, one OK each in chat, per `AGENTS.md` gate 3.
   Reviewer name goes in `--reviewer`. Not started.
5. **AIFI and TS donor overlap.** Study-level comparison only warns (Q5). Check
   before any run.
6. **tusoskill scope.** Iterations change corpus building and engine settings on the
   AIFI dev split only. Any change found there is a new row the user approves. TS is
   read once, after the freeze commit.
7. **Predicted effect.** Not stated. Add before freezing, as reported-either-way.
8. **Wall-clock cap and download plan.** About 66 GB of H5AD. Slurm node, memory and
   cap not set.
9. **Gate A** assumptions approval stays the user's step.

## Not done

No plot, table of results or accuracy number exists for v8. Any later conclusion
needs an inspected saved figure (`plot-first` rule) or reads
`unsupported: no inspected figure`.

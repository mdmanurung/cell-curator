# v8 consensus mode: design plan (not implemented)

Status: draft for the user's approval. No `src/` change has been made.
Decisions by the user, 2026-11-12: provenance = per cluster plus per-cell flags; scope now =
this design plan only. The role of the marker evidence was first set to referee, then
replaced in round 2 (below) by equal co-voter plus ontology tie-break.

## Goal

cell-curator builds a consensus label from the annotations of the standalone tools used
on a dataset (CellTypist, popV, SingleR, scANVI, Azimuth), breaks ties between them, and
records why each cluster, and each cell that departs from its cluster, got its label.
Accuracy is not the pitch. The pitch is an ontology-checked label with a reason.

## What exists today

- `levels/L1/labels.tsv` has `rationale`, `rejected_alternatives` (with scores),
  `confidence`, `margin`, `agreement`, `reviewer` per cluster. Example, AIDA run 003
  cluster 0: "vote share below threshold", rejected: natural killer cell (score 1.339),
  mature gamma-delta T cell (1.014).
- `decision.manifest.json` and `evidence.manifest.json` record inputs and artifacts.
- `annotation.py:604-650`: classifier predictions are mapped to the corpus vocabulary
  per cluster, an agreement share is computed, and a share below
  `min_classifier_agreement` (0.5) is an advisory flag (`classifier_disagreement`), or a
  gate failure if the classifier is `required`. A classifier never supplies a label.
- The v8 dev runs set `classifier_resources: {}`, so no classifier was used (AIDA run
  checked; AIFI config not opened).

## Gaps

1. No record of what each tool called for each cluster, or which one decided.
2. No per-cell record. Cells inherit the cluster label (`cluster_propagated_to_cells`).
3. No mode where tools propose labels and the marker evidence chooses.

## Design

Role of the marker evidence: referee.

1. Candidate set per cluster: each tool's majority label over the cluster's cells, mapped
   to a CL term and collapsed along ontology lineage (a call and its ancestor or
   descendant count as one candidate, kept at the most specific term the corpus holds).
2. Referee rule: among the candidates, pick the one with the highest marker vote share in
   the cluster. Call it only if that share passes the vote gate and the existing evidence
   gates pass (parent remainder, closest sibling, negative markers). Otherwise `Unknown`.
3. Decision rule recorded per cluster: `all_agree`, `referee_pick`, `tie_break`
   (top two candidates within the lead gate; ontology specificity decides), or `abstain`.
4. Cell level: a cell keeps its cluster label. A cell is flagged only if it takes a
   per-cell fallback call or its own tool labels disagree with the cluster label.

Open, to be settled by the benchmark below and by the user:
- vote gate value, and whether votes are smoothed over the kNN graph (smoothing is under
  test in the 0.95-coverage search; the current best smoothed setting is share 0.5,
  fallback margin 0.40);
- how a missing or failed tool is treated;
- whether tools get unequal weight.

## Provenance schema

`levels/L1/consensus.tsv`, one row per cluster:

| column | content |
| --- | --- |
| cluster_id, label, cl_id | the consensus call, or Unknown |
| decision_rule | all_agree / referee_pick / tie_break / abstain |
| tool_calls | JSON: tool -> {label, cl_id, cluster_share} |
| marker_votes | JSON: top terms with vote share |
| gate_failures, gate_flags | as in scores.tsv |
| rationale, rejected_alternatives | text, as in labels.tsv |
| n_cells, n_cells_flagged | counts |

`levels/L1/cell_flags.tsv`: only flagged cells (barcode, cluster_id, reason, tool_labels).
The manifest gains, per tool: prediction file sha256, tool name and version, model name.
Size: AIFI has 339,530 cells; the flag file is a small share of them (not measured yet).

## Benchmark before any code

Offline harmonizer test (`cc-v8-work/tusoskill/cov95-2026-11-11/candidates/c012_harmonizer/`):
each tool alone, plain majority vote, marker evidence alone, and referee, co-voter and
tie-break designs; gate, weight and tie-break rank chosen on one dataset and applied frozen
to the other. Tools: CellTypist, popV, SingleR on both datasets (scANVI AIFI only; it
failed on AIDA with a gene-count mismatch).

Measured (descendant accuracy among calls, coverage in brackets; figures
`figs/riskcov_aida_split_a.png` and `figs/riskcov_aifi_dev.png`, both inspected):

| design | AIDA frozen | AIFI frozen |
| --- | --- | --- |
| best single tool (popV) | 0.835 (0.973) | 0.888 (0.998) |
| plain majority, no gate | 0.818 (1.000) | 0.894 (1.000) |
| majority + vote-share gate | 0.841 (0.945) | 0.894 (1.000) |
| referee per cell (H2) | 0.811 (1.000) | 0.872 (0.933) |
| referee per cluster (H1) | 0.845 (0.844) | 0.763 (0.970) |
| marker evidence as equal co-voter (H4) | 0.849 (0.907) | 0.900 (0.988) |
| marker evidence alone | 0.730 (0.95) | 0.695 (0.95) |

Findings: the co-voter role is best at coverage >= 0.95 and adds about 1 point over plain
majority on AIFI and nothing on AIDA; the referee role is worse and, per cluster, unstable.
No design reaches 0.85 at coverage >= 0.95 on AIDA. All three tools agree on 70% (AIDA) and
83% (AIFI) of cells and are wrong on 7.6% and 4.4% of all cells, so agreement is not proof;
marker top-1 disagrees with the agreed wrong term on 21% (AIDA) and 39% (AIFI) of those
cells. Caveats: AIDA truth is coarse; one weight grid and one ontology-snap rule tried;
the comparators were trained on AIFI, so AIFI is in-domain for them.

Consequence for this plan: the measured evidence favours the co-voter role over the referee
role chosen on 2026-11-12. The user decides whether to change it.

Option 3 (co-voter plus ontology tie-break: more specific term, then marker share, then best
tool; `candidates/c014_tiebreak/`, figures opened) ties the co-voter: per cell it is identical
(AIDA 0.8441 at 0.907, AIFI 0.8920 at 0.9885); per cluster it is +0.0024 [-0.0085, 0.0140]
(AIDA) and 0.0000 [-0.0018, 0.0015] (AIFI), so no gain. With three tools and one marker vote a
winning share of 0.5 needs two tool votes, so the marker vote cannot flip a majority. Cells that
reach a tie-break are the hard ones: accuracy before the gate is 0.885 (AIDA) and 0.941 (AIFI)
for all_agree, 0.70 and 0.65 for majority, 0.49 and 0.32 for the marker tie-break. Specificity
is right more often than wrong but covers only 31 to 727 cells; marker flips against a plain
majority are net wrong on AIFI (1,235 right vs 3,244 wrong, per cell, delta 0.5).

## Decisions, round 1 (user, 2026-11-12)

| # | Question | Answer |
| --- | --- | --- |
| 1 | Pass bar | Non-inferior to the best single tool (within 1 point at coverage >= 0.95) and every call has a recorded reason |
| 2 | Headline dataset | AIDA and AIFI equal. Report both, with AIFI marked in-domain for the AIFI-trained tools |
| 3 | Unresolved clusters | Two columns: strict label (`Unknown` kept) and a best-guess label flagged low-confidence |
| 4 | Tool set | Matched (AIFI-trained) CellTypist, popV, SingleR as primary; public pretrained models as a sensitivity run |

## Decisions, round 2 (user, 2026-11-12)

| # | Question | Answer |
| --- | --- | --- |
| 5 | Role of the marker evidence | Equal co-voter plus ontology tie-break. Supersedes the referee choice above. Measured: identical to the co-voter per cell; the tie-break step "marker share" is net wrong on AIFI (1,235 right vs 3,244 wrong), the step "more specific term" is right more often (231 vs 56). Whether to keep the marker-share step is open |
| 6 | Gate for the strict column | Frozen from the other dataset (chosen on one, applied unchanged to the other) |
| 7 | Last 3 iterations | Public-model sensitivity run (cached public CellTypist models, plus any public popV and SingleR models found), run through the new `annotation-consensus` skill |
| 8 | Commit | Now, on local main (nothing pushed); the skill directory is committed after its workflow finishes |

## Marker analysis against known markers (`candidates/c015_markers`, figures opened)

Wilcoxon and logistic-regression markers (top 50 by both against the rest of the parent) tested
against PanglaoDB, the Azimuth PBMC file and the package corpus (hypergeometric, BH within run
and source; celltype.info gave no usable entries). Top hit matches truth on called vs Unknown
clusters: AIDA 003 Azimuth 4/8 vs 5/15, corpus 5/8 vs 8/15, PanglaoDB 0/8 vs 5/15; AIFI 002
Azimuth 3/12 vs 1/10, corpus 10/12 vs 2/10 (circular, learned from AIFI), PanglaoDB 3/12 vs
2/10. Almost every top hit has FDR < 0.05, so significance does not separate right from wrong
hits. Known-marker enrichment resolves 5 (Azimuth) to 8 (corpus) of the 15 Unknown AIDA
clusters (B cell, plasma cell, hematopoietic stem cell, CD141+ DC, NK subsets). Caveats found by
the skeptic: BH was applied per run and source, not run-wide as planned (recomputing run-wide
leaves the significant-hit counts unchanged); the AIFI corpus circularity flag is in the heatmap
title and README only; AIDA truth is coarse, so descendant credit is generous.

Still open: whether marker enrichment becomes a voter; whether TF activity, pathway and
cytokine-response scores (CollecTRI, PROGENy, huCIRA) enter as provenance, state calls or
votes (round 3); when the novel sets are read.

## Package changes (not made)

- `config.py`: classifier role `candidate` beside `required` and advisory.
- `annotation.py` (604-650) and `hierarchy.py`: the candidate path and referee rule.
- `reconciliation.py`, `review.py`: new columns, `consensus.tsv`, `cell_flags.tsv`.
- `provenance.py`: per-tool hashes and versions in the manifests.
- Tests for each rule and for the schema.

## Gates that do not change

Gate A, critics, corpus review and write-back are untouched. A consensus label is
evidence like any other: it never writes back by itself. LLM evidence stays limited to
breaking a near-tie the markers left open.

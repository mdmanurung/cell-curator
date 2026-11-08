# Annotation corpus implementation validation

Validated locally on 2026-11-08. Changes remain in the working tree; no release,
push, production annotation, or production write-back was performed.

## Implemented scope

- Strict, versioned corpus snapshots with descriptions, context, source provenance,
  signed markers, ontology relations, source mappings, state/study scope, and
  reviewed proposals. Context cannot override dataset biology, and colliding
  sources require explicit mapping. Each plan freezes the corpus bytes.
- Bounded imports of saved normalized CAP, Allen, and generic records. Imports
  remain unreviewed. The packaged examples are not a callable biological corpus.
- CPU Wilcoxon and balanced logistic-regression markers with actual convergence
  receipts, parent-local comparisons, positive/negative checks, explicit method
  failures, and deterministic abstention or retention of a supported parent.
- Raw CellTypist predictions and a local pretrained popV CPU-fast adapter.
  Classifier evidence cannot override measured marker gates. Optional isolated
  interpreters avoid popV/workflow dependency conflicts.
- A Snakemake DAG for classifier evidence, each hierarchy level, final assembly,
  and strict review. L1–L3 identity labels stay separate from the L4 study/state
  display. Existing assumptions, critic reconciliation, and human write-back
  gates remain enforced.

## Execution evidence

| Check | Observed result | What the result establishes |
| --- | --- | --- |
| Complete regression suite | 668 passed, 3 skipped, 55 warnings; 2055.06 seconds | Existing and new behavior passed the collected suite, including the real PBMC3k tutorial |
| Actual synthetic Snakemake annotation | All three identity levels and L4 matched the constructed truth | The compiled computational DAG executes the shared implementation |
| Snakemake replay | No work scheduled; final manifest unchanged | Completed artifacts resume without recomputation |
| Actual Snakemake review target | Strict report, cards, colors, and critic packet created | Review compilation works without approving findings or applying labels |
| Final focused checks | State approval, worker isolation, progress, review preview, and frozen packet checks passed | Later integration changes retain their declared gates and invariants |
| Ruff, mypy, schema, lock, and diff checks | Passed | Static checks and shipped configuration/dependency metadata agree |

The complete-suite output is retained in
`plans/annotation-corpus-validation/pytest-full.log`. Later focused checks follow
integration edits made while that longer run was executing; the full-suite
result is not a claim that every later edit was exercised by that same process.

The three skipped tests require actual GPU hardware and CuPy. Native GPU
execution was not qualified here. CellTypist/popV adapters were tested through
controlled fixtures and mocks; no pretrained classifier accuracy study ran.

The public PBMC3k tutorial exercises the existing workflow with the updated CPU
marker implementation. It is not a comparative benchmark or an independent
biological validation of the new corpus.

## Remaining scientific qualification

Review and curate the corpus for the intended species, tissue, assay, stage,
and condition. Provision suitable local models and references before a study
run. Evaluate model performance and new corpus revisions on independent data.
No source is designated a gold standard, and repeated predictions do not
automatically become accepted knowledge.

The corpus workflow currently scores one declared CPU assay. Native GPU and
multimodal lanes remain in the existing workflow. Source thresholds are stored
as assertions; a reviewer must choose executable study rules. Corpus imports do
not scrape live websites or synchronize CAP's API.

See the [user guide](../skills/cell-curator/references/annotation-corpus.md) and
[compiled DAG](../workflow/annotation.smk) for commands and output contracts.

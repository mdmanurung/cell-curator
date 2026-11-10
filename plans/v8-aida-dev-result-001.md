# v8 AIDA split A dev result 001 (provisional, pre-validate)

Run `v8-dev-aida-split-a-003`, pinned worktree at a12664d, job 25991468 (16 min).
Predictions come from `levels/L1/labels.tsv`, not `run export-predictions` (blocked on
critic findings). Not for write-back.

| Method | Coverage | Descendant accuracy among calls | Strict accuracy (all cells) |
| --- | --- | --- | --- |
| cell-curator, 34 records | 0.386 | 0.892 | 0.175 |
| CellTypist, AIFI-trained | 1.000 | 0.763 | 0.360 |
| popV, AIFI-trained | 0.973 | 0.835 | 0.397 |

At the matched coverage 0.386: descendant 0.892 / 0.894 / 0.917 and strict among calls
0.452 / 0.441 / 0.446 (cell-curator / CellTypist / popV). Table:
`cc-v8-work/mvp-aida-2026-11-10/compare/matched_coverage.tsv`. Figure inspected:
`cc-v8-work/mvp-aida-2026-11-10/compare/risk_coverage.png`.

Why 15 of 23 clusters are `Unknown`: no term reaches 60% of cell votes with a 20-point
lead. Votes split between related terms (cluster 0: CD56-dim NK 0.53, NK 0.31; cluster 16:
naive B 0.52, transitional B 0.48). Source: `levels/L1/scores.tsv` (`table fallback:
per-cluster vote fractions`).

Negative result: back-off and looser gates raise coverage only by lowering accuracy among
calls (`fig_backoff_floor_sim.png`). Not adopted.

Not covered: SingleR, scANVI (mapping failed: gene count mismatch), AIFI dev, novel sets.
One dataset, dev split, one seed.

# Protocol v7 (draft): the corpus workflow on AIDA split B

**Status: DRAFT. Not frozen.** Written 2026-11-08. Nothing here is declared until you
review it and it moves into `benchmarks/README.md` under a dated "Declared" line,
before any v7 run. Items marked `[user decides]` are open. Items marked `[proposed]`
are my suggestion and carry no authority until you accept them.

Source of the design: tracker rows 21 and 22 in `plans/annotation-corpus-decisions.md`.

## What v7 claims, and what it does not

The claim stays the v6 claim: **at matched coverage, cell-curator's calls are
comparable to the named comparators, and its abstention and confidence are usefully
calibrated.** v6 withdrew that claim (both stop rules fired on v14). v7 tests whether
the corpus workflow changes the outcome. It does not claim a published-leaderboard
result. No accuracy claim exists for the corpus workflow today.

## What changes from v6

Three things change at once, so v7 is not a single-change protocol.

1. **Engine.** Tracker rows 2 to 14, 16 to 18 and 20: unmeasured negatives flagged
   not failed, family BH, at-least-2-marker gate, advisory classifiers, and the rest.
2. **Knowledge.** The v5/v6 vocabulary (derived from Azimuth) is replaced by a corpus
   imported from scGate_models v0.16, human `PBMC`, 13 labels.
3. **Isolation.** `annotation.isolation` is set, so the contract refuses a corpus
   that names a split B donor.

**Attribution problem.** With all three changed, a v6-to-v7 difference cannot be put
down to one of them. `[user decides]` Q1.

## Decisions (round 1, approved by the user 2026-11-08)

| Q | Decision |
| --- | --- |
| Q1 | Two runs. **v7a**: new engine plus the v6 vocabulary converted into a corpus. **v7b**: new engine plus the scGate PBMC corpus. A converter from `pbmc_markers.json` to corpus records is needed first. |
| Q2 | Time the CPU methods on a Slurm node before declaring v7. Job 25978636, `benchmarks/scripts/time_marker_methods.py`, output `cc-v5-work/timing-row19-2026-11-08/timing.tsv`. Any gene prefilter becomes its own tracker row. |
| Q3 | The user reviews each of the 13 scGate PBMC records from a review packet, one OK each in chat. The reviewed corpus sha256 is then frozen here. |
| Q4 | Accept the flat L1 hierarchy for v7b. Deeper truth types cannot be called; they score only through descendant and relaxed agreement. |
| Q5 | Reuse the v6 comparator results unchanged, citing their result file hashes. |
| Q6 | Report a coverage-matched comparison beside the primary curves. |
| Q7 | A wall-clock cap per run, set after the Q2 timing. Number not chosen yet. |

## Still open

- The wall-clock cap (Q7), after the timing result.
- How the v6 vocabulary maps to a corpus for v7a: levels, parents, tissue. The v6 file has its own hierarchy; the converter must not invent one. `[user decides]` once the converter's output is shown.
- The predicted-effect text, once v7a and v7b are both defined.

## Fixed from v6 (unchanged)

- **Dataset.** AIDA Phase 1 Data Freeze v2, dataset `c838aec3-03ef-4398-b882-0e3912abfff0`,
  the whole-transcriptome slice (17,394 cells, 61,497 genes, sha256 `1e74a9d6…`),
  same cell order, same `bench_cluster`, `annotation_scope`, `X_pca_bench`.
- **Split.** `data/aida-pbmc-protocol-v5/donor_split.tsv`, seed 0. **Scored on split B**
  (8,547 cells); full-slice numbers reported beside it.
- **Ontology and scoring.** Pinned `cl.obo`; exact, descendant, relaxed, coverage,
  selective, macro F1; strict and descendant risk-coverage sweeps.
- **Thresholds.** Every `evidence` threshold as in `data/aida-pbmc-protocol-v6/config.yaml`.
  New v7 parameters take the approved defaults: `min_state_positive_fraction` 0.5,
  `min_state_fraction` 0.5, marker gate of 2 positives, top 100 by effect after the
  adjusted-p filter.
- **Gate A.** Interactive, approved by you.

## Isolation

- Config sets `annotation.isolation`: `donor_split` = the v5 split file, `test_split` =
  `B`, `dataset_id` = the AIDA dataset id above.
- Run `cell-curator corpus check-isolation` first and keep its output with the run.
- **What this proves.** Only what records declare. scGate records carry no
  `derived_from`, so the check passes trivially.
- **What it does not prove.** scGate_models were built from public data. I did not
  check whether any of it overlaps AIDA donors. `[unknown]` State this as a limit in
  the result.

## Pre-run checklist (all before the first v7 run)

1. Q1 to Q6 answered and this file moved into `benchmarks/README.md` as "Declared".
2. Row 19 timing done; any prefilter approved as its own row.
3. scGate_models clone at tag v0.16 (`git describe --tags --exact-match`), path and
   file hashes recorded. Nothing from it is committed to this repository.
4. Corpus imported with `--source-kind scgate`, reviewed by you, hash frozen.
5. `corpus check-isolation` output saved.
6. Clean worktree; full test suite passing at the run commit; commit SHA recorded.
7. Predictions exported with `export benchmark-predictions` and scored with
   `benchmark evaluate` into a new `results/aida-pbmc-v15*` set. Run id `aida-pbmc-v15`
   (v14 was v6's). `publish_json` refuses a rewrite, so the first result is final.

## Predicted effect, stated before measuring `[proposed]`

- Coverage rises above v14's 0.298, because corpus labels carry both positives and
  negatives and the family-BH and 2-marker rules differ from v14's gate. **Not
  predicted:** whether calls are right.
- Calls are limited to the 13 corpus labels. Truth types below them (for example naive
  CD4 T cell) can score only through descendant and relaxed agreement.
- Whether the v14 error "natural killer cell called platelet at 0.976" repeats is
  **not predicted.** Row 30 exclusions make it testable: the corpus NK and Platelet
  labels exclude each other's signatures. Report cluster 12 explicitly either way.

## Stop rules (unchanged, quoted from `benchmarks/README.md`)

- **Harness invalid.** If the synthetic re-derivation in `tests/test_benchmark_instrument.py`
  fails, no number is read.
- **Not comparable.** If `descendant_accuracy` against AIDA is below 0.60, the problem is
  label taxonomy. Reported as the finding; the corpus is not retuned to chase the number.
- **Dominated.** If the risk-coverage curve is dominated by a comparator across its whole
  range, the selective-prediction claim is withdrawn. Not re-tuned.
- **Blocked.** If a bounded slice cannot be read within 4 GB, or the run cannot finish
  within the budget set after Q2, the run stops. No silent downsampling.

## What gets reported

- Both slices (full and split B) in one table with v14 and the v6 comparators beside them.
- Flag counts from the new rules: `unmeasured_negatives`, `no_declared_negatives`,
  `comparator_unavailable:*`, `exclusion_unmeasured:*`, `classifier_disagreement:*`.
- Per-cluster calls for clusters 12 and 26 (the two confident errors in v14).
- A figure for every conclusion, opened and inspected before it is cited. Until then a
  conclusion is written `unsupported: no inspected figure`.
- Threats: row 4 does not neutralise labels with no declared negatives; the isolation
  check is declarative; scGate overlap with AIDA is unchecked.

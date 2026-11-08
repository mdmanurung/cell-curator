# Annotation corpus: decision tracker

Opened 2026-11-08. Source: Opus review of `plans/annotation-corpus.md` and
`plans/annotation-corpus-validation.md` against commits `90eab77`..`61b3d89`.

Rules for this file:
- Every row is an impactful choice. It is the user's, never a default.
- Status moves `open` -> `approved: <choice>` -> `done: <commit>`.
- Only the user sets `approved`. An agent never edits that column.
- Code locations come from the review and are not re-verified (`unverified`).
- No row has a supporting figure yet: `unsupported: no inspected figure`.

## Gate-related

| # | Decision | Where (unverified) | Current | Alternatives | Source | Status |
|---|---|---|---|---|---|---|
| 1 | Who approves corpus records | `corpus.py:503-562`, `484-487` | Free-text `--reviewer`; accept makes records callable | Typed terminal token as in write-back; list as an `AGENTS.md` gate | inherited: plans/annotation-corpus.md | open |

## Label-call rules

| # | Decision | Where (unverified) | Current | Alternatives | Source | Status |
|---|---|---|---|---|---|---|
| 2 | Unmeasured negative marker | `annotation.py:384` | NaN fails the gate | Score separately | inherited: plans/annotation-corpus.md line 143 | open |
| 3 | Missing comparator table | `annotation.py:395-396` | Fails `{role}:{method}` | Mark unavailable, skip vote | inherited | open |
| 4 | Label with no declared negatives | `evidence.py:690`, `718` | Coverage 1.0, violation 0.0 | Neutral, or require negatives | inherited | open |
| 5 | Parent with no differing sibling | `annotation.py:357-361` | No sibling rows, all labels fail | Use a corpus sibling | inherited | open |
| 6 | Positive-marker gate | `annotation.py:395` | 1 marker in top 100 | Minimum count or fraction | inherited | open |
| 7 | Top-N ranking | `evidence.py:224` | Mean difference, N=100 | Rank by p-value; other N | inherited | open |
| 8 | Multiple-testing family | `evidence.py` | BH per cluster per comparison | One family across siblings | inherited | open |
| 9 | Exact-mapping merge | `corpus.py:373` | Union of positives | Intersection | inherited | open |
| 14 | Margin definition | `annotation.py:417` | Failed runner-up gets vote 0 | Keep its real vote | inherited | open |
| 15 | LLM near-tie hook | corpus DAG | `evidence.llm.proposals` not read | Wire in, or declare out of scope | inherited: plan line 154 | open |
| 16 | Corpus context scope | `corpus.py:125`, `293` | Empty tissue/stage/condition matches anywhere | Require tissue and assay | inherited: plan line 71 | open |
| 17 | State-call rule | `annotation.py:435` | Every positive gene in every cell | Fraction or score | inherited | open |

## Classifier rules

| # | Decision | Where (unverified) | Current | Alternatives | Source | Status |
|---|---|---|---|---|---|---|
| 10 | Classifier role | `config.py:367` | `required` by default | Advisory | inherited | open |
| 11 | popV and CellTypist dependence | `crosschecks.py:966`, `reconciliation.py` | Both can veto; `reconcile_predictions` never called | One vote; call reconciliation | inherited | open |
| 12 | Coarser classifier label at finer level | `annotation.py:262-267` | Counts as disagreement | Counts as compatible | inherited | open |
| 13 | Unmapped raw labels | `annotation.py:406` | Unknown, lowers agreement | Drop from denominator | inherited | open |

## Engine and corpus lifecycle

| # | Decision | Where (unverified) | Current | Alternatives | Source | Status |
|---|---|---|---|---|---|---|
| 18 | CPU marker methods | `decisions.py:311-332` | Wilcoxon and logistic regression both needed | Wilcoxon only (earlier behaviour) | unplanned: Codex edit, commit `85b8077` | open |
| 19 | Logistic regression setup | `evidence.py:304-340` | Balanced L2, L-BFGS, all genes | Gene prefilter, other penalty | unplanned: Codex edit, commit `85b8077` | open |
| 20 | Corpus update rule | `corpus.py:479` | Same `label_id` replaces record | Mark old `superseded` | inherited | open |
| 21 | Learning loop | `experience` unused | No held-out set; run output could feed its own test | Hold out a split | inherited | open |
| 22 | Validation target | `tests/test_annotation_pipeline.py:189` | Circular synthetic truth; no real-data accuracy test | Held-out real data plus a comparator | inherited | open |

## Not decisions, still to do

- Re-run the full suite on the tip (the last full run predates the final edits).
- Check whether the all-gene CPU logistic fits explain the slow PBMC3k run.
- Native CellTypist and popV installs; real-data validation. Not started.

## Suggested order

1. Row 1 (the only gate).
2. Rows 2-4, 6, 7, 14 (change calls most directly).
3. Rows 5, 8-13, 15-17.
4. Rows 18-22.

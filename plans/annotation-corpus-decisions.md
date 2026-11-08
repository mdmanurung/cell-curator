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
| 1 | Who approves corpus records | `corpus.py:503-562`, `484-487` | Free-text `--reviewer`; accept makes records callable | Typed terminal token as in write-back; list as an `AGENTS.md` gate | inherited: plans/annotation-corpus.md | approved: user's own OK in chat per proposal, no token; rule added to `AGENTS.md` (2026-11-08) |

## Label-call rules

| # | Decision | Where (unverified) | Current | Alternatives | Source | Status |
|---|---|---|---|---|---|---|
| 2 | Unmeasured negative marker | `annotation.py:384` | NaN fails the gate | Score separately | inherited: plans/annotation-corpus.md line 143 | approved: pass with `unmeasured_negatives` flag (user, 2026-11-08); done in working tree, uncommitted |
| 3 | Missing comparator table | `annotation.py:395-396` | Fails `{role}:{method}` | Mark unavailable, skip vote | inherited | approved: mark `comparator_unavailable`, skip that vote, show the flag (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |
| 4 | Label with no declared negatives | `evidence.py:690`, `718` | Coverage 1.0, violation 0.0 | Neutral, or require negatives | inherited | approved: pass with `no_declared_negatives` flag; no score bonus (neutral coverage) (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |
| 5 | Parent with no differing sibling | `annotation.py:357-361` | No sibling rows, all labels fail | Use a corpus sibling | inherited | approved: use a corpus sibling of the same parent; nearest cluster as fallback (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |
| 6 | Positive-marker gate | `annotation.py:395` | 1 marker in top 100 | Minimum count or fraction | inherited | approved: at least 2 positive markers, or all if the label has only 1 (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |
| 7 | Top-N ranking | `evidence.py:224` | Mean difference, N=100 | Rank by p-value; other N | inherited | approved: filter by adjusted p first, rank by effect, N=100 (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |
| 8 | Multiple-testing family | `evidence.py` | BH per cluster per comparison | One family across siblings | inherited | approved: one BH family per parent scope and method, across all sibling tests (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |
| 9 | Exact-mapping merge | `corpus.py:373` | Union of positives | Intersection | inherited | approved: intersection of positives; each source's own list kept in the record (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |
| 14 | Margin definition | `annotation.py:417` | Failed runner-up gets vote 0 | Keep its real vote | inherited | approved: keep the failed runner-up's real vote so near-ties stay visible (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |
| 15 | LLM near-tie hook | corpus DAG | `evidence.llm.proposals` not read | Wire in, or declare out of scope | inherited: plan line 154 | approved: out of scope for now; a near-tie goes to the parent or Unknown (user, 2026-11-08); no code change needed |
| 16 | Corpus context scope | `corpus.py:125`, `293` | Empty tissue/stage/condition matches anywhere | Require tissue and assay | inherited: plan line 71 | approved: require species and tissue; assay and condition optional but recorded (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |
| 17 | State-call rule | `annotation.py:435` | Every positive gene in every cell | Fraction or score | inherited | approved: at least half the positive genes per cell and no negative gene; threshold in config with its reason (user, 2026-11-08); implemented in the working tree 2026-11-08 (see Status) |

## Classifier rules

| # | Decision | Where (unverified) | Current | Alternatives | Source | Status |
|---|---|---|---|---|---|---|
| 10 | Classifier role | `config.py:367` | `required` by default | Advisory | inherited | approved: advisory by default; disagreement shown as a flag and in the review packet; `required` settable per classifier in config (user, 2026-11-08); implemented in the working tree 2026-11-08 |
| 11 | popV and CellTypist dependence | `crosschecks.py:966`, `reconciliation.py` | Both can veto; `reconcile_predictions` never called | One vote; call reconciliation | inherited | approved: popV and its embedded CellTypist count as one vote, via `reconcile_predictions` (user, 2026-11-08); implemented in the working tree 2026-11-08 |
| 12 | Coarser classifier label at finer level | `annotation.py:262-267` | Counts as disagreement | Counts as compatible | inherited | approved: abstain: leave out of the agreement rate and report the count (user, 2026-11-08); implemented in the working tree 2026-11-08 |
| 13 | Unmapped raw labels | `annotation.py:406` | Unknown, lowers agreement | Drop from denominator | inherited | approved: drop from the agreement rate and report the unmapped count (user, 2026-11-08); implemented in the working tree 2026-11-08 |

## Engine and corpus lifecycle

| # | Decision | Where (unverified) | Current | Alternatives | Source | Status |
|---|---|---|---|---|---|---|
| 18 | CPU marker methods | `decisions.py:311-332` | Wilcoxon and logistic regression both needed | Wilcoxon only (earlier behaviour) | unplanned: Codex edit, commit `85b8077` | approved: Wilcoxon is the default; logistic regression only in the corpus workflow (user, 2026-11-08); implemented in the working tree 2026-11-08 |
| 19 | Logistic regression setup | `evidence.py:304-340` | Balanced L2, L-BFGS, all genes | Gene prefilter, other penalty | unplanned: Codex edit, commit `85b8077` | approved: keep balanced L2 on all genes; time one profile first; any prefilter gets its own row (user, 2026-11-08); not yet implemented |
| 20 | Corpus update rule | `corpus.py:479` | Same `label_id` replaces record | Mark old `superseded` | inherited | approved: mark the old record `superseded` with a link; old runs keep their pinned version (user, 2026-11-08); implemented 2026-11-08 |
| 21 | Learning loop | `experience` unused | No held-out set; run output could feed its own test | Hold out a split | inherited | approved: learn only from split A donors (`benchmarks/data/aida-pbmc-protocol-v5/donor_split.tsv`); split B never proposes or accepts a record; the corpus pins the donor list it learned from. No record may cite AIDA or its labels as a source for the benchmark run: records get `derived_from`, and the evaluator refuses any link to test donors (user, 2026-11-08); not yet implemented |
| 22 | Validation target | `tests/test_annotation_pipeline.py:189` | Circular synthetic truth; no real-data accuracy test | Held-out real data plus a comparator | inherited | approved: new frozen protocol v7 declared before any run: AIDA split B through the corpus workflow; comparators CellTypist, SingleR, scANVI, popV plus v14 cell-curator as baseline (no scGate); corpus = scGate_models `generic` and `PBMC` at tag v0.16, every record reviewed by the user before the run; headline claim stays selective prediction at matched coverage and calibrated abstention, plus unknown rate and the row 2, 3, 4 flag counts reported; v6 stop rules reused unchanged; synthetic test kept as a plumbing check only (user, 2026-11-08); not yet implemented |

## External methods added 2026-11-08 (user request)

Facts read with `gh` on 2026-11-08, not from the repos' docs sites (those timed out):

- scGate v1.7.0 (2025-04-24), R, GPL-3 in `DESCRIPTION`. Scores cells with UCell, smooths over kNN, runs a hierarchical gating model, calls `Pure`/`Impure`.
- scGate_models v0.16 (2023-12-18), no licence file seen. Human and mouse. Each model is a TSV with columns `levels`, `use_as` (`positive` or `negative`), `name`, `signature`. `master_table.tsv` holds `name` and `signature` (genes joined by `;`, a trailing `-` marks a negative gene). Collections: PBMC, TME_broad, TME_HiRes, HiTME, CD4_TIL, CD8_TIL. `CellOntology_dictionary.tsv` maps names to Cell Ontology.
- MMoCHi, Python, GPL-3.0, no release tag. Hierarchy of random forests trained on high-confidence thresholded cells, built for CITE-seq. Uses ADT, RNA and metadata.

| # | Decision | Alternatives | What it changes | Source | Status |
|---|---|---|---|---|---|
| 23 | Role of scGate_models (first import: human `generic` and `PBMC`; TME and mouse after the first review round, approved 2026-11-08) | Import as corpus source (via `corpus import`); reference only; runtime download | Which labels and negatives the corpus holds | user request 2026-11-08 | approved: import as corpus source, records start as proposals needing the user's OK (user, 2026-11-08); licence (row 24) must clear before any file is shipped |
| 24 | Licence and bundling | Ship copies; fetch at run time and pin; do not ship | Whether GPL-3 or unlicensed files enter the package. Needs a licence check first (scGate_models has none) | user request 2026-11-08 | approved: (a) ship nothing; import at run time from a clone the user points to; records carry source URL, tag and file hash. (b) asking carmonalab for a licence is the user's to send (user, 2026-11-08) |
| 25 | Map scGate negatives | Level-wise negative = a whole population signature; `GENE-` = a single gene. Map both to corpus `negative`, or only the single genes | Negative-marker gate (rows 2-4) | user request 2026-11-08 | approved: level negatives become a "must not match" link to another label, measured through that label's markers; single `GENE-` negatives map to corpus `negative` (user, 2026-11-08) |
| 26 | Run scGate as a classifier | Add as a crosscheck beside CellTypist and popV (R dependency, cell-level, own UCell/kNN); or markers only | Adds an evidence source and a possible veto (rows 10-11) | user request 2026-11-08 | approved: no scGate crosscheck for now; markers only; revisit after real-data validation (user, 2026-11-08) |
| 27 | Use of MMoCHi | Design reference for hierarchy and high-confidence thresholds; or a crosscheck trained on the user's own data (needs labelled cells, CITE-seq suits it best) | Hierarchy rules; a possible training-set circularity (row 21) | user request 2026-11-08 | approved: design reference only; no MMoCHi classifier (user, 2026-11-08) |
| 28 | Version pins | scGate_models tag v0.16, scGate v1.7.0, MMoCHi commit SHA (no tag) | Reproducibility of imported records | user request 2026-11-08 | approved: scGate_models tag v0.16 by default; a newer tag or commit needs an explicit config value and a `corpus diff` first. scGate and MMoCHi pins moot while unused at run time (user, 2026-11-08) |
| 29 | Ontology mapping of imported labels | Use `CellOntology_dictionary.tsv` as-is; or re-map all through the package check | Which terms enter the corpus | user request 2026-11-08 | approved: dictionary is a first guess only; package ontology check runs on every label and a failing term blocks the proposal (user, 2026-11-08) |

| 30 | How a level negative (`scgate_excludes`) acts as a gate | Names refer to signatures in the collection's `master_table.tsv`, and most have no model. Options: (a) import every referenced signature as a helper label and require the cluster not to pass it; (b) test the signature genes directly with the row 17 state rule; (c) keep as reviewer-only text | Whether exclusion has any effect on calls; adds labels to the vocabulary under (a) | arises from rows 23 and 25 | approved: (b) test the signature genes with the row 17 state rule (user, 2026-11-08); implemented 2026-11-08 |

## Not decisions, still to do

- Re-run the full suite on the tip (the last full run predates the final edits).
- Check whether the all-gene CPU logistic fits explain the slow PBMC3k run.
- Native CellTypist and popV installs; real-data validation. Not started.

## Status 2026-11-08

Rows 1 to 30 are approved.
Implemented: rows 1 (`AGENTS.md` rule), 2 to 14, 16 to 18, 20, 23, 24, 28, 29. Rows 25 and 30 are done. Row 15 needs no code. Rows 26 and 27 were decided "no", so no code.
Not yet implemented: row 19 (needs a timing run on real data), rows 21 and 22 (learning rule and frozen protocol v7).

Known limits of what was built:
- Row 4: the label gets the `no_declared_negatives` flag and still has violation 0.0. The score has no penalty for it, so it is not made neutral. A neutral value would need a number you have not chosen.
- Row 8: the family is every kept Wilcoxon test in a parent scope. Sibling-pair tests keep only the target cluster's side.
- Row 9: negatives are still a union across mapped sources. Intersection of positives can make an exact mapping fail with `share no positive marker`.
- Row 18: `compute_bottom_up_markers(methods=("wilcoxon",))` is the default and fits no logistic model. The corpus workflow passes both methods. The CPU receipt accepts Wilcoxon only; GPU lanes still need both. A fixture check against commit `69e1deb` (random counts, 3 groups, 60 genes) gives identical shared columns. No real-data benchmark rerun was done.
- Row 20: a replaced record is kept as `<label_id>~<base revision>` with status `superseded` and `superseded_by` set. `corpus diff` lists it under `added`. A `dispute` decision also replaces a reviewed record, as before; the old one is now kept.
- Rows 23, 24, 28, 29: `--source-kind scgate` imports a local clone at the named tag (default v0.16); nothing is shipped; models become unreviewed labels with the file hash kept; the dictionary gives a first-guess `cl_id`. Rows 25 and 30: `GENE-` genes are negatives; each excluded population with a signature in `master_table.tsv` becomes an `exclusions` entry, and `collect-level` fails the label with `excluded:<name>` when that signature passes the row 17 state rule in the cluster. Excluded names with no signature are listed in `scgate_excludes_without_signature` and do not gate.
- Row 11: the one-vote rule is in `annotation.py` (members of a dependency group share one mean agreement). `reconcile_predictions` is still not called; it reports per-cluster composition and does not combine votes.
- Row 16: reviewed labels now need `tissues`, and the dataset context needs `tissue`. Existing reviewed corpora without tissue will no longer load.

## Suggested order

1. Row 1 (the only gate).
2. Rows 2-4, 6, 7, 14 (change calls most directly).
3. Rows 5, 8-13, 15-17.
4. Rows 18-22.

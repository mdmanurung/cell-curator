# Tutorial: PBMC 3k

This tutorial labels the 10x Genomics PBMC 3k dataset with three public annotation
tools and a marker vote, and records why each cluster got its label. cell-curator chooses
no label itself: the `annotation-consensus` skill does, and `annotation import-consensus`
checks its output against the run and writes the labels.

The worked example is in `skills/annotation-consensus/examples/pbmc3k/` in the
repository. Its README lists the scripts in run order, the two consensus runs (marker
vote on and off), and five figures. PBMC 3k has one donor and no truth labels, so the
example shows the workflow, not accuracy. `Unknown` is a result: the `rationale`
column of the consensus table says which gate stopped each call.

The order is fixed. Run `run scene` and `run plan`; run the consensus skill with
`plan/vocabulary.tsv`; set `annotation.consensus_directory`; then run
`annotation collect-level`, `annotation import-consensus` and `annotation finalize`.

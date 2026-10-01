# Benchmark protocol

**This document is frozen before any result is published.** `benchmark evaluate`
publishes through `publish_json`, which refuses a non-identical rewrite, so the
first result written for a dataset is a permanent commitment. Anything decided
here after a number has been seen is a result that was tuned, and the only honest
record of that is a second dataset and a second protocol version.

`results/synthetic.json` is the sole exception already on record. It predates
`ontology_distance_error` becoming a real measurement, and the value it carries
for that field came from the alias the metric replaced. It is superseded, not
corrected, and `test_the_committed_artifact_still_awaiting_supersession` is a
strict xfail that fires the moment that happens.

## What is and is not being claimed

No accuracy or performance claim has been published for this package, and
neither has one for the two public QuadBio projects this is compared against. The
stated claim is a **selective-prediction** claim: *at matched coverage, this
package's calls are comparable to the named comparators, and its abstention and
confidence are usefully calibrated.* That is a narrower claim than being more
accurate, and it is the one the design can support, because the package gates
labels on evidence and preserves `Unknown` rather than answering every cell.

**A fresh comparison of CellTypist, SingleR, scANVI, and Azimuth on any public
PBMC dataset is confirmatory, not novel.** At least twelve published papers
compare multiple automatic annotators on one dataset against one ground truth;
the closest are popV (*Nat Genet* 2024), scTab (*Nat Commun* 2024), the
SCimilarity paper (*Nature* 2025), and Hu et al. 2025 on HLCA. A scanpy pipeline
adding a data point to that is not a contribution, and this protocol does not
pretend otherwise. The two things that could make it worth doing are recorded
under "Where a novel result could still come from".

## Ground truth, and what it is not

**Dataset.** AIDA Phase 1 Data Freeze v2, CELLxGENE collection
`ced320a1-29f3-47c1-a735-513c7084d508`, dataset
`c838aec3-03ef-4398-b882-0e3912abfff0`, read from Census release `2025-11-08`.
Acquired slice: **17,394 cells x 2,000 highly variable genes, 394 donors, 28
CL types, 75.6 MB**, sha256 `e7833ebb…cee09`, pinned in
`manifests/rna.json`. The two acquisition scripts sit beside the data so the
slice is re-derivable rather than asserted. The 2,000-gene restriction is a
memory constraint, not a design choice: Census peak memory scales with tile
scatter, and a contiguous 2,000-gene read of this cohort measures ~2.6 GB peak
against a 4 GB budget, so the slice is taken as contiguous windows straddling
the joinid axis plus a small scattered backfill for rare types. The whole
transcriptome at the same cell count measures ~5.5 GB and does not fit.

**Composition, recorded before any annotation was run.** The slice is T/NK-heavy
because the joinid axis is not ordered by cell type, so window placement
determines the mix. Counts: CD16-dim NK 2,854; naive CD4 1,729; T cell 1,720;
CD14 mono 1,614; central memory CD4 1,500; naive CD8 1,023; CD8 1,002; CD8
cytotoxic 823; NK 789; gamma-delta T 685; myeloid 530; CD4 470; CD14-low
CD16 mono 438; CD8 memory 343; effector memory CD4 295; iNKT 271; Treg 152; CD4
cytotoxic 151; **B cell 147**; CD1c+ mDC 115; CD16-bright NK 109; double
negative Treg 108; ILC 107; pre-cDC 106; **plasma cell 101**; CD141+ mDC 100;
HSPC 100; **platelet 12**.

The B-cell and rare-type tail is thin, and the slice was **not** re-cut to
improve it. Re-acquiring after inspecting the composition, with no result in
hand, is still tuning the dataset toward the annotation, and the protocol is
worth less if that is available. The cost is stated instead: per-cell-type F1 for
platelet rests on 12 cells and must be read with its support count or not at
all, and the B-cell axis is where this evaluation is weakest. The `scenarios`
field in the manifest records which partitions the slice is meant to exercise and
`rare-types-under-200-cells` names the ones it covers only thinly.

**The caveat that governs how this result may be read.** CELLxGENE
`cell_type_ontology_term_id` is curator-standardized but not
expert-curated. The CZI dataset card states that annotations are "provided by
the original authors according to their chosen annotation methods" which "may
involve expert knowledge using canonical markers **or automated classification
tools such as CellTypist**", and that Discover "does not include any alignment,
cluster or annotation pipeline". The Census schema carries **no per-dataset
provenance field** saying whether a given dataset's labels were expert or
tool-derived. There is no way, from the data, to exclude that AIDA's labels were
produced by a tool.

AIDA is better than the average case in one respect that was checked rather than
assumed: the labelling code is public and readable, uses hand-written marker
panels and two rounds of manual sub-clustering, and calls no automated classifier
for the final labels. Its quality control did use Cell BLAST against
Monaco/ImmGen, which is a reference mapping but an FACS-sorted one, outside the
three tools excluded here, and it gated QC rather than identity.

**Consequence.** AIDA is legitimate ontology-consistent ground truth and is
independent of Azimuth, CellTypist, and scANVI as *training* data, which is what
the circularity requirement is about. It is **not** orthogonal ground truth in the
sense of protein-measured or FACS-sorted labels, and no result computed on it may
be described as measuring cell-type accuracy. It measures agreement with
expert-reviewed CL-mapped author labels. Anyone wanting the stronger claim needs
the bead-sorted PBMC arm of Zheng et al. 2017, which Xu et al. used for scANVI
and which nothing in this table benchmarks.

**Circularity, declared and non-overlapping.** cell-curator cross-checks with
CellTypist, Azimuth, scANVI/scArches, and SCimilarity, so truth derived from any
of them would be circular. Rejected outright: 10x `cell_type_annotations.csv`
for the dataset in `manifests/rna.json` (unobtainable *and* Azimuth-derived, per
10x's own documentation of its Cell Ranger Arc annotation output); the nygc
multimodal PBMC set, which *is* the Azimuth reference at 161,764 cells; Tabula
Sapiens, since Pan-human Azimuth was trained on it.

## Retraction: the 8,000-gene input was row-scrambled (recorded 2026-10-01)

`acquire_slice_8000.py` assigned `source.obs_names` positionally to the rows Census
returns, which come back in `soma_joinid` order; `source.obs` is not in that order,
so every cell's expression was paired with another cell's labels and cluster. Runs
v9 and v10 used that file (input sha256 `a8c5e332…5c8380`), and the committed
`aida-pbmc-v10` result, its predictions, and its receipt are retracted by
`results/aida-pbmc-v10.retraction.json`. They stay in the repository, unedited,
because a published result is a commitment and a retraction is a record, not an
erasure.

What this withdraws in the text below: every 8,000-gene number, including "the
8,000-gene slice clusters *worse* ... ceiling 0.44 against 0.64" and the
conclusion that "clustering on the 2,000 most variable of the measured genes is
therefore not a compromise, it is the better input". Both are what a scrambled
matrix produces and say nothing about panel width. The v2 row of the v1/v2 table
came from the v10 run and is withdrawn with it. Runs on the 2,000-gene file did not
use the 8,000-gene acquisition path and are not affected. The script now reorders
each block by `soma_joinid` and refuses a block missing any source cell.

## Protocol v2: the negative-marker bar, changed after v1 and derived rather than fitted

**v1 is not withdrawn and not hidden.** It labelled 1 cluster of 28 at a
0.2 violation threshold, and that number stands as recorded above. Everything from
here is v2, and v2 changed after v1's outcome was seen. That is the one thing this
document exists to make visible, so it is stated before any v2 number exists.

**What changed.** `max_negative_violation` goes from 0.2 to 0.5, and every marker
program carries three discriminative negative markers rather than two. Nothing
else changed: the cell set, the clustering and its recipe, the measured panel, the
other three thresholds, the marker positives, and the CL identifiers are as v1.

**The threshold is derived, not chosen to raise coverage.** The evidence records
`negative_violation` as the *fraction of a program's measured negative markers
that are flagged*, not the fraction of cells expressing one. A clustering whose
ceiling is 0.636 admits up to 0.364 of a cluster belonging to a neighbouring
lineage, so a neighbour's marker is expected to be flagged in many clusters
through contamination alone. Two markers cannot express that distinction: the
reachable values are 0.0, 0.5, and 1.0, so one marker flagged by contamination
and one marker flagged by a genuine contradiction are the same number, and any
threshold above 0.0 admits every row and one at or below 0.0 rejects every row.
Three markers give 0.0, 0.333, 0.667, and 1.0, so **0.5 separates one flag from
two**: contamination produces one, contradiction produces two. The measured
violation distribution bears this out -- on the 2-marker panels 225 of 280 rows sit
at exactly 0.5, and on the 3-marker panels 101 of 112 sit at 0.333.

**What a v2 number may and may not mean.** With a 0.636 clustering ceiling, no
system can exceed 0.636 per-cell accuracy on this input, and cell-curator's own
negative gate is the thing that detects the ceiling. A v2 result therefore
measures whether naming the dominant type of each cluster is a useful summary of
these clusters, not whether the clusters are pure. If it is, coverage should rise
substantially and accuracy should land near the ceiling, and that agreement is the
result: the tool is doing what it is designed to do and the input is what limits
it. If coverage stays low, the gate is not the cause and the diagnosis above is
wrong.

The honest reading of either outcome is bounded by this section. No v2 number is
comparable to a published figure, and none is a measurement of cell-type accuracy.

### What v2 measured, and what it cost

| | v1 (0.2, 2 negatives) | v2 (0.5, 3 negatives) |
| --- | --- | --- |
| clusters labelled | 1 of 28 | 2 of 28 |
| cells called | 153 | 279 |
| per-cell coverage | 0.0088 | 0.0160 |
| selective accuracy on called cells | 0.961 (147/153) | 0.570 (159/279) |
| cells called and exactly wrong | 6 | 120 |

The derivation was right about its own terms. The negative gate stopped being the
binding constraint exactly as predicted, blocking 12 clusters instead of 24, and
the marker-panel fix held. Coverage roughly doubled.

**The result is still not a benchmark, and the change made precision worse.** The
126 additional cells are one cluster of 1,103 whose dominant true type is CD16-dim
natural killer cell, called `platelet`, of which 12 are platelet. So v2 bought 120
confidently wrong cells to double the denominator, and selective accuracy fell by
39 points. v1's bar was better on both counts: it called fewer cells and they were
nearly all right.

That is a real result about the gate rather than about the tool, and it points the
other way from the change. The derivation assumed contamination flags *one*
negative marker, which held for the 3-marker panels; it did not hold for the
winner, because a wrong winner on a contaminated cluster can also flag only one.
One flag is ambiguous between "this cluster is a third neighbour" and "this cluster
is not what the winner says it is", and 0.5 cannot separate two from one, since
no two-marker program exists any more. **A third state is needed -- panels of four
or five negatives, where one flag is 0.2 or 0.25, two is 0.4 or 0.5, and the gate
can sit between them.** That is the direction the threshold work should have gone,
and it was not found by lowering a bound.

The remaining blockers are the margin and entropy gates, at 26 of 28 and 21 of 28
respectively, and those were left alone deliberately. They are the package's
uncertainty calibration, which is the core claim; a system with those lowered
labels whatever scores highest, which is not a selective-prediction system and
would make every other number here decorative. Median top-2 margin on this input
is 0.099 with 14 competing programs, so the evidence genuinely does not separate
them, and that is a property of the data.

**Two measurements are absent for understood reasons and are not defects.**
`ontology_distance_error`, `descendant_accuracy`, and `relaxed_accuracy` are `None`
because the run configured no ontology, so the pipeline cleared every marker
program's accession and refused to let an unvalidated one reach a proposal. The
exported table therefore carries no predicted identifier and there is nothing to
measure. That is the designed refusal, and it is why a real CL snapshot was
supplied only at scoring time. `unknown_precision` and `unknown_recall` are 0.0
because `Unknown` is spelled `Unknown` on the truth side only for a handful of
cells, so the binary problem is nearly one-sided; the coverage figure is the one to
read, and the protocol says so.

## Protocol v3: per-cell votes, an ontology, and the unscrambled input

**Frozen 2026-10-01, before any v3 run.** v1 and v2 stand as recorded. Nothing in
this section has been run on AIDA, and nothing in it may change once it has: a
change to a threshold, a marker, or the engine after the v3 result is seen is v4.
Run config: `data/aida-pbmc-protocol-v3/config.yaml`, run id `aida-pbmc-v11`.

**Engine.** Commit `6f11904`, 502 tests passing. Four changes since v2 decide the
result:

1. Genes are standardised across cells, not across the 28 cluster means. Each
   cell votes for its best program among siblings, a cell no program fits does
   not vote, and confidence is the share of the cluster voting for the winner.
   v2's confidence was a softmax over cluster scores, which fell as the
   vocabulary grew.
2. The Cell Ontology parser keeps `is_a` edges that carry `{...}` qualifiers
   (28,860 edges; T cell at depth 8, where v2's graph had it at 0).
3. A cluster that fails only on confidence, margin, or entropy is labelled with
   the unique deepest CL ancestor of its near-tied candidates, if that ancestor is
   not a generic term (cell, leukocyte, mononuclear cell, and similar) and the
   candidates under it hold 60% of the votes. The run plan states this as an
   assumption.
4. Calibration is computed over calls only, and `predicted_path` no longer
   carries level and `ROOT` tokens.

**Input.** `aida_pbmc_clustered.h5ad`, sha256 `440ba41a…97c159`: the v1 input,
17,394 cells, 2,000 genes, 28 clusters, ceiling 0.636. The 8,000-gene file is not
used. Alignment was checked before freezing, without looking at any annotation
output. Truth labels match `obs` on 17,394 of 17,394 cells, and marker means follow
the truth: CD14 is 1.99 in CD14 monocytes and at most 0.21 in five other types,
PPBP is 6.23 in platelets.

**Ontology.** `cl.obo`, `data-version: releases/2026-06-08`, sha256 `7e19b4ae…8fc51f`.
v2 ran without one.

**Vocabulary.** `data/aida-pbmc-protocol-v3/pbmc_markers.json`, sha256
`8147b66e…d0e355`. It makes four edits to the v2 file. Each can be decided from the
vocabulary and the ontology alone, and none comes from AIDA labels or output:

1. Four labels take their CL names (for example "CD4-positive T cell" becomes
   "CD4-positive, alpha-beta T cell"), so a label and its accession validate
   together. Markers are unchanged, and repeated entries are deduplicated.
2. The plasmacytoid dendritic cell program pointed at CL:0002010 (pre-conventional
   dendritic cell). It now points at CL:0000784. Its markers (LILRA4, CLEC4C,
   IL3RA) are pDC markers. Consequence, stated now: AIDA's 106 cells labelled
   pre-conventional dendritic cell can no longer match exactly.
3. CD8-positive T cell listed NKG7 as both positive and negative, which the
   package now refuses. NKG7 stays positive, as a cytotoxic effector gene. The
   negative slot takes CD14, the lineage negative the other programs use.
4. "myeloid cell" is dropped. It is the CL ancestor of five other labels, and 6
   of its 11 positives are monocyte markers, so per-cell votes cannot tell it from
   its descendants. Its 530 AIDA cells are reachable only through backoff. Innate
   lymphoid cell (ancestor of NK, distinct markers) and the CD14/CD16 monocyte pair
   (a real CL `is_a`) stay, and the plan records both as assumptions.

**Thresholds.** These are the v2 numbers, re-declared in vote-share units:
- `min_confidence` 0.6: 60% of a cluster's cells vote for the call.
- `min_margin` 0.2: a 20-point lead over the runner-up.
- `max_normalized_entropy` 0.85.
- `min_positive_coverage` 0.5.
- `max_negative_violation` 0.5.

**Lanes and arms.** Local markers only: no Census- or CellGuide-derived lane,
because AIDA is in Census. This freeze covers the arm without an LLM lane. An
LLM arm is an amendment pinned before it runs. It is reported separately and
flagged for possible training leakage.

**How backoff calls are scored.** The comparison with v1 is at the exact tier:
`accuracy`, `selective_accuracy`, `coverage`. A backoff call is coarser than its
truth by construction, so exact scoring and `descendant_accuracy` give it no
credit. Backoff calls are counted on their own row, with `relaxed_accuracy`
and `ontology_distance_error` reported alongside.

**Expected, stated before measuring.** v1 called 1 of 28 clusters (coverage
0.016). The scoring diagnosis predicts that v3 calls more than 10 of 28 clusters,
exact or backoff, with exact-tier selective accuracy at or below the 0.636 ceiling.
If v3 calls 3 of 28 or fewer, the diagnosis is wrong and is reported that way.
The stop rules below apply unchanged.

### What v3 measured

Recorded 2026-10-01 from run `aida-pbmc-v11`, executed in a clean checkout of the
freeze commit `bad7532`. A first execution overlapped uncommitted edits to the
network modules. It was set aside, and its label table and proposals are
byte-identical to the clean run's. Result: `results/aida-pbmc-v11.cell-curator.json`.

| | v1 | v3 |
| --- | --- | --- |
| clusters labelled | 1 of 28 | 18 of 28 (3 by backoff) |
| cells called | 153 | 10,634 |
| per-cell coverage | 0.0088 | 0.611 |
| selective accuracy, exact | 0.961 (147/153) | 0.277 (2,942/10,634) |
| `relaxed_accuracy` (called cells; a coarser ancestor counts) | — | 0.874 |
| `descendant_accuracy` | — | 0.326 |
| `ontology_distance_error` | — | 0.165 |
| Brier / ECE, exact tier | — | 0.532 / 0.585 |

**The expectation held.** v3 called more than 10 of 28 clusters, and exact
selective accuracy is below the 0.636 ceiling. The diagnosis was right: the
scoring and the plumbing caused the 1-of-28 coverage, not caution.

**A stop rule fires, and it is reported rather than worked around.**
`descendant_accuracy` is 0.326, under 0.60, so the main limit is label taxonomy.
The vocabulary is a 13-way PBMC split. AIDA labels most T cells one or two levels
finer: the vocabulary says CD4-positive αβ T cell, the truth says naive or
central-memory CD4-positive αβ T cell. Every such call is exact-wrong and
relaxed-right, which is why the two tiers are 0.277 and 0.874. The exact-tier
calibration is poor for the same reason. Vote share measures how pure a cluster is
at the vocabulary's granularity, not whether the label matches AIDA's finer one.

**Errors, by cluster.**
- Cluster 14 (374 cells, truth naive CD8 T) was labelled "anatomical structure"
  (UBERON:0000061). This is a bug. HSC, T cell, and ILC programs tied, and their
  shared ancestors are all on the floor list, but the snapshot imports UBERON terms
  and backoff climbed into one.
- Cluster 12 (126 cells, truth NK) is called platelet at 0.94, as in v2.
- Cluster 25 (100 cells, truth CD141-positive mDC) is called CD1c-positive mDC, a
  sibling.

## Protocol v4: backoff stays inside the Cell Ontology

**Declared 2026-10-01, after v3 was seen and before v4 runs.** The single change:
a backoff ancestor must be a `CL:` term. Every other part of v3 is unchanged: input,
vocabulary, thresholds, lanes, and scoring. The change follows from the cluster 14
bug, not from tuning, and its only predicted effect is that cluster 14 becomes
Unknown. Any other difference between v3 and v4 means the change did more than it
says, and is reported as such. Run id `aida-pbmc-v12`.

## Why the pre-declared evidence bar is not met, measured rather than guessed

Two runs, differing only in measured panel width. Both label one cluster of 28.

| | 2,000 genes | 8,000 genes |
| --- | --- | --- |
| cells / clusters | 17,394 / 28 | 17,394 / 28 (clusters reused) |
| clustering ceiling | 0.636 | 0.636 (same clusters) |
| marker programs below the 0.5 coverage bar | 3 | 1 |
| clusters labelled | 1 of 28 | 1 of 28 |
| per-cell coverage | 0.009 | 0.009 |

Widening the panel helped exactly where it was predicted to and nowhere else.
Programs whose markers the panel could not measure became measurable, so the
coverage bar was met by more of the vocabulary. The number of labels did not move,
because the bar they are failing is not the coverage bar.

The gate they fail is the negative-marker one, and it is failing correctly. On 24
of 28 clusters the winning program's negative markers are flagged, and the
evidence shows why: the winning row for a *correct* call carries
`negative_signal = 0.0` with `negative_violation = 0.5`, and 0.5 is one violated
marker out of two. A cluster whose dominant true type is natural killer cell has
half its cells expressing CD3D, because the cluster is not half natural killer
cells. The ceiling says exactly this -- roughly a third of every cluster belongs to
a neighbouring lineage, so a neighbouring lineage's marker clears the 20% violation
threshold in most clusters. The gate is measuring cluster impurity and declining to
name a cluster that does not warrant the name.

**Resolution does not fix it, and this was measured rather than assumed.** Louvain
from 17 to 53 clusters moves the ceiling from 0.423 to 0.448. It is not a
granularity problem: the cells of different AIDA subtypes overlap in expression
space on this representation, so splitting further splits an impure cluster into
smaller impure clusters. The missing ingredient is donor integration across 394
donors, and this machine has no `scanpy`, `leidenalg`, `harmony`, or `scvi-tools`
available. Batch correction is explicitly outside this package's scope, so it
cannot be done by the tool under test and would have to be declared as input
preparation rather than performed silently.

One finding separates from the failure. The 8,000-gene slice clusters *worse*
than the 2,000-gene one on its own gene set: ceiling 0.44 against 0.64. The cause
is that log-normalisation depth depends on panel size, so a wider panel is a
different representation, not strictly more information. Clustering on the 2,000
most variable of the measured genes is therefore not a compromise, it is the
better input, and the wider panel's value is confined to marker measurement.

A second finding is about the vocabulary rather than the input.
Dispersion-based HVG selection systematically drops lineage-restricted genes --
CD34, TRDC, TRGC1, TRGC2, RORC, KIT, CRHBP -- so the cells it defines are exactly
the ones whose programs become uncallable. On this dataset hematopoietic stem cell
reaches 5 of 11 markers even at 8,000 genes. Widening a panel of highly variable
genes is the wrong lever for rare lineages; a panel not selected by dispersion
would be the right one.

## The input caps this evaluation before any annotator runs

The clustering is fixed a priori, not chosen: unit-variance scaling, 50 PCs, 15
nearest neighbours on the first 30, Louvain at resolution 1.0, seed 0, via the
`networkx` the package already depends on. It yields **28 clusters, median
purity 0.718 against AIDA truth, and a hard ceiling of 0.636 on per-cell
accuracy contributed by the clustering alone** — no annotator, however good, can
exceed it on this input.

The cause is a narrow gene panel plus 394 donors with no batch correction, on a
machine where `scanpy`, `leidenalg`, and `harmony` are unavailable. Resolution was
deliberately not tuned against purity, because tuning a clustering by its
agreement with the ground truth fits the input to the answer.

**Two corrections to the first draft of this section, both because a stated
constraint turned out to be false.** It said the gene panel could not be widened
because it would not fit in memory. That was doubly wrong: the 4 GB figure was a
budget the acquisition script set for itself, not a limit of the machine, which
has 3.9 TB of RAM and no cgroup cap. Nothing constrained the panel width. The
benchmark runs on 8,000 genes, with the cell set, the clustering recipe, the
marker vocabulary, and all four thresholds unchanged, so the 2,000-gene and
8,000-gene runs differ only in how many genes were measured and the comparison
between them is clean.

**This is recorded before the run and it governs how every number may be read.**
The literature's fine-PBMC figures of 52–75% sit *at or above* this ceiling, so
an absolute accuracy from this input is not comparable to any published number
and no such comparison is made. What survives is the *relative* comparison
between systems on identical cells, which is what the claim below is about and
which a shared ceiling depresses equally. Two consequences follow and both are
accepted rather than worked around: the per-cell-type numbers are properties of
this input as much as of the systems, and a system that abstains on clusters
below a usable purity may look better here than it would on full-transcriptome
data, because abstention is the correct response to an input that cannot support
the call.

Recording this is the protocol doing its job. The alternative — running the
annotation, seeing a number near 0.6, and explaining it afterwards as "PBMC is
hard" — is the failure this document exists to prevent.

## Scoring, frozen

Three granularity tiers are reported side by side, never one in place of
another:

| Field | Rule | Established precedent |
| --- | --- | --- |
| `accuracy`, `macro_f1`, `micro_f1` | exact string equality | — |
| `descendant_accuracy` | credit when the prediction is the truth or a **more specific** `is_a` descendant | scTab (Fischer et al. 2024) |
| `relaxed_accuracy` | also credit a **coarser** ancestor | — |
| `ontology_distance_error` | `(a+b)/(a+b+d(S))` over the `is_a` graph | — |

popV's exact/parent/child/sibling/no-match tiers are not implemented; the two
rules here cover its middle tiers and the distance metric separates its ends.

**Granularity is a first-order confound, not a footnote.** Published evidence:
CellTypist and Azimuth agree on 98.4% of lineage for 9,327 PBMC cells while
scoring ARI 0.661, attributed explicitly to "CellTypist classification labels
occupying an intermediate granularity"; one report saw Cohen's kappa fall from
0.895 to 0.27 through an unfixed label map alone. Any single number here is
close to meaningless without its tier.

**Always reported together, never separately:** `coverage`, `selective_accuracy`,
`coverage_risk`, `unknown_precision`, `unknown_recall`, `unknown_f1`, and the
full risk-coverage curve. This is not ceremony. Hu et al. found CellTypist with
an abstention threshold reached the *highest* cell-level agreement (0.917) and
the *worst* macro-F1 (0.540) because it answered on roughly one cell in five.
A selective-prediction claim without its coverage is the same defect.

`runtime_seconds` and `peak_memory_mb` are only published when measured by the
harness. They are absent from a result otherwise, and the shipped synthetic
result's constant values are a property of a hand-written fixture.

## Per-cell honesty

Annotation in this package is per cluster. A cell inherits its cluster's label,
confidence, and CL accession. The exported table carries the cluster key and the
receipt names the granularity, because a per-cell confidence column is otherwise
a per-cluster column wearing a per-cell name, and comparing it to a system that
varies within a cluster compares two different quantities.

A single-level run emits a real `predicted_path` and therefore
`hierarchy_error` is measured. A recursive run emits none, the column is absent,
and `hierarchy_error` is `None`. It is never reconstructed by inference.

## Expected results, stated before measuring

Recorded so the outcome cannot be rationalised afterwards. Fine-granularity human
PBMC is the hard case in the literature, and a macro-F1 of 0.5–0.6 on fine
subtypes is a normal respectable published result, not a failure.

Anchors on 22 immune subsets (SCimilarity, *Nature* 2025): SCimilarity 75.3%,
CellTypist 59.1%, scANVI 52.2%, TOSICA 44.4%. On HLCA at 61 types (Hu et al.
2025), per-cluster accuracy: Azimuth 69.3%, CellTypist 60.9%, scArches 71.7% —
with pooled F1 above 0.8 for all four tools. On scTab, 164 types cross-tissue,
CellTypist macro-F1 0.730.

The pooled-versus-per-cluster gap is the single largest trap in this field, so
per-cell-type F1 with support counts is mandatory, not optional. **T-cell
subtypes are the weakest class across every published method**, and the CD4 T
partition specifically is where Azimuth collapses 38.9% of cells into CD4 TCM
while leaving CD4 TEM nearly vacant.

## Comparators

Frozen as CellTypist, SingleR, and scANVI, each at a recorded version, each run
on the same slice through its own documented route, with `scANVI`'s reference
trained on the *same* data and seed and that fact stated. **Azimuth is excluded
from the primary comparison** for a reason that is about reproducibility rather
than fairness: no clean published Azimuth PBMC annotation accuracy exists, so
including it would make the headline number unreviewable.

## Stop rules

Decided now, so the evaluation terminates.

- **Harness invalid.** If the synthetic re-derivation in
  `tests/test_benchmark_instrument.py` fails, the instrument is wrong and no
  number from it is read.
- **Not comparable.** If `descendant_accuracy` against AIDA is below 0.60, the
  problem is label taxonomy rather than model quality. That is reported as the
  finding; the vocabulary is not retuned to chase the number.
- **Dominated.** If the risk-coverage curve is dominated by a comparator across
  its whole range, the selective-prediction claim is withdrawn. Not re-tuned.
- **Blocked.** If a bounded slice cannot be read within 4 GB, the run stops. No
  silent downsampling presented as the dataset.

## Where a novel result could still come from

Recorded so the work is aimed rather than merely done. Novelty is not in
reproducing a leaderboard; it would have to come from one of:

1. **A three-tier granularity stack with abstention reported beside it** on one
   cohort, with per-cell-type support counts. Hu et al. showed pooled F1 above
   0.8 while per-cluster accuracy was 60–70%; reporting both is the cheapest
   available improvement over existing comparisons.
2. **A disagreement analysis on the CD4 T partition**, where the published
   evidence says every method fails and where the granularity argument is
   strongest. This is where a real result is most likely to be informative.
3. **Coverage-matched comparison**, which no cited paper performs: hold coverage
   equal across systems and compare risk, rather than comparing numbers taken at
   each system's own preferred operating point.

None of these is claimed in advance.

## What this protocol does not do

No GPU result. The RAPIDS lane is skipped without a contract-satisfying device
and CI states plainly that a green run says nothing about it; the record lives in
a gitignored directory and its last entry is a pre-fix traceback. No GPU number
is published until it is captured on a real node and committed somewhere
durable.

No multi-dataset claim. One dataset, one cohort, one protocol version.

No superiority claim over any tool. The comparators are measured on the same
data by the same harness, and the claim is about matching at matched coverage.

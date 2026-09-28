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

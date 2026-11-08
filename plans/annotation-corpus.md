# Annotation corpus and evidence workflow

Status: implementation contract, started 2026-11-07.
The corpus core and compiled workflow are implemented in this checkout. See the
[user guide](../skills/cell-curator/references/annotation-corpus.md) for executable
commands and current limits. This document does not approve a biological corpus
or a study analysis.

## Intended behavior

Make the analyses needed for careful manual cell annotation repeatable. Collect
marker and classifier evidence, compare it with documented label definitions,
and resolve labels through explicit rules. Preserve uncertainty and evidence for
rejected labels. Reuse reviewed knowledge across studies without making an
earlier prediction into its own validation.

The corpus is a collection of source-backed label definitions. It is more than
a gene-list database. Allen's Immune Health Atlas is an example of the desired
depth, not a gold standard, universal taxonomy, or default preferred reference.
No source gains priority merely because it was imported first.

## What a label record needs

Each record needs a stable identifier and revision. Names alone are not keys.
Keep the following information as structured fields with readable descriptions:

| Area | Required information |
| --- | --- |
| Meaning | Preferred label, expanded name, biological description, identity or state axis, synonyms with their source |
| Hierarchy | Source hierarchy and source level; canonical parent links; explicit mappings to the study hierarchy |
| Scope | Species, tissue, age or developmental stage, condition, assay, feature namespace, applicable populations and exclusions |
| Positive evidence | Defining markers, supporting markers, programs, combinations, and the assay that measures each feature |
| Contrary evidence | Exclusion markers, nearest competing identities, discriminating features, and known failure cases |
| Decision criteria | Source-reported thresholds and their units; separate reviewed rules for the current study; required versus optional evidence |
| Provenance | Source record, URL or DOI, source version, retrieval date, content digest, reuse terms, and original label |
| Review | Imported, proposed, reviewed, disputed, superseded, or retired status; reviewer, rationale, revision links |
| Experience | Supporting and contradictory run receipts, sample and donor support, context, validation status, and known limitations |

Separate established marker knowledge from the markers observed in a particular
dataset. A source that supplies positives only has unknown negatives. Do not
invent exclusions, thresholds, or ontology identifiers to complete a record.
Store missing information explicitly.

Keep gene expression, measured protein, inferred activity, and transcript isoform
claims distinct. A decision rule must name its required measurement. A display
name such as `PD1+TIGIT+` needs declared features and positivity criteria; the
plugin must not infer the measurement from the text alone.

## Multiple sources and mappings

Support ontology definitions, atlas descriptions, curated marker resources,
published studies, reference model metadata, and reviewed local knowledge.
Preserve each source's assertion before proposing a common interpretation.

An adapter imports records with provenance and a receipt. It must expose partial
imports, omitted fields, retrieval failures, and source reuse limits. Importing
a label does not make it callable or reviewed. Keep downloaded content under
the existing network policy; annotation reads a pinned local snapshot.

Record mappings as exact, broader, narrower, related, conflicting, or unresolved.
Exact equivalence needs evidence. A similar label string or the same broad Cell
Ontology ancestor does not establish equivalence. Keep disagreements visible.

Source resolution levels are local to that source. A source's L3 must not become
the plugin's L3 automatically. Some sources combine state and identity in one
label. Preserve their label, then review the decomposition into identity and
state. Do not discard such records or silently reinterpret them.

Select applicable records from species, tissue, assay, stage, condition, and
the study's biological question. Record every selected source and exclusion.
Unknown context must remain unknown. Do not select a healthy PBMC reference as
a universal reference for other tissues or conditions.

Two websites copying one study are one underlying evidence source. Multiple
models trained on overlapping data do not create independent biological support.
Track source and training-reference dependencies when that information exists;
otherwise record independence as unknown.

## Evidence collected for an annotation

1. Inspect the input and freeze its identity, cluster membership, assay meanings,
   count and ranking matrices, sample structure, and available embeddings.
   Profile measured features before selecting applicable corpus records.
2. Compute Wilcoxon and logistic-regression marker evidence for each relevant
   partition. At finer levels, compare within the frozen parent. Include targeted
   sibling comparisons when candidate identities need them. Preserve effects,
   detection fractions, ranks, adjusted P values where available, and model
   convergence. Logistic-regression coefficients are not P values.
3. Profile positive and negative programs using the existing marker machinery.
   Evaluate competing identities, mixed populations, cell-cycle and other state
   programs, QC, and donor or batch concentration. Missing features remain
   unmeasured. These are descriptive cluster analyses, not condition-level DE.
4. Run applicable automated methods, including the requested CellTypist and popV
   lanes. Pin their models, references, preprocessing, versions, seeds, and
   hardware. Preserve per-cell predictions, cluster composition, uncertainty,
   and disagreements. A missing requested method is a visible unmet requirement,
   not a substitute result marked complete.
5. Assemble evidence cards and measured plots. Resolve known labels through the
   frozen rule set, record alternatives and failures, validate hierarchy nesting,
   and prepare independent critic review and the write-back preview.

The configured evidence contract determines which methods are required for a
run. Comprehensive coverage means the relevant analyses are accounted for;
unavailable or unsuitable methods need explicit reasons. It does not mean every
available classifier must run on every dataset.

popV can include CellTypist internally. Retain its constituent predictions and
method list. Do not count an internal CellTypist prediction and an identical
external prediction as independent votes. Wilcoxon, logistic regression, and
program scores also use overlapping data; their agreement is corroboration,
not independent replication.

Cluster-derived evidence stays identified as cluster evidence when propagated
to cells. Per-cell assignments need their own evidence and membership checks.
Mixed clusters must not receive a fabricated homogeneous cell-level result.
Subclustering is conditional on measured heterogeneity and existing refinement
gates; a target number of hierarchy levels is not a reason to manufacture splits.

## Hierarchy and deterministic decisions

The user wants L1 as broad lineage, L2 as an intermediate identity, and L3 as a
finer known identity. One illustrative path is:

```text
L1: T cell
L2: CD8 T cell
L3: CD8 T EMRA
L4 display: CD8 T EMRA PD1+TIGIT+
```

This is a requested output example, not a validated marker definition.

Maintain structured identity labels, parallel state calls, and a study-specific
extension. L4 may display their composition without changing the known parent
identity. A study extension can describe a supported local subtype or state;
it must have a stable local identifier and an accepted parent. Do not assign a
precise Cell Ontology identifier unless a matching term is validated. A broader
ontology mapping must be recorded as broader.

The rule set must define admissible candidates, measurement coverage, positive
support, negative contradictions, sibling separation, classifier reconciliation,
and abstention. Bind its version and thresholds before final label assignment.
Score missing evidence separately from evidence against a candidate. Preserve
the highest supported parent when a child cannot be resolved; use `Unknown`
when no identity is supported. Ties remain unresolved unless an authorized,
recorded tie-break rule succeeds.

Deterministic resolution means identical frozen evidence, corpus, mappings, and
rules produce identical labels and reasons. Record stochastic upstream methods
and their execution details separately. Pinning a seed does not establish
identical results across arbitrary hardware or package versions.

Existing LLM restrictions continue to apply. LLMs can help propose corpus
additions and explanations. LLM label proposals can only break qualified
near-ties through `evidence.llm.proposals`; they cannot rewrite measured calls.

## Learning across studies

Learning means accumulating reviewed, traceable knowledge. It does not imply
automatic classifier retraining or acceptance of the plugin's predictions.

```text
new source or reviewed run
  -> proposed addition or correction
  -> provenance, scope, and conflict checks
  -> reviewer accepts, rejects, or retains a dispute
  -> new corpus revision
  -> later runs may explicitly pin that revision
```

Keep candidate changes separate from the released corpus. Store reviewer
corrections and negative cases as well as accepted calls. Repeated labels from
the same donor, dataset, or reference are not new independent confirmations.
A study-specific observation remains local unless wider applicability receives
review. The workflow can suggest additions; the default is that a suggestion is
pending review and cannot affect production calls.

Every run pins a corpus revision and the exact records used. Corpus updates
leave completed runs unchanged. Reannotation produces a new run with a label
diff and reasons. Support superseding or retiring a bad rule without deleting
the original record or its history.

Test corpus changes on fixed, independently annotated fixtures or held-out
datasets. Keep evidence used to build the corpus out of its independent test
set. A deterministic replay test establishes repeatability only; it does not
establish biological accuracy. Track changed calls, new abstentions, finer
resolution, and regressions before accepting a revised release.

## Starting package coverage

Checked against the local source on 2026-11-07:

| Area | Current foundation | Work still needed |
| --- | --- | --- |
| Knowledge | `MarkerProgram`, local and cached marker providers, checksum-pinned remote tables | Rich label records, per-assertion provenance, source mappings, corpus revisions and review lifecycle |
| Markers | CPU Mann-Whitney evidence; native GPU Wilcoxon and logistic regression; signed programs | A uniform requested-method contract and reviewed execution path for both requested marker methods on supported hardware |
| Classifiers | CellTypist and several reference cross-check adapters | popV adapter, method-dependency tracking, and explicit reconciliation across requested classifiers |
| Hierarchy | Parent-scoped execution, retained parents, parallel states, critic and write-back gates | Corpus-driven candidate selection, reviewed source-to-study mappings, and a study extension presentation |
| Learning | Immutable runs and evidence receipts | Proposed corpus changes, reviewer actions, release history, and independent regression evaluation |

These foundations should be extended in place. Do not replace the existing
evidence engine, add a second approval system, or change the current source
preservation and human write-back contracts.

## Implementation order and checks

1. Add and validate the corpus format, local loading, explicit mapping relations,
   revision binding, and export to existing marker programs. Check duplicate IDs,
   bad parent links, cycles, ambiguous aliases, contradictory rules, unsupported
   measurements, and unreviewed records entering a callable vocabulary.
2. Add bounded source imports for a small reviewed slice from multiple sources.
   Start with CAP and Allen examples without privileging either. Preserve source
   identities, incomplete records, and conflicting assertions. Test imports with
   saved fixtures; keep network access out of normal annotation and unit tests.
3. Complete the requested marker-method contract and classifier adapters. Check
   both methods and parent scopes, logreg convergence, model preprocessing,
   feature overlap, barcode alignment, and duplicate-source accounting.
4. Connect the corpus to deterministic resolution and study extensions. Check
   exact replay, parent nesting, unknown labels, classifier disagreement,
   unmeasured negatives, and multi-state composition. Keep cluster propagation
   distinguishable from per-cell inference in every output.
5. Add the reviewed learning lifecycle and independent evaluation. Check that
   proposals cannot alter existing runs, corrections create new revisions, and
   a new revision produces a readable change report before adoption.

The first corpus slice is an implementation and review exercise, not a claim
that all immune labels are covered. Wider coverage needs source-specific work
and qualification. Per-run Gate A and human write-back remain mandatory under
the repository's existing rules.

## Sources inspected

- [CAP annotation metadata](https://celltype.info/docs/cell-annotation-metadata-terms):
  distinguishes author labels, ontology mappings, canonical markers, measured
  marker evidence, and rationales.
- [CAP API client](https://celltype.info/docs/python-client-for-cap-api):
  documents programmatic access to annotation metadata. API integration has not
  been exercised here.
- [Allen cell type descriptions](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/cell-type-descriptions/):
  one example of detailed hierarchical definitions and associated marker plots.
- [Allen CD8 descriptions](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/cell-type-descriptions/cd8-t-cells-gdt-cells-and-mait-cells/):
  illustrates source-specific subdivision and state-bearing labels. Apparent
  inconsistencies between table and prose need review during import.
- [Allen labeling models](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/downloads/models/):
  supplies classifiers at several levels with assay-specific variants.
- [QuadBio marker recipes](https://github.com/quadbio/cell-type-annotation/blob/main/skills/cell-type-annotation/references/marker_tools.md):
  documents parent-local marker analysis and multivariate marker checks.
- [QuadBio method selection](https://github.com/quadbio/cell-type-annotation/blob/main/skills/cell-type-annotation/references/method_selection.md):
  its one-primary-cross-check policy is distinct from this user's requested
  CellTypist-plus-popV evidence workflow.
- [CellAnnotator](https://github.com/quadbio/cell-annotator/blob/main/README.md):
  an example of LLM annotation proposals; it is not the deterministic resolver.
- [popV documentation](https://popv.readthedocs.io/en/latest/):
  identifies constituent classifiers, including CellTypist, and method outputs.

Source reuse terms and exact import endpoints still need verification. No atlas
corpus, reference model, or study data was downloaded by writing this plan.

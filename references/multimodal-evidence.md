# Assay-aware evidence contract

## Binary marker lanes

For every candidate cross each configured evidence lane with candidate versus
parent remainder and candidate versus closest biological/sibling comparator.
Group 0 is candidate and group 1 is comparator.

Each RAPIDS lane must run native Wilcoxon and balanced binary logistic
regression. Do not substitute a multiclass one-versus-rest coefficient.
Non-RAPIDS lanes require a declared backend and immutable receipt.

## RAPIDS receipt

Validate exactly one allowed GPU, model/compute capability, RAPIDS distribution,
`rapids-singlecell` version, CUDA/driver, native API, both methods, cell/output
hashes, and explicit false CPU-fallback booleans. Validate before data read or
output creation.

## Absolute and relative evidence

Join candidate/comparator means, valid detection fractions when available,
signed effects with semantics, adjusted P values, method concordance,
positive-program coverage, and expected-negative violations. A signed
representation may contribute means/effects. Detection requires a declared
nonnegative matrix; otherwise report it unavailable rather than thresholding.

## Program and cell-level evidence

Score identity, subtype-positive, expected-negative, and incompatible programs
only on compatible feature spaces. Aggregate markers cannot distinguish separate
singlets from within-cell co-occurrence, so report per-cell distributions and a
predeclared separation statistic. Integrate available doublet evidence without
inventing it for studies lacking scores.

## Technical evidence

Report configured sample/capture representation, dominance/effect sizes, and
technical metrics. These may include mitochondrial fraction, RNA/ATAC/protein
complexity, TSS enrichment, FRiP, nucleosome signal, antibody complexity,
sex-linked programs, batch, or spatial/capture effects.

## Source corroboration

- PanglaoDB ORA: gene or gene-activity rankings only, with measured background.
- Peak ATAC: declared peak, peak-to-gene, or motif enrichment where available.
- Protein: curated panel/gate evidence with antibody and panel limitations.
- Reference classifiers: hypotheses only; record confidence and mismatch.

Source corroboration never replaces direct positive and negative evidence.
Panel-missing features are unavailable, not negative.

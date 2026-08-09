# Validation and completion

## Structural and scientific

- skill validator, package compilation/build, tests, and `--help` checks pass;
- reusable logic contains no study paths, IDs, marker panels, or scheduler resources;
- RNA-only, ATAC-only, and multimodal input contracts are tested;
- absent/mismatched GPUs fail RAPIDS lanes with no CPU path;
- arbitrary evidence-lane counts create both binary comparators;
- optional detection is unavailable and required signed detection blocks;
- stored/local candidate, one-cluster rejection, technical/doublet, positive and
  expected-negative, and evidence-only contracts are tested;
- every non-`Unknown` label has an executable evidence route; CL-free labels remain
  valid when ontology evidence is unavailable and that limitation is recorded.

## Mapping

Require one row per input cell, no cross-parent reassignment, parent count
conservation, unchanged unsplit membership, disjoint accepted children, failed
candidate fold-back, accepted-candidate leaf IDs, and correct scope/global totals.

## Immutability and resume

Input, canonical partition, completed packet, and reviewed-annotation hashes
remain unchanged. Partial valid outputs resume; changed configuration creates a
new run. Manifests contain input/config/code/output and backend provenance.

## Forward tests

Use realistic prompts without expected labels on RNA-only, peak-ATAC-only,
RNA-plus-protein or RNA-plus-ATAC, a local biological split, and a mixed/doublet
parent. Compare artifacts against invariants rather than expected labels.

Complete only when every parent has one outcome, accepted children pass all
applicable gates, every cell maps once, totals conserve, critics reconcile,
hashes remain invariant, and human approval is recorded.

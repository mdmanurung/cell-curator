# Artifact schemas

TSV files use UTF-8, tab separators, one header, stable row order, explicit
empty strings, and no implicit index. JSON manifests are sorted and newline-
terminated. Every artifact carries `schema_version`, `run_id`, and lineage or
scope when applicable.

## Control documents

- `context.md`: dataset purpose, inputs, assays/capabilities, canonical partition,
  environment, and declared limitations.
- `assumptions.md`: Gate A living ledger and sign-off.
- `plan.md`: phase DAG, resources, review order, and stop conditions.
- `run_state.json`: mode, phase state, blockers, immutable hashes, and timestamps.

## Parent and candidate tables

- `adaptive_parent_audit.tsv`: one row per canonical parent; counts, guide
  entropy/minorities, incompatible programs, doublets, donor/capture dominance,
  QC/complexity/sex/batch metrics, flags, reasons.
- `candidate_registry.tsv`: one row per candidate; candidate/parent IDs, source
  kind, source resolution/cluster, membership hash, counts, nesting, stability,
  closest comparator, triage action, and blank label/CL fields.
- `candidate_membership.tsv`: one row per candidate per parent cell with Boolean
  membership. Unique key is `(candidate_id, barcode)`.
- `parent_decision_summary.tsv`: exactly one row per parent with outcome,
  accepted child IDs, reasons, alternatives, confidence, reviewer state.
- `parent_child_mapping.tsv`: exactly one row per input cell with canonical
  parent, final leaf ID, membership role, evidence outcome, and review status.
- `child_evidence_summary.tsv`: one row per tested candidate with every gate,
  failed reasons, evidence paths, confidence, and blank/final label fields.

## Comparator and marker artifacts

- comparator manifest: one row per candidate x comparator role x evidence lane;
  expected groups, cell hashes, input path, output paths, and status.
- marker tables: candidate/comparator, assay kind, feature space, representation,
  backend/method,
  feature, rank, score/effect semantics, P/adjusted P, target/reference levels
  and valid detection fractions.
- RAPIDS receipt: validated hardware/software contract, method parameters,
  convergence, hashes, and explicit fallback booleans.

## Recursive hierarchy invocation

Pass `run execute --score-manifest FILE` a strict JSON object with
`schema_version: 3` and a non-empty `score_tables` array. Each entry contains
exactly `level`, `parent_scope`, and `path`, plus optional `state_path`. Paths are
absolute or relative to the manifest directory. Optional top-level fields are
`vocabulary_path`, `minimum_confidence`, and `minimum_margin`.

Every score TSV contains `cluster_id`, `label`, and finite numeric `score`; optional
`specificity`, `agreement`, and `coverage` columns feed advisory confidence. Every
state JSON maps cluster IDs to structured parallel-state values. Duplicate scopes,
unknown fields, invalid thresholds, absent files, cross-parent labels, and stale
artifacts are hard errors.

The executor copies scoped inputs into `hierarchy/scores/` and `hierarchy/states/`,
writes scoped labels below `labels/scopes/`, consolidates sibling parents into one
durable `labels/L*.tsv`, writes retained leaves, and hashes each scope checkpoint.
Missing child scopes retain their accepted parent leaf and never reuse sibling
evidence.

## Review packet

Include evidence cards, parent decision cards, parent/full embeddings,
assay-appropriate plots, sample/capture/QC/doublet summaries,
persistence/alluvial tables, critic findings,
reconciliation, blank review template, mapping, leaf totals, and provenance.

An evidence card states hypothesis, direct evidence, negative evidence,
alternatives, technical/doublet assessment, stability, guide divergence,
confidence margin, gate failures, and links to raw tables.

## Leaf totals

Report canonical parents, accepted candidates, remainder leaves, effective
leaves by lineage, and global effective leaves. The effective leaf count is the
number of unique final leaf IDs, not canonical parents plus accepted-parent
count.

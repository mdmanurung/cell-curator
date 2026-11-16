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

## Mapping table

- `parent_child_mapping.tsv`: exactly one row per input cell with canonical
  parent, final leaf ID, membership role, evidence outcome, and review status.

## Comparator and marker artifacts

- comparator manifest: one row per candidate x comparator role x evidence lane;
  expected groups, cell hashes, input path, output paths, and status.
- marker tables: candidate/comparator, assay kind, feature space, representation,
  backend/method,
  feature, rank, score/effect semantics, P/adjusted P, target/reference levels
  and valid detection fractions.
- RAPIDS receipt: validated hardware/software contract, method parameters,
  convergence, hashes, and explicit fallback booleans.

## Consensus import and hierarchy execution

`annotation import-consensus` reads the run named by `annotation.consensus_directory`:
`consensus.tsv`, `run_manifest.json` and `work/cells.tsv`. `consensus.tsv` needs the
columns `cluster_id`, `label`, `cl_id`, `decision_rule`, `rationale`,
`rejected_alternatives`, `cluster_cell_share` and `cluster_cell_lead`. It writes
`annotation/levels/L1/labels.tsv`, copies the three files to `annotation/levels/L1/consensus/`,
and publishes `decision.manifest.json` with every input hash and `decided_by:
annotation-consensus`.

`finalize` passes those labels to the hierarchy executor as one decision table per
`(level, parent_scope)`. The executor chooses no label. It copies scoped inputs into
`hierarchy/scores/` and `hierarchy/decisions/`, writes scoped labels below
`labels/scopes/`, consolidates sibling parents into one durable `labels/L*.tsv`,
writes retained leaves, and hashes each scope checkpoint, so a changed decision
re-executes its scope. A decision label must be `Unknown` or an identity label the
vocabulary lists for that scope; a cell state program is refused. Only the vocabulary
rows of the configured levels need an approved assumption.

## Review packet

`review_packet/` is the immutable snapshot critics review, frozen when execution
completes and checked only against its own manifest. It includes evidence cards,
parent decision cards, parent/full embeddings, assay-appropriate plots,
sample/capture/QC/doublet summaries, persistence/alluvial tables, the critic packet,
blank review template, mapping, leaf totals, and provenance. The critic findings and
reconciliation are added afterwards under `review/`; the write-back approval token
binds them, together with the rebuilt report.

An evidence card states hypothesis, direct evidence, negative evidence,
alternatives, technical/doublet assessment, stability, guide divergence,
confidence margin, gate failures, and links to raw tables.

## Leaf totals

Report canonical parents, accepted candidates, remainder leaves, effective
leaves by lineage, and global effective leaves. The effective leaf count is the
number of unique final leaf IDs, not canonical parents plus accepted-parent
count.

## Critic ledgers and write-back receipt

`review/critic_findings.tsv` is append-only: an import adds findings, and re-importing
an existing `finding_id` with different content is refused. Every pass named in
`review.critic_passes` must contribute at least one row; a pass that found nothing
records one row with severity `NONE`. Every row needs a terminal reconciliation
(`RESOLVED` or `DISMISSED`) with a reviewer and resolution.

`writeback/receipt.json` is the preview plus `output_sha256` and an `approval`
object: `approval_mode` (`interactive_confirmed` or
`preauthorized_non_interactive`), `reviewer`, `os_user`, `host`, and
`approved_utc`. The preview's `approval_token` is a 12-character digest over every
other preview field. An interactive `write-back apply` writes only after that token
is typed back, so a change to labels, assumptions, report, critic ledgers, diff, or
target invalidates a confirmation given earlier.

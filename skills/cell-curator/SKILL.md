---
name: cell-curator
description: Run independent, local-first, evidence-driven annotation of clustered AnnData or MuData with guided context confirmation, hierarchical parent-scoped labels, multimodal marker and reference cross-checks, adaptive refinement, human gates, critic reconciliation, self-contained reports, and guarded write-back. Use for offline or optionally enhanced single-cell and spatial cell-type annotation workflows.
---

# Cell Curator

Run `cell-curator` independently. Keep local AnnData or MuData, marker files,
ontology files, and labeled references sufficient for an offline run. Treat optional
classifiers, atlases, foundation models, downloads, and web services as corroborating
providers. Never require QuadBio, an MCP, or a particular agent platform.

## Preserve the invariants

1. Freeze the source hash, canonical working copy, and one parent membership per
   cell. Never mutate the source during computation.
2. Constrain every child vocabulary, score, refinement, and checkpoint to its
   parent. Preserve frozen parents and fold rejected children into the correct
   remainder.
3. Keep `cell_type` and `cell_state` separate. Preserve marker, reference,
   classifier, modality, and spatial disagreement. Accept `Unknown` as success.
4. Treat missing panel features as unmeasured. Use a declared nonnegative matrix
   for detection; use signed or corrected matrices only for valid means and effects.
5. Validate a real GPU before every GPU route. Fail closed without CPU substitution.
6. Require approved assumptions before biological labels. Record the same ledger,
   decisions, checkpoints, and provenance in interactive and non-interactive runs.
7. Keep labels tables durable. Reconcile critics and preview the complete diff before
   approved atomic write-back.
8. Let Snakemake own dependencies, staleness, and resume when executing the bundled
   workflow.

## Configure the run

Copy [config.template.yaml](assets/config.template.yaml). Declare the input, cluster
and scope keys, embedding and visualization keys, assays, layers, representations,
detection sources, organism, tissue, developmental stage, condition, technology,
experimental context, composition covariates, vocabulary, references, compute
limits, hierarchy levels, and write-back prefix.

Use `source: raw` for a declared `.raw` representation and `gene_mapping_path` for
a local Ensembl-to-symbol mapping when automatic symbol detection is insufficient.

Register local signed marker JSON/TSV, local ontology, labeled AnnData/MuData
references, or cached exports under `knowledge.providers`. Permit CL-free labels
when ontology evidence is absent; validate every supplied CL identifier when it is
present.

Read only the direct reference needed:

- Input, IDs, matrices, and QC: [input-contract.md](references/input-contract.md)
- Assay and provider routing: [assay-routing.md](references/assay-routing.md)
- Parent-local refinement: [adaptive-subclustering.md](references/adaptive-subclustering.md)
- Signed and multimodal evidence: [multimodal-evidence.md](references/multimodal-evidence.md)
- Assumptions, decisions, and critics: [decision-gates.md](references/decision-gates.md)
- Durable artifacts: [artifact-schemas.md](references/artifact-schemas.md)
- Compute and resume: [compute-and-snakemake.md](references/compute-and-snakemake.md)
- Strict completion: [validation.md](references/validation.md)
- Neutral capability mapping: [quadbio-parity.md](references/quadbio-parity.md)

Validate the configuration before loading data:

```bash
cell-curator config validate --config config.yaml
cell-curator environment detect --config config.yaml
cell-curator data inspect --config config.yaml
cell-curator ontology provision --config config.yaml
```

Stop on an ambiguous representation, invalid normalization fingerprint, duplicate
feature identity, missing cluster/modality/embedding contract, mismatched cells,
source-hash drift, scheduler ambiguity, or an unvalidated requested GPU.

## Stage 1: Set the scene

Create the read-only inspection, canonical working copy, QC metrics, environment
record, context, assumptions ledger, progress state, and compute profile:

```bash
cell-curator data canonicalize --config config.yaml
cell-curator data qc --config config.yaml
cell-curator run scene --config config.yaml
cell-curator progress status --config config.yaml
```

For an interactive run, review `context.json`, correct every unresolved biological
or write-back fact in a JSON response, and record the reviewer:

```bash
cell-curator run confirm-context --config config.yaml \
  --input context-confirmation.json --reviewer REVIEWER
```

For a pre-authorized non-interactive run, set `guidance.interactive: false`,
`guidance.preauthorized: true`, a named reviewer, complete `guidance.context`, and
recorded `guidance.assumptions`. Do not use non-interactive mode to bypass a missing
context or gate.

## Stage 2: Build and approve the plan

Assemble a parent-scoped vocabulary, profile measured panels, compute cluster
markers, score signed programs, and select one primary cross-check per hierarchy
level:

```bash
cell-curator markers assemble --config config.yaml
cell-curator markers compute --config config.yaml --assay transcriptome \
  --representation logcounts --group-key canonical_cluster
cell-curator markers profile --config config.yaml --assay transcriptome \
  --representation logcounts --group-key canonical_cluster
cell-curator markers score --config config.yaml --assay transcriptome \
  --representation logcounts --group-key canonical_cluster
cell-curator route select --request route-request.json
cell-curator run plan --config config.yaml
cell-curator assumptions list --config config.yaml
```

Resolve assumptions before labeling. Add newly discovered material assumptions and
approve only the relevant scope:

```bash
cell-curator assumptions add --config config.yaml --scope SCOPE \
  --statement "MATERIAL ASSUMPTION"
cell-curator assumptions approve --config config.yaml \
  --scope SCOPE --reviewer REVIEWER
```

Let a valid non-interactive configuration record pre-authorization automatically.
Require renewed approval only when a later level introduces a materially new
assumption.

## Stage 3: Execute recursively

Execute one level within one parent, checkpoint it, and continue only into accepted
children. Inspect stored higher resolutions before one bounded parent-local
subclustering pass. Profile ambiguity, technical association, low quality, cycling,
stress, gradients, mixed/doublet state, composition, and spatial context before
assigning a child or `Unknown`.

```bash
cell-curator run execute --config config.yaml
cell-curator cluster subcluster --config config.yaml
cell-curator cluster profile --config config.yaml
cell-curator cluster confidence --margin 0.3 --specificity 0.8 \
  --agreement 0.7 --coverage 0.9
cell-curator hierarchy assign --config config.yaml --scores scores.tsv \
  --level L1 --parent ROOT
cell-curator hierarchy validate --labels results/RUN/labels/L1.tsv
cell-curator run status --config config.yaml
```

For an arbitrary-depth vocabulary, list each real parent-scoped score table in a
schema-v3 JSON manifest and execute it through the same gate and checkpoint engine:

```json
{
  "schema_version": 3,
  "score_tables": [
    {"level": "L1", "parent_scope": "ROOT", "path": "scores/L1-ROOT.tsv"},
    {"level": "L2", "parent_scope": "Immune", "path": "scores/L2-Immune.tsv"}
  ]
}
```

```bash
cell-curator run execute --config config.yaml \
  --score-manifest hierarchy-scores.json
```

Resolve manifest paths relative to the manifest. Retain an accepted parent as a
leaf when its child score table is absent. Never borrow a sibling parent's scores.
See [artifact-schemas.md](references/artifact-schemas.md) for the strict contract.

Use `route select` to choose the primary cross-check. Invoke only applicable routes;
inspect every `unavailable` receipt and remediation. Discover or run optional routes
with their own help contracts:

```bash
cell-curator crosscheck obs-pred --help
cell-curator crosscheck celltypist-discover --help
cell-curator crosscheck celltypist-predict --help
cell-curator crosscheck reference-knn --help
cell-curator crosscheck azimuth --help
cell-curator crosscheck scimilarity --help
cell-curator reference acquire-census --help
cell-curator reference train-scanvi --help
cell-curator reference map-scarches --help
```

Treat all automated outputs as hypotheses. Keep CPU reference kNN available; permit
GPU kNN, scVI/scANVI training, or scArches mapping only under their validated GPU
contracts. Bound and checksum every acquisition.

## Review, reconcile, and validate

Scaffold evidence cards from raw values, preserve manual reasoning in labels tables,
build a platform-neutral critic packet, import findings, and reconcile each finding
as `RESOLVED` or `DISMISSED` with reviewer and rationale:

```bash
cell-curator evidence-card scaffold --config config.yaml
cell-curator critic packet --config config.yaml
cell-curator critic import-findings --config config.yaml --input critic-findings.tsv
cell-curator critic reconcile --config config.yaml --input critic-reconciliation.tsv
cell-curator critic validate --config config.yaml
cell-curator report build --config config.yaml --strict
cell-curator report check --config config.yaml --strict
cell-curator run validate --config config.yaml
```

Require the self-contained report to include embeddings, marker dotplots, hierarchy,
uncertainty, rejected alternatives, QC, spatial context when available, composition,
categorical colors, provenance, and all approved levels. Use report checking as a
read-only post-approval validation.

## Preview and write back

Generate the complete diff first:

```bash
cell-curator write-back preview --config config.yaml --output annotated.h5ad
```

After assumptions, hierarchy, reports, and critics pass, create the signed approval
for the current run plus the source, mapping, durable-label, assumptions, report,
critic, preview-diff, exact-output-path, and in-place-intent fields.
Write atomically to a new file:

```bash
cell-curator write-back apply --config config.yaml \
  --approval approval.json --output annotated.h5ad
```

Preserve category order and report colors. Mutate the source only through the
explicit in-place option. Finish by checking both status surfaces:

```bash
cell-curator progress status --config config.yaml
cell-curator run status --config config.yaml
```

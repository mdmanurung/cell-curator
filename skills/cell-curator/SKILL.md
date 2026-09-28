---
name: cell-curator
description: Annotate pre-clustered AnnData or MuData from local marker, ontology, reference, and multimodal evidence. Use for review-gated cell-type annotation, uncertainty assessment, parent-scoped refinement, and guarded write-back.
---

# Cell Curator

Use cell-curator to annotate an existing clustering. Local inputs are sufficient;
optional classifiers, atlases, and web resources provide corroborating evidence.

## Analytical constraints

1. Freeze the source hash, working copy, and parent membership before analysis.
2. Scope every child vocabulary, score, and refinement to its parent.
3. Keep cell type and cell state separate. Preserve disagreement and `Unknown`.
4. Treat missing features as unmeasured. Declare separate ranking and detection
   representations when their semantics differ.
5. Run GPU methods only after the declared hardware and software contract passes.
6. Approve material assumptions before assigning biological labels.
7. Reconcile critic findings and preview the complete diff before write-back.
8. Use Snakemake to manage dependencies and resume for workflow execution.

## References

Read only the contract needed for the current decision.

| Topic | Contract |
| --- | --- |
| Inputs, matrices, and QC | [input-contract.md](references/input-contract.md) |
| Assays and providers | [assay-routing.md](references/assay-routing.md) |
| Parent-local refinement | [adaptive-subclustering.md](references/adaptive-subclustering.md) |
| Multimodal evidence | [multimodal-evidence.md](references/multimodal-evidence.md) |
| Decisions and critics | [decision-gates.md](references/decision-gates.md) |
| Durable artifacts | [artifact-schemas.md](references/artifact-schemas.md) |
| Compute and resume | [compute-and-snakemake.md](references/compute-and-snakemake.md) |
| Completion criteria | [validation.md](references/validation.md) |
| Capability mapping | [quadbio-parity.md](references/quadbio-parity.md) |

Copy [config.template.yaml](assets/config.template.yaml). Declare the input, cluster
key, assays, representations, detection source, biological context, knowledge
providers, hierarchy, compute limits, and output path.

## 1. Validate and freeze

Validate the configuration and inspect all requested resources before analysis:

```bash
cell-curator config validate --config config.yaml
cell-curator environment detect --config config.yaml
cell-curator data inspect --config config.yaml
cell-curator ontology provision --config config.yaml
cell-curator data canonicalize --config config.yaml
cell-curator data qc --config config.yaml
cell-curator run scene --config config.yaml
cell-curator progress status --config config.yaml
```

Resolve ambiguous matrices, duplicate features, missing keys, cell mismatches,
source-hash drift, scheduler ambiguity, and invalid GPU contracts before proceeding.
Review `context.json`, then record corrections and the reviewer:

```bash
cell-curator run confirm-context --config config.yaml \
  --input context-confirmation.json --reviewer REVIEWER
```

A non-interactive run requires `guidance.interactive: false`,
`guidance.preauthorized: true`, a reviewer, complete context, and recorded
assumptions. It uses the same gates and ledger as an interactive run.

## 2. Build and approve the plan

Construct the parent-scoped vocabulary, profile the measured panel, compute markers,
score signed programs, and select a cross-check:

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

Record and approve material assumptions for the relevant scope:

```bash
cell-curator assumptions add --config config.yaml --scope SCOPE \
  --statement "MATERIAL ASSUMPTION"
cell-curator assumptions approve --config config.yaml \
  --scope SCOPE --reviewer REVIEWER
```

Renew approval when a later hierarchy level introduces a material assumption.

## 3. Generate evidence and assign labels

Execute each level within one parent. Assess technical association, quality, stress,
cycling, gradients, mixed states, composition, and spatial context before accepting
a child label or `Unknown`.

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

For multiple levels, pass `run execute` a schema-v3 score manifest. Each entry must
identify one `(level, parent_scope)` and its score table. Resolve paths relative to
the manifest, retain an accepted parent when no child table exists, and do not reuse
scores across parents.

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

Use only applicable cross-checks. Inspect `unavailable` receipts and their stated
remediation. Automated predictions remain hypotheses.

```bash
cell-curator crosscheck obs-pred --input canonical.h5ad \
  --prediction-column predicted_label --cluster-column canonical_cluster \
  --output-dir results/RUN/crosschecks
cell-curator crosscheck celltypist-discover --output-dir results/RUN/crosschecks
cell-curator crosscheck celltypist-predict --input canonical.h5ad \
  --model MODEL --cluster-column canonical_cluster \
  --output-dir results/RUN/crosschecks
cell-curator crosscheck reference-knn --query canonical.h5ad \
  --reference reference.h5ad --reference-label cell_type --backend cpu \
  --output-dir results/RUN/crosschecks
cell-curator crosscheck azimuth --query canonical.h5ad \
  --reference reference.h5ad --reference-label cell_type \
  --output-dir results/RUN/crosschecks
cell-curator crosscheck scimilarity --input canonical.h5ad \
  --model-path models/scimilarity --output-dir results/RUN/crosschecks
cell-curator reference acquire-census --organism homo_sapiens --tissue blood \
  --output references/census.h5ad --version CENSUS_VERSION --max-cells 50000 \
  --expected-sha256 SHA256
cell-curator reference acquire-azimuth --url https://HOST/PATH \
  --output references/azimuth.h5ad --expected-sha256 SHA256
cell-curator reference acquire-scimilarity --url https://HOST/PATH \
  --output models/scimilarity.tar.gz --expected-sha256 SHA256
cell-curator reference train-scanvi --reference references/census.h5ad \
  --labels-key cell_type --counts-layer counts --output-dir models/scanvi
cell-curator reference map-scarches --query canonical.h5ad \
  --model-dir models/scanvi --output-dir results/RUN/crosschecks
```

Add `--help` to a route for its remaining flags. A `.h5mu` input additionally
requires `--input-modality`, `--query-modality`, or `--reference-modality`.
`--backend gpu` on `reference-knn` runs only after GPU validation passes.

Each `acquire-*` route is bounded by `--max-bytes` and `--timeout-seconds`, is
pinned by a required `--expected-sha256`, and writes a receipt recording the
cache hit and checksum result. A bundle is stored as published; extract it
yourself before passing it to `--reference` or `--model-path`.

## 4. Review and validate

Create evidence cards, import independent critic findings, reconcile each finding,
and build the self-contained report:

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

The local ontology decides whether a CL identifier is valid. When a
`live_ontology` provider is configured and remote enhancement is authorized,
corroborate the accessions you assigned against it before write-back:

```bash
cell-curator ontology corroborate --config config.yaml \
  --term CL:0000084 --expected-name "T cell"
```

It reports agreement, disagreement, and unreachability, stores each answer it
received as a hashed artifact, and changes no label. An absent network is an
expected result, not a failure.

The report must show embeddings, marker evidence, hierarchy, uncertainty, rejected
alternatives, QC, composition, provenance, and spatial context when available.
Reconciliation records a reviewer and rationale for every resolved or dismissed
finding.

To score this run against an external ground truth, export its labels as a
neutral prediction table:

```bash
cell-curator run export-predictions --config config.yaml \
  --truth ground-truth.tsv --output predictions.tsv
```

The truth table must name exactly the barcodes the run labelled, carrying
`barcode` and `true_label` and optionally `true_cl_id` and `true_path`. Any
overlap that is not an exact match is refused rather than scored, so the cell set
being measured is one that was chosen. An immutable receipt beside the table
records the truth digest, the run's mapping digest, and every label table it was
built from. Annotation is per cluster, so the per-cell confidence column is a
cluster value carried to each cell, and the receipt says so.

## 5. Preview and write back

Preview the complete mapping:

```bash
cell-curator write-back preview --config config.yaml --output annotated.h5ad
```

After assumptions, hierarchy, report, and critic validation pass, create an approval
bound to the current run, source, labels, assumptions, report, findings, preview
diff, output path, and in-place intent. Apply it atomically to a new file:

```bash
cell-curator write-back apply --config config.yaml \
  --approval approval.json --output annotated.h5ad
cell-curator progress status --config config.yaml
cell-curator run status --config config.yaml
```

Preserve category order and report colours. Modify the source only through the
explicit in-place option.

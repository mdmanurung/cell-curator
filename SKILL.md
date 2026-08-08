---
name: adaptive-multimodal-cell-annotation
description: Orchestrate evidence-gated annotation and adaptive nested subclustering of pre-clustered single-cell or spatial AnnData/MuData. Use for RNA-only, protein-only, ATAC-only, gene-activity, CITE-seq, multiome, or other assay combinations when canonical clusters must stay frozen, heterogeneous parents need bounded refinement, evidence must be assay-aware, or immutable Snakemake/HPC review packets and ontology-guided human approval are required.
---

# Adaptive Multimodal Cell Annotation

Use this skill as the production-grade orchestration layer over
`$cell-type-annotation`. The agent remains the annotator. Tools, classifiers,
gates, and enrichment produce hypotheses and evidence.

## Hold these contracts throughout

1. Freeze the canonical parent partition and input hashes before analysis.
2. Create versioned outputs. Preserve source objects, canonical membership,
   completed packets, and reviewed annotations.
3. Make Snakemake authoritative for dependencies, staleness, phases, and resume.
   Run heavy computation in deterministic scripts; keep notebooks read-only.
4. Treat assays as configured capabilities, not a fixed RNA-plus-protein bundle.
   Require only scientifically applicable evidence and record unavailable
   capabilities as limitations or blockers.
5. Require a real validated GPU for each evidence lane whose backend is
   `rapids_singlecell`. Fail before data read or output creation when its
   GPU/software contract fails. Provide no CPU replacement for that lane.
6. Treat automated annotations and gates as overridable guides. Guide consensus
   is not biological evidence.
7. Keep `cell_type` and `cell_state` separate. `Unknown` is a successful outcome.
8. Write no label or CL ID until Gate A, evidence completion, critic
   reconciliation, Knowledgebase semantic validation, and explicit human
   approval pass.

## Select the operating mode

Verify the QuadBio Knowledgebase MCP before biological reasoning.

- `annotation`: the Knowledgebase is callable and at least one declared lane is
  ontology-compatible. Labels may enter review after all gates.
- `evidence-only`: compute label-free audits, markers, stability, technical
  summaries, evidence cards, and blank review templates. Record missing
  capabilities as blockers; leave labels and CL IDs blank and block write-back.

An assay can also force an annotation run into an effectively evidence-only
state. For example, peak-level ATAC without gene activity, peak-to-gene mapping,
or a validated reference may reveal stable subpopulations but cannot alone
justify gene-marker ontology labels.

## Route references by branch

- Assay selection and support matrix: [assay-routing.md](references/assay-routing.md)
- Input and freeze checks: [input-contract.md](references/input-contract.md)
- Candidate discovery/local splitting: [adaptive-subclustering.md](references/adaptive-subclustering.md)
- Relative, absolute, and source evidence: [multimodal-evidence.md](references/multimodal-evidence.md)
- Parent/child decisions and critics: [decision-gates.md](references/decision-gates.md)
- Snakemake, HPC, backends, and immutability: [compute-and-snakemake.md](references/compute-and-snakemake.md)
- Artifact column contracts: [artifact-schemas.md](references/artifact-schemas.md)
- Completion and invariance checks: [validation.md](references/validation.md)
- QuadBio reuse and extensions: [quadbio-integration.md](references/quadbio-integration.md)
- Development parity and implementation backlog: [quadbio-parity.md](references/quadbio-parity.md)

## Workflow

### 1. Establish Gate A

Invoke `$cell-type-annotation` and use its living assumptions ledger. Copy
[assumptions.template.md](assets/assumptions.template.md) into the versioned run.
Record biology, vocabulary, assay limitations, feature spaces, missing
capabilities, canonical partition, covariates, QC/doublet keys, guide roles,
thresholds, compute environment, and operating mode. Obtain explicit sign-off
before proposing biological labels.

Complete when assumptions are signed and inputs/configuration are frozen.

### 2. Declare and inspect assay capabilities

Start from [config.template.yaml](assets/config.template.yaml). Register each
assay and each evidence representation. For every representation declare:

- assay kind and feature space;
- source matrix/layer and feature identifiers;
- ranking backend;
- whether it supports ontology-marker reasoning;
- an optional nonnegative detection source and whether detection is required;
- whether the lane is required for child acceptance.

Run:

```bash
python scripts/inspect_input_contract.py \
  --config config.yaml --output-root RUN_ROOT
```

Stop on missing or ambiguous matrices. Do not infer layers from names. A signed
representation may support means/effects, but detection must come from a
separately declared nonnegative source. Detection may be unavailable when the
scientific gate does not require it.

Complete when `context.md`, `run_state.json`, and the input manifest describe
every assay, lane, limitation, hash, and blocker.

### 3. Audit every canonical parent

Export one label-free metadata table and run:

```bash
python scripts/audit_parents.py \
  --config config.yaml --cells cells.tsv \
  --output adaptive_parent_audit.tsv
```

Audit guide disagreement/entropy, meaningful minorities, incompatible programs,
available doublet evidence, donor/capture dominance when present, and configured
technical metrics. Missing donor, capture, RNA complexity, mitochondrial, or
other study-specific fields are declared limitations rather than fabricated
values. Flagging opens investigation; it never proves a split.

Complete when every parent has one audit row.

### 4. Discover bounded candidate children

Screen stored higher-resolution membership first. Do not assume Leiden
partitions are nested. Collapse equivalent children into persistent families by
parent-restricted membership and require parent purity plus persistence.

If stored partitions fail, permit one bounded parent-local pass using the frozen
preprocessing, graph inputs, batch covariates, and assay representations. The
local method may be WNN, a single-assay neighbor graph, or another declared
study graph; reproduce the canonical construction rather than imposing WNN.

```bash
python scripts/discover_candidates.py \
  --config config.yaml --parents parent_membership.tsv \
  --stored-memberships stored_membership_long.tsv \
  --local-candidates parent_local_candidates.tsv \
  --audit adaptive_parent_audit.tsv --output-dir discovery
python scripts/validate_candidate_registry.py \
  --config config.yaml --audit adaptive_parent_audit.tsv \
  --registry discovery/candidate_registry.tsv \
  --membership discovery/candidate_membership.tsv
```

Select the lowest stable biologically resolving partition. Reject stable
one-cluster results, continuous gradients, and technical partitions. Keep local
discoveries label-free until full evidence testing.

Complete when every flagged parent is routed to testing or documented rejection.

### 5. Build binary comparator lanes

For each credible candidate construct two independent comparisons:

1. candidate versus canonical-parent remainder;
2. candidate versus closest biological or sibling comparator.

Expand both comparisons over every configured evidence lane:

```bash
python scripts/build_candidate_contrasts.py \
  --config config.yaml --cells cells.tsv \
  --registry discovery/candidate_registry.tsv \
  --membership discovery/candidate_membership.tsv \
  --output-dir contrasts
```

The lane count is `candidates x 2 x configured_lanes`; it is never hard-coded.
A multinomial one-versus-rest coefficient satisfies neither binary comparison.

### 6. Run assay-aware evidence through Snakemake

Adapt [workflow-template](assets/workflow-template/) only for paths, registered
backends, and resources. For each `rapids_singlecell` lane run native Wilcoxon
and balanced binary logistic regression, then validate receipts:

```bash
python scripts/rank_markers_gpu.py --help
python scripts/validate_gpu_receipts.py \
  --config config.yaml --lane-registry contrasts/comparator_manifest.tsv \
  --receipts-root gpu_markers --output gpu_receipts.validated.json
```

Non-RAPIDS lanes use a declared project/backend plugin and their own immutable
receipt contract. Never relabel a CPU result as RAPIDS evidence.

For every applicable lane report means, effect sizes, adjusted P values,
positive-program coverage, expected-negative violations, and detection fractions
when valid. Add source corroboration only where compatible: PanglaoDB ORA for
gene or gene-activity spaces, peak/gene or motif methods for ATAC when declared,
and panel-aware protein evidence for antibody assays. Also compute cell-level
program co-occurrence, stability, available QC/doublet evidence, and technical
composition summaries.

Summarize with `scripts/summarize_multimodal_evidence.py`, then consolidate
configured candidate gates with `scripts/evaluate_candidate_evidence.py`.

Complete when every required lane and source gate has a validated artifact, or
the candidate is failed closed with a specific missing capability.

### 7. Make one outcome per parent

Assign exactly one reviewed outcome:

- `ACCEPT SPLIT`
- `RETAIN PARENT`
- `DOUBLET/MIXED`
- `TECHNICAL/UNRESOLVED`

Accept a child only when all configured biological, assay, technical, diversity,
doublet, stability, and comparator gates pass. Do not require an assay the study
does not contain. Do require every lane declared necessary for the claims being
made. Failed candidates fold into the correct remainder even when a sibling
candidate passes.

### 8. Critic review and interactive annotation

Run independent critic passes over raw tables. Test biological coherence,
positive and negative evidence, doublet explanations, technical confounding,
assay limitations, and semantic hierarchy. Record each finding and resolution.

Present one scope, parent, or candidate at a time. Show rejected alternatives,
uncertainty, guide divergence, and type/state separately. Use the Knowledgebase
to validate labels and CL IDs. Computational completion is not approval.

### 9. Build the immutable packet and validate mapping

```bash
python scripts/build_review_packet.py --config config.yaml --run-root RUN_ROOT
python scripts/validate_final_mapping.py \
  --config config.yaml --canonical parent_membership.tsv \
  --decisions parent_decision_summary.tsv \
  --registry discovery/candidate_registry.tsv \
  --candidate-membership discovery/candidate_membership.tsv \
  --mapping parent_child_mapping.tsv --build-mapping \
  --output mapping.validation.json
```

Build leaf IDs from accepted candidate IDs. Assign accepted children first and
fold all remaining cells into a deterministic remainder. Preserve unsplit parent
membership byte-for-byte.

Complete only when each cell maps once, parent counts conserve, leaf totals are
correct, input/packet hashes are invariant, critics are reconciled, ontology
semantics are valid, and human approval is recorded. Keep write-back separate;
prohibit it in `evidence-only` mode.

## Bundled utilities

Treat scripts as deterministic workflow components, not substitutes for
biological review. Put marker panels, guide columns, assay names, feature spaces,
thresholds, backend plugins, and HPC resources in `config.yaml`.

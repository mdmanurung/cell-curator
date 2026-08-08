# Parity with QuadBio cell-type-annotation

This is a development assessment, not a claim of feature parity. It compares
this repository with QuadBio
[`cell-type-annotation`](https://github.com/quadbio/cell-type-annotation) at
commit `b4381198b85b7fd422f6640c1d7a04c2978ace1d` (`cta` 0.1.0), the version
pinned in [quadbio-integration.md](quadbio-integration.md).

## Verdict

QuadBio is the stronger executable product for transcriptomic annotation today.
This skill has the stronger design for immutable canonical partitions, adaptive
nested subclustering, multimodal evidence, and strict HPC/GPU provenance. It is
currently an orchestration and contract layer over QuadBio, not a feature-complete
replacement.

Use QuadBio as the annotation engine. Develop this package as its assay-aware,
adaptive, production-workflow extension.

## Current comparison

| Capability | Stronger implementation | Assessment |
|---|---|---|
| RNA and spatial-transcriptomic annotation | QuadBio | Installable `cta` CLI, canonicalization, QC, marker computation, cross-check routing, reports, and guarded write-back are implemented. |
| CITE-seq and multimodal decision design | This skill | Evidence lanes, representation-aware detection, nested children, binary comparators, and immutable mappings are explicit contracts. |
| ATAC-only or protein-only annotation | Neither | This skill can declare these assays, but it does not yet ship differential-accessibility, motif, peak-to-gene, or protein-specific annotation backends. |
| Frozen canonical partitions | This skill | Parents remain fixed and accepted children receive hierarchical leaf identities. |
| Adaptive subclustering | This skill in design; QuadBio in executable availability | This skill specifies persistence and parent-local stability; QuadBio ships a simpler bounded Leiden command. The full adaptive engine is not yet implemented here. |
| Automated annotation cross-checks | QuadBio | CellTypist, CellMapper, scANVI/scArches, Azimuth, SCimilarity, and environment-aware routing are concrete commands. |
| HPC, immutable artifacts, and GPU provenance | This skill | Snakemake authority, versioned packets, RAPIDS receipts, and GPU fail-closed behavior are stronger contracts. |
| Packaging and reproducibility | QuadBio | `pyproject.toml`, `pixi.toml`, committed lockfile, unified CLI, and broader tests are present. |
| Interactive evidence review | QuadBio | Evidence-card scaffolding, HTML reports, completeness checks, hierarchy validation, and dry-run write-back are implemented. |
| Parent/child accounting | This skill | Exactly one decision per parent, one leaf per cell, count conservation, and lineage/global leaf totals are explicit. |

At the pinned versions, QuadBio contains about 8,022 Python source lines under
`src/cta`, 12 test modules, and 64 test functions. This skill contains about
2,853 Python lines under `scripts`, five test modules, and 22 test functions.
These counts measure implementation surface, not scientific quality.

## Capabilities this skill adds

- Freeze and hash canonical cluster membership before refinement.
- Represent accepted children with stable hierarchical IDs and conserve every
  parent cell in an immutable one-row-per-cell mapping.
- Distinguish `ACCEPT SPLIT`, `RETAIN PARENT`, `DOUBLET/MIXED`, and
  `TECHNICAL/UNRESOLVED` as mutually exclusive reviewed parent outcomes.
- Screen stored partitions for parent purity and membership persistence before
  permitting one bounded parent-local clustering pass.
- Require independent candidate-versus-parent-remainder and
  candidate-versus-sibling binary comparisons.
- Register arbitrary assays and representations without assuming every study
  has RNA and ADT.
- Separate signed corrected values used for means/effects from nonnegative
  matrices used for detection fractions.
- Treat donor, capture, doublet, QC, and technical structure as acceptance
  gates rather than report-only metadata.
- Fail closed when a configured `rapids_singlecell` lane lacks a validated GPU;
  no CPU marker result may impersonate that lane.
- Keep computational completion, biological review, ontology validation, human
  approval, and write-back as separate states.

## Material gaps and correctness risks

### Evidence gates are not yet provenance-derived

`evaluate_candidate_evidence.py` accepts boolean fields such as positive-marker,
negative-marker, source-corroboration, and lane-pass gates from an input table.
It does not derive those decisions from hashed marker, enrichment, QC, and
stability artifacts. A malformed or optimistic upstream table can therefore
produce a passing decision without a cryptographic and semantic path back to raw
evidence.

Replace free gate booleans with deterministic evaluators whose inputs are
validated receipts and whose outputs include source hashes and threshold traces.

### The Snakemake workflow is a skeleton

The bundled workflow does not produce prepared cell tables, parent membership,
stored-resolution membership, parent-local candidates, dynamic comparator jobs,
source corroboration, program co-expression, technical summaries, critics,
parent decisions, final mappings, or leaf totals. It also makes the review packet
the `all_evidence` target even though the written contract says review packets
and finalizers are outside default aggregate targets.

Implement a checkpoint- or manifest-driven DAG and make evidence completion,
review-packet construction, and finalization separate explicit targets.

### Candidate discovery trusts external decisions

Parent-local candidate input currently supplies `stable`,
`biologically_resolving`, `technical_dominated`, and `eligible` booleans. The
bundled code neither builds the local graph nor calculates resampling stability.
The stored-resolution family matcher is greedy, is not mutual-best or one-to-one,
and selects the lowest-resolution family member without proving it is the best
representative. `closest_sibling` is also supplied rather than derived.

Implement graph construction, adjacent-resolution persistence, resampling ARI,
technical association tests, comparator selection, and deterministic candidate
promotion inside the workflow.

### Declared modality support exceeds executable support

RNA-only, CITE-seq, ATAC-only, protein-only, gene-activity, and multiome inputs
can be described by the configuration schema. Only RAPIDS marker lanes have a
bundled ranking implementation, and that runner consumes prepared AnnData rather
than extracting arbitrary MuData modalities. No ATAC accessibility, motif,
peak-to-gene, or protein-specific backend is included.

Report modality support at four levels:

1. `annotation-capable`: complete identity evidence and ontology validation;
2. `corroboration-only`: can support or refute another assay's identity claim;
3. `qc-only`: usable only for technical assessment;
4. `unavailable`: declared but not implemented for the run.

Do not describe a modality as supported merely because it can be registered.

### Parent decisions have a nullable-threshold bug

`decide_parents` converts `adaptive.audit.doublet_fraction_max` directly with
`float(...)`. The template permits that threshold to be null, causing
`float(None)` in a default-config decision run. Tests currently avoid the bug by
supplying a threshold.

Define null as a disabled gate consistently and add a default-template regression
test.

### Review-packet completion is weaker than the contract

Mapping validation is optional to the packet builder even though one-row-per-cell
mapping and count conservation are completion requirements. The packet builder
also expects a GPU receipt aggregate even when a study has no RAPIDS lanes, while
the workflow does not define the zero-lane artifact path.

Make mapping validation mandatory and define explicit zero-applicable-lane
receipts for every optional backend family.

### GPU receipts need stronger runtime checks

The RAPIDS runner requires one visible validated device, but records the
configured distribution string without confirming the installed distribution.
It does not validate logistic-regression convergence or implement the documented
retry behavior. Prepared matrix identity is not fully tied to comparator and
feature manifests.

Record observed package distribution, solver convergence, retry history, input
matrix hash, feature order hash, group-membership hash, and comparator-manifest
hash in each receipt.

### Packaging and user-facing execution lag QuadBio

This repository has no unified CLI, package metadata, locked environment, method
router, HTML report renderer, canonicalization command, or ontology-aware
write-back command. The scripts are useful workflow components but require
project-specific glue.

Create one installable CLI that wraps the deterministic scripts and delegates
transcriptomic annotation functions to QuadBio rather than duplicating them.

## Generalization target

The portable abstraction is an evidence-capability graph, not an RNA-plus-ADT
checklist:

```text
study input
  -> assay representations
  -> executable evidence capabilities
  -> claims each capability may support or refute
  -> candidate-specific required gates
  -> reviewed type and state decisions
```

Examples:

| Study | Valid route |
|---|---|
| RNA only | QuadBio RNA evidence plus adaptive parent/child, stability, technical, and mapping gates. |
| CITE-seq | RNA identity evidence plus panel-aware protein corroboration; corrected protein detection requires a separate nonnegative source. |
| Protein only | Annotation only when a sufficiently discriminative validated panel and ontology-compatible rules exist; otherwise evidence-only. |
| ATAC only | Use accessibility, motif, gene-activity, and reference-mapping evidence when implemented; peak-level clustering alone cannot justify gene-marker labels. |
| RNA plus ATAC | RNA supplies direct identity evidence; ATAC supplies regulatory corroboration or independently supported subtype evidence. |
| New assay | Add a backend plugin, receipt validator, evidence schema, claim scope, and tests before marking it annotation-capable. |

Missing modalities are not failures. Missing capabilities required for a claim
are blockers. Acceptance gates must be assembled from the study's available,
validated capabilities rather than from a universal modality checklist.

## Development order toward genuine parity

### P0: make decisions trustworthy

1. Derive all gates from hashed raw artifacts.
2. Fix nullable thresholds and require mapping validation.
3. Separate evidence, packet, and finalization workflow targets.
4. Define zero-lane behavior and fail closed on every required capability.

### P1: make the workflow executable

1. Implement input preparation for AnnData and MuData.
2. Implement stored-child persistence and parent-local graph/resampling stability.
3. Expand comparator lanes dynamically through Snakemake.
4. Implement program co-expression, source corroboration, technical tests,
   critics, mappings, and leaf totals as producing rules.
5. Add a unified CLI, package metadata, and locked environments.

### P2: become broadly multimodal

1. Ship reference RNA, protein, ATAC, gene-activity, and multiome adapters.
2. Add backend-specific receipt validators and capability routing.
3. Reuse QuadBio cross-checks, reports, hierarchy checks, and write-back instead
   of rebuilding them.
4. Run end-to-end tests on RNA-only, CITE-seq, ATAC-only, multiome, real GPU,
   evidence-only, doublet/mixed, and successful adaptive-split fixtures.

## Readiness bar for the term "supercharged QuadBio"

Use that description only after:

- QuadBio's executable annotation paths remain available through this skill;
- every advertised modality has at least one complete reference backend;
- all acceptance gates are reproducible from immutable raw evidence;
- the Snakemake workflow resumes end to end without project-specific glue;
- the package has a locked environment and unified CLI;
- real-data and synthetic forward tests cover accepted splits, retained
  gradients, technical partitions, and mixed/doublet parents;
- downstream invariance checks demonstrate that canonical inputs and unsplit
  memberships remain unchanged.

Until then, describe this repository as an adaptive multimodal orchestration
layer over QuadBio.

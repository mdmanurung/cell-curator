# Neutral capability-parity matrix

Use the two public QuadBio surfaces below only as a capability checklist. The
mappings identify independent `cell-curator` implementations; they do not import,
invoke, vendor, or require QuadBio, its CLI, package, skill, MCP configuration, or
artifacts. They make no comparative accuracy or performance claim.

**This is reviewer orientation, not a routing input.** An agent reading this file
is deciding how to annotate data, not which project to prefer. Nothing here should
change what a run does, and a gap recorded below is a statement about this
package's surface, not an instruction to route around it.

Checklist sources, both documentation only and never used at runtime:

- QuadBio [`cell-type-annotation`](https://github.com/quadbio/cell-type-annotation/blob/main/skills/cell-type-annotation/SKILL.md) — an agent skill plus a `cta` CLI, last commit `b438119` (2026-08-06).
- QuadBio [`cell-annotator`](https://github.com/quadbio/cell-annotator) — a PyPI package, `cell-annotator` 0.3.0, head `0b17fc3` (2026-07-07).

### The two QuadBio projects are unrelated to each other

They share a lab and a problem, and nothing else. Neither repository contains a
cross-reference to the other, they share no package or namespace, and
`cell-annotator` was created roughly eighteen months *before* the skill, which
rules out the skill being a port of it. They also make opposite architectural
bets: `cell-annotator` uses a language model as the classifier, while
`cell-type-annotation` uses a language model as the reasoning agent over evidence
computed classically. The matrix below is therefore drawn from the skill, and
`cell-annotator` is compared separately further down.

## The matrix, against `cell-type-annotation`

| QuadBio capability | `cell-curator` equivalent | Implementation file | CLI / API entry point | Test coverage |
| --- | --- | --- | --- | --- |
| Independently installable annotation skill | Build a typed wheel containing the CLI, skill, schema, workflow, references, and local-first Python package without QuadBio artifacts | [`pyproject.toml`](../../../pyproject.toml), [`__init__.py`](../../../src/cell_curator/__init__.py) | `cell-curator --version`, `cell-curator --help`; package imports | [`test_clean_wheel.py`](../../../tests/test_clean_wheel.py), [`test_package_core.py`](../../../tests/test_package_core.py) |
| Three-stage guided workflow | Persist scene context, construct an executable plan, then execute parent-scoped levels | [`guided.py`](../../../src/cell_curator/guided.py), [`pipeline.py`](../../../src/cell_curator/pipeline.py) | `run scene`, `run plan`, `run execute`; `prepare_scene()`, `construct_plan()`, `run_annotation()` | [`test_end_to_end.py`](../../../tests/test_end_to_end.py), [`test_environment_routing.py`](../../../tests/test_environment_routing.py) |
| Interactive context confirmation | Confirm or correct a complete biological and write-back context with a named reviewer | [`guided.py`](../../../src/cell_curator/guided.py) | `run confirm-context`; `confirm_context()` | [`test_end_to_end.py`](../../../tests/test_end_to_end.py) |
| Pre-authorized non-interactive operation | Validate complete context, reviewer, assumptions, and pre-authorization; record the same gate trail | [`config.py`](../../../src/cell_curator/config.py), [`guided.py`](../../../src/cell_curator/guided.py) | `run scene`, `run plan`, `run execute`; `GuidanceConfig`, `approve_assumptions()` | [`test_environment_routing.py`](../../../tests/test_environment_routing.py), [`test_end_to_end.py`](../../../tests/test_end_to_end.py) |
| Living assumptions and renewed gates | Add, list, hash, scope, approve, and re-check material assumptions before labels | [`guided.py`](../../../src/cell_curator/guided.py), [`models.py`](../../../src/cell_curator/models.py) | `assumptions list`, `assumptions add`, `assumptions approve`; `add_assumptions()`, `approve_assumptions()`, `require_label_gate()` | [`test_end_to_end.py`](../../../tests/test_end_to_end.py) |
| CPU, RAM, Apple Silicon, CUDA, scheduler, and resource-limit detection | Read-only structured environment detection with login-node ambiguity and real-GPU validation | [`environment.py`](../../../src/cell_curator/environment.py) | `environment detect`; `detect_environment()`, `validate_real_gpu()` | [`test_environment_routing.py`](../../../tests/test_environment_routing.py) |
| Deterministic method routing per hierarchy level | Select one primary cross-check from biological context, technology, resources, memory, models, and references; retain every skip reason | [`routing.py`](../../../src/cell_curator/routing.py) | `route select`; `route_crosscheck()`, `route_hierarchy()` | [`test_environment_routing.py`](../../../tests/test_environment_routing.py) |
| AnnData and MuData inspection | Validate cells, clusters, modalities, layers, raw, embeddings, visualization, spatial coordinates, representations, and colors | [`data.py`](../../../src/cell_curator/data.py) | `data inspect`; `inspect()` | [`test_data_parity.py`](../../../tests/test_data_parity.py), [`test_input_and_audit.py`](../../../tests/test_input_and_audit.py) |
| Canonical resumable working copy | Freeze a source-preserving canonical copy and manifests with drift checks | [`data.py`](../../../src/cell_curator/data.py), [`provenance.py`](../../../src/cell_curator/provenance.py) | `data canonicalize`; `canonicalize()` | [`test_data_parity.py`](../../../tests/test_data_parity.py), [`test_input_and_audit.py`](../../../tests/test_input_and_audit.py) |
| Read-only QC and cell-cycle metrics | Compute depth, detected features, mitochondrial, ribosomal, cell-cycle, phase, and feasible assay metrics without source mutation | [`data.py`](../../../src/cell_curator/data.py) | `data qc`; `compute_qc()` | [`test_data_parity.py`](../../../tests/test_data_parity.py) |
| Representation and normalization safeguards | Fingerprint counts, log-normalized, corrected/signed, raw, layered, and detection-safe matrices; hard-stop ambiguity | [`data.py`](../../../src/cell_curator/data.py) | `data inspect`, `data canonicalize`, `data qc`; `normalization_fingerprint()`, `matrix_for()` | [`test_data_parity.py`](../../../tests/test_data_parity.py) |
| Gene-ID canonicalization | Detect symbols and Ensembl IDs, remove Ensembl versions, apply supplied mappings, and reject ambiguous identifiers | [`data.py`](../../../src/cell_curator/data.py), [`crosschecks.py`](../../../src/cell_curator/crosschecks.py) | `data canonicalize`; `gene_symbols()`, `load_gene_mapping()`, `harmonize_reference_genes()` | [`test_data_parity.py`](../../../tests/test_data_parity.py), [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py) |
| Targeted-panel and covariate awareness | Summarize panel capability and missingness, discover categorical covariates, and preserve categorical colors | [`data.py`](../../../src/cell_curator/data.py), [`review.py`](../../../src/cell_curator/review.py) | `data inspect`, `data qc`; `targeted_panel_summary()`, `categorical_covariate_inventory()` | [`test_data_parity.py`](../../../tests/test_data_parity.py) |
| Marker-program assembly | Merge local signed lists, model-export tables, cached tables, and labeled references with strict parent filtering | [`markers.py`](../../../src/cell_curator/markers.py), [`knowledge.py`](../../../src/cell_curator/knowledge.py) | `markers assemble`; `assemble_marker_programs()`, `collect_programs()` | [`test_marker_parity.py`](../../../tests/test_marker_parity.py), [`test_package_core.py`](../../../tests/test_package_core.py) |
| Bottom-up cluster-marker discovery | Compute one-cluster-versus-remainder markers and specificity filters without conflating condition-level DE | [`evidence.py`](../../../src/cell_curator/evidence.py), [`markers.py`](../../../src/cell_curator/markers.py) | `markers compute`; `compute_bottom_up_markers()`, `compute_de_marker_overlap()` | [`test_marker_parity.py`](../../../tests/test_marker_parity.py), [`test_multimodal_engine.py`](../../../tests/test_multimodal_engine.py) |
| Top-down positive and negative panel profiling | Emit raw mean/detection matrices, panel gates, missing features, uncallable programs, and expected-negative violations | [`markers.py`](../../../src/cell_curator/markers.py) | `markers profile`; `profile_marker_panel()`, `negative_marker_violations()` | [`test_marker_parity.py`](../../../tests/test_marker_parity.py) |
| AUCell, UCell, and DE-overlap cross-checks | Compute deterministic cell-rank enrichment and cluster-marker overlap for signed programs | [`markers.py`](../../../src/cell_curator/markers.py) | `markers score`; `score_aucell()`, `score_ucell()`, `compute_de_marker_overlap()` | [`test_marker_parity.py`](../../../tests/test_marker_parity.py) |
| Consolidated marker evidence and plots | Preserve route components and disagreement in durable tables plus a data-backed self-contained SVG dotplot | [`markers.py`](../../../src/cell_curator/markers.py) | `markers score`; `consolidate_marker_crosschecks()`, `write_marker_artifacts()` | [`test_marker_parity.py`](../../../tests/test_marker_parity.py) |
| Existing prediction columns | Convert existing labels into cell and cluster hypotheses with uncertainty rather than verdicts | [`crosschecks.py`](../../../src/cell_curator/crosschecks.py) | `crosscheck obs-pred`; `existing_prediction_hypothesis()` | [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py) |
| CellTypist discovery and prediction | Discover installed/cached models and run optional prediction with actionable unavailability | [`crosschecks.py`](../../../src/cell_curator/crosschecks.py) | `crosscheck celltypist-discover`, `crosscheck celltypist-predict`; `discover_celltypist_models()`, `celltypist_hypothesis()` | [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py) |
| Labeled-reference kNN | Harmonize genes, map by CPU kNN with uncertainty, and allow GPU kNN only after real-GPU validation | [`crosschecks.py`](../../../src/cell_curator/crosschecks.py), [`environment.py`](../../../src/cell_curator/environment.py) | `crosscheck reference-knn`; `reference_knn_hypothesis()` | [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py), [`test_environment_routing.py`](../../../tests/test_environment_routing.py) |
| scVI/scANVI reference training | Train an optional reference model under an explicit validated GPU contract | [`crosschecks.py`](../../../src/cell_curator/crosschecks.py) | `reference train-scanvi`; `train_scanvi_reference()` | [`test_optional_provider_contracts.py`](../../../tests/test_optional_provider_contracts.py), [`test_environment_routing.py`](../../../tests/test_environment_routing.py) |
| scArches-style query mapping | Map queries through an optional compatible model and retain uncertainty | [`crosschecks.py`](../../../src/cell_curator/crosschecks.py) | `reference map-scarches`; `map_scarches_query()` | [`test_optional_provider_contracts.py`](../../../tests/test_optional_provider_contracts.py), [`test_environment_routing.py`](../../../tests/test_environment_routing.py) |
| Curated human reference mapping | Run an optional Azimuth-compatible or equivalent curated human-reference hypothesis route | [`crosschecks.py`](../../../src/cell_curator/crosschecks.py) | `crosscheck azimuth`; `curated_human_reference_hypothesis()` | [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py) |
| Foundation-model cross-check | Run an optional SCimilarity-compatible route or return exact dependency/model remediation | [`crosschecks.py`](../../../src/cell_curator/crosschecks.py) | `crosscheck scimilarity`; `scimilarity_hypothesis()` | [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py) |
| Bounded CELLxGENE Census acquisition | Acquire a bounded optional subset with cache, version, timeout, checksum, and explicit failure state | [`crosschecks.py`](../../../src/cell_curator/crosschecks.py), [`knowledge.py`](../../../src/cell_curator/knowledge.py) | `reference acquire-census`; `fetch_census_reference()` | [`test_optional_provider_contracts.py`](../../../tests/test_optional_provider_contracts.py), [`test_environment_routing.py`](../../../tests/test_environment_routing.py) |
| Pinned Azimuth and SCimilarity bundle provisioning | Acquire one pinned bundle over https under a required digest, a streamed size bound, and a timeout, with an explicit failure receipt that never deletes an operator's cache | [`crosschecks.py`](../../../src/cell_curator/crosschecks.py) | `reference acquire-azimuth`, `reference acquire-scimilarity`; `fetch_pinned_reference_bundle()`, `fetch_azimuth_reference()`, `fetch_scimilarity_model()` | [`test_optional_provider_contracts.py`](../../../tests/test_optional_provider_contracts.py) |
| Multimodal evidence preservation | Score RNA, protein, gene activity, ATAC, and spatial routes according to declared capabilities and retain disagreement | [`evidence.py`](../../../src/cell_curator/evidence.py), [`capabilities.py`](../../../src/cell_curator/capabilities.py) | `run execute`; `run_lane()`, `aggregate_multimodal_evidence()` | [`test_multimodal_engine.py`](../../../tests/test_multimodal_engine.py) |
| Parent-local refinement and stored partitions | Inspect stored resolutions first, run one bounded local graph pass, and emit child markers, QC, composition, and spatial summaries | [`refinement.py`](../../../src/cell_curator/refinement.py), [`hierarchy.py`](../../../src/cell_curator/hierarchy.py) | `cluster subcluster`, `cluster profile`; `discover_refinement_candidates()`, `profile_refined_clusters()` | [`test_ontology_and_refinement.py`](../../../tests/test_ontology_and_refinement.py) |
| Ambiguity and parallel state profiling | Detect ambiguous, retained-parent, mixed/doublet, technical, cycling, stressed, gradient, and low-quality states independently of identity | [`hierarchy.py`](../../../src/cell_curator/hierarchy.py), [`audit_parents.py`](../../../src/cell_curator/audit_parents.py) | `cluster profile`; `profile_ambiguity()`, `classify_parallel_state()` | [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py), [`test_input_and_audit.py`](../../../tests/test_input_and_audit.py) |
| Advisory confidence | Combine evidence margin, specificity, agreement, and coverage without converting automation into a verdict | [`hierarchy.py`](../../../src/cell_curator/hierarchy.py) | `cluster confidence`; `advisory_confidence()` | [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py) |
| Recursive parent-scoped assignment | Execute arbitrary-depth exact parent scopes, preserve rejected alternatives and rationale, retain missing-child parents, hash checkpoints, and resume only valid scopes | [`pipeline.py`](../../../src/cell_curator/pipeline.py), [`hierarchy.py`](../../../src/cell_curator/hierarchy.py), [`guided.py`](../../../src/cell_curator/guided.py) | `run execute --score-manifest`; `execute_recursive_hierarchy()`, `execute_recursive_hierarchy_manifest()`, `checkpoint_level()` | [`test_recursive_hierarchy_execution.py`](../../../tests/test_recursive_hierarchy_execution.py), [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py) |
| Strict hierarchy and one-cell-one-leaf validation | Reject child reuse across parents, invalid nesting, overlapping children, duplicate cells, and blank leaves; permit `Unknown` | [`hierarchy.py`](../../../src/cell_curator/hierarchy.py), [`decisions.py`](../../../src/cell_curator/decisions.py) | `hierarchy validate`; `validate_hierarchy_nesting()`, `build_mapping()` | [`test_crosschecks_hierarchy_focused.py`](../../../tests/test_crosschecks_hierarchy_focused.py) |
| Durable labels and evidence cards | Use labels TSVs as source of truth, preserve manual reasoning, and scaffold cards from raw evidence values | [`hierarchy.py`](../../../src/cell_curator/hierarchy.py), [`review.py`](../../../src/cell_curator/review.py) | `evidence-card scaffold`; `scaffold_labels_table()`, `build_evidence_cards()` | [`test_end_to_end.py`](../../../tests/test_end_to_end.py) |
| Platform-neutral critic workflow | Build structured packets, import validated findings, and reconcile `RESOLVED` or `DISMISSED` dispositions with reviewer and rationale | [`review.py`](../../../src/cell_curator/review.py) | `critic packet`, `critic import-findings`, `critic reconcile`, `critic validate`; `build_critic_packet()`, `import_critic_findings()`, `import_reconciliation()` | [`test_package_core.py`](../../../tests/test_package_core.py) |
| Per-level and full-run HTML reporting | Build self-contained embeddings, dotplots, hierarchy, uncertainty, alternatives, QC, spatial, composition, color, and provenance panels | [`review.py`](../../../src/cell_curator/review.py) | `report build`; `build_html_report()` | [`test_end_to_end.py`](../../../tests/test_end_to_end.py), [`test_multimodal_engine.py`](../../../tests/test_multimodal_engine.py) |
| Strict read-only report checking | Validate required labels, reasoning, panels, assumptions, and existing report without rebuilding approved artifacts | [`review.py`](../../../src/cell_curator/review.py) | `report check`; `check_review_report()`, `check_report_completeness()` | [`test_end_to_end.py`](../../../tests/test_end_to_end.py) |
| Progress and next-action reporting | Report current stage, hierarchy level, parent, blockers, completed checkpoints, artifacts, and next action | [`guided.py`](../../../src/cell_curator/guided.py) | `progress status`, `run status`; `status()` | [`test_end_to_end.py`](../../../tests/test_end_to_end.py) |
| Local Cell Ontology validation | Resolve names, synonyms, accessions, obsolete replacements, ancestors, and compatibility; reject invalid supplied mappings | [`ontology.py`](../../../src/cell_curator/ontology.py) | `ontology validate`, `hierarchy validate`; `CellOntology`, `validate_program_ontology()` | [`test_ontology_and_refinement.py`](../../../tests/test_ontology_and_refinement.py) |
| Optional remote ontology enhancement | Provision a versioned, timeout-bounded, checksum-pinned ontology cache and record explicit unavailability while retaining CL-free operation | [`ontology.py`](../../../src/cell_curator/ontology.py), [`knowledge.py`](../../../src/cell_curator/knowledge.py) | `ontology provision`; `ontology_from_config()`, `RemoteTableProvider.provision()` | [`test_ontology_and_refinement.py`](../../../tests/test_ontology_and_refinement.py) |
| Optional live ontology corroboration | Query a pinned https endpoint for CL terms, store each answer as canonical hashed JSON with its endpoint, query, digest, and retrieval time, and report agreement, obsolescence, name disagreement, or unreachability without grounding any label | [`ontology.py`](../../../src/cell_curator/ontology.py) | `ontology corroborate`; `corroborate_terms_live()`, `live_ontology_provider()` | [`test_ontology_and_refinement.py`](../../../tests/test_ontology_and_refinement.py) |
| CL-free annotation | Keep labels valid without ontology evidence while preventing unvalidated CL identifiers from entering write-back | [`ontology.py`](../../../src/cell_curator/ontology.py), [`pipeline.py`](../../../src/cell_curator/pipeline.py) | `run plan`, `run validate`; `ontology_from_config()` | [`test_ontology_and_refinement.py`](../../../tests/test_ontology_and_refinement.py), [`test_end_to_end.py`](../../../tests/test_end_to_end.py) |
| Complete guarded write-back | Validate hierarchy and ontology, preview a complete diff, preserve category order and colors, require approval, and write atomically to a new file by default | [`writeback.py`](../../../src/cell_curator/writeback.py) | `write-back preview`, `write-back apply`; `validate_writeback()`, `write_back()` | [`test_end_to_end.py`](../../../tests/test_end_to_end.py) |
| Immutable provenance and resume | Hash inputs, configuration, mappings, receipts, and artifacts; resume valid checkpoints through the authoritative DAG | [`provenance.py`](../../../src/cell_curator/provenance.py), [`Snakefile`](../../../workflow/Snakefile) | `run status`, `run validate`; Snakemake targets `all_evidence`, `review_packet`, `finalization`, `writeback` | [`test_end_to_end.py`](../../../tests/test_end_to_end.py) |
| Discoverable unified CLI and shared Python logic | Expose command help while keeping commands thin over the same typed package APIs used by tests and workflows | [`cli.py`](../../../src/cell_curator/cli.py), [`api.py`](../../../src/cell_curator/api.py) | `cell-curator --help`; all command groups listed above | [`test_cli_and_manifests.py`](../../../tests/test_cli_and_manifests.py), [`test_package_core.py`](../../../tests/test_package_core.py) |

## Where `cell-annotator` differs

The skill is the checklist; the package is a different project with a different
bet, and the comparison is worth recording because it is not one-sided in either
direction.

| Capability | `cell-annotator` | `cell-curator` |
| --- | --- | --- |
| Where the call is made | A language model classifies each cluster from its marker list | Scored locally; no model call, no provider key |
| Cell Ontology | None. Labels are free text, de-duplicated only against each other | Local CL snapshot, semantic validation, optional live tier |
| Hierarchy | Flat. A single `cluster_key`, no recursion or sub-clustering | Arbitrary depth, every child scoped to its parent |
| State axis | `cell_state` is a field the agent may fill, and is not written to `obs` | Orthogonal to identity and refused as an identity candidate |
| Reference mapping | None. No CellTypist, Azimuth, scANVI, or kNN route | Routed, with bounded acquisition and skip reasons |
| Uncertainty | A self-reported `Low`/`Medium`/`High` string, also not written to `obs` | Margin-anchored confidence; `Unknown` is machine-gated |
| Cross-sample vocabulary consolidation | Yes, local labels are mapped onto one global vocabulary | Absent. See the gap recorded below |
| Category colour assignment | Model-proposed hex colours validated with a ΔE threshold and re-queried on failure | Preserves existing category order and colours |
| Network | Mandatory. Every annotation is an outbound call to a commercial provider | Not required. Runs with no network |
| Provenance | None. Re-running re-issues every call | Hash-bound, immutable artifacts, resumable DAG |
| Orchestration | None | Snakemake DAG with a scheduler profile |
| Published evaluation | None | None |

The two shared refusals are worth stating plainly. Both projects disclaim any
accuracy comparison, and both require a network service that `cell-curator` does
not: `cell-annotator` requires a commercial LLM provider, `cell-type-annotation`
requires the BioContextAI knowledge MCP.

## Reverse direction: gaps and deliberate divergences

The matrix above records what `cell-curator` implements independently. It says
nothing about the reverse direction, which anyone comparing the projects needs.
Reviewed against both public sources on 2026-09-28.

### Open gaps

| QuadBio capability | `cell-curator` status | What closing it requires |
| --- | --- | --- |
| Live query-time grounding against the Human Protein Atlas and Europe PMC | Absent. PanglaoDB is reachable as a checksum-pinned `remote_table`, and the Cell Ontology has a live tier (below), but expression-atlas and literature queries have no provider. | Two further `live_*` provider kinds on the contract the live ontology tier already establishes. Lower value than the ontology tier: neither is direct marker evidence for a label. |
| Curated literature grounding through a Zotero MCP | Not a gap of equal standing, and previously recorded here as one. The skill treats Zotero as strongly preferred rather than required and documents a fallback to its own Europe PMC tools, so a run continues without it. | A local-first provider reading a BibTeX or CSL-JSON export would be a strict convenience over reading papers directly, which the labels table already carries. A documented non-goal until there is a reason to prefer it. |

### Gaps this package has that the review exposed

These are absent from `cell-curator` and are recorded here because a checklist
that only lists another project's gaps is not a neutral document.

| Gap | Status | What closing it requires |
| --- | --- | --- |
| Cross-sample label-vocabulary consolidation | Absent. One clustering over one vocabulary is assumed, so two independently clustered samples can carry `Monocyte` at different granularity with nothing reconciling them. `cell-annotator` maps per-sample labels onto one global vocabulary. | A consolidation pass over the per-parent vocabularies, writing a mapping table and refusing an unresolvable label rather than picking a winner. Note the same machinery already exists for gene identifiers in `harmonize_reference_genes()`; this is the label-axis analogue. |
| A neutral evaluation harness that has never been pointed at real data | The harness exists and is independent of the annotation engine, but five of six dataset manifests are `pending-external-resource` and the only result scores an eight-row fixture. | One real dataset with independent ground truth, run through the shipped CLI, with the label-harmonization and threshold grid frozen before the first result is published. The frozen decisions are recorded in [`benchmarks/README.md`](../../../benchmarks/README.md), because the first published result is immutable and a protocol fixed afterwards is not a protocol. |
| No stated empirical basis for the negative-marker violation threshold | The threshold that decides whether an expected-absent marker counts as violated is asserted. The skill documents a chosen value and the reasoning behind it, and reports the rates behind that choice in a pull request body rather than in the repository, so the asymmetry is one of documented reasoning rather than of published evidence. | Record the rate at which real positive markers exceed the threshold across several reference datasets, and either justify the value or calibrate it. Until then it is a tuned constant with no recorded provenance. |
| The matrix above can only be checked for existence, not for truth | The skill generates its method-selection document from the router and fails CI on drift. `tests/test_parity_matrix_claims.py` closes the mechanical half of the difference by resolving every implementation file, entry point, and symbol this document names. | A reader still has to judge whether the named implementation actually provides the capability, and nothing here measures that. Existence is checkable; equivalence is not. |

### Closed already, listed because a surface comparison suggests otherwise

| QuadBio capability | `cell-curator` equivalent |
| --- | --- |
| Semantic CL identifier validation before write-back: the term exists, is not obsolete, and its name matches the label | `CellOntology.validate_term()` resolves the accession, rejects `is_obsolete` while reporting `replaced_by`, and checks the label against the term name and synonyms; `validate_durable_label_ontology()` runs it at write-back. |
| Live Cell Ontology lookup, so a stale local copy cannot answer silently | `ontology corroborate` queries a configured `live_ontology` endpoint, stores each answer under `knowledge/ontology_live/` as canonical JSON, hashes it, and records the endpoint, query, digest, and retrieval time. It reports agreement, disagreement, obsolescence, and unreachability, and changes no label. Unlike the MCP it parallels, it is optional in both directions: absent, disabled, or unreachable all yield a receipt, and `require_local_grounding_for_labels` stays pinned true so no label can depend on it. |
| Self-bootstrapping isolated environment through pixi | `uv.lock` with `uv sync --frozen --extra workflow`, and `bin/cell-curator` puts the CLI on PATH for an installed plugin without a project install. Parity by other means, not a gap. |
| Sized reference provisioning for Azimuth and SCimilarity | `reference acquire-azimuth` and `reference acquire-scimilarity` mirror the Census route through `fetch_pinned_reference_bundle()`: https-only URL, required SHA-256, `--max-bytes` enforced while the body streams, `--timeout-seconds`, and a receipt recording cache hit and checksum result. A failed digest never deletes the operator's file, and a refused transfer leaves no partial. |

### Deliberate divergence, not a gap

A hard MCP requirement is the one QuadBio design decision this package rejects on
purpose. `cell-curator` must run from local AnnData or MuData, marker files, and
ontology files with no agent platform, no MCP, and no network. Adopting that
requirement would contradict the invariant, so the live provider above is
specified as optional corroboration that degrades to the local cache. It closes
the capability gap without reversing the trade.

The trade is sharper than a checklist alone conveys, and the reason it was made
that way is on the record. In the skill the knowledge MCP is a strict stop rather
than a preference: the run halts instead of degrading to web-only, on the stated
ground that a web-only fallback produces exactly the confident-but-unfounded
annotation the project exists to prevent. That same server also supplies the
Ensembl-to-symbol map, so without it an Ensembl-indexed query cannot be
canonicalized at all. The cost of the requirement is therefore not that a lookup
is slower; it is that the package cannot be used at all in a firewalled or
air-gapped setting, and that a third-party service sits in the path of every
label. `cell-curator` takes the opposite position and pays for it by carrying a
local ontology snapshot that must be provisioned before use.

Against the same checklist, `cell-curator` additionally carries MuData and
multimodal capability routing, a Snakemake DAG with a scheduler profile, a
write-back that a person confirms by typing a token bound to the whole preview,
and a critic loop that runs as a file packet rather than as agent subagents.

The language model is also placed differently. `cell-annotator` uses it as the
classifier. Here it is an optional evidence lane that may only break a near-tie
the measured markers left open: every model must agree, the candidate must pass
its own marker gates, and every cited gene must be a positive marker of the
cluster. The deterministic engine runs without it, and without a network.

## LLM annotators: what leaves the machine, and what checks the answer

Three further projects use a language model to name clusters. Every entry below
was read from the source at the pinned commit, and the path cited is relative to
that repository. None of them checks a model's answer against the measured
evidence.

| | [CyteType](https://github.com/NygenAnalytics/CyteType) `42c4c0e` | [mLLMCelltype](https://github.com/cafferychen777/mLLMCelltype) `26bda56` | [CellTypeAI](https://github.com/rhdaw/CellTypeAI) `5507d62` | `cell-curator` |
| --- | --- | --- | --- | --- |
| What is sent | expression matrix, raw counts, the whole `obs` table with barcodes and any donor columns, embedding coordinates (`cytetype/main.py:301-309`, `core/artifacts.py:515`) | marker gene symbols per cluster, species, tissue (`python/mllmcelltype/prompts.py:29`) | top Wilcoxon gene names per cluster, species, tissue (`src/celltypeai/h5ad_to_json.py:31-32`) | near-tied candidates, top markers, vote shares, pseudobulk means for clusters of at least 20 cells; never per cell |
| Where it goes | Nygen's closed server, `https://cytetype.nygen.io` (`config.py:16`); the user's LLM key goes too (`api/schemas.py:35`, `core/payload.py:90-91`) | a cloud provider; a `base_url` override exists, but an API key is still required (`annotate.py:54-63`) | Ollama at a hard-coded `http://localhost:11434` (`prompt_submit.py:104`) | `knowledge.network`: `offline` (default), `institution` (localhost and allowed hosts), or `online` |
| Local only | no; local Ollama must be exposed to Nygen through ngrok (`docs/ollama.md:1`) | not as a supported mode | yes, only | yes, and the default |
| Abstains | not in the client; confidence is server-side and "model-assessed" (`docs/results.md:24`) | no; without consensus it returns the last majority (`consensus.py:2042`) | no; the mode over repeats (`main.py:106-108`) | yes; Unknown unless the gates pass, or a shared CL ancestor holds the votes |
| Checks the answer | not in the client; the server's code is closed | no; cited genes are stored as text (`utils.py:794`) | no; it checks only the reply's JSON shape (`prompt_submit.py:184-187`) | yes; the gates above |
| Cell Ontology | output fields, mapped on the server (`core/results.py:67-76`) | none | none | the label and its CL id are validated together against a pinned snapshot |

mLLMCelltype's consensus proportion and entropy are computed in code first
(`consensus.py:1143`). In harder cases a model's own figures are used instead,
but only if they pass a plausibility check (`:1194-1199`).

### No published comparison exists for any of the three projects

This is worth stating separately because it is the most consequential fact in
the document, and it is checkable rather than a matter of opinion.

None of these projects publishes a measured accuracy figure for cell-type
annotation. The skill states this three times, and disclaims it most explicitly
where a reader might mistake its method ranking for a result: the ordering "is
not a claim that any tool annotates better than another", and the internal
sanity-checking it describes is "preliminary and unpublished". `cell-annotator`
publishes no benchmark either; the only figure in its tutorials is a confusion
matrix computed against clusters that are the ground-truth labels relabelled, so
each cluster is by construction exactly one true type and the comparison cannot
detect a tool failing to find types in real clusters. Its own README recommends
the output as "an initial, coarse set of annotations that must be further
refined".

`cell-curator` is in the same position and says so: comparative claims require
measured benchmark results, and no such result has been published. No accuracy or
performance claim is made or implied by any row above.

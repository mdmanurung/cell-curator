# Acceptance run, AIDA split A (2026-11-17)

Run `v8-accept-aida-split-a-001` under `cc-v8-work/acceptance-2026-11-17/`. Code: branch `drop-labelling-src` at `a7083bb`, plus the receipt-path fix `5c79dd2` for the review step. Gate A was approved by the user's yes on 2026-11-17 (reviewer Mikhael Manurung), under the standing dev-set confirmation of 2026-11-10. The run stops before critics, `run validate` and write-back apply.

Evidence type: `table fallback: each verdict is exact equality of two label tables`. The only figure, `acceptance-2026-11-17/figs/aida_marker_evidence.png` (from `annotation/levels/L1/marker_evidence.svg`, inspected), is a 23-row marker dotplot that renders. It is too dense at 7856 px wide to read a label from, so it supports "the plot builds", not any label.

## Steps and results

| Step | Result |
| --- | --- |
| `run scene`, `run plan` | pass. The plan's `vocabulary.tsv` has sha256 `25f65f5a…`, the same as run `v8-dev-aida-split-a-003` |
| Gate A | 43 of 44 assumptions identical to run 003. The 44th records the annotation rules, which now include `consensus_directory`, so its text differs |
| `annotation collect-level` | pass, 1 h 20 min on the login node (no ranking-cache hits for this run) |
| annotation-consensus (h5, snap on) | 19 of 23 clusters labelled, 1702 cells flagged, `verify` OK (17 stored files byte-identical to the rerun). Values: the frozen demo values the user chose on 2026-11-13, not tuned here |
| `annotation import-consensus`, `finalize` | pass; hierarchy nesting accepted |
| Imported `labels.tsv` against `consensus.tsv` | 23 of 23 clusters equal in `label` and `cl_id`; 4 `Unknown` (clusters 1, 11, 14, 22) kept |
| `cell_labels.tsv` | 8847 cells; each cell's `L1` equals its cluster's consensus label |
| Review build | 23 evidence cards, critic packet, `report.html` 4.8 MB with 3 inline SVG and 0 `<img>` |

## Findings

- **Fixed:** a receipt written from a relative config path was stale once the working directory changed. The Snakemake review step runs from the run directory, so it refused the evidence receipt. Receipt paths are now absolute (`5c79dd2`, test `test_a_receipt_written_from_a_relative_config_path_is_valid_from_any_directory`). This run's `evidence.manifest.json` still holds the relative key, so its review step was built by direct calls from the acceptance directory, not by Snakemake.
- **Cost:** a second Snakemake attempt deleted the report it was asked to rebuild; it was rebuilt by hand.
- **Container:** rebuilt after the fix (SIF `01703be8d325…`); passes on three CPU types, see `container/README.md`.
- **Not checked:** accuracy. The run scores against split-A truth only if someone runs the skill's `score` step.

## Still the user's

Critic findings, `run validate`, `write-back preview` review and `write-back apply`.

# Acceptance run, PBMC3K (2026-11-17)

Run `v8-accept-pbmc3k-001` under `cc-v8-work/acceptance-2026-11-17/`, config `pbmc3k.yaml`, code at `6ab8d75` (local `main`). Consensus run: `cc-v8-work/pbmc3k-2026-11-13/out_h5` (rule h5, the frozen values the user approved for the showcase on 2026-11-13). Gate A was approved by the user ("approve", reviewer Mikhael Manurung) on 2026-11-17; PBMC3K is not a dev set, so the standing dev-set confirmation does not cover it.

**Corpus review was not needed.** The plan runs on `benchmarks/data/v8-corpus/reviewed_v3.json`, which the user reviewed on 2026-11-09 (34 reviewed labels, 4 disputed). The plan's `vocabulary.tsv` has sha256 `25f65f5a…`, equal to `examples/pbmc3k/vocabulary.sha256`. The earlier plan claim that this corpus is unreviewed was wrong.

Evidence: `table fallback: each verdict is exact equality of two label tables`. Figure inspected for the reading below: `skills/annotation-consensus/examples/pbmc3k/figures/fig1_umap.png` (consensus strict panel).

| Step | Result |
| --- | --- |
| `run scene`, `run plan` | pass, no blockers. The input is a copy with one added column, `annotation_scope = ROOT` (`pbmc3k_package_input.h5ad`); the original `pbmc3k.h5ad` is unchanged |
| Gate A | 43 of 44 assumptions identical to the AIDA run; the 44th differs only by the `consensus_directory` path |
| `collect-level`, `import-consensus`, `finalize` | pass |
| Imported `labels.tsv` against `consensus.tsv` | 8 of 8 clusters equal in `label` and `cl_id`; clusters 3, 4, 6 `Unknown` (rule `abstain`) |
| `cell_labels.tsv` | 2638 cells, each equals its cluster's consensus label |
| Review build and `pipeline packet` | pass; `report.html` 1.9 MB, 3 inline SVG, 0 `<img>` |

Reading of the figure: the strict panel follows the scanpy cluster names for B cells, both monocyte groups and megakaryocytes. The CD8, NK and dendritic-cell regions carry `Unknown`. PBMC3K has one donor and no truth labels, so this shows the workflow, not accuracy.

Status after the packet: `cell-curator critic import-findings --config pbmc3k.yaml --input FINDINGS.tsv`. Critic findings, `run validate`, the write-back preview and `write-back apply` stay with the user.

# Reruns with rule h4 (2026-11-17)

After the user dropped the tie-break (`plans/v8-agent-decisions.md`, 2026-11-17), both datasets were rerun with rule h4 and snap on as new package runs `v8-accept-aida-split-a-002` and `v8-accept-pbmc3k-002` (configs `aida_split_a_002.yaml`, `pbmc3k_002.yaml`). Gate A was approved by the user for both on 2026-11-17 (reviewer Mikhael Manurung); 43 of 44 assumptions equal run 001's, and the 44th is the annotation-rules text, which names the new consensus directory.

Evidence type: `table fallback: each verdict is exact equality of two label tables`.

| Check | AIDA split A | PBMC3K |
| --- | --- | --- |
| consensus under h4 (verify OK) | 19 of 23 clusters labelled, coverage 0.9072, 1702 cells flagged | 5 of 8 labelled, coverage 0.636088, 958 flagged |
| cluster `label`, `cl_id` against the h5 run | 23 of 23 equal | 8 of 8 equal |
| per-cell strict label against the h5 run | 8847 of 8847 equal | all cells equal (`cell_labels.tsv` identical) |
| `cell_flags.tsv` against the h5 run | identical | identical |
| per-cell best-guess label against the h5 run | 11 of 8847 differ (abstained cells in cluster 11 with 4 tied terms: h5 picked `platelet`, h4 picks a CD16-positive NK term by tool order) | not compared |
| `labels.tsv` and `cell_labels.tsv` of the package run 002 against 001 | equal | equal |
| `collect-level` | minutes (ranking cache hit) | 21 s |

Snap off was tested on AIDA and rejected: rule h4 with `snap_to_vocabulary: false` gave the cluster label `memory regulatory T cell`, which is outside the approved vocabulary and which `import-consensus` refuses.

Critic step (agent drafts, user approves; `plans/v8-agent-decisions.md`, 2026-11-17): the agent drafted 6 findings for PBMC3K and 13 for AIDA from `scores.tsv`, the consensus tool calls and donor counts (`critic_summary_*.tsv`, `critic-drafts/`). On the user's yes ("Keep all labels, record findings", then "Record as drafted") the findings were imported and reconciled with reviewer Mikhael Manurung: PBMC3K 3 RESOLVED and 3 DISMISSED, AIDA 7 RESOLVED and 6 DISMISSED; every RESOLVED line keeps the label and names the alternative. The strict report was rebuilt after reconcile (`report build --strict`) and `run validate` returned `complete` for both runs.

Cost: the stale `consensus_aida_config.yaml` pointed at the deleted worktree for the ontology; the path now points at the main checkout, and `verify` of the original h5 run still passes.

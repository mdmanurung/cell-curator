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
- **Open:** the container image (SIF `97febc139e1c…`, from `a7083bb`) predates `5c79dd2`. It must be rebuilt before merge.
- **Not checked:** accuracy. The run scores against split-A truth only if someone runs the skill's `score` step.

## Still the user's

Critic findings, `run validate`, `write-back preview` review and `write-back apply`.

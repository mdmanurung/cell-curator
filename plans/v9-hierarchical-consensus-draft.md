# Hierarchical consensus: design draft (future version, dev sets only)

Status: draft. Nothing here touches the frozen novel-set method. Every AIDA number below comes from scratch-script tables and has no inspected figure: `unsupported: no inspected figure`. Code anchors come from the reader maps; ac.py = skills/annotation-consensus/annotation_consensus.py. The anchors were not re-read for this draft.

## Answer

Yes, it can work, and the package already has most of the plumbing. But a hard top-down restriction fixes few errors on the one dev set measured (AIDA split A). The data show a different gain: a coarse back-off label for cells that abstain today. So the recommended design keeps the flat h4 call and rolls uncalled cells and clusters up the Cell Ontology to a coarse term. It does not narrow the search space, which is what you asked about. A true parent-restricted pass stays a later option, and only if a dev pilot justifies it. The package keeps assigning no label.

## Evidence from dev data

All numbers are AIDA split A. Inputs: truth.tsv (8847 cells), calls.tsv from acceptance-2026-11-17/consensus_aida_h4, and cl.obo via `CellOntology.is_descendant` (ontology.py:306). Scripts: /tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung-cell-curator/8e80509a-206c-44f5-b1b6-e34297fdb855/scratchpad/hier/01_lineage_errors.py and 02_backoff.py. Truth and calls join 1:1 on barcode.

Each strict call was compared with its truth CL term.
- Strict-called cells: 8026 (0.9072). Exact 3610 (0.4498), finer than truth 3165 (0.3943), coarser 56 (0.0070), divergent 1195 (0.1489). Abstained: 821 of 8847.
- Lineage-wrong calls: 104 cells (1.30% of strict calls) in 11 of 22 strict-called clusters. 97 of 104 are T vs NK (65 T call on NK truth, 32 NK call on T truth). Anchor scheme A, 5 anchors: T, NK/ILC, B, myeloid, precursor.
- Scheme B (7 anchors, monocyte, DC and platelet split out) gives 153 lineage-wrong cells (0.0191). The anchor set is my choice and needs your confirm.
- Fine-divergent but lineage-right: 1091 cells (0.1359). Only 191 of them have any tool naming a truth-consistent term. 988 (90.6%) already have all three tools in one lineage. 403 cells are a truth-vocabulary mismatch (effector memory CD8 call vs cytotoxic CD8 truth). A parent restriction cannot choose between siblings inside one lineage.
- Of the 104 lineage-wrong cells, 55 have all three tools unanimously in one wrong lineage. A tool-agreement check cannot catch them.
- Back-off: of 821 abstained cells, 442 have all three tools in one lineage, and 439 of those are in the truth lineage. That is 5.0% of all cells at 99% lineage accuracy. The best-guess lineage is right for 702 of 821 (0.8551).
- Search space: 38 corpus terms, 3 to 17 per lineage. The three tools already agree on lineage in most cases, so narrowing changes few votes.

Caveats. Labels are mostly cluster-level, so the effective n is about 22 clusters. AIDA PBMC has few lineages, so lineage errors will look rarer than in tissue. Truth is mixed-depth. I tested no hierarchical vote, only the existing flat calls.

Literature: I read five PubMed abstracts (scHPL PMID 33990598, popV PMID 39567746, OnClass PMID 34548483, CellTypist PMID 35549406, Azimuth reference PMID 34062119). None gives a hierarchical vs flat accuracy number, and I found no source on error propagation. This is background, not support.

## Recommended design: flat h4 plus ontology rollup, opt-in

The flat h4 call runs first and is unchanged. Every cell and cluster that passes keeps its fine label. Rollup acts only on cells and clusters that fail their gates.

Grafts: the Step 0 pilot and per-cluster reporting from the strict design; the `outside_parent` idea (a tool term outside the anchor lineage abstains instead of snapping); the `parent_conflict` flag as a diagnostic; mask-not-renormalise (keep M+w); the stay-at-parent idea for coarse tool terms.

Steps:
1. Pilot, no code, about 0.5 day. On existing flat runs for AIDA split A, AIFI dev and PBMC3K, count abstained cells whose tools agree on one lineage. Make the figure and inspect it (needs a grill plan). Stop if the back-off gain is small.
2. Rollup table, a new input `inputs.rollup_terms` (label, cl_id, tier). Rows are the only coarse terms a rolled label may take. Each row must resolve in cl.obo and have a vocabulary descendant-or-equal (check in step_validate, ac.py:521-591, pattern of load_vocabulary ac.py:457-467). Dev runs use unreviewed rows under `annotation.development_candidates` and can never be written back (AGENTS.md gate 3).
3. Config block `consensus.rollup` {enabled, terms, cell_gate, cluster_gate{share, lead}, anchor, tiers}, absent by default. Add it to the schema (ac.py:127-166) and the IMPACTFUL set (ac.py:109-124). canonical_config (ac.py:369-385) must omit the block when absent, so frozen config hashes do not move.
4. Tool mapping (step_map_tools, ac.py:802-865). With rollup on, a tool term that is in the rollup table keeps its identity (new status `kept_rollup_term`). It is not snapped down. Today a coarse term snaps down to the nearest descendant (ac.py:821-826), which invents depth. A tool term outside the flat best guess's lineage abstains (new status `outside_parent`) and is not snapped. Test on dev data that this does not change cells that passed the flat gate. Preferred form: apply it only to cells that would abstain, so the flat path stays identical.
5. Term universe (step_consensus, ac.py:868-870, :982-986). Rollup terms join the sorted list with S=0. They carry tool votes only.
6. Flat call, unchanged (ac.py:994-1008).
7. Cell rollup, uncalled cells only. support(a) = sum of vote mass over a and its descendants, via the DESC matrix (ac.py:880-886), divided by M+w (denominator fixed, ac.py:900, :934). Take the deepest a with support at or above `rollup.cell_gate`. Depth is the shortest is_a path (up_dist, ac.py:782-794). Ties go to deeper, then name order. With `anchor` on, only ancestors of the flat best guess are tried. Otherwise the cell stays Unknown.
8. Cluster rollup, only for clusters that fail the cluster gate (ac.py:1045-1056). share(a) is the fraction of cluster cells whose final term is a or below. lead(a) is share(a) minus the best share outside a's subtree. Take the deepest a above `cluster_gate.share` and `.lead`. New gate_failures strings: `rollup_share`, `rollup_lead`, `no_rollup_term`. decision_rule is `rolled_h4`.
9. Diagnostics. Flag `rolled_off_anchor` when the fine best guess is not under the rolled label. Flag `parent_conflict` when the flat cluster call passes the cluster gates but lies outside the rolled subtree. Both need a critic finding.
10. Outputs. consensus.tsv gains resolution (fine or rolled), depth, label_fine_best_guess, cl_id_fine_best_guess, rolled_from, and tier columns, all written by the skill. The import's required columns (consensus_import.py:20-29) stay a subset, so extra columns are safe. calls.tsv and report.md get a depth table. Rollup files join hashed_outputs (ac.py:1338-1348) only when enabled, so expected_hashes.json (18 files) stays valid.
11. Keep `hierarchy_levels = [L1]` (config.py:494, :558-562). Tiers are extra columns in a sidecar levels/L1/rollup.tsv, because labels.tsv has a fixed column set (hierarchy.py:15-30) and new levels would need per-level scores and evidence (annotation.py:687-696). Whether finalize can read the sidecar is not yet inspected.
12. Package side. Import accepts a label in the L1 vocabulary slice (consensus_import.py:123-127) or in the Gate A approved rollup table. It checks tier ancestry with `is_descendant` (ontology.py:306-312) and `ancestors` (ontology.py:222-233). It assigns nothing.
13. Tests: a rolled label accepted, a label outside vocabulary and table refused, a non-ancestor tier refused. Run `.venv/bin/python -m pytest -q`, `ruff check .`, `mypy`.

## Skill changes

- ac.py:109-166: rollup config keys, all in the IMPACTFUL set.
- ac.py:457-467: rollup-table loader beside load_vocabulary.
- ac.py:521-591: step_validate checks (rows resolve, each has a vocabulary descendant, tiers are ancestor-or-self).
- ac.py:802-865: `kept_rollup_term` and `outside_parent` statuses in step_map_tools.
- After ac.py:1008 and :1056: cell and cluster rollup.
- ac.py:1180-1203: new columns; ac.py:1338-1348 and :1614-1653: hash and verify coverage.
- skills/annotation-consensus/SKILL.md: document the block, table format, new strings, depth table, and that the flat path is unchanged. Label every default as demo, not approved.
- A separate demo config with its own hashes; examples/demo/config.yaml and its 18 pinned hashes stay untouched.

## Package changes

- consensus_import.py:123-127 (label filter), :20-29 (columns unchanged): accept rollup-table labels, write levels/L1/rollup.tsv, copy the table into the receipt.
- config.py: new `annotation.rollup_terms` path and optional tiers. The L1-only guards stay (config.py:558-562, pipeline.py:1015-1018, consensus_import.py:94-99).
- hierarchy.py:191-242: reuse validate_hierarchy_nesting if the tier table fits; else add a small ancestor-or-self helper.
- annotation.py near :712-717 and :763: carry tier columns into the cell table and report (inspect first).
- guided.py:942-960 and :1076-1092: add the table to the Gate A scope hash and to the next-command text.
- tests/test_consensus_import.py: three new cases.

If a true parent-restricted pass is wanted later, the anchors are: lift the guards at config.py:558-562, pipeline.py:1015-1018 and consensus_import.py:94-99; one consensus directory per (level, parent); replace the equal-cells check (consensus_import.py:47-86) with a subset check; write the real parent_scope (consensus_import.py:146); restrict `vset` (ac.py:808); build kNN within the parent (ac.py:771-775); fix the per-level state read (annotation.py:720, :739). That needs an L2 vocabulary that does not exist (reviewed_v3.json holds 110 labels, all L1; corpus schema supports it, corpus.py:104-135, :253-260). Estimate 10 to 13 days.

## Impactful choices you must decide

Each is a plan row for you. No default exists unless stated.
1. Build at all after the Step 0 pilot. Evidence: lineage errors 1.30%, back-off 5.0% of cells (no inspected figure).
2. Lineage anchor set: scheme A (5) or B (7). Inherited from my scratch scripts; one-line confirm.
3. Rollup term set and depths. No default. Needs your corpus review.
4. Cell rollup gate. New knob; no default. The flat gate 0.53 was chosen on L1 only (examples/demo/config.yaml:44).
5. Cluster rollup gates. No default. Flat values 0.5 / 0.2 / 0.40 (config.yaml:51-53) were set on L1 only.
6. Anchor rule: ancestors of the flat best guess only, or any ancestor of the vote mass. No default.
7. Coarse tool terms unsnapped, or snapped as today. No default; it changes the tool matrix.
8. Marker mass in rollup support (marker_weight*S included, or tools only). No default.
9. Support denominator: fixed M+w (inherited, ac.py:900, :934; one-line confirm) or renormalised.
10. Tier definitions, and Unknown vs root-tier fallback when no rollup term passes.
11. Per-level gates and marker settings for any later L2 pass: z-score reference (all cells, ac.py:648-653, vs parent only), kNN scope (ac.py:771-775), share denominator after masking, tie rule (h4 vs h5), abstained-parent handling. No package default for a parent pass; each needs a named dev yardstick.
12. Dev yardstick and pass bars (below). Any seeds or sensitivity grids: no recommended value; the method has no seed, since the kNN is deterministic.

## Risks

- Error propagation. Rollup cannot fix the 55 of 104 lineage-wrong cells where all tools agree on the wrong lineage. A hard top-down gate would send all 104 to the wrong subtree with no recovery. This is why the winner leaves passing calls untouched.
- Unsnapped coarse terms change the tool matrix and can flip winners for cells that passed the flat gate. Report this. Prefer applying it only to abstaining cells.
- Rolled labels raise coverage and lower information. Never quote coverage without depth.
- The ontology is a DAG, so the lowest common ancestor is not unique. The support rule and tie order decide, and that is an impactful choice.
- Rolled terms have no marker program or scores.tsv row. The evidence gate cannot score them, and critics must attest that the label rests on descendants' evidence. finalize and the sidecar path are not inspected.
- Widening the import filter (consensus_import.py:123-127) loosens one check. Tie it to a Gate A approved table.
- AIDA has few lineages, and truth is mixed-depth (403 cells of vocabulary mismatch), so exact-match scores mislead for any arm.
- A config default injected into canonical_config would break every frozen hash. Test first.
- Any package-side projection of tiers would break the rule that the package assigns no label.
- Unverified: Snakemake DAG (`annotation_workflow_path`), finalize sidecar handling, and h5 behaviour under a mask.

## Dev-set validation and pass criterion

Dev sets only: AIDA split A (654-cell demo and the 8847-cell set), AIFI dev split, PBMC3K. No file under cc-v8-work/novel-2026-11-17 and no novel truth.
1. Regression, rollup off: the 18 hashes in expected_hashes.json and the lock hash verify. SKILL.md:305-308 values match exactly (AIDA h4 snap off 0.848617 at 0.9072; AIFI h4 0.899894 at 0.988437).
2. Rollup on, with your approved table and gates. Report fine accuracy on flat-called cells, coverage by depth, rolled-cell lineage accuracy, lineage-wrong count among rolled cells, and the number of flat-passing cells whose winner changed. Report per cluster (about 22 on AIDA), not only per cell.
3. Sensitivity runs only on grids you name, one choice at a time, each recorded.
4. Error-propagation check on the AIDA T vs NK clusters 3 and 6.
5. Every conclusion cites a saved, inspected figure by path (plot-first). Until then it is written `unsupported: no inspected figure`. Making the plots needs a grill plan.
6. `verify` with the scrubbed-env rerun (ac.py:1582-1586, :1614-1653) on the rollup demo.

Proposed pass criterion (you set the actual bars): fine accuracy on flat-called cells is unchanged on all three sets; rolled-cell lineage accuracy is at or above a bar you set (AIDA baseline 99% on the 439 cells, unsupported: no inspected figure); measurable coverage at coarse depth; zero flat-passing cells with a changed winner, or a count you accept.

## Effort

About 6 days: pilot 0.5, skill 3, package 1.5, tests and dev validation 1. This omits the uninspected finalize and sidecar change. A true parent-restricted pass adds 10 to 13 days and needs a new L2 vocabulary and your corpus review.

## What stays unchanged

- Rule h4, marker_weight 1.0, gate 0.53, cluster gate 0.5 / 0.2 / 0.40 and snap_to_vocabulary true stay as frozen for the running novel-set benchmark. The new block is absent by default, so the flat path stays byte-identical.
- The package imports, checks and finalizes only. It assigns, edits and overrules no label.
- Gate A (you approve assumptions, including the rollup table), critics, corpus review (your OK, your name as `--reviewer`) and the write-back token (you type it) stay as they are. No agent runs them.
- `hierarchy_levels` stays `["L1"]`.
- Any novel-set use is a separate future version with a new frozen protocol and your decision.

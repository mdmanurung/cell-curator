---
name: annotation-consensus
description: >-
  Deterministic, rerunnable consensus of standalone cell-type annotation tools (CellTypist,
  popV, SingleR, scANVI) with cell-curator marker-program evidence. Writes a strict label,
  a best-guess label flagged low-confidence, and the reason for each cluster and each
  flagged cell. Local only; truth labels are read by the optional score step alone.
version: 0.1.0
author: Mikhael Manurung
domain: single-cell annotation
license: see repository
tags:
  - single-cell
  - cell-type-annotation
  - consensus
  - provenance
  - reproducibility
  - cell-curator
inputs:
  - name: config
    type: file
    format: [yaml]
    description: ONE strict config naming every input and every impactful choice
outputs:
  - name: report
    type: file
    format: [md]
  - name: result
    type: file
    format: [json]
  - name: provenance_tables
    type: file
    format: [tsv]
demo_data:
  - path: examples/demo/config.yaml
    description: 654 cells from AIDA split A, 23 clusters, three tools (demo values, not approved defaults)
---

# annotation-consensus

## 1. Why this exists

Standalone annotation tools disagree with each other and give no reason for a label.
All three tools agree on 70% of AIDA cells and 83% of AIFI cells, yet they are wrong on
7.6% and 4.4% of all cells (plan `plans/v8-consensus-mode.md`), so agreement does not
prove a call. This pipeline takes the per-cell calls of several tools and the marker-program
votes of cell-curator, builds one label per cell and per cluster, and records why.

The pitch is an ontology-checked label with a recorded reason. Accuracy is not the pitch:
the measured gain over plain majority is about 1 point (section 7, results).

The run is deterministic. The same inputs and config give byte-identical hashed outputs, and
`verify` checks that.

## 2. Core capabilities

- Maps each tool's labels to the corpus vocabulary along the Cell Ontology (`snap_to_vocabulary`).
- Votes per cell: each tool counts one vote, the marker evidence counts `marker_weight`.
- Three rules: `majority`, `h4` (marker co-voter), `h5` (co-voter plus ontology tie-break).
- Two label columns: strict (`Unknown` kept) and best guess (flagged low-confidence where strict is `Unknown`).
- Cluster label from the plurality of strict per-cell calls, behind a share and a lead gate.
- `consensus.tsv`: one row per cluster with the decision rule, each tool's call, marker votes, gate
  failures, rationale and rejected alternatives.
- `cell_flags.tsv`: only the cells that depart from their cluster label, abstain, take a fallback call
  or reach a tie-break, each with a reason.
- Optional per-cell marker votes computed from an expression matrix (a port of the
  `score_marker_programs` loop in `src/cell_curator/evidence.py`, with a built-in check).
- `verify` hashes the stored outputs, reruns in a fresh subprocess with a scrubbed environment and compares all hashed outputs.

## 3. Input formats

All inputs are named in one YAML config. Relative paths resolve against the config file.

| key | format | content |
| --- | --- | --- |
| `inputs.cells` | h5ad | `obs[cluster_key]`, `obsm[smoothing.embedding_key]`; `X` only when votes are computed |
| `inputs.vocabulary` | TSV | `label`, `cl_id` (plus marker columns when votes are computed) |
| `inputs.ontology` | `cl.obo` | referenced by path, never copied |
| `inputs.tools.<name>` | TSV | `barcode`, `predicted_label`, `predicted_cl_id`, one file per tool |
| `votes.path` | TSV | `barcode`, `cluster_id`, `top1`, `top2`, `score1`, `score2`, `voted` (when `votes.mode: precomputed`) |
| `inputs.markers` | TSV | optional, `levels/L1/markers.tsv` of a cell-curator run |
| `scoring.truth` | TSV | optional, `barcode`, `cell_type`, `cell_type_ontology_term_id`; the score step only |

Barcodes are aligned by name, never by position. A barcode mismatch fails `validate`.
The tools themselves are not run here. Their predictions are inputs, and their sha256 is recorded.

## 4. Workflow

1. `validate`: inputs exist, barcodes match by name, labels map, required keys present, no unknown keys.
2. `votes`: per-cell marker votes, computed (with the built-in check) or read from the supplied file.
3. `smooth`: kNN smoothing of the one-hot top-1 votes.
4. `map-tools`: tool terms to the corpus vocabulary, written to `tool_term_map.tsv`.
5. `consensus`: vote, gate, tie-break; writes `calls.tsv`, `consensus.tsv`, `cell_flags.tsv`.
6. `markers`: tidy table of genes in the top N under both Wilcoxon and logistic regression.
7. `score` (optional): the only step that reads truth.
8. `report`: `report.md`, `result.json` and the `reproducibility/` bundle.

`run` does all steps into an empty directory. Without `--force` it refuses to overwrite. `--force` removes outputs
only from a directory this pipeline created (it holds a `.annotation_consensus_output` marker file) and refuses
any other directory. An `--output` inside a cell-curator run directory (one holding `run_state.json` or
`annotation/levels`) is refused.

## 5. CLI reference

Checked against `annotation_consensus.py --help` and each subcommand's `--help`.

```
usage: annotation_consensus.py [-h] [--input INPUT] [--output OUTPUT] [--demo]
                               [--force] [--score]
                               {validate,votes,smooth,map-tools,consensus,markers,score,report,run,verify}
                               ...
```

Options, the same on the top-level command and on every subcommand:

| option | meaning |
| --- | --- |
| `--input CONFIG`, `--config CONFIG` | the one strict config YAML |
| `--output DIR` | output directory; default is `run.output_dir` of the config |
| `--demo` | use the bundled demo config `examples/demo/config.yaml` |
| `--force` | overwrite existing outputs |
| `--score` | also run the score step (needs the `scoring` block and truth) |
| `--allow-version-drift` | run on versions that differ from `requirements.lock`; the drift is recorded in `run_manifest.json` and `result.json` |

Subcommands:

| subcommand | does |
| --- | --- |
| `validate` | check inputs exist, barcodes match by name, labels map, required config keys |
| `votes` | port of the marker-program vote loop; skipped when votes are supplied precomputed |
| `smooth` | kNN smoothing of one-hot top-1 votes over the embedding |
| `map-tools` | snap tool labels to the corpus vocabulary along the ontology |
| `consensus` | tool plus marker vote, gate, tie-break; writes calls, consensus, cell_flags |
| `markers` | tidy table of genes in the top N under wilcoxon and logreg (enrichment NOT built yet) |
| `score` | optional: score against truth with `benchmarks/scripts/score_v8.py` (the only step that reads truth) |
| `report` | `report.md`, `result.json` and the reproducibility bundle |
| `run` | run every step into an empty output directory (the default when no subcommand is given) |
| `verify` | hash the stored outputs, rerun in a fresh scrubbed subprocess, compare all hashed outputs |
| `expected-hashes` | write `expected_hashes.json` from a finished demo run in `--output` (`--expected-file` picks the target) |

Exit codes: 0 on success; 2 on a `PipelineError` (a message starting `error:` goes to stderr).
`verify` exits 2 and lists each file that differs, including a stored file that was edited after the run.
Every step except `validate` needs `validation.json` in the output directory, so run `validate` first
when you run steps by hand.

Examples:

```bash
python skills/annotation-consensus/annotation_consensus.py --demo --output out/demo
python skills/annotation-consensus/annotation_consensus.py run --input CONFIG.yaml --output OUT [--score] [--force]
python skills/annotation-consensus/annotation_consensus.py verify --input CONFIG.yaml --output OUT
```

Snakemake runs the same steps through the same CLI. Exactly one `--configfile` is allowed, and the
config must hold `run.output_dir`:

```bash
snakemake -s skills/annotation-consensus/Snakefile --configfile CONFIG.yaml -c1
snakemake -s skills/annotation-consensus/Snakefile --configfile CONFIG.yaml -c1 score   # optional, reads truth
```

The Snakefile fails on any config overlay (`--config key=value`). To show scores in the report, run the
`score` target first, then the default target.

### Strict config

`config.example.yaml` lists every key. A missing required key stops the run, and so does an unknown key.
These choices are impactful, are the user's, and have no default:

| key | meaning |
| --- | --- |
| `consensus.rule` | `majority`, `h4` or `h5` |
| `consensus.marker_weight` | marker voter weight, in units of one tool vote |
| `consensus.gate` | winning share needed for a strict label |
| `consensus.tie_break_delta` | vote-split width in tool votes (rule `h5`) |
| `consensus.snap_to_vocabulary` | move tool terms to the nearest vocabulary term (`true` or `false`) |
| `consensus.marker_scope` | `cell` or `cluster`: marker share per cell or cluster mean |
| `consensus.tool_order` | best tool first; breaks exact ties |
| `consensus.values_status` | who approved the values, or the demo label |
| `smoothing.k`, `smoothing.embedding_key` | kNN size and embedding |
| `cluster_gate.share`, `.lead`, `.fallback_margin` | cluster label gate and marker-only fallback |
| `markers_table.top_n`, `.comparator_role` | marker table (required when `inputs.markers` is set) |

The demo config carries values copied from the frozen results and labels them
`demo values, not approved defaults`. Choose your own values and say who approved them in
`consensus.values_status`. `determinism.threads` and `determinism.float_digits` are also required keys.

## 6. Demo

```bash
.venv/bin/python skills/annotation-consensus/annotation_consensus.py --demo --output out/demo
.venv/bin/python skills/annotation-consensus/annotation_consensus.py verify --demo --output out/demo
```

The demo is 654 cells from AIDA split A that cover all 23 `bench_cluster` values, drawn with a seeded
stratified draw by `examples/make_demo.py`. It has three tool prediction files, precomputed marker
votes, a marker table and a truth file. Truth is read only when you add `--score`.

A run on this machine labelled 19 of 23 clusters and flagged 117 cells (rule `h5`, coverage 0.899083).
With `--score` the same run gave consensus_strict 0.82483 at 0.899083 and 19 hashed files (checked here).
That fixture is small, so do not compare its accuracy to the frozen numbers in section 7.

## 7. Algorithm

The logic is a port of the offline scripts under
`cc-v8-work/tusoskill/cov95-2026-11-11/candidates/` (`c004`, `c009`, `c010`, `c012`, `c014`).

**Marker votes per cell (`votes`).** A copy of the `score_marker_programs` loop. With `votes.mode: compute`
the step compares its `raw_vote_fraction` with the run's `levels/L1/scores.tsv` (`votes.compute.check_scores`)
and fails when the maximum absolute difference is 1e-9 or more. The builder reports a maximum of 1.11e-16 on
AIDA. With `mode: precomputed` the supplied file is read and no expression matrix is needed.

**Smoothing (`smooth`).** The one-hot top-1 vote of each cell (weighted by `voted`) is averaged over its
`k` nearest neighbours in `obsm[embedding_key]`, with the cell itself included (brute-force kNN, `k + 1`
points). The result is the marker share of each vocabulary term per cell.

**Tool mapping (`map-tools`).** Each tool term moves to the nearest vocabulary term on its ontology lineage
(fewest `is_a` steps up or down; the first sorted term wins a tie). Terms with no vocabulary term on their
lineage stay as they are. With `snap_to_vocabulary: false` no term moves. Every move is listed in
`tool_term_map.tsv`.

**Per-cell vote (`consensus`).**

- `majority`: each tool gives one vote to its term. The best tool in `tool_order` breaks exact ties.
- `h4`: the marker share (per cell or per cluster mean, by `marker_scope`) is added with weight `marker_weight`.
- `h5`: as `h4`, but where the top total and the runner-up differ by `tie_break_delta` or less, break the
  tie in this order: (1) a tied term that descends from another tied term, `tie_break_specificity`;
  (2) a unique top marker share, `tie_break_marker`; (3) the best tool, `tie_break_tool`.
  A tie nobody backs becomes `abstain`.
- Winning share is the winner's total over (number of tools + `marker_weight`). A share below
  `consensus.gate` gives `Unknown` in the strict column; the best-guess column keeps the winner and
  is flagged low-confidence.
- `decision_rule` per cell: `all_agree`, `majority`, `tie_break_specificity`, `tie_break_marker`,
  `tie_break_tool`, `abstain`.

**Cluster label.** The plurality of strict per-cell calls, if its share is at least `cluster_gate.share`
and its lead over the runner-up is at least `cluster_gate.lead`. Otherwise the cluster is unresolved and
its strict label is `Unknown`. The cluster `decision_rule` is `all_agree` when every tool's cluster majority
agrees with the label, `abstain` when the gate fails, else the most common cell rule among its carriers.

**Marker-only call (comparator).** The c009 rule: the cluster takes the top term of the mean marker share
when its share and lead pass `cluster_gate`; otherwise a cell takes its own top term when its smoothed
margin is at least `cluster_gate.fallback_margin`. It is written per cell and scored as a comparator.
It never feeds the consensus.

**Flags.** A cell is flagged when its call departs from its cluster label, abstains in a labelled
cluster, takes a fallback call in an unresolved cluster, or reaches a tie-break.

**Determinism.** Single process. `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` are set from
`determinism.threads` before numpy loads (the demo uses 1). The thread count is in no hashed file and not in
the config hash: it lives in `run_manifest.json`. Floats in every table are written with
`determinism.float_digits` significant digits. Hashed outputs hold no timestamp, host name, user name or path,
only file names and sha256; those live in `run_manifest.json`, which is excluded from the hash comparison.
It lists the sha256 of every input (config, cells, votes, tool predictions, vocabulary, ontology, markers, and
truth when scored), the version lock and any allowed version drift.

- **Pinned environment.** `requirements.lock` lists python and the package versions (numpy, pandas, scipy,
  scikit-learn, anndata, pyyaml, pyarrow, cell-curator). A run on any other version stops with an error that
  names it, unless you pass `--allow-version-drift`, which is recorded in `run_manifest.json` and `result.json`.
  `reproducibility/environment.txt` is hashed, so `verify` and `expected_hashes.json` pass only on the same
  pinned environment. It also holds the line `requirements.lock-sha256 <hex>` (the same hash is in
  `run_manifest.json` as `lock_file_sha256`). `verify` stops with "requirements.lock changed since the run" when
  the lock file on disk differs from that line, before it reruns anything, unless you pass `--allow-version-drift`.
- **Ties.** Scores are rounded to 1e-9 before every argmax and argsort. A tie goes to the first label name
  (ascending), because terms are indexed in label-name order. kNN neighbours are ordered by distance, then by
  barcode rank: `smooth` searches on the rows sorted by barcode, so an exact distance tie goes to the lower barcode
  and the neighbour set does not depend on the input row order or the thread count (tables keep the input order).
  Distances are recomputed in float64 from the pair of points alone, so two points tie only when their float64
  squared distances are equal bit for bit; there is no tolerance. The candidate window widens until it ends more
  than 1e-10 of the largest squared norm beyond the k-th neighbour. The kNN search itself runs on one thread.
- **`verify`.** It hashes the stored output files, then reruns in a fresh subprocess whose environment holds only
  `PATH` and the pinned thread variables, and exits 2 on any difference in either.
- **`expected_hashes.json`.** The committed sha256 of every hashed output of the demo run, written by
  `annotation_consensus.py expected-hashes`. A test fails when an output changes until the file is regenerated on
  purpose.

### Results (measured offline, frozen settings)

Descendant accuracy among called cells, coverage in brackets, scored with `benchmarks/scripts/score_v8.py`.
Gate and tool order were chosen on one dataset and applied unchanged to the other.
Sources: `c012_harmonizer/frozen.tsv`, `c014_tiebreak/frozen.tsv`, `c009_smooth_cluster/sweep_summary.tsv`
(AIDA marker-only row), `c010/summary.tsv` (AIFI marker-only row) and the design plan `plans/v8-consensus-mode.md`.
The two marker-only rows sit at different coverage levels, so read each against its own coverage.

| design | AIDA | AIFI |
| --- | --- | --- |
| best single tool (popV) | 0.835 (0.973) | 0.888 (0.998) |
| plain majority, no gate | 0.818 (1.000) | 0.894 (1.000) |
| majority + vote-share gate | 0.841 (0.945) | 0.894 (1.000) |
| marker evidence as equal co-voter (H4, per cell, snap off) | 0.8486 (0.9072) | 0.8999 (0.9884) |
| co-voter plus tie-break (H5, per cell, delta 0, snap on) | 0.8441 (0.9072) | 0.8920 (0.9885) |
| marker evidence alone, per cluster, share 0.5, lead 0.2, fallback 0.4 (c009, c010; not at 0.95 coverage) | 0.850055 (0.724426) | 0.852398 (0.701222) |
| marker evidence alone, per cell, at coverage 0.95 (from the plan; a different coverage level from the row above) | 0.730 (0.95) | 0.695 (0.95) |

Figures inspected for these rows:
`cc-v8-work/tusoskill/cov95-2026-11-11/candidates/c012_harmonizer/figs/riskcov_aida_split_a.png`,
`.../c012_harmonizer/figs/riskcov_aifi_dev.png` and `.../c014_tiebreak/figs/riskcov_aida_split_a.png`.
In the first two, the H4 per-cell and per-cluster curves lie at or above the plain-majority and
single-tool points at coverage near 1. The referee curves are lower, and the per-cluster referee is
the lowest on AIFI. In the third, the H5 curves and the H4 per-cell curve overlap.

The pipeline reproduces the H4 and H5 per-cell rows, scored with `score_v8.py` on the full inputs
(builder's report; the full-data regression tests are marked slow and were not run here):

| setting | AIDA | AIFI |
| --- | --- | --- |
| rule `h4`, snap off | 0.848617 at 0.9072 (file: 0.8486 at 0.9072) | 0.899894 at 0.988437 (file: 0.8999 at 0.9884) |
| rule `h5`, delta 0, snap on | 0.844132 at 0.9072 (file: 0.8441 at 0.9072) | 0.892044 at 0.988502 (file: 0.8920 at 0.9885) |

The two rows differ because c012 did not snap tool labels and c014 did, so `snap_to_vocabulary` is a required key.

Statements the measurements support:

- Marker evidence as an equal co-voter is the best design at coverage 0.95 or higher. It adds about 1 point
  over plain majority on AIFI and nothing on AIDA.
- The referee role is lower, and per cluster it is unstable (AIFI 0.763 at 0.970).
- The tie-break adds nothing over the co-voter. Per cell it is identical. Per cluster the paired difference is
  +0.0024 [-0.0085, 0.0140] on AIDA and 0.0000 [-0.0018, 0.0015] on AIFI (`c014_tiebreak/boot.tsv`).
  With three tools and one marker vote, a winning share of 0.5 needs two tool votes, so the marker vote cannot flip a majority.
- Within the tie-break, the specificity step is right more often than wrong but covers only 31 to 727 cells.
  The marker-share step is net wrong on AIFI (1,235 right, 3,244 wrong, per cell, delta 0.5; from the plan).
  Whether to keep that step is open.
- No design reaches 0.85 accuracy at coverage 0.95 or higher on AIDA. The target of 0.95 coverage at 0.85 accuracy was not reached there.
- Agreement is not proof: the three tools agree and are wrong on 7.6% (AIDA) and 4.4% (AIFI) of all cells.

Caveats:

- The comparators were trained on AIFI. AIFI is in-domain for them and must be read that way.
- AIDA truth is coarse, so descendant credit is generous.
- One weight grid and one ontology-snap rule were tried.
- scANVI ran on AIFI only; it failed on AIDA with a gene-count mismatch. The benchmark tool set is CellTypist, popV and SingleR.
- No public-pretrained model was run.
- The demo fixture (654 cells) is too small to compare with these numbers.

## 8. Example queries

- "Build a consensus annotation for this dataset from my CellTypist, popV and SingleR outputs."
- "Why did cluster 7 get Unknown, and what would the best guess be?"
- "Which cells disagree with their cluster label, and which tools disagreed?"
- "Does the result change if I turn `snap_to_vocabulary` off?" (a new config, a new output directory)
- "Verify that the last run reproduces byte for byte."

An impactful choice is the user's. Agents that run this skill ask for the config values and do not pick them.

## 9. Output structure

```
OUT/
  report.md                     human-readable summary, parameters, tables
  result.json                   machine-readable summary (valid JSON); has "scores" only after --score
  calls.tsv                     one row per cell: strict and best-guess label, decision_rule, tool labels, marker-only call
  consensus.tsv                 one row per cluster: label, decision_rule, tool_calls, marker_votes,
                                gate_failures, gate_flags, rationale, rejected_alternatives, n_cells_flagged
  cell_flags.tsv                only flagged cells: barcode, cluster_id, cell_label, cell_best_guess, cluster_label,
                                reason, decision_rule, winner_share, tool_labels
  tool_term_map.tsv             every tool term and where it was mapped
  marker_table.tsv              genes in the top N under both methods (when inputs.markers is set)
  scores.tsv                    only after --score
  validation.json
  work/                         intermediate files (votes status, smoothed votes, mapped tools, summaries)
  reproducibility/
    commands.sh                 rerun commands; set CONFIG to the strict config
    environment.txt             Python and package versions (hashed: verify passes only on the pinned environment)
    input_checksums.sha256      sha256 of every input (file names, not paths) and of the canonical config
    checksums.sha256            sha256 of every hashed output
  run_manifest.json             timestamp, host, argv, paths, thread count, input hashes, version drift;
                                excluded from the hash comparison
  .annotation_consensus_output  marker: this pipeline created the directory (--force needs it); not hashed
  score/                        score_v8.py working files; excluded from the hash comparison
```

The demo run wrote 18 hashed files in `checksums.sha256` without `--score`; the score step adds `scores.tsv`.
`consensus.tsv` columns are `cluster_id, label, cl_id, label_best_guess, cl_id_best_guess,
best_guess_low_confidence, decision_rule, tool_calls, marker_votes, marker_only_call, marker_gate_pass,
gate_failures, gate_flags, rationale, rejected_alternatives, n_cells, n_cells_flagged, cluster_cell_share,
cluster_cell_lead`.

## 10. Dependencies

Python 3.11 with numpy, pandas, scipy, scikit-learn, anndata, pyarrow, PyYAML. Snakemake for the optional DAG.
The repository `.venv` has all of them; versions used in the checked run are in `reproducibility/environment.txt`
(numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, scikit-learn 1.9.0, anndata 0.12.19, pyarrow 25.0.0, PyYAML 6.0.3);
`requirements.lock` pins them with python 3.11.12 and cell-curator 0.2.0, and a run refuses other versions.
The ontology file `cl.obo` and, for `score`, `benchmarks/scripts/score_v8.py` are referenced by path.
The pipeline installs nothing. It imports `cell_curator.ontology` (the `cl.obo` reader) and, when it computes votes, `cell_curator.contracts`, so the `cell-curator` package (0.2.0 in the checked run) must be importable. `src/` is read, not changed.

Not built:

- Enrichment of marker genes against known markers (PanglaoDB, Azimuth PBMC, package corpus). The marker table lists genes only.
  The offline analysis in `candidates/c015_markers` is not part of this pipeline.
- The public-pretrained sensitivity run (cached public CellTypist models and any public popV and SingleR models).
- Any change to `src/`. The package changes listed in the plan (`config.py`, `annotation.py`, `hierarchy.py`,
  `reconciliation.py`, `review.py`, `provenance.py`) are not made.
- Determinism rows 8 (login node vs compute node) and 12 (network blocked, resource-file hash checks) of the plan,
  and the Apptainer image. Rows 1 to 7 and 9 to 11 have tests (`tests/test_annotation_consensus.py`, section 7).
- Running the annotation tools themselves. Their predictions are inputs.

## 11. Safety

- Offline. The pipeline makes no network call and sends nothing off the machine.
- No write-back. The outputs are evidence tables. A consensus label never writes back by itself.
  `cell-curator write-back apply` stays a person's step, with the approval token typed in a terminal.
- No gate bypass. The skill does not run `assumptions approve`, `corpus review`, critic attestation,
  `run validate` or `run execute`. Gate A, the critics and corpus review are untouched.
- Strict config. Every impactful choice is a required key. A missing key or an unknown key stops the run.
  The demo values are labelled `demo values, not approved defaults`.
- Truth only in the score step. Truth labels enter no other step, and scoring does not change the consensus outputs.
- No overwrite. The run refuses an existing output directory without `--force`, and `--force` removes only what a
  directory marked by this pipeline holds. It refuses an `--output` inside a cell-curator run directory.
- LLM evidence stays evidence. This pipeline uses no LLM. It does not read `cell-curator llm packet` output.
- Outputs hold cluster and cell labels and per-cell barcodes. Treat them as the input data's sensitivity level.

## 12. Integration

- Reads: the h5ad and the `levels/L1/` outputs of a cell-curator run (`markers.tsv`, and `scores.tsv` for the vote check).
- Feeds: the consensus tables are evidence for review. They do not enter a run on their own.
  Moving them into `consensus.tsv` and `cell_flags.tsv` under `levels/L1/` of a run needs the package changes in
  section 10, which are not made.
- Scoring uses `benchmarks/scripts/score_v8.py` from this repository.
- Design plan: `plans/v8-consensus-mode.md`. Tests: `tests/test_annotation_consensus.py`
  (`.venv/bin/python -m pytest skills/annotation-consensus/tests -q` runs the fast suite, about 6 min;
  `ANNOTATION_CONSENSUS_FULL=1` adds the full-data AIDA regression, about 15 min in total, so submit it with `sbatch`).
- Snakemake: `Snakefile`, one rule per step, same CLI.

## 13. Citations

Tool references below were written from memory and not checked against the sources.

- Domínguez Conde C et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans (CellTypist). Science 2022. [CITATION NEEDED: confirm volume and pages]
- Ergen C et al. Consensus prediction of cell type labels in single-cell data with popV. Nature Genetics 2024. [CITATION NEEDED: confirm]
- Aran D et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage (SingleR). Nature Immunology 2019. [CITATION NEEDED: confirm]
- Xu C et al. Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models (scANVI). Molecular Systems Biology 2021. [CITATION NEEDED: confirm]
- AIDA and AIFI datasets: [CITATION NEEDED: the data sources used in `benchmarks/`]
- Cell Ontology (`cl.obo`): [CITATION NEEDED: version used in `benchmarks/data/aida-phase1-v2-pbmc-hvg2000/`]
- Offline measurements: `cc-v8-work/tusoskill/cov95-2026-11-11/candidates/` and `plans/v8-consensus-mode.md` in this repository.

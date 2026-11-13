# annotation-consensus

A deterministic consensus of standalone annotation tools (CellTypist, popV, SingleR, scANVI) plus
cell-curator marker evidence. It writes a strict label, a best-guess label, and the reason for each
cluster and each flagged cell. It runs offline. See `SKILL.md` for the full description.

## Run the demo

From the repository root:

```bash
.venv/bin/python skills/annotation-consensus/annotation_consensus.py --demo --output out/demo
```

The demo is 654 cells and three tool prediction files. Its config values are demo values, not
approved defaults. The run refuses to overwrite `out/demo`; use a new `--output` or add `--force`.

Look at `out/demo/report.md`, `consensus.tsv` (one row per cluster) and `cell_flags.tsv`
(only flagged cells).

## Verify

```bash
.venv/bin/python skills/annotation-consensus/annotation_consensus.py verify --demo --output out/demo
```

`verify` hashes the stored outputs, reruns in a fresh subprocess with a scrubbed environment and
compares the sha256 of every hashed output. It exits 2 and names each file that differs.
`run_manifest.json` is not compared. `environment.txt` is hashed, so `verify` passes only on the
environment pinned in `requirements.lock`; a run on other versions stops unless you pass
`--allow-version-drift`. `environment.txt` also records the sha256 of `requirements.lock`; `verify` fails
with "requirements.lock changed since the run" when the lock file on disk differs, unless you pass
`--allow-version-drift`. Exact kNN distance ties are broken by barcode rank, so the input row order does not
change the neighbours.

## Score the demo (optional)

```bash
.venv/bin/python skills/annotation-consensus/annotation_consensus.py run --demo --score --output out/demo_scored
```

The score step is the only step that reads `examples/demo/truth.tsv`. It calls
`benchmarks/scripts/score_v8.py` once per method and takes a few minutes.

## Your own data

Copy `config.example.yaml`, fill in every key, and run:

```bash
.venv/bin/python skills/annotation-consensus/annotation_consensus.py run --input CONFIG.yaml --output OUT
```

Every impactful choice (rule, marker weight, gate, tie-break delta, smoothing k, cluster-gate share,
lead and fallback margin, snap option, marker scope, tool order) is a required key. A missing or
unknown key stops the run. Choose these values yourself and record who approved them in
`consensus.values_status`.

## Snakemake

```bash
snakemake -s skills/annotation-consensus/Snakefile --configfile CONFIG.yaml -c1
```

The config needs `run.output_dir`. Exactly one `--configfile` is allowed and `--config` overlays fail.

## Tests

```bash
.venv/bin/python -m pytest skills/annotation-consensus/tests -q -m "not slow"
```

## Not built

Enrichment of marker genes against known markers, the public-pretrained sensitivity run, and any
change to `src/`.

# Changelog

## 0.3.0 (unreleased)

Correctness and reproducibility release; 57 commits since 0.2.0.

- **Benchmark results** carry `per_type` precision, recall, F1, and support, and store file
  names rather than absolute paths, so a published result is the same on every machine.
  `results/synthetic.json` is superseded by `results/synthetic.v2.json`; the original stays.
- **Removed** the ATAC, CITE-seq, multiome, and spatial benchmark manifests and the five
  `benchmarks/public_*.pending.json` files every run wrote. None pointed at data with truth labels.
- **Docs:** Installation, Troubleshooting, a generated command-line reference, and a
  PBMC 3k tutorial (`examples/pbmc3k_tutorial.py`) that runs from raw counts to labels.
- **`cell-curator benchmark sweep`** writes a risk-coverage curve; `--abstain-unknown`
  stops an `Unknown` label from counting as a wrong call.
- **Benchmark protocol v5** is declared: the Azimuth PBMC level-2 vocabulary, a human
  Gate A, and CellTypist, SingleR, scANVI, popV, and four LLM annotators as comparators.
- **Benchmark protocols v5 and v6 measured.** On AIDA PBMC, v5 (2,000 genes) called
  nothing; v6 (whole transcriptome) calls 11 of 28 clusters at 0.48 descendant accuracy.
  Both "Not comparable" and "Dominated" stop rules fire; SingleR, scANVI, and popV lead.
- **Negative markers** are judged over the measured ones when enough are in the panel;
  one unmeasured negative no longer voids the check.
- **Gene matching** is case-insensitive when unambiguous. The wheel ships a human PBMC marker
  panel, and `references/markers.md` documents the marker file format.
- **Evidence.** Score cells, not cluster means; confidence is a vote share. Cells with no
  marker signal do not vote. Near-ties back off to their shared Cell Ontology ancestor, and
  backoff stays inside the Cell Ontology.
- **Ontology.** Parse the real Cell Ontology as one connected hierarchy.
- **Network.** One guarded door to the network, offline by default.
- **LLM lane.** A gated LLM lane that may only break measured near-ties.
- **Gates.** Critic and write-back gates refuse what they claim to; write-back completion is
  durable and retry-safe.
- **Benchmark.** Neutral prediction export, ontology distance, risk-coverage sweep, and AIDA
  PBMC protocol v3 and v4 results. The v10 result on a row-scrambled input is retracted.
- **Harness.** One skill and one approval gate for Claude Code, Codex, and OpenCode.
- **Input.** Five crashes on real objects fixed; Snakefile paths resolve against the config file.

## 0.2.0 (2026-08-11)

First tagged release.

# Changelog

## 0.3.0 (unreleased)

Correctness and reproducibility release; 57 commits since 0.2.0.

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

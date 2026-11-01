# cell-curator

This package annotates **pre-clustered** AnnData and MuData from signed marker
programs, bottom-up markers, ontology checks, reference mappings, and multimodal
evidence. It records uncertainty and preserves `Unknown` when evidence is
insufficient.

Runs are local and review-gated, and the source object remains unchanged. Writing an
annotated output requires a person to confirm the exact reviewed state. Running cell-curator
requires Python 3.11 or newer and an existing clustering.

Documentation: <https://mdmanurung.github.io/cell-curator/>

## Install

Clone the repository, then install the workflow dependencies:

```bash
git clone https://github.com/mdmanurung/cell-curator.git
cd cell-curator
uv sync --frozen --extra workflow
uv run cell-curator --version
```

Alternatively, install the current wheel from
[GitHub Releases](https://github.com/mdmanurung/cell-curator/releases/latest):

```bash
uv pip install https://github.com/mdmanurung/cell-curator/releases/download/v0.3.0/cell_curator-0.3.0-py3-none-any.whl
```

The wheel installs the workflow, schemas, templates, and agent skill under
`<prefix>/share/cell-curator/`. Available extras are `workflow`, `celltypist`,
`reference-models`, `census`, and `scimilarity`.

### Agent harnesses

The skill and the CLI work the same in Claude Code, Codex, and OpenCode. Install
the CLI once, so every harness can run it:

```bash
uv tool install "cell-curator[workflow] @ git+https://github.com/mdmanurung/cell-curator"
```

Inside this repository nothing else is needed. Codex and OpenCode find the skill
at `.agents/skills/cell-curator` and read `AGENTS.md`, and Claude Code reads
`CLAUDE.md`. To use the skill from any project:

| Harness | Skill | Ask before approvals |
| --- | --- | --- |
| Claude Code | `/plugin marketplace add https://github.com/mdmanurung/cell-curator`, then `/plugin install cell-curator@cell-curator` | the plugin's `hooks/hooks.json` |
| Codex | `mkdir -p ~/.agents/skills && ln -s "$PWD/skills/cell-curator" ~/.agents/skills/` | copy `.codex/rules/cell-curator.rules` to `~/.codex/rules/` |
| OpenCode | the same `~/.agents/skills` link | merge `opencode.json`'s `permission.bash` into your config |

The approval rules make the harness ask a person before `cell-curator
assumptions approve` (Gate A) and `cell-curator write-back apply`, the two commands
that stand in for a human decision.

## Agent-guided run

Provide the input path, cluster key, assay representations, biological context,
and available marker programs. For example:

> Annotate `data/pancreas_clustered.h5ad`. The samples are healthy adult human
> pancreas measured with 10x 3'. Clusters are in `obs["leiden_1"]`, log-normalized
> values in `layers["logcounts"]`, and counts in `layers["counts"]`. Use
> `markers/pancreas.json`. Begin in evidence-only mode and report marker coherence,
> uncertainty, and unresolved clusters before proposing labels.

The agent then follows [SKILL.md](skills/cell-curator/SKILL.md) and the linked analytical
contracts. Interactive and pre-authorized runs produce the same decision ledger.

## Python run

For programmatic use, `CellCurator` exposes the implementation shared by the CLI
and Snakemake workflow:

```python
from cell_curator import CellCurator

cur = CellCurator(
    adata=adata,  # or path="data/pancreas_clustered.h5ad"
    organism="human",
    tissue="pancreas",
    cluster_key="leiden_1",
    markers="markers/pancreas.json",
    run_id="pancreas-v1",
    output_root="/abs/path/results",
    mode="annotation",  # or "evidence-only"
    developmental_stage="adult",
    condition="healthy",
    experimental_context="pancreas atlas v1",
)

cur.inspect()
cur.prepare()
cur.plan()
cur.approve_assumptions(reviewer="your-name")  # Gate A: you approve the plan's assumptions
scores, proposals = cur.evidence()
packet = cur.annotate()
cur.validate_run()
```

`annotate()` resumes completed phases. Write-back requires a person to confirm it. A
tested example is available in
[quickstart_synthetic.ipynb](examples/quickstart_synthetic.ipynb).

## Configure the input

For file-based CLI or Snakemake runs, copy and validate the configuration template:

```bash
cp skills/cell-curator/assets/config.template.yaml config.yaml
uv run cell-curator config validate --config config.yaml
```

Declare each assay and representation explicitly. The table summarizes the main
distinctions:

| Data | Required distinction |
| --- | --- |
| scRNA-seq | Use log-normalized values for ranking and counts for detection. |
| CITE-seq or MuData | Keep RNA and protein as separate evidence lanes. |
| scATAC-seq or multiome | Distinguish peak accessibility from inferred gene activity. |
| Spatial transcriptomics | Treat proximity as context rather than cell identity. |
| Hierarchical labels | Scope every child decision to its frozen parent. |

Set `input.keys.canonical_parent` to the cluster column. Register markers,
ontologies, or labeled references under `knowledge.providers`. Use a new
`run.run_id` when the input or configuration changes. See
[assay-routing.md](skills/cell-curator/references/assay-routing.md) for the full
representation contract.

## Network: offline, institution, or online

Every HTTP request goes through one guard, set by `knowledge.network`:

| `network` | Reaches | Use it for |
| --- | --- | --- |
| `offline` (default) | nothing | local markers and a local Cell Ontology |
| `institution` | localhost and `knowledge.allowed_hosts` | a local or institutional model server, ontology mirror, or table server; nothing leaves the institution |
| `online` | any host | remote tables, live ontology lookups, Census, CellTypist models |

An `allowed_hosts` entry that starts with `.` allows that domain's subdomains.
Redirects are checked too, so an allowed host cannot send a request outside.
CellTypist routes use cached models only and never download. Plan offers a
Census route only in `online` mode. `reference acquire-azimuth` and
`acquire-scimilarity` are downloads a person runs by name, each one GET for a
public, checksum-pinned file. `reference acquire-census` queries Census through tiledb, outside the guard.

LLM tie-breaks follow the same rule: `cell-curator llm propose` reaches only what
`network` allows, so `institution` with a local Ollama or vLLM server keeps every
payload inside the institution. A host agent that reads `llm packet` output sends
it to its own model, so that command refuses unless `network` is `online` or
`evidence.llm.host_agent_local` is set.

**What `network` cannot govern: the agent itself.** In Claude Code or Codex with a cloud model, every file the agent reads goes to that vendor, and the evidence cards, reports, and score tables carry the same numbers as the LLM packet. Nothing leaves the institution only when the harness model is local too (for example OpenCode with Ollama), or when you run the CLI or Python without a cloud agent.

## Execute the workflow

With `config.yaml` prepared, use the guided CLI to move through each review boundary:

```bash
uv run cell-curator config validate --config config.yaml
uv run cell-curator data inspect --config config.yaml
uv run cell-curator run scene --config config.yaml
uv run cell-curator run confirm-context --config config.yaml \
  --input context-confirmation.json --reviewer REVIEWER
uv run cell-curator run plan --config config.yaml
uv run cell-curator assumptions approve --config config.yaml --reviewer REVIEWER
uv run cell-curator run execute --config config.yaml
uv run cell-curator run validate --config config.yaml
```

For additional commands, `cell-curator --help` lists all 22 command groups.

For resumable execution, Snakemake manages dependencies and stale artifacts. The
targets below preserve the evidence, review, and finalization boundaries. Write-back
remains a separate step:

```bash
uv run snakemake --snakefile workflow/Snakefile all_evidence \
  --configfile config.yaml --cores 4
uv run snakemake --snakefile workflow/Snakefile review_packet \
  --configfile config.yaml --cores 4
uv run snakemake --snakefile workflow/Snakefile finalization \
  --configfile config.yaml --cores 4
```

Before executing a new configuration, use `--dry-run`. The Slurm profile is in
[workflow/profiles/slurm](workflow/profiles/slurm/).

## Review and write back

After the run completes and before write-back, review the artifacts under
`run.output_root/run.run_id`. They include evidence receipts, hierarchy decisions,
critic reconciliation, the review packet, and provenance. Every configured critic
pass must record findings, or a severity `NONE` attestation, before write-back
proceeds. Preview the complete mapping, then apply it from your own terminal:

```bash
uv run cell-curator write-back preview --config config.yaml --output annotated.h5ad
uv run cell-curator write-back apply --config config.yaml --output annotated.h5ad
```

`apply` prints a summary and an approval token covering the run, source, labels,
assumptions, report, critic ledgers, diff, and output path, and writes only after you
type the token. A non-interactive run writes only when its configuration set
`guidance.preauthorized: true` before scene setting. See
[artifact-schemas.md](skills/cell-curator/references/artifact-schemas.md) for the
receipt.

## Analytical constraints

- Cell type and cell state are separate annotation axes.
- Missing panel features are unmeasured, not negative evidence.
- ATAC peaks provide direct accessibility evidence; gene activity is inferred.
- A RAPIDS route requires its declared GPU and software contract. CPU output is
  not substituted for a required GPU result.
- `Unknown` is retained when coverage, coherence, margin, or capability is
  inadequate.
- Supplied Cell Ontology identifiers require local semantic validation. Labels may
  remain ontology-free when no ontology is available.
- Comparative claims require measured benchmark results.

## Known limits

- **No comparative result yet.** The only published benchmark scores cell-curator alone
  on one AIDA PBMC slice, and it fires its own "not comparable" stop rule. See
  [benchmarks/README.md](benchmarks/README.md).
- **Confidence is a vote share, not a calibrated probability.**
- **Human-first.** Marker genes match case-insensitively, but there is no ortholog
  mapping, and only a human PBMC marker set ships.
- **CPU scale.** The CPU bottom-up marker lane densifies each cluster against the rest,
  and stability refinement runs networkx Louvain. Neither has been measured at atlas
  scale. For hundreds of thousands of cells, use the RAPIDS route.
- **scANVI and scArches need a CUDA GPU.** The Azimuth-compatible lane is a kNN on
  harmonized genes, not Azimuth itself.

## Documentation

| Topic | Document |
| --- | --- |
| Agent workflow | [SKILL.md](skills/cell-curator/SKILL.md) |
| Inputs and matrices | [input-contract.md](skills/cell-curator/references/input-contract.md) |
| Marker file format | [markers.md](skills/cell-curator/references/markers.md) |
| Assays and providers | [assay-routing.md](skills/cell-curator/references/assay-routing.md) |
| Parent-local refinement | [adaptive-subclustering.md](skills/cell-curator/references/adaptive-subclustering.md) |
| Multimodal evidence | [multimodal-evidence.md](skills/cell-curator/references/multimodal-evidence.md) |
| Decisions and critics | [decision-gates.md](skills/cell-curator/references/decision-gates.md) |
| Artifact schemas | [artifact-schemas.md](skills/cell-curator/references/artifact-schemas.md) |
| Compute and resume | [compute-and-snakemake.md](skills/cell-curator/references/compute-and-snakemake.md) |
| Validation | [validation.md](skills/cell-curator/references/validation.md) |
| Capability mapping | [quadbio-parity.md](skills/cell-curator/references/quadbio-parity.md) |

# QuadBio integration and deliberate extensions

Source: <https://github.com/quadbio/cell-type-annotation>

Pinned source: commit `b4381198b85b7fd422f6640c1d7a04c2978ace1d`, package
`cta` version `0.1.0`. The locally installed skill and routed references were
byte-identical to this commit when this package was authored.

## Reused unchanged

| QuadBio behavior | Use here |
|---|---|
| The agent is the annotator | Tools generate hypotheses and evidence; no classifier writes the verdict. |
| Artifact-centric, resumable work | Context, assumptions, plan, run state, raw evidence, cards, critics, and decisions are durable artifacts. |
| Gate A | A living assumptions ledger requires explicit sign-off before labels are reviewable. |
| Gate B | `Unknown` is a valid, non-failure outcome. |
| Evidence hierarchy | Positive and negative programs, panel-aware interpretation, rejected alternatives, and explicit confidence are required. |
| Type/state separation | Cell type and state are separate annotation axes. |
| Bounded subclustering | Only heterogeneous parents are split, with a limited local pass. |
| Independent critics | Critics inspect raw tables before reconciliation. |
| Hierarchy and ontology | Labels must nest semantically and CL IDs require Knowledgebase validation. |
| Explicit write-back approval | Evidence completion never authorizes mutation of source data. |

## Extended for assay-aware production workflows

| Extension | Reason |
|---|---|
| AnnData and MuData assay registry | RNA-only, protein-only, ATAC-only, gene-activity, CITE-seq, multiome, spatial, and future assays declare their own capabilities. |
| Representation-aware detection | Signed values provide means/effects; detection uses a separate nonnegative source or is explicitly unavailable. |
| Frozen canonical parents and nested leaf IDs | Global reclustering would discard the stability-selected parent partition. |
| Stored-resolution families plus one local pass | Existing partitions are screened before bounded recomputation; Leiden nesting is verified rather than assumed. |
| Two native binary comparators per candidate | Candidate-versus-remainder and candidate-versus-sibling answer different biological questions. |
| Configurable sample/QC/doublet gates | Discrete-looking children can be technical or multiplet artifacts; absent metadata is not fabricated. |
| Feature-space-aware corroboration | PanglaoDB applies to gene-like spaces; peak/motif/protein sources route separately. |
| Immutable per-cell mappings and leaf totals | Nested identities must conserve parent membership and support downstream invariance checks. |
| Snakemake-authoritative HPC phases | Heavy multimodal evidence needs reliable dependency, staleness, resume, and resource boundaries. |

## Intentionally overridden

| QuadBio behavior | Override | Reason |
|---|---|---|
| Marker tooling may route to a CPU implementation when GPU execution is unavailable | Native `rapids-singlecell` marker calculations require one real validated GPU and fail closed. No CPU marker fallback exists or may be scheduled. | CPU/GPU method substitution changes the frozen evidence contract and hides environment failures. |
| General GPU routing can choose among available devices | Allowed GPU model, compute capability, CUDA/RAPIDS distribution, and visible-device count are explicit configuration contracts. | Reproducible HPC evidence requires validated runtime provenance. |
| Transcriptomic annotation is the main substrate | A capability registry routes arbitrary unimodal or multimodal evidence; unsupported claims remain blocked. | Studies may contain RNA only, ATAC only, protein only, or different combinations. |
| A flat subclustering artifact can be sufficient | Every child remains nested under a frozen canonical parent, with a one-row-per-cell immutable mapping. | Downstream totals and parent membership must remain invariant. |

## Knowledgebase dependency

The QuadBio Knowledgebase MCP is required before proposing, assigning, or
semantically validating biological labels or CL IDs. If it is unavailable, use
`evidence-only` mode: leave labels and CL IDs blank, record the blocker, and
prohibit final annotation and write-back.

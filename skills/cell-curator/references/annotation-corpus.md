# Annotation corpus

Use a corpus to keep label definitions, their sources, and their revisions next
to the measured annotation evidence. The corpus has no preferred atlas. Allen,
CAP, other published references, and reviewed local studies remain separate
sources. Importing a label does not approve its definition or establish that two
sources mean the same thing.

## Inspect the example

The bundled example contains short normalized documentation excerpts from two
sources. All four labels are `imported`, so none is callable:

```bash
cell-curator corpus validate --input skills/cell-curator/assets/corpus/example.json
```

It is a format example, not a complete immune-cell reference or a validated
marker panel. Source reuse terms are not established. The snapshot points to
the original pages and hashes the saved normalized excerpts. No reference model
is bundled, downloaded, or selected automatically.

## Label records

A snapshot contains `schema_version: 3`, `corpus_id`, `revision`, `sources`,
`labels`, and optional `mappings`, `parent_revision`, and `review`.

Each source has a stable `source_id`, name, URL, version, content SHA-256,
retrieval date, reuse information, and optional `depends_on` source IDs. Source
dependencies record shared provenance; source counts are not independent votes.

A label records its source and biological definition. For example, this is a
synthetic proposed child record, not a biological marker recommendation:

```json
{
  "label_id": "study:subtype-a",
  "name": "Subtype A",
  "description": "Study-specific subdivision supported by measured evidence.",
  "source_id": "my-study",
  "source_label": "original subtype A",
  "source_level": "fine",
  "level": "L3",
  "parent_id": "known:parent-a",
  "axis": "cell_type",
  "species": ["Homo sapiens"],
  "tissues": ["blood"],
  "positive": ["FEATURE_A"],
  "negative": [],
  "feature_space": "gene",
  "measurement": "gene_expression",
  "study_id": "my-study",
  "status": "proposed"
}
```

The referenced source and parent must already exist or be included in the same
proposal. A reviewer must supply a rationale before the record becomes
`reviewed`. Defining a label requires more than replacing the example feature.
Reviewed identity levels must agree with parent depth: roots are L1 and each
child advances one level. State labels can use L4 beneath an identity at any depth.

Other fields include `full_name`, `synonyms`, `alternatives`, `limitations`,
`assertions`, and `experience`. Assertions distinguish canonical markers,
source-reported marker evidence, and markers observed during a recorded run.
Source thresholds retain their measurement and units; they are source assertions,
not silently adopted study decision thresholds.

Cell Ontology mappings record both `cl_id` and `ontology_relation`. Only an exact
mapping validated against a supplied local ontology is exported as an exact CL
ID. Broader and unresolved mappings stay on the rich record. Export without an
ontology keeps marker labels usable and omits exact CL claims.

When exporting CL IDs, the display label must also be a name or synonym accepted
by the existing local ontology checks. A full ontology name in `full_name` does
not make an unrelated display abbreviation valid.

Mappings between corpus records require source and target label IDs, a relation,
and a rationale. Relations are `exact`, `broader`, `narrower`, `related`,
`conflicting`, or `unresolved`. Exact mappings require a named reviewer. Similar
names do not establish equivalence; ambiguous exported names or synonyms fail.

## Import saved records

The import commands read bounded local JSON. They do not scrape HTML or query
CAP's GraphQL API. Normalize a saved source into a list, or this wrapper:

```json
{
  "context": {"species": "Homo sapiens", "tissue": "blood"},
  "labels": [
    {
      "id": "source-type-a",
      "cell_label": "Source type A",
      "rationale": "The source's biological description.",
      "canonical_marker_genes": ["FEATURE_A"],
      "marker_gene_evidence": ["FEATURE_B"],
      "level": "source-fine-level"
    }
  ]
}
```

CAP-style fields include `cell_label`, `cell_fullname`, `rationale`,
`canonical_marker_genes`, `marker_gene_evidence`, and `ontology_term_id`.
Allen-style normalized records use `name`, `description`, `markers`, `level`,
and `parent_id`. Generic normalized records can use `positive`, `negative`,
`source_label`, `source_level`, and an explicit `canonical_level`.

The importer preserves source level and parent identifiers separately. It does
not turn a source's L3 into the study's L3. Missing canonical levels use an
unreviewed L1 placeholder. The import receipt reports this omission. Reviewers
must map the hierarchy before those records become callable.

```bash
cell-curator corpus import --input corpus-v1.json --records saved-records.json \
  --output corpus-imported.json --source-kind cap --source-id source-version-a \
  --name "Source A" --url https://example.org/source-a \
  --source-version VERSION --retrieved-at YYYY-MM-DD
```

Use `--source-kind allen` or `generic` for the other normalized formats. The
default limits are 10 MiB and 10,000 records. The Python import API can set
smaller limits. Input values are preserved in the adjacent `.source.json` file.
The `.import-receipt.json` records failed rows, omitted fields, digests, and the
partial-import status. An imported row stays unreviewed even if its source JSON
claims an approved status. Negative markers are never fabricated from positives.

### scGate_models

`--source-kind scgate` reads a local clone of `carmonalab/scGate_models`. The
clone must sit exactly at `--scgate-tag` (default `v0.16`). Nothing is downloaded
and nothing is shipped, because the repository has no licence file. `--records` is
the clone folder.

```bash
cell-curator corpus import --input corpus-v1.json --records /path/to/scGate_models \
  --output corpus-imported.json --source-kind scgate --scgate-collection PBMC \
  --species human --tissue blood --source-id scgate-pbmc \
  --name "scGate_models human PBMC" --url https://github.com/carmonalab/scGate_models \
  --source-version v0.16 --retrieved-at YYYY-MM-DD
```

Each model file becomes one label. Its positive genes are the signature of the
model's last level. A gene written `GENE-` is a negative gene. The gating path,
the populations the model excludes (`scgate_excludes`) and the model file hash
stay in the `.source.json` copy for the reviewer. Each excluded population that has
a signature in `master_table.tsv` becomes an `exclusions` entry on the label. At
`collect-level`, a cluster fails the label (`excluded:<name>`) when that population's
own genes pass the state rule (`min_state_positive_fraction`, `min_state_fraction`).
An exclusion the panel does not measure is flagged `exclusion_unmeasured:<name>`.
`CellOntology_dictionary.tsv` gives a first guess for `cl_id`; the ontology check
still runs at review. Import one collection per source id: two collections that
define the same label name need an explicit exact mapping.

## Review and release corpus knowledge

Create `changes.json` with `labels`, `sources`, and/or `mappings`. Supply full
replacement records for corrected IDs. Keep unrelated records out of the change.

```bash
cell-curator corpus propose --input corpus-v1.json --changes changes.json \
  --output proposal.json --reason "Source-backed correction" --proposer ANALYST
```

A person reviews the proposal and records their decision:

```bash
cell-curator corpus review --input corpus-v1.json --proposal proposal.json \
  --output corpus-v2.json --decision accept --revision v2 \
  --reviewer REVIEWER --rationale "Reviewed provenance, scope, and criteria"
cell-curator corpus diff --input corpus-v1.json --other corpus-v2.json
```

Use `reject` to preserve the old vocabulary in a new review revision, or `dispute`
to retain the proposed records outside the callable vocabulary. A disputed
replacement suspends the affected definition until resolved. The original
snapshot remains available.

Review checks the proposal's base digest. It refuses changed bases and input
overwrites. A new revision gets an immutable receipt with reviewer, rationale,
proposal hash, and changed fields. Existing runs retain their pinned corpus.
Corpus review does not approve a run's assumptions and does not approve
write-back. Agents must not claim a human review on their own behalf.

Export reviewed marker programs for a specified context:

```bash
cell-curator corpus export --input corpus-v2.json --context context.json \
  --measurement gene_expression --output reviewed-markers.json
```

`context.json` can name `species`, `tissue`, `assay`, `stage`, `condition`, and
`study_id`. Missing context cannot satisfy a constrained record. Human and mouse
species names accept their common names, Latin names, and taxonomy IDs. Pass
`--ontology ontology.obo` to validate and export exact CL mappings. Transcript
isoform rules currently have no executable marker backend and are refused.

## Configure an annotation

Extend the existing strict run configuration. Get the corpus file digest with
`sha256sum corpus-v2.json`; use that actual digest and snapshot revision below:

```yaml
knowledge:
  network: offline
  providers:
    - kind: corpus
      path: corpus-v2.json
      sha256: REPLACE_WITH_FILE_SHA256
      version: v2
      context:
        species: Homo sapiens
        tissue: blood
guidance:
  hierarchy_levels: [L1, L2, L3]
annotation:
  enabled: true
  assay: transcriptome
  corpus_assay: rna
  representation: logcounts
  study_id: my-study
  seed: 0
  min_classifier_agreement: 0.5
  min_state_fraction: 0.5
  classifiers:
    - name: matched-celltypist
      method: celltypist
      model: models/matched.pkl
      counts_layer: counts
      required: true
      label_map:
        reference_alias_a: reviewed-corpus-label-a
    - name: matched-popv
      method: popv
      model: models/popv
      reference: references/matched.h5ad
      reference_label: cell_type
      counts_layer: counts
      reference_counts_layer: counts
      methods: [CELLTYPIST, Support_Vector, XGboost]
      required: true
      label_map:
        reference_alias_a: reviewed-corpus-label-a
```

Model and reference paths must point to locally provisioned resources. Matching
requires species, tissue, stage, condition, and assay assessment; a filename is
not evidence of a suitable reference. Declare the rest of the biological and
input context in the normal run configuration.

Corpus assay scopes use a stable context string, not the arbitrary input assay
key. By default, the selected assay's declared kind supplies `rna`, `protein`,
or another configured kind. `annotation.corpus_assay` can instead declare a
specific platform such as `10x-3prime-v3`. Source record `assays` values must match
that declared string. The example's input key is `transcriptome`, while its
corpus assay scope is `rna`.

Provider context cannot contradict known dataset biology or supply missing
species, tissue, stage, or condition. Establish those facts through the normal
context review. Provider selection must not turn an unknown dataset species
into a confidently selected human reference.

The corpus DAG currently uses the explicit CPU marker backend. It collects
Wilcoxon and logistic-regression markers within each accepted parent, evaluates
signed programs, and retains contradictory evidence. Native GPU lanes remain
in the existing workflow; this DAG never substitutes CPU output for a required
GPU route.

This DAG scores one explicitly selected assay. Its corpus programs must use
that assay's measurement and feature space. A protein rule cannot run on an RNA
matrix because another modality happens to be present. Unsupported measurements
fail explicitly; the existing multimodal workflow remains available separately.

The popV adapter currently supports CPU `fast` inference from a local model and
reference. It does not train models or download resources. The constituent
method list is recorded. Shared CellTypist evidence is accounted for rather
than counted as an independent confirmation. Preserve classifier provenance
when deciding which model evidence is independent.

popV 0.6 requires Python 3.11 or 3.12. Its TensorFlow/numpy dependencies conflict
with the current workflow/Slurm extra. Keep the environments separate:

```bash
uv sync --frozen --extra workflow
UV_PROJECT_ENVIRONMENT=.venv-popv uv sync --frozen --extra popv --python 3.11
```

Set `python_executable: .venv-popv/bin/python` on the popV classifier. Paths
resolve relative to the configuration file. The main Snakemake environment
schedules the worker in that interpreter. The worker runs the same configured
classifier stage and publishes its receipts; it never falls back to a different
environment. The interpreter file is included among the planned resources.

Classifier labels that exactly match a reviewed corpus name are recognized.
Other reference names require an explicit `label_map`. Unmapped names remain
`Unknown`. Mapping a name does not override measured positive, negative, or
hierarchy checks. Required classifier failures block completion; optional
unavailable classifiers remain visible.

## Execute and inspect

Complete the normal input inspection, scene confirmation, and plan first. A
person must approve material run assumptions under Gate A. Neither Snakemake
nor the corpus commands grant that approval.

```bash
cell-curator config validate --config config.yaml
cell-curator run scene --config config.yaml
cell-curator run confirm-context --config config.yaml \
  --input context-confirmation.json --reviewer REVIEWER
cell-curator run plan --config config.yaml
cell-curator assumptions list --config config.yaml
```

After the person approves those assumptions, inspect the DAG before running it:

```bash
snakemake --snakefile workflow/annotation.smk --configfile config.yaml --cores 4 --dry-run
snakemake --snakefile workflow/annotation.smk --configfile config.yaml --cores 4
snakemake --snakefile workflow/annotation.smk --configfile config.yaml --cores 4 annotation_review
```

Snakemake schedules classifiers, per-level marker evidence, per-level decisions,
and final assembly. The next level depends on the previous level's accepted
parent decisions. Config overlays are refused. The plan pins corpus, reference,
classifier resources, and rules; changed resources require a new run.

`cell-curator run execute --config config.yaml` and the Python `annotate()` API
also invoke this compiled DAG when `annotation.enabled` is true. The direct
Snakemake command lets you choose cores and inspect a dry-run. The separate
`annotation_review` target creates the normal report, cards, and critic packet;
it does not record critic findings, approve assumptions, or write annotations.

Inspect `annotation/report.html`, each level's `marker_evidence.svg`, marker and
score tables, and `annotation/cell_labels.tsv` under the run directory. Calls
retain their alternatives and failed gates. An unresolved child keeps its
accepted parent. The per-cell table declares `cluster_propagated_to_cells`;
cluster-derived assignments are not presented as independent cell-level calls.

L1 is lineage, L2 an intermediate identity, and L3 a finer supported identity.
State records use `axis: cell_state`, a parent identity, and an explicit
measurement. Accepted states remain separate in `cell_state`. L4 displays the
accepted identity plus states; it does not replace L3 or create an ontology
term. Inspect this study display in the annotation output and write-back preview.

Finish the existing report, hierarchy, and critic checks before write-back.
Every configured critic needs findings or a `NONE` attestation. Preview the
complete diff, then ask the person to run `write-back apply` and enter its token
in their own terminal. See [decision gates](decision-gates.md) and
[artifact schemas](artifact-schemas.md).

## Limits of learning

New studies can supply proposed definitions, corrections, source observations,
and contradictions. They do not retrain classifiers or promote their own
predictions automatically. Validate new corpus revisions against independent
held-out evidence before treating them as scientifically qualified. Repeated
calls from the same source or donor do not establish new independent support.
Deterministic replay checks reproducibility of a decision, not its biological
correctness.

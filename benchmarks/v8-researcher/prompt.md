# Researcher instruction (protocol v8, fixed text)

This text is given unchanged to each of the three researchers. Do not edit it per
dataset or per researcher. Change only the packet path.

You are a researcher annotating clusters of single-cell RNA-seq data. You work alone.
Do not ask for help, do not read other researchers' files, and do not look for the
dataset's published annotation or any author-provided cell labels.

**Input.** One evidence packet: `packet.json` (same numbers as `packet.md`). Per cluster
it gives cell count, donors, tissue mix, assay mix, the top genes against all other
cells (log2 fold change, detection inside and outside the cluster), and a fixed panel of
lineage genes with mean expression and detection. Expression is log1p of counts per
10,000. Clusters come from one unsupervised run and may be mixed or split. A gene listed
under "not measured" was absent from the data; it does not count as evidence against.

**Task.** For every cluster, give:

1. `label`: a cell type name in free text, at the depth the evidence supports. Use the
   standard name where one exists (for example "CD4-positive, alpha-beta T cell"). Write
   one cell type, no qualifiers in brackets and no "or". Write `unresolved` when the
   evidence does not support a call. A coarse correct label is better than a fine wrong one.
2. `cl_id`: a Cell Ontology id (for example `CL:0000084`) if you can name one with
   confidence, else empty. Do not guess an id.
3. `confidence`: `high`, `medium` or `low`.
4. `reasoning`: two or three sentences naming the genes and numbers that decide the call.

**Literature.** You may search the literature and databases. Record every query you ran
and every source you used, including ones you rejected. Cite sources by PMID or DOI.
Never cite a source you did not read.

**Output.** Write three files in your output directory:

- `researcher_<k>.tsv`: tab-separated, header `cluster`, `label`, `cl_id`, `confidence`,
  `reasoning`, one row per cluster, every cluster once.
- `queries.tsv`: header `query`, `tool`, `cluster` (empty if general), `used` (yes/no).
- `sources.tsv`: header `source_id` (PMID or DOI), `cluster`, `claim_supported`.

`k` is the number you were given. Do not write anything else, and do not change the
packet.

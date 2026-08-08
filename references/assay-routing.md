# Assay routing and capability matrix

Use this reference while writing `config.yaml` and when an evidence gate is not
obviously applicable. Register what the study actually measures; never create a
placeholder assay merely to satisfy a template.

| Assay / feature space | Direct evidence | Typical corroboration | Important limitation |
|---|---|---|---|
| RNA / genes | identity, subtype, state, programs | PanglaoDB ORA, ontology markers, reference mapping | dropout and panel coverage affect negative evidence |
| Protein / antibody features | surface/intracellular identity and state | cytometry gates, curated antibody panels | antibody specificity, background, and panel omissions |
| ATAC / peaks | differential accessibility, regulatory state, discrete chromatin populations | peak enrichment, peak-to-gene links, motifs | raw peaks do not directly support gene-marker ontology claims |
| Gene activity / genes | gene-level regulatory accessibility | gene-marker sets, PanglaoDB as qualified corroboration | activity is inferred accessibility, not RNA expression |
| Spatial / genes or proteins | identity plus spatial context | neighborhood and anatomical priors | proximity is not identity; panel limitations apply |
| Other | only declared capabilities | project-specific validated backend | unsupported claims remain blocked |

## Routing rules

1. Declare each assay's `kind`, `feature_space`, feature identifier, graph inputs,
   evidence representations, and capabilities.
2. Mark a representation `ontology_compatible: true` only when its feature IDs
   can be validly compared with the marker/ontology source being used.
3. Use PanglaoDB ORA only for `gene` or `gene_activity` rankings with a measured
   background. Qualify gene-activity conclusions as accessibility evidence.
4. For peak ATAC, prefer a registered differential-accessibility backend. Add
   peak-to-gene, motif, or reference-label evidence only when its method and
   receipt are declared. RAPIDS ranking is allowed only when scientifically
   appropriate for the configured numerical representation.
5. For protein, define detection from a nonnegative count/logcount source. A
   signed corrected layer supplies means/effects, not `>0` detection.
6. If no ontology-compatible lane or validated reference exists, remain
   label-free: audit parents, discover stable children, and report evidence, but
   block biological labels and CL IDs.

## Common configurations

- RNA-only h5ad: one `rna/gene` assay and one or more RNA evidence lanes.
- ATAC-only h5ad: one `atac/peak` assay, optionally with a valid gene-activity
  representation and feature mapping.
- CITE-seq h5mu: separate RNA and protein assays; corrected and nonnegative
  protein representations may share one detection source.
- Multiome h5mu: RNA plus ATAC, optionally gene activity; choose the graph used
  by the frozen canonical workflow.

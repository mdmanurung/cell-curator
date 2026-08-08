# Input and freeze contract

## Required declarations

Declare input kind/path, canonical parent and scope keys, stored resolutions,
versioned output root, and operating mode. Register every actual assay with:

- logical assay ID and, for MuData, modality name;
- assay kind, feature space, and feature ID mapping;
- frozen graph representations;
- evidence representations, matrix sources, ranking backends, and capabilities;
- optional nonnegative detection source and whether detection is required;
- optional donor, capture, doublet, sex, batch, guide, and technical-metric keys.

An h5ad declares one root assay. An h5mu may declare any number of modalities.
Names such as `rna`, `adt`, or `atac` are configuration, not logic.

## Freeze receipt

Hash input, configuration, scripts/rules, and external guide/reference snapshots.
Record object shape, ordered and sorted barcode hashes, parent counts, assay
dimensions, representations, layers, stored resolutions, packages, and compute
environment. Use the ordered hash for exact joins and the sorted hash for
membership invariance.

## Checks

1. Barcodes are unique and aligned across declared MuData modalities.
2. Canonical parent/scope keys are complete and declared parent totals match.
3. Every representation and detection source exists and has cell-feature shape.
4. A required detection source is nonnegative; signed sources block
   detection-dependent acceptance.
5. At least one ontology-compatible lane or validated reference exists before
   biological labels become reviewable.
6. Optional metadata is either present or explicitly unavailable. Configure no
   donor/capture gate when the study lacks that axis.
7. Guides remain hypotheses and source/completed packets remain read-only.

Stop on ambiguity. Record missing optional capabilities as limitations.

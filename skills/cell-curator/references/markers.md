# Marker files

A `local_markers` provider reads one JSON, TSV, or CSV file. Each entry is one
marker program: a label with the genes expected on and the genes expected off.

## JSON

A mapping from label to program:

```json
{
  "B cell": {
    "cl_id": "CL:0000236",
    "positive": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "negative": ["CD3E", "NKG7", "LYZ", "CD14"]
  }
}
```

A list of objects with a `label` field also works.

## TSV or CSV

One row per program. Gene lists split on `;`, `,`, or `|`.

```text
label	positive	negative	cl_id
B cell	CD79A;CD79B;MS4A1;CD19	CD3E;NKG7;LYZ;CD14	CL:0000236
```

## Fields

| Field | Required | Meaning |
| --- | --- | --- |
| `label` | yes | Cell type name. Use the Cell Ontology label when `cl_id` is set. |
| `positive` (or `pos`) | yes | Genes expected to be detected. |
| `negative` (or `neg`) | no | Genes expected to be absent. A gene may not be both. |
| `cl_id` | no | Cell Ontology ID, such as `CL:0000236`. Needed for ontology backoff. |
| `parent` | no | Label of the parent program, for hierarchical runs. |
| `parent_cl_id` | no | Cell Ontology ID of the parent. |
| `axis` | no | `cell_type` (default) or `cell_state`. A state program is never an identity call. |
| `feature_space` | no | `gene` (default), or the assay's feature space for non-RNA programs. |

## How genes are matched

Genes match `var_names`, or the symbol column the input contract names. Matching is
exact first, then case-insensitive, so a human `CD3E` panel finds mouse `Cd3e`. That
is not ortholog mapping: genes whose names differ between species are not found.

Both signs need `evidence.min_positive_coverage` of their genes measured. A program
with fewer measured positives is left `Unknown`. Negatives are judged over the
measured ones; with fewer than that share measured, the negative check is
unavailable and the program is left `Unknown` too. On a highly-variable-gene panel,
check both lists against `adata.var_names` before running.

## Writing a panel

- Use 5 to 10 positive genes and at least 4 negative genes. Negatives are what stop a
  confident wrong call between siblings.
- Pick negatives from the closest siblings, not from distant lineages.
- Record where the panel came from. A panel written by the person scoring the result
  inflates that result.

## Shipped panel

`assets/markers/human-pbmc.json` is a 13-type human PBMC vocabulary with Cell Ontology
IDs. It was written from standard PBMC immunology for the AIDA benchmark and carries
3 negatives per type. It has not been checked against an independent panel. Treat it as
a starting point.

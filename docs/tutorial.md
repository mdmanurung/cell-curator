# Tutorial: PBMC 3k

This tutorial annotates the 10x Genomics PBMC 3k dataset with the human PBMC marker
panel that ships with the package. It needs scanpy and a network connection for the
two downloads.

```bash
uv run python examples/pbmc3k_tutorial.py
```

The script keeps the clusters from the scanpy tutorial but gives cell-curator only
their numeric codes. At the end it prints each cluster's cell-curator label beside the
scanpy name, so you can see which calls agree and which clusters stay `Unknown`.

`Unknown` is a result, not a failure. The `rationale` column of
`labels/L1.tsv` in the run directory says which gate stopped each call. A common one on
this dataset is "expected-negative markers were violated", where a cluster expresses a
gene the panel says it should not.

If every label lacks a `cl_id`, the ontology download failed and the run continued
without it. See [Troubleshooting](troubleshooting.md).

```{literalinclude} ../examples/pbmc3k_tutorial.py
:language: python
```

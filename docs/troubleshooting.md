# Troubleshooting

Each entry is an error cell-curator prints, what causes it, and the fix.

## Input

`declared layer is absent: logcounts`
: The config names a layer the object does not have. Add it, for example
  `adata.layers["logcounts"] = adata.X.copy()` after `sc.pp.normalize_total` and
  `sc.pp.log1p`, or point the representation at the layer you have.

`feature ID column is absent: gene_symbol`
: The config expects gene symbols in a `.var` column. Add
  `adata.var["gene_symbol"] = adata.var_names`, or name the column you have.

`features use Ensembl identifiers but no gene-symbol column or mapping was supplied`
: `var_names` are Ensembl IDs. Supply a symbol column or `gene_mapping_path`.
  Marker panels are written in symbols.

`feature IDs are not unique: ...`
: Two features share a symbol. Run `adata.var_names_make_unique()` or pick a unique column.

`declared canonical parent total (N) differs from input (M)`
: `expected_parent_count` in the config does not match the number of clusters.
  Set it to the real count.

## Markers

`marker source does not exist: ...`
: The marker path is wrong. Paths in a config resolve against the config file.

`marker program 'X' has contradictory signed features: [...]`
: A gene is listed as both positive and negative for one type. Remove it from one list.

Every cluster is `Unknown` with "too little of the signed marker program was measured"
: Fewer than `evidence.min_positive_coverage` of a panel's genes are in the data. This
  is common on a highly-variable-gene subset. Check which panel genes are in
  `adata.var_names`, or run on all genes.

## Ontology and network

`network is offline (knowledge.network); refusing 'host'`
: A provider needs the network and the run is offline, which is the default. Use a
  local ontology file, or set `knowledge.network: online` for the first fetch; later runs
  read the cache.

`'host' is outside the institution (knowledge.allowed_hosts)`
: The run is in `institution` mode and the host is not listed. Add it to
  `knowledge.allowed_hosts`.

Labels have no `cl_id`, and `knowledge/ontology_provider_status.json` says
`CERTIFICATE_VERIFY_FAILED`
: Python could not verify the ontology server's certificate, so the run continued
  without the ontology. uv-managed Python uses its own OpenSSL, which may not find the
  system certificates. Point it at them, for example
  `export SSL_CERT_FILE=/etc/pki/tls/certs/ca-bundle.crt` on RHEL or
  `/etc/ssl/certs/ca-certificates.crt` on Debian and Ubuntu, then use a new `run_id`.

`CL IDs cannot be written without a validated local ontology`
: Labels carry Cell Ontology IDs but no ontology provider is configured. Add one, or
  drop the IDs.

`ontology term is obsolete: CL:...`
: The marker panel uses a retired Cell Ontology term. Replace it with the term the
  message names.

## Gates

`Gate A blocks biological labels until N assumption(s) are ...`
: The plan's assumptions need a person's approval. Read them, then run
  `cell-curator assumptions approve --config CONFIG --reviewer NAME`.

`Gate A blocks biological labels while the run is evidence-only`
: The config sets `run.mode: evidence-only`, which never assigns labels. Set
  `run.mode: annotation` for a new run.

`refusing to overwrite write-back output: ...`
: Write-back never replaces a file. Choose a new `--output` path.

`refusing to replace non-identical immutable artifact: ...`
: A run directory is write-once. Use a new `run.run_id` instead of editing a finished run.

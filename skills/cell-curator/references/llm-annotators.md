# LLM annotators: what leaves the machine, and what checks the answer

Three projects use a language model to name clusters. Every entry below
was read from the source at the pinned commit, and the path cited is relative to
that repository. None of them checks a model's answer against the measured
evidence.

| | [CyteType](https://github.com/NygenAnalytics/CyteType) `42c4c0e` | [mLLMCelltype](https://github.com/cafferychen777/mLLMCelltype) `26bda56` | [CellTypeAI](https://github.com/rhdaw/CellTypeAI) `5507d62` | `cell-curator` |
| --- | --- | --- | --- | --- |
| What is sent | expression matrix, raw counts, the whole `obs` table with barcodes and any donor columns, embedding coordinates (`cytetype/main.py:301-309`, `core/artifacts.py:515`) | marker gene symbols per cluster, species, tissue (`python/mllmcelltype/prompts.py:29`) | top Wilcoxon gene names per cluster, species, tissue (`src/celltypeai/h5ad_to_json.py:31-32`) | near-tied candidates, top markers, vote shares, pseudobulk means for clusters of at least 20 cells; never per cell |
| Where it goes | Nygen's closed server, `https://cytetype.nygen.io` (`config.py:16`); the user's LLM key goes too (`api/schemas.py:35`, `core/payload.py:90-91`) | a cloud provider; a `base_url` override exists, but an API key is still required (`annotate.py:54-63`) | Ollama at a hard-coded `http://localhost:11434` (`prompt_submit.py:104`) | `knowledge.network`: `offline` (default), `institution` (localhost and allowed hosts), or `online` |
| Local only | no; local Ollama must be exposed to Nygen through ngrok (`docs/ollama.md:1`) | not as a supported mode | yes, only | yes, and the default |
| Abstains | not in the client; confidence is server-side and "model-assessed" (`docs/results.md:24`) | no; without consensus it returns the last majority (`consensus.py:2042`) | no; the mode over repeats (`main.py:106-108`) | yes; Unknown unless the gates pass, or a shared CL ancestor holds the votes |
| Checks the answer | not in the client; the server's code is closed | no; cited genes are stored as text (`utils.py:794`) | no; it checks only the reply's JSON shape (`prompt_submit.py:184-187`) | yes; the gates above |
| Cell Ontology | output fields, mapped on the server (`core/results.py:67-76`) | none | none | the label and its CL id are validated together against a pinned snapshot |

mLLMCelltype's consensus proportion and entropy are computed in code first
(`consensus.py:1143`). In harder cases a model's own figures are used instead,
but only if they pass a plausibility check (`:1194-1199`).

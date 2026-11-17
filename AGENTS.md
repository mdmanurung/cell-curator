# cell-curator: instructions for coding agents

Annotate pre-clustered AnnData/MuData with evidence-gated, ontology-checked
labels. The workflow is in `skills/cell-curator/SKILL.md`. Read it before running
anything.

## Commands

- `cell-curator run scene|plan|execute|validate --config CONFIG`: the guided run.
- `cell-curator annotation collect-level|import-consensus --config CONFIG --level L1`, then `annotation finalize --config CONFIG`: evidence, the consensus labels, the durable tables.
- `cell-curator progress status --config CONFIG`: where the run is, and what is next.
- `cell-curator write-back preview|apply --config CONFIG --output OUT.h5ad`.
- Tests: `.venv/bin/python -m pytest -q`, plus `ruff check .` and `mypy`.

## Gates you must not route around

1. **Gate A.** Assumptions are approved by a person. Do not run
   `cell-curator assumptions approve` on your own initiative. Ask the user.
2. **Critics.** Every configured critic pass needs a finding or an attestation
   before write-back.
3. **Corpus review.** `cell-curator corpus review` makes records callable. Run it
   only after the user says OK in the conversation for that proposal. Pass the
   user's name as `--reviewer`. Never review on your own initiative.
   Exception (user, 2026-11-10): `annotation.development_candidates` lets unreviewed candidate records run on development splits only; such runs can never be written back.
4. **Write-back.** `cell-curator write-back apply` asks a person to type the
   approval token in a terminal. Stop and ask the user to run it. Never pipe the
   token in.

## Labels come from the consensus run, not from you

An `annotation-consensus` run (`skills/annotation-consensus`) chooses every cluster
label. `cell-curator annotation import-consensus` only checks that run against this
run's cells, clusters and vocabulary and writes its labels. You never assign, edit or
overrule a label, and no LLM lane exists in the package.

## The network setting does not cover you

`knowledge.network` decides what the package may send: `offline`, `institution`
(localhost and allowed hosts), or `online`. It does not cover you. If your model is not local, everything you read from a run (cards, reports, score tables) leaves the institution. When the user asks for local-only work and your model is remote, say so before reading run outputs.

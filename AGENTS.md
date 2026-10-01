# cell-curator: instructions for coding agents

Annotate pre-clustered AnnData/MuData with evidence-gated, ontology-checked
labels. The workflow is in `skills/cell-curator/SKILL.md`. Read it before running
anything.

## Commands

- `cell-curator run scene|plan|execute|validate --config CONFIG`: the guided run.
- `cell-curator progress status --config CONFIG`: where the run is, and what is next.
- `cell-curator write-back preview|apply --config CONFIG --output OUT.h5ad`.
- Tests: `.venv/bin/python -m pytest -q`, plus `ruff check .` and `mypy`.

## Gates you must not route around

1. **Gate A.** Assumptions are approved by a person. Do not run
   `cell-curator assumptions approve` on your own initiative. Ask the user.
2. **Critics.** Every configured critic pass needs a finding or an attestation
   before write-back.
3. **Write-back.** `cell-curator write-back apply` asks a person to type the
   approval token in a terminal. Stop and ask the user to run it. Never pipe the
   token in.

## LLM evidence is evidence, not a verdict

An LLM, including you, may only break a near-tie that the measured markers left
open, and only through `evidence.llm.proposals`. It never overrules a call.

`knowledge.network` decides what may leave the machine: `offline`, `institution`
(localhost and allowed hosts), or `online`. Reading `cell-curator llm packet` output
sends it to your model provider. Do not read it unless the network is `online` or
the user set `evidence.llm.host_agent_local: true`.

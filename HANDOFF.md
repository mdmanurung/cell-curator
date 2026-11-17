# Handoff: cell-curator after the consensus removal

**Written:** 2026-11-17 · **Branch:** `main` (local only, nothing pushed, no tag) · **Status:** package labelling removed, consensus importer in, AIDA split A reviewed up to critics, PBMC3K waiting on you

## Check first
Run `git log --oneline -12`, `git status --short` and `squeue -u $USER -o "%i %j %T %M"`. If any shows work this file does not mention, say so first. The previous handoff (2026-11-10, v8 dev runs) is saved at `/exports/para-lipg-hpc/mdmanurung/cc-v8-work/HANDOFF-2026-11-10.md`.

## Goal
cell-curator checks, reviews and writes back labels. The `annotation-consensus` skill chooses every label. The package assigns none.

## What exists
- Flow: `run scene`, `run plan`, annotation-consensus with `plan/vocabulary.tsv`, `annotation collect-level`, `annotation import-consensus`, `annotation finalize`, review build, `run validate`, `write-back preview|apply`. Written in `README.md`, `skills/cell-curator/SKILL.md`, `docs/tutorial.md`.
- Config: `annotation.consensus_directory`; `guidance.hierarchy_levels: [L1]` only (a corpus config with a second level is refused at load).
- Removed: non-corpus pipeline, `decisions`, `audit_parents`, `refinement`, evidence lanes, LLM lane, `workflow/Snakefile`. The `markers` group stays.
- Tests at the merge: 678 passed, 2 GPU skipped; ruff and mypy clean.
- Container: `skills/annotation-consensus/container/`; SIF `01703be8d325…` (wheel `1083cee6…`) reproduces `expected_hashes.json` on three CPU types, with no network, and under `verify --demo`.
- Records: `plans/v8-consensus-audit.md` (12-row determinism audit), `plans/v8-acceptance-2026-11-17.md` (AIDA run), `plans/v8-agent-decisions.md`.

## Next action
Record critic findings for the AIDA and PBMC3K runs (`critic import-findings`, `critic reconcile`), then `run validate` and `write-back preview`. Both runs sit at the critic step.

## State
- AIDA split A, `cc-v8-work/acceptance-2026-11-17/runs/v8-accept-aida-split-a-001`: imported labels equal the consensus for 23 of 23 clusters; 4 `Unknown`; report built. Run commands from `cc-v8-work/acceptance-2026-11-17/` because this run's evidence receipt holds a relative config path (fixed for new runs in `5c79dd2`).
- PBMC3K, `.../runs/v8-accept-pbmc3k-001` (config `pbmc3k.yaml`): Gate A approved by the user 2026-11-17; imported labels equal the consensus for 8 of 8 clusters (3 `Unknown`); report and review packet built. No corpus review was needed: the plan uses `reviewed_v3.json`, which you reviewed on 2026-11-09. Record: `plans/v8-acceptance-2026-11-17.md` (merged to `main`).

## Locked decisions
- Per-choice asks off for v8; choices are logged in `plans/v8-agent-decisions.md` (user, 2026-11-09).
- Standing Gate A for v8 dev-set runs, reviewer Mikhael Manurung; novel sets need the user (user, 2026-11-10; renewed by a yes on 2026-11-17 for the AIDA run).
- Consensus values are the frozen demo values the user chose on 2026-11-13. Do not tune them on a run.

## Dead ends: do not redo
- Ontology back-off and looser vote gates for the L1 path; they only trade coverage for accuracy (`plans/v8-aida-dev-result-001.md`).
- Running from the live repo: `.venv` is an editable install of the main checkout's `src`. Use `PYTHONPATH=<worktree>/src` without `python -I`, which ignores it.
- Snakemake `annotation_review` after a failed first try: it deletes the outputs it was asked to build.

## Open items (yours)
- Critic findings for `run validate` (config is `network: offline`; the agent is a remote model, so no run output goes to a critic).
- `write-back apply` for any run, typed by a person.
- Per-cell flag rule (19.2% of AIDA, 32.6% of AIFI cells flagged), marker-share tie-break, state-call rule, when the novel sets are read.
- Push and tag.

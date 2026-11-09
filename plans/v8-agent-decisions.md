# v8 decisions made by the agent (user turned off per-choice asks, 2026-11-09)

The user chose "stop asking me per choice" for the v8 effort. Each row is a choice the
agent made alone, with the alternatives and the reason, for the user to review later.
The AGENTS.md gates (Gate A, critics, corpus review, write-back) are not covered.

| Date | Choice | Alternatives | Reason |
| --- | --- | --- | --- |
| 2026-11-10 | Tabula Sapiens and Cross-tissue run the package on their blood cells only, keeping the frozen `bench_cluster` labels (no re-clustering); every non-blood cell is scored as an abstention, because the blood-scoped corpus selects no record for any other tissue (checked: `snapshot.select` gives 34 records for blood, 0 for spleen). | (a) Pooled run with context tissue `blood`: gives blood records to spleen and lymph-node cells, against the user's blood-scope choice. (b) One run per tissue: non-blood runs select 0 records, so they equal abstention at higher cost. | The package reads one tissue context per run and has no per-cell tissue; (a) breaks the user's scope decision and (b) adds runs with an identical outcome. |
| 2026-11-10 | (user decision, recorded here) Standing Gate A approval for v8 dev-set runs only, reviewer Mikhael Manurung; a dev-only package flag lets unreviewed candidate corpora run on dev splits. | Ask per run; review every candidate corpus | Up to 100 approvals over 50 iterations. |

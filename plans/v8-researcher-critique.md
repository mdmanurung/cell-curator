# v8 strategy critique by three independent Opus researchers (2026-11-09)

Three Opus agents with literature search each read the v8 text and `AGENTS.md` alone
and wrote a critique. They saw no run output and no novel-set data. Their sources were
checked by PubMed metadata only, no full texts, and the marker genes they name come from
their own knowledge. Treat every point as a hypothesis until measured.
`unsupported: no inspected figure` applies to every claim about v8 performance.

## Agreement (counts of 3 researchers who raised the point)

| Weak point | Raised by | Checked here? |
| --- | --- | --- |
| Blood-only corpus has no tissue types (GC B, Tfh, resident T, tissue macrophages); cells get the nearest blood label | 3 of 3 | The corpus is blood only: yes (API) |
| Negatives learned in blood can veto true calls in tissue (CD69, dissociation genes) | 3 of 3 | No |
| Coarse truth labels make descendant scoring easy; the benchmark may mostly test lineage | 3 of 3 | Partly: 15% of TS, 36% of Indonesia, 53% of Cross-tissue match a corpus term by name |
| Cell-level Wilcoxon pseudo-replicates; top 10 may separate lineages, not siblings (Squair 2021, PMID 34584091) | 3 of 3 | No |
| One Leiden run at 1.0, no batch correction, can split by donor or tissue or merge siblings | 3 of 3 | No |
| Assay shift: Indonesia is 10x 5' v2 | 3 of 3, all said "unchecked" | Yes: the atlas is 10x 3' v3 only, so the corpus never saw 5' |
| No doublet or ambient-RNA step | 3 of 3 | No |
| Cross-tissue truth labels came from CellTypist | 2 of 3 | Stated in amendment 1 |
| The dev split is in-distribution, so 50 iterations cannot see tissue or assay shift | 1 of 3 | Reasoning holds: the dev donors are the same blood 3' v3 study |

## Agreement on what the package should do

- All three: name a type only at the depth the evidence reaches, else the parent or
  unresolved; require lineage before fine labels; require support in most donors; run
  QC, doublet and ambient flags the same way every time.
- All three: do not name states outside the corpus; do not merge, split or re-cluster to
  fit a label.
- One: corpus approval stays a person's gate (already `AGENTS.md` gate 3).

## Where they differ

Rank order of the weak points differs. Two put the blood-only vocabulary first; one put
the scoring granularity first. No researcher disagreed with another on a fact.

## What this changes

Nothing in the frozen v8 text. Proposed, each needing the user's decision:
1. A development set that is not in-distribution (tissue or 5' data), so the 50
   iterations can see the shifts the novel sets contain.
2. Report-only strata for the scoring: truth is itself a corpus term versus a coarse
   parent. The pass bar stays the v6 rules.
3. Diagnostics added to the learning output, no decision: per-donor support of each
   marker, and each negative's detection in 5' data where available.
4. Roadmap, not v8: doublet and ambient flags per cluster in the package.

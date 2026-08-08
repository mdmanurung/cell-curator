# Adaptive nested subclustering

## Parent audit routing

Flag a canonical parent for investigation when cell-level evidence shows guide
disagreement/entropy, a meaningful minority, incompatible programs, doublet
enrichment, donor/capture dominance, QC/complexity/sex structure, or a plausible
subtype mixture. A flag does not imply `ACCEPT SPLIT`.

## Stored-resolution families

For each higher stored resolution:

1. intersect each child with the canonical parent;
2. calculate child-to-parent purity and fraction of the parent;
3. enforce configured cell/donor/capture eligibility;
4. match children across adjacent resolutions by Jaccard on parent-restricted
   barcode sets;
5. collapse persistent matches into one family;
6. reject singleton/transient families unless a documented rare exception is
   approved.

Leiden labels are arbitrary and partitions are not formally nested. Compare
memberships, not child IDs.

## Parent-local fallback

Run at most one local pass per flagged parent after stored screening fails.
Reuse frozen preprocessing, transformation, batch covariate, declared assay
representations, graph construction, neighbor settings, clustering backend, and
stability protocol. The graph may be single-assay, WNN, or another frozen study
graph. Do not reselect global preprocessing.

For every local resolution record:

- cluster count and size distribution;
- configured cell/sample/capture eligibility per child;
- resampling stability across reference seeds;
- adjacent-resolution persistence;
- guide distributions only as descriptive hypotheses;
- identity, positive, negative, and incompatible program separation;
- configured sample/capture association and maximum dominance;
- assay-appropriate QC, complexity, sex, and batch effect sizes;
- doublet fraction and within-cell incompatible-program co-expression.

Select the lowest resolution that is stable **and** resolves at least two
eligible, biologically distinguishable groups. A stable one-cluster solution is
uninformative. A stable technical partition is not a biological split.

## Candidate registry

Stored and local children enter one schema and one evidence workflow. Each row
has a unique candidate ID, parent ID, source kind, frozen membership hash,
eligibility/stability receipts, closest comparator, blank biological labels,
and `FULL_TESTING` routing. Membership contains exactly one Boolean row per
parent cell per candidate.

Local discovery is not accepted directly. It must pass both binary comparisons
for every required evidence lane and the same gates as a stored child.

## Boundedness

Do not recursively split accepted children in the same run. If evidence supports
a further split, register it as a new reviewed run with a new frozen parent
contract.

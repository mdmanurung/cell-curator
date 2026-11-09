# Mapping rule: researcher free text to Cell Ontology ids (protocol v8)

Fixed before any researcher output exists. `map_and_agree.py` applies it; nobody edits
a label by hand.

**Ontology.** `benchmarks/data/aida-phase1-v2-pbmc-hvg2000/cl.obo`, `data-version:
releases/2026-06-08`, sha256 `7e19b4aef8e7fe7720b0c4474c6a5b8713c904ea9ab7b7d5d96f2cb28b8fc51f`
(the pin in the README).

**Rule, in order.**

1. Trim and collapse spaces, drop a final period. Case is ignored.
2. If the label is `unresolved`, `unknown`, `unassigned`, `na`, `n/a` or empty: `UNRESOLVED`
   (status `no_call`).
3. Look up the label as an exact Cell Ontology name or an EXACT synonym
   (`CellOntology.resolve`). If it ends in " cells" and fails, try the singular.
4. One match to a non-obsolete term: that CL id (status `exact`).
5. Otherwise `UNRESOLVED`, with status `unmatched` (no match), `ambiguous` (the name or
   synonym belongs to more than one term) or `obsolete`.

No fuzzy matching, no stemming, no BROAD, NARROW or RELATED synonyms, no manual fixes.

**The id the researcher states is not used.** It is kept in `stated_cl_id` for inspection.
Free text decides, so a wrong remembered id cannot change a result. A mismatch between
the stated id and the mapped id is counted and reported, not corrected.

**Agreement.** Per dataset, over clusters, with `UNRESOLVED` as its own category (all
three saying "unresolved" is agreement):

- `exact`: pairwise agreement (mean over the 3 researcher pairs of the fraction of
  clusters with the same CL id) and Fleiss' kappa.
- `coarse_lineage`: the same two measures after replacing each CL id with its deepest
  ancestor-or-self in the coarse list in `map_and_agree.py` (T cell, B cell, natural
  killer cell, innate lymphoid cell, plasma cell, monocyte, macrophage, dendritic
  cell, granulocyte, platelet, erythrocyte, hematopoietic stem cell, epithelial cell,
  endothelial cell, fibroblast). An id under none of them stays itself. This is the
  "parent level".
- `ancestor_or_equal`: pairwise only. Two labels agree if one id equals the other or is
  its ancestor in the pinned ontology.
- `exact_all_resolved` and `coarse_all_resolved`: the same measures on clusters that
  all three resolved, so abstention does not inflate them.
- `fraction_unresolved`: share of researcher-cluster calls that are `UNRESOLVED`.

Fleiss' kappa is `NaN` when chance agreement is 1 (every call is one category).
Unmapped labels are listed in `unmapped.tsv` and reported with the result.

## Open questions (proposed value, for the user to confirm before any packet is made)

Each choice below can change a headline number. None is settled.

1. **Genes per cluster in the packet.** Proposed 20 (the user's figure), ranked by
   Wilcoxon score, positive direction, cell level, no cell cap. More genes give more
   evidence and more room for a researcher to be led by noise.
2. **Lineage gene panel.** Proposed the 13 groups in `LINEAGE_PANEL` of
   `benchmarks/scripts/researcher_packets_v8.py` (74 genes, immune-heavy, with a
   few tissue and non-immune genes). Tabula Sapiens and Cross-tissue contain non-immune
   cells only if the subset keeps them; the subset is immune, so the non-immune group is small.
3. **Mapping source.** Proposed: free text decides, the stated CL id is ignored. The
   alternative is to trust a stated id that resolves.
4. **Coarse (parent) level.** Proposed the 15 terms in `COARSE_LINEAGE`. A direct
   `is_a` parent is not used because terms have several parents and kappa needs one category.
5. **Abstention in agreement.** Proposed: `UNRESOLVED` is a category, and the
   all-resolved measures are reported beside it.
6. **Unit of agreement.** Proposed: clusters, per dataset, equal weight. Weighting by
   cluster size is the alternative.
7. **Wilcoxon on X as given.** log1p CP10k, no scaling, as the packet text says. Ranking
   by score, not log fold change.

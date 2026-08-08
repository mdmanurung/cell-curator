# Decision gates and critic reconciliation

## Parent outcomes

Every canonical parent receives exactly one reviewed outcome:

- `ACCEPT SPLIT`: reproducible biological children pass all applicable gates;
- `RETAIN PARENT`: heterogeneity is continuous, weak, or unsupported;
- `DOUBLET/MIXED`: incompatible programs co-occur or multiplets dominate;
- `TECHNICAL/UNRESOLVED`: technical structure or missing evidence dominates.

## Capability-aware child gate

Require configured size/replication eligibility, parent nesting, stability,
positive identity/subtype evidence, expected-negative evidence, both binary
comparators, cell-level separation, applicable source corroboration, technical
checks, doublet/mixed assessment, alternatives, and confidence margin.

Each required evidence lane contributes `lane__<lane_key>__evidence_pass`.
`evidence.required_candidate_gates` supplies study-level gates. Do not require an
absent assay. Do fail a candidate when a claim requires a capability that was
not declared or validated.

Failure of one child does not invalidate a passing sibling; success of one child
does not rescue failed candidates. Build mappings from accepted candidate IDs.

## Critic and label gates

Critics inspect raw tables for hierarchy, negative evidence, doublet/mixed
explanations, technical/sample confounding, assay limitations, and ontology
semantics. Before labels require signed Gate A, annotation mode, callable
Knowledgebase, an ontology-compatible evidence path, completed gates,
reconciled critics, semantic validity, separate type/state, and human approval.
Keep label/CL fields blank in evidence-only mode.

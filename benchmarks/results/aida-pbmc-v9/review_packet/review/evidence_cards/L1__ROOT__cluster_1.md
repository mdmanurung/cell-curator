# L1 · ROOT · cluster 1

- Parent scope: `ROOT`
- Proposed cell type: `Unknown`
- Cell state: `not_assessed`
- Confidence: `0.242`
- Margin: `0.146`
- Normalized entropy: `0.933`
- Abstained: `True`
- Rationale: no direct annotation-capable lane supported the winner; required lanes lacked valid detection-dependent program evidence; expected-negative markers were violated; confidence or winner margin was insufficient; candidate evidence remained too entropic
- Rejected alternatives: hematopoietic stem cell | CD8-positive T cell | plasmacytoid dendritic cell

## Actual multimodal evidence

| label | score | probability | positive_coverage | negative_violation | lane_winners | multimodal_disagreement |
| --- | --- | --- | --- | --- | --- | --- |
| CD4-positive T cell | 0.2224623411893844 | 0.2419205299828843 | 1.0 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |
| hematopoietic stem cell | -0.7044624735911686 | 0.0957446229155912 | 0.2 | 0.3333333333333333 | transcriptome__lognorm=CD4-positive T cell | False |
| CD8-positive T cell | -0.7654334902763367 | 0.0900813771416578 | 0.6 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |
| plasmacytoid dendritic cell | -0.9902414580186208 | 0.0719452551391146 | 0.0 | 0.3333333333333333 | transcriptome__lognorm=CD4-positive T cell | False |
| platelet | -1.0419609745343523 | 0.0683188671627363 | 0.0 | 0.3333333333333333 | transcriptome__lognorm=CD4-positive T cell | False |
| innate lymphoid cell | -1.173765391111374 | 0.0598823348788583 | 0.7142857142857143 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |
| CD14-positive monocyte | -1.238837629556656 | 0.0561097342196779 | 0.0909090909090909 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |
| myeloid dendritic cell | -1.271268626054128 | 0.0543192304665558 | 0.1 | 0.3333333333333333 | transcriptome__lognorm=CD4-positive T cell | False |
| B cell | -1.271339386701584 | 0.0543153869386255 | 0.1111111111111111 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |
| CD16-positive monocyte | -1.364245980978012 | 0.0494964505804138 | 0.2727272727272727 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |

# L1 · ROOT · cluster 1

- Parent scope: `ROOT`
- Proposed cell type: `Unknown`
- Cell state: `not_assessed`
- Confidence: `0.224`
- Margin: `0.135`
- Normalized entropy: `0.946`
- Abstained: `True`
- Rationale: no direct annotation-capable lane supported the winner; required lanes lacked valid detection-dependent program evidence; expected-negative markers were violated; confidence or winner margin was insufficient; candidate evidence remained too entropic
- Rejected alternatives: CD8-positive T cell | innate lymphoid cell | platelet

## Actual multimodal evidence

| label | score | probability | positive_coverage | negative_violation | lane_winners | multimodal_disagreement |
| --- | --- | --- | --- | --- | --- | --- |
| CD4-positive T cell | 0.2396078556776046 | 0.2240769124876578 | 1.0 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |
| CD8-positive T cell | -0.6840944886207581 | 0.0889689122314558 | 0.6 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |
| innate lymphoid cell | -0.7425936609506607 | 0.0839136115414952 | 1.0 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |
| platelet | -0.8943019409974415 | 0.0721018393863661 | 0.0 | 0.3333333333333333 | transcriptome__lognorm=CD4-positive T cell | False |
| hematopoietic stem cell | -0.9275032778580984 | 0.0697472657143351 | 0.0 | 0.3333333333333333 | transcriptome__lognorm=CD4-positive T cell | False |
| plasmacytoid dendritic cell | -0.9832293291886648 | 0.0659668384067491 | 0.0 | 0.3333333333333333 | transcriptome__lognorm=CD4-positive T cell | False |
| gamma-delta T cell | -1.0065525323152542 | 0.0644460838058775 | 0.6666666666666666 | 1.0 | transcriptome__lognorm=CD4-positive T cell | False |
| myeloid dendritic cell | -1.0823666354020436 | 0.0597407795148891 | 0.1 | 0.3333333333333333 | transcriptome__lognorm=CD4-positive T cell | False |
| CD14-positive monocyte | -1.263086199760437 | 0.0498638004095686 | 0.0909090909090909 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |
| plasma cell | -1.2889749705791471 | 0.0485894547049944 | 0.375 | 0.5 | transcriptome__lognorm=CD4-positive T cell | False |

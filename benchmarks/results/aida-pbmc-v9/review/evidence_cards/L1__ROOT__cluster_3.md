# L1 · ROOT · cluster 3

- Parent scope: `ROOT`
- Proposed cell type: `Unknown`
- Cell state: `not_assessed`
- Confidence: `0.197`
- Margin: `0.045`
- Normalized entropy: `0.922`
- Abstained: `True`
- Rationale: no direct annotation-capable lane supported the winner; required lanes lacked valid detection-dependent program evidence; expected-negative markers were violated; confidence or winner margin was insufficient; candidate evidence remained too entropic
- Rejected alternatives: myeloid dendritic cell | CD14-positive monocyte | natural killer cell

## Actual multimodal evidence

| label | score | probability | positive_coverage | negative_violation | lane_winners | multimodal_disagreement |
| --- | --- | --- | --- | --- | --- | --- |
| platelet | 0.2878572046756744 | 0.1973440389211554 | 0.6666666666666666 | 1.0 | transcriptome__lognorm=platelet | False |
| myeloid dendritic cell | 0.0307463357845942 | 0.1526026870491002 | 0.5 | 0.3333333333333333 | transcriptome__lognorm=platelet | False |
| CD14-positive monocyte | -0.4121636748313904 | 0.0979961033873916 | 0.9090909090909092 | 1.0 | transcriptome__lognorm=platelet | False |
| natural killer cell | -0.4124849662184715 | 0.0979646231408426 | 0.9090909090909092 | 0.5 | transcriptome__lognorm=platelet | False |
| B cell | -0.6630835831165314 | 0.0762492674510312 | 0.5555555555555556 | 1.0 | transcriptome__lognorm=platelet | False |
| plasmacytoid dendritic cell | -0.7303529183069865 | 0.0712887458717757 | 0.1818181818181818 | 0.6666666666666666 | transcriptome__lognorm=platelet | False |
| myeloid cell | -0.7693160176277161 | 0.0685645319274571 | 0.8181818181818182 | 1.0 | transcriptome__lognorm=platelet | False |
| gamma-delta T cell | -1.210777536034584 | 0.0440935650831338 | 1.0 | 1.0 | transcriptome__lognorm=platelet | False |
| CD16-positive monocyte | -1.3371188640594482 | 0.0388602752599391 | 0.8181818181818182 | 1.0 | transcriptome__lognorm=platelet | False |
| hematopoietic stem cell | -1.3877474069595337 | 0.0369418104600113 | 0.4 | 1.0 | transcriptome__lognorm=platelet | False |

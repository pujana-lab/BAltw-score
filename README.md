# BAltw-score

R implementation for computing the **βAltw score**, a weighted gene expression signature that captures the anti-correlated relationship between TGFβ signaling and alternative end-joining (alt-EJ) DNA repair activity.

## Background

Loss of TGFβ signaling shifts DNA double-strand break repair toward alt-EJ, an error-prone pathway associated with genomic instability. TGFβ and alt-EJ gene expression signatures are anti-correlated across cancer types, and this relationship predicts patient outcome after genotoxic therapy (radiation, platinum-based chemotherapy).

The βAltw score refines this relationship by weighting (also unweighted version available) each signature gene according to its association with functional response, improving prediction of treatment response.

This work builds on:

- Liu Q, Palomero L, Moore J, et al. *Loss of TGFβ signaling increases alternative end-joining DNA repair that sensitizes to genotoxic therapies across cancer types.* Sci Transl Med. 2021;13(580):eabc4465. [PMID: 33568520](https://pubmed.ncbi.nlm.nih.gov/33568520)
- Guix I, Liu Q, Pujana MA, et al. *Validation of Anticorrelated TGFβ Signaling and Alternative End-Joining DNA Repair Signatures that Predict Response to Genotoxic Cancer Therapy.* Clin Cancer Res. 2022;28(7):1372-1382. [PMID: 35022323](https://pubmed.ncbi.nlm.nih.gov/35022323)

## Contents

- `compute_bAltw_score.R` — computes the βAltw score from gene expression data
- `data/` — supporting data files

## Usage

```r
source("compute_bAltw_score.R")
# add specific usage example once function signature/inputs are confirmed
```

## Citation

If you use this code, please cite the publications above.

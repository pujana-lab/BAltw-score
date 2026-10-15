# BAltw-score

R implementation for computing the **βAltw score**, a weighted gene expression signature that captures the anti-correlated relationship between TGFβ signaling and alternative end-joining (alt-EJ) DNA repair activity.

> **Note:** This repository contains a standalone implementation of the βAltw scoring method. For the full analysis pipeline used to derive and validate the score across all datasets in the paper (HNSC explants, TCGA, GSE41613, GSE26712, NCI-60, etc.), see [bAltw-score-validation](https://github.com/pujana-lab/bAltw-score-validation).

## Background

Loss of TGFβ signaling shifts DNA double-strand break repair toward alt-EJ, an error-prone pathway associated with genomic instability. TGFβ and alt-EJ gene expression signatures are anti-correlated across cancer types, and this relationship predicts patient outcome after genotoxic therapy (radiation, platinum-based chemotherapy).

The βAltw score refines this relationship by weighting each signature gene according to its association with functional response, improving prediction of treatment response. The script also provides an optional unweighted *sum* (`use_weights = FALSE`), which is **not** the original ssGSEA-based βAlt method of Liu et al. 2021.

This work builds on:
- Liu Q, Palomero L, Moore J, et al. *Loss of TGFβ signaling increases alternative end-joining DNA repair that sensitizes to genotoxic therapies across cancer types.* Sci Transl Med. 2021;13(580):eabc4465. [PMID: 33568520](https://pubmed.ncbi.nlm.nih.gov/33568520)
- Guix I, Liu Q, Pujana MA, et al. *Validation of Anticorrelated TGFβ Signaling and Alternative End-Joining DNA Repair Signatures that Predict Response to Genotoxic Cancer Therapy.* Clin Cancer Res. 2022;28(7):1372-1382. [PMID: 35022323](https://pubmed.ncbi.nlm.nih.gov/35022323)

## Contents

- `compute_bAltw_score.R` — calculates weighted βAltw and reports missing genes
- `data/` — TGFβ and alt-EJ gene sets, weights, Ensembl mapping, and `alias_map.txt` (symbol aliases)

## Usage

1. Install the R packages `dplyr`, `data.table`, and `readxl` (`readxl` is needed for Excel input).
2. Prepare a CSV, TSV, semicolon-delimited text, or Excel `.xlsx`/`.xls` expression table. Genes are **rows**, samples are **columns**, and the **first column** holds the gene symbol or unversioned Ensembl ID. Expression values must be log2-transformed and **z-scored per gene across the cohort's samples**. A full-transcriptome or signature-only matrix can be used.
3. Edit these settings near the top of `compute_bAltw_score.R`:

```r
path_to_input_files  <- "./data/"
path_to_dataset_file <- "./my_zscored_expression.csv"
path_to_output_file  <- NA  # or "./bAltw_scores.csv"
use_ensembl_ids      <- FALSE  # TRUE for unversioned ENSG... IDs
use_weights          <- TRUE   # FALSE: unweighted sum, not ssGSEA
```

4. From the repository directory, run:

```r
source("compute_bAltw_score.R")
head(weighted_scores)
missing_genes  # genes absent from the input or with entirely missing values
```

`weighted_scores` contains `SampleID`, `TGFBUPgeneset_ws`, `ALTEJgeneset_ws`, `bAlt_wscore`, separate `tertileW`/`TertileW` (three bins) and `quartileW`/`QuartileW` (four bins; only extreme quartiles receive text labels). If `path_to_output_file` is set, scores are also saved as a CSV.

The script reports signature coverage and missing genes; it stops if **either entire signature is absent**. Dataset rows corresponding to the same signature gene are averaged **after** alias mapping, with a warning (for array data, collapsing probes **before** z-scoring is preferable). Isolated missing z-scores contribute zero to the weighted sum and are reported. Missing signature genes are omitted from the sum, so coverage should be documented.

**Cohort dependence:** the geometric score depends on sample-level signature extrema; both the z-scoring and tertile/quartile assignments depend on the cohort. Scores calculated on different sample populations are not directly interchangeable.

## Citation

If you use this code, please cite the publications above.

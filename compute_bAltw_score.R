# =============================================================================
# compute_bAltw_score.R
#
# Calculation of the weighted beta-Alt (betaAltw) score described in:
# Guix et al., Clinical Cancer Research (2022), 28:1372-1382.
#
# INPUT:
# Gene expression matrix with genes as rows and samples as columns.
# The first column must contain gene identifiers (HGNC symbols or Ensembl IDs).
# Expression values must be log2-transformed and z-scored per gene across
# samples. Supported formats: CSV, TSV, delimited text, and Excel.
#
# OUTPUT:
# weighted_scores : Per-sample TGFbeta and alt-EJ weighted signature scores,
#                   betaAltw scores, and tertile/quartile classifications.
# missing_genes   : Signature genes absent or entirely missing in the dataset.
#
# Higher betaAltw scores indicate higher alt-EJ and lower TGFbeta
# signature expression.
#
# See README.md for usage instructions and additional information.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(readxl)
})

# ---- Settings ---------------------------------------------------------------
path_to_input_files  <- "./data/"     # folder containing the bAlt data files
path_to_dataset_file <- "./my_zscored_expression.csv"   # expression file
path_to_output_file  <- NA            # e.g. "./bAltw_scores.csv"; NA = do not write

use_ensembl_ids <- FALSE  # TRUE if row identifiers are Ensembl gene IDs (ENSG...)
use_weights     <- TRUE   # FALSE = all gene factors set to 1 (unweighted weighted-sum
                          # score). This is NOT the ssGSEA-based βAlt method.

# ---- Reference files --------------------------------------------------------
TGFBlist  <- unique(scan(file.path(path_to_input_files, "TGFBUPgeneset.txt"), character(), quiet = TRUE))
ALTEJlist <- unique(scan(file.path(path_to_input_files, "ALTEJgeneset.txt"),  character(), quiet = TRUE))

tgfb_weighted   <- data.table::fread(file.path(path_to_input_files, "tgfb_weighted.txt"))
Alt_ej_Weighted <- data.table::fread(file.path(path_to_input_files, "altej_weighted.txt"))
weights_factors <- rbind(tgfb_weighted, Alt_ej_Weighted)
tgfb_keys  <- tgfb_weighted$geneID
altej_keys <- Alt_ej_Weighted$geneID

if (use_ensembl_ids) {
  gene_map <- data.table::fread(file.path(path_to_input_files, "AltEJ_TGFb_mapping.txt"))
  symbol_to_ens <- setNames(trimws(gene_map$EnsemblID), trimws(gene_map$Symbol))
  TGFBlist  <- unique(unname(symbol_to_ens[TGFBlist]))
  ALTEJlist <- unique(unname(symbol_to_ens[ALTEJlist]))
  weights_factors$geneID <- unname(symbol_to_ens[weights_factors$geneID])
  tgfb_keys  <- unname(symbol_to_ens[tgfb_keys])
  altej_keys <- unname(symbol_to_ens[altej_keys])
  message("Ensembl mode: reference files converted to Ensembl IDs.")
} else {
  message("Symbol mode: using gene symbols as identifiers.")
}

# ---- Expression matrix ------------------------------------------------------
# Preserve gene IDs in a separate vector. Assigning them as data-frame row
# names at import would reject duplicate IDs *before* they can be averaged.
ext <- tolower(tools::file_ext(path_to_dataset_file))
if (ext == "xlsb") stop(".xlsb is not supported by readxl; please save as .xlsx, CSV or TSV.")
if (ext %in% c("xlsx", "xls", "xlsm")) {
  header_info <- readxl::read_excel(path_to_dataset_file, n_max = 0)
  if (ncol(header_info) < 2) stop("Input must contain a gene-ID column and at least one sample column.")
  types <- c("text", rep("numeric", ncol(header_info) - 1))
  input_table <- as.data.frame(readxl::read_excel(path_to_dataset_file, col_types = types), check.names = FALSE)
} else {
  input_table <- as.data.frame(data.table::fread(path_to_dataset_file, check.names = FALSE),
                               check.names = FALSE)
}
if (ncol(input_table) < 2 || nrow(input_table) == 0)
  stop("Input must contain gene IDs and expression values for at least one sample.")
gene_ids <- trimws(as.character(input_table[[1]]))
if (anyNA(gene_ids) || any(!nzchar(gene_ids)))
  stop("The first column contains missing or empty gene identifiers.")
genie <- input_table[, -1, drop = FALSE]
if (anyDuplicated(colnames(genie))) stop("Sample names (column headers) must be unique.")
if (any(colnames(genie) %in% c("Gene", "geneID", "Weight", "Factor", "factor")))
  stop("Sample names cannot be Gene, geneID, Weight, Factor or factor.")
if (!all(vapply(genie, is.numeric, logical(1))))
  stop("Expression columns must be numeric; check delimiters and input values.")
if (!any(is.finite(as.matrix(genie)))) stop("No finite expression values found.")

# Input sanity checks (warnings only)
if (any(genie > 80, na.rm = TRUE)) warning("Data has values > 80. Likely not z-scored.")
if (min(as.matrix(genie), na.rm = TRUE) >= 0) warning("No negative values found. Likely not z-scored.")
if (any(abs(rowMeans(genie, na.rm = TRUE)) > 1e-5, na.rm = TRUE))
  warning("Gene (row) means are not ~0. Not z-scored per gene, or samples removed after z-scoring.")

# ---- Select signature genes and harmonize symbols ---------------------------
selected <- which(gene_ids %in% TGFBlist | gene_ids %in% ALTEJlist)
assay <- genie[selected, , drop = FALSE]
assay$Gene <- gene_ids[selected]

if (!use_ensembl_ids) {
  am <- data.table::fread(file.path(path_to_input_files, "alias_map.txt"))
  alias_map <- setNames(am$weight_key, am$dataset_symbol)
  assay$Gene <- ifelse(assay$Gene %in% names(alias_map), alias_map[assay$Gene], assay$Gene)
}

# Rows that are entirely NA (e.g. zero-variance genes after z-scoring) carry no information.
sample_cols_assay <- setdiff(names(assay), "Gene")
all_na <- rowSums(!is.na(assay[, sample_cols_assay, drop = FALSE])) == 0
if (any(all_na)) {
  warning("Signature gene rows with no non-missing values were removed: ",
          paste(assay$Gene[all_na], collapse = ", "))
  assay <- assay[!all_na, , drop = FALSE]
}
# Partial NA: rowSums(na.rm = TRUE) below treats a missing z-score as 0 (the gene mean).
n_na <- sum(is.na(assay[, sample_cols_assay]))
if (n_na > 0) warning(n_na, " missing values in signature genes; each is treated as z = 0 (gene mean).")

# Two dataset rows mapping to the same signature gene (NUDT1 and MTH1, two probes or IDs): average.
dup_genes <- unique(assay$Gene[duplicated(assay$Gene)])
if (length(dup_genes) > 0) {
  warning("Multiple dataset rows map to the same signature gene; rows averaged for: ",
          paste(dup_genes, collapse = ", "))
  assay <- aggregate(assay[, sample_cols_assay, drop = FALSE], by = list(Gene = assay$Gene),
                     FUN = function(v) mean(v, na.rm = TRUE))
}

# ---- Missing-gene report ----------------------------------------------------
missing_tgfb  <- setdiff(tgfb_keys,  assay$Gene)
missing_altej <- setdiff(altej_keys, assay$Gene)
message(sprintf("TGFb genes matched: %d/%d; alt-EJ genes matched: %d/%d",
                length(tgfb_keys) - length(missing_tgfb), length(tgfb_keys),
                length(altej_keys) - length(missing_altej), length(altej_keys)))
if (length(missing_tgfb)  > 0) message("Missing TGFb genes: ",  paste(missing_tgfb,  collapse = ", "))
if (length(missing_altej) > 0) message("Missing alt-EJ genes: ", paste(missing_altej, collapse = ", "))
missing_genes <- data.frame(geneID = c(missing_tgfb, missing_altej),
                            signature = c(rep("TGFb", length(missing_tgfb)), rep("altEJ", length(missing_altej))))

# Each signature must contain at least one usable gene; otherwise a numeric
# score would misleadingly treat the entirely absent signature as zero.
if (length(missing_tgfb) == length(tgfb_keys) || length(missing_altej) == length(altej_keys)) {
  absent <- c(if (length(missing_tgfb) == length(tgfb_keys)) "TGFb",
              if (length(missing_altej) == length(altej_keys)) "alt-EJ")
  stop("Cannot compute bAltw score: no usable ", paste(absent, collapse = " or "),
       " signature genes matched. Check gene identifiers, use_ensembl_ids, and missing values.")
}

# ---- Weighted signature scores ----------------------------------------------
genes2 <- merge(weights_factors, assay, by.x = "geneID", by.y = "Gene", all = TRUE)
# factor = 1 + Weight if Weight > 0, else 1 (Guix et al. 2022, Methods)
genes2$factor <- ifelse(genes2$Weight > 0, 1 + genes2$Weight, 1)
if (!use_weights) genes2$factor <- 1
genes2 <- as.data.frame(genes2)
rownames(genes2) <- genes2$geneID

sample_cols <- setdiff(names(genes2), c("geneID", "Weight", "factor", "Factor", "Gene"))
scores <- sweep(genes2[, sample_cols, drop = FALSE], 1, genes2$factor, `*`)
rownames(scores) <- genes2$geneID
scores <- as.data.frame(t(scores))

scores$TGFBUPgeneset_ws <- rowSums(scores[, which(colnames(scores) %in% tgfb_keys),  drop = FALSE], na.rm = TRUE)
scores$ALTEJgeneset_ws  <- rowSums(scores[, which(colnames(scores) %in% altej_keys), drop = FALSE], na.rm = TRUE)

# ---- βAltw (geometric formula; sign flipped so high = high alt-EJ / low TGFβ)
scores$bAlt_wscore <- -1 * (
  sqrt((max(scores$ALTEJgeneset_ws) - scores$ALTEJgeneset_ws)^2 +
       (min(scores$TGFBUPgeneset_ws) - scores$TGFBUPgeneset_ws)^2) -
  sqrt((min(scores$ALTEJgeneset_ws) - scores$ALTEJgeneset_ws)^2 +
       (max(scores$TGFBUPgeneset_ws) - scores$TGFBUPgeneset_ws)^2))

# ---- Tertiles and quartiles (kept in separate columns) ----------------------
scores <- scores %>%
  dplyr::mutate(tertileW  = ntile(bAlt_wscore, 3),
                TertileW  = dplyr::case_when(tertileW == 1 ~ "low βalt", tertileW == 2 ~ "mid βalt",
                                             tertileW == 3 ~ "high βalt"),
                quartileW = ntile(bAlt_wscore, 4),
                QuartileW = dplyr::case_when(quartileW == 1 ~ "low βalt", quartileW == 4 ~ "high βalt",
                                             TRUE ~ NA_character_))
scores$SampleID <- rownames(scores)
weighted_scores <- scores %>%
  dplyr::select(SampleID, TGFBUPgeneset_ws, ALTEJgeneset_ws, bAlt_wscore,
                tertileW, TertileW, quartileW, QuartileW)

if (!is.na(path_to_output_file)) {
  write.csv(weighted_scores, path_to_output_file, row.names = FALSE)
  message("Scores written to ", path_to_output_file)
}

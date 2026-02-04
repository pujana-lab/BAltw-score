# in order to run this script review the next steps after loading the libraries.

library("dplyr")
library("readr")
library("data.table")
library("readxl")

# This script requires normalized count data that has been log-transformed and Z-scored
# gene tables (.txt) and weighting factors (.txt) should be in the path_to_input_files directory.
# normalized_counts %>% log2() %>% scale(count_matrix, center = T, scale = T) # to scale, genes are columns samples are rows
# transpose t() genes to rows
# convert to a data.frame()
# This file should be in the path_to_dataset_file directory. Separated by tab, comma, or semi colon.

# Set this to TRUE if your input datasetuses Ensembl IDs (ENSG...)
# Set to FALSE if it uses Gene Symbols (TP53, MRE11...)
use_ensembl_ids <- F

path_to_input_files  = "./data/"
path_to_dataset_file = "./"

# Load gene tables, .txt gene lists contain complete gene-sets as well as gene alias where relevant.
TGFBlist  = scan(paste0(path_to_input_files, "TGFBUPgeneset.txt"), as.character(), quiet = TRUE)
ALTEJlist = scan(paste0(path_to_input_files, "ALTEJgeneset.txt"), as.character(), quiet = TRUE)

# Load weights
tgfb_weighted   <- data.table::fread(paste0(path_to_input_files, "tgfb_weighted.txt"))
Alt_ej_Weighted <- data.table::fread(paste0(path_to_input_files, "altej_weighted.txt"))
weights_factors <- rbind(tgfb_weighted, Alt_ej_Weighted)

# Bifurcation based on execution mode (Ensembl vs Symbol)
if (use_ensembl_ids) {
  #Convert the genesets and weights to Ensembl ID
  path_to_ensembl_map  = paste0(path_to_input_files,"AltEJ_TGFb_mapping.txt")
  if (!file.exists(path_to_ensembl_map)) {
    stop("Ensembl mode is ON, but mapping file not found at: ", path_to_ensembl_map)
  }
  # Load Map (Expects header: Symbol, EnsemblID)
  gene_map <- data.table::fread(path_to_ensembl_map)
  
  # Create look-up vector
  symbol_to_ens <- setNames(gene_map$EnsemblID, gene_map$Symbol)
  
  # Translate genesets
  TGFBlist  <- symbol_to_ens[TGFBlist]
  ALTEJlist <- symbol_to_ens[ALTEJlist]
  
  # Translate weights table
  weights_factors$geneID <- symbol_to_ens[weights_factors$geneID]
  message("Ensembl Mode: Reference files converted to Ensembl IDs.")
} else {
  message("Symbol Mode: Using Gene Symbols as identifiers.")
}

# Load normalized expression matrix
ext <- tolower(tools::file_ext(path_to_dataset_file))
if(ext %in% c("xlsx", "xls", "xlsm", "xlsb")){
  header_info <- read_excel(path_to_dataset_file, n_max = 0)
  num_cols <- ncol(header_info)
  types <- c("text", rep("numeric", num_cols - 1))
  genie <- read_excel(path_to_dataset_file, col_types = types)
  genie <- column_to_rownames(genie, var = colnames(genie)[1])
} else {
  genie <- data.frame(data.table::fread(path_to_dataset_file), row.names = 1)  
}

# 1. Check for high values (unlikely in Z-scores)
if (any(genie > 80, na.rm = T)) warning("Data has values > 80. Likely not Z-scored.")
# 2. Check for non-negativity (Z-scores should center at 0)
if (min(genie, na.rm = T) >= 0) warning("No negative values found. Likely not Z-scored.")
# 3. Check if column means are ~0.
if (any(abs(rowMeans(genie, na.rm = TRUE)) > 1e-5, na.rm = TRUE)) warning("Column means are not ~0. Not Z-scored or samples removed afterwards")

# Select genes in the BAlt gene sets
assay <- genie[which(rownames(genie) %in% TGFBlist | rownames(genie) %in% ALTEJlist),] 
assay$Gene <- rownames(assay)

#  Only run Alias Mapping if we are in SYMBOL mode.
# If in Ensembl mode, IDs are unique and don't need aliasing.
if (!use_ensembl_ids) {
  
  alias_map <- c(
    "APEX2"    = "APE2",
    "HRAS"     = "HRAS1",
    "PRPF19"   = "PRP19",
    "KAT5"     = "TIP60",
    "ARHGAP32" = "RICS",
    "CYTH1"    = "PSCD1",
    "NUDT1"    = "MTH1",
    "NABP2"    = "NABP2",
    "CTGF"     = "CCN2",
    "CCDC99"   = "SPDL1",
    "TMEPAI"   = "PMEPA1",
    "C19orf40" = "FAAP24",
    "RAD51L3"  = "RAD51D",
    "MRE11"   = "MRE11A"
  )
  
  # If gene is in map, swap it; otherwise keep original
  assay$Gene <- ifelse(assay$Gene %in% names(alias_map), alias_map[assay$Gene], assay$Gene)
}

# Merge weights with assay data
genes2 <- merge(weights_factors, assay, by.x="geneID", by.y="Gene", 
                all = TRUE) # record missing genes here.  all = T to identify missing genes for accounting
# Calculate the TGFB and ALTEJ refined scores (refined from HNSC explants) in each sample, by multiplying each gene by its factor (1+W or 1 if W<0).
genes2$factor <- ifelse(genes2$Weight < 0, 1, genes2$Weight)
genes2$factor <- ifelse(genes2$factor == 1, 1, 1 + genes2$factor)

# rename rows as the weights_factors geneIDs
rownames(genes2) <- genes2$geneID
genes2 <- as.data.frame(genes2)

# Sync column names and perform multiplication (Sweep)
metadata_cols <- c("geneID", "Weight", "factor", "Factor", "Gene") 
sample_cols <- setdiff(names(genes2), metadata_cols)

# If merge(all=T) introduced NAs (missing genes), sweep/sapply will create NAs.

# Multiply factor by each sample's gene expression to get weighted gene scores
scores <- sweep(genes2[, sample_cols, drop=FALSE], 1, genes2$factor, `*`)
# Ensure rownames are set for the output
rownames(scores) <- genes2$geneID

# transpose and convert to data.frame to add the gene-set score variables (samples to rows)
scores <- as.data.frame(t(scores))

# sum adjusted gene scores in each gene-set for each sample, removing genes not represented with na.rm argument.
scores$TGFBUPgeneset_ws <- rowSums(scores[ , which(colnames(scores) %in% TGFBlist)], na.rm=TRUE)
scores$ALTEJgeneset_ws  <- rowSums(scores[ , which(colnames(scores) %in% ALTEJlist)], na.rm=TRUE)

# calculate bAlt score from TGFB and ALTEJ
scores$bAlt_wscore <- sqrt((max(scores$ALTEJgeneset_ws) - scores$ALTEJgeneset_ws)^2 +
                             (min(scores$TGFBUPgeneset_ws) - scores$TGFBUPgeneset_ws)^2) - 
  sqrt((min(scores$ALTEJgeneset_ws) - scores$ALTEJgeneset_ws)^2 +
         (max(scores$TGFBUPgeneset_ws) - scores$TGFBUPgeneset_ws)^2)

scores$bAlt_wscore <- scores$bAlt_wscore * -1 

# generate BAlt tertiles
scores <- scores %>% dplyr::mutate(quartileW = ntile(bAlt_wscore, 3))
scores$QuartileW <- ifelse(scores$quartileW == 1, "low βalt", ifelse(scores$quartileW == 3, "high βalt", NA))
# generate BAlt quartiles
scores <- scores %>% dplyr::mutate(quartileW = ntile(bAlt_wscore, 4))
scores$QuartileW <- ifelse(scores$quartileW == 1, "low βalt", ifelse(scores$quartileW == 4, "high βalt", NA))

# this dataframe can be merged to clinical or experimenal data for downstream analysis
scores$SampleID <- rownames(scores)
weighted_scores <- scores %>% dplyr::select(TGFBUPgeneset_ws, ALTEJgeneset_ws, bAlt_wscore, 
                                            quartileW, QuartileW, SampleID)



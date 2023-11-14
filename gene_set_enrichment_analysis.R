# Script Name: gene_set_enrichment_analysis.R
#
# Description:
# This script conducts Gene Set Enrichment Analysis (GSEA) on differential gene expression results.
#
# Inputs:
# 1. Combined differential gene expression file with symbol annotations.
# 2. GMT files (gene sets).
#
# Outputs:
# 1. Combined GSEA results file.
#


# Libraries & file/dir names
library("fgsea")
library("dplyr")
library("data.table")

deseq_in_fn <- "data/differential_expression_annotated_combined.csv"
gmt_dir <- "data/gmt_files"
gsea_out_fn <- "data/fgsea_combined_new.csv"

set.seed(42)

# Check input data exists
if (!file.exists(deseq_in_fn) || !dir.exists(gmt_dir)) {
  stop("Input file or directory not found")
}

# Load gmt files
gene_sets <-
  lapply(list.files(gmt_dir, full.names = T),
         fgsea::gmtPathways)
names(gene_sets) <- list.files(gmt_dir)

# Load deseq results file
deseq_df <- read.csv(deseq_in_fn, header = T)

# Drop all columns other than HGNC, t-score and subset
deseq_df <- deseq_df[, c("HGNC.symbol", "stat", "subset")]

# Remove unmapped genes
deseq_df <- deseq_df[deseq_df$HGNC.symbol != "", ]

# Split into seperate dfs on subset
deseq_dfs <-
  split(deseq_df,
        factor(deseq_df$subset, levels = (unique(deseq_df$subset))))

# Sort by t-score, descending and average t-score for duplicates
sorted_deseq_dfs <- lapply(deseq_dfs, function(x) {
  x %>%
    group_by(HGNC.symbol) %>%
    summarise(stat = mean(stat, na.rm = TRUE)) %>%
    arrange(desc(stat))
})

# Create t-score rank vectors
rank_vectors <-
  lapply(sorted_deseq_dfs, function(x)
    setNames(x$stat, x$HGNC.symbol))

# Perform GSEA on rank vectors (all combos of rank vectors and pathways)
fgsea_res <- lapply(gene_sets,
                    function(x)
                      lapply(rank_vectors, function(y)
                        fgsea(
                          pathways = x,
                          stats = y,
                          minSize = 15,
                          maxSize = 300
                        )))

# Condense to a combined df per category
category_dfs <-
  lapply(fgsea_res, function(x)
    dplyr::bind_rows(x, .id = "subset"))

# Condense to one combined df
combined_df <- dplyr::bind_rows(category_dfs, .id = "category")
combined_df$category <- gsub(".gmt", "", combined_df$category)

# Save combined df
fwrite(combined_df, gsea_out_fn)
# Script Name: differential_gene_expression_analysis.R
#
# Description:
# This script conducts differential gene expression analysis.
#
# Inputs:
# 1. Estimated counts files (downloaded from GeneLab).
# 2. Metadata files (downloaded from GeneLab).
# 3. Manually curated metadata table.
#
# Outputs:
# 1. Raw counts files.
# 2. Normalized counts files.
# 3. Differential gene expression files.
# 4. Single combined differential gene expression file with symbol annotations.
#


################
### PREAMBLE ###
################

# Libs
library(purrr)
library(dplyr)
library(tximport)
library(DESeq2)
library(stringr)
library(data.table)
library(readr)
library(biomaRt)

# Define GeneLab RNA-Seq datasets to be analysed
dataset_ids <- c('238', '239', '240', '241', '254')

# Define exclusion criteria for dataset conditions
exclusion_criteria <- "VIV|BSL|FLT_1G"

# Define file paths and names
data_dir <- 'data'
data_subdirs <- file.path(data_dir, "downloaded_counts", dataset_ids)
subsets_fn <- file.path(data_dir, "subset_table.csv")
deseq_out_fn <- file.path(data_dir, "deseq.csv")
all_data_files <- list.files(data_subdirs, full.names = TRUE)
mdata_file_paths <-
  all_data_files[grepl("*metadata_GLDS*", all_data_files)]

# Set seed and save session info
set.seed(42)
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(data_dir, "session_info.txt"))


#####################
### LOAD METADATA ###
#####################

# Read in metadata for each dataset
metadata <- map(set_names(mdata_file_paths, dataset_ids), ~ {
  unzip(.x, exdir = dirname(.x), overwrite = TRUE)
  
  metadata_file <-
    list.files(dirname(.x), pattern = "s_.*txt$", full.names = TRUE)
  
  read.delim(metadata_file, stringsAsFactors = TRUE) %>%
    mutate(Condition = gsub("^Mmus_|-|*_Rep.*$", "", Sample.Name))
})

# Drop unused conditions
if (!is.null(exclusion_criteria)) {
  metadata <-
    map(metadata, function(x)
      dplyr::filter(x,!grepl(exclusion_criteria, Condition)))
}
metadata <- map(metadata, droplevels)


########################
### LOAD COUNTS DATA ###
########################

# Helper function to load counts files for a dataset
import_counts_files <- function(dataset_idx) {
  # Get names of counts files for the dataset id
  data_folders <- list.dirs('data/downloaded_counts', recursive = FALSE)
  files <-
    list.files(data_folders[[dataset_idx]], full.names = TRUE)
  counts_files <- files[grepl("*.genes.results$", files)]
  
  # Filter down to only counts files for samples in metadata
  filtered_files <- c()
  for (sample_name in metadata[[dataset_idx]]$Sample.Name) {
    sample_name <-
      gsub(" ", "", sample_name) # Remove any extra whitespaces from sample name
    
    file <- counts_files[grepl(sample_name, counts_files)]
    names(file) <- sample_name
    filtered_files <- c(filtered_files, file)
  }
  
  # Import counts data
  imported <-
    tximport(filtered_files,
             type = "rsem",
             txIn = FALSE,
             txOut = FALSE)
  imported$length[imported$length == 0] <-
    1 # Correct for zero gene length issue (https://support.bioconductor.org/p/92763/)
  
  return(imported)
}

# Call helper function to load estimated counts
estimated_counts <- vector("list", length(dataset_ids))
for (idx in 1:length(dataset_ids)) {
  estimated_counts[[idx]] <- import_counts_files(idx)
}
names(estimated_counts) <- dataset_ids


############################
### SAVE RAW COUNTS DATA ###
############################

raw_counts <-
  lapply(estimated_counts, function(x)
    as.data.frame(x$counts))

if (!file.exists("data/raw_counts")) {
  dir.create("data/raw_counts")
}

raw_counts <-
  lapply(raw_counts, tibble::rownames_to_column, "ENSEMBL")

mapply(function(x, y)
  fwrite(x, file = paste0("data/raw_counts/", y, "_raw_counts.csv")), x = raw_counts, y = dataset_ids)


###############
### RUN DGE ###
###############

# Create dds objects
dds_list <- vector("list", length(dataset_ids))
dds_list <-
  mapply(
    function(x, y)
      DESeqDataSetFromTximport(
        txi = x,
        colData = y,
        design = ~ Condition
      ),
    x = estimated_counts,
    y = metadata
  )
names(dds_list) <- dataset_ids

# Reduce down to only ENSEMBL genes (remove ERCC)
dds_list <-
  lapply(dds_list, function(x)
    x[grepl("ENSMUSG", row.names(x)),])

# Filter out low count, counts must be > 0 in at least the minimum number of samples corresponding to a single condition
dds_list <-
  lapply(dds_list, function(x)
    x[rowSums(counts(x) > 0) >= (min(table(x$Condition))),])

# Run differential expression analysis
dds_list <- lapply(dds_list, DESeq)


##############################
### SAVE NORMALIZED COUNTS ###
##############################

normalized_counts <-
  lapply(dds_list, function(x)
    as.data.frame(counts(x, normalized = TRUE)))

# Create subdir for normalized counts if non-existent
if (!file.exists("data/normalized_counts")) {
  dir.create("data/normalized_counts")
}

normalized_counts <-
  lapply(normalized_counts, tibble::rownames_to_column, "ENSEMBL")

# Write normalized counts to files
mapply(function(x, y)
  fwrite(
    x,
    file = paste0("data/normalized_counts/", y, "_normalized_counts.csv")
  ), x = normalized_counts, y = dataset_ids)


###################
### DGE RESULTS ###
###################

# Read in custom metadata file/sample table
data_subsets <- read.csv(file = subsets_fn)

# Calculate differential expression results
diff_express_results <- lapply(1:nrow(data_subsets), function(i) {
  results(
    object = dds_list[[data_subsets$dataset_id[i]]],
    contrast = c(
      "Condition",
      data_subsets$flight_condition[i],
      data_subsets$ground_condition[i]
    ),
    cooksCutoff = FALSE,
    independentFiltering = FALSE
  )
})
names(diff_express_results) <- data_subsets$subset_name

# Convert deseq results into data frames
diff_express_results_dfs <-
  mapply(
    function(x, y)
      data.frame(row.names(x), x),
    x = diff_express_results,
    y = names(diff_express_results),
    SIMPLIFY = F
  )

# Create subdir for deseq results if non-existent
if (!file.exists("data/differential_expression_results")) {
  dir.create("data/differential_expression_results")
}

# write deseq results to files, one per data subset
mapply(
  function(x, y)
    fwrite(
      x,
      file = paste0(
        "data/differential_expression_results/",
        y,
        "_differential_expression_results.csv"
      )
    ),
  x = diff_express_results_dfs,
  y = names(diff_express_results_dfs)
)


######################
### SYMBOL MAPPING ###
######################

# Helper function for ortholog mapping
gene_symbol_mapping <- function(df, mart_mouse, mart_human) {

  # Rename column
  colnames(df)[1] <- "mouse_ensembl_gene_id"
  
  # Get mouse gene symbols
  genes_mouse_symbols <- getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'),
                               filters='ensembl_gene_id',
                               values=df$mouse_ensembl_gene_id,
                               mart=mart_mouse)
  
  
  # Get human orthologs one-to-one
  human_orthologs <- getLDS(
    attributes = c("mgi_symbol"), 
    filters = "mgi_symbol",
    values = genes_mouse_symbols$mgi_symbol,
    mart = mart_mouse,
    attributesL = c("hgnc_symbol"), 
    martL = mart_human
  )
  
  # Join together
  merged_df <- left_join(genes_mouse_symbols, human_orthologs, by = c("mgi_symbol" = "MGI.symbol"), relationship = "many-to-many")
  
  df <- left_join(df, merged_df, by=c("mouse_ensembl_gene_id" = "ensembl_gene_id"))
  return(df)

}

# ensembl marts
mart_mouse <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org//")
mart_human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org//")

diff_express_results_mapped <- mapply(gene_symbol_mapping, diff_express_results_dfs, MoreArgs = list(mart_mouse = mart_mouse, mart_human = mart_human), SIMPLIFY = FALSE)
 

###################################################
### COMBINED RESULTS FILE & MULTI-MISSION GENES ###
###################################################

# combine each diff_express_results dataframe into one
combined_df <- imap_dfr(diff_express_results_mapped, ~ .x %>% 
                          mutate(subset = .y) %>% 
                          dplyr::select(-baseMean, -lfcSE))

# connect results with metadata
combined_df_with_meta <- inner_join(combined_df, data_subsets, by = c("subset" = "subset_name"))

# find multi-mission genes
mission_genes <- combined_df_with_meta %>%
  filter(padj <= 0.1) %>% # filter to significant genes
  filter(mgi_symbol != "") %>% # ensure mapped to mgi symbol
  filter(HGNC.symbol != "") %>% # ensure mapped to hgnc symbol
  group_by(mgi_symbol) %>% 
  summarise(n_missions = n_distinct(mission)) %>%
  ungroup() %>%
  filter(n_missions >= 2) # must be present(significant & mapped) in multiple missions

# Add the isKeyGene column
combined_df <- combined_df %>%
  mutate(isKeyGene = if_else(mgi_symbol %in% mission_genes$mgi_symbol, TRUE, FALSE))

# write to csv
fwrite(combined_df, "data/differential_expression_annotated_combined.csv")

# write unique multi-mission gene to .txt
filtered_df <- combined_df[combined_df$isKeyGene == TRUE, ]
mgi_key_genes <- filtered_df[!duplicated(filtered_df$mgi_symbol), ]$mgi_symbol
hgnc_key_genes <- filtered_df[!duplicated(filtered_df$HGNC.symbol), ]$HGNC.symbol

writeLines(mgi_key_genes, "key_genes_mgi.txt")
writeLines(hgnc_key_genes, "key_genes_hgnc.txt")




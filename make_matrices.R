# Libs
library(tidyr)
library(dplyr)

#####################
## Pathway results ##
#####################

# Load GSEA results as dataframe
fgsea_df <- read.csv("data/fgsea_combined.csv")
fgsea_df$subset <- factor(fgsea_df$subset, levels = unique(fgsea_df$subset))

# 1. NES matrix

# Convert dataframe into matrix format
fgsea_nes_matrix <- fgsea_df %>%
  dplyr::select(pathway, subset, NES) %>%          
  spread(key = subset, value = NES)  

# Set rownames as pathways
rownames(fgsea_nes_matrix) <- fgsea_nes_matrix$pathway
fgsea_nes_matrix$pathway <- NULL

# convert to matrix
fgsea_nes_matrix <- as.matrix(fgsea_nes_matrix)

# 2. Pathway significance matrix

# Convert dataframe into matrix format
fgsea_sig_matrix <- fgsea_df %>%
  dplyr::select(pathway, subset, padj) %>%         
  spread(key = subset, value = padj)         

# Set rownames as pathways
rownames(fgsea_sig_matrix) <- fgsea_sig_matrix$pathway
fgsea_sig_matrix$pathway <- NULL

# convert to matrix
fgsea_sig_matrix <- as.matrix(fgsea_sig_matrix)


###################
## Deseq results - HGNC ##
###################

# Load deseq results as dataframe
deseq_df <- read.csv("data/differential_expression_annotated_combined.csv")
deseq_df$subset <- factor(deseq_df$subset, levels = unique(deseq_df$subset))

## HGNC

# Remove unmapped HGNC symbols, i.e. ""
deseq_df <- deseq_df %>%
  filter(HGNC.symbol != "")

# Average tscore and significance for duplicate HGNC symbols within the same subset
deseq_df <- deseq_df %>%
  group_by(HGNC.symbol, subset) %>%
  summarise(stat = mean(stat, na.rm = TRUE), 
            padj = mean(padj, na.rm = TRUE)) %>%
  ungroup()

# 1. tscore matrix

# Convert dataframe into matrix format
deseq_tscore_matrix <- deseq_df %>%
  dplyr::select(HGNC.symbol, subset, stat) %>%          
  spread(key = subset, value = stat)  

# Set rownames as pathways
deseq_tscore_matrix <- as.data.frame(deseq_tscore_matrix)
rownames(deseq_tscore_matrix) <- deseq_tscore_matrix$HGNC.symbol
deseq_tscore_matrix$HGNC.symbol <- NULL

# convert to matrix
deseq_tscore_matrix_hgnc <- as.matrix(deseq_tscore_matrix)

# 2. gene significance matrix 

# Convert dataframe into matrix format
deseq_sig_matrix <- deseq_df %>%
  dplyr::select(HGNC.symbol, subset, padj) %>%         
  spread(key = subset, value = padj)         

# Set rownames as pathways
deseq_sig_matrix <- as.data.frame(deseq_sig_matrix)
rownames(deseq_sig_matrix) <- deseq_sig_matrix$HGNC.symbol
deseq_sig_matrix$HGNC.symbol <- NULL

# convert to matrix
deseq_sig_matrix_hgnc <- as.matrix(deseq_sig_matrix)

###################
## Deseq results - MGI ##
###################

# Load deseq results as dataframe
deseq_df <- read.csv("data/differential_expression_annotated_combined.csv")
deseq_df$subset <- factor(deseq_df$subset, levels = unique(deseq_df$subset))

## MGI

# Remove unmapped mgi symbols, i.e. ""
deseq_df <- deseq_df %>%
  filter(mgi_symbol != "")

# Average tscore and significance for duplicate mgi symbols within the same subset
deseq_df <- deseq_df %>%
  group_by(mgi_symbol, subset) %>%
  summarise(stat = mean(stat, na.rm = TRUE), 
            padj = mean(padj, na.rm = TRUE)) %>%
  ungroup()

# 1. tscore matrix

# Convert dataframe into matrix format
deseq_tscore_matrix <- deseq_df %>%
  dplyr::select(mgi_symbol, subset, stat) %>%          
  spread(key = subset, value = stat)  

# Set rownames as pathways
deseq_tscore_matrix <- as.data.frame(deseq_tscore_matrix)
rownames(deseq_tscore_matrix) <- deseq_tscore_matrix$mgi_symbol
deseq_tscore_matrix$mgi_symbol <- NULL

# convert to matrix
deseq_tscore_matrix_mgi <- as.matrix(deseq_tscore_matrix)

# 2. gene significance matrix 

# Convert dataframe into matrix format
deseq_sig_matrix <- deseq_df %>%
  dplyr::select(mgi_symbol, subset, padj) %>%         
  spread(key = subset, value = padj)         

# Set rownames as pathways
deseq_sig_matrix <- as.data.frame(deseq_sig_matrix)
rownames(deseq_sig_matrix) <- deseq_sig_matrix$mgi_symbol
deseq_sig_matrix$mgi_symbol <- NULL

# convert to matrix
deseq_sig_matrix_mgi <- as.matrix(deseq_sig_matrix)



#################
## Write Files ##
#################

# Create directory if it doesn't exist
dir.create("matrices", showWarnings = FALSE)


write.csv(fgsea_nes_matrix, file = "matrices/fgsea_nes_matrix.csv", row.names = TRUE)

write.csv(fgsea_sig_matrix, file = "matrices/fgsea_sig_matrix.csv", row.names = TRUE)

write.csv(deseq_tscore_matrix_hgnc, file = "matrices/deseq_tscore_matrix_hgnc.csv", row.names = TRUE)

write.csv(deseq_sig_matrix_hgnc, file = "matrices/deseq_sig_matrix_hgnc.csv", row.names = TRUE)

write.csv(deseq_tscore_matrix_mgi, file = "matrices/deseq_tscore_matrix_mgi.csv", row.names = TRUE)

write.csv(deseq_sig_matrix_mgi, file = "matrices/deseq_sig_matrix_mgi.csv", row.names = TRUE)


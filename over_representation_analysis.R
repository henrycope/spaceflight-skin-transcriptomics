# Libs 
library("clusterProfiler")
library("org.Mm.eg.db")

# Load deseq results
deseq_tscore_matrix <- read.csv("matrices/deseq_tscore_matrix_mgi.csv", row.names = 1, check.names = F)

# Get all gene names (rownames of matrix)
all_genes_mgi <- row.names(deseq_tscore_matrix)

# Load key genes symbols
key_genes_mgi <- read.csv("key_genes_mgi.txt", header = F)[[1]]

# Perform ORA on key genes - using all genes as background
ora_results_mgi <- enrichGO(gene         = key_genes_mgi,
                            universe     = all_genes_mgi, # using all genes as background
                            OrgDb        = org.Mm.eg.db,
                            keyType      = "SYMBOL",
                            ont          = "BP", 
                            pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)

simplify_ora_results_mgi <- as.data.frame(simplify(ora_results_mgi, cutoff = 0.5))


simplify_ora_results_mgi$Percentage <- (simplify_ora_results_mgi$Count / 189) * 100

simplify_ora_results_mgi <- simplify_ora_results_mgi %>%
  mutate(Label = paste0(tools::toTitleCase(Description)))


###################################################
###################################################
###################################################
###################################################

# # Open deseq results file
df_deseq <- read.csv("data/differential_expression_annotated_combined.csv",
                     header = T)


sig_deseq_df <- df_deseq[df_deseq$padj <= 0.1, ]


################# MHU2


mhu2_only <- sig_deseq_df[sig_deseq_df$subset %in% c("238_GCvsFLT_JC", "238_GCvsFLT_JCwFOS", "239_GCvsFLT_JC", "239_GCvsFLT_JCwFOS"), ]


unique_mgi_238_GCvsFLT_JC <- unique(mhu2_only[mhu2_only$subset == "238_GCvsFLT_JC", ]$mgi_symbol)
unique_mgi_238_GCvsFLT_JCwFOS <- unique(mhu2_only[mhu2_only$subset == "238_GCvsFLT_JCwFOS", ]$mgi_symbol)
unique_mgi_239_GCvsFLT_JCC <- unique(mhu2_only[mhu2_only$subset == "239_GCvsFLT_JC", ]$mgi_symbol)
unique_mgi_239_GCvsFLT_JCwFOS <- unique(mhu2_only[mhu2_only$subset == "239_GCvsFLT_JCwFOS", ]$mgi_symbol)

# Finding intersections pairwise
intersect_1_2 <- intersect(unique_mgi_238_GCvsFLT_JC, unique_mgi_238_GCvsFLT_JCwFOS)
intersect_1_3 <- intersect(unique_mgi_238_GCvsFLT_JC, unique_mgi_239_GCvsFLT_JCC)
intersect_1_4 <- intersect(unique_mgi_238_GCvsFLT_JC, unique_mgi_239_GCvsFLT_JCwFOS)
intersect_2_3 <- intersect(unique_mgi_238_GCvsFLT_JCwFOS, unique_mgi_239_GCvsFLT_JCC)
intersect_2_4 <- intersect(unique_mgi_238_GCvsFLT_JCwFOS, unique_mgi_239_GCvsFLT_JCwFOS)
intersect_3_4 <- intersect(unique_mgi_239_GCvsFLT_JCC, unique_mgi_239_GCvsFLT_JCwFOS)

# Combin the intersections
combined_intersections <- c(intersect_1_2, intersect_1_3, intersect_1_4, intersect_2_3, intersect_2_4, intersect_3_4)

# Get unique elements
unique_combined_intersections <- unique(combined_intersections)

# Perform ORA on key genes - using all genes as background
ora_results_mgi_mhu2 <- enrichGO(gene         = unique_combined_intersections,
                            universe     = all_genes_mgi, # using all genes as background
                            OrgDb        = org.Mm.eg.db,
                            keyType      = "SYMBOL",
                            ont          = "BP", 
                            pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)

simplify_ora_results_mgi_mhu2 <- as.data.frame(simplify(ora_results_mgi_mhu2, cutoff = 0.5))

simplify_ora_results_mgi_mhu2 <- simplify_ora_results_mgi_mhu2 %>%
  mutate(Label = paste0(tools::toTitleCase(Description)))


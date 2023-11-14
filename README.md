# spaceflight-skin-transcriptomics

Execution Order:

1. download_datasets.R - downloads skin datasets from GeneLab
2. differential_gene_expression_analysis.R - Performs differential gene expression analysis from raw counts data
3. make_matrices.R - helper script to make matrices used for heatmap plots
4. over_representation_analysis.R - over representation analysis on cross-mission and MHU-2 gene lists
5. gene_set_enrichment_analysis.R - perform gene set enrichment analysis
6. gene_set_enrichment_analysis_mitochondrial.R - gene set enrichment analysis mitochondrial pathways only

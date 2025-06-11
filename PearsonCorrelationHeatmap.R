# Load required packages
# if (!requireNamespace("CEMiTool", quietly = TRUE)) {
#     install.packages("BiocManager")
#     BiocManager::install("CEMiTool")
# }
library(CEMiTool)

# Load expression data (genes as rows, samples as columns)
# expr_matrix should be a data.frame or matrix object
# pheno_data should be a data.frame with sample info (sample names should match column names of expr_matrix)

# Example:
# expr_matrix <- read.csv("expression_matrix.csv", row.names = 1)
# pheno_data <- read.csv("phenotype_data.csv")

# Run CEMiTool
cem <- cemitool(expr=batch_corrected, annot=pheno_hnVScopd_clean, filter=TRUE, plot=TRUE, verbose=TRUE)

# Generate report and plots
generate_report(cem, directory="CEMiTool_Report")


BiocManager::install(c(
  "DOSE",
  "fgsea",
  "clusterProfiler",
  "enrichplot",
  "CEMiTool"
))


# Unable to install CEMiTool directly from Bioconductor thus alternativly using alter method
# Compute Pearson correlation between samples
cor_matrix <- cor(batch_corrected, method = "pearson")
cor_matrix_genes <- cor(t(batch_corrected), method = "pearson")

# Load and plot
library(pheatmap)
# Save to PNG
png("temp/sample_correlation_heatmap.png", width = 1200, height = 1000, res = 150)
pheatmap::pheatmap(cor_matrix,
         main = "Sample Correlation Heatmap (Pearson)",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         display_numbers = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

annotation_col <- data.frame(
  Condition = pheno_hnVScopd$phenotype
)
rownames(annotation_col) <- colnames(cor_matrix)
# Save to PNG
png("temp/sample_annotated_correlation_heatmap.png", width = 1200, height = 1000, res = 150)
pheatmap::pheatmap(cor_matrix,
         annotation_col = annotation_col,
         main = "Sample Correlation with Condition",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         display_numbers = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()
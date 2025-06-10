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

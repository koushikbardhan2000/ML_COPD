# Differential gene expression analysis with limma
# 
library(limma)

# Step 1: Log2 transformation
expr_hnVScopd_log2 <- log2(expr_hnVScopd + 1)

# Step 2: Quantile normalization (optional but common for microarray)
expr_norm <- normalizeBetweenArrays(expr_hnVScopd_log2, method = "quantile")

# Step 3: Build phenotype design matrix
group <- factor(pheno_hnVScopd$phenotype)
design <- model.matrix(~ 0 + group)
colnames(design) <- make.names(levels(group))

# Step 4: Linear model fitting and contrast setup
fit <- lmFit(expr_norm, design)
contrast_matrix <- makeContrasts(SmokerwithCOPD_vs_Healthy = Smoker.with.COPD - Healthy.Non.Smoker, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Step 5: Get top DEGs
top_degs <- topTable(fit2, number = Inf, adjust.method = "BH")


# Filter significant DEGs
deg_filtered <- top_degs[top_degs$adj.P.Val < 0.05 & abs(top_degs$logFC) > 1, ]

# View results
head(deg_filtered)

# Load required package
library(ggplot2)

# Transpose expression matrix: samples as rows, genes as columns
expr_t <- t(expr_hnVScopd)

# Perform PCA
pca_result <- prcomp(expr_t, scale. = TRUE)

# Create data frame for plotting
pca_df <- data.frame(pca_result$x)
pca_df$Phenotype <- pheno_hnVScopd$phenotype  # Assuming this matches the sample order

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA: Smoker with COPD vs Healthy Non-Smoker",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
  theme_minimal() +
  scale_color_manual(values = c("Healthy Non-Smoker" = "forestgreen", "Smoker with COPD" = "firebrick"))



library(ggplot2)
library(pheatmap)

# === Volcano Plot ===
# Add gene names as a column if not already present
deg_filtered$Gene <- rownames(deg_filtered)

# Create significance labels
deg_filtered$threshold <- with(deg_filtered,
                               ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                                      ifelse(logFC > 1, "Upregulated", "Downregulated"), "Not Significant"))

# Volcano plot
volcano_plot <- ggplot(deg_filtered, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Smoker with COPD vs Healthy Non-Smoker",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(legend.title = element_blank())

print(volcano_plot)

# === Heatmap ===
# Select top DEGs (e.g., top 50)
top_genes <- rownames(deg_filtered[deg_filtered$adj.P.Val < 0.05 & abs(deg_filtered$logFC) > 1, ][1:50, ])

# Subset expression data for these genes
expr_subset <- expr_hnVScopd[top_genes, ]

# Scale rows (genes) for better visualization
expr_scaled <- t(scale(t(expr_subset)))

# Sample annotations for columns
annotation_col <- data.frame(Phenotype = factor(pheno_hnVScopd$phenotype))
rownames(annotation_col) <- pheno_hnVScopd$GSM_IDs

# Heatmap
pheatmap(expr_scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         main = "Top 50 DEGs Heatmap")


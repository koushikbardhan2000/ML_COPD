# Load required packages
library(limma)
library(ggplot2)
library(ggrepel)
library(ggfortify)

# STEP 1: Log2 transform and normalize expression data
# expr_log2 <- log2(expr_hnVScopd + 1)
# expr_norm <- normalizeBetweenArrays(expr_log2, method = "quantile")
keepSamples <- pheno_hnVScopd$GSM_IDs
batch_corrected <- batch_corrected[,keepSamples]
# STEP 2: Prepare phenotype and design matrix
group <- factor(pheno_hnVScopd$phenotype)
design <- model.matrix(~ 0 + group)
colnames(design) <- make.names(levels(group))  # syntactically valid

# STEP 3: Linear model and contrasts
fit <- lmFit(batch_corrected, design)
contrast_matrix <- makeContrasts(SmokerwithCOPD_vs_Healthy = Smoker.with.COPD - Healthy.Non.Smoker, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# STEP 4: Get DEGs
deg_table <- topTable(fit2, number = Inf, adjust.method = "BH")
deg_table$Gene <- rownames(deg_table)

# STEP 5: Volcano Plot
deg_table$Significance <- "Not Significant"
deg_table$Significance[deg_table$logFC > 1 & deg_table$adj.P.Val < 0.05] <- "Upregulated"
deg_table$Significance[deg_table$logFC < -1 & deg_table$adj.P.Val < 0.05] <- "Downregulated"

ggplot(deg_table, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Smoker with COPD vs Healthy Non-Smoker",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(plot.title = element_text(hjust = 0.5))

# STEP 6: PCA Plot
pca_data <- t(batch_corrected)  # transpose to samples x genes
autoplot(prcomp(pca_data, scale. = TRUE),
         data = pheno_hnVScopd,
         colour = 'phenotype',
         shape = 'phenotype') +
  theme_minimal() +
  ggtitle("PCA Plot: Samples by Phenotype")



# Load package
library(pheatmap)

# STEP 1: Select top 50 DEGs by adjusted p-value
top_degs <- deg_table[order(deg_table$adj.P.Val), ][1:50, ]
top_gene_ids <- top_degs$Gene

# STEP 2: Extract expression data for top DEGs
expr_top_degs <- batch_corrected[top_gene_ids, ]

# STEP 3: Prepare annotation for columns (samples)
annotation_col <- data.frame(Phenotype = pheno_hnVScopd$phenotype)
rownames(annotation_col) <- colnames(expr_top_degs)

# STEP 4: Generate heatmap
pheatmap(expr_top_degs,
         annotation_col = annotation_col,
         scale = "row",                    # z-score normalization by row
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 50 DEGs Heatmap: Smoker with COPD vs Healthy Non-Smoker")

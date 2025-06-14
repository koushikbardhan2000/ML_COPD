# Create a sample group mapping
sample_groups <- data.frame(
  SampleID = colnames(expr_hnVScopd),
  Group = c(rep("HN", 161), rep("COPD", 92))  # Adjust based on your actual group sizes
)
# Extract gene symbols
gene_symbols <- sapply(strsplit(rownames(expr_hnVScopd), " /// "), `[`, 1)

# Identify valid (non-NA and non-duplicated) entries
valid <- !is.na(gene_symbols) & !duplicated(gene_symbols)

# Subset the expression matrix and assign new rownames
expr_hnVScopd <- expr_hnVScopd[valid, ]
rownames(expr_hnVScopd) <- gene_symbols[valid]



expr_hnVScopd_log <- log2(expr_hnVScopd + 1)


library(limma)

# Design matrix
group <- factor(sample_groups$Group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit model
fit <- lmFit(expr_hnVScopd_log, design)
contrast.matrix <- makeContrasts(COPDvsHN = COPD - HN, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Top DE genes
topTable(fit2, adjust="fdr", number=20)


pca <- prcomp(t(batch_corrected), scale. = TRUE)
plot(pca$x[,1:2], col=as.factor(sample_groups$Group), pch=16)




expr_log2_norm <- normalizeBetweenArrays(expr_hnVScopd_log, method = "quantile")
# Transpose: samples as rows
expr_t <- t(expr_log2_norm)

# Scale the data
expr_scaled <- scale(expr_t)
pca <- prcomp(expr_scaled, center = TRUE, scale. = TRUE)
metadata <- pheno_hnVScopd


library(ggplot2)

# Create PCA data frame
pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     Batch = metadata$GSE_ID,
                     phenotype = metadata$phenotype)

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch, shape = phenotype)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Merged GEO Data (Batch Effect Visualization)",
       x = "PC1", y = "PC2")




# Batch effect removal
library(sva)

# Batch correction
batch_corrected <- ComBat(dat = as.matrix(expr_log2_norm), batch = metadata$phenotype, par.prior = TRUE)

# Redo PCA
pca_corrected <- prcomp(t(batch_corrected), center = TRUE, scale. = TRUE)


library(ggplot2)

# Create PCA data frame
pca_df <- data.frame(PC1 = pca_corrected$x[,1],
                     PC2 = pca_corrected$x[,2],
                     Batch = metadata$GSE_ID,
                     phenotype = metadata$phenotype)

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch, shape = phenotype)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Merged GEO Data (Batch Effect Visualization)",
       x = "PC1", y = "PC2")





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
expr_top_degs <- expr_norm[top_gene_ids, ]

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

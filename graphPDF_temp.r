# This file is 1st part of the COPD analysis using LASSO and various ML models.
# After this file, you should have run the COPD_new_LASSO_ML.R script to generate the required datasets.
# To clear the global environment in R
rm(list = ls())
########################
# Pre-processing starts
########################
matrixes <- list.files(path = "data_required/", full.names = T, pattern = ".gz")


# Read and merge all files by "ID_REF"
merged_data <- lapply(matrixes, function(file) {
  # Read file
  df <- read.table(file, header = TRUE, sep = "\t", comment.char = "!", fill = TRUE, quote = "", stringsAsFactors = FALSE)
  # Keep only relevant columns
  df
})

# View merged result
# head(merged_data)


# Rename the first column in each data frame to "ID" for uniform merging
merged_data <- lapply(merged_data, function(df) {
  colnames(df)[1] <- "ID"
  return(df)
})

# Now merge all data frames by "ID" using Reduce
final_merged_data <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), merged_data)

# View merged data
# head(final_merged_data)
dim(final_merged_data)


# Clean up column names
clean_names <- function(names_vec) {
  names_vec <- gsub("^X\\.", "", names_vec)         # Remove leading "X."
  names_vec <- gsub("\\.\\.x$|\\.\\.y$|\\.$", "", names_vec)  # Remove trailing "..x", "..y", or "."
  return(names_vec)
}

# Apply to your data
colnames(final_merged_data) <- clean_names(colnames(final_merged_data))
names(final_merged_data)

# Check for duplicates
# dup_cols <- colnames(final_merged_data)[duplicated(colnames(final_merged_data))]
# print(dup_cols)

# Remove duplicated columns, keeping the first occurrence
# final_merged_data <- final_merged_data[, !duplicated(colnames(final_merged_data))]

# Remove double quotes from 'ID' column
final_merged_data$ID <- gsub('"', '', final_merged_data$ID)

# soft file for annotation
soft <- read.delim("data_required/GPL570-55999.txt")
soft <- as.data.frame(soft[,c(1,11)])


# Merge gene symbols from soft into final_merged_data based on ID
final_annotated_data <- merge(soft, final_merged_data, by = "ID")

final_annotated_data$cleaned_names <- sub(" ///.*", "", final_annotated_data$Gene.Symbol)

# Remove rows with duplicated Gene.Symbols, keeping the first occurrence
final_annotated_data_unique <- final_annotated_data[!duplicated(final_annotated_data$cleaned_names) & 
                                                      !is.na(final_annotated_data$cleaned_names), ]


# Assign 'ID' column as rownames
rownames(final_annotated_data_unique) <- final_annotated_data_unique$cleaned_names

# Drop the 'ID' column
final_annotated_data_unique$ID <- NULL
final_annotated_data_unique$Gene.Symbol <- NULL
final_annotated_data_unique$cleaned_names <- NULL

dim(final_annotated_data_unique)
# write.csv(final_annotated_data_unique,"output/final_annotated_data_unique.csv")
# final_annotated_data_unique <- read.csv("output/final_annotated_data_unique.csv", row.names = 1)
head(final_annotated_data_unique[1:5, 1:5])
dim(final_annotated_data_unique)

# Preprocessing phenotype data
# Load required libraries
library(GEOquery)
library(Biobase)
library(dplyr)
library(stringr)

# Function to clean and set row names
clean_and_set_rownames <- function(df) {
  df$ID <- str_replace_all(df$ID, '"', '')  # Remove quotation marks
  rownames(df) <- df$ID
  df$ID <- NULL  # Optionally remove the ID column after setting row names
  return(df)
}

# List of GEO datasets
gse_ids <- c("GSE11784", "GSE13896", "GSE11906", "GSE37768", "GSE130928")

# Function to annotate samples based on dataset
annotate_samples <- function(gse_id) {
  gse <- getGEO(gse_id, GSEMatrix = TRUE)
  if (length(gse) > 1) idx <- grep("GPL", attr(gse, "names")) else idx <- 1
  gse_data <- gse[[idx]]
  pheno_data <- pData(gse_data)
  
  # Annotate phenotype
  pheno_data$phenotype <- case_when(
    grepl("non-smoker", pheno_data$title, ignore.case = TRUE) |
      grepl("tissue_NS", pheno_data$title, ignore.case = TRUE) |
      grepl("\\bNS\\b", pheno_data$title, ignore.case = TRUE) ~ "Healthy Non-Smoker",
    
    grepl("COPD", pheno_data$title, ignore.case = TRUE) ~ "Smoker with COPD",
    
    grepl("smoker", pheno_data$title, ignore.case = TRUE) & 
      !grepl("COPD", pheno_data$title, ignore.case = TRUE) ~ "Smoker",
    
    grepl("SYMs", pheno_data$title, ignore.case = TRUE) |
      grepl("\\bS\\b", pheno_data$title, ignore.case = TRUE) |
      grepl("tissue_S", pheno_data$title, ignore.case = TRUE) |
      grepl("LODs", pheno_data$title, ignore.case = TRUE) ~ "Smoker",
    
    TRUE ~ "Unknown"
  )
  
#   pheno_data$phenotype <- case_when(
#   # Healthy non-smokers
#   grepl("non-smoker", pheno_data$title, ignore.case = TRUE) |
#     grepl("tissue_NS", pheno_data$title, ignore.case = TRUE) |
#     grepl("\\bNS\\b", pheno_data$title, ignore.case = TRUE) ~ "Healthy Non-Smoker",
  
#   # Early COPD samples are to be treated as Smokers
#   grepl("early-COPD", pheno_data$title, ignore.case = TRUE) ~ "Smoker",
  
#   # Smoker with COPD
#   grepl("COPD", pheno_data$title, ignore.case = TRUE) ~ "Smoker with COPD",
  
#   # Generic smokers
#   grepl("smoker", pheno_data$title, ignore.case = TRUE) & 
#     !grepl("COPD", pheno_data$title, ignore.case = TRUE) ~ "Smoker",
  
#   # S label implies Smoker
#   grepl("\\bS\\b", pheno_data$title, ignore.case = TRUE) |
#     grepl("tissue_S", pheno_data$title, ignore.case = TRUE) ~ "Smoker",
  
#   # Default
#   TRUE ~ "Unknown"
# )

  # Add GSE_ID column
  pheno_data$GSE_ID <- gse_id
  
  return(pheno_data[, c("GSE_ID", "title", "phenotype")])
}

# Annotate samples for all datasets
annotated_samples_list <- lapply(gse_ids, annotate_samples)

# Combine all annotations into a single data frame
all_annotated_samples <- bind_rows(annotated_samples_list)


# View result
# View(all_annotated_samples)
dim(all_annotated_samples)
table(all_annotated_samples$phenotype)

# List of Unknown samples
Unknown <- data.frame(all_annotated_samples[all_annotated_samples$phenotype == "Unknown", ])
# View(Unknown)
# write.csv(Unknown, "output/Unknown samples removed.csv", row.names = F)


# Remove unknown phenotypes
# all_annotated_samples_clean <- all_annotated_samples %>%
#   filter(phenotype != "Unknown")

all_annotated_samples_withIDs <- data.frame("GSM_IDs" = row.names(all_annotated_samples), all_annotated_samples)


# Remove '...numbers' from GSM IDs
all_annotated_samples_withIDs$GSM_IDs <- gsub("\\.\\.\\..*", "", all_annotated_samples_withIDs$GSM_IDs)
# write.csv(all_annotated_samples_withIDs, "output/all_annotated_samples_withIDs.csv", row.names = T)

# View the combined annotations
# View(all_annotated_samples_withIDs)


# Remove rows with duplicates
all_annotated_samples_withIDs_unique <- all_annotated_samples_withIDs[!duplicated(all_annotated_samples_withIDs$GSM_IDs) & 
                                                         !is.na(all_annotated_samples_withIDs$GSM_IDs), ]

# Sort the dataframe by phenotype
all_annotated_samples_withIDs_unique_sorted <- all_annotated_samples_withIDs_unique[order(all_annotated_samples_withIDs_unique$phenotype), ]
row.names(all_annotated_samples_withIDs_unique_sorted) <- seq_len(nrow(all_annotated_samples_withIDs_unique_sorted))

# Final phenotype dataframe
phenotype <- all_annotated_samples_withIDs_unique_sorted[, c(1, 2, 4)]

# remove Unknown samples
phenotype <- phenotype[phenotype$phenotype != "Unknown", ]
head(phenotype[1:3, 1:3])
dim(phenotype)
table(phenotype$phenotype)
# write.csv(phenotype, "final/final_phenotype_final.csv", row.names = F) # nolint




# 1. Keep only columns in final_annotated_data_unique that are in phenotype$GSM_IDs
common_samples <- intersect(phenotype$GSM_IDs, colnames(final_annotated_data_unique))
filtered_expression <- final_annotated_data_unique[, common_samples]


# 2. Sort phenotype by phenotype column (e.g., alphabetical order) (already done before but just to be sure)
phenotype_sorted <- phenotype[order(phenotype$phenotype), ]
table(phenotype_sorted$phenotype)

# 3. Reorder the columns of expression data to match sorted phenotype
filtered_expression <- filtered_expression[, phenotype_sorted$GSM_IDs]
# write.csv(filtered_expression,"filtered_expression.csv")

# Subset GSM IDs by phenotype group
hn_ids <- phenotype_sorted$GSM_IDs[phenotype_sorted$phenotype == "Healthy Non-Smoker"]
smoker_ids <- phenotype_sorted$GSM_IDs[phenotype_sorted$phenotype == "Smoker"]
copd_ids <- phenotype_sorted$GSM_IDs[phenotype_sorted$phenotype == "Smoker with COPD"]

# Create subset expression matrices
expr_hn <- filtered_expression[, hn_ids]
expr_smoker <- filtered_expression[, smoker_ids]
expr_copd <- filtered_expression[, copd_ids]
# Healthy vs COPD expression matrix
expr_hnVScopd <- filtered_expression[,c(hn_ids,copd_ids)]
# write.csv(expr_hnVScopd,"final/expr_hnVScopd.csv")
# expr_hnVScopd <- read.csv("final/expr_hnVScopd.csv", row.names = 1)
# Smoker vs COPD expression matrix
expr_smokerVScopd <- filtered_expression[,c(smoker_ids,copd_ids)]
# write.csv(expr_smokerVScopd,"output/expr_smokerVScopd.csv")



# Subset phenotype information for the selected samples
pheno_hnVScopd <- phenotype_sorted[phenotype_sorted$GSM_IDs %in% colnames(expr_hnVScopd), ]

# Ensure phenotype vector matches the column order
pheno_hnVScopd <- pheno_hnVScopd[match(colnames(expr_hnVScopd), pheno_hnVScopd$GSM_IDs), ]
dim(pheno_hnVScopd)
dim(expr_hnVScopd)
table(pheno_hnVScopd$GSE_ID)
table(pheno_hnVScopd$phenotype)

# save the phenotype data
# write.csv(pheno_hnVScopd, "final/pheno_hnVScopd.csv", row.names = FALSE)
# pheno_hnVScopd1 <- read.csv("final/pheno_hnVScopd.csv")
#####################
# Pre-processing end
#####################


#############################
# Batch effect removal starts
#############################
# Create a sample group mapping
sample_groups <- data.frame(
  SampleID = colnames(expr_hnVScopd),
  Group = factor(pheno_hnVScopd$phenotype) #c(rep("HN", 161), rep("COPD", 92))  # Adjust based on your actual group sizes
)
# Extract gene symbols
gene_symbols <- sapply(strsplit(rownames(expr_hnVScopd), " /// "), `[`, 1)

# Identify valid (non-NA and non-duplicated) entries
valid <- !is.na(gene_symbols) & !duplicated(gene_symbols)

# Subset the expression matrix and assign new rownames
expr_hnVScopd <- expr_hnVScopd[valid, ]
rownames(expr_hnVScopd) <- gene_symbols[valid]


# log normalization
expr_hnVScopd_log <- log2(expr_hnVScopd + 1)

# Load the limma package for normalization
library(limma)
expr_log2_norm <- normalizeBetweenArrays(expr_hnVScopd_log, method = "quantile")



# Batch effect removal
library(sva)

# Batch correction
batch_corrected <- ComBat(dat = as.matrix(expr_log2_norm), batch = pheno_hnVScopd$GSE_ID, par.prior = TRUE)
# write.csv(batch_corrected, "final/batch_corrected.csv", row.names = TRUE)
head(batch_corrected[1:5, 1:5])
dim(batch_corrected)
# batch_corrected <- read.csv("final/batch_corrected.csv", row.names = 1)

# PCA Plot before normalization
library(ggplot2)
library(ggfortify)
pca_data <- t(expr_hnVScopd)  # transpose to samples x genes
pca_plot <- autoplot(prcomp(pca_data, scale. = TRUE),
         data = pheno_hnVScopd,
         colour = 'GSE_ID',
         shape = 'phenotype') +
  theme_minimal() +
  ggtitle("PCA Plot: Samples by Phenotype")
ggsave("final/final_PCA_beforeNormalization.png", plot = pca_plot, width = 8, height = 8, dpi = 300)

# PCA Plot before batch effect removal
library(ggplot2)
library(ggfortify)
pca_data <- t(expr_log2_norm)  # transpose to samples x genes
pca_plot <- autoplot(prcomp(pca_data, scale. = TRUE),
         data = pheno_hnVScopd,
         colour = 'GSE_ID',
         shape = 'phenotype') +
  theme_minimal() +
  ggtitle("PCA Plot: Samples by Phenotype")
ggsave("final/final_PCA_beforeBatchCorrectionWithlogNormalized.png", plot = pca_plot, width = 8, height = 8, dpi = 300)

# PCA Plot After batch effect removal
library(ggplot2)
library(ggfortify)
pca_data <- t(batch_corrected)  # transpose to samples x genes
pca_plot <- autoplot(prcomp(pca_data, scale. = TRUE),
         data = pheno_hnVScopd,
         colour = 'GSE_ID',
         shape = 'phenotype') +
  theme_minimal() +
  ggtitle("PCA Plot: Samples by Phenotype")
ggsave("final/final_PCA_afterBatchCorrection.png", plot = pca_plot, width = 8, height = 8, dpi = 300)



#############################
# Batch effect removal done
#############################



#############
# DGEA starts
#############

# Load required packages
library(limma)
library(ggplot2)
library(ggrepel)
library(ggfortify)

# STEP 1: Log2 transform and normalize expression data
expr_log2 <- log2(expr_hnVScopd + 1)
expr_norm <- normalizeBetweenArrays(expr_log2, method = "quantile")
# keepSamples <- pheno_hnVScopd$GSM_IDs
# batch_corrected <- batch_corrected[,keepSamples]
# STEP 2: Prepare phenotype and design matrix
group <- factor(pheno_hnVScopd$phenotype)
design <- model.matrix(~ 0 + group)
colnames(design) <- make.names(levels(group))  # syntactically valid

# STEP 3: Linear model and contrasts
fit <- lmFit(expr_norm, design)
contrast_matrix <- makeContrasts(SmokerwithCOPD_vs_Healthy = Smoker.with.COPD - Healthy.Non.Smoker, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# STEP 4: Get DEGs
deg_table <- topTable(fit2, number = Inf, adjust.method = "BH")
deg_table$Gene <- rownames(deg_table)
head(deg_table[1:5, 1:5])
length(deg_table$Gene)
# STEP 5: Volcano Plot
deg_table$Significance <- "Not Significant"
deg_table$Significance[deg_table$logFC > 1 & deg_table$adj.P.Val < 0.05] <- "Upregulated"
deg_table$Significance[deg_table$logFC < -1 & deg_table$adj.P.Val < 0.05] <- "Downregulated"
cat("Total No. of DEGs: ", sum(deg_table$Significance == "Upregulated") + sum(deg_table$Significance == "Downregulated"), "\n",
    "Upregulated: ", sum(deg_table$Significance == "Upregulated"), "\n",
    "Downregulated: ", sum(deg_table$Significance == "Downregulated"), "\n")

# write.csv(deg_table[deg_table$Significance != "Not Significant", ], "final/final_DEG_table_full.csv", row.names = FALSE)
# write.csv(deg_table[deg_table$Significance == "Upregulated", ], "final/final_DEG_table_Upregulated.csv", row.names = FALSE)
# write.csv(deg_table[deg_table$Significance == "Downregulated", ], "final/final_DEG_table_Downregulated.csv", row.names = FALSE)

png("final/final_volcano_plot.png", width = 4000, height = 4000, res = 300)
ggplot(deg_table, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "#af0c0c", "Downregulated" = "#161697", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Smoker with COPD vs Healthy Non-Smoker",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


# Load package
library(pheatmap)

# STEP 1: Select top 50 DEGs by adjusted p-value
top_degs <- deg_table[deg_table$Significance != "Not Significant", ]
top_degs <- top_degs[order(top_degs$adj.P.Val), ][1:50, ]
top_gene_ids <- top_degs$Gene

# STEP 2: Extract expression data for top DEGs
expr_top_degs <- batch_corrected[top_gene_ids, ]

# STEP 3: Prepare annotation for columns (samples)
annotation_col <- data.frame(Phenotype = pheno_hnVScopd$phenotype)
rownames(annotation_col) <- colnames(expr_top_degs)

# STEP 4: Generate heatmap
library(pheatmap)
# Save to PNG
png("final/final_Heatmap_top_50_DEGs.png", width = 4000, height = 4000, res = 300)
pheatmap(expr_top_degs,
         annotation_col = annotation_col,
         scale = "row",                    # z-score normalization by row
         ## Improved clustering settings:
         clustering_distance_rows = "euclidean",           # magnitude-based
         clustering_distance_cols = "correlation",         # retains pattern similarities
         clustering_method = "ward.D2",                    # better cluster structure
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 50 DEGs Heatmap: Smoker with COPD vs Healthy Non-Smoker")
dev.off()


# Gene Ontology (GO) and KEGG enrichment analysis
library(ggplot2)

# Barchart for UP regulated
barchart <- read.csv("final/UP/barchart.csv")
Category <- barchart$Category
Count <- barchart$Count
Term <- barchart$Term

png("final/UP/Up_barchart.png", width = 4000, height = 4000, res = 300)
ggplot(barchart,aes(Count, reorder(Term,Count), fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  facet_wrap( ~Category, nrow = 2, scales = "free") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank())
dev.off()


# Barchart for Down regulated
barchart <- read.csv("final/DOWN/barchart-down.csv")
Category <- barchart$Category
Count <- barchart$Count
Term <- barchart$Term

png("final/DOWN/Down_barchart.png", width = 4000, height = 4000, res = 300)
ggplot(barchart,aes(Count, reorder(Term,Count), fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  facet_wrap( ~Category, nrow = 2, scales = "free") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank())
dev.off()

#############
# DGEA ends
#############

############################
# Pearson correlation starts
############################

# Compute Pearson correlation between samples (Top 50 DEGs)
cor_matrix <- cor(expr_top_degs, method = "pearson")
# cor_matrix_genes <- cor(t(expr_top_degs), method = "pearson")

# Load and plot
library(pheatmap)
# Save to PNG
annotation_col <- data.frame(
  Condition = pheno_hnVScopd$phenotype
)
rownames(annotation_col) <- colnames(cor_matrix)
# Save to PNG
png("final/Annotated_correlation_heatmap.png", width = 4000, height = 4000, res = 300)
pheatmap::pheatmap(cor_matrix,
         annotation_col = annotation_col,
         main = "Correlation Heatmap (Pearson) Top 50 DEGs",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         display_numbers = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()
############################
# Pearson correlation ends
############################

####################################
# Network Centrality Analysis starts
####################################

# install.packages("igraph")
# install.packages("CINNA")

library(igraph)
library(CINNA)

# Load the PPI data
ppi <- read.delim("final/final_204DEGs_string_interactions_short.tsv", header = TRUE, sep = "\t")
colnames(ppi)

# Create an igraph object from the edge list
ppi_graph <- graph_from_data_frame(ppi[, c("X.node1", "node2")], directed = FALSE)


# Extract the largest connected component (giant component)
# `clusters()` was deprecated in igraph 2.0.0. Thus use `components()` instead
components <- components(ppi_graph)
giant_component <- induced_subgraph(ppi_graph, which(components$membership == which.max(components$csize)))


# Identify available centrality measures for the network
available_measures <- proper_centralities(giant_component)

# Calculate centralities for all valid measures
centrality_scores <- calculate_centralities(giant_component)

# Filter out zero-length centrality vectors
centrality_scores_clean_list <- centrality_scores[sapply(centrality_scores, function(x) length(x) > 0)]

# Convert remaining valid ones to data frame
centrality_df <- as.data.frame(centrality_scores_clean_list)
# write.csv(centrality_df, "final/final_centrality_scores.csv", row.names = TRUE)
# centrality_df <- read.csv("final/final_centrality_scores.csv", row.names = 1)
# Remove columns with NA values
centrality_df_clean <- centrality_df[, colSums(is.na(centrality_df)) == 0]
# write.csv(centrality_df_clean, "final/final_centrality_scores_clean.csv", row.names = TRUE)

# Load library
if (!require("factoextra")) install.packages("factoextra")
library(factoextra)

# Perform PCA
pca_res <- prcomp(centrality_df_clean, scale. = TRUE)

# Plot variable contributions
png("final/final_PCA_Variable_Contributions_with_direction.png", width = 4000, height = 4000, res = 300)
fviz_pca_var(pca_res, col.var = "contrib", repel = TRUE, axes = c(1, 2))
dev.off()





# Calculate contributions (squared loadings) of each centrality measure to PC1
contributions <- abs(pca_res$rotation[, 1])^2
contributions <- contributions / sum(contributions) * 100  # Convert to percentage

# Create a data frame for plotting
contrib_df <- data.frame(
  Centrality = names(contributions),
  Contribution = contributions
)

# Sort by contribution descending
contrib_df <- contrib_df[order(-contrib_df$Contribution), ]

# Bar chart
# Mark top 5 contributors
contrib_df$Color <- ifelse(rank(-contrib_df$Contribution) <= 5, "Top5", "Other")

# Plot with color grouping
png("final/final_PCA_Centrality_Contributions_to_PC1_barchart.png", width = 4000, height = 4000, res = 300)
ggplot(contrib_df, aes(x = reorder(Centrality, -Contribution), y = Contribution, fill = Color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Top5" = "#a90f0f", "Other" = "#276aa0")) +
  labs(x = "Centrality Measure",
       y = "Contribution (%)",
       fill = "Legend") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



# Get top 5 contributors to PC1
top5_contributors <- head(contrib_df, 5)
print(top5_contributors$Centrality)

# selection pf Hub Genes
# Step 1: Select top 5 centrality measures from the full centrality dataframe
top5_measures <- top5_contributors$Centrality
top5_scores <- centrality_df_clean[, top5_measures]

# Step 2: Scale the centrality scores
top5_scores_scaled <- scale(top5_scores)

# Step 3: Calculate composite centrality score (mean of scaled values)
composite_score <- rowMeans(top5_scores_scaled)

# Step 4: Select hub genes with composite score > mean(composite_score)
threshold <- mean(composite_score)
hub_genes <- names(composite_score[composite_score > threshold])

# Optionally: Create a data frame of hub genes and their composite scores
hub_gene_df <- data.frame(
  Gene = hub_genes,
  CompositeScore = composite_score[hub_genes],
  top5_scores = top5_scores[hub_genes, ]
)

# View or export
names(hub_gene_df)
length(hub_gene_df$Gene)
dim(hub_gene_df)
head(hub_gene_df[1:5, 1:5])
# Order by CompositeScore descending
hub_gene_df_ordered <- hub_gene_df[order(-hub_gene_df$CompositeScore), ]
# Select top 20 hub genes
top20_hub_genes <- head(hub_gene_df_ordered, 20)
names(top20_hub_genes)
length(top20_hub_genes$Gene)

# Save the hub gene data frame
# write.csv(hub_gene_df_ordered, "final/final_Hub_genes_fromCentralityMeasuresOnPPI-WithCompositeScoreAndIndividualCentralityScores.csv", row.names = FALSE)
# read.csv("final/final_Hub_genes_fromCentralityMeasuresOnPPI-WithCompositeScoreAndIndividualCentralityScores.csv", row.names = 1)

# Save the top 20 hub genes
# write.csv(top20_hub_genes, "final/final_Hub_genes_Top_20_fromCentralityMeasuresOnPPI.csv", row.names = FALSE)
# read.csv("final/final_Hub_genes_Top_20_fromCentralityMeasuresOnPPI.csv", row.names = 1)


######################################## Don't remove below commented codes

# # with only the top centrality measure
# # Step 1: Select the top 1 centrality measure
# top5_measures <- top5_contributors$Centrality
# top5_scores <- centrality_df_clean[, top5_measures]
# top1_measures <- top5_contributors$Centrality[1]  # Select the top measure
# top1_scores <- data.frame("Genes" = rownames(top5_scores),"Top_measure" = top5_scores[, 1])  # Select the first column (top measure)
# rownames(top1_scores) <- top1_scores$Genes  # Set rownames to gene names
# # Step 2: Compute threshold (mean)
# threshold <- mean(top1_scores$Top_measure)

# # Step 3: Identify hub genes with score > mean
# hub_genes <- row.names.data.frame(top1_scores[top1_scores$Top_measure > threshold,])

# # Step 4: Create data frame of hub genes and their centrality scores
# hub_gene_df <- data.frame(
#   Gene = top1_scores$Genes[top1_scores$Top_measure > threshold],
#   CentralityScore = top1_scores$Top_measure[top1_scores$Top_measure > threshold]
# )

# # View first few hub genes and their count
# head(hub_gene_df)
# length(hub_gene_df$Gene)


# hub_gene_df_compositeScore <- read.csv("final/final_Hub_genes_composite_score_fromCentralityMeasuresOnPPI.csv", row.names = 1)
# head(hub_gene_df_compositeScore)
# length(hub_gene_df_compositeScore$Gene)

######################################## Don't remove above commented codes



# # Get top contributing measure
# var_contrib <- get_pca_var(pca_res)$contrib
# mean_contrib <- rowMeans(var_contrib[, 1:2])
# top_5_measure <- names(sort(mean_contrib, decreasing = TRUE)[1:5])
# cat("Top 5 contributing centrality:", top_5_measure, "\n")
# top_measure <- names(sort(mean_contrib, decreasing = TRUE)[1])

# # Extract and filter hub genes (HUB genes not identified in this step)
# # Print column names
# colnames(centrality_df_clean)

# # Print the top_measure value
# print(top_measure)





# # Top 20 genes ranked from centrality measures were saved as CSV files in the temp directory
# # these files are generated by cytoscape using the cytohubban plugin and starts with "20"
# # Merge top 20 genes from each centrality measure
# top20 <- list.files("temp/", pattern = ".csv", full.names = T)
# top20 <- top20[1:5]

# merged_data <- lapply(top20, function(file) {
#   df <- read.csv(file, header = TRUE)
#   df
# })

# final_merged_data <- Reduce(function(x, y) merge(x, y, by = "Name", all = TRUE), merged_data)
# final_merged_data <- final_merged_data[!is.na(final_merged_data$Rank.x),]
# final_merged_data <- final_merged_data[!is.na(final_merged_data$Score),]
# dim(final_merged_data)
# # Save the final merged data to a CSV file
# write.csv(final_merged_data$Name, "temp/listofHubGenes.csv", row.names = FALSE, quote = FALSE)

##################################
# Network Centrality Analysis Ends
##################################


##############
# LASSO starts
##############
# Install if not already installed
# packages <- c("glmnet", "randomForest", "caret", "doParallel")
# installed <- packages %in% rownames(installed.packages())
# if (any(!installed)) install.packages(packages[!installed])

# Load required libraries
library(glmnet)
library(randomForest)
library(caret)
library(doParallel)

# # STEP 1: Select top 50 DEGs by adjusted p-value
# top_degs <- deg_table[deg_table$Significance != "Not Significant", ]
# top_gene_ids <- top_degs$Gene

# STEP 1: Select the hub genes from the centrality analysis
# top20_hub_genes <- read.csv("final/final_Hub_genes_Top_20_fromCentralityMeasuresOnPPI.csv")
hub_genes <- top20_hub_genes$Gene

# STEP 2: Extract expression data for top DEGs
expr_top_hub <- batch_corrected[hub_genes, ]

# Transpose expression data for glmnet (samples x genes)
x <- t(expr_top_hub)
y <- factor(pheno_hnVScopd$phenotype)

# Binary classification: 0 = Healthy Non-Smoker, 1 = Smoker with COPD
y_bin <- ifelse(y == "Smoker with COPD", 1, 0)

# Standardize features
# x_scaled <- scale(x)

# Create cross-validation folds
set.seed(123)
cv_lasso <- cv.glmnet(x, y_bin, alpha = 1, family = "binomial", nfolds = 10) # x_scaled not used here

# Plot CV results
# Plot with lambda.min and lambda.1se indicated
png("final/final_LASSO_CV_plot.png", width = 4000, height = 4000, res = 300)
plot(cv_lasso)

# Add vertical lines for lambda.min and lambda.1se
abline(v = log(cv_lasso$lambda.min), col = "red", lty = 2, lwd = 2)
abline(v = log(cv_lasso$lambda.1se), col = "blue", lty = 2, lwd = 2)

# Add text labels
legend("bottomright",
       legend = c(paste0("lambda.min = ", signif(cv_lasso$lambda.min, 1)),
                  paste0("lambda.1se = ", signif(cv_lasso$lambda.1se, 1))),
       col = c("red", "blue"),
       lty = 2, lwd = 2, bty = "n")
dev.off()

# Fit LASSO model (alpha = 1 means LASSO)
lasso_fit <- glmnet(x, y_bin, alpha = 1, family = "binomial")

# Plot coefficient profiles
png("final/final_LASSO_Coefficient_Profiles.png", width = 4000, height = 4000, res = 300)
plot(lasso_fit, xvar = "lambda", label = TRUE, lwd = 1.5)
dev.off()


# Best lambda
best_lambda <- cv_lasso$lambda.min
cat("Best lambda (min):", best_lambda, "\n")
# Best lambda (1se)
best_lambda_1se <- cv_lasso$lambda.1se
cat("Best lambda (1se):", best_lambda_1se, "\n")

# Coefficients at best lambda
lasso_coef <- coef(cv_lasso, s = best_lambda)
selected_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
selected_genes <- selected_genes[!selected_genes %in% "(Intercept)"]

# Expression matrix with selected genes
lasso_expr <- x[, selected_genes]
head(lasso_expr[1:5, 1:5])

# write.csv(lasso_expr, "final/final_lasso_selected_genes_expression.csv", row.names = TRUE)
# lasso_expr <- read.csv("final/final_lasso_selected_genes_expression.csv", row.names = 1)


# ================================ preprocessing for the network analysis with 13 selected genes
# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)

# Read the CSV file
# Replace "your_file.csv" with your actual file path
data <- read_csv("final/Experimental_functional_annotation_table_on_13_genes.csv")


# Transform: split the Genes column and unnest rows
transformed_data <- data %>%
  separate_rows(Genes, sep = ",\\s*") %>%   # expand rows
  mutate(Genes = trimws(Genes)) %>%         # clean whitespace
  select(Genes, Term, Category)             # reorder columns

# View the transformed data in R console
print(transformed_data)

# Write the transformed data to a new CSV file
write_csv(transformed_data, "final/Experimental_transformed_table_for_network_analysis.csv")



lasso_expr <- read.csv("final/final_lasso_selected_genes_expression.csv", row.names = 1)
cat(names(lasso_expr))
full_DEG_table <- read.csv("final/final_DEG_table_full.csv")
selected <- subset(full_DEG_table, Gene %in% c("CCL2","SPP1","CD163","LTF","MMP12", 
                                               "TREM2","PPARG","KRT5","SCGB1A1", 
                                               "TEKT1","DNAI1","FCN1","LPL"))

selected <- subset(selected, c("logFC","Gene", "Significance"))
write.csv(selected, "final/logFC_of_13_selected_genes.csv", row.names = F)


read.csv("final/Experimental_transformed_table_for_network_analysis.csv") -> transformed_table


merged_transformed_table <- merge(transformed_table, selected, by.x = "Genes", by.y = "Gene")
View(merged_transformed_table)
write.csv(merged_transformed_table, "final/Experimental_merged_transformed_table_for_network_analysis.csv",row.names = F)

# ================================ preprocessing for the network analysis with 13 selected genes



# ================================ ML medels starts here
###############
# Random Forest
###############

# Convert phenotype to factor with two levels
rf_df <- data.frame(lasso_expr, phenotype = factor(y))

# Set up repeated RF modeling
set.seed(123)
n_vars <- seq_along(selected_genes) # 1:length(selected_genes)
error_rates <- numeric(length(n_vars))

for (i in n_vars) {
  rf_model <- randomForest(phenotype ~ ., 
                           data = rf_df[, c(1:i, ncol(rf_df))],
                           ntree = 500,
                           importance = TRUE)
  error_rates[i] <- 1 - rf_model$confusion["Smoker with COPD", "Smoker with COPD"] / sum(rf_model$confusion["Smoker with COPD", ])
}

# Plot error rate vs number of variables
png("final/final_error_rate_vs_genes1.png", width = 1200, height = 1200, res = 150)
plot(n_vars, error_rates, type = "b", col = "darkgreen", 
     xlab = "Number of Top LASSO-selected Genes",
     ylab = "Random Forest Error Rate",
     main = "Error Rate vs. Number of Genes")
dev.off()
# Optional: Identify minimum error configuration
min_error_index <- which.min(error_rates)
cat("Minimum error rate:", error_rates[min_error_index], "using", min_error_index, "genes.\n")


# Train final model with optimal number of genes
best_genes <- selected_genes[1:min_error_index]
final_rf <- randomForest(phenotype ~ ., 
                         data = rf_df[, c(best_genes, "phenotype")],
                         ntree = 500,
                         importance = TRUE)

# Plot variable importance
png("final/final_RF_Variable_Importance_forBest9Genes.png", width = 1200, height = 1200, res = 150)
varImpPlot(final_rf, type = 1, main = "Top Gene Importance (MeanDecreaseAccuracy)")
dev.off()

# # Train final model with all 13 selected genes
# final_rf <- randomForest(phenotype ~ ., 
#                          data = rf_df[, c(selected_genes, "phenotype")],
#                          ntree = 500,
#                          importance = TRUE)

# # Plot variable importance
# png("final/final_RF_Variable_Importance_forSelected13Genes.png", width = 1200, height = 1200, res = 150)
# varImpPlot(final_rf, type = 1, main = "Top Gene Importance (MeanDecreaseAccuracy)")
# dev.off()


library(caret)

ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
rf_cv <- train(phenotype ~ ., data = rf_df[, c(best_genes, "phenotype")],
               method = "rf", trControl = ctrl)
print(rf_cv)



################################################### performancs evaluation on only training data
# This is not not presentable as the model is trained on the same data
# but it is useful to check the performance of the model on the training data
# Load required libraries

library(caret)
library(pROC)
library(e1071)

# Predict on the training data using the final model
rf_pred <- predict(final_rf, rf_df[, best_genes])
rf_prob <- predict(final_rf, rf_df[, best_genes], type = "prob")

# Create confusion matrix
conf_mat <- confusionMatrix(rf_pred, rf_df$phenotype, positive = "Smoker with COPD")
print(conf_mat)

# Extract metrics
accuracy     <- conf_mat$overall["Accuracy"]
sensitivity  <- conf_mat$byClass["Sensitivity"]
specificity  <- conf_mat$byClass["Specificity"]
f1_score     <- conf_mat$byClass["F1"]

# Print metrics
cat("Accuracy   :", round(accuracy, 4), "\n")
cat("Sensitivity:", round(sensitivity, 4), "\n")
cat("Specificity:", round(specificity, 4), "\n")
cat("F1 Score   :", round(f1_score, 4), "\n")

# Compute AUROC (needs numeric predictions)
# Create binary vector for true class
y_true_bin <- ifelse(rf_df$phenotype == "Smoker with COPD", 1, 0)

# Probabilities for the positive class
rf_prob_copd <- rf_prob[, "Smoker with COPD"]

# AUROC
roc_obj <- roc(y_true_bin, rf_prob_copd)
auc_score <- auc(roc_obj)
cat("AUROC      :", round(auc_score, 4), "\n")

# Plot ROC curve
png("final/final_AUROC_RF1.png", width = 1200, height = 1200, res = 150)
plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve - Random Forest", legacy.axes = FALSE)
abline(a = 0, b = 1, lty = 2, col = "gray")
text(0.6, 0.2, labels = paste0("AUC = ", round(auc_score, 4)), col = "black", cex = 1.5)
dev.off()



###############
# SVM Model
###############

# Prepare data for SVM
svm_df <- rf_df[, c(best_genes, "phenotype")]
svm_df$phenotype <- factor(svm_df$phenotype, levels = c("Healthy Non-Smoker", "Smoker with COPD"))

# Train SVM model with radial kernel
svm_model <- svm(phenotype ~ ., data = svm_df, kernel = "radial", probability = TRUE)

# Predict on training data
svm_pred <- predict(svm_model, svm_df[, best_genes])
svm_prob <- attr(predict(svm_model, svm_df[, best_genes], probability = TRUE), "probabilities")

# Confusion Matrix
svm_conf <- confusionMatrix(svm_pred, svm_df$phenotype, positive = "Smoker with COPD")
print(svm_conf)

# Extract metrics
svm_accuracy     <- svm_conf$overall["Accuracy"]
svm_sensitivity  <- svm_conf$byClass["Sensitivity"]
svm_specificity  <- svm_conf$byClass["Specificity"]
svm_f1_score     <- svm_conf$byClass["F1"]

# Print metrics
cat("\n=== SVM Performance ===\n")
cat("Accuracy   :", round(svm_accuracy, 4), "\n")
cat("Sensitivity:", round(svm_sensitivity, 4), "\n")
cat("Specificity:", round(svm_specificity, 4), "\n")
cat("F1 Score   :", round(svm_f1_score, 4), "\n")

# AUROC Calculation
svm_prob_copd <- svm_prob[, "Smoker with COPD"]
svm_roc_obj <- roc(y_true_bin, svm_prob_copd)
svm_auc <- auc(svm_roc_obj)
cat("AUROC      :", round(svm_auc, 4), "\n")

# Plot ROC curve
png("final/final_AUROC_SVM1.png", width = 1200, height = 1200, res = 150)
plot(svm_roc_obj, col = "red", lwd = 2, main = "ROC Curve - SVM")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()


# # List of required packages
# required_packages <- c("caret", "pROC", "e1071", "glmnet", "xgboost", 
#                        "ggplot2", "dplyr", "tibble", "forcats", "tidyr")

# # Install missing packages
# new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
# if(length(new_packages)) install.packages(new_packages)

# # Load all required packages
# lapply(required_packages, require, character.only = TRUE)


###############
# Logistic Regression
###############
library(caret)
library(pROC)

# Fit logistic regression model
glm_model <- glm(phenotype ~ ., data = rf_df[, c(best_genes, "phenotype")], family = binomial)

# Predict class probabilities and labels
glm_prob <- predict(glm_model, rf_df[, best_genes], type = "response")
glm_pred <- factor(ifelse(glm_prob > 0.5, "Smoker with COPD", "Healthy Non-Smoker"), 
                   levels = levels(rf_df$phenotype))

# Evaluate
conf_glm <- confusionMatrix(glm_pred, rf_df$phenotype, positive = "Smoker with COPD")
print(conf_glm)

# Extract performance
acc_glm <- conf_glm$overall["Accuracy"]
sens_glm <- conf_glm$byClass["Sensitivity"]
spec_glm <- conf_glm$byClass["Specificity"]
f1_glm   <- conf_glm$byClass["F1"]

# AUROC
glm_roc <- roc(y_true_bin, glm_prob)
auc_glm <- auc(glm_roc)

# Output
cat("Logistic Regression Results:\n")
cat("Accuracy   :", round(acc_glm, 4), "\n")
cat("Sensitivity:", round(sens_glm, 4), "\n")
cat("Specificity:", round(spec_glm, 4), "\n")
cat("F1 Score   :", round(f1_glm, 4), "\n")
cat("AUROC      :", round(auc_glm, 4), "\n")

# Plot
png("final/final_AUROC_GLM1.png", width = 1200, height = 1200, res = 150)
plot(glm_roc, col = "green", lwd = 2, main = "ROC Curve - Logistic Regression")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()



###############
# XGBoost
###############
# Load necessary libraries
library(xgboost)
library(caret)
library(pROC)

# Prepare data
xgb_df <- rf_df[, c(best_genes, "phenotype")]
xgb_df$phenotype <- factor(ifelse(xgb_df$phenotype == "Smoker with COPD", "COPD", "Healthy"))

# Ensure phenotype has two levels: "Healthy", "COPD" (required for twoClassSummary)
xgb_df$phenotype <- factor(xgb_df$phenotype, levels = c("Healthy", "COPD"))

# Set up trainControl with ROC summary
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = "final")

# Train XGBoost using caret
set.seed(123)
xgb_model <- train(phenotype ~ ., data = xgb_df, method = "xgbTree",
                   trControl = ctrl, metric = "ROC")

# Predictions
xgb_pred <- predict(xgb_model, xgb_df[, best_genes])
xgb_prob <- predict(xgb_model, xgb_df[, best_genes], type = "prob")[, "COPD"]

# Confusion Matrix
conf_xgb <- confusionMatrix(xgb_pred, xgb_df$phenotype, positive = "COPD")
print(conf_xgb)

# Extract metrics
acc_xgb <- conf_xgb$overall["Accuracy"]
sens_xgb <- conf_xgb$byClass["Sensitivity"]
spec_xgb <- conf_xgb$byClass["Specificity"]
f1_xgb   <- conf_xgb$byClass["F1"]

# AUROC
xgb_roc <- roc(response = ifelse(xgb_df$phenotype == "COPD", 1, 0), predictor = xgb_prob)
auc_xgb <- auc(xgb_roc)

# Print performance
cat("XGBoost Results:\n")
cat("Accuracy   :", round(acc_xgb, 4), "\n")
cat("Sensitivity:", round(sens_xgb, 4), "\n")
cat("Specificity:", round(spec_xgb, 4), "\n")
cat("F1 Score   :", round(f1_xgb, 4), "\n")
cat("AUROC      :", round(auc_xgb, 4), "\n")

# Save ROC Plot
png("final/final_AUROC_XGB1.png", width = 1200, height = 1200, res = 150)
plot(xgb_roc, col = "darkorange", lwd = 2, main = "ROC Curve - XGBoost")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()
# ================================ ML medels ends here


###############
# Dependencies
###############
# packages <- c("caret", "pROC", "e1071", "glmnet", "xgboost", "ggplot2", "dplyr", "tibble", "forcats", "tidyr")
# installed <- packages %in% rownames(installed.packages())
# if (any(!installed)) install.packages(packages[!installed])

# Load libraries
library(caret)
library(pROC)
library(e1071)
library(glmnet)
library(xgboost)
library(ggplot2)
library(dplyr)
library(tibble)
library(forcats)
library(tidyr)

#####################
# Input Preprocessing
#####################
lasso_expr <- read.csv("final/final_lasso_selected_genes_expression.csv", row.names = 1)
pheno_hnVScopd <- read.csv("final/pheno_hnVScopd.csv")
y <- factor(pheno_hnVScopd$phenotype)
# make y binary
y_bin <- ifelse(y == "Smoker with COPD", 1, 0)
# Convert phenotype to factor with two levels
rf_df <- data.frame(lasso_expr, phenotype = factor(y))



x_data <- rf_df[, 1:13]  # Select only the gene expression columns instead of best_genes
y_data <- rf_df$phenotype

train_index <- createDataPartition(y_bin, p = 0.7, list = FALSE)
x_train <- lasso_expr[train_index, ]
x_test  <- lasso_expr[-train_index, ]
y_train <- y[train_index]
y_test  <- y[-train_index]

y_train <- factor(y_train)
levels(y_train) <- make.names(levels(y_train))

y_test <- factor(y_test)
levels(y_test) <- make.names(levels(y_test))
########################
# Train ML Models
########################

ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5,
                     classProbs = TRUE, summaryFunction = twoClassSummary,
                     savePredictions = "final")

models <- list()

set.seed(123)
models$RF <- train(x = x_train, y = y_train, method = "rf", trControl = ctrl, metric = "ROC")

set.seed(123)
models$SVM <- train(x = x_train, y = y_train, method = "svmRadial", trControl = ctrl, metric = "ROC")

set.seed(123)
models$LogReg <- train(x = x_train, y = y_train, method = "glm", family = "binomial", trControl = ctrl, metric = "ROC")

set.seed(123)
models$XGB <- train(x = x_train, y = y_train, method = "xgbTree", trControl = ctrl, metric = "ROC")

########################
# Collect Performance Metrics
########################

results <- data.frame(Model = character(), Accuracy = numeric(), Sensitivity = numeric(),
                      Specificity = numeric(), F1 = numeric(), AUROC = numeric())

for (model_name in names(models)) {
  model <- models[[model_name]]
  preds <- model$pred

  # Filter predictions by best tune parameters if available
  if (!is.null(model$bestTune)) {
    for (param in names(model$bestTune)) {
      preds <- preds[preds[[param]] == model$bestTune[[param]], ]
    }
  }

  pred_labels <- preds$pred
  true_labels <- preds$obs
  prob_copd <- preds$Smoker.with.COPD

  cm <- confusionMatrix(pred_labels, true_labels, positive = "Smoker.with.COPD")

  if (length(unique(true_labels)) == 2) {
    roc_obj <- roc(response = true_labels, predictor = prob_copd,
                   levels = c("Healthy.Non.Smoker", "Smoker.with.COPD"), direction = "<")
    auc_val <- auc(roc_obj)
  } else {
    auc_val <- NA
  }

  results <- rbind(results, data.frame(
    Model = model_name,
    Accuracy = cm$overall["Accuracy"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    F1 = cm$byClass["F1"],
    AUROC = as.numeric(auc_val)
  ))
}
print(results)

# ============================ Test the models on the test set
########################
# Evaluate on Test Set
########################

results_test <- data.frame(Model = character(), Accuracy = numeric(), Sensitivity = numeric(),
                           Specificity = numeric(), F1 = numeric(), AUROC = numeric())

for (model_name in names(models)) {
  model <- models[[model_name]]
  
  # Predict class labels
  pred_labels <- predict(model, x_test)
  
  # Predict probabilities
  prob_preds <- predict(model, x_test, type = "prob")
  
  # Probabilities for positive class ("Smoker.with.COPD")
  prob_copd <- prob_preds$Smoker.with.COPD
  
  # Confusion matrix on test set
  cm <- confusionMatrix(pred_labels, y_test, positive = "Smoker.with.COPD")
  
  # ROC/AUC
  if (length(unique(y_test)) == 2) {
    roc_obj <- roc(response = y_test, predictor = prob_copd,
                   levels = c("Healthy.Non.Smoker", "Smoker.with.COPD"), direction = "<")
    auc_val <- auc(roc_obj)
  } else {
    auc_val <- NA
  }
  
  # Store metrics
  results_test <- rbind(results_test, data.frame(
    Model = model_name,
    Accuracy = cm$overall["Accuracy"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    F1 = cm$byClass["F1"],
    AUROC = as.numeric(auc_val)
  ))
}

cat("\n===== TEST SET PERFORMANCE =====\n")
print(results_test)


########################
# Plot Performance Metrics
########################
results$Set <- "CV"
results_test$Set <- "Test"

combined_results <- bind_rows(results, results_test)
combined_results_long <- combined_results %>%
  pivot_longer(cols = c(Accuracy, Sensitivity, Specificity, F1, AUROC),
               names_to = "Metric", values_to = "Score")

p_combined <- ggplot(combined_results_long, aes(x = fct_reorder(Model, Score), y = Score, fill = Set)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Metric, scales = "free_y", labeller = labeller(Metric = label_value)) +
  theme_minimal(base_size = 14) +
  labs(title = "Performance Comparison: Cross-Validation vs Test Set",
       x = "Model", y = "Score", fill = "Dataset") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 14))

print(p_combined)

ggsave("final/final_integrated_model_comparison_with_train_and_test.png", plot = p_combined, width = 12, height = 12, dpi = 300)
################################################### performancs evaluation on only training data



################################################### performancs evaluation on both training & test data
# This is presentable as the model is trained on the training data and evaluated on the test data
# this generate AUROC plots for each model on both training and test data
# ================================ Load All Required Libraries

# Machine Learning models
library(randomForest)  # For Random Forest
library(e1071)         # For SVM
library(glmnet)        # For LASSO
library(xgboost)       # For XGBoost

# Model training and cross-validation
library(caret)
library(doParallel)

# ROC and AUC
library(pROC)

# Data manipulation and visualization
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(forcats)

# ================================ Load the dataset
lasso_expr <- read.csv("final/final_lasso_selected_genes_expression.csv", row.names = 1)
pheno_hnVScopd <- read.csv("final/pheno_hnVScopd.csv")

# ================================ Pre-process the data
# Ensure the phenotype is a factor
y <- factor(pheno_hnVScopd$phenotype)
# make y binary
y_bin <- ifelse(y == "Smoker with COPD", 1, 0)
set.seed(42)  # for reproducibility
# Prepare the feature matrix
train_index <- createDataPartition(y_bin, p = 0.7, list = FALSE)
x_train <- lasso_expr[train_index, ]
x_test  <- lasso_expr[-train_index, ]
y_train <- y[train_index]
y_test  <- y[-train_index]




# RF model
# Prepare training and test data
rf_train_df <- data.frame(x_train, phenotype = y_train)
rf_test_df  <- data.frame(x_test, phenotype = y_test)

# Train RF model
set.seed(123)
rf_model <- randomForest(phenotype ~ ., data = rf_train_df, ntree = 500, importance = TRUE)

# Predictions
rf_train_pred <- predict(rf_model, rf_train_df[, -ncol(rf_train_df)], type = "prob")[, "Smoker with COPD"]
rf_test_pred  <- predict(rf_model, rf_test_df[, -ncol(rf_test_df)], type = "prob")[, "Smoker with COPD"]

# ================================
# TRAIN SET CONFUSION MATRIX
# Convert probabilities to class labels using 0.5 threshold
rf_train_class <- ifelse(rf_train_pred >= 0.5, "Smoker with COPD", "Healthy Never Smoker")
rf_train_class <- factor(rf_train_class, levels = levels(y_train))

# Create confusion matrix
conf_matrix_train <- confusionMatrix(rf_train_class, y_train)
print("Confusion Matrix - Training Set")
print(conf_matrix_train)

# ================================
# TEST SET CONFUSION MATRIX
rf_test_class <- ifelse(rf_test_pred >= 0.5, "Smoker with COPD", "Healthy Never Smoker")
rf_test_class <- factor(rf_test_class, levels = levels(y_test))

conf_matrix_test <- confusionMatrix(rf_test_class, y_test)
print("Confusion Matrix - Test Set")
print(conf_matrix_test)

# ROC
rf_train_roc <- roc(ifelse(y_train == "Smoker with COPD", 1, 0), rf_train_pred)
rf_test_roc  <- roc(ifelse(y_test == "Smoker with COPD", 1, 0), rf_test_pred)

# Plot
png("temp/temp_AUROC_RF_Train_Test.png", width = 1200, height = 1200, res = 150)
plot(rf_train_roc, col = "darkblue", lwd = 2, main = "Random Forest ROC - Train vs Test")
plot(rf_test_roc, col = "#000000", lwd = 2, add = TRUE)
abline(a = 0, b = 1, col = "gray", lty = 2)
legend("bottomright", legend = c(
  paste0("Train AUC = ", round(auc(rf_train_roc), 4)),
  paste0("Test AUC = ", round(auc(rf_test_roc), 4))),
  col = c("darkblue", "#000000"), lwd = 2)
dev.off()


# SVM model
library(e1071)
svm_train_df <- data.frame(x_train, phenotype = y_train)
svm_test_df  <- data.frame(x_test, phenotype = y_test)

# Train SVM
svm_model <- svm(phenotype ~ ., data = svm_train_df, kernel = "radial", probability = TRUE)

# Predict
svm_train_pred <- predict(svm_model, svm_train_df[, -ncol(svm_train_df)])
svm_test_pred  <- predict(svm_model, svm_test_df[, -ncol(svm_test_df)])
svm_train_prob <- attr(predict(svm_model, svm_train_df[, -ncol(svm_train_df)], probability = TRUE), "probabilities")[, "Smoker with COPD"]
svm_test_prob  <- attr(predict(svm_model, svm_test_df[, -ncol(svm_test_df)], probability = TRUE), "probabilities")[, "Smoker with COPD"]


# ================================
# TRAIN SET CONFUSION MATRIX - SVM
svm_train_class <- ifelse(svm_train_prob >= 0.5, "Smoker with COPD", "Healthy Never Smoker")
svm_train_class <- factor(svm_train_class, levels = levels(y_train))

conf_matrix_svm_train <- confusionMatrix(svm_train_class, y_train)
print("Confusion Matrix - SVM Training Set")
print(conf_matrix_svm_train)

# ================================
# TEST SET CONFUSION MATRIX - SVM
svm_test_class <- ifelse(svm_test_prob >= 0.5, "Smoker with COPD", "Healthy Never Smoker")
svm_test_class <- factor(svm_test_class, levels = levels(y_test))

conf_matrix_svm_test <- confusionMatrix(svm_test_class, y_test)
print("Confusion Matrix - SVM Test Set")
print(conf_matrix_svm_test)


# ROC
svm_train_roc <- roc(ifelse(y_train == "Smoker with COPD", 1, 0), svm_train_prob)
svm_test_roc  <- roc(ifelse(y_test == "Smoker with COPD", 1, 0), svm_test_prob)

# Plot
png("temp/temp_AUROC_SVM_Train_Test.png", width = 1200, height = 1200, res = 150)
plot(svm_train_roc, col = "red", lwd = 2, main = "SVM ROC Curve - Train vs Test")
plot(svm_test_roc, col = "#000000", lwd = 2, add = TRUE)
abline(a = 0, b = 1, col = "gray", lty = 2)
legend("bottomright", legend = c(
  paste0("Train AUC = ", round(auc(svm_train_roc), 4)),
  paste0("Test AUC = ", round(auc(svm_test_roc), 4))),
  col = c("red", "#000000"), lwd = 2)
dev.off()



# LR model
glm_train_df <- data.frame(x_train, phenotype = y_train)
glm_test_df  <- data.frame(x_test, phenotype = y_test)

glm_model <- glm(phenotype ~ ., data = glm_train_df, family = "binomial")

# Probabilities
glm_train_prob <- predict(glm_model, glm_train_df[, -ncol(glm_train_df)], type = "response")
glm_test_prob  <- predict(glm_model, glm_test_df[, -ncol(glm_test_df)], type = "response")


# ================================
# TRAIN SET CONFUSION MATRIX - LR
glm_train_class <- ifelse(glm_train_prob >= 0.5, "Smoker with COPD", "Healthy Never Smoker")
glm_train_class <- factor(glm_train_class, levels = levels(y_train))

conf_matrix_glm_train <- confusionMatrix(glm_train_class, y_train)
print("Confusion Matrix - LR Training Set")
print(conf_matrix_glm_train)

# ================================
# TEST SET CONFUSION MATRIX - LR
glm_test_class <- ifelse(glm_test_prob >= 0.5, "Smoker with COPD", "Healthy Never Smoker")
glm_test_class <- factor(glm_test_class, levels = levels(y_test))

conf_matrix_glm_test <- confusionMatrix(glm_test_class, y_test)
print("Confusion Matrix - LR Test Set")
print(conf_matrix_glm_test)



# ROC
glm_train_roc <- roc(ifelse(y_train == "Smoker with COPD", 1, 0), glm_train_prob)
glm_test_roc  <- roc(ifelse(y_test == "Smoker with COPD", 1, 0), glm_test_prob)

# Plot
png("temp/temp_AUROC_GLM_Train_Test.png", width = 1200, height = 1200, res = 150)
plot(glm_train_roc, col = "darkgreen", lwd = 2, main = "Logistic Regression ROC - Train vs Test")
plot(glm_test_roc, col = "#000000", lwd = 2, add = TRUE)
abline(a = 0, b = 1, col = "gray", lty = 2)
legend("bottomright", legend = c(
  paste0("Train AUC = ", round(auc(glm_train_roc), 4)),
  paste0("Test AUC = ", round(auc(glm_test_roc), 4))),
  col = c("darkgreen", "#000000"), lwd = 2)
dev.off()



# XGBoost model
library(xgboost)
library(caret)

xgb_train_df <- data.frame(x_train, phenotype = y_train)
xgb_test_df  <- data.frame(x_test, phenotype = y_test)

xgb_train_df$phenotype <- factor(ifelse(xgb_train_df$phenotype == "Smoker with COPD", "COPD", "Healthy"))
xgb_test_df$phenotype  <- factor(ifelse(xgb_test_df$phenotype == "Smoker with COPD", "COPD", "Healthy"))

# Ensure correct factor levels
xgb_train_df$phenotype <- factor(xgb_train_df$phenotype, levels = c("Healthy", "COPD"))
xgb_test_df$phenotype  <- factor(xgb_test_df$phenotype, levels = c("Healthy", "COPD"))

ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5,
                     classProbs = TRUE, summaryFunction = twoClassSummary)

# Train model
xgb_model <- train(phenotype ~ ., data = xgb_train_df, method = "xgbTree",
                   trControl = ctrl, metric = "ROC")

# Probabilities
xgb_train_prob <- predict(xgb_model, xgb_train_df[, -ncol(xgb_train_df)], type = "prob")[, "COPD"]
xgb_test_prob  <- predict(xgb_model, xgb_test_df[, -ncol(xgb_test_df)], type = "prob")[, "COPD"]


# ================================
# TRAIN SET CONFUSION MATRIX - XGBoost
xgb_train_class <- ifelse(xgb_train_prob >= 0.5, "COPD", "Healthy")
xgb_train_class <- factor(xgb_train_class, levels = levels(xgb_train_df$phenotype))

conf_matrix_xgb_train <- confusionMatrix(xgb_train_class, xgb_train_df$phenotype)
print("Confusion Matrix - XGBoost Training Set")
print(conf_matrix_xgb_train)

# ================================
# TEST SET CONFUSION MATRIX - XGBoost
xgb_test_class <- ifelse(xgb_test_prob >= 0.5, "COPD", "Healthy")
xgb_test_class <- factor(xgb_test_class, levels = levels(xgb_test_df$phenotype))

conf_matrix_xgb_test <- confusionMatrix(xgb_test_class, xgb_test_df$phenotype)
print("Confusion Matrix - XGBoost Test Set")
print(conf_matrix_xgb_test)



# ROC
xgb_train_roc <- roc(response = ifelse(xgb_train_df$phenotype == "COPD", 1, 0), predictor = xgb_train_prob)
xgb_test_roc  <- roc(response = ifelse(xgb_test_df$phenotype == "COPD", 1, 0), predictor = xgb_test_prob)

# Plot
png("temp/temp_AUROC_XGB_Train_Test.png", width = 1200, height = 1200, res = 150)
plot(xgb_train_roc, col = "orange", lwd = 2, main = "XGBoost ROC - Train vs Test")
plot(xgb_test_roc, col = "#000000", lwd = 2, add = TRUE)
abline(a = 0, b = 1, col = "gray", lty = 2)
legend("bottomright", legend = c(
  paste0("Train AUC = ", round(auc(xgb_train_roc), 4)),
  paste0("Test AUC = ", round(auc(xgb_test_roc), 4))),
  col = c("orange", "#000000"), lwd = 2)
dev.off()
###################################
# Combined ROC plots
# Plot all test ROC curves
png("temp/temp_Combined_Test_AUROC_AllModels.png", width = 1200, height = 1200, res = 150)
plot(glm_test_roc, col = "green", lwd = 2, main = "Test ROC Curves - All Models")
plot(svm_test_roc, col = "red", lwd = 2, add = TRUE)
plot(xgb_test_roc, col = "darkorange", lwd = 2, add = TRUE)
plot(rf_test_roc, col = "blue", lwd = 2, add = TRUE)

abline(a = 0, b = 1, lty = 2, col = "gray")

legend("bottomright", legend = c(
  paste("Logistic Regression (AUC =", round(auc(glm_test_roc), 3), ")"),
  paste("SVM (AUC =", round(auc(svm_test_roc), 3), ")"),
  paste("XGBoost (AUC =", round(auc(xgb_test_roc), 3), ")"),
  paste("Random Forest (AUC =", round(auc(rf_test_roc), 3), ")")
), col = c("green", "red", "darkorange", "blue"), lwd = 2)
dev.off()

# Plot all train ROC curves
png("temp/temp_Combined_Train_AUROC_AllModels.png", width = 1200, height = 1200, res = 150)
plot(glm_train_roc, col = "green", lwd = 2, main = "Train ROC Curves - All Models")
plot(svm_train_roc, col = "red", lwd = 2, add = TRUE)
plot(xgb_train_roc, col = "darkorange", lwd = 2, add = TRUE)
plot(rf_train_roc, col = "blue", lwd = 2, add = TRUE)

abline(a = 0, b = 1, lty = 2, col = "gray")

legend("bottomright", legend = c(
  paste("Logistic Regression (AUC =", round(auc(glm_train_roc), 3), ")"),
  paste("SVM (AUC =", round(auc(svm_train_roc), 3), ")"),
  paste("XGBoost (AUC =", round(auc(xgb_train_roc), 3), ")"),
  paste("Random Forest (AUC =", round(auc(rf_train_roc), 3), ")")
), col = c("green", "red", "darkorange", "blue"), lwd = 2)
dev.off()






# ================================ Performance Comparison on all the metrics of all ML Models
get_metrics <- function(conf_mat, model_name, set_name) {
  # Helper: safely extract metric or NA if missing
  safe_extract <- function(x, name) if (name %in% names(x)) x[name] else NA
  
  data.frame(
    Model = model_name,
    Set = set_name,
    Accuracy = safe_extract(conf_mat$overall, "Accuracy"),
    BalancedAccuracy = safe_extract(conf_mat$byClass, "Balanced Accuracy"),
    Sensitivity = safe_extract(conf_mat$byClass, "Sensitivity"),
    Specificity = safe_extract(conf_mat$byClass, "Specificity"),
    Precision = safe_extract(conf_mat$byClass, "Pos Pred Value"),
    F1 = safe_extract(conf_mat$byClass, "F1")
  )
}

# Combine all into one big data frame
metrics_df <- rbind(
  get_metrics(conf_matrix_train, "RF", "Train"),
  get_metrics(conf_matrix_test, "RF", "Test"),
  get_metrics(conf_matrix_svm_train, "SVM", "Train"),
  get_metrics(conf_matrix_svm_test, "SVM", "Test"),
  get_metrics(conf_matrix_glm_train, "LR", "Train"),
  get_metrics(conf_matrix_glm_test, "LR", "Test"),
  get_metrics(conf_matrix_xgb_train, "XGB", "Train"),
  get_metrics(conf_matrix_xgb_test, "XGB", "Test")
)

print(metrics_df)



library(tidyr)

metrics_long <- pivot_longer(metrics_df,
  cols = c(Accuracy, BalancedAccuracy, Sensitivity, Specificity, Precision, F1),
  names_to = "Metric", values_to = "Value"
)



ggplot(metrics_long, aes(x = Model, y = Value, fill = Set)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = c("Train" = "#32289c", "Test" = "#c62525")) +
  labs(title = "Performance Metrics Comparison Across Models",
       y = "Metric Value", x = "ML Model") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")



ggsave("final/final_Model_Performance_Metrics_Comparison_for_all_metrics.png", width = 8, height = 8, dpi = 300)






# ================================ Performance Comparison on Accuracy of all ML Models
# Gather accuracies from TRAIN confusion matrices
acc_train <- c(
  RF = conf_matrix_train$overall["Accuracy"],
  SVM = conf_matrix_svm_train$overall["Accuracy"],
  LR = conf_matrix_glm_train$overall["Accuracy"],
  XGB = conf_matrix_xgb_train$overall["Accuracy"]
)

# Gather accuracies from TEST confusion matrices
acc_test <- c(
  RF = conf_matrix_test$overall["Accuracy"],
  SVM = conf_matrix_svm_test$overall["Accuracy"],
  LR = conf_matrix_glm_test$overall["Accuracy"],
  XGB = conf_matrix_xgb_test$overall["Accuracy"]
)

# Create summary data frame
perf_df <- data.frame(
  Model = rep(names(acc_train), times = 2),
  Set = rep(c("Train", "Test"), each = length(acc_train)),
  Accuracy = c(acc_train, acc_test)
)

print(perf_df)


library(ggplot2)

ggplot(perf_df, aes(x = Model, y = Accuracy, fill = Set)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("Train" = "#32289c", "Test" = "#c62525")) +
  ylim(0, 1) +
  labs(title = "Model Performance Comparison",
       y = "Accuracy",
       x = "ML Model") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")
ggsave("final/final_Model_Performance_Comparison_on_accuarcy.png", width = 8, height = 8, dpi = 300)


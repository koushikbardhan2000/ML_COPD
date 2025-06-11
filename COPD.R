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

# List of Unknown samples
Unknown <- data.frame(all_annotated_samples[all_annotated_samples$phenotype == "Unknown", ])
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

# write.csv(phenotype, "output/phenotype_final.csv", row.names = F) # nolint





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
# write.csv(expr_hnVScopd,"output/expr_hnVScopd.csv")
# Smoker vs COPD expression matrix
expr_smokerVScopd <- filtered_expression[,c(smoker_ids,copd_ids)]
# write.csv(expr_smokerVScopd,"output/expr_smokerVScopd.csv")



# Subset phenotype information for the selected samples
pheno_hnVScopd <- phenotype_sorted[phenotype_sorted$GSM_IDs %in% colnames(expr_hnVScopd), ]

# Ensure phenotype vector matches the column order
pheno_hnVScopd <- pheno_hnVScopd[match(colnames(expr_hnVScopd), pheno_hnVScopd$GSM_IDs), ]
#####################
# Pre-processing end
#####################


#############################
# Batch effect removal starts
#############################
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


# log normalization
expr_hnVScopd_log <- log2(expr_hnVScopd + 1)

# Load the limma package for normalization
library(limma)
expr_log2_norm <- normalizeBetweenArrays(expr_hnVScopd_log, method = "quantile")



# Batch effect removal
library(sva)

# Batch correction
batch_corrected <- ComBat(dat = as.matrix(expr_log2_norm), batch = pheno_hnVScopd$phenotype, par.prior = TRUE)
write.csv(batch_corrected, "output/batch_corrected.csv")
head(batch_corrected[1:5, 1:5])
batch_corrected <- read.csv("output/batch_corrected.csv", row.names = 1)


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

#############
# DGEA ends
#############



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


# Transpose expression data for glmnet (samples x genes)
x <- t(batch_corrected)
y <- factor(pheno_hnVScopd$phenotype)

# Binary classification: 0 = Healthy Non-Smoker, 1 = Smoker with COPD
y_bin <- ifelse(y == "Smoker with COPD", 1, 0)

# Standardize features
x_scaled <- scale(x)

# Create cross-validation folds
set.seed(123)
cv_lasso <- cv.glmnet(x_scaled, y_bin, alpha = 1, family = "binomial", nfolds = 10)

# Plot CV results
plot(cv_lasso)

# Best lambda
best_lambda <- 0.05

# Coefficients at best lambda
lasso_coef <- coef(cv_lasso, s = best_lambda)
selected_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
selected_genes <- selected_genes[!selected_genes %in% "(Intercept)"]

# Expression matrix with selected genes
lasso_expr <- x[, selected_genes]

###############
# Random Forest
###############

# Convert phenotype to factor with two levels
rf_df <- data.frame(lasso_expr, phenotype = factor(y))

# Set up repeated RF modeling
set.seed(123)
n_vars <- 1:length(selected_genes)
error_rates <- numeric(length(n_vars))

for (i in n_vars) {
  rf_model <- randomForest(phenotype ~ ., 
                           data = rf_df[, c(1:i, ncol(rf_df))],
                           ntree = 500,
                           importance = TRUE)
  error_rates[i] <- 1 - rf_model$confusion["Smoker with COPD", "Smoker with COPD"] / sum(rf_model$confusion["Smoker with COPD", ])
}

# Plot error rate vs number of variables
png("temp/temp_error_rate_vs_genes.png", width = 1200, height = 1000, res = 150)
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
varImpPlot(final_rf, type = 1, main = "Top Gene Importance (MeanDecreaseAccuracy)")


library(caret)

ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
rf_cv <- train(phenotype ~ ., data = rf_df[, c(best_genes, "phenotype")],
               method = "rf", trControl = ctrl)
print(rf_cv)



# performancs evaluation
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
png("temp/temp_AUROC_RF.png", width = 1200, height = 1000, res = 150)
plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve - Random Forest")
abline(a = 0, b = 1, lty = 2, col = "gray")
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
png("temp/temp_AUROC_SVM.png", width = 1200, height = 1000, res = 150)
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
png("temp/temp_AUROC_GLM.png", width = 1200, height = 1000, res = 150)
plot(glm_roc, col = "red", lwd = 2, main = "ROC Curve - Logistic Regression")
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
png("temp/temp_AUROC_XGB.png", width = 1200, height = 1000, res = 150)
plot(xgb_roc, col = "darkorange", lwd = 2, main = "ROC Curve - XGBoost")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()




# ============================ this portion is moved to the final script
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

# Assuming `rf_df` and `best_genes` are already available in environment

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
# ============================ this portion is moved to the final script





# ================================ this portion is made for Maitree and has been sent to her by WP 
####################################
# Network Centrality Analysis starts
####################################

# install.packages("igraph")
# install.packages("CINNA")

library(igraph)
library(CINNA)

# Load the PPI data
ppi <- read.delim("DEGs_string_interactions_short.tsv", header = TRUE, sep = "\t")
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
png("PCA_Variable_Contributions_with_direction.png", width = 4000, height = 4000, res = 300)
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
png("PCA_Centrality_Contributions_to_PC1_barchart.png", width = 4000, height = 4000, res = 300)
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


##############
# LASSO starts
##############
# Load required libraries
library(glmnet)
library(randomForest)
library(caret)
library(doParallel)


# STEP 1: Select the hub genes from the centrality analysis
hub_genes <- top20_hub_genes$Gene

# STEP 2: Extract expression data for top DEGs
expr_top_hub <- batch_corrected[hub_genes, ]

# Transpose expression data for glmnet (samples x genes)
x <- t(expr_top_hub)
y <- factor(pheno_hnVScopd$phenotype)

# Binary classification: 0 = Healthy Non-Smoker, 1 = Smoker with COPD
y_bin <- ifelse(y == "AM", 1, 0)

# Standardize features
# x_scaled <- scale(x)

# Create cross-validation folds
set.seed(123)
cv_lasso <- cv.glmnet(x, y_bin, alpha = 1, family = "binomial", nfolds = 10) # x_scaled not used here

# Plot CV results
# Plot with lambda.min and lambda.1se indicated
png("LASSO_CV_plot.png", width = 4000, height = 4000, res = 300)
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
png("LASSO_Coefficient_Profiles.png", width = 4000, height = 4000, res = 300)
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
cat("Selected genes at lambda.min:", length(selected_genes), "\n")

# Expression matrix with selected genes
lasso_expr <- x[, selected_genes]
head(lasso_expr[1:5, 1:5])

###############
# Random Forest to find the best genes
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
  error_rates[i] <- 1 - rf_model$confusion["AM", "AM"] / sum(rf_model$confusion["AM", ])
}

# Plot error rate vs number of variables
png("error_rate_vs_genes.png", width = 1200, height = 1200, res = 150)
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

# after this point run the ML models on the best_genes and validate the results
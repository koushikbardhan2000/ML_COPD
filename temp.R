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





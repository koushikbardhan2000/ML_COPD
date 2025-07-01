# This file is 2nd part of the COPD analysis using LASSO and various ML models.
# Before this file, you should have run the final.R script to generate the required datasets.
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
plot(rf_test_roc, col = "blue", lwd = 2, add = TRUE)
abline(a = 0, b = 1, col = "gray", lty = 2)
legend("bottomright", legend = c(
  paste0("Train AUC = ", round(auc(rf_train_roc), 4)),
  paste0("Test AUC = ", round(auc(rf_test_roc), 4))),
  col = c("darkblue", "blue"), lwd = 2)
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
plot(svm_test_roc, col = "blue", lwd = 2, add = TRUE)
abline(a = 0, b = 1, col = "gray", lty = 2)
legend("bottomright", legend = c(
  paste0("Train AUC = ", round(auc(svm_train_roc), 4)),
  paste0("Test AUC = ", round(auc(svm_test_roc), 4))),
  col = c("red", "blue"), lwd = 2)
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
plot(glm_test_roc, col = "blue", lwd = 2, add = TRUE)
abline(a = 0, b = 1, col = "gray", lty = 2)
legend("bottomright", legend = c(
  paste0("Train AUC = ", round(auc(glm_train_roc), 4)),
  paste0("Test AUC = ", round(auc(glm_test_roc), 4))),
  col = c("darkgreen", "blue"), lwd = 2)
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
plot(xgb_test_roc, col = "blue", lwd = 2, add = TRUE)
abline(a = 0, b = 1, col = "gray", lty = 2)
legend("bottomright", legend = c(
  paste0("Train AUC = ", round(auc(xgb_train_roc), 4)),
  paste0("Test AUC = ", round(auc(xgb_test_roc), 4))),
  col = c("orange", "blue"), lwd = 2)
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



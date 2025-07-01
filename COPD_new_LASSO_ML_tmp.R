


set.seed(123)
train_idx <- createDataPartition(rf_df$phenotype, p = 0.7, list = FALSE)
train_df <- rf_df[train_idx, ]
test_df <- rf_df[-train_idx, ]

y_train_bin <- ifelse(train_df$phenotype == "Smoker with COPD", 1, 0)
y_test_bin <- ifelse(test_df$phenotype == "Smoker with COPD", 1, 0)



# Train
rf_model <- randomForest(phenotype ~ ., data = train_df[, c(best_genes, "phenotype")],
                         ntree = 500, importance = TRUE)

# Predict
rf_train_pred <- predict(rf_model, train_df[, best_genes])
rf_test_pred <- predict(rf_model, test_df[, best_genes])

rf_train_prob <- predict(rf_model, train_df[, best_genes], type = "prob")[, "Smoker with COPD"]
rf_test_prob <- predict(rf_model, test_df[, best_genes], type = "prob")[, "Smoker with COPD"]

# Confusion Matrices
rf_train_conf <- confusionMatrix(rf_train_pred, train_df$phenotype, positive = "Smoker with COPD")
rf_test_conf <- confusionMatrix(rf_test_pred, test_df$phenotype, positive = "Smoker with COPD")

# AUROC
rf_train_auc <- auc(roc(y_train_bin, rf_train_prob))
rf_test_auc  <- auc(roc(y_test_bin, rf_test_prob))

# Print
cat("\n=== Random Forest Metrics ===\n")
cat("Train - Accuracy:", rf_train_conf$overall["Accuracy"], 
    " | AUROC:", round(rf_train_auc, 4), "\n")
cat("Test  - Accuracy:", rf_test_conf$overall["Accuracy"], 
    " | AUROC:", round(rf_test_auc, 4), "\n")




# Train SVM with radial kernel
svm_model <- svm(phenotype ~ ., data = train_df[, c(best_genes, "phenotype")],
                 kernel = "radial", probability = TRUE)

# Predict
svm_train_pred <- predict(svm_model, train_df[, best_genes])
svm_test_pred <- predict(svm_model, test_df[, best_genes])

svm_train_prob <- attr(predict(svm_model, train_df[, best_genes], probability = TRUE), "probabilities")[, "Smoker with COPD"]
svm_test_prob <- attr(predict(svm_model, test_df[, best_genes], probability = TRUE), "probabilities")[, "Smoker with COPD"]

# Confusion Matrices
svm_train_conf <- confusionMatrix(svm_train_pred, train_df$phenotype, positive = "Smoker with COPD")
svm_test_conf <- confusionMatrix(svm_test_pred, test_df$phenotype, positive = "Smoker with COPD")

# AUROC
svm_train_auc <- auc(roc(y_train_bin, svm_train_prob))
svm_test_auc <- auc(roc(y_test_bin, svm_test_prob))

# Print
cat("\n=== SVM Metrics ===\n")
cat("Train - Accuracy:", svm_train_conf$overall["Accuracy"],
    " | AUROC:", round(svm_train_auc, 4), "\n")
cat("Test  - Accuracy:", svm_test_conf$overall["Accuracy"],
    " | AUROC:", round(svm_test_auc, 4), "\n")




# Train
glm_model <- glm(phenotype ~ ., data = train_df[, c(best_genes, "phenotype")], family = "binomial")

# Predict
glm_train_prob <- predict(glm_model, train_df[, best_genes], type = "response")
glm_test_prob  <- predict(glm_model, test_df[, best_genes], type = "response")

glm_train_pred <- factor(ifelse(glm_train_prob > 0.5, "Smoker with COPD", "Healthy Non-Smoker"),
                         levels = levels(train_df$phenotype))
glm_test_pred <- factor(ifelse(glm_test_prob > 0.5, "Smoker with COPD", "Healthy Non-Smoker"),
                        levels = levels(test_df$phenotype))

# Confusion Matrices
glm_train_conf <- confusionMatrix(glm_train_pred, train_df$phenotype, positive = "Smoker with COPD")
glm_test_conf <- confusionMatrix(glm_test_pred, test_df$phenotype, positive = "Smoker with COPD")

# AUROC
glm_train_auc <- auc(roc(y_train_bin, glm_train_prob))
glm_test_auc  <- auc(roc(y_test_bin, glm_test_prob))

# Print
cat("\n=== Logistic Regression Metrics ===\n")
cat("Train - Accuracy:", glm_train_conf$overall["Accuracy"],
    " | AUROC:", round(glm_train_auc, 4), "\n")
cat("Test  - Accuracy:", glm_test_conf$overall["Accuracy"],
    " | AUROC:", round(glm_test_auc, 4), "\n")





# XGBoost format
train_xgb <- train_df
test_xgb <- test_df
train_xgb$phenotype <- factor(ifelse(train_xgb$phenotype == "Smoker with COPD", "COPD", "Healthy"))
test_xgb$phenotype  <- factor(ifelse(test_xgb$phenotype == "Smoker with COPD", "COPD", "Healthy"))

# Train
xgb_ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5,
                         classProbs = TRUE, summaryFunction = twoClassSummary)
xgb_model <- train(phenotype ~ ., data = train_xgb[, c(best_genes, "phenotype")],
                   method = "xgbTree", trControl = xgb_ctrl, metric = "ROC")

# Predict
xgb_train_pred <- predict(xgb_model, train_xgb[, best_genes])
xgb_test_pred <- predict(xgb_model, test_xgb[, best_genes])

xgb_train_prob <- predict(xgb_model, train_xgb[, best_genes], type = "prob")[, "COPD"]
xgb_test_prob  <- predict(xgb_model, test_xgb[, best_genes], type = "prob")[, "COPD"]

# Confusion Matrices
xgb_train_conf <- confusionMatrix(xgb_train_pred, train_xgb$phenotype, positive = "COPD")
xgb_test_conf <- confusionMatrix(xgb_test_pred, test_xgb$phenotype, positive = "COPD")

# AUROC
xgb_train_auc <- auc(roc(ifelse(train_xgb$phenotype == "COPD", 1, 0), xgb_train_prob))
xgb_test_auc <- auc(roc(ifelse(test_xgb$phenotype == "COPD", 1, 0), xgb_test_prob))

# Print
cat("\n=== XGBoost Metrics ===\n")
cat("Train - Accuracy:", xgb_train_conf$overall["Accuracy"],
    " | AUROC:", round(xgb_train_auc, 4), "\n")
cat("Test  - Accuracy:", xgb_test_conf$overall["Accuracy"],
    " | AUROC:", round(xgb_test_auc, 4), "\n")






# ----------------------------------------
# 1. Compute Train/Test Metrics for Each Model
# ----------------------------------------

# Logistic Regression
glm_train_acc <- glm_train_conf$overall["Accuracy"]
glm_test_acc  <- glm_test_conf$overall["Accuracy"]
glm_train_sens <- glm_train_conf$byClass["Sensitivity"]
glm_test_sens  <- glm_test_conf$byClass["Sensitivity"]
glm_train_spec <- glm_train_conf$byClass["Specificity"]
glm_test_spec  <- glm_test_conf$byClass["Specificity"]
glm_train_f1   <- glm_train_conf$byClass["F1"]
glm_test_f1    <- glm_test_conf$byClass["F1"]

# SVM
svm_train_conf <- confusionMatrix(svm_train_pred, y_train, positive = "Smoker with COPD")
svm_test_conf  <- confusionMatrix(svm_test_pred,  y_test,  positive = "Smoker with COPD")
svm_train_acc  <- svm_train_conf$overall["Accuracy"]
svm_test_acc   <- svm_test_conf$overall["Accuracy"]
svm_train_sens <- svm_train_conf$byClass["Sensitivity"]
svm_test_sens  <- svm_test_conf$byClass["Sensitivity"]
svm_train_spec <- svm_train_conf$byClass["Specificity"]
svm_test_spec  <- svm_test_conf$byClass["Specificity"]
svm_train_f1   <- svm_train_conf$byClass["F1"]
svm_test_f1    <- svm_test_conf$byClass["F1"]

# XGBoost
xgb_train_conf <- confusionMatrix(xgb_train_pred, train_xgb$phenotype, positive = "COPD")
xgb_test_conf  <- confusionMatrix(xgb_test_pred,  test_xgb$phenotype,  positive = "COPD")
xgb_train_acc  <- xgb_train_conf$overall["Accuracy"]
xgb_test_acc   <- xgb_test_conf$overall["Accuracy"]
xgb_train_sens <- xgb_train_conf$byClass["Sensitivity"]
xgb_test_sens  <- xgb_test_conf$byClass["Sensitivity"]
xgb_train_spec <- xgb_train_conf$byClass["Specificity"]
xgb_test_spec  <- xgb_test_conf$byClass["Specificity"]
xgb_train_f1   <- xgb_train_conf$byClass["F1"]
xgb_test_f1    <- xgb_test_conf$byClass["F1"]

# Random Forest
rf_train_conf <- confusionMatrix(predict(rf_model, rf_train_df[, -ncol(rf_train_df)]), y_train, positive = "Smoker with COPD")
rf_test_conf  <- confusionMatrix(predict(rf_model, rf_test_df[, -ncol(rf_test_df)]),  y_test,  positive = "Smoker with COPD")
rf_train_acc  <- rf_train_conf$overall["Accuracy"]
rf_test_acc   <- rf_test_conf$overall["Accuracy"]
rf_train_sens <- rf_train_conf$byClass["Sensitivity"]
rf_test_sens  <- rf_test_conf$byClass["Sensitivity"]
rf_train_spec <- rf_train_conf$byClass["Specificity"]
rf_test_spec  <- rf_test_conf$byClass["Specificity"]
rf_train_f1   <- rf_train_conf$byClass["F1"]
rf_test_f1    <- rf_test_conf$byClass["F1"]

# ----------------------------------------
# 2. Assemble Summary Performance Table
# ----------------------------------------
performance_summary <- data.frame(
  Model            = c("Logistic Regression", "SVM", "XGBoost", "Random Forest"),
  Train_AUC        = c(auc(glm_train_roc), auc(svm_train_roc), auc(xgb_train_roc), auc(rf_train_roc)),
  Test_AUC         = c(auc(glm_test_roc),  auc(svm_test_roc),  auc(xgb_test_roc),  auc(rf_test_roc)),
  Train_Accuracy   = c(glm_train_acc,  svm_train_acc,  xgb_train_acc,  rf_train_acc),
  Test_Accuracy    = c(glm_test_acc,   svm_test_acc,   xgb_test_acc,   rf_test_acc),
  Train_Sensitivity= c(glm_train_sens, svm_train_sens, xgb_train_sens, rf_train_sens),
  Test_Sensitivity = c(glm_test_sens,  svm_test_sens,  xgb_test_sens,  rf_test_sens),
  Train_Specificity= c(glm_train_spec, svm_train_spec, xgb_train_spec, rf_train_spec),
  Test_Specificity = c(glm_test_spec,  svm_test_spec,  xgb_test_spec,  rf_test_spec),
  Train_F1         = c(glm_train_f1,   svm_train_f1,   xgb_train_f1,   rf_train_f1),
  Test_F1          = c(glm_test_f1,    svm_test_f1,    xgb_test_f1,    rf_test_f1)
)

# Round and display
performance_summary[-1] <- round(performance_summary[-1], 4)
print(performance_summary)

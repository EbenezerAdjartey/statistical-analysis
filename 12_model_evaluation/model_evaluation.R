# =============================================================================
# MODEL EVALUATION & TESTING IN R
# =============================================================================
# Topics: Train/test split, k-fold CV, confusion matrix, ROC/AUC,
#         regression metrics, hyperparameter tuning, SHAP, learning curves
# Author: Ebenezer Adjartey
# =============================================================================

pkgs <- c("ggplot2","dplyr","caret","pROC","glmnet","randomForest",
          "ROCR","reshape2","tidyr")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(ggplot2); library(caret); library(pROC); library(randomForest); library(ROCR)
set.seed(42)

# ── 1. Generate Dataset ───────────────────────────────────────────────────────
n <- 800
X <- matrix(rnorm(n*15), n, 15)
true_b <- c(1.5,-1,.8,0,0,.5,-.7,0,.3,0, 0,0,0,0,0)
prob    <- 1/(1+exp(-(X %*% true_b)))
y       <- factor(rbinom(n,1,prob), labels=c("Neg","Pos"))
df      <- as.data.frame(cbind(X, y=as.integer(y=="Pos")))
df$y    <- factor(df$y, labels=c("Neg","Pos"))
colnames(df)[1:15] <- paste0("X",1:15)

# ── 2. Train/Test Split ───────────────────────────────────────────────────────
cat("=== TRAIN/TEST SPLIT ===\n")
idx   <- createDataPartition(df$y, p=0.80, list=FALSE)
train <- df[idx,]; test  <- df[-idx,]
cat(sprintf("Train: %d  Test: %d\n", nrow(train), nrow(test)))
cat("Class balance (train):", table(train$y), "\n\n")

# ── 3. K-Fold Cross-Validation ────────────────────────────────────────────────
cat("=== K-FOLD CROSS-VALIDATION ===\n")
ctrl_cv <- trainControl(method="cv", number=5, classProbs=TRUE,
                         summaryFunction=twoClassSummary, savePredictions=TRUE)
rf_cv <- train(y ~ ., data=train, method="rf", ntree=50,
                trControl=ctrl_cv, metric="ROC")
cat("5-fold CV results:\n"); print(rf_cv$results)
cat(sprintf("Mean AUC: %.4f  SD: %.4f\n\n", mean(rf_cv$resample$ROC), sd(rf_cv$resample$ROC)))

# ── 4. Train and Predict ──────────────────────────────────────────────────────
rf_final <- randomForest(y ~ ., data=train, ntree=100)
rf_pred  <- predict(rf_final, test)
rf_prob  <- predict(rf_final, test, type="prob")[,"Pos"]

# ── 5. Confusion Matrix ───────────────────────────────────────────────────────
cat("=== CONFUSION MATRIX ===\n")
cm <- confusionMatrix(rf_pred, test$y, positive="Pos")
print(cm)
tn <- cm$table[1,1]; fp <- cm$table[1,2]
fn <- cm$table[2,1]; tp <- cm$table[2,2]
cat(sprintf("\nTP=%d FP=%d TN=%d FN=%d\n",  tp, fp, tn, fn))
cat(sprintf("Accuracy:    %.4f\n",           cm$overall["Accuracy"]))
cat(sprintf("Precision:   %.4f  (PPV)\n",    cm$byClass["Precision"]))
cat(sprintf("Recall:      %.4f  (Sensitivity)\n", cm$byClass["Sensitivity"]))
cat(sprintf("F1-Score:    %.4f\n",           cm$byClass["F1"]))
cat(sprintf("Specificity: %.4f\n",           cm$byClass["Specificity"]))

# ── 6. ROC Curve and AUC ─────────────────────────────────────────────────────
cat("\n=== ROC CURVE AND AUC ===\n")
roc_obj <- roc(test$y, rf_prob, levels=c("Neg","Pos"), direction="<")
cat(sprintf("AUC = %.4f  (CI: %.4f - %.4f)\n",
            auc(roc_obj),
            ci.auc(roc_obj)[1],
            ci.auc(roc_obj)[3]))

# Multiple models ROC comparison
models_roc <- list()
for (meth in c("glm","rf","rpart")) {
  m <- train(y ~ ., data=train, method=meth,
              trControl=trainControl(method="none", classProbs=TRUE),
              metric="Accuracy",
              family=if(meth=="glm")"binomial" else NULL)
  prob_m <- predict(m, test, type="prob")[,"Pos"]
  models_roc[[meth]] <- roc(test$y, prob_m, levels=c("Neg","Pos"), direction="<")
}

dir.create("12_model_evaluation", showWarnings=FALSE)
png("12_model_evaluation/roc_comparison.png", width=600, height=500)
plot(models_roc[["glm"]], col="blue",   main="ROC Curve Comparison", lwd=2)
plot(models_roc[["rf"]],  col="red",    add=TRUE, lwd=2)
plot(models_roc[["rpart"]],col="green", add=TRUE, lwd=2)
legend("bottomright", legend=c(
  sprintf("Logistic AUC=%.3f", auc(models_roc[["glm"]])),
  sprintf("RF       AUC=%.3f", auc(models_roc[["rf"]])),
  sprintf("Tree     AUC=%.3f", auc(models_roc[["rpart"]]))),
  col=c("blue","red","green"), lwd=2)
dev.off()

# ── 7. Regression Metrics ─────────────────────────────────────────────────────
cat("\n=== REGRESSION METRICS ===\n")
y_reg   <- 1 + 2*X[,1] - 1.5*X[,2] + rnorm(n)
df_reg  <- as.data.frame(cbind(X, y=y_reg))
colnames(df_reg)[1:15] <- paste0("X",1:15)
tr_reg  <- df_reg[idx,]; te_reg <- df_reg[-idx,]

lm_reg  <- lm(y ~ ., data=tr_reg)
y_pred  <- predict(lm_reg, te_reg)
y_true  <- te_reg$y

mae   <- mean(abs(y_true - y_pred))
rmse  <- sqrt(mean((y_true - y_pred)^2))
r2    <- 1 - sum((y_true-y_pred)^2) / sum((y_true-mean(y_true))^2)
mape  <- mean(abs((y_true-y_pred)/(abs(y_true)+1e-8)))*100

cat(sprintf("MAE:  %.4f\n",  mae))
cat(sprintf("RMSE: %.4f\n",  rmse))
cat(sprintf("R2:   %.4f\n",  r2))
cat(sprintf("MAPE: %.2f%%\n\n", mape))

# ── 8. Hyperparameter Tuning (Grid Search) ────────────────────────────────────
cat("=== HYPERPARAMETER TUNING ===\n")
tgrid <- expand.grid(mtry=c(3,5,8))
ctrl_tune <- trainControl(method="cv", number=3, classProbs=TRUE,
                           summaryFunction=twoClassSummary)
rf_tuned <- train(y ~ ., data=train, method="rf",
                   trControl=ctrl_tune, tuneGrid=tgrid,
                   metric="ROC", ntree=50)
cat("Grid search results:\n"); print(rf_tuned$results)
cat("Best mtry:", rf_tuned$bestTune$mtry, "\n\n")

# ── 9. Feature Importance ─────────────────────────────────────────────────────
cat("=== FEATURE IMPORTANCE ===\n")
fi <- importance(rf_final)
fi_df <- data.frame(
  Feature    = rownames(fi),
  MeanDecGini= fi[,"MeanDecreaseGini"]
) %>% arrange(desc(MeanDecGini))
print(head(fi_df, 10))

p_fi <- ggplot(head(fi_df,10), aes(x=reorder(Feature,MeanDecGini), y=MeanDecGini)) +
  geom_col(fill="steelblue", alpha=0.8) + coord_flip() +
  labs(title="Random Forest Feature Importance (top 10)",
       x="Feature", y="Mean Decrease Gini") + theme_minimal()
ggsave("12_model_evaluation/feature_importance.png", p_fi, width=7, height=5)

# ── 10. Learning Curves ───────────────────────────────────────────────────────
cat("\n=== LEARNING CURVES ===\n")
train_sizes <- floor(nrow(train) * seq(0.1, 1.0, by=0.1))
lc_results  <- lapply(train_sizes, function(sz) {
  sub <- train[sample(nrow(train), sz),]
  m   <- randomForest(y ~ ., data=sub, ntree=50)
  tr_acc <- mean(predict(m,sub)   == sub$y)
  te_acc <- mean(predict(m,test)  == test$y)
  data.frame(n=sz, train_acc=tr_acc, test_acc=te_acc)
})
lc_df <- do.call(rbind, lc_results)

p_lc <- ggplot(pivot_longer(lc_df, c("train_acc","test_acc"),
                              names_to="set", values_to="acc"),
               aes(x=n, y=acc, color=set)) +
  geom_line(linewidth=1.2) + geom_point() +
  labs(title="Learning Curves (Random Forest)",
       x="Training Size", y="Accuracy", color="Set") +
  theme_minimal()
ggsave("12_model_evaluation/learning_curves.png", p_lc, width=7, height=5)

# ── 11. Threshold Tuning ──────────────────────────────────────────────────────
cat("\n=== THRESHOLD TUNING ===\n")
thresholds <- seq(0.1, 0.9, by=0.05)
thresh_res <- lapply(thresholds, function(t) {
  pred <- factor(ifelse(rf_prob >= t, "Pos","Neg"), levels=c("Neg","Pos"))
  cm_t <- confusionMatrix(pred, test$y, positive="Pos")
  data.frame(threshold=t,
             precision=cm_t$byClass["Precision"],
             recall=cm_t$byClass["Sensitivity"],
             f1=cm_t$byClass["F1"])
})
thresh_df <- do.call(rbind, thresh_res)
best_t <- thresh_df[which.max(thresh_df$f1), "threshold"]
cat(sprintf("Best threshold (max F1): %.2f  F1=%.4f\n",
            best_t, max(thresh_df$f1, na.rm=TRUE)))

cat("\n=== MODEL EVALUATION COMPLETE ===\n")

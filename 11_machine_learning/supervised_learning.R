# =============================================================================
# SUPERVISED MACHINE LEARNING IN R
# =============================================================================
# Topics: Linear/Ridge/Lasso, Logistic, Decision Trees, Random Forest,
#         Gradient Boosting (XGBoost), SVM, KNN, Neural Network
# Author: Ebenezer Adjartey
# =============================================================================

pkgs <- c("ggplot2","dplyr","caret","glmnet","rpart","rpart.plot",
          "randomForest","gbm","e1071","class","nnet","pROC","xgboost")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(ggplot2); library(caret); library(glmnet); library(rpart)
library(randomForest); library(gbm); library(e1071); library(pROC)
set.seed(42)

# ── 1. Generate Datasets ──────────────────────────────────────────────────────
n <- 600
# Classification
X_cls <- matrix(rnorm(n*10), n, 10)
true_b <- c(1.5, -1, 0.8, 0, 0, 0.5, -0.7, 0, 0.3, 0)
prob    <- 1/(1+exp(-(X_cls %*% true_b)))
y_cls   <- rbinom(n, 1, prob)
df_cls  <- as.data.frame(cbind(X_cls, y=y_cls))
colnames(df_cls)[1:10] <- paste0("X",1:10)
df_cls$y <- factor(df_cls$y, labels=c("No","Yes"))

# Regression
y_reg <- 1 + 2*X_cls[,1] - 1.5*X_cls[,2] + rnorm(n)
df_reg <- as.data.frame(cbind(X_cls, y=y_reg))
colnames(df_reg)[1:10] <- paste0("X",1:10)

# Train/test split
train_idx <- createDataPartition(df_cls$y, p=0.75, list=FALSE)
train_c <- df_cls[train_idx,];  test_c <- df_cls[-train_idx,]
train_r <- df_reg[train_idx,];  test_r <- df_reg[-train_idx,]
cat("Train:", nrow(train_c), "  Test:", nrow(test_c), "\n\n")

# ── 2. Linear/Ridge/Lasso Regression ─────────────────────────────────────────
cat("=== PENALIZED REGRESSION (OLS, Ridge, Lasso) ===\n")
X_tr <- as.matrix(train_r[,1:10])
X_te <- as.matrix(test_r[,1:10])
y_tr <- train_r$y; y_te <- test_r$y

# OLS
ols_r <- lm(y ~ ., data=train_r)
rmse_ols <- sqrt(mean((predict(ols_r, test_r) - y_te)^2))
cat(sprintf("OLS  RMSE=%.4f  R2=%.4f\n", rmse_ols, cor(predict(ols_r,test_r),y_te)^2))

# Ridge (alpha=0)
ridge_cv <- cv.glmnet(X_tr, y_tr, alpha=0)
pred_ridge <- predict(ridge_cv, X_te, s="lambda.min")
rmse_ridge <- sqrt(mean((pred_ridge-y_te)^2))
cat(sprintf("Ridge RMSE=%.4f  lambda.min=%.4f\n", rmse_ridge, ridge_cv$lambda.min))

# Lasso (alpha=1)
lasso_cv <- cv.glmnet(X_tr, y_tr, alpha=1)
pred_lasso <- predict(lasso_cv, X_te, s="lambda.min")
rmse_lasso <- sqrt(mean((pred_lasso-y_te)^2))
nonzero    <- sum(coef(lasso_cv, s="lambda.min") != 0) - 1
cat(sprintf("Lasso RMSE=%.4f  nonzero_coefs=%d\n\n", rmse_lasso, nonzero))

# ── 3. Logistic Regression ────────────────────────────────────────────────────
cat("=== LOGISTIC REGRESSION ===\n")
logit_m <- glm(y ~ ., data=train_c, family=binomial)
logit_prob <- predict(logit_m, test_c, type="response")
logit_pred <- factor(ifelse(logit_prob>0.5,"Yes","No"), levels=c("No","Yes"))
cm_logit <- confusionMatrix(logit_pred, test_c$y, positive="Yes")
auc_logit <- auc(roc(as.numeric(test_c$y)-1, logit_prob))
cat(sprintf("Logistic: Acc=%.4f  AUC=%.4f\n\n", cm_logit$overall["Accuracy"], auc_logit))

# ── 4. Decision Tree ──────────────────────────────────────────────────────────
cat("=== DECISION TREE ===\n")
dt <- rpart(y ~ ., data=train_c, method="class", cp=0.01)
dt_pred <- predict(dt, test_c, type="class")
dt_prob <- predict(dt, test_c, type="prob")[,"Yes"]
cm_dt   <- confusionMatrix(dt_pred, test_c$y, positive="Yes")
auc_dt  <- auc(roc(as.numeric(test_c$y)-1, dt_prob))
cat(sprintf("Decision Tree: Acc=%.4f  AUC=%.4f  Depth=%d\n\n",
            cm_dt$overall["Accuracy"], auc_dt, max(dt$frame$depth)))

# ── 5. Random Forest ──────────────────────────────────────────────────────────
cat("=== RANDOM FOREST ===\n")
rf <- randomForest(y ~ ., data=train_c, ntree=100, importance=TRUE)
rf_pred <- predict(rf, test_c)
rf_prob <- predict(rf, test_c, type="prob")[,"Yes"]
cm_rf   <- confusionMatrix(rf_pred, test_c$y, positive="Yes")
auc_rf  <- auc(roc(as.numeric(test_c$y)-1, rf_prob))
cat(sprintf("Random Forest: Acc=%.4f  AUC=%.4f\n\n", cm_rf$overall["Accuracy"], auc_rf))
cat("Variable importance:\n"); print(importance(rf)[order(-importance(rf)[,"MeanDecreaseGini"]),])

# ── 6. Gradient Boosting ──────────────────────────────────────────────────────
cat("=== GRADIENT BOOSTING ===\n")
train_gb <- train_c; test_gb <- test_c
train_gb$y <- as.integer(train_gb$y=="Yes")

gb <- gbm(y ~ ., data=train_gb, distribution="bernoulli",
          n.trees=100, interaction.depth=3, shrinkage=0.1, verbose=FALSE)
gb_prob <- predict(gb, test_c, n.trees=100, type="response")
gb_pred <- factor(ifelse(gb_prob>0.5,"Yes","No"), levels=c("No","Yes"))
auc_gb  <- auc(roc(as.numeric(test_c$y)-1, gb_prob))
cat(sprintf("GBM: AUC=%.4f\n\n", auc_gb))

# ── 7. SVM ────────────────────────────────────────────────────────────────────
cat("=== SUPPORT VECTOR MACHINE ===\n")
svm_m <- svm(y ~ ., data=train_c, kernel="radial", probability=TRUE)
svm_pred <- predict(svm_m, test_c)
svm_prob <- attr(predict(svm_m, test_c, probability=TRUE),"probabilities")[,"Yes"]
cm_svm   <- confusionMatrix(svm_pred, test_c$y, positive="Yes")
auc_svm  <- auc(roc(as.numeric(test_c$y)-1, svm_prob))
cat(sprintf("SVM: Acc=%.4f  AUC=%.4f\n\n", cm_svm$overall["Accuracy"], auc_svm))

# ── 8. KNN ────────────────────────────────────────────────────────────────────
cat("=== KNN (k=5) ===\n")
pre_proc <- preProcess(train_c[,1:10], method=c("center","scale"))
tr_sc    <- predict(pre_proc, train_c[,1:10])
te_sc    <- predict(pre_proc, test_c[,1:10])
knn_pred <- class::knn(tr_sc, te_sc, train_c$y, k=5)
cat(sprintf("KNN: Acc=%.4f\n\n", mean(knn_pred==test_c$y)))

# ── 9. Neural Network ────────────────────────────────────────────────────────
cat("=== NEURAL NETWORK (MLP) ===\n")
nn_m <- nnet::nnet(y ~ ., data=train_c, size=10, decay=0.01,
                    maxit=500, trace=FALSE)
nn_prob <- predict(nn_m, test_c, type="raw")[,1]
nn_pred <- factor(ifelse(nn_prob>0.5,"Yes","No"), levels=c("No","Yes"))
auc_nn  <- auc(roc(as.numeric(test_c$y)-1, nn_prob))
cat(sprintf("Neural Net: AUC=%.4f\n\n", auc_nn))

# ── 10. Model Comparison ─────────────────────────────────────────────────────
cat("=== MODEL COMPARISON ===\n")
comp <- data.frame(
  Model    = c("Logistic","Decision Tree","Random Forest","GBM","SVM","KNN","NNet"),
  Accuracy = c(cm_logit$overall["Accuracy"], cm_dt$overall["Accuracy"],
               cm_rf$overall["Accuracy"],    mean(gb_pred==test_c$y),
               cm_svm$overall["Accuracy"],   mean(knn_pred==test_c$y),
               mean(nn_pred==test_c$y)),
  AUC      = c(auc_logit, auc_dt, auc_rf, auc_gb, auc_svm, NA, auc_nn)
)
comp <- comp[order(-comp$AUC, na.last=TRUE),]
print(round(comp, 4))

# ROC curves
p_roc <- ggroc(list(
  "Random Forest"  = roc(as.numeric(test_c$y)-1, rf_prob),
  "GBM"            = roc(as.numeric(test_c$y)-1, gb_prob),
  "Logistic"       = roc(as.numeric(test_c$y)-1, logit_prob),
  "SVM"            = roc(as.numeric(test_c$y)-1, svm_prob)
)) +
  geom_abline(slope=1, intercept=1, linetype="dashed", color="grey") +
  labs(title="ROC Curves", x="1-Specificity", y="Sensitivity") +
  theme_minimal()

dir.create("11_machine_learning", showWarnings=FALSE)
ggsave("11_machine_learning/roc_curves.png", p_roc, width=7, height=6)

cat("\n=== SUPERVISED LEARNING COMPLETE ===\n")

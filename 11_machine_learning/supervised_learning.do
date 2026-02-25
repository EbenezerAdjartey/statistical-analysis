/* ============================================================================
   SUPERVISED MACHINE LEARNING IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   Note: Stata has limited ML capabilities. This file covers what's available
         natively and via SSC packages (svm, boost, lasso).
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. Generate Dataset ─────────────────────────────────────────────────── */
set obs 600
forvalues i = 1/10 {
    gen x`i' = rnormal()
}
gen latent = 1.5*x1 - x2 + 0.8*x3 + 0.5*x6 - 0.7*x7 + 0.3*x9 + rnormal()
gen y_bin  = (latent > 0)
gen y_cont = 1 + 2*x1 - 1.5*x2 + rnormal()

* Train/test split (75/25)
gen train = (runiform() < 0.75)
di "Train: " (train==1) "  Test: " (train==0)

/* ── 2. OLS Regression ───────────────────────────────────────────────────── */
di _n "=== OLS REGRESSION ==="
regress y_cont x1-x10 if train
predict y_hat_ols
gen resid_ols = y_cont - y_hat_ols if !train
quietly summarize resid_ols if !train
di "OLS Test RMSE = " sqrt(r(Var))  /* approximate */

/* ── 3. LASSO (via Stata 16+ lasso command) ─────────────────────────────── */
di _n "=== LASSO REGRESSION ==="
lasso linear y_cont x1-x10 if train, rseed(42)
cvplot   /* cross-validation plot */
lasso linear y_cont x1-x10 if train, selection(cv) rseed(42)
di "Selected variables:"
lassoinfo
di "Nonzero coefficients:"
coef, display

/* ── 4. Ridge Regression (via lasso with ridge option) ────────────────────── */
di _n "=== RIDGE REGRESSION ==="
lasso linear y_cont x1-x10 if train, alpha(0) rseed(42)  /* alpha=0 = ridge */

/* ── 5. Elastic Net ──────────────────────────────────────────────────────── */
di _n "=== ELASTIC NET ==="
lasso linear y_cont x1-x10 if train, alpha(0.5) rseed(42)  /* alpha=0.5 = elastic net */

/* ── 6. Logistic Regression ──────────────────────────────────────────────── */
di _n "=== LOGISTIC REGRESSION ==="
logit y_bin x1-x10 if train
predict p_logit
classify p_logit y_bin if !train, cutoff(0.5)
estat classification if !train
roctab y_bin p_logit if !train, graph

/* ── 7. LASSO Logistic Regression ────────────────────────────────────────── */
di _n "=== LASSO LOGISTIC ==="
lasso logit y_bin x1-x10 if train, rseed(42) selection(cv)
predict p_lasso_logit
estat classification, cutoff(0.5)

/* ── 8. Decision Tree (via parfit or similar - manual threshold approach) ── */
di _n "=== NOTE: DECISION TREES IN STATA ==="
di "Stata does not have native decision tree commands."
di "For trees, use: ssc install randomforest (if available)"
di "Or use Stata's integration with Python: python: from sklearn..."
di ""
di "Below demonstrates a simple logistic decision boundary approximation:"

* Simple threshold-based classification
gen pred_tree = (x1 > 0) & (x2 < 0)   /* simple rule */
tabulate pred_tree y_bin if !train

/* ── 9. SVM (via svmlight Stata plugin or lsvm) ──────────────────────────── */
di _n "=== SVM (via svmachines if installed) ==="
* ssc install svmachines
* svmachines y_bin x1-x10 if train, type(csvc) kernel(rbf)
di "Requires: ssc install svmachines"
di "Alternative: use Python bridge: python script using sklearn"

/* ── 10. Post-Estimation: ROC and AUC ───────────────────────────────────── */
di _n "=== MODEL EVALUATION: ROC CURVE ==="
logit y_bin x1-x10 if train
predict p_logit2 if !train

roctab y_bin p_logit2 if !train, graph ///
    title("ROC Curve: Logistic Regression") ///
    note("AUC shown in table above")
graph export "11_machine_learning/roc_logit.png", replace

/* ── 11. Cross-Validation Framework ─────────────────────────────────────── */
di _n "=== CROSS-VALIDATION ==="
* Stata's lasso has built-in CV
lasso linear y_cont x1-x10 if train, selection(cv) nfolds(5) rseed(42)
di "5-fold CV RMSE: " e(cvmse)^0.5

/* ── 12. Gradient Boosting (via boost package) ───────────────────────────── */
di _n "=== GRADIENT BOOSTING ==="
* ssc install boost
* boost y_bin x1-x10 if train, dist(bernoulli) shrinkage(0.1) ///
*       trees(100) interaction(3) bag(0.5)
di "Requires: ssc install boost"

/* ── 13. Model Comparison Table ──────────────────────────────────────────── */
di _n "=== MODEL COMPARISON ==="
di "Compare models using AUC and classification accuracy:"
di "(Run each model above and record metrics)"

foreach model in "OLS" "Lasso" "Ridge" "Logistic" "Lasso-Logit" {
    di "`model': see model-specific output above"
}

di _n "=== SUPERVISED LEARNING COMPLETE ==="
di "For full ML capabilities in Stata, use Python integration:"
di "  python: import sklearn, xgboost, keras"
di "  Or use R via rsource/rcall packages"

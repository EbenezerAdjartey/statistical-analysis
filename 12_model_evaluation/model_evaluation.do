/* ============================================================================
   MODEL EVALUATION & TESTING IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. Generate Dataset ─────────────────────────────────────────────────── */
set obs 800
forvalues i = 1/15 {
    gen x`i' = rnormal()
}
gen latent = 1.5*x1 - x2 + 0.8*x3 + 0.5*x6 - 0.7*x7 + 0.3*x9 + rnormal()
gen y = (latent > 0)   /* binary outcome */
gen y_cont = 1 + 2*x1 - 1.5*x2 + rnormal()

* Train/test split (80/20)
gen train = (runiform() < 0.80)
di "Train: " (train==1) "  Test: " (train==0)

/* ── 2. Logistic Regression Model ────────────────────────────────────────── */
di _n "=== LOGISTIC REGRESSION ==="
logit y x1-x10 if train
predict p_train if  train
predict p_test  if !train

/* ── 3. Classification Metrics ───────────────────────────────────────────── */
di _n "=== CLASSIFICATION METRICS ==="
estat classification if !train   /* confusion matrix + accuracy */

* Sensitivity, specificity at cutoff=0.5
di _n "Classification at threshold = 0.5:"
gen yhat = (p_test >= 0.5) if !train

* Manual confusion matrix
count if yhat==1 & y==1 & !train; local TP = r(N)
count if yhat==1 & y==0 & !train; local FP = r(N)
count if yhat==0 & y==1 & !train; local FN = r(N)
count if yhat==0 & y==0 & !train; local TN = r(N)

di "TP=" `TP' "  FP=" `FP' "  TN=" `TN' "  FN=" `FN'
di "Accuracy    = " (`TP'+`TN')/(`TP'+`FP'+`TN'+`FN')
di "Precision   = " `TP'/(`TP'+`FP')
di "Recall      = " `TP'/(`TP'+`FN')
di "Specificity = " `TN'/(`TN'+`FP')
di "F1-Score    = " 2*`TP'/(2*`TP'+`FP'+`FN')

/* ── 4. ROC Curve and AUC ─────────────────────────────────────────────────── */
di _n "=== ROC CURVE AND AUC ==="
roctab y p_test if !train, graph ///
    title("ROC Curve: Logistic Regression") ///
    xtitle("1 - Specificity") ytitle("Sensitivity")
graph export "12_model_evaluation/roc_curve.png", replace

di "AUC = " r(area)

/* ── 5. Alternative Threshold Analysis ──────────────────────────────────── */
di _n "=== THRESHOLD SENSITIVITY ==="
forvalues t = 1(1)9 {
    local thresh = `t'/10
    quietly gen yhat_`t' = (p_test >= `thresh') if !train
    quietly count if yhat_`t'==1 & y==1 & !train; local tp=r(N)
    quietly count if yhat_`t'==1 & y==0 & !train; local fp=r(N)
    quietly count if yhat_`t'==0 & y==1 & !train; local fn=r(N)
    local prec = cond(`tp'+`fp'>0, `tp'/(`tp'+`fp'), 0)
    local rec  = cond(`tp'+`fn'>0, `tp'/(`tp'+`fn'), 0)
    local f1   = cond(`prec'+`rec'>0, 2*`prec'*`rec'/(`prec'+`rec'), 0)
    di "Threshold=" `thresh' "  Precision=" %5.3f `prec' "  Recall=" %5.3f `rec' "  F1=" %5.3f `f1'
    drop yhat_`t'
}

/* ── 6. Cross-Validation via Caret-style Bootstrap ───────────────────────── */
di _n "=== K-FOLD CROSS-VALIDATION (5-fold) ==="
* Stata's crossfold package (ssc install crossfold)
* crossfold logit y x1-x10, k(5) stub(cv)
* Alternative: Manual 5-fold CV
gen fold = mod(_n-1,5)+1

local cv_auc = 0
forvalues k = 1/5 {
    quietly logit y x1-x10 if fold!=`k' & train
    quietly predict p_fold if fold==`k' & train
    quietly roctab y p_fold if fold==`k' & train
    local cv_auc = `cv_auc' + r(area)/5
    quietly drop p_fold
}
di "5-fold CV AUC (train) = " %6.4f `cv_auc'

/* ── 7. Regression Metrics ────────────────────────────────────────────────── */
di _n "=== REGRESSION METRICS ==="
regress y_cont x1-x10 if train
predict yhat_cont if !train
gen resid_cont = y_cont - yhat_cont if !train

quietly summarize resid_cont if !train
local mae  = r(mean)   /* note: this is mean, not MAE */
local rmse = sqrt(r(Var) + r(mean)^2)

* Proper MAE
gen abs_resid = abs(resid_cont) if !train
quietly summarize abs_resid if !train
local mae = r(mean)

* RMSE
quietly gen sq_resid = resid_cont^2 if !train
quietly summarize sq_resid if !train
local rmse = sqrt(r(mean))

* R-squared
quietly correlate y_cont yhat_cont if !train
local r2_test = r(rho)^2

di "Test MAE  = " %8.4f `mae'
di "Test RMSE = " %8.4f `rmse'
di "Test R2   = " %8.4f `r2_test'

/* ── 8. LASSO Cross-Validation ───────────────────────────────────────────── */
di _n "=== LASSO WITH CV (binary) ==="
lasso logit y x1-x15 if train, selection(cv) nfolds(5) rseed(42)
lassoinfo
predict p_lasso if !train
roctab y p_lasso if !train
di "LASSO Logit AUC = " r(area)

/* ── 9. Model Comparison ─────────────────────────────────────────────────── */
di _n "=== MODEL COMPARISON ==="
di "Logistic Regression:"
quietly logit y x1-x10 if train
quietly predict p_logit_te if !train
quietly roctab y p_logit_te if !train
di "  AUC = " %6.4f r(area)

di "Probit Model:"
quietly probit y x1-x10 if train
quietly predict p_probit_te if !train
quietly roctab y p_probit_te if !train
di "  AUC = " %6.4f r(area)

/* ── 10. Feature Importance (via standardized coefficients) ──────────────── */
di _n "=== FEATURE IMPORTANCE (via std. coefficients) ==="
* Standardize variables
foreach v of varlist x1-x10 {
    quietly summarize `v' if train
    gen `v'_std = (`v' - r(mean)) / r(sd)
}
quietly logit y x1_std-x10_std if train
di "Standardized logit coefficients (larger |coef| = more important):"
di "(sorted by |coefficient|)"
matrix B = e(b)
matrix list B

/* ── 11. Calibration Plot ─────────────────────────────────────────────────── */
di _n "=== CALIBRATION: Hosmer-Lemeshow Test ==="
quietly logit y x1-x10 if train
quietly predict p_cal if !train
estat gof, group(10)  /* Hosmer-Lemeshow goodness-of-fit */

di _n "=== MODEL EVALUATION COMPLETE ==="

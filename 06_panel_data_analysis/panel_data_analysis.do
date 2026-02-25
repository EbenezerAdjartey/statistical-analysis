/* ============================================================================
   PANEL DATA ANALYSIS IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. Generate Balanced Panel ─────────────────────────────────────────── */
local N = 100
local T = 10
local n = `N'*`T'
set obs `n'

gen id   = ceil(_n/`T')
gen time = mod(_n-1, `T') + 1
xtset id time

* Individual fixed effects
bysort id: gen alpha_i = rnormal(0,2) if _n==1
bysort id: replace alpha_i = alpha_i[1]

gen x1 = rnormal() + alpha_i*0.3
gen x2 = rnormal()
gen y  = alpha_i + 2*x1 + 1.5*x2 + rnormal()

label var y  "Outcome"
label var x1 "Regressor 1 (correlated with FE)"
label var x2 "Regressor 2"

di _n "Panel: `N' individuals x `T' periods = `n' observations"
xtdescribe

/* ── 2. Pooled OLS ────────────────────────────────────────────────────────── */
di _n "=== POOLED OLS ==="
regress y x1 x2
di "True: x1=2.0, x2=1.5"

/* ── 3. Fixed Effects (Within Estimator) ─────────────────────────────────── */
di _n "=== FIXED EFFECTS (WITHIN) ==="
xtreg y x1 x2, fe
estimates store fe_est
di "FE: coef(x1) = " _b[x1] "  coef(x2) = " _b[x2]

/* ── 4. Random Effects ────────────────────────────────────────────────────── */
di _n "=== RANDOM EFFECTS ==="
xtreg y x1 x2, re
estimates store re_est
di "RE: coef(x1) = " _b[x1] "  coef(x2) = " _b[x2]

/* ── 5. Hausman Test ──────────────────────────────────────────────────────── */
di _n "=== HAUSMAN TEST ==="
hausman fe_est re_est
di "(p<0.05: use Fixed Effects; p>0.05: Random Effects preferred)"

/* ── 6. First-Difference Estimator ──────────────────────────────────────── */
di _n "=== FIRST-DIFFERENCE ESTIMATOR ==="
* First-difference manually
gen dy  = d.y
gen dx1 = d.x1
gen dx2 = d.x2
regress dy dx1 dx2, noconstant

/* ── 7. Two-Way Fixed Effects ────────────────────────────────────────────── */
di _n "=== TWO-WAY FIXED EFFECTS ==="
* Absorb both entity and time FE
reghdfe y x1 x2, absorb(id time)
* ssc install reghdfe if not installed

/* ── 8. Cluster-Robust Standard Errors ───────────────────────────────────── */
di _n "=== CLUSTER-ROBUST SEs (clustered by entity) ==="
xtreg y x1 x2, fe vce(cluster id)

/* ── 9. Panel Diagnostic Tests ───────────────────────────────────────────── */
di _n "=== PANEL DIAGNOSTICS ==="

* F-test for FE (significance of individual effects)
xtreg y x1 x2, fe
xttest3   /* modified Wald test for groupwise heteroskedasticity */

* Breusch-Pagan LM test for RE vs Pooled
xtreg y x1 x2, re
xttest0

* Wooldridge test for serial correlation in panels
xtserial y x1 x2

/* ── 10. Dynamic Panel (Arellano-Bond GMM) ───────────────────────────────── */
di _n "=== DYNAMIC PANEL GMM (Arellano-Bond) ==="
* Requires: ssc install xtabond2
* xtabond2 y l.y x1 x2, gmm(l.y, lag(2 4)) iv(x1 x2) twostep robust
* Below uses built-in xtabond (less flexible)
xtabond y x1 x2, lags(1)
estat sargan   /* Sargan-Hansen overidentification test */

/* ── 11. Visualization ────────────────────────────────────────────────────── */
* Coefficient comparison
estimates clear
quietly xtreg y x1 x2, fe
estimates store fe
quietly xtreg y x1 x2, re
estimates store re
quietly regress y x1 x2
estimates store ols

coefplot fe re ols, ///
    keep(x1 x2) ///
    xline(2, lp(dash) lc(blue)) xline(1.5, lp(dash) lc(red)) ///
    title("Panel Estimators (dashed=true)") ///
    legend(order(1 "FE" 2 "RE" 3 "OLS"))
graph export "06_panel_data_analysis/panel_comparison.png", replace

di _n "=== PANEL DATA ANALYSIS COMPLETE ==="

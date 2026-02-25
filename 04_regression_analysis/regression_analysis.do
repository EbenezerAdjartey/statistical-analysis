/* ============================================================================
   REGRESSION ANALYSIS IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. Generate Data ─────────────────────────────────────────────────────── */
set obs 300
gen educ  = round(8 + runiform()*12)
gen exper = round(runiform()*30)
gen iq    = rnormal(100, 15)
gen fem   = rbinomial(1, 0.5)
gen e     = rnormal(0, 20)

gen wage_log = 10 + 2*educ + 0.5*exper + 0.05*iq - 5*fem + e
gen wage     = exp(wage_log/40)

gen latent   = -5 + 0.3*educ + 0.05*exper + rnormal()
gen employed = (latent > 0)

gen pub_rate = exp(0.5 + 0.1*educ + rnormal(0,.3))
gen pubs     = rpoisson(pub_rate)

gen college_prox = 0.5*educ + rnormal()

label var wage  "Hourly wage"
label var educ  "Years education"
label var exper "Years experience"
label var iq    "IQ score"
label var fem   "Female (=1)"

/* ── 2. Simple OLS ────────────────────────────────────────────────────────── */
di _n "=== SIMPLE OLS: wage ~ educ ==="
regress wage educ
di "R-squared = " %6.4f e(r2)

/* ── 3. Multiple OLS ──────────────────────────────────────────────────────── */
di _n "=== MULTIPLE OLS ==="
regress wage educ exper iq fem
estat ic       /* AIC/BIC */

/* ── 4. OLS Diagnostics ───────────────────────────────────────────────────── */
di _n "=== OLS DIAGNOSTICS ==="

* VIF (multicollinearity)
regress wage educ exper iq fem
estat vif

* Breusch-Pagan / Cook-Weisberg heteroskedasticity test
estat hettest

* Ramsey RESET test (functional form misspecification)
estat ovtest

* Robust standard errors (HC1)
di _n "--- Robust SEs ---"
regress wage educ exper iq fem, robust

* Clustered robust SEs (illustrative)
* regress wage educ exper iq fem, vce(cluster id)

/* ── 5. Logistic Regression ───────────────────────────────────────────────── */
di _n "=== LOGISTIC REGRESSION ==="
logit employed educ exper fem
* Odds ratios
logit employed educ exper fem, or
* Marginal effects
margins, dydx(*)

/* ── 6. Probit Model ──────────────────────────────────────────────────────── */
di _n "=== PROBIT MODEL ==="
probit employed educ exper fem
margins, dydx(*)

/* ── 7. Tobit Model ───────────────────────────────────────────────────────── */
di _n "=== TOBIT MODEL (censored at 0) ==="
tobit wage educ exper fem, ll(0)

/* ── 8. Poisson Regression ────────────────────────────────────────────────── */
di _n "=== POISSON REGRESSION ==="
poisson pubs educ exper
* Incidence rate ratios
poisson pubs educ exper, irr
* Check for overdispersion
quietly poisson pubs educ exper
di "Mean: " r(N) "  Var: "

/* ── 9. Negative Binomial ─────────────────────────────────────────────────── */
di _n "=== NEGATIVE BINOMIAL ==="
nbreg pubs educ exper
nbreg pubs educ exper, irr
di "Alpha (overdispersion) > 0 => overdispersed"

/* ── 10. IV/2SLS Regression ──────────────────────────────────────────────── */
di _n "=== IV/2SLS: educ instrumented by college_prox ==="
ivregress 2sls wage exper (educ = college_prox)
estat endogenous      /* endogeneity test */
estat firststage      /* first-stage F-statistic */
estat overid          /* overidentification test (needs >1 instrument) */

/* ── 11. Quantile Regression ─────────────────────────────────────────────── */
di _n "=== QUANTILE REGRESSION ==="
foreach q of numlist 0.25 0.50 0.75 {
    qreg wage educ exper fem, quantile(`q')
    di "Q(" `q' "): educ=" _b[educ] "  exper=" _b[exper]
}

/* ── 12. Regression table output ─────────────────────────────────────────── */
di _n "=== COMPARING OLS SPECIFICATIONS ==="
quietly regress wage educ
estimates store m1
quietly regress wage educ exper
estimates store m2
quietly regress wage educ exper iq fem
estimates store m3

estimates table m1 m2 m3, b(%8.4f) se stats(r2 N)

/* ── 13. Visualization ────────────────────────────────────────────────────── */
twoway (scatter wage educ, msize(vsmall) mcolor(blue%30)) ///
       (lfit wage educ, lcolor(red) lwidth(medthick)), ///
    title("OLS: Wage vs Education") xtitle("Education (years)") ytitle("Wage")
graph export "04_regression_analysis/ols_scatter.png", replace

di _n "=== REGRESSION ANALYSIS COMPLETE ==="

/* ============================================================================
   SURVIVAL ANALYSIS IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. Generate Survival Data ────────────────────────────────────────────── */
set obs 300

* Simulate Weibull survival times: T ~ Weibull(shape=1.5, scale=5)
gen u         = runiform()
gen true_time = 5 * (-ln(u))^(1/1.5)   /* inverse CDF of Weibull */
gen cens_time = runiform() * 10
gen obs_time  = min(true_time, cens_time)
gen event     = (true_time <= cens_time)

gen age     = rnormal(55, 12)
gen treated = rbinomial(1, 0.5)
gen stage   = cond(runiform() < .6, "Early", "Late")

* Treatment improves survival
replace obs_time = obs_time + treated * rexponential(1)

* Declare survival data
stset obs_time, failure(event==1)

di _n "Survival summary:"
di "N=" _N "  Events=" r(N_fail)

/* ── 2. Kaplan-Meier Estimator ────────────────────────────────────────────── */
di _n "=== KAPLAN-MEIER ==="
sts graph, by(treated) ///
    title("Kaplan-Meier: by Treatment") ///
    xtitle("Time") ytitle("Survival Probability") ///
    legend(label(1 "Control") label(2 "Treated"))
graph export "07_survival_analysis/km_treatment.png", replace

sts list, by(treated) at(1 2 5 8 10)

/* ── 3. Log-Rank Test ─────────────────────────────────────────────────────── */
di _n "=== LOG-RANK TEST ==="
sts test treated, logrank
sts test stage,   logrank

/* ── 4. Cox PH Model ──────────────────────────────────────────────────────── */
di _n "=== COX PH MODEL ==="
stcox age treated i.stage, efron nolog
di _n "Hazard Ratios:"
stcox age treated i.stage, efron nolog hr

* PH assumption test (Schoenfeld residuals)
estat phtest, detail

* Baseline hazard plot
stcurve, survival at1(treated=0) at2(treated=1) ///
    title("Cox PH: Predicted Survival by Treatment") ///
    xtitle("Time") ytitle("Survival")
graph export "07_survival_analysis/cox_survival.png", replace

/* ── 5. Parametric Survival Models ───────────────────────────────────────── */
di _n "=== PARAMETRIC MODELS ==="

* Weibull model
streg age treated i.stage, distribution(weibull) nolog
estat ic

* Exponential model
streg age treated i.stage, distribution(exponential) nolog
estat ic

* Log-Normal model
streg age treated i.stage, distribution(lognormal) nolog
estat ic

* Log-Logistic
streg age treated i.stage, distribution(loglogistic) nolog
estat ic

di "Compare AIC across models; lowest = best fit"

/* ── 6. Stratified Cox Model ──────────────────────────────────────────────── */
di _n "=== STRATIFIED COX MODEL ==="
stcox age treated, strata(stage) nolog
di "Stratified by stage: separate baseline hazard per stratum"

/* ── 7. Time-Varying Covariates (Split episodes) ────────────────────────── */
di _n "=== TIME-VARYING COVARIATES ==="
* Expand dataset to create time-varying treatment variable
gen start = 0
gen stop  = obs_time
* In a real analysis, use st_ep_split or stsplit for TVC
* Example: treatment changes at time 3
stsplit period, at(3)
gen tvc_treat = treated * (period >= 3)
stcox age tvc_treat, nolog efron
di "Coefficient on TVC treatment: " _b[tvc_treat]

/* ── 8. Competing Risks ───────────────────────────────────────────────────── */
di _n "=== COMPETING RISKS ==="
* Create a second event type
clear
set obs 300
set seed 42
gen u = runiform()
gen time2 = 5*(-ln(u))^(1/1.5)
gen cens2  = runiform()*10
gen time_obs = min(time2, cens2)
gen cause    = cond(runiform()<.6, 1, 2) if time2<=cens2
replace cause = 0 if missing(cause)
gen treated2 = rbinomial(1, .5)

stset time_obs, failure(cause==1)
sts graph, failure by(treated2) ///
    title("Cumulative Incidence: Event of Interest") ///
    xtitle("Time") ytitle("Cumulative Incidence")
graph export "07_survival_analysis/cumulative_incidence.png", replace

* Competing risks regression (Fine-Gray model)
* Requires stcomprisk (ssc install stcomprisk)
* stcomprisk, compet1(2) covariates(treated2) ...

/* ── 9. Nelson-Aalen Cumulative Hazard ───────────────────────────────────── */
clear
set obs 300; set seed 42
gen u=runiform(); gen t=5*(-ln(u))^(1/1.5)
gen c=runiform()*10; gen obs_t=min(t,c); gen ev=(t<=c)
gen treat=rbinomial(1,.5)
stset obs_t, failure(ev==1)

sts graph, cumhaz by(treat) ///
    title("Nelson-Aalen Cumulative Hazard") ///
    xtitle("Time") ytitle("Cumulative Hazard")
graph export "07_survival_analysis/cumhaz.png", replace

di _n "=== SURVIVAL ANALYSIS COMPLETE ==="

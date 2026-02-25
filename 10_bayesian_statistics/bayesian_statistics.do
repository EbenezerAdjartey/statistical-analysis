/* ============================================================================
   BAYESIAN STATISTICS IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   Note: Uses Stata's bayes: prefix (Stata 15+) and manual calculations
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. Bayes Theorem Illustration ──────────────────────────────────────── */
di _n "=== 1. BAYES THEOREM: MEDICAL TEST ==="
local P_D     = 0.01
local P_pos_D  = 0.95
local P_pos_nD = 0.10
local P_nD    = 1 - `P_D'
local P_pos   = `P_pos_D'*`P_D' + `P_pos_nD'*`P_nD'
local P_D_pos = (`P_pos_D' * `P_D') / `P_pos'

di "P(disease)          = " `P_D'
di "P(+|disease)        = " `P_pos_D'  " (sensitivity)"
di "P(+|no disease)     = " `P_pos_nD' " (false positive)"
di "P(+) total          = " %6.4f `P_pos'
di "P(disease|+) POST   = " %6.4f `P_D_pos'
di "Only " %4.1f `P_D_pos'*100 "% of positives actually have disease"

/* ── 2. Beta-Binomial Conjugate Prior ────────────────────────────────────── */
di _n "=== 2. BETA-BINOMIAL ==="
local a0=2; local b0=2   /* prior */
local k=12; local n=20   /* data: 12 successes in 20 trials */
local a1 = `a0' + `k'
local b1 = `b0' + `n' - `k'

di "Prior: Beta(" `a0' "," `b0' ")  mean=" `a0'/(`a0'+`b0')
di "Data: " `k' "/" `n' "  MLE=" `k'/`n'
di "Posterior: Beta(" `a1' "," `b1' ")"
di "  Posterior mean = " `a1'/(`a1'+`b1')
di "  Posterior mode = " (`a1'-1)/(`a1'+`b1'-2)
di "  95% CrI = (" invbeta(`a1',`b1',.025) "," invbeta(`a1',`b1',.975) ")"

* Plot prior, likelihood, posterior
clear; set obs 200
gen theta = _n/200
gen prior_dens = exp(lngamma(`a0'+`b0') - lngamma(`a0') - lngamma(`b0')) * ///
                 theta^(`a0'-1) * (1-theta)^(`b0'-1)
gen post_dens  = exp(lngamma(`a1'+`b1') - lngamma(`a1') - lngamma(`b1')) * ///
                 theta^(`a1'-1) * (1-theta)^(`b1'-1)

twoway (line prior_dens theta, lp(dash) lc(blue)) ///
       (line post_dens  theta, lc(red)), ///
    title("Bayesian Updating: Prior -> Posterior") ///
    xtitle("theta") ytitle("Density") ///
    legend(order(1 "Prior Beta(2,2)" 2 "Posterior Beta(14,10)"))
graph export "10_bayesian_statistics/bayesian_updating.png", replace

/* ── 3. Normal-Normal Conjugate ──────────────────────────────────────────── */
di _n "=== 3. NORMAL-NORMAL CONJUGATE ==="
local mu0=5; local tau2=4; local sigma2=9
local xbar=7.45; local n_obs=8

local tau_n = 1/(1/`tau2' + `n_obs'/`sigma2')
local mu_n  = `tau_n' * (`mu0'/`tau2' + `n_obs'*`xbar'/`sigma2')
local ci_lo = `mu_n' - 1.96*sqrt(`tau_n')
local ci_hi = `mu_n' + 1.96*sqrt(`tau_n')

di "Prior: N(" `mu0' "," `tau2' ")"
di "Data: n=" `n_obs' "  xbar=" `xbar'
di "Posterior mean  = " %7.4f `mu_n'
di "Posterior var   = " %7.4f `tau_n'
di "95% CrI:          (" %6.4f `ci_lo' "," %6.4f `ci_hi' ")"

local freq_lo = `xbar' - 1.96*sqrt(`sigma2'/`n_obs')
local freq_hi = `xbar' + 1.96*sqrt(`sigma2'/`n_obs')
di "Frequentist 95% CI: (" %6.4f `freq_lo' "," %6.4f `freq_hi' ")"

/* ── 4. Bayesian Linear Regression via Stata's bayes: prefix ─────────────── */
di _n "=== 4. BAYESIAN LINEAR REGRESSION ==="
clear; set obs 100; set seed 42
gen x = rnormal()
gen y = 2 + 1.5*x + rnormal()

* OLS for comparison
regress y x
di "OLS: intercept=" _b[_cons] "  slope=" _b[x]

* Bayesian regression (flat prior by default)
bayes, rseed(42) mcmcsize(5000): regress y x
di "Bayesian: posterior means shown above"

* Custom priors: normal(0,10) on coefficients
bayes, rseed(42) mcmcsize(5000) ///
    prior({y:x},    normal(0,10)) ///
    prior({y:_cons},normal(0,10)) ///
    prior({y:sigma2},igamma(0.001,0.001)): ///
    regress y x
bayesstats summary

/* ── 5. MCMC Diagnostics ─────────────────────────────────────────────────── */
di _n "=== 5. MCMC DIAGNOSTICS ==="
bayesgraph diagnostics {y:x} {y:_cons}
graph export "10_bayesian_statistics/mcmc_diagnostics.png", replace

bayesstats grubin     /* Gelman-Rubin convergence */
bayesstats ess        /* Effective sample size */
bayesstats ic         /* DIC information criterion */

/* ── 6. Bayesian Logistic Regression ─────────────────────────────────────── */
di _n "=== 6. BAYESIAN LOGISTIC REGRESSION ==="
gen p   = invlogit(-2 + 0.5*x)
gen grp = rbinomial(1, p)

bayes, rseed(42) mcmcsize(5000): logit grp x
bayesstats summary
bayesgraph diagnostics {grp:x}
graph export "10_bayesian_statistics/bayes_logit_trace.png", replace

/* ── 7. Bayes Factor (manual computation) ────────────────────────────────── */
di _n "=== 7. BAYES FACTOR ==="
di "H0: p=0.5  vs H1: p~U(0,1)"
di "Data: 15 successes in 20 trials"
local k=15; local n_bf=20
* Marginal likelihood under H0: Binomial(20,0.5)
local m0 = exp(lncomb(`n_bf',`k') + `k'*log(0.5) + (`n_bf'-`k')*log(0.5))
* Marginal likelihood under H1: integral = 1/(n+1) (Beta-Binomial with a=b=1)
local m1 = 1/(`n_bf'+1)   /* exact for uniform prior */
local BF = `m1'/`m0'
di "M0=" %8.6f `m0' "  M1=" %8.6f `m1'
di "Bayes Factor BF_10 = " %6.3f `BF'
di "Interpretation: " cond(`BF'>100,"Decisive",cond(`BF'>10,"Strong",cond(`BF'>3,"Moderate","Anecdotal"))) " evidence for H1"

di _n "=== BAYESIAN STATISTICS COMPLETE ==="

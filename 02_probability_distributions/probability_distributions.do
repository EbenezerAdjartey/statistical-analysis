/* ============================================================================
   PROBABILITY DISTRIBUTIONS IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   Topics: Discrete & continuous distributions, PDF/PMF, CDF,
           quantile functions, distribution fitting, GOF tests
   ============================================================================ */

clear all
set more off
set seed 42

di "========================================"
di "   DISCRETE PROBABILITY DISTRIBUTIONS   "
di "========================================"

/* --------------------------------------------------------------------------
   1. Binomial Distribution B(10, 0.5)
   -------------------------------------------------------------------------- */
di _n "--- 1. BINOMIAL B(n=10, p=0.5) ---"
di "Scenario: 10 fair coin flips"

clear; set obs 11
gen k   = _n - 1
gen pmf = binomialp(10, k, 0.5)
gen cdf = binomial(10, k, 0.5)

list k pmf cdf, noobs divider

di "Mean    = " 10 * 0.5
di "Variance= " 10 * 0.5 * 0.5
di "P(X<=5) = " binomial(10, 5, 0.5)
di "P(X>=7) = " 1 - binomial(10, 6, 0.5)
di "90th pct= " invbinomial(10, 0.5, 0.90)

/* --------------------------------------------------------------------------
   2. Poisson Distribution (lambda=3)
   -------------------------------------------------------------------------- */
di _n "--- 2. POISSON(lambda=3) ---"
di "Scenario: 3 customers/minute"

clear; set obs 13
gen k   = _n - 1
gen pmf = poissonp(3, k)
gen cdf = poisson(3, k)
list k pmf cdf, noobs divider

di "P(X=0)=" poissonp(3,0) "  P(X=3)=" poissonp(3,3) "  P(X>5)=" 1-poisson(3,5)
di "Mean=Variance=3"

/* --------------------------------------------------------------------------
   3. Geometric Distribution (p=0.3)
   -------------------------------------------------------------------------- */
di _n "--- 3. GEOMETRIC(p=0.3) ---"
di "Trials until first success"
di "P(success on trial 1) = " 0.3
di "P(success on trial 3) = " 0.7^2 * 0.3
di "P(X<=5)               = " 1 - 0.7^5
di "Mean (expected trials)= " 1/0.3

/* --------------------------------------------------------------------------
   4. Hypergeometric (N=25, K=10, n=8)
   -------------------------------------------------------------------------- */
di _n "--- 4. HYPERGEOMETRIC(N=25, K=10, n=8) ---"
di "Urn: 10 red + 15 blue balls; draw 8 without replacement"
di "P(exactly 3 red) = " hypergeometricp(25, 8, 10, 3)
di "P(at least 4 red)= " 1 - hypergeometric(25, 8, 10, 3)
di "P(exactly 4 red) = " hypergeometricp(25, 8, 10, 4)
di "Mean = " 8 * 10/25

/* --------------------------------------------------------------------------
   5. Negative Binomial (r=3, p=0.4)
   -------------------------------------------------------------------------- */
di _n "--- 5. NEGATIVE BINOMIAL(r=3, p=0.4) ---"
di "Number of failures before 3rd success"
local r = 3; local p = 0.4
* P(k failures) = C(k+r-1,k) * p^r * (1-p)^k
di "P(5 failures before 3rd success) = " exp(lncomb(5+3-1,5)) * 0.4^3 * 0.6^5
di "Mean failures = " `r' * (1-`p') / `p'

di _n "========================================"
di "  CONTINUOUS PROBABILITY DISTRIBUTIONS  "
di "========================================"

/* --------------------------------------------------------------------------
   6. Normal Distribution N(100, 15^2)
   -------------------------------------------------------------------------- */
di _n "--- 6. NORMAL N(100, 15^2) ---"
di "P(X<115)       = " normal((115-100)/15)
di "P(85<X<115)    = " normal((115-100)/15) - normal((85-100)/15)
di "P(X>130)       = " 1 - normal((130-100)/15)
di "95th percentile= " 100 + 15*invnormal(0.95)
di "z-score for 120= " (120-100)/15

di _n "68-95-99.7 Rule:"
di "  P(mu+-1sd)=" normal(1) - normal(-1)
di "  P(mu+-2sd)=" normal(2) - normal(-2)
di "  P(mu+-3sd)=" normal(3) - normal(-3)

/* --------------------------------------------------------------------------
   7. t-Distribution
   -------------------------------------------------------------------------- */
di _n "--- 7. t-DISTRIBUTION (alpha=0.05, two-tail) ---"
foreach df of numlist 1 5 10 30 100 {
    di "  df=" `df' ":  t_crit=" invttail(`df', 0.025)
}

/* --------------------------------------------------------------------------
   8. Chi-square Distribution
   -------------------------------------------------------------------------- */
di _n "--- 8. CHI-SQUARE (95th percentiles) ---"
foreach df of numlist 1 3 5 10 20 30 {
    di "  df=" `df' ":  chi2_crit=" invchi2(`df', 0.95)
}

/* --------------------------------------------------------------------------
   9. F-Distribution
   -------------------------------------------------------------------------- */
di _n "--- 9. F-DISTRIBUTION (alpha=0.05) ---"
foreach d1 of numlist 2 3 5 {
    foreach d2 of numlist 10 20 60 {
        di "  F(" `d1' "," `d2' ") = " invFtail(`d1', `d2', 0.05)
    }
}

/* --------------------------------------------------------------------------
   10. Exponential Distribution (rate=0.5, mean=2)
   -------------------------------------------------------------------------- */
di _n "--- 10. EXPONENTIAL(rate=0.5, mean=2) ---"
di "P(X>3)  = " exp(-0.5*3)
di "P(X<1)  = " 1 - exp(-0.5*1)
di "Median  = " ln(2)/0.5
di "Note: memoryless property: P(X>s+t|X>s) = P(X>t)"

/* --------------------------------------------------------------------------
   11. Gamma Distribution (shape=3, rate=0.5)
   -------------------------------------------------------------------------- */
di _n "--- 11. GAMMA(shape=3, rate=0.5) ---"
di "Mean    = " 3/0.5
di "Variance= " 3/0.5^2
di "Note: Sum of 3 exponential(rate=0.5) r.v.s"

/* --------------------------------------------------------------------------
   12. Weibull Distribution (shape=2, scale=10)
   -------------------------------------------------------------------------- */
di _n "--- 12. WEIBULL(shape=2, scale=10) ---"
di "P(failure before t=8) = " 1 - exp(-(8/10)^2)
di "P(failure before t=10)= " 1 - exp(-(10/10)^2)
di "Median survival time  = " 10 * (-ln(0.5))^(1/2)

/* --------------------------------------------------------------------------
   13. Beta Distribution (alpha=2, beta=5) -- approximate with simulation
   -------------------------------------------------------------------------- */
di _n "--- 13. BETA(alpha=2, beta=5) ---"
clear; set obs 10000
gen x = rbeta(2, 5)
summarize x, detail
di "Theoretical mean  = " 2/7
di "Theoretical var   = " (2*5)/((7^2)*8)

/* --------------------------------------------------------------------------
   14. Uniform Distribution U(0, 10)
   -------------------------------------------------------------------------- */
di _n "--- 14. UNIFORM U(0, 10) ---"
di "P(2<X<7) = " (7-2)/(10-0)
di "Mean     = " (0+10)/2
di "Variance = " (10-0)^2/12

/* --------------------------------------------------------------------------
   15. Distribution Fitting & Goodness-of-Fit
   -------------------------------------------------------------------------- */
di _n "========================================"
di "  DISTRIBUTION FITTING & GOF TESTS      "
di "========================================"

clear; set obs 200
* Generate Gamma(shape=2, scale=2) sample: rgamma(shape, scale)
gen x = rgamma(2, 2)
label var x "Gamma(2,2) sample"

di _n "--- Sample Statistics ---"
summarize x, detail

di _n "--- Shapiro-Wilk Normality Test ---"
swilk x

di _n "--- Skewness & Kurtosis Test ---"
sktest x

* Histogram with normal overlay
histogram x, normal kdensity ///
    title("Gamma Sample: Normal & KDE Overlay") xtitle("x") ytitle("Density")
graph export "02_probability_distributions/histogram_fit.png", replace

* Q-Q plot
qnorm x, title("Q-Q Plot vs Normal Distribution")
graph export "02_probability_distributions/qqplot_normal.png", replace

* Binomial PMF bar chart
clear; set obs 11; gen k=_n-1; gen pmf=binomialp(10,k,0.5)
label var pmf "P(X=k)"
graph bar pmf, over(k) ///
    title("Binomial PMF: B(10, 0.5)") ytitle("P(X=k)")
graph export "02_probability_distributions/binomial_pmf.png", replace

* Poisson PMF
clear; set obs 13; gen k=_n-1; gen pmf=poissonp(3,k)
graph bar pmf, over(k) ///
    title("Poisson PMF: lambda=3") ytitle("P(X=k)")
graph export "02_probability_distributions/poisson_pmf.png", replace

di _n "=== PROBABILITY DISTRIBUTIONS COMPLETE ==="

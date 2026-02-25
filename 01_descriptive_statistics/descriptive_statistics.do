/* ============================================================================
   DESCRIPTIVE STATISTICS IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   Topics: Central tendency, dispersion, skewness, kurtosis,
           frequency tables, cross-tabulations, visualizations
   ============================================================================ */

clear all
set more off
set seed 42

/* --------------------------------------------------------------------------
   1. Generate Synthetic Dataset
   -------------------------------------------------------------------------- */
set obs 200
gen age    = max(18, round(rnormal(35, 10)))
gen income = round(exp(rnormal(10, 0.5)))
gen score  = round(rnormal(70, 15))

gen edu_r  = runiform()
gen education = "Primary"   if edu_r < 0.20
replace education = "Secondary" if edu_r >= 0.20 & edu_r < 0.60
replace education = "Tertiary"  if edu_r >= 0.60

gen gender = cond(runiform() < 0.5, "Male", "Female")

label var age       "Age (years)"
label var income    "Monthly income"
label var score     "Test score"
label var education "Highest education level"
label var gender    "Gender"

di "=== DATASET: `=_N' observations ==="
list in 1/6

/* --------------------------------------------------------------------------
   2. Central Tendency
   -------------------------------------------------------------------------- */
di _n "=== CENTRAL TENDENCY ==="
summarize age income score

di _n "Detailed measures:"
foreach v of varlist age income score {
    quietly summarize `v', detail
    di "`v': mean=" %8.2f r(mean) "  median=" %8.2f r(p50)
}
di "(Mode = most frequent value -- see tabulate below)"

/* --------------------------------------------------------------------------
   3. Dispersion
   -------------------------------------------------------------------------- */
di _n "=== DISPERSION ==="
foreach v of varlist age income score {
    quietly summarize `v', detail
    di "`v':"
    di "  Variance:  " %10.2f r(Var)
    di "  Std Dev:   " %10.2f r(sd)
    di "  Range:     " %10.2f r(max) - r(min)
    di "  IQR:       " %10.2f r(p75) - r(p25)
    di "  CV (%):    " %10.2f 100*r(sd)/r(mean)
    di "  Q1/Q3:     " %6.2f r(p25) "/" %6.2f r(p75)
}

/* --------------------------------------------------------------------------
   4. Skewness & Kurtosis
   -------------------------------------------------------------------------- */
di _n "=== SKEWNESS & KURTOSIS ==="
summarize age income score, detail
foreach v of varlist age income score {
    quietly summarize `v', detail
    di sprintf("`v': skewness=%6.3f  kurtosis=%6.3f  (excess=%6.3f)",
               r(skewness), r(kurtosis), r(kurtosis)-3)
}

/* --------------------------------------------------------------------------
   5. Frequency Tables
   -------------------------------------------------------------------------- */
di _n "=== FREQUENCY TABLE: EDUCATION ==="
tabulate education, sort

di _n "=== FREQUENCY TABLE: GENDER ==="
tabulate gender

/* --------------------------------------------------------------------------
   6. Cross-Tabulation
   -------------------------------------------------------------------------- */
di _n "=== CROSS-TAB: GENDER x EDUCATION ==="
tabulate gender education, chi2 row col

/* --------------------------------------------------------------------------
   7. Grouped Statistics
   -------------------------------------------------------------------------- */
di _n "=== SCORE BY EDUCATION ==="
tabstat score, by(education) stat(n mean sd median min max) nototal

di _n "=== INCOME BY GENDER ==="
tabstat income, by(gender) stat(n mean sd median p25 p75) nototal

/* --------------------------------------------------------------------------
   8. Correlations
   -------------------------------------------------------------------------- */
di _n "=== PEARSON CORRELATIONS ==="
correlate age income score

di _n "=== SPEARMAN CORRELATIONS ==="
spearman age income score

/* --------------------------------------------------------------------------
   9. Visualizations
   -------------------------------------------------------------------------- */
* Histogram with normal + KDE overlay
histogram age, bin(20) normal kdensity ///
    title("Age Distribution") xtitle("Age") ytitle("Density") ///
    note("Normal curve and kernel density estimate shown")
graph export "01_descriptive_statistics/histogram_age.png", replace

* Boxplot of score by education
graph box score, over(education, sort(1)) ///
    title("Test Score by Education Level") ///
    ytitle("Score") note("Boxes = IQR; whiskers = 1.5xIQR")
graph export "01_descriptive_statistics/boxplot_score_edu.png", replace

* Boxplot of income by gender
graph box income, over(gender) ///
    title("Income Distribution by Gender") ytitle("Income")
graph export "01_descriptive_statistics/boxplot_income_gender.png", replace

* Q-Q plot for normality assessment
qnorm score, title("Q-Q Plot: Score vs Normal Distribution")
graph export "01_descriptive_statistics/qqplot_score.png", replace

* Scatter matrix
graph matrix age income score, half ///
    title("Scatter Plot Matrix")
graph export "01_descriptive_statistics/scatter_matrix.png", replace

di _n "=== DESCRIPTIVE STATISTICS COMPLETE ==="

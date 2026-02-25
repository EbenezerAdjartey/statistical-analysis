/* ============================================================================
   HYPOTHESIS TESTING IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. One-Sample t-Test ──────────────────────────────────────────────────── */
di _n "=== 1. ONE-SAMPLE t-TEST ==="
di "H0: mu = 70  vs  H1: mu != 70"
clear; set obs 30
gen scores = rnormal(72, 12)
ttest scores == 70
di "Verdict: " cond(r(p)<0.05,"Reject H0","Fail to reject H0")

/* ── 2. Two-Sample t-Test ──────────────────────────────────────────────────── */
di _n "=== 2. TWO-SAMPLE t-TEST ==="
clear; set obs 80
gen group = cond(_n<=40, 0, 1)
gen score = rnormal(75,10) if group==0
replace score = rnormal(70,12) if group==1

* Levene's test
robvar score, by(group)

* Student's t (equal var)
ttest score, by(group)

* Welch's t (unequal var)
ttest score, by(group) welch

/* ── 3. Paired t-Test ──────────────────────────────────────────────────────── */
di _n "=== 3. PAIRED t-TEST ==="
clear; set obs 25
gen before = rnormal(120, 15)
gen after  = before - rnormal(8, 5)
gen diff   = before - after
di "Mean difference: " %6.3f r(mean)
ttest before == after
* Or: ttest diff == 0

/* ── 4. One-Way ANOVA ──────────────────────────────────────────────────────── */
di _n "=== 4. ONE-WAY ANOVA ==="
clear; set obs 90
gen group_id = ceil(_n/30)
label define grp 1 "G1" 2 "G2" 3 "G3"
label values group_id grp
gen y = rnormal(70,10) if group_id==1
replace y = rnormal(75,10) if group_id==2
replace y = rnormal(80,10) if group_id==3

oneway y group_id, tabulate

/* ── 5. Tukey HSD Post-Hoc ─────────────────────────────────────────────────── */
di _n "=== 5. TUKEY HSD POST-HOC ==="
* In Stata, use pwmean after oneway
pwmean y, over(group_id) mcompare(tukey) effects

/* ── 6. Two-Way ANOVA ──────────────────────────────────────────────────────── */
di _n "=== 6. TWO-WAY ANOVA ==="
clear; set obs 120
gen method = mod(_n-1, 3) + 1  /* 1=A, 2=B, 3=C */
gen gender = cond(mod(_n-1, 2)==0, 1, 2)
gen score  = rnormal(70,10)
replace score = score + 5 if method==2
replace score = score + 10 if method==3

label define meth 1 "A" 2 "B" 3 "C"
label define gend 1 "M" 2 "F"
label values method meth
label values gender gend

anova score method##gender

/* ── 7. Chi-Square Goodness-of-Fit ─────────────────────────────────────────── */
di _n "=== 7. CHI-SQUARE GOODNESS-OF-FIT ==="
clear; set obs 200
gen category = ceil(runiform()*4)
tabulate category

* Compute chi-square manually vs uniform
foreach k of numlist 1/4 {
    quietly count if category==`k'
    local obs_`k' = r(N)
}
di "Expected per category: " 200/4 " = 50"
local chi2_stat = 0
foreach k of numlist 1/4 {
    local chi2_stat = `chi2_stat' + (`obs_`k'' - 50)^2 / 50
}
di "Chi2 statistic = " `chi2_stat'
di "p-value (df=3) = " 1 - chi2(`chi2_stat', 3)

/* ── 8. Chi-Square Test of Independence ────────────────────────────────────── */
di _n "=== 8. CHI-SQUARE INDEPENDENCE TEST ==="
clear; set obs 100
gen gender2 = cond(runiform()<0.5, "Male", "Female")
gen response = cond(runiform()<0.6, "Yes", "No")
tabulate gender2 response, chi2 row col expected

/* ── 9. Fisher's Exact Test ─────────────────────────────────────────────────── */
di _n "=== 9. FISHER'S EXACT TEST ==="
tabulate gender2 response, exact

/* ── 10. Z-Test for Proportion ──────────────────────────────────────────────── */
di _n "=== 10. Z-TEST FOR PROPORTIONS ==="
di "H0: p=0.5  vs  H1: p!=0.5  (60 successes in 100 trials)"
local n_obs=100; local x_obs=60; local p0=0.5
local p_hat = `x_obs'/`n_obs'
local se_z  = sqrt(`p0'*(1-`p0')/`n_obs')
local z_stat= (`p_hat'-`p0')/`se_z'
local p_val = 2*(1-normal(abs(`z_stat')))
local ci_lo  = `p_hat' - 1.96*sqrt(`p_hat'*(1-`p_hat')/`n_obs')
local ci_hi  = `p_hat' + 1.96*sqrt(`p_hat'*(1-`p_hat')/`n_obs')
di "p_hat=" `p_hat' "  z=" %6.4f `z_stat' "  p=" %6.4f `p_val'
di "95% CI: (" %6.4f `ci_lo' ", " %6.4f `ci_hi' ")"
di "Verdict: " cond(`p_val'<0.05,"Reject H0","Fail to reject H0")

/* ── 11. F-Test for Variance Equality ──────────────────────────────────────── */
di _n "=== 11. F-TEST FOR VARIANCE EQUALITY ==="
clear; set obs 60
gen grp2 = cond(_n<=30, 1, 2)
gen x2   = rnormal(50,8)  if grp2==1
replace x2 = rnormal(50,12) if grp2==2

sdtest x2, by(grp2)
robvar x2, by(grp2)  /* Levene's test */

/* ── 12. Visualization ──────────────────────────────────────────────────────── */
clear; set obs 90
gen group_v = ceil(_n/30)
gen y_v = rnormal(70,10) if group_v==1
replace y_v = rnormal(75,10) if group_v==2
replace y_v = rnormal(80,10) if group_v==3
label values group_v grp

graph box y_v, over(group_v) ///
    title("One-Way ANOVA: Score by Group") ytitle("Score")
graph export "03_hypothesis_testing/anova_boxplot.png", replace

di _n "=== HYPOTHESIS TESTING COMPLETE ==="

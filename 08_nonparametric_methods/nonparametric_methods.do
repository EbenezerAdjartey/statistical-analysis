/* ============================================================================
   NONPARAMETRIC METHODS IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. Generate Skewed Data ─────────────────────────────────────────────── */
set obs 200
gen group_id = cond(_n<=100, 1, 2)  /* 1=A, 2=B */
gen value    = rexponential(3) if group_id==1
replace value = rexponential(4) if group_id==2

gen x_rank = rnormal(50,10)
gen y_rank = 0.7*x_rank + rnormal(0,8)

gen before = rgamma(3, 0.5)     /* Gamma(shape=3, scale=2) */
gen after  = before*0.8 + rnormal(0,.5)

/* ── 2. Mann-Whitney U Test (Rank-Sum) ───────────────────────────────────── */
di _n "=== MANN-WHITNEY U TEST (Rank-Sum) ==="
ranksum value, by(group_id)
di "(Z statistic and p-value reported above)"

/* ── 3. Wilcoxon Signed-Rank Test (Paired) ───────────────────────────────── */
di _n "=== WILCOXON SIGNED-RANK TEST ==="
signrank before = after
di "Median difference: " %6.3f r(z)

/* ── 4. Kruskal-Wallis Test ──────────────────────────────────────────────── */
di _n "=== KRUSKAL-WALLIS TEST ==="
clear; set obs 90; set seed 42
gen grp  = ceil(_n/30)
gen val  = rexponential(2) if grp==1
replace val = rexponential(3) if grp==2
replace val = rexponential(5) if grp==3

kwallis val, by(grp)
di "(H statistic and p-value above; reject H0 if p<0.05)"

/* ── 5. Post-Hoc Pairwise Tests with Bonferroni ──────────────────────────── */
di _n "=== PAIRWISE RANK-SUM (Bonferroni) ==="
* Pairwise between groups 1-2, 1-3, 2-3
foreach g1 of numlist 1 2 {
    foreach g2 of numlist 2 3 {
        if `g1' < `g2' {
            quietly ranksum val if inlist(grp,`g1',`g2'), by(grp)
            di "G`g1' vs G`g2': z=" r(z) "  p=" 2*normal(-abs(r(z)))
        }
    }
}

/* ── 6. Spearman Rank Correlation ────────────────────────────────────────── */
di _n "=== RANK CORRELATIONS ==="
clear; set obs 60; set seed 42
gen x_r = rnormal(50,10)
gen y_r = 0.7*x_r + rnormal(0,8)

di "Pearson correlation:"
correlate x_r y_r

di _n "Spearman rank correlation:"
spearman x_r y_r

di _n "Kendall's tau-b:"
ktau x_r y_r

/* ── 7. Kolmogorov-Smirnov Tests ─────────────────────────────────────────── */
di _n "=== KOLMOGOROV-SMIRNOV TESTS ==="
clear; set obs 100; set seed 42
gen ga = rexponential(3)
gen gb = rexponential(4)

* One-sample KS vs normal
ksmirnov ga = normal((ga - r(mean))/r(sd))
di "KS vs Normal (using standardized values)"

* Two-sample KS
ksmirnov ga, by(cond(_n<=50,0,1))

/* ── 8. Kernel Density Estimation ───────────────────────────────────────── */
di _n "=== KERNEL DENSITY ESTIMATION ==="
clear; set obs 100; set seed 42
gen ga2 = rexponential(3)
gen gb2 = rexponential(4)

* KDE with default bandwidth
kdensity ga2, title("KDE: Group A") xtitle("Value") ytitle("Density")
graph export "08_nonparametric_methods/kde_groupA.png", replace

* Compare two groups
twoway (kdensity ga2, lcolor(blue)) ///
       (kdensity gb2, lcolor(red) lp(dash)), ///
    title("KDE: Group A vs B") legend(order(1 "A" 2 "B"))
graph export "08_nonparametric_methods/kde_comparison.png", replace

/* ── 9. Nonparametric Tests: Mood's Median ───────────────────────────────── */
di _n "=== MOOD'S MEDIAN TEST ==="
clear; set obs 100; set seed 42
gen g  = cond(_n<=50,1,2)
gen v  = rexponential(3) if g==1
replace v = rexponential(4) if g==2
median v, by(g)

/* ── 10. Bootstrap CI (Stata's bsample) ──────────────────────────────────── */
di _n "=== BOOTSTRAP CONFIDENCE INTERVALS ==="
clear; set obs 100; set seed 42
gen z = rexponential(3)

* Bootstrap median
bootstrap r(p50), reps(1000) seed(42): summarize z, detail
estat bootstrap, all
di "Bootstrap 95% CI for median shown above"

/* ── 11. Permutation Test ────────────────────────────────────────────────── */
di _n "=== PERMUTATION TEST ==="
clear; set obs 80; set seed 42
gen grp2 = cond(_n<=40,0,1)
gen val2 = rnormal(0,1) if grp2==0
replace val2 = rnormal(0.5,1) if grp2==1

* Permutation test using permtest (ssc install permtest if needed)
* permtest val2, by(grp2) reps(1000)
* Alternative: two-sample t-test as reference
ttest val2, by(grp2)
di "(Note: for non-normal data, prefer rank-based or permutation tests)"

di _n "=== NONPARAMETRIC METHODS COMPLETE ==="

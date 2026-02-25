/* ============================================================================
   MULTIVARIATE ANALYSIS IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. Generate Correlated Dataset ─────────────────────────────────────── */
set obs 200
gen f1 = rnormal(); gen f2 = rnormal(); gen f3 = rnormal()

gen reading    = 0.8*f1 + 0.1*f2 + rnormal(0,.3)
gen vocabulary = 0.7*f1 + 0.2*f2 + rnormal(0,.3)
gen math       = 0.1*f1 + 0.8*f2 + rnormal(0,.3)
gen statistics = 0.2*f1 + 0.9*f2 + rnormal(0,.3)
gen memory1    = 0.1*f1 + 0.1*f2 + 0.8*f3 + rnormal(0,.3)
gen memory2    = 0.2*f1 + 0.1*f2 + 0.7*f3 + rnormal(0,.3)
gen noise1     = rnormal(); gen noise2 = rnormal()

* Group variable
gen group = 1 + (f1+f2 > 0.5) + (f1+f2 > 1.5)
label define grp 1 "Low" 2 "Med" 3 "High"
label values group grp

drop f1 f2 f3

correlate reading vocabulary math statistics memory1 memory2
di "Note: reading/vocabulary correlate strongly; so do math/statistics"

/* ── 2. Principal Component Analysis ─────────────────────────────────────── */
di _n "=== PRINCIPAL COMPONENT ANALYSIS ==="
pca reading vocabulary math statistics memory1 memory2 noise1 noise2, mineigen(1) blanks(.3)
estat loadings, cnorm blanks(.3)

* Scree plot
screeplot, yline(1) title("Scree Plot") ci(het)
graph export "09_multivariate_analysis/screeplot.png", replace

* Biplot
biplot reading vocabulary math statistics memory1 memory2, ///
    stretch(1.5) title("PCA Biplot")
graph export "09_multivariate_analysis/pca_biplot.png", replace

* Score the observations
predict pc1 pc2, score
scatter pc1 pc2, by(group) title("PCA Scores by Group")
graph export "09_multivariate_analysis/pca_scores.png", replace

/* ── 3. Exploratory Factor Analysis ─────────────────────────────────────── */
di _n "=== EXPLORATORY FACTOR ANALYSIS ==="
factor reading vocabulary math statistics memory1 memory2 noise1 noise2, ///
    factors(3) ml
rotate, varimax blanks(.3)
di "Factor loadings above (blanked at |0.3|)"

estat common    /* communalities */

/* ── 4. K-Means Cluster Analysis ─────────────────────────────────────────── */
di _n "=== K-MEANS CLUSTER ANALYSIS ==="
* Try k=2,3,4
forvalues k = 2/4 {
    cluster kmeans reading vocabulary math statistics, k(`k') name(km`k') start(random)
    di "k=`k': within-group SS =  (stored in cluster solution)"
}

* Use k=3
cluster kmeans reading vocabulary math statistics memory1 memory2, k(3) name(km3) start(random(42))
cluster list km3
tabulate km3

* Profile clusters
tabstat reading vocabulary math statistics, by(km3) stat(mean)

/* ── 5. Hierarchical Cluster Analysis ────────────────────────────────────── */
di _n "=== HIERARCHICAL CLUSTERING ==="
* Ward's linkage
cluster ward reading vocabulary math statistics memory1 memory2, name(hc_ward) measure(L2)
cluster dendrogram hc_ward, title("Ward Linkage Dendrogram") ///
    cutoff(3) showcount
graph export "09_multivariate_analysis/dendrogram.png", replace

* Cut into 3 clusters
cluster generate hc3 = groups(3), name(hc_ward)
tabulate hc3

/* ── 6. Discriminant Analysis ────────────────────────────────────────────── */
di _n "=== LINEAR DISCRIMINANT ANALYSIS ==="
discrim lda reading vocabulary math statistics, group(group)
estat grsummarize
estat classtable    /* classification table */
estat errorrate     /* error rate */

/* ── 7. MANOVA ───────────────────────────────────────────────────────────── */
di _n "=== MANOVA ==="
manova reading vocabulary math statistics = i.group
di "Pillai's trace test:"
manova reading vocabulary math statistics = i.group, test(pillai)

/* ── 8. Canonical Correlation Analysis ───────────────────────────────────── */
di _n "=== CANONICAL CORRELATION ANALYSIS ==="
canon (reading vocabulary) (math statistics), test(lr)
di "Canonical correlations and tests shown above"

/* ── 9. Summary Plots ─────────────────────────────────────────────────────── */
* Heatmap of correlation matrix
correlate reading vocabulary math statistics memory1 memory2
matrix C = r(C)
heatplot C, color(hcl sequential, reverse) ///
    title("Correlation Heatmap") ///
    aspectratio(1)
graph export "09_multivariate_analysis/corr_heatmap.png", replace

di _n "=== MULTIVARIATE ANALYSIS COMPLETE ==="

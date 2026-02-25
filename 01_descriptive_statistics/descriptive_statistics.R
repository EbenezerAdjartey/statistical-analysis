# =============================================================================
# DESCRIPTIVE STATISTICS IN R
# =============================================================================
# Topics: Central tendency, dispersion, skewness, kurtosis,
#         frequency tables, cross-tabulations, visualizations
# Author: Ebenezer Adjartey
# =============================================================================

if (!requireNamespace("moments", quietly = TRUE)) install.packages("moments")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr",   quietly = TRUE)) install.packages("tidyr")

library(moments); library(ggplot2); library(dplyr); library(tidyr)
set.seed(42)

# ── 1. Synthetic Dataset ─────────────────────────────────────────────────────
n <- 200
df <- data.frame(
  age       = pmax(18, round(rnorm(n, 35, 10))),
  income    = round(rlnorm(n, 10, 0.5)),
  score     = round(rnorm(n, 70, 15)),
  education = sample(c("Primary","Secondary","Tertiary"), n, TRUE, c(.2,.4,.4)),
  gender    = sample(c("Male","Female"), n, TRUE)
)
cat("=== DATASET (first 6 rows) ===\n"); print(head(df))
cat("Dimensions:", nrow(df), "x", ncol(df), "\n\n")

# ── 2. Central Tendency ───────────────────────────────────────────────────────
cat("=== CENTRAL TENDENCY ===\n")
get_mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
for (v in c("age","income","score")) {
  cat(sprintf("%-8s mean=%-10.2f median=%-8.1f mode=%s\n",
      v, mean(df[[v]]), median(df[[v]]), get_mode(df[[v]])))
}

# ── 3. Dispersion ─────────────────────────────────────────────────────────────
cat("\n=== DISPERSION ===\n")
for (v in c("age","income","score")) {
  x <- df[[v]]
  cat(sprintf("%-8s var=%-12.2f sd=%-10.2f IQR=%-8.2f range=%-6.0f CV=%.2f%%\n",
      v, var(x), sd(x), IQR(x), diff(range(x)), sd(x)/mean(x)*100))
}
cat("\nQuantiles (age):\n"); print(quantile(df$age))

# ── 4. Skewness & Kurtosis ────────────────────────────────────────────────────
cat("\n=== SKEWNESS & KURTOSIS ===\n")
for (v in c("age","income","score")) {
  x <- df[[v]]
  cat(sprintf("%-8s skewness=%7.4f  excess_kurtosis=%7.4f\n",
      v, skewness(x), kurtosis(x) - 3))
}
cat("Interpretation: |skew|>1=highly skewed; excess kurtosis>0=leptokurtic\n")

# ── 5. Five-Number Summary ────────────────────────────────────────────────────
cat("\n=== FIVE-NUMBER SUMMARY ===\n")
print(summary(df[, c("age","income","score")]))

# ── 6. Frequency Table ────────────────────────────────────────────────────────
cat("\n=== FREQUENCY TABLE: EDUCATION ===\n")
ft  <- table(df$education)
pct <- prop.table(ft) * 100
print(data.frame(
  Category       = names(ft),
  Frequency      = as.integer(ft),
  Percent        = round(as.numeric(pct), 1),
  Cumulative_Pct = round(cumsum(as.numeric(pct)), 1)
))

# ── 7. Cross-Tabulation ───────────────────────────────────────────────────────
cat("\n=== CROSS-TAB: GENDER x EDUCATION ===\n")
ct <- table(df$gender, df$education)
print(ct)
cat("\nRow proportions (%):\n")
print(round(prop.table(ct, 1) * 100, 1))
chi2 <- chisq.test(ct)
cat("\nChi-square test:\n"); print(chi2)

# ── 8. Grouped Statistics ─────────────────────────────────────────────────────
cat("\n=== SCORE BY EDUCATION ===\n")
print(df %>% group_by(education) %>%
  summarise(n=n(), mean=round(mean(score),2), median=median(score),
            sd=round(sd(score),2), min=min(score), max=max(score), .groups="drop"))

# ── 9. Correlations ───────────────────────────────────────────────────────────
cat("\n=== PEARSON CORRELATION MATRIX ===\n")
print(round(cor(df[, c("age","income","score")]), 3))
cat("\n=== SPEARMAN CORRELATION MATRIX ===\n")
print(round(cor(df[, c("age","income","score")], method="spearman"), 3))

# ── 10. Visualizations ────────────────────────────────────────────────────────
p1 <- ggplot(df, aes(x=age)) +
  geom_histogram(bins=20, fill="steelblue", color="white", alpha=0.8) +
  geom_vline(xintercept=mean(df$age),   color="red",       linetype="dashed", linewidth=1) +
  geom_vline(xintercept=median(df$age), color="darkgreen", linetype="dashed", linewidth=1) +
  labs(title="Age Distribution", subtitle="Red=Mean, Green=Median",
       x="Age", y="Count") + theme_minimal()

p2 <- ggplot(df, aes(x=education, y=score, fill=education)) +
  geom_boxplot(alpha=0.7) +
  geom_jitter(width=0.2, alpha=0.3, size=0.8) +
  labs(title="Score by Education", x="Education", y="Score") +
  theme_minimal() + theme(legend.position="none")

p3 <- ggplot(df, aes(x=score, fill=gender, color=gender)) +
  geom_density(alpha=0.4) +
  labs(title="Score Density by Gender", x="Score", y="Density") +
  theme_minimal()

p4 <- ggplot(df, aes(sample=score)) +
  stat_qq() + stat_qq_line(color="red") +
  labs(title="Q-Q Plot: Score vs Normal",
       x="Theoretical Quantiles", y="Sample Quantiles") +
  theme_minimal()

dir.create("01_descriptive_statistics", showWarnings=FALSE)
ggsave("01_descriptive_statistics/histogram_age.png",        p1, width=7, height=5)
ggsave("01_descriptive_statistics/boxplot_score_edu.png",    p2, width=7, height=5)
ggsave("01_descriptive_statistics/density_score_gender.png", p3, width=7, height=5)
ggsave("01_descriptive_statistics/qqplot_score.png",         p4, width=7, height=5)
cat("\nPlots saved to 01_descriptive_statistics/\n")
cat("=== DESCRIPTIVE STATISTICS COMPLETE ===\n")

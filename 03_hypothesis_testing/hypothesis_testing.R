# =============================================================================
# HYPOTHESIS TESTING IN R
# =============================================================================
# Topics: t-tests, ANOVA, chi-square, z-tests, F-test, multiple comparisons
# Author: Ebenezer Adjartey
# =============================================================================

if (!requireNamespace("ggplot2",       quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",         quietly=TRUE)) install.packages("dplyr")
if (!requireNamespace("multcomp",      quietly=TRUE)) install.packages("multcomp")
if (!requireNamespace("car",           quietly=TRUE)) install.packages("car")
if (!requireNamespace("PropCIs",       quietly=TRUE)) install.packages("PropCIs")

library(ggplot2); library(dplyr); library(multcomp); library(car)
set.seed(42)

# ── 1. One-Sample t-Test ──────────────────────────────────────────────────────
cat("=== 1. ONE-SAMPLE t-TEST ===\n")
cat("H0: mu = 70  vs  H1: mu != 70\n")
scores <- rnorm(30, mean=72, sd=12)
t_test1 <- t.test(scores, mu=70, alternative="two.sided", conf.level=0.95)
print(t_test1)
cat("Verdict:", if(t_test1$p.value<0.05) "Reject H0" else "Fail to reject H0", "\n\n")

# ── 2. Two-Sample t-Test (Independent) ───────────────────────────────────────
cat("=== 2. TWO-SAMPLE t-TEST ===\n")
group_a <- rnorm(40, mean=75, sd=10)
group_b <- rnorm(40, mean=70, sd=12)

# Levene's test for equal variances
lev <- leveneTest(c(group_a, group_b),
                  factor(rep(c("A","B"), each=40)))
cat("Levene's test (equal variance):\n"); print(lev)

# Equal variance (Student's t)
t_eq <- t.test(group_a, group_b, var.equal=TRUE)
cat("\nStudent's t-test (equal var):\n"); print(t_eq)

# Welch's t (unequal variance)
t_w <- t.test(group_a, group_b, var.equal=FALSE)
cat("\nWelch's t-test (unequal var):\n"); print(t_w)

# Effect size (Cohen's d)
pooled_sd <- sqrt((var(group_a)+var(group_b))/2)
cohens_d  <- (mean(group_a)-mean(group_b)) / pooled_sd
cat(sprintf("Cohen's d = %.4f\n\n", cohens_d))

# ── 3. Paired t-Test ─────────────────────────────────────────────────────────
cat("=== 3. PAIRED t-TEST ===\n")
before <- rnorm(25, mean=120, sd=15)
after  <- before - rnorm(25, mean=8, sd=5)
t_pair <- t.test(before, after, paired=TRUE)
print(t_pair)
cat(sprintf("Mean difference: %.3f\n\n", mean(before-after)))

# ── 4. One-Way ANOVA ─────────────────────────────────────────────────────────
cat("=== 4. ONE-WAY ANOVA ===\n")
g1 <- rnorm(30, 70, 10); g2 <- rnorm(30, 75, 10); g3 <- rnorm(30, 80, 10)
y  <- c(g1, g2, g3)
grp <- factor(rep(c("G1","G2","G3"), each=30))

anova1 <- aov(y ~ grp)
print(summary(anova1))

# Eta-squared
ss <- summary(anova1)[[1]]$`Sum Sq`
eta_sq <- ss[1] / sum(ss)
cat(sprintf("Eta-squared (effect size) = %.4f\n\n", eta_sq))

# ── 5. Tukey HSD Post-Hoc ────────────────────────────────────────────────────
cat("=== 5. TUKEY HSD POST-HOC ===\n")
tukey <- TukeyHSD(anova1)
print(tukey)
plot(tukey, las=1, main="Tukey HSD: 95% Family-Wise CI")

# ── 6. Bonferroni Correction ─────────────────────────────────────────────────
cat("\n=== 6. MULTIPLE COMPARISONS: BONFERRONI ===\n")
raw_p <- c(0.01, 0.04, 0.06, 0.12, 0.20)
adj_bonf <- p.adjust(raw_p, method="bonferroni")
adj_holm <- p.adjust(raw_p, method="holm")
adj_bh   <- p.adjust(raw_p, method="BH")  # Benjamini-Hochberg FDR
result <- data.frame(raw_p=raw_p, bonferroni=adj_bonf, holm=adj_holm, BH_FDR=adj_bh)
print(result)

# ── 7. Two-Way ANOVA ─────────────────────────────────────────────────────────
cat("\n=== 7. TWO-WAY ANOVA ===\n")
n2 <- 120
df2 <- data.frame(
  score  = rnorm(n2, 70, 10),
  method = factor(rep(c("A","B","C"), n2/3)),
  gender = factor(rep(c("M","F"), n2/2))
)
df2$score[df2$method=="B"] <- df2$score[df2$method=="B"] + 5
df2$score[df2$method=="C"] <- df2$score[df2$method=="C"] + 10

anova2 <- aov(score ~ method * gender, data=df2)
print(summary(anova2))

# Type III SS (balanced design)
cat("\nType III SS:\n")
print(Anova(anova2, type="III"))

# ── 8. Chi-Square Tests ───────────────────────────────────────────────────────
cat("\n=== 8. CHI-SQUARE TESTS ===\n")
cat("--- Goodness-of-fit ---\n")
observed <- c(45, 60, 55, 40)
expected_p <- c(0.25, 0.25, 0.25, 0.25)  # expected: uniform
gof <- chisq.test(observed, p=expected_p)
print(gof)

cat("\n--- Independence test ---\n")
cont_table <- matrix(c(30,20,15,35), nrow=2,
                      dimnames=list(c("Male","Female"), c("Yes","No")))
ind <- chisq.test(cont_table)
print(ind)
cat("Expected frequencies:\n"); print(ind$expected)

# Fisher's exact test (for small samples)
cat("\nFisher's exact test:\n")
print(fisher.test(cont_table))

# ── 9. Z-Test for Proportions ────────────────────────────────────────────────
cat("\n=== 9. Z-TEST FOR PROPORTIONS ===\n")
# H0: p = 0.50  (two-tailed)
n_z  <- 100; x_z  <- 60; p0 <- 0.5
p_hat <- x_z / n_z
se_z  <- sqrt(p0*(1-p0)/n_z)
z_stat <- (p_hat - p0) / se_z
p_val_z <- 2*(1-pnorm(abs(z_stat)))
ci_z <- c(p_hat - 1.96*sqrt(p_hat*(1-p_hat)/n_z),
          p_hat + 1.96*sqrt(p_hat*(1-p_hat)/n_z))
cat(sprintf("p_hat=%.2f  z=%.4f  p=%.4f  95%%CI=(%.4f,%.4f)\n",
            p_hat, z_stat, p_val_z, ci_z[1], ci_z[2]))

# ── 10. F-Test for Equal Variances ────────────────────────────────────────────
cat("\n=== 10. F-TEST FOR VARIANCE EQUALITY ===\n")
s1 <- rnorm(30, 50, 8); s2 <- rnorm(30, 50, 12)
f_result <- var.test(s1, s2)
print(f_result)
cat(sprintf("F-ratio = %.4f  p = %.4f\n", f_result$statistic, f_result$p.value))

# ── 11. Visualizations ────────────────────────────────────────────────────────
df_plot <- data.frame(score=y, group=grp)
p1 <- ggplot(df_plot, aes(x=group, y=score, fill=group)) +
  geom_boxplot(alpha=0.7) +
  geom_jitter(width=0.2, alpha=0.3, size=0.8) +
  labs(title=sprintf("One-Way ANOVA\nF=%.2f, p=%.4f",
       summary(anova1)[[1]]$`F value`[1],
       summary(anova1)[[1]]$`Pr(>F)`[1]),
       x="Group", y="Score") +
  theme_minimal() + theme(legend.position="none")

dir.create("03_hypothesis_testing", showWarnings=FALSE)
ggsave("03_hypothesis_testing/anova_boxplot.png", p1, width=7, height=5)

cat("\n=== HYPOTHESIS TESTING COMPLETE ===\n")

# =============================================================================
# NONPARAMETRIC METHODS IN R
# =============================================================================
# Topics: Mann-Whitney, Wilcoxon signed-rank, Kruskal-Wallis,
#         Spearman/Kendall, KS test, KDE, bootstrap CI
# Author: Ebenezer Adjartey
# =============================================================================

pkgs <- c("ggplot2","dplyr","tidyr","coin","boot")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(ggplot2); library(dplyr); library(boot)
set.seed(42)

# ── 1. Generate Skewed Data ────────────────────────────────────────────────────
group_a <- rexp(50, rate=1/3)      # Exp(mean=3), right-skewed
group_b <- rexp(50, rate=1/4)      # Exp(mean=4), shifted right
before  <- rgamma(40, shape=3, rate=0.5)
after   <- before * 0.8 + rnorm(40, 0, 0.5)
g1 <- rexp(30, 0.5); g2 <- rexp(30, 0.33); g3 <- rexp(30, 0.2)
x_rank <- rnorm(60, 50, 10)
y_rank <- x_rank * 0.7 + rnorm(60, 0, 8)

cat("Group A: mean=", round(mean(group_a),3), " sd=", round(sd(group_a),3), "\n")
cat("Group B: mean=", round(mean(group_b),3), " sd=", round(sd(group_b),3), "\n\n")

# ── 2. Mann-Whitney U Test (Wilcoxon Rank-Sum) ───────────────────────────────
cat("=== 2. MANN-WHITNEY U TEST ===\n")
mw <- wilcox.test(group_a, group_b, alternative="two.sided", conf.int=TRUE)
print(mw)

# Rank-biserial correlation (effect size)
U <- mw$statistic
r_rb <- 1 - 2*U / (length(group_a)*length(group_b))
cat(sprintf("Rank-biserial r = %.4f  (|r|>0.5 = large effect)\n\n", r_rb))

# ── 3. Wilcoxon Signed-Rank Test (Paired) ───────────────────────────────────
cat("=== 3. WILCOXON SIGNED-RANK TEST ===\n")
wsr <- wilcox.test(before, after, paired=TRUE, conf.int=TRUE)
print(wsr)
cat(sprintf("Median difference: %.4f\n\n", median(before-after)))

# ── 4. Kruskal-Wallis Test ────────────────────────────────────────────────────
cat("=== 4. KRUSKAL-WALLIS TEST ===\n")
all_vals <- c(g1, g2, g3)
all_grps <- factor(rep(c("G1","G2","G3"), each=30))
kw <- kruskal.test(all_vals ~ all_grps)
print(kw)
cat("Verdict:", if(kw$p.value<0.05)"At least one group differs" else "No significant diff", "\n\n")

# Dunn's post-hoc test (pairwise)
cat("Pairwise Mann-Whitney (Holm correction):\n")
pairs_list <- combn(c("G1","G2","G3"), 2, simplify=FALSE)
p_vals <- sapply(pairs_list, function(pair) {
  x <- list(G1=g1, G2=g2, G3=g3)
  wilcox.test(x[[pair[1]]], x[[pair[2]]])$p.value
})
adj_p <- p.adjust(p_vals, method="holm")
for (i in seq_along(pairs_list)) {
  cat(sprintf("  %s vs %s: p_raw=%.4f  p_adj=%.4f\n",
              pairs_list[[i]][1], pairs_list[[i]][2], p_vals[i], adj_p[i]))
}

# ── 5. Rank Correlations ──────────────────────────────────────────────────────
cat("\n=== 5. RANK CORRELATIONS ===\n")
r_pear  <- cor.test(x_rank, y_rank, method="pearson")
r_spear <- cor.test(x_rank, y_rank, method="spearman")
r_kend  <- cor.test(x_rank, y_rank, method="kendall")

cat(sprintf("Pearson  r   = %.4f  p = %.4f\n", r_pear$estimate,  r_pear$p.value))
cat(sprintf("Spearman rho = %.4f  p = %.4f\n", r_spear$estimate, r_spear$p.value))
cat(sprintf("Kendall  tau = %.4f  p = %.4f\n", r_kend$estimate,  r_kend$p.value))
cat("Spearman/Kendall: robust to outliers and monotone non-linear relationships\n\n")

# ── 6. Kolmogorov-Smirnov Test ────────────────────────────────────────────────
cat("=== 6. KOLMOGOROV-SMIRNOV TESTS ===\n")
# One-sample KS: test vs exponential
ks_exp  <- ks.test(group_a, "pexp", rate=1/mean(group_a))
ks_norm <- ks.test(group_a, "pnorm", mean=mean(group_a), sd=sd(group_a))
cat(sprintf("KS vs Exp(rate=1/mean):  D=%.4f  p=%.4f\n", ks_exp$statistic,  ks_exp$p.value))
cat(sprintf("KS vs Normal:            D=%.4f  p=%.4f\n", ks_norm$statistic, ks_norm$p.value))

# Two-sample KS
ks2 <- ks.test(group_a, group_b)
cat(sprintf("2-sample KS (A vs B):    D=%.4f  p=%.4f\n\n", ks2$statistic, ks2$p.value))

# ── 7. Kernel Density Estimation ─────────────────────────────────────────────
cat("=== 7. KERNEL DENSITY ESTIMATION ===\n")
bw_vals <- c(0.2, 0.5, 1.0)
x_eval  <- seq(min(group_a)-1, max(group_a)+1, length.out=300)

kde_list <- lapply(bw_vals, function(bw) {
  d  <- density(group_a, bw=bw, from=min(x_eval), to=max(x_eval), n=300)
  data.frame(x=d$x, y=d$y, bw=paste0("bw=",bw))
})
kde_df <- do.call(rbind, kde_list)

p1 <- ggplot(kde_df, aes(x=x, y=y, color=bw)) +
  geom_line(linewidth=1.2) +
  geom_histogram(data=data.frame(x=group_a), aes(x=x, y=after_stat(density)),
                 inherit.aes=FALSE, bins=15, fill="grey80", alpha=.5) +
  labs(title="KDE: Effect of Bandwidth", x="Value", y="Density") +
  theme_minimal()

p2_df <- rbind(
  data.frame(x=group_a, grp="Group A"),
  data.frame(x=group_b, grp="Group B")
)
p2 <- ggplot(p2_df, aes(x=x, fill=grp, color=grp)) +
  geom_density(alpha=0.4) +
  labs(title="KDE Comparison: Group A vs B", x="Value", y="Density") +
  theme_minimal()

dir.create("08_nonparametric_methods", showWarnings=FALSE)
ggsave("08_nonparametric_methods/kde_bandwidth.png",   p1, width=7, height=5)
ggsave("08_nonparametric_methods/kde_comparison.png",  p2, width=7, height=5)

# ── 8. Bootstrap Confidence Intervals ────────────────────────────────────────
cat("=== 8. BOOTSTRAP CONFIDENCE INTERVALS ===\n")

# Bootstrap CI for median
n_boot <- 5000
boot_med <- replicate(n_boot, median(sample(group_a, replace=TRUE)))
ci_med   <- quantile(boot_med, c(0.025, 0.975))
cat(sprintf("Sample median = %.4f\n", median(group_a)))
cat(sprintf("Bootstrap 95%% CI: (%.4f, %.4f)\n", ci_med[1], ci_med[2]))
cat(sprintf("Bootstrap SE    = %.4f\n\n", sd(boot_med)))

# Bootstrap CI for Spearman correlation
boot_rho <- replicate(n_boot, {
  idx <- sample(length(x_rank), replace=TRUE)
  cor(x_rank[idx], y_rank[idx], method="spearman")
})
ci_rho <- quantile(boot_rho, c(0.025, 0.975))
cat(sprintf("Spearman rho    = %.4f\n", r_spear$estimate))
cat(sprintf("Bootstrap 95%% CI: (%.4f, %.4f)\n\n", ci_rho[1], ci_rho[2]))

# BCa bootstrap using boot package
boot_stat <- function(data, i) median(data[i])
boot_res  <- boot(group_a, boot_stat, R=n_boot)
bca_ci    <- boot.ci(boot_res, type="bca")
cat("BCa Bootstrap CI for median:\n"); print(bca_ci)

cat("\n=== NONPARAMETRIC METHODS COMPLETE ===\n")

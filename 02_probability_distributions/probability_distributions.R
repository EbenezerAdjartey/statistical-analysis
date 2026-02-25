# =============================================================================
# PROBABILITY DISTRIBUTIONS IN R
# =============================================================================
# Topics: Discrete & continuous distributions, PDF/PMF, CDF,
#         quantile functions, distribution fitting, goodness-of-fit
# Author: Ebenezer Adjartey
# =============================================================================

if (!requireNamespace("ggplot2",      quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("MASS",         quietly=TRUE)) install.packages("MASS")
if (!requireNamespace("fitdistrplus", quietly=TRUE)) install.packages("fitdistrplus")
if (!requireNamespace("tidyr",        quietly=TRUE)) install.packages("tidyr")
library(ggplot2); library(MASS); library(fitdistrplus); library(tidyr)
set.seed(42)

# =============================================================================
# PART A: DISCRETE DISTRIBUTIONS
# =============================================================================
cat("========================================\n")
cat("   DISCRETE PROBABILITY DISTRIBUTIONS  \n")
cat("========================================\n\n")

# ── 1. Binomial B(n=10, p=0.5) ───────────────────────────────────────────────
cat("--- 1. BINOMIAL B(n=10, p=0.5) ---\n")
cat("Scenario: 10 fair coin flips\n")
n <- 10; p <- 0.5; k <- 0:n
df_b <- data.frame(k=k, PMF=round(dbinom(k,n,p),4), CDF=round(pbinom(k,n,p),4))
print(df_b, row.names=FALSE)
cat(sprintf("Mean=%.1f  Var=%.2f\n", n*p, n*p*(1-p)))
cat(sprintf("P(X<=5)=%.4f  P(X>=7)=%.4f  90th pct=%d\n\n",
    pbinom(5,n,p), 1-pbinom(6,n,p), qbinom(0.9,n,p)))

# ── 2. Poisson(lambda=3) ──────────────────────────────────────────────────────
cat("--- 2. POISSON(lambda=3) ---\n")
lam <- 3; k <- 0:12
df_p <- data.frame(k=k, PMF=round(dpois(k,lam),4), CDF=round(ppois(k,lam),4))
print(df_p, row.names=FALSE)
cat(sprintf("P(X=0)=%.4f  P(X=3)=%.4f  P(X>5)=%.4f\n",
    dpois(0,lam), dpois(3,lam), 1-ppois(5,lam)))
cat(sprintf("Mean=Var=%.0f\n\n", lam))

# ── 3. Geometric(p=0.3) ───────────────────────────────────────────────────────
cat("--- 3. GEOMETRIC(p=0.3) ---\n")
cat("Trials until first success\n")
p_g <- 0.3; k <- 1:10
# R's dgeom uses number of failures before success: k-1 failures
df_g <- data.frame(Trial_k=k, PMF=round(dgeom(k-1,p_g),4), CDF=round(pgeom(k-1,p_g),4))
print(df_g, row.names=FALSE)
cat(sprintf("Mean=%.2f  P(X<=5)=%.4f\n\n", 1/p_g, pgeom(4,p_g)))

# ── 4. Hypergeometric(N=25, K=10, n=8) ───────────────────────────────────────
cat("--- 4. HYPERGEOMETRIC(N=25, K=10, n=8) ---\n")
cat("Urn: 10 red + 15 blue balls; draw 8 without replacement\n")
k <- 0:8
df_h <- data.frame(k=k, PMF=round(dhyper(k,10,15,8),4), CDF=round(phyper(k,10,15,8),4))
print(df_h, row.names=FALSE)
cat(sprintf("P(exactly 3 red)=%.4f\nP(at least 4 red)=%.4f\nMean=%.3f\n\n",
    dhyper(3,10,15,8), 1-phyper(3,10,15,8), 8*10/25))

# ── 5. Negative Binomial(r=3, p=0.4) ─────────────────────────────────────────
cat("--- 5. NEGATIVE BINOMIAL(r=3, p=0.4) ---\n")
cat("Number of failures before 3rd success\n")
r_nb <- 3; p_nb <- 0.4; k <- 0:15
df_nb <- data.frame(failures_k=k, PMF=round(dnbinom(k,r_nb,p_nb),4))
print(head(df_nb,10), row.names=FALSE)
cat(sprintf("P(5 failures before 3rd success)=%.4f\nMean failures=%.2f\n\n",
    dnbinom(5,r_nb,p_nb), r_nb*(1-p_nb)/p_nb))

# =============================================================================
# PART B: CONTINUOUS DISTRIBUTIONS
# =============================================================================
cat("========================================\n")
cat("  CONTINUOUS PROBABILITY DISTRIBUTIONS  \n")
cat("========================================\n\n")

# ── 6. Normal N(100, 15) ──────────────────────────────────────────────────────
cat("--- 6. NORMAL N(mu=100, sigma=15) ---\n")
cat(sprintf("P(X<115)       = %.4f\n", pnorm(115,100,15)))
cat(sprintf("P(85<X<115)    = %.4f\n", pnorm(115,100,15)-pnorm(85,100,15)))
cat(sprintf("P(X>130)       = %.4f\n", 1-pnorm(130,100,15)))
cat(sprintf("95th percentile= %.2f\n",  qnorm(0.95,100,15)))
cat(sprintf("z-score for 120= %.3f\n\n",(120-100)/15))
cat("68-95-99.7 rule:\n")
for (z in 1:3) cat(sprintf("  P(mu+-%dsd) = %.4f\n", z, pnorm(z)-pnorm(-z)))

# ── 7. t-Distribution ────────────────────────────────────────────────────────
cat("\n--- 7. t-DISTRIBUTION ---\n")
cat("Critical values (alpha=0.05, two-tail):\n")
for (df in c(1,5,10,30,100))
  cat(sprintf("  df=%4d: t_crit=%.4f\n", df, qt(0.975,df)))

# ── 8. Chi-square ─────────────────────────────────────────────────────────────
cat("\n--- 8. CHI-SQUARE (95th percentiles) ---\n")
for (df in c(1,3,5,10,20,30))
  cat(sprintf("  df=%2d: chi2_crit=%.3f\n", df, qchisq(0.95,df)))

# ── 9. F-Distribution ────────────────────────────────────────────────────────
cat("\n--- 9. F-DISTRIBUTION (alpha=0.05) ---\n")
for (d1 in c(2,3,5)) for (d2 in c(10,20,60))
  cat(sprintf("  F(%d,%2d): F_crit=%.3f\n", d1, d2, qf(0.95,d1,d2)))

# ── 10. Exponential(rate=0.5) ────────────────────────────────────────────────
cat("\n--- 10. EXPONENTIAL(rate=0.5, mean=2) ---\n")
cat(sprintf("P(X>3)=%.4f  P(X<1)=%.4f  Median=%.4f\n",
    1-pexp(3,0.5), pexp(1,0.5), qexp(0.5,0.5)))

# ── 11. Gamma(shape=3, rate=0.5) ─────────────────────────────────────────────
cat("\n--- 11. GAMMA(shape=3, rate=0.5) ---\n")
cat(sprintf("Mean=%.1f  Var=%.1f  P(X<5)=%.4f\n",
    3/0.5, 3/0.5^2, pgamma(5,3,0.5)))

# ── 12. Beta(alpha=2, beta=5) ────────────────────────────────────────────────
cat("\n--- 12. BETA(alpha=2, beta=5) ---\n")
cat(sprintf("Mean=%.4f  P(X<0.3)=%.4f  P(0.2<X<0.5)=%.4f\n",
    2/7, pbeta(0.3,2,5), pbeta(0.5,2,5)-pbeta(0.2,2,5)))

# ── 13. Weibull(shape=2, scale=10) ───────────────────────────────────────────
cat("\n--- 13. WEIBULL(shape=2, scale=10) ---\n")
cat(sprintf("P(failure before t=8)=%.4f  Median=%.4f\n",
    pweibull(8,2,10), qweibull(0.5,2,10)))

# ── 14. Uniform U(0, 10) ─────────────────────────────────────────────────────
cat("\n--- 14. UNIFORM U(0,10) ---\n")
cat(sprintf("P(2<X<7)=%.4f  Mean=%.1f  Var=%.4f\n",
    punif(7,0,10)-punif(2,0,10), 5, 100/12))

# =============================================================================
# PART C: DISTRIBUTION FITTING & GOODNESS-OF-FIT
# =============================================================================
cat("\n========================================\n")
cat("  DISTRIBUTION FITTING & GOF TESTS      \n")
cat("========================================\n\n")

# Generate Gamma(shape=2, rate=0.5) sample
samp <- rgamma(200, shape=2, rate=0.5)
cat("Sample: 200 obs from Gamma(shape=2, rate=0.5)\n")
cat(sprintf("mean=%.3f  sd=%.3f  skew=%.3f\n\n",
    mean(samp), sd(samp), moments::skewness(samp)))

cat("--- Fitting Normal ---\n")
fit_n <- fitdist(samp, "norm")
print(summary(fit_n))

cat("\n--- Fitting Gamma ---\n")
fit_g <- fitdist(samp, "gamma")
print(summary(fit_g))

cat(sprintf("\nAIC: Normal=%.2f  Gamma=%.2f  (lower=better)\n\n",
    fit_n$aic, fit_g$aic))

cat("--- Kolmogorov-Smirnov test vs Normal ---\n")
print(ks.test(scale(samp), "pnorm"))

cat("\n--- Shapiro-Wilk normality test (n=50 subsample) ---\n")
print(shapiro.test(samp[1:50]))

# =============================================================================
# VISUALIZATIONS
# =============================================================================

# 1. Discrete PMF comparison
x_k <- 0:12
df_disc <- data.frame(
  k     = rep(x_k, 2),
  PMF   = c(dbinom(x_k,10,0.5), dpois(x_k,3)),
  Dist  = rep(c("Binomial(10,0.5)","Poisson(3)"), each=length(x_k))
)
p1 <- ggplot(df_disc, aes(x=k, y=PMF, fill=Dist)) +
  geom_col(position="dodge", alpha=0.8) +
  labs(title="Discrete PMF Comparison", x="k", y="P(X=k)") + theme_minimal()

# 2. Normal family PDFs
x <- seq(-4,4,0.01)
df_norm <- data.frame(
  x=rep(x,3),
  PDF=c(dnorm(x,0,1), dnorm(x,0,2), dnorm(x,1,0.5)),
  Dist=rep(c("N(0,1)","N(0,2)","N(1,0.5)"), each=length(x))
)
p2 <- ggplot(df_norm, aes(x=x, y=PDF, color=Dist)) +
  geom_line(linewidth=1.2) +
  labs(title="Normal Distribution Family", x="x", y="PDF f(x)") + theme_minimal()

# 3. Distribution fitting
x_plot <- seq(min(samp), max(samp), length.out=300)
df_fit <- data.frame(
  x=rep(x_plot,2),
  Density=c(dnorm(x_plot,fit_n$estimate[1],fit_n$estimate[2]),
            dgamma(x_plot,fit_g$estimate[1],fit_g$estimate[2])),
  Fit=rep(c("Normal","Gamma"), each=length(x_plot))
)
p3 <- ggplot() +
  geom_histogram(data=data.frame(x=samp), aes(x=x, y=after_stat(density)),
                 bins=30, fill="grey80", color="white") +
  geom_line(data=df_fit, aes(x=x, y=Density, color=Fit), linewidth=1.2) +
  labs(title="Distribution Fitting: Gamma-generated Data",
       x="Value", y="Density") + theme_minimal()

# 4. CDF comparison: Normal vs t
x <- seq(-4,4,0.01)
df_cdf <- data.frame(
  x=rep(x,4),
  CDF=c(pnorm(x), pt(x,3), pt(x,10), pt(x,30)),
  Dist=rep(c("Normal","t(3)","t(10)","t(30)"), each=length(x))
)
p4 <- ggplot(df_cdf, aes(x=x, y=CDF, color=Dist)) +
  geom_line(linewidth=1.2) +
  labs(title="CDF: Normal vs t-Distributions", x="x", y="F(x)") + theme_minimal()

dir.create("02_probability_distributions", showWarnings=FALSE)
ggsave("02_probability_distributions/discrete_pmf.png",    p1, width=8, height=5)
ggsave("02_probability_distributions/normal_family.png",   p2, width=8, height=5)
ggsave("02_probability_distributions/dist_fitting.png",    p3, width=8, height=5)
ggsave("02_probability_distributions/cdf_comparison.png",  p4, width=8, height=5)
cat("Plots saved to 02_probability_distributions/\n")
cat("=== PROBABILITY DISTRIBUTIONS COMPLETE ===\n")

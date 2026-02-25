# =============================================================================
# BAYESIAN STATISTICS IN R
# =============================================================================
# Topics: Bayes theorem, conjugate priors, Bayesian regression,
#         MCMC (Metropolis-Hastings), credible intervals, Bayes Factor
# Author: Ebenezer Adjartey
# =============================================================================

pkgs <- c("ggplot2","dplyr","tidyr","coda","BayesFactor","bayesplot")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(ggplot2); library(dplyr); library(tidyr)
set.seed(42)

# ── 1. Bayes Theorem Illustration ─────────────────────────────────────────────
cat("=== 1. BAYES THEOREM: MEDICAL TEST ===\n")
P_D   <- 0.01    # P(disease) - prevalence
P_pos_D  <- 0.95 # sensitivity
P_pos_nD <- 0.10 # false positive rate
P_nD  <- 1 - P_D
P_pos <- P_pos_D*P_D + P_pos_nD*P_nD  # total probability of positive

P_D_pos <- (P_pos_D * P_D) / P_pos    # posterior

cat(sprintf("P(disease)         = %.3f (prior)\n",        P_D))
cat(sprintf("P(+|disease)       = %.3f (sensitivity)\n",  P_pos_D))
cat(sprintf("P(+|no disease)    = %.3f (false pos rate)\n",P_pos_nD))
cat(sprintf("P(+)               = %.4f\n",                 P_pos))
cat(sprintf("P(disease|+)       = %.4f (posterior)\n",     P_D_pos))
cat(sprintf("\nDespite 95%% sensitivity, only %.1f%% of positives have disease!\n", P_D_pos*100))
cat("This is the base-rate fallacy in action.\n\n")

# ── 2. Beta-Binomial Conjugate Prior ─────────────────────────────────────────
cat("=== 2. BETA-BINOMIAL (Bayesian Proportion) ===\n")
alpha_prior <- 2; beta_prior <- 2
successes <- 12; n_trials <- 20

# Posterior: Beta(alpha+k, beta+n-k)
alpha_post <- alpha_prior + successes
beta_post  <- beta_prior  + (n_trials - successes)

# Summaries
post_mean <- alpha_post / (alpha_post + beta_post)
post_mode <- (alpha_post-1) / (alpha_post+beta_post-2)
ci        <- qbeta(c(0.025,0.975), alpha_post, beta_post)

cat(sprintf("Prior: Beta(%d,%d)  mean=%.3f\n", alpha_prior, beta_prior, alpha_prior/(alpha_prior+beta_prior)))
cat(sprintf("Data:  %d/%d  MLE=%.3f\n",         successes, n_trials, successes/n_trials))
cat(sprintf("Posterior: Beta(%d,%d)\n",           alpha_post, beta_post))
cat(sprintf("  Mean  = %.4f\n", post_mean))
cat(sprintf("  Mode  = %.4f\n", post_mode))
cat(sprintf("  95%% Credible Interval: (%.4f, %.4f)\n\n", ci[1], ci[2]))

# Plot
theta <- seq(0,1,length.out=300)
df_prior   <- data.frame(theta=theta, density=dbeta(theta, alpha_prior, beta_prior), curve="Prior")
df_like    <- data.frame(theta=theta, density=dbinom(successes,n_trials,theta)/
                           integrate(function(p)dbinom(successes,n_trials,p),0,1)$value, curve="Likelihood")
df_post    <- data.frame(theta=theta, density=dbeta(theta, alpha_post, beta_post), curve="Posterior")
df_plot    <- bind_rows(df_prior, df_like, df_post)

p1 <- ggplot(df_plot, aes(x=theta, y=density, color=curve)) +
  geom_line(linewidth=1.5) +
  geom_vline(xintercept=ci, linetype="dotted", color="orange", linewidth=1) +
  labs(title="Bayesian Updating: Prior -> Posterior",
       x="theta (proportion)", y="Density", color="") +
  theme_minimal()

dir.create("10_bayesian_statistics", showWarnings=FALSE)
ggsave("10_bayesian_statistics/bayesian_updating.png", p1, width=8, height=5)

# ── 3. Normal-Normal Conjugate ────────────────────────────────────────────────
cat("=== 3. NORMAL-NORMAL CONJUGATE ===\n")
mu0 <- 5.0; tau2 <- 4.0; sigma2 <- 9.0
data_obs <- c(7.2, 6.8, 8.1, 7.5, 6.9, 8.3, 7.0, 7.8)
n_obs <- length(data_obs); x_bar <- mean(data_obs)

tau_n2 <- 1 / (1/tau2 + n_obs/sigma2)
mu_n   <- tau_n2 * (mu0/tau2 + n_obs*x_bar/sigma2)
ci_n   <- qnorm(c(0.025,0.975), mu_n, sqrt(tau_n2))

cat(sprintf("Prior: N(%.1f, %.1f)\n",    mu0, tau2))
cat(sprintf("Data:  n=%d, x_bar=%.4f\n", n_obs, x_bar))
cat(sprintf("Posterior: N(%.4f, %.4f)\n", mu_n, tau_n2))
cat(sprintf("95%% Credible Interval: (%.4f, %.4f)\n\n", ci_n[1], ci_n[2]))

# Compare with frequentist CI
freq_ci <- x_bar + c(-1,1)*1.96*sqrt(sigma2/n_obs)
cat(sprintf("Frequentist 95%% CI:    (%.4f, %.4f)\n", freq_ci[1], freq_ci[2]))
cat("Bayesian CI is shrunk toward prior; frequentist CI centered on sample mean\n\n")

# ── 4. MCMC: Metropolis-Hastings ──────────────────────────────────────────────
cat("=== 4. MCMC: METROPOLIS-HASTINGS ===\n")
true_mu <- 3.0; sigma_mc <- 2.0
data_mc <- rnorm(50, true_mu, sigma_mc)

log_posterior <- function(mu, data, sigma=2, mu_prior=0, sigma_prior=10) {
  sum(dnorm(data, mu, sigma, log=TRUE)) + dnorm(mu, mu_prior, sigma_prior, log=TRUE)
}

n_iter     <- 10000; proposal_sd <- 0.2
chain      <- numeric(n_iter); chain[1] <- 0
accepted   <- 0

for (i in 2:n_iter) {
  proposal  <- chain[i-1] + rnorm(1, 0, proposal_sd)
  log_ratio <- log_posterior(proposal, data_mc) - log_posterior(chain[i-1], data_mc)
  if (log(runif(1)) < log_ratio) { chain[i] <- proposal; accepted <- accepted+1 }
  else chain[i] <- chain[i-1]
}

burn_in <- 2000; samples <- chain[(burn_in+1):n_iter]
acc_rate <- accepted / n_iter

cat(sprintf("n_iter=%d  burn_in=%d  acceptance=%.3f\n", n_iter, burn_in, acc_rate))
cat(sprintf("Posterior mean: %.4f  (true: %.1f)\n",  mean(samples), true_mu))
cat(sprintf("Posterior SD:   %.4f\n",                  sd(samples)))
ci_mc <- quantile(samples, c(0.025,0.975))
cat(sprintf("95%% Credible Interval: (%.4f, %.4f)\n\n", ci_mc[1], ci_mc[2]))

# MCMC diagnostics
cat("Geweke diagnostic (|z|<2 = good):", geweke.diag(as.mcmc(samples))$z, "\n")
cat("Effective sample size:", round(effectiveSize(as.mcmc(samples))), "\n\n")

# Trace + posterior plot
p_trace <- ggplot(data.frame(iter=1:1000, chain=chain[1:1000]),
                   aes(x=iter, y=chain)) +
  geom_line(alpha=0.7) +
  geom_vline(xintercept=burn_in, color="red", linetype="dashed") +
  labs(title="MCMC Trace Plot", x="Iteration", y="mu") +
  theme_minimal()
ggsave("10_bayesian_statistics/mcmc_trace.png", p_trace, width=7, height=4)

p_post <- ggplot(data.frame(mu=samples), aes(x=mu)) +
  geom_histogram(aes(y=after_stat(density)), bins=50, fill="steelblue", alpha=.7) +
  stat_function(fun=dnorm, args=list(mean=mean(samples), sd=sd(samples)),
                color="red", linewidth=1.2) +
  geom_vline(xintercept=true_mu, color="green", linetype="dashed", linewidth=1.2) +
  labs(title="Posterior Distribution (MCMC)", x="mu", y="Density") +
  theme_minimal()
ggsave("10_bayesian_statistics/mcmc_posterior.png", p_post, width=7, height=4)

# ── 5. Bayes Factor ───────────────────────────────────────────────────────────
cat("=== 5. BAYES FACTOR ===\n")
k_bf <- 15; n_bf <- 20  # 15 successes in 20 trials
p0   <- 0.5

# BF10: H1 (p~U[0,1]) vs H0 (p=0.5)
marginal_H0 <- dbinom(k_bf, n_bf, p0)
marginal_H1 <- integrate(function(p) dbinom(k_bf, n_bf, p), 0, 1)$value
BF_10 <- marginal_H1 / marginal_H0
BF_01 <- 1 / BF_10

cat(sprintf("Marginal likelihood H0: %.6f\n", marginal_H0))
cat(sprintf("Marginal likelihood H1: %.6f\n", marginal_H1))
cat(sprintf("Bayes Factor BF_10 = %.4f  (H1 vs H0)\n", BF_10))
cat("Jeffreys scale: >100=Decisive; 10-100=Strong; 3-10=Moderate; 1-3=Anecdotal\n")
cat(sprintf("Interpretation: %s\n\n",
            if(BF_10>100)"Decisive evidence for H1"
            else if(BF_10>10)"Strong evidence for H1"
            else if(BF_10>3)"Moderate evidence for H1"
            else "Anecdotal/no evidence for H1"))

# ── 6. Bayesian Linear Regression (manual) ───────────────────────────────────
cat("=== 6. BAYESIAN LINEAR REGRESSION ===\n")
n_reg <- 100
x_reg <- rnorm(n_reg)
y_reg <- 2 + 1.5*x_reg + rnorm(n_reg)

# Conjugate Normal-Inverse-Gamma prior
X_mat  <- cbind(1, x_reg)
# Posterior mean (with flat prior ~ standard OLS)
beta_post_mean <- solve(t(X_mat)%*%X_mat) %*% t(X_mat) %*% y_reg
beta_post_cov  <- solve(t(X_mat)%*%X_mat) * var(y_reg - X_mat%*%beta_post_mean)

cat("Posterior means (approximately equal to OLS with flat prior):\n")
cat(sprintf("  Intercept: %.4f  (true: 2.0)\n",  beta_post_mean[1]))
cat(sprintf("  Slope:     %.4f  (true: 1.5)\n",  beta_post_mean[2]))

# Compare with OLS
ols_lm <- lm(y_reg ~ x_reg)
cat(sprintf("OLS Intercept: %.4f  Slope: %.4f\n",
            coef(ols_lm)[1], coef(ols_lm)[2]))

cat("\n=== BAYESIAN STATISTICS COMPLETE ===\n")

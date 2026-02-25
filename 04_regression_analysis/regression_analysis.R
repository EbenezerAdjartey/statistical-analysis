# =============================================================================
# REGRESSION ANALYSIS IN R
# =============================================================================
# Topics: Simple/Multiple OLS, diagnostics, logistic, probit, Tobit,
#         Poisson, NegBin, IV/2SLS, quantile regression
# Author: Ebenezer Adjartey
# =============================================================================

pkgs <- c("ggplot2","dplyr","lmtest","sandwich","car","MASS","AER",
          "quantreg","sampleSelection","pscl","stargazer","corrplot")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(ggplot2); library(lmtest); library(sandwich); library(car)
library(MASS); library(AER); library(quantreg)
set.seed(42)

# ── 1. Data Generation ────────────────────────────────────────────────────────
n    <- 300
educ <- sample(8:20, n, replace=TRUE)
exper<- sample(0:30, n, replace=TRUE)
iq   <- rnorm(n, 100, 15)
fem  <- rbinom(n, 1, 0.5)

# True DGP
wage_log <- 10 + 2*educ + 0.5*exper + 0.05*iq - 5*fem + rnorm(n,0,20)
wage     <- exp(wage_log/40)
latent   <- -5 + 0.3*educ + 0.05*exper + rnorm(n)
employed <- as.integer(latent > 0)
pub_rate <- exp(0.5 + 0.1*educ + rnorm(n,0,.3))
pubs     <- rpois(n, pub_rate)

df <- data.frame(wage=wage, educ=educ, exper=exper, iq=iq,
                 female=fem, employed=employed, publications=pubs)
cat("Dataset head:\n"); print(head(df))

# ── 2. Simple OLS ─────────────────────────────────────────────────────────────
cat("\n=== SIMPLE OLS: wage ~ educ ===\n")
ols1 <- lm(wage ~ educ, data=df)
print(summary(ols1))

# ── 3. Multiple OLS ───────────────────────────────────────────────────────────
cat("\n=== MULTIPLE OLS ===\n")
ols2 <- lm(wage ~ educ + exper + iq + female, data=df)
print(summary(ols2))
cat(sprintf("Adjusted R2: %.4f\n", summary(ols2)$adj.r.squared))

# ── 4. OLS Diagnostics ────────────────────────────────────────────────────────
cat("\n=== OLS DIAGNOSTICS ===\n")

# VIF (multicollinearity)
cat("VIF values:\n"); print(vif(ols2))

# Breusch-Pagan heteroskedasticity test
cat("\nBreusch-Pagan test:\n"); print(bptest(ols2))

# Robust standard errors (HC3)
cat("\nHC3 Robust SEs:\n")
print(coeftest(ols2, vcov=vcovHC(ols2, type="HC3")))

# Ramsey RESET test (functional form)
cat("\nRamsey RESET test:\n"); print(resettest(ols2))

# Durbin-Watson (autocorrelation)
cat("\nDurbin-Watson:\n"); print(dwtest(ols2))

# ── 5. Logistic Regression ────────────────────────────────────────────────────
cat("\n=== LOGISTIC REGRESSION ===\n")
logit <- glm(employed ~ educ + exper + female, data=df, family=binomial(link="logit"))
print(summary(logit))

# Odds ratios
cat("\nOdds Ratios (exp(coef)):\n")
print(exp(coef(logit)))
cat("\nOR with 95% CI:\n")
print(exp(confint(logit)))

# McFadden's pseudo-R2
null_ll <- logLik(glm(employed~1, data=df, family=binomial))
mcfadden <- 1 - logLik(logit)/null_ll
cat(sprintf("\nMcFadden R2 = %.4f\n", mcfadden))

# ── 6. Probit Model ────────────────────────────────────────────────────────────
cat("\n=== PROBIT MODEL ===\n")
probit <- glm(employed ~ educ + exper + female, data=df, family=binomial(link="probit"))
print(summary(probit))

# ── 7. Poisson Regression ─────────────────────────────────────────────────────
cat("\n=== POISSON REGRESSION ===\n")
pois <- glm(publications ~ educ + exper, data=df, family=poisson)
print(summary(pois))
cat("\nIncidence Rate Ratios:\n"); print(exp(coef(pois)))

# Test for overdispersion
cat("\nOverdispersion test:\n")
print(dispersiontest(pois))

# ── 8. Negative Binomial Regression ──────────────────────────────────────────
cat("\n=== NEGATIVE BINOMIAL ===\n")
nb <- glm.nb(publications ~ educ + exper, data=df)
print(summary(nb))

# LR test: Poisson vs NB
lr_stat <- 2*(logLik(nb) - logLik(pois))
cat(sprintf("LR test stat=%.4f, p=%.4f (prefer NB if p<0.05)\n",
            lr_stat, pchisq(lr_stat, 1, lower.tail=FALSE)))

# ── 9. IV/2SLS Regression ────────────────────────────────────────────────────
cat("\n=== IV/2SLS REGRESSION ===\n")
# Create an instrument: proximity to college (random)
df$college_prox <- 0.5*educ + rnorm(n)
iv_model <- ivreg(wage ~ educ + exper | college_prox + exper, data=df)
print(summary(iv_model, diagnostics=TRUE))

# ── 10. Quantile Regression ───────────────────────────────────────────────────
cat("\n=== QUANTILE REGRESSION ===\n")
for (tau in c(0.25, 0.50, 0.75)) {
  qr <- rq(wage ~ educ + exper + female, data=df, tau=tau)
  coefs <- coef(qr)
  cat(sprintf("Q(%s): educ=%.4f  exper=%.4f  female=%.4f\n",
              tau, coefs["educ"], coefs["exper"], coefs["female"]))
}

# ── 11. Diagnostic Plots ──────────────────────────────────────────────────────
par(mfrow=c(2,2))
plot(ols2)
par(mfrow=c(1,1))

# Scatter plot
p1 <- ggplot(df, aes(x=educ, y=wage)) +
  geom_point(alpha=0.4, size=1.5) +
  geom_smooth(method="lm", color="red", se=TRUE) +
  labs(title=sprintf("OLS: wage ~ educ  (R2=%.3f)", summary(ols1)$r.squared),
       x="Education", y="Wage") + theme_minimal()

dir.create("04_regression_analysis", showWarnings=FALSE)
ggsave("04_regression_analysis/ols_scatter.png", p1, width=7, height=5)

cat("\n=== REGRESSION ANALYSIS COMPLETE ===\n")

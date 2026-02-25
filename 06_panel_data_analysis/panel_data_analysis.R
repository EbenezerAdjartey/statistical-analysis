# =============================================================================
# PANEL DATA ANALYSIS IN R
# =============================================================================
# Topics: Pooled OLS, Fixed Effects, Random Effects, Hausman, FD, GMM
# Author: Ebenezer Adjartey
# =============================================================================

pkgs <- c("plm","lmtest","sandwich","ggplot2","dplyr","stargazer")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(plm); library(lmtest); library(sandwich); library(ggplot2)
set.seed(42)

# ── 1. Generate Balanced Panel ────────────────────────────────────────────────
N <- 100; T <- 10; n <- N*T
id   <- rep(1:N, each=T)
time <- rep(1:T, N)

alpha_i <- rep(rnorm(N, 0, 2), each=T)  # individual FE
x1 <- rnorm(n) + alpha_i*0.3             # correlated with FE
x2 <- rnorm(n)
y  <- alpha_i + 2*x1 + 1.5*x2 + rnorm(n)

panel <- data.frame(id=id, time=time, y=y, x1=x1, x2=x2)
pdata <- pdata.frame(panel, index=c("id","time"))
cat("Panel dimensions:", N, "individuals x", T, "periods =", n, "obs\n\n")

# ── 2. Pooled OLS ─────────────────────────────────────────────────────────────
cat("=== POOLED OLS ===\n")
pooled <- plm(y ~ x1 + x2, data=pdata, model="pooling")
print(summary(pooled))
cat(sprintf("Pooled OLS: coef(x1)=%.4f (true=2.0)\n\n", coef(pooled)["x1"]))

# ── 3. Fixed Effects (Within Estimator) ──────────────────────────────────────
cat("=== FIXED EFFECTS (WITHIN) ===\n")
fe <- plm(y ~ x1 + x2, data=pdata, model="within")
print(summary(fe))
cat(sprintf("FE: coef(x1)=%.4f  coef(x2)=%.4f\n\n", coef(fe)["x1"], coef(fe)["x2"]))

# ── 4. Random Effects ─────────────────────────────────────────────────────────
cat("=== RANDOM EFFECTS ===\n")
re <- plm(y ~ x1 + x2, data=pdata, model="random")
print(summary(re))
cat(sprintf("RE: coef(x1)=%.4f  coef(x2)=%.4f\n\n", coef(re)["x1"], coef(re)["x2"]))

# ── 5. Hausman Test ───────────────────────────────────────────────────────────
cat("=== HAUSMAN TEST (FE vs RE) ===\n")
hausman <- phtest(fe, re)
print(hausman)
cat("Verdict:", if(hausman$p.value<0.05) "Use Fixed Effects (RE inconsistent)"
               else "Random Effects preferred", "\n\n")

# ── 6. First-Difference Estimator ─────────────────────────────────────────────
cat("=== FIRST-DIFFERENCE ESTIMATOR ===\n")
fd <- plm(y ~ x1 + x2, data=pdata, model="fd")
print(summary(fd))
cat(sprintf("FD: coef(x1)=%.4f  coef(x2)=%.4f\n\n", coef(fd)["x1"], coef(fd)["x2"]))

# ── 7. Time Fixed Effects ─────────────────────────────────────────────────────
cat("=== TWO-WAY FIXED EFFECTS (entity + time) ===\n")
twfe <- plm(y ~ x1 + x2, data=pdata, model="within", effect="twoways")
print(summary(twfe))

# ── 8. Cluster-Robust Standard Errors ────────────────────────────────────────
cat("\n=== CLUSTER-ROBUST SEs ===\n")
cat("Fixed Effects with cluster-robust SEs (by entity):\n")
fe_robust <- coeftest(fe, vcov=vcovHC(fe, type="HC1", cluster="group"))
print(fe_robust)

# ── 9. Tests ──────────────────────────────────────────────────────────────────
cat("\n=== PANEL DIAGNOSTIC TESTS ===\n")
# F-test for individual effects (FE vs Pooled)
cat("F-test (FE vs Pooled):\n")
print(pFtest(fe, pooled))

# Breusch-Pagan LM test (RE vs Pooled)
cat("\nBreusch-Pagan LM test (RE vs Pooled):\n")
print(plmtest(pooled, type="bp"))

# Test for serial correlation
cat("\nBreusch-Godfrey serial correlation test:\n")
print(pbgtest(fe, order=2))

# ── 10. Coefficient Comparison Plot ──────────────────────────────────────────
comp <- data.frame(
  Model    = c("Pooled OLS","Fixed Effects","Random Effects","First Diff"),
  x1_coef  = c(coef(pooled)["x1"], coef(fe)["x1"], coef(re)["x1"], coef(fd)["x1"]),
  x2_coef  = c(coef(pooled)["x2"], coef(fe)["x2"], coef(re)["x2"], coef(fd)["x2"])
)
cat("\nCoefficient Comparison (True: x1=2.0, x2=1.5):\n")
print(round(comp, 4))

library(tidyr)
comp_long <- pivot_longer(comp, c("x1_coef","x2_coef"), names_to="Var", values_to="Estimate")
p1 <- ggplot(comp_long, aes(x=Model, y=Estimate, fill=Var)) +
  geom_col(position="dodge", alpha=0.8) +
  geom_hline(yintercept=2.0, linetype="dashed", color="blue",  alpha=0.6) +
  geom_hline(yintercept=1.5, linetype="dashed", color="red",   alpha=0.6) +
  labs(title="Panel Estimators (dashed lines = true values)",
       x="Model", y="Coefficient") +
  theme_minimal() + theme(axis.text.x=element_text(angle=30, hjust=1))

dir.create("06_panel_data_analysis", showWarnings=FALSE)
ggsave("06_panel_data_analysis/panel_comparison.png", p1, width=9, height=5)

cat("\n=== PANEL DATA ANALYSIS COMPLETE ===\n")

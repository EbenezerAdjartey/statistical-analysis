# =============================================================================
# SURVIVAL ANALYSIS IN R
# =============================================================================
# Topics: Kaplan-Meier, log-rank test, Cox PH model, parametric models,
#         competing risks, time-varying covariates
# Author: Ebenezer Adjartey
# =============================================================================

pkgs <- c("survival","survminer","ggplot2","dplyr","flexsurv","cmprsk")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(survival); library(survminer); library(ggplot2); library(flexsurv)
set.seed(42)

# ── 1. Generate Survival Data ─────────────────────────────────────────────────
n <- 300
true_time <- rweibull(n, shape=1.5, scale=5)
cens_time <- runif(n, 0, 10)
obs_time  <- pmin(true_time, cens_time)
event     <- as.integer(true_time <= cens_time)

age     <- rnorm(n, 55, 12)
treated <- rbinom(n, 1, 0.5)
stage   <- sample(c("Early","Late"), n, replace=TRUE, prob=c(.6,.4))

# Treatment increases survival time
obs_time <- obs_time + treated * rexp(n, rate=1)

df <- data.frame(time=obs_time, event=event, age=age,
                 treated=treated, stage=stage)
df$surv_obj <- with(df, Surv(time, event))

cat("N =", n, "| Events =", sum(event), "| Censoring rate =",
    round((1-mean(event))*100, 1), "%\n\n")

# ── 2. Kaplan-Meier Estimator ─────────────────────────────────────────────────
cat("=== KAPLAN-MEIER ESTIMATOR ===\n")
km_all <- survfit(Surv(time, event) ~ 1, data=df)
cat("Overall KM summary:\n"); print(summary(km_all, times=c(1,2,5,8,10)))
cat("Median survival:", km_all$time[min(which(km_all$surv<=0.5))], "\n\n")

km_trt <- survfit(Surv(time, event) ~ treated, data=df)
print(km_trt)

# Plot
dir.create("07_survival_analysis", showWarnings=FALSE)
p1 <- ggsurvplot(km_trt, data=df,
                  risk.table=TRUE, pval=TRUE, conf.int=TRUE,
                  palette=c("dodgerblue","red"),
                  legend.labs=c("Control","Treated"),
                  title="KM Curves by Treatment",
                  xlab="Time", ylab="Survival Probability")
ggsave("07_survival_analysis/km_curves.png", print(p1), width=8, height=7)

# ── 3. Log-Rank Test ──────────────────────────────────────────────────────────
cat("=== LOG-RANK TEST ===\n")
lr_trt   <- survdiff(Surv(time,event) ~ treated, data=df)
lr_stage <- survdiff(Surv(time,event) ~ stage,   data=df)
print(lr_trt)
print(lr_stage)

# ── 4. Cox Proportional Hazards Model ────────────────────────────────────────
cat("\n=== COX PH MODEL ===\n")
cox <- coxph(Surv(time,event) ~ age + treated + stage, data=df)
print(summary(cox))
cat("\nHazard Ratios:\n"); print(exp(coef(cox)))
cat("\nHR with 95% CI:\n"); print(exp(confint(cox)))

# Test PH assumption (Schoenfeld residuals)
cat("\nPH assumption test (Schoenfeld residuals):\n")
print(cox.zph(cox))

# Plot Schoenfeld residuals
png("07_survival_analysis/schoenfeld_resid.png", width=800, height=400)
par(mfrow=c(1,2))
plot(cox.zph(cox)["age"],   main="Schoenfeld: age")
plot(cox.zph(cox)["treated"],main="Schoenfeld: treated")
par(mfrow=c(1,1))
dev.off()

# ── 5. Baseline Survival Function ────────────────────────────────────────────
cat("\n=== BASELINE SURVIVAL ===\n")
bh <- basehaz(cox, centered=FALSE)
cat("Baseline cumulative hazard (first 5):\n"); print(head(bh))

# Predicted survival for specific profiles
new_pts <- data.frame(
  age     = c(50, 70),
  treated = c(1,  0),
  stage   = c("Early","Late")
)
pred_surv <- survfit(cox, newdata=new_pts)
cat("\nPredicted median survival:\n")
print(pred_surv)

# ── 6. Parametric Survival Models ────────────────────────────────────────────
cat("\n=== PARAMETRIC SURVIVAL MODELS ===\n")

# Weibull
weibull_fit <- flexsurvreg(Surv(time,event) ~ age + treated + stage,
                            data=df, dist="weibull")
cat("Weibull AFT model:\n"); print(weibull_fit)

# Exponential
exp_fit <- flexsurvreg(Surv(time,event) ~ age + treated,
                        data=df, dist="exponential")
cat("\nExponential model AIC:", exp_fit$AIC, "\n")
cat("Weibull model AIC:    ", weibull_fit$AIC, "\n")

# Log-Normal
lnorm_fit <- flexsurvreg(Surv(time,event) ~ age + treated + stage,
                          data=df, dist="lognormal")
cat("Log-Normal model AIC: ", lnorm_fit$AIC, "\n")

# ── 7. Stratified Cox Model ───────────────────────────────────────────────────
cat("\n=== STRATIFIED COX MODEL ===\n")
cox_strat <- coxph(Surv(time,event) ~ age + treated + strata(stage), data=df)
print(summary(cox_strat))

# ── 8. Competing Risks ────────────────────────────────────────────────────────
cat("\n=== COMPETING RISKS ===\n")
# Create competing events: 1=event of interest, 2=competing event
cause <- sample(c(0,1,2), n, replace=TRUE, prob=c(.2,.5,.3))
cause[event==0] <- 0  # censored
df$cause <- cause

library(cmprsk)
cg <- cuminc(df$time, df$cause, group=df$treated, cencode=0)
cat("Cumulative incidence (competing risks):\n")
print(cg$Tests)

# ── 9. Summary Visualization ──────────────────────────────────────────────────
p2 <- ggsurvplot(km_all, data=df,
                  fun="cumhaz", conf.int=TRUE,
                  title="Nelson-Aalen Cumulative Hazard",
                  xlab="Time", ylab="Cumulative Hazard")
ggsave("07_survival_analysis/cumhaz.png", print(p2), width=7, height=5)

cat("\n=== SURVIVAL ANALYSIS COMPLETE ===\n")

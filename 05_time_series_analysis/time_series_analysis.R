# =============================================================================
# TIME SERIES ANALYSIS IN R
# =============================================================================
# Topics: Plotting, ACF/PACF, ADF/KPSS stationarity, ARIMA, SARIMA,
#         VAR, VECM (cointegration), ARCH/GARCH
# Author: Ebenezer Adjartey
# =============================================================================

pkgs <- c("ggplot2","forecast","tseries","vars","urca","rugarch","lmtest","zoo")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(ggplot2); library(forecast); library(tseries)
library(vars);    library(urca);     library(rugarch)
set.seed(42)

# ── 1. Generate Time Series ───────────────────────────────────────────────────
n <- 200
t_idx <- 1:n

# AR(2): y_t = 0.6*y_{t-1} - 0.2*y_{t-2} + e_t
e  <- rnorm(n)
y  <- numeric(n)
for (i in 3:n) y[i] <- 0.6*y[i-1] - 0.2*y[i-2] + e[i]
ar2_ts <- ts(y, start=c(2005,1), frequency=12)

# Random walk
rw_ts <- ts(cumsum(rnorm(n)), start=c(2005,1), frequency=12)

# Seasonal series
seas_ts <- ts(0.1*t_idx + 5*sin(2*pi*t_idx/12) + rnorm(n,0,.5),
              start=c(2005,1), frequency=12)

# Financial returns (GARCH-like)
ret <- rnorm(n)
for (i in 3:n) {
  vol <- sqrt(0.01 + 0.1*ret[i-1]^2 + 0.85*var(ret[1:(i-1)]))
  ret[i] <- vol * rnorm(1)
}
ret_ts <- ts(ret, start=c(2005,1), frequency=12)

cat("Series created: AR(2), RandomWalk, Seasonal, Returns\n\n")

# ── 2. Plots ──────────────────────────────────────────────────────────────────
dir.create("05_time_series_analysis", showWarnings=FALSE)
png("05_time_series_analysis/ts_overview.png", width=900, height=700)
par(mfrow=c(2,2))
plot(ar2_ts,  main="AR(2) Process",       ylab="Value")
plot(rw_ts,   main="Random Walk",          ylab="Value")
plot(seas_ts, main="Seasonal Series",      ylab="Value")
plot(ret_ts,  main="Financial Returns",    ylab="Return")
par(mfrow=c(1,1))
dev.off()

# ACF/PACF for AR(2)
png("05_time_series_analysis/acf_pacf.png", width=900, height=400)
par(mfrow=c(1,2))
acf(ar2_ts,  lag.max=24, main="ACF: AR(2)")
pacf(ar2_ts, lag.max=24, main="PACF: AR(2)")
par(mfrow=c(1,1))
dev.off()

# ── 3. Stationarity Tests ─────────────────────────────────────────────────────
cat("=== STATIONARITY TESTS ===\n")

run_tests <- function(x, name) {
  cat(sprintf("\n--- %s ---\n", name))
  adf <- adf.test(x, alternative="stationary")
  cat(sprintf("ADF:  stat=%.4f  p=%.4f  => %s\n",
              adf$statistic, adf$p.value,
              if(adf$p.value<0.05)"Stationary" else "Non-stationary"))
  kp  <- kpss.test(x, null="Level")
  cat(sprintf("KPSS: stat=%.4f  p=%.4f  => %s\n",
              kp$statistic, kp$p.value,
              if(kp$p.value>0.05)"Stationary" else "Non-stationary"))
  pp  <- pp.test(x)
  cat(sprintf("PP:   stat=%.4f  p=%.4f  => %s\n",
              pp$statistic, pp$p.value,
              if(pp$p.value<0.05)"Stationary" else "Non-stationary"))
}

run_tests(ar2_ts,          "AR(2)")
run_tests(rw_ts,           "Random Walk")
run_tests(diff(rw_ts),     "Differenced Random Walk")
run_tests(seas_ts,         "Seasonal")

# ── 4. ARIMA Modeling ─────────────────────────────────────────────────────────
cat("\n=== ARIMA MODELING ===\n")

# Manual: ARIMA(2,0,0)
arima_manual <- Arima(ar2_ts, order=c(2,0,0))
cat("ARIMA(2,0,0):\n"); print(arima_manual)

# Auto ARIMA
cat("\nAuto ARIMA:\n")
auto_fit <- auto.arima(ar2_ts, ic="aic", stepwise=FALSE, approximation=FALSE)
print(summary(auto_fit))

# Forecasting
fc <- forecast(auto_fit, h=12)
cat("\n12-step forecast:\n"); print(fc)

png("05_time_series_analysis/arima_forecast.png", width=800, height=500)
plot(fc, main="ARIMA Forecast (12 months)")
dev.off()

# Ljung-Box test on residuals
cat("\nLjung-Box test on residuals:\n")
print(Box.test(auto_fit$residuals, lag=10, type="Ljung-Box"))

# ── 5. SARIMA ─────────────────────────────────────────────────────────────────
cat("\n=== SARIMA ===\n")
sarima_fit <- auto.arima(seas_ts, seasonal=TRUE, ic="aic")
print(summary(sarima_fit))
cat(sprintf("SARIMA order: (%d,%d,%d)(%d,%d,%d)[%d]\n",
            sarima_fit$arma[1], sarima_fit$arma[6], sarima_fit$arma[2],
            sarima_fit$arma[3], sarima_fit$arma[7], sarima_fit$arma[4],
            sarima_fit$arma[5]))

# ── 6. VAR (Vector Autoregression) ────────────────────────────────────────────
cat("\n=== VAR ===\n")
# Two correlated AR series
e12 <- MASS::mvrnorm(n, c(0,0), matrix(c(1,.7,.7,1),2,2))
y1  <- numeric(n); y2 <- numeric(n)
for (i in 3:n) {
  y1[i] <- 0.5*y1[i-1] + 0.2*y2[i-1] + e12[i,1]
  y2[i] <- 0.1*y1[i-1] + 0.4*y2[i-1] + e12[i,2]
}
var_data <- ts(cbind(y1, y2), start=c(2005,1), frequency=12)

var_select <- VARselect(var_data, lag.max=8, type="const")
cat("Optimal lag by AIC:", var_select$selection["AIC(n)"], "\n")

var_fit <- VAR(var_data, p=var_select$selection["AIC(n)"], type="const")
cat("VAR summary:\n"); print(summary(var_fit))

# Granger causality
cat("\nGranger causality (y2 -> y1?):\n")
print(causality(var_fit, cause="y2")$Granger)

# Impulse Response Function
irf_result <- irf(var_fit, response="y1", impulse="y2", n.ahead=10)
png("05_time_series_analysis/var_irf.png", width=700, height=500)
plot(irf_result, main="IRF: Response of y1 to Shock in y2")
dev.off()

# ── 7. ARCH/GARCH ─────────────────────────────────────────────────────────────
cat("\n=== ARCH/GARCH VOLATILITY ===\n")

# ARCH test (Engle's LM test)
arch_test <- ArchTest(ret_ts, lags=5)
cat("ARCH-LM test:\n"); print(arch_test)

# Fit GARCH(1,1)
spec <- ugarchspec(
  variance.model = list(model="sGARCH", garchOrder=c(1,1)),
  mean.model     = list(armaOrder=c(0,0), include.mean=TRUE),
  distribution.model = "norm"
)
garch_fit <- ugarchfit(spec=spec, data=ret_ts)
cat("\nGARCH(1,1) results:\n"); print(garch_fit)

# Conditional volatility forecast
garch_fc <- ugarchforecast(garch_fit, n.ahead=10)
cat("\nVolatility forecast (next 10 periods):\n")
print(sigma(garch_fc))

cat("\n=== TIME SERIES ANALYSIS COMPLETE ===\n")

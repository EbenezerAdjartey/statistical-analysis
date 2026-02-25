/* ============================================================================
   TIME SERIES ANALYSIS IN STATA
   Author: Ebenezer Adjartey | Stata 16+
   Packages: tsset, arima, var, arch (built-in)
   ============================================================================ */
clear all
set more off
set seed 42

/* ── 1. Generate AR(2) series ─────────────────────────────────────────────── */
local n = 200
set obs `n'
gen t = _n
tsset t

* AR(2): y_t = 0.6*y_{t-1} - 0.2*y_{t-2} + e_t
gen e = rnormal()
gen y = 0
replace y = 0.6*l.y - 0.2*l2.y + e in 3/l

* Random Walk
gen rw = 0
replace rw = l.rw + rnormal() in 2/l

* Seasonal series
gen seas = 0.1*t + 5*sin(2*_pi*t/12) + rnormal(0,.5)

* Financial returns
gen ret = rnormal()
replace ret = sqrt(0.01 + 0.1*l.ret^2 + 0.85*0.5) * rnormal() in 3/l
label var y    "AR(2) series"
label var rw   "Random walk"
label var seas "Seasonal series"
label var ret  "Financial returns"

/* ── 2. Visualizations ────────────────────────────────────────────────────── */
tsline y,   title("AR(2) Process") ytitle("Value")
graph export "05_time_series_analysis/ar2_plot.png", replace

tsline rw,  title("Random Walk") ytitle("Value")
graph export "05_time_series_analysis/rw_plot.png", replace

tsline seas, title("Seasonal Series") ytitle("Value")
graph export "05_time_series_analysis/seasonal_plot.png", replace

/* ── 3. ACF / PACF ────────────────────────────────────────────────────────── */
di _n "=== ACF and PACF ==="
ac y,  lags(24) title("ACF: AR(2)")
graph export "05_time_series_analysis/acf_ar2.png", replace

pac y, lags(24) title("PACF: AR(2)")
graph export "05_time_series_analysis/pacf_ar2.png", replace

/* ── 4. Stationarity Tests ────────────────────────────────────────────────── */
di _n "=== STATIONARITY TESTS ==="

* ADF test
di _n "--- ADF test: AR(2) ---"
dfuller y, lags(2) regress
di _n "--- ADF test: Random Walk ---"
dfuller rw, lags(1) regress
di _n "--- ADF test: Differenced Random Walk ---"
gen d_rw = d.rw
dfuller d_rw, lags(1) regress

* KPSS test (requires dfgls from SSC or use pperron)
* ssc install kpss if needed
* kpss y

* Phillips-Perron test
di _n "--- Phillips-Perron: AR(2) ---"
pperron y, lags(2)
di _n "--- Phillips-Perron: Random Walk ---"
pperron rw, lags(1)

/* ── 5. ARIMA Modeling ────────────────────────────────────────────────────── */
di _n "=== ARIMA(2,0,0) ==="
arima y, ar(1/2)
estat ic

di _n "=== ARIMA(1,1,1) on Random Walk ==="
arima rw, arima(1,1,1)
estat ic

di _n "=== ARIMA with Forecast ==="
arima y, ar(1/2)
predict y_hat, xb
predict y_se,  stdp

di "12-step forecast:"
arima y, ar(1/2)
predict yf, dynamic(190)
list t y yf in 190/200

/* ── 6. SARIMA ────────────────────────────────────────────────────────────── */
di _n "=== SARIMA ==="
arima seas, arima(1,1,1) sarima(1,1,0,12)
estat ic

/* ── 7. VAR ────────────────────────────────────────────────────────────────── */
di _n "=== VAR (Vector Autoregression) ==="
gen y2 = 0.1*y + rnormal()   /* a second series correlated with y */

* Lag selection
varsoc y y2, maxlag(8)

* Fit VAR(2)
var y y2, lags(1/2)

* Granger causality
vargranger

* Impulse Response Functions
irf create var_irf, step(10) set(myirf) replace
irf graph oirf, impulse(y2) response(y) ///
    title("IRF: Response of y to Shock in y2")
graph export "05_time_series_analysis/var_irf.png", replace

/* ── 8. Johansen Cointegration ────────────────────────────────────────────── */
di _n "=== JOHANSEN COINTEGRATION TEST ==="
* Create two cointegrated series
gen z  = rw + rnormal()   /* z is cointegrated with rw */
vecrank rw z, lags(2) trend(constant)

/* ── 9. ARCH/GARCH ────────────────────────────────────────────────────────── */
di _n "=== ARCH/GARCH ==="

* Test for ARCH effects
regress ret l.ret
estat archlm, lags(5)

* ARCH(1) model
arch ret, arch(1)

* GARCH(1,1)
arch ret, arch(1) garch(1)
predict vol, variance
tsline vol, title("Conditional Variance: GARCH(1,1)") ytitle("Variance")
graph export "05_time_series_analysis/garch_volatility.png", replace

di _n "=== TIME SERIES ANALYSIS COMPLETE ==="

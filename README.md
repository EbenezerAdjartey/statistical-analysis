# Statistical Analysis Portfolio

A comprehensive collection of statistical analysis implementations across three platforms:
**Python (Jupyter Notebooks)**, **R scripts**, and **Stata do-files**.

Each topic is self-contained, fully documented, and uses synthetic datasets so all code runs
without external data dependencies.

## Repository Structure

| Folder | Topic | Methods Covered |
|--------|-------|-----------------|
| `01_descriptive_statistics/` | Descriptive Statistics | Central tendency, dispersion, skewness, kurtosis, visualizations |
| `02_probability_distributions/` | Probability Distributions | Discrete & continuous distributions, PDF/CDF, goodness-of-fit |
| `03_hypothesis_testing/` | Hypothesis Testing | t-tests, ANOVA, chi-square, z-tests, multiple comparisons |
| `04_regression_analysis/` | Regression Analysis | OLS, logistic, probit, Tobit, Poisson, IV/2SLS, quantile |
| `05_time_series_analysis/` | Time Series Analysis | ARIMA, SARIMA, VAR, VECM, ARCH/GARCH, stationarity tests |
| `06_panel_data_analysis/` | Panel Data Analysis | Pooled OLS, fixed effects, random effects, Hausman, GMM |
| `07_survival_analysis/` | Survival Analysis | Kaplan-Meier, log-rank, Cox PH, parametric models, competing risks |
| `08_nonparametric_methods/` | Nonparametric Methods | Mann-Whitney, Kruskal-Wallis, KDE, bootstrap, rank correlations |
| `09_multivariate_analysis/` | Multivariate Analysis | PCA, EFA, clustering, LDA, MANOVA, canonical correlation |
| `10_bayesian_statistics/` | Bayesian Statistics | Bayes' theorem, MCMC, Bayesian regression, credible intervals |
| `11_machine_learning/` | Supervised Learning | Regression, classification, trees, ensembles, SVM, neural nets |
| `12_model_evaluation/` | Model Evaluation | Cross-validation, ROC/AUC, confusion matrix, SHAP, tuning |

## File Formats

Each topic folder contains three equivalent implementations:

- **`.ipynb`** — Jupyter Notebook (Python): uses `numpy`, `pandas`, `scipy`, `statsmodels`, `sklearn`, `matplotlib`, `seaborn`
- **`.R`** — R Script: uses base R plus `tidyverse`, `lmtest`, `sandwich`, `survival`, `caret`, and topic-specific packages
- **`.do`** — Stata Do-File: compatible with **Stata 16+** using standard built-in commands

## Getting Started

### Python (Jupyter)
```bash
pip install numpy pandas scipy statsmodels scikit-learn matplotlib seaborn lifelines pymc arviz xgboost lightgbm shap
jupyter notebook
```

### R
```r
install.packages(c("tidyverse", "lmtest", "sandwich", "car", "MASS", "survival",
                   "survminer", "plm", "tseries", "forecast", "vars", "rugarch",
                   "FactoMineR", "factoextra", "cluster", "caret", "randomForest",
                   "e1071", "pROC", "xgboost", "ggplot2", "reshape2"))
```

### Stata
Open any `.do` file in Stata 16+ and run with `Do`. External packages (`outreg2`, `estout`,
`xtabond2`, `stcomprisk`) can be installed with `ssc install <package>`.

## Topics Overview

### 01 — Descriptive Statistics
Comprehensive summary statistics including mean, median, mode, variance, standard deviation,
IQR, range, skewness, and kurtosis. Frequency tables, cross-tabulations, and visualizations
(histograms, boxplots, density plots, Q-Q plots).

### 02 — Probability Distributions
All major discrete distributions (Binomial, Poisson, Geometric, Hypergeometric, Negative Binomial)
and continuous distributions (Normal, t, Chi-square, F, Exponential, Gamma, Beta, Uniform, Weibull).
Covers PDF/PMF, CDF, quantile functions, distribution fitting, and goodness-of-fit tests.

### 03 — Hypothesis Testing
One-sample and two-sample t-tests, paired t-test, one-way and two-way ANOVA, chi-square tests
(goodness-of-fit and independence), z-test for proportions, F-test for variance equality, and
multiple comparison corrections (Bonferroni, Tukey, Scheffé).

### 04 — Regression Analysis
Simple and multiple OLS regression with diagnostics (VIF, heteroskedasticity tests). Logistic
regression (binary, multinomial, ordinal), probit and Tobit models, Poisson and Negative Binomial
regression, IV/2SLS regression, and quantile regression.

### 05 — Time Series Analysis
Time series visualization, ACF/PACF, stationarity tests (ADF, KPSS, PP), ARIMA/SARIMA modeling
and forecasting, VAR (Vector Autoregression), VECM (cointegration/error correction), and
ARCH/GARCH volatility modeling.

### 06 — Panel Data Analysis
Pooled OLS, fixed effects (within estimator), random effects (GLS), Hausman specification test,
first-difference estimator, dynamic panel GMM (Arellano-Bond), and cluster-robust standard errors.

### 07 — Survival Analysis
Kaplan-Meier survival curves, log-rank test, Cox Proportional Hazards model with diagnostics,
parametric survival models (Weibull, exponential, log-normal), competing risks analysis, and
time-varying covariates.

### 08 — Nonparametric Methods
Mann-Whitney U / Wilcoxon rank-sum test, Wilcoxon signed-rank test, Kruskal-Wallis test,
Spearman and Kendall rank correlations, Kolmogorov-Smirnov test, kernel density estimation,
and bootstrap resampling with confidence intervals.

### 09 — Multivariate Analysis
Principal Component Analysis (PCA), Exploratory Factor Analysis (EFA), cluster analysis
(k-means and hierarchical), Linear Discriminant Analysis (LDA), MANOVA, and Canonical
Correlation Analysis.

### 10 — Bayesian Statistics
Bayes' theorem illustration, Bayesian inference with conjugate priors, Bayesian linear regression,
MCMC sampling (Metropolis-Hastings), credible intervals vs confidence intervals, and Bayesian
hypothesis testing (Bayes Factor).

### 11 — Machine Learning (Supervised)
Linear, Ridge, and Lasso regression; logistic regression for classification; decision trees;
random forests and gradient boosting (XGBoost/LightGBM); SVM; KNN; and basic MLP neural networks.

### 12 — Model Evaluation & Testing
Train/test split, k-fold cross-validation, confusion matrix metrics (accuracy, precision, recall,
F1), ROC/AUC, regression metrics (RMSE, MAE, R²), hyperparameter tuning (GridSearch/RandomSearch),
SHAP values, and learning curves.

## Author

**Ebenezer Adjartey**
GitHub: [@EbenezerAdjartey](https://github.com/EbenezerAdjartey)

---
*All analyses use synthetic or built-in datasets. Code is written for educational and portfolio purposes.*

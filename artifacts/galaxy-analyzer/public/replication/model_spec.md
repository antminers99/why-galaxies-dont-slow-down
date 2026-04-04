# Model Specifications — Replication Package
# Galaxy Rotation Curve a0 Variation Analysis

## Fitting Method (all models)

- **Algorithm**: Ordinary Least Squares (OLS)
- **Target**: logA0 (log10 of per-galaxy MOND acceleration scale)
- **No regularization**: All models are standard multivariate linear regression
- **Software**: Any OLS implementation will reproduce results (R lm(), Python
  statsmodels, numpy.linalg.lstsq, etc.)

## Model Hierarchy

### M0 — Null (intercept only)
```
logA0 = beta_0
```
- **k** = 1 parameter
- **Purpose**: Baseline — no galaxy properties predict a0 variation

### A' — Gas + Environment (2 predictors)
```
logA0 = beta_0 + beta_1 * logMHI + beta_2 * logMhost
```
- **k** = 3 parameters
- **Expected signs**: beta_1 < 0, beta_2 < 0

### B' — 4-variable model
```
logA0 = beta_0 + beta_1 * logMHI + beta_2 * rcWiggliness
        + beta_3 * logMhost + beta_4 * logSigma0
```
- **k** = 5 parameters
- **Expected signs**: beta_1 < 0, beta_2 > 0, beta_3 < 0, beta_4 > 0

### B″ — Final 6-variable model (FROZEN BASELINE)
```
logA0 = beta_0 + beta_1 * logMHI + beta_2 * rcWiggliness
        + beta_3 * logMhost + beta_4 * logSigma0
        + beta_5 * logMeanRun + beta_6 * Upsilon_perp
```
- **k** = 7 parameters (intercept + 6 predictors)
- **Expected coefficients**:
  - beta_0 (intercept) = 4.60 +/- ...
  - beta_1 (logMHI) = -0.243 (negative)
  - beta_2 (rcWiggliness) = +1.279 (positive)
  - beta_3 (logMhost) = -0.150 (negative)
  - beta_4 (logSigma0) = +0.168 (positive)
  - beta_5 (logMeanRun) = +0.446 (positive)
  - beta_6 (Upsilon_perp) = +0.659 (positive)

### IMPORTANT: Upsilon_perp construction
Before fitting B″, you must first construct Upsilon_perp:
1. Compute logUpsilon = log10(Upsilon_disk) for all 45 galaxies
2. Fit auxiliary OLS: logUpsilon ~ logMHI + logSigma0 + morphT
3. Take residuals = Upsilon_perp
4. Use these residuals as the 6th predictor in B″

This is NOT optional — using raw logUpsilon instead of Upsilon_perp
will give different (worse) results because raw Upsilon is confounded
with gas mass and surface density.

## Evaluation Metrics

### LOO gap% (primary metric)
- **Definition**: Leave-one-out cross-validated variance explained,
  expressed as a percentage
- **Computation**:
  1. For each galaxy i = 1..N:
     a. Remove galaxy i from the sample
     b. Fit OLS on remaining N-1 galaxies
     c. Predict logA0 for galaxy i using the held-out coefficients
     d. Record squared prediction error: e_i^2 = (logA0_i - predicted_i)^2
  2. LOO RMS = sqrt(mean(e_i^2))
  3. SD_Y = standard deviation of all N logA0 values
  4. gap% = 100 * (1 - LOO_RMS^2 / SD_Y^2)
- **Interpretation**: Percentage of logA0 variance explained on truly
  unseen data. Higher = better. Can be negative if model is worse
  than predicting the mean.
- **B″ expected value**: 49.7%

### Residual scatter (tau)
- **Definition**: Standard deviation of OLS residuals (in-sample)
- **Computation**: tau = SD(logA0_obs - logA0_predicted) using all N
- **Units**: dex
- **B″ expected value**: 0.183 dex (full-sample), 0.156 dex (LOO-adjusted
  with df correction)

### tau progression across models
- M0: tau ~ 0.258 dex
- A': tau ~ 0.215 dex
- B': tau ~ 0.190 dex
- B″: tau ~ 0.183 dex

### Intrinsic scatter (tau_int)
- **Definition**: Residual scatter after subtracting estimated
  measurement noise in quadrature
- **Computation**: tau_int = sqrt(max(0, tau^2 - mean(sigma_meas_i^2)))
- **B″ expected value**: ~0.105 dex (95% CI: [0.041, 0.143])
- **Note**: Does NOT touch zero — genuine intrinsic scatter remains

## Expected Outputs for Verification

A successful replication should produce:

| Quantity | Expected Value | Tolerance |
|----------|---------------|-----------|
| B″ intercept | ~4.60 | +/- 0.5 |
| beta(logMHI) | -0.243 | sign must be negative |
| beta(rcWig) | +1.279 | sign must be positive |
| beta(logMhost) | -0.150 | sign must be negative |
| beta(logSigma0) | +0.168 | sign must be positive |
| beta(logMeanRun) | +0.446 | sign must be positive |
| beta(Ups_perp) | +0.659 | sign must be positive |
| LOO gap% | ~49.7% | within [45%, 55%] |
| In-sample R^2 | ~0.73 | within [0.65, 0.80] |
| Residual SD | ~0.156 dex | within [0.13, 0.18] |

All 6 coefficient signs are the critical test. If all signs match
and LOO gap% is in range, the replication is confirmed.

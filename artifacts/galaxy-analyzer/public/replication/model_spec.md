# Model Specification — Final M5 and M3 Laws

## M5: Five-Axis Predictive Law (Final)

### Equation
```
log10(a0) = 4.978
           - 0.236 * log10(MHI / 1e9 Msun)
           - 0.172 * log10(Mhost / Msun)
           + 0.145 * log10(Sigma0 / [Msun/pc^2])
           + 0.452 * log10(MeanRun)
           + 0.658 * Upsilon_perp
```

### Coefficients
| Variable | Coefficient | Rational | Sign | Physical axis |
|----------|------------|----------|------|---------------|
| intercept | +4.978 | - | + | - |
| logMHI | -0.236 | -1/4 | - | Gas reservoir |
| logMhost | -0.172 | -1/6 | - | Environmental depth |
| logSigma0 | +0.145 | +1/7 | + | Baryon concentration |
| logMeanRun | +0.452 | +3/7 | + | Dynamical coherence |
| Ups_perp | +0.658 | +2/3 | + | Stellar structure |

### Performance
- LOO gap% = 51.0%
- RMS = 0.157 dex
- In-sample R^2 ~ 0.73
- N = 45
- k = 6 (intercept + 5 predictors)

## M3: Three-Axis Compressed State Law

### Equation
```
log10(a0) = 5.182
           - 0.198 * log10(MHI / 1e9 Msun)
           - 0.155 * log10(Mhost / Msun)
           + 0.459 * log10(MeanRun)
```

### Performance
- LOO gap% = 44.1%
- RMS = 0.193 dex
- Retains 86% of M5 signal
- N = 45
- k = 4 (intercept + 3 predictors)

## Upsilon_perp Construction

```
Step 1: Fit OLS on N=45
  logUpsilon_disk = a0 + a1*logMHI + a2*logSigma0 + a3*morphT + eps

Step 2: This regression explains R^2 ~ 0.68

Step 3: Upsilon_perp_i = eps_i (the OLS residual for galaxy i)
```

## Dropped Variable: rcWiggliness

rcWiggliness (rotation-curve irregularity) was part of an earlier 6-variable
model (B''). It was dropped from the final M5 model because:
- Incremental F-test = 1.50 (not significant at any standard level)
- 14.5% bootstrap sign-flip rate (worst of all candidates)
- Removing it improves LOO gap% from ~49.7% to 51.0%

## Fitting Method

Ordinary Least Squares (OLS) on all N=45 galaxies.
Verified identical results under Huber IRLS and bootstrap mean.

## Evaluation Metrics

| Metric | Formula |
|--------|---------|
| LOO gap% | 100 * (1 - LOO_RMS^2 / SD(y)^2) |
| LOO_RMS | sqrt(mean(e_LOO_i^2)) where e_LOO_i = y_i - y_hat_{-i} |
| SD(y) | standard deviation of logA0, ddof=1 |

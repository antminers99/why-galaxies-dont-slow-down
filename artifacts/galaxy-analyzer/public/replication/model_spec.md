# Model Specification — Final M5 and M3 Laws

**STATUS: SUPERSEDED.** These models are from Phases 1-100. The underlying
Vflat/InnerVmax relation was shown to be a geometric artifact (Phase 101).
The current active result (Phases 101-122) uses logOMD as target and
log(MHI/L3.6) as predictor. See REPRODUCIBLE_RESULT.md.

---

## M5: Five-Axis Predictive Law (Corrected)

### Equation
```
log10(a0) = 5.016
           - 0.235 * log10(MHI / 1e9 Msun)
           - 0.175 * log10(Mhost / Msun)
           + 0.146 * log10(Sigma0 / [Msun/pc^2])
           + 0.447 * log10(MeanRun)
           + 0.372 * Upsilon_perp
```

### Coefficients
| Variable | Corrected | Old (buggy) | Sign | Physical axis |
|----------|-----------|-------------|------|---------------|
| intercept | +5.016 | +4.978 | + | - |
| logMHI | -0.235 | -0.236 | - | Gas reservoir |
| logMhost | -0.175 | -0.172 | - | Environmental depth |
| logSigma0 | +0.146 | +0.145 | + | Baryon concentration |
| logMeanRun | +0.447 | +0.452 | + | Dynamical coherence |
| Ups_perp | +0.372 | +0.658 | + | Stellar structure |

### Performance
- LOO gap% = 46.6% (was 50.9% with bug)
- LOO RMS = 0.189 dex (was 0.181 dex)
- In-sample R^2 = 0.588
- N = 45
- k = 6 (intercept + 5 predictors)

### Errata: T||5 Bug
The old coefficients and performance numbers were inflated by a JavaScript
bug: `sparcMap[g.name]?.T || 5` treated morphological type T=0 as falsy,
replacing it with T=5 for NGC4138 and UGC06786. This distorted Upsilon_perp
construction and inflated its coefficient from 0.372 to 0.658.
Fixed by using `?? 5` (nullish coalescing).

## M3: Three-Axis Compressed State Law (Unaffected)

### Equation
```
log10(a0) = 5.182
           - 0.198 * log10(MHI / 1e9 Msun)
           - 0.155 * log10(Mhost / Msun)
           + 0.459 * log10(MeanRun)
```

### Performance
- LOO gap% = 44.1%
- LOO RMS = 0.193 dex
- Retains 95% of corrected M5 signal (was 86% relative to buggy M5)
- N = 45
- k = 4 (intercept + 3 predictors)

## Upsilon_perp Construction

```
Step 1: Fit OLS on N=45
  logUpsilon_disk = a0 + a1*logMHI + a2*logSigma0 + a3*morphT + eps
  CRITICAL: morphT=0 is valid (S0 galaxies). Use ?? 5, NOT || 5.

Step 2: This regression explains R^2 ~ 0.68

Step 3: Upsilon_perp_i = eps_i (the OLS residual for galaxy i)
```

## Dropped Variable: rcWiggliness

rcWiggliness (rotation-curve irregularity) was part of an earlier 6-variable
model (B''). It was dropped from the final M5 model because:
- Incremental F-test = 1.50 (not significant at any standard level)
- 14.5% bootstrap sign-flip rate (worst of all candidates)
- Removing it improves LOO

## Fitting Method

Ordinary Least Squares (OLS) on all N=45 galaxies.

## Evaluation Metrics

| Metric | Formula |
|--------|---------|
| LOO gap% | 100 * (1 - LOO_RMS^2 / SD(y)^2) |
| LOO_RMS | sqrt(mean(e_LOO_i^2)) where e_LOO_i = y_i - y_hat_{-i} |
| SD(y) | standard deviation of logA0, ddof=1 |

## Units Note

logA0 values in the CSV are in **log10((km/s)^2/kpc)**, NOT log10(m/s^2).
To convert: log10(m/s^2) = logA0 - 13.489

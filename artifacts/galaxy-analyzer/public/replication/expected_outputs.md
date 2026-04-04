# Expected Outputs — Replication Verification
# Run replicate_from_scratch.py and compare against these values

## B″ Coefficients (from CSV dataset)

| Variable | Expected Coeff | Expected Sign | Tolerance |
|----------|---------------|---------------|-----------|
| intercept | +4.55 | (any) | +/- 1.0 |
| logMHI | -0.243 | NEGATIVE | must be negative |
| rcWiggliness | +1.50 | POSITIVE | must be positive |
| logMhost | -0.148 | NEGATIVE | must be negative |
| logSigma0 | +0.172 | POSITIVE | must be positive |
| logMeanRun | +0.435 | POSITIVE | must be positive |
| Upsilon_perp | +0.458 | POSITIVE | must be positive |

## Key Metrics

| Metric | Expected Value | Acceptable Range |
|--------|---------------|-----------------|
| LOO gap% | 45.7% | [40%, 55%] |
| In-sample R^2 | ~0.73 | [0.65, 0.80] |
| Residual SD | 0.162 dex | [0.13, 0.19] |
| N galaxies | 45 | exactly 45 |

## Upsilon_perp Construction Check

| Quantity | Expected |
|----------|---------|
| Auxiliary R^2 (logUps ~ logMHI + logSigma0 + morphT) | ~0.68 |
| Upsilon_perp mean | ~0.000 (by construction) |

## Tau Progression

| Model | Residual SD (dex) |
|-------|------------------|
| M0 (null) | ~0.258 |
| A' (2 vars) | ~0.215 |
| B' (4 vars) | ~0.190 |
| B″ (6 vars) | ~0.162 |

## Critical Verification Checks

A replication is CONFIRMED if ALL of:
1. All 6 predictor coefficient signs match expected signs
2. LOO gap% falls within [40%, 55%]
3. R^2 falls within [0.65, 0.80]
4. N = 45 galaxies (no missing data)

## Residual Diagnostics (informational)

| Statistic | Expected | Notes |
|-----------|---------|-------|
| Skewness | near 0 | [-1, +1] acceptable |
| Excess kurtosis | near 0 | [-2, +3] acceptable |
| Runs test z | near 0 | [-2, +2] acceptable |

## Notes on Numerical Precision

Small differences (< 10%) in coefficient magnitudes are expected
due to:
- Floating-point rounding in CSV (6 decimal places)
- Different OLS implementations (QR vs normal equations vs SVD)
- Upsilon_perp re-computation from rounded logUpsilon_disk values

The SIGNS are the robust, implementation-independent result.
The magnitudes carry the model's physical content.
The LOO gap% measures out-of-sample predictive power.

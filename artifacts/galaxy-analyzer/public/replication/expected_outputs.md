# Expected Outputs — Replication Verification
# Run replicate_from_scratch.py and compare against these values

## M5 Coefficients (from CSV dataset)

| Variable | Expected | Tolerance |
|----------|----------|-----------|
| intercept | +4.978 | +/- 0.10 |
| logMHI | -0.236 | +/- 0.02 |
| logMhost | -0.172 | +/- 0.02 |
| logSigma0 | +0.145 | +/- 0.02 |
| logMeanRun | +0.452 | +/- 0.03 |
| Ups_perp | +0.658 | +/- 0.05 |

## M3 Coefficients

| Variable | Expected | Tolerance |
|----------|----------|-----------|
| intercept | +5.182 | +/- 0.10 |
| logMHI | -0.198 | +/- 0.02 |
| logMhost | -0.155 | +/- 0.02 |
| logMeanRun | +0.459 | +/- 0.03 |

## M5 Performance Metrics

| Metric | Expected | Pass range |
|--------|----------|------------|
| LOO gap% | 51.0% | [45%, 55%] |
| RMS | 0.157 dex | [0.14, 0.17] |
| In-sample R^2 | ~0.73 | [0.65, 0.80] |
| All 5 signs correct | 5/5 | 5/5 |

## M3 Performance Metrics

| Metric | Expected | Pass range |
|--------|----------|------------|
| LOO gap% | 44.1% | [40%, 50%] |
| RMS | 0.193 dex | [0.17, 0.22] |
| M3/M5 retention | 86% | [80%, 95%] |

## Upsilon_perp Auxiliary Regression

| Metric | Expected | Pass range |
|--------|----------|------------|
| R^2 (confounders) | 0.68 | [0.60, 0.75] |
| Ups_perp mean | ~0.000 | [-0.01, 0.01] |

## Residual Diagnostics (M5)

| Metric | Expected | Conclusion |
|--------|----------|------------|
| Skewness | -0.124 | Near-symmetric |
| Jarque-Bera | 1.22 | Cannot reject normality (critical 5.99) |
| Breusch-Pagan nR^2 | 4.63 | No heteroscedasticity (critical 11.1) |
| Stratification biases | 0/7 significant | No subpopulation bias |

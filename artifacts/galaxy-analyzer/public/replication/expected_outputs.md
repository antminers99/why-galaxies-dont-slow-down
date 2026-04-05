# Expected Outputs — Replication Verification
# Run replicate_from_scratch.py and compare against these values

**STATUS: SUPERSEDED.** These are for the Phases 1-100 M5/M3 models.
The current active result is in REPRODUCIBLE_RESULT.md (Phases 101-122).

**ERRATA (T||5 bug fix):** All M5 numbers updated to reflect corrected
morphological type handling. M3 numbers are unchanged.

---

## M5 Coefficients (CORRECTED — from CSV dataset)

| Variable | Expected | Tolerance |
|----------|----------|-----------|
| intercept | +5.016 | +/- 0.10 |
| logMHI | -0.235 | +/- 0.02 |
| logMhost | -0.175 | +/- 0.02 |
| logSigma0 | +0.146 | +/- 0.02 |
| logMeanRun | +0.447 | +/- 0.03 |
| Ups_perp | +0.372 | +/- 0.05 |

## M3 Coefficients (unchanged)

| Variable | Expected | Tolerance |
|----------|----------|-----------|
| intercept | +5.182 | +/- 0.10 |
| logMHI | -0.198 | +/- 0.02 |
| logMhost | -0.155 | +/- 0.02 |
| logMeanRun | +0.459 | +/- 0.03 |

## M5 Performance Metrics (CORRECTED)

| Metric | Expected | Pass range |
|--------|----------|------------|
| LOO gap% | 46.6% | [42%, 52%] |
| LOO RMS | 0.189 dex | [0.17, 0.21] |
| In-sample R^2 | 0.588 | [0.55, 0.65] |
| All 5 signs correct | 5/5 | 5/5 |

## M3 Performance Metrics (unchanged)

| Metric | Expected | Pass range |
|--------|----------|------------|
| LOO gap% | 44.1% | [40%, 50%] |
| LOO RMS | 0.193 dex | [0.17, 0.22] |
| M3/M5 retention | 95% | [85%, 100%] |

## Upsilon_perp Auxiliary Regression

| Metric | Expected | Pass range |
|--------|----------|------------|
| R^2 (confounders) | 0.68 | [0.60, 0.75] |
| Ups_perp mean | ~0.000 | [-0.01, 0.01] |

## Residual Diagnostics (M5)

| Metric | Expected | Conclusion |
|--------|----------|------------|
| Skewness | ~-0.1 | Near-symmetric |
| Jarque-Bera | < 5.99 | Cannot reject normality |
| Stratification biases | 0/7 significant | No subpopulation bias |

## Units Note

logA0 in the CSV is in **log10((km/s)^2/kpc)**, NOT log10(m/s^2).
Range: [3.1, 4.1]. To convert to m/s^2: subtract 13.489.

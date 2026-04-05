# Reproduction Guide

**STATUS: SUPERSEDED.** This guide covers the Phases 1-100 M5/M3 models.
The current active analysis (Phases 101-122) uses a different framework.
See REPRODUCIBLE_RESULT.md for the current result and reproduction steps.

**ERRATA (T||5 bug):** All M5 numbers have been corrected. See below.

---

## How to Reproduce All Results from Scratch

### Requirements
- Python 3.7 or later
- NumPy (any recent version)
- No other dependencies required

### Files in This Package

| File | Description |
|------|-------------|
| `N45_final_dataset.csv` | Complete N=45 sample with all variables |
| `replicate_from_scratch.py` | Self-contained replication script |
| `variable_definitions.md` | Precise definitions of every variable |
| `model_spec.md` | Model equations and expected outputs |
| `expected_outputs.md` | Numerical targets for verification |

### Quick Start

```bash
cd replication/
python replicate_from_scratch.py
```

This will:
1. Load the N=45 dataset
2. Construct Upsilon_perp via orthogonalization
3. Fit M0 (universal constant), M3 (3-axis), and M5 (5-axis)
4. Compute leave-one-out cross-validation for each
5. Run residual diagnostics
6. Print verification checks (PASS/FAIL)

### Expected Output Summary (CORRECTED)

| Check | Expected |
|-------|----------|
| M5 LOO gap% | 46.6% (range: 42-52%) |
| M5 LOO RMS | 0.189 dex |
| M5 In-sample R^2 | 0.588 |
| M3 LOO gap% | 44.1% (range: 40-50%) |
| All M5 coefficient signs correct | 5/5 |
| M3 retains ~95% of corrected M5 signal | Yes |

### The Two Final Laws (CORRECTED)

**M5 (Predictive Law)**:
```
log(a0) = 5.016 - 0.235*logMHI - 0.175*logMhost
        + 0.146*logSigma0 + 0.447*logMeanRun + 0.372*Upsilon_perp
```

**M3 (Compressed State Law)**:
```
log(a0) = 5.182 - 0.198*logMHI - 0.155*logMhost + 0.459*logMeanRun
```

### Dataset Columns

| Column | Description | Units |
|--------|-------------|-------|
| galaxy_name | SPARC galaxy identifier | - |
| logA0 | log10(a0) target | **log10((km/s)^2/kpc)** |
| logMHI | log10(M_HI / 10^9 Msun) | dex |
| rcWiggliness | RC irregularity (dropped from M5, kept for reference) | dimensionless |
| logMhost | log10(M_host / Msun) | dex |
| logSigma0 | log10(central surface mass density) | log10(Msun/pc^2) |
| logMeanRun | log10(mean run length in MOND residuals) | dex |
| logUpsilon_disk | log10(stellar M/L at 3.6um) | dex |
| morphT | de Vaucouleurs morphological type (T=0 is valid!) | - |
| Vflat_km_s | Flat rotation velocity (metadata only) | km/s |
| D_Mpc | Distance | Mpc |
| inc_deg | Inclination | degrees |
| RHI_kpc | HI radius | kpc |
| Rdisk_kpc | Disk scale length | kpc |
| L36_1e9Lsun | 3.6um luminosity | 10^9 Lsun |
| Q_SPARC | SPARC quality flag | 1 or 2 |

**Units warning**: logA0 is in log10((km/s)^2/kpc), NOT log10(m/s^2).
To convert: log10(m/s^2) = logA0 - 13.489. Values in CSV range [3.1, 4.1].

### Upsilon_perp Construction

Upsilon_perp is NOT in the CSV directly. It is constructed by the replication script:

1. Regress logUpsilon_disk on [logMHI, logSigma0, morphT]
2. Take the OLS residual = Upsilon_perp
3. This removes 68% of M/L variance, isolating the independent component

**CRITICAL**: morphT=0 is a valid morphological type (S0 galaxies). When
providing a default for missing T values, use nullish coalescing (?? 5 in JS,
or explicit None check in Python), NOT logical OR (|| 5), which treats T=0
as falsy and replaces it with T=5.

### Errata: T||5 Bug

A bug in the original JavaScript scripts used `sparcMap[g.name]?.T || 5`,
which incorrectly replaced T=0 with T=5 for NGC4138 and UGC06786.
This distorted the Upsilon_perp construction and inflated the M5 Upsilon_perp
coefficient from 0.372 to 0.658, and LOO gap% from 46.6% to 50.9%.

The Python replication script (replicate_from_scratch.py) reads T directly
from the CSV and handles T=0 correctly, which is why it produces different
numbers from the original JS scripts. The corrected JS scripts now match.

### Verification

If all checks print "PASS", you have successfully replicated the analysis.
Note that the replication script's expected values should match the corrected
numbers in this document, not the original (buggy) values.

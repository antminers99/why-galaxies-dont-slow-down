# Reproduction Guide

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
| `paper.md` | Full paper text |
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

### Expected Output Summary

| Check | Expected |
|-------|----------|
| M5 LOO gap% | 51.0% (range: 45-55%) |
| M5 RMS | 0.157 dex |
| M3 LOO gap% | 44.1% (range: 40-50%) |
| All M5 coefficient signs correct | 5/5 |
| M3 retains ~86% of M5 signal | Yes |

### The Two Final Laws

**M5 (Predictive Law)**:
```
log(a0) = 4.978 - 0.236*logMHI - 0.172*logMhost
        + 0.145*logSigma0 + 0.452*logMeanRun + 0.658*Upsilon_perp
```

**M3 (Compressed State Law)**:
```
log(a0) = 5.182 - 0.198*logMHI - 0.155*logMhost + 0.459*logMeanRun
```

### Dataset Columns

| Column | Description | Units |
|--------|-------------|-------|
| galaxy_name | SPARC galaxy identifier | - |
| logA0 | log10(a0) target | log10(m/s^2) |
| logMHI | log10(M_HI / 10^9 Msun) | dex |
| rcWiggliness | RC irregularity (dropped from M5, kept for reference) | dimensionless |
| logMhost | log10(M_host / Msun) | dex |
| logSigma0 | log10(central surface mass density) | log10(Msun/pc^2) |
| logMeanRun | log10(mean run length in MOND residuals) | dex |
| logUpsilon_disk | log10(stellar M/L at 3.6um) | dex |
| morphT | de Vaucouleurs morphological type | - |
| Vflat_km_s | Flat rotation velocity (metadata only) | km/s |
| D_Mpc | Distance | Mpc |
| inc_deg | Inclination | degrees |
| RHI_kpc | HI radius | kpc |
| Rdisk_kpc | Disk scale length | kpc |
| L36_1e9Lsun | 3.6um luminosity | 10^9 Lsun |
| Q_SPARC | SPARC quality flag | 1 or 2 |

### Upsilon_perp Construction

Upsilon_perp is NOT in the CSV directly. It is constructed by the replication script:

1. Regress logUpsilon_disk on [logMHI, logSigma0, morphT]
2. Take the OLS residual = Upsilon_perp
3. This removes 68% of M/L variance, isolating the independent component

### Verification

If all checks print "PASS", you have successfully replicated the analysis.
If any check prints "FAIL", please verify:
- CSV file is unmodified (46 lines: 1 header + 45 data rows)
- Python version >= 3.7
- NumPy is installed correctly

### Contact

For questions about this replication package, see the accompanying paper
and supplementary materials.

---
name: galaxy-data
description: Reference data, constants, and computed results for the Galaxy Rotation Curve Analyzer project. Use whenever working on galaxy analysis, modifying the report, running scripts, or needing physical constants and dataset summaries.
---

# Galaxy Data Reference

## Physical Constants

| Constant | Value | Units |
|----------|-------|-------|
| a₀ (canonical) | 1.2×10⁻¹⁰ | m/s² |
| a₀ in galaxy units | 3702 | (km/s)²/kpc |
| a₀ (global fit, combined) | 1730.4 | (km/s)²/kpc |
| c | 299,792,458 | m/s |
| H₀ | 67.4 | km/s/Mpc |
| cH₀ | 6.548×10⁻¹⁰ | m/s² |
| cH₀/2π | 1.042×10⁻¹⁰ | m/s² |
| a₀/(cH₀/2π) | 0.87 (per-galaxy median) or 0.56–0.62 (global fit) | — |
| κ_dS = c²√(Λ/3) | 5.47×10⁻¹⁰ | m/s² |
| κ_dS/2π | 8.7×10⁻¹¹ | m/s² |
| 1 (km/s)²/kpc → m/s² | 3.241×10⁻¹⁴ | conversion |
| Υ_disk (fixed) | 0.5 | M☉/L☉ |
| Υ_bulge (fixed) | 0.7 | M☉/L☉ |

## Dataset Summary

| Property | SPARC | LITTLE THINGS | Combined |
|----------|-------|---------------|----------|
| Galaxies | 175 (171 fitted) | 22 | 197 (193 fitted) |
| Source | Lelli+2016 | Oh+2015 | — |
| Photometry | Spitzer 3.6μm | VLA HI + Spitzer | — |
| Total points | ~3500 | ~600 | 4,123 |
| Global a₀ | 3702 (median) | — | 1730.4 (global fit) |
| RMS scatter | 0.198 dex | — | 0.2364 dex |
| Well-constrained | — | — | 128 galaxies |

## Per-Galaxy Metadata Fields (rar-analysis-real.json)

Each of the 175 SPARC galaxies has:
- `name`: Galaxy identifier
- `Vmax`: Maximum rotation velocity (km/s)
- `Rmax`: Radius at Vmax (kpc)
- `n`: Number of data points
- `distance`: Distance (Mpc)
- `inc`: Inclination (degrees)
- `L36`: 3.6μm luminosity (10⁹ L☉)
- `Rdisk`: Disk scale length (kpc)
- `MHI`: HI gas mass (10⁹ M☉)
- `sigma_bar`: Baryonic surface density
- `k_ratio`, `r_fid`, `eta_rot`, `eta_bar`, `S_out`, `Q_kin`
- `meanDeltaInner`, `meanDeltaOuter`: Inner/outer RAR residuals

## Diagnostic Test Results (April 2026)

### Test 1: Interpolation Function Audit

| Function | Best-fit a₀ | RMS (dex) | AIC | BIC |
|----------|-------------|-----------|-----|-----|
| McGaugh RAR: 1/(1-e^(-√y)) | 1809 | 0.2385 | -5910 | -5904 |
| Simple MOND: ½+√(¼+1/y) | 1988 | 0.2385 | -5910 | -5904 |
| Standard MOND: (1+√(1+4/y))/2 | 1988 | 0.2385 | -5910 | -5904 |
| Our ν(x)=x/(1+x) | 1972 | 0.2443 | -5811 | -5805 |

a₀ variation: 4.1% — **ROBUST** to interpolation choice.
Our framework fits worst (ΔAIC=99 vs McGaugh RAR).

### Test 2: Baryonic Calibration
- Best f_bar: 1.64 (but only 0.9% RMS improvement)
- f_bar and a₀ are degenerate
- Calibration is NOT the dominant scatter source

### Test 3: Residual Correlations (156 galaxies)

| Parameter | Pearson r | Significant? |
|-----------|-----------|-------------|
| Inclination | -0.159 | No |
| Distance | -0.139 | No |
| V_max | -0.313 | YES |
| Gas fraction (MHI/L36) | +0.641 | YES (STRONG) |
| Surface brightness | -0.644 | YES (STRONG) |
| Disk scale | -0.275 | No |
| N points | -0.157 | No |
| Luminosity | -0.524 | YES |

### Sub-Population a₀ Stability

| Sub-population | n_gal | a₀ (fit) | RMS | Δa₀ |
|---------------|-------|----------|-----|------|
| High SB | 144 | 3118 | 0.188 | -16% |
| Low SB | 31 | 1322 | 0.264 | -64% |
| High mass (V>150) | 56 | 3720 | 0.128 | +0.5% |
| Low mass (V<100) | 88 | 46170 | 0.265 | +1147% |
| Gas-rich | 73 | 1137 | 0.241 | -69% |
| Gas-poor | 102 | 3106 | 0.178 | -16% |

## Stress Test Results

| Test | Result |
|------|--------|
| Null (random g_bar) | 6.2σ rejection |
| Scramble (shuffled labels) | 6.1σ rejection |
| M/L perturbation (±0.15 dex) | 21/21 galaxies retain signal |
| Distance MC (±30%, 1000 iter) | 100% negative slopes |
| CV = 0.3% under distance perturbation | a₀ stable |

## Key Insight

The transition scale EXISTS and is robust. But:
- Precise a₀ value is method-dependent (1730 global vs 3702 median)
- Residuals correlate with gas fraction and surface brightness
- Low-mass galaxies are poorly constrained
- The cosmological ratio a₀ ≈ cH₀/2π depends on which a₀ you use

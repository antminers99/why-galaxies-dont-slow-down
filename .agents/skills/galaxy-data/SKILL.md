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
| a₀ (GOLD+i45 global fit) | 5275 | (km/s)²/kpc |
| a₀ (GOLD global fit) | 4285 | (km/s)²/kpc |
| a₀ (ALL global fit) | 1929 | (km/s)²/kpc |
| a₀ uncertainty | ±0.185 | dex |
| a₀ 1σ range | [3443, 8082] | (km/s)²/kpc |
| c | 299,792,458 | m/s |
| H₀ | 67.4 | km/s/Mpc |
| cH₀ | 6.548×10⁻¹⁰ | m/s² |
| cH₀/2π | 1.042×10⁻¹⁰ | m/s² |
| a₀/(cH₀/2π) | 1.641 (GOLD+i45) or 1.333 (GOLD) or 0.600 (ALL) | — |
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

### Test 3: Residual Correlations (CORRECTED — Phase 2)

Phase 1 had a code bug (gbar multiplied by a₀). Corrected values:

| Parameter | Pearson r | Significant? |
|-----------|-----------|-------------|
| Inclination | -0.159 | No |
| Distance | -0.139 | No |
| V_max | +0.414 | YES |
| Gas fraction (MHI/L36) | -0.079 | No |
| Surface brightness | +0.084 | No |
| Luminosity | +0.158 | Weak |

Gas fraction and SB correlations from Phase 1 were ARTIFACTS.
Vmax is the only real correlate, and it vanishes under quality cuts.

### Phase 2 Results (Partial Correlations + Quality Cuts)

- Vmax is primary residual driver (r=+0.414), strengthens controlling for SB (r=+0.610)
- SB becomes significant only when controlling for Vmax (r=-0.498)
- Gas fraction is never significant in any partial correlation
- Full 4-variable model: R²=0.638, dominated by Vmax (coeff=+1.56) and L (coeff=-0.30)
- Gas+SB alone: R²=0.007 (explain NOTHING)
- Under STRICT quality cuts (V>=80, inc>=30, n>=5, range>=0.5): ALL correlations vanish
- 54 galaxies with >1 dex g_bar range: residuals essentially clean

### Sub-Population a₀ Stability

| Sub-population | n_gal | n_pts | a₀ (fit) | RMS |
|---------------|-------|-------|----------|-----|
| ALL COMBINED | 175 | 1698 | 2670 | 0.200 |
| SPARC only | 173 | 1689 | 2688 | 0.200 |
| High mass (V>150) | 56 | 820 | 3720 | 0.128 |
| Low mass (V<100) | 88 | 555 | 46170 | 0.265 |
| Gas-poor | 102 | 1189 | 3106 | 0.178 |
| Gas-rich | 73 | 509 | 1137 | 0.241 |
| High SB | 144 | 1501 | 3118 | 0.188 |
| Low SB | 31 | 197 | 1322 | 0.264 |
| High quality (n>=10) | 124 | 1536 | 3024 | 0.183 |
| Low quality (n<10) | 51 | 162 | 3029 | 0.291 |

Low-mass divergence EXPLAINED: 73% of V<80 galaxies have g_bar range <0.5 dex.
They cannot constrain a₀. The huge values are fitting artifacts.
High vs low quality gives SAME a₀ (3024 vs 3029) — reassuring.

## Stress Test Results

| Test | Result |
|------|--------|
| Null (random g_bar) | 6.2σ rejection |
| Scramble (shuffled labels) | 6.1σ rejection |
| M/L perturbation (±0.15 dex) | 21/21 galaxies retain signal |
| Distance MC (±30%, 1000 iter) | 100% negative slopes |
| CV = 0.3% under distance perturbation | a₀ stable |

### Gold-Standard Pipeline v2.0.0 Results (April 2026)

Script: `scripts/gold-standard-pipeline.cjs`
Results: `public/gold-standard-results.json`

**5-Layer Global Fits (mean of 3 interpolation functions):**

| Sample | n_gal | n_pts | a₀ (km/s)²/kpc | a₀ (m/s²) | a₀/(cH₀/2π) |
|--------|-------|-------|----------------|-----------|-------------|
| ALL | 175 | 2062 | 1929 | 6.25×10⁻¹¹ | 0.600 |
| CLEAN | 78 | 1141 | 4883 | 1.58×10⁻¹⁰ | 1.519 |
| GOLD | 52 | 872 | 4285 | 1.39×10⁻¹⁰ | 1.333 |
| GOLD+i45 | 47 | 767 | 5275 | 1.71×10⁻¹⁰ | 1.641 |
| HIGH-MASS | 56 | 820 | 4035 | 1.31×10⁻¹⁰ | 1.255 |

**Per-Galaxy Medians:**

| Sample | n | median a₀ | MAD (dex) |
|--------|---|----------|-----------|
| ALL | 73 | 2270 | 0.558 |
| CLEAN | 58 | 3277 | 0.399 |
| GOLD | 44 | 3479 | 0.347 |
| GOLD+i45 | 39 | 4074 | 0.304 |
| HIGH-MASS | 42 | 3277 | 0.383 |

**Blind Residual Audit — GOLD (52 galaxies):**

| Variable | r | Verdict |
|----------|---|---------|
| V_max | +0.143 | CLEAN |
| Surface brightness | -0.085 | CLEAN |
| Gas fraction | -0.056 | CLEAN |
| Luminosity | -0.072 | CLEAN |
| Inclination | -0.307 | MARGINAL ⚠ |
| Distance | -0.151 | CLEAN |

**Blind Residual Audit — GOLD+i45 (47 galaxies):**

| Variable | r | Verdict |
|----------|---|---------|
| V_max | +0.095 | CLEAN |
| Surface brightness | -0.091 | CLEAN |
| Gas fraction | +0.003 | CLEAN |
| Luminosity | -0.104 | CLEAN |
| Inclination | -0.285 | CLEAN ✓ |
| Distance | -0.191 | CLEAN |

ALL 6 CLEAN. Inclination resolved by tightening to inc≥45°.

**Uncertainty Budget (6 components, quadrature):**

| Source | ±σ (dex) |
|--------|----------|
| Interpolation function | ±0.028 |
| Sample selection | ±0.058 |
| Global vs per-galaxy | ±0.112 |
| Υ★ (±15%) | ±0.066 |
| Inclination cut | ±0.106 |
| Dynamic range cut | ±0.044 |
| **TOTAL** | **±0.185** |

**Low-mass stress test:** a₀ = 104,713 (at fit boundary). Cannot constrain a₀.

## Key Insight

The transition scale EXISTS and is robust. GOLD+i45 is the recommended sample:
- ALL 6 residual audit variables clean (inclination resolved with inc≥45°)
- a₀ = 5275 ± 0.185 dex [3443, 8082] (km/s)²/kpc
- Per-galaxy MAD tightest of all samples (0.304 dex)
- Gas fraction and SB do NOT correlate with residuals
- Low-mass a₀ divergence is a fitting artifact (limited dynamic range)
- The cosmological ratio a₀/(cH₀/2π) = 1.641 for GOLD+i45 global fit
- Literature value (1.20×10⁻¹⁰) falls within our 1σ band
- Dominant uncertainties: global-vs-median and inclination-cut sensitivity

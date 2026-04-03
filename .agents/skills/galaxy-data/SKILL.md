---
name: galaxy-data
description: Reference data, constants, and computed results for the Galaxy Rotation Curve Analyzer project. Use whenever working on galaxy analysis, modifying the report, running scripts, or needing physical constants and dataset summaries.
---

# Galaxy Data Reference

## Physical Constants

| Constant | Value | Units |
|----------|-------|-------|
| a₀ (canonical, McGaugh 2016) | 1.2×10⁻¹⁰ | m/s² |
| a₀ in galaxy units | 3703 | (km/s)²/kpc |
| a₀ (GOLD+i45 unweighted) | 4837 | (km/s)²/kpc |
| a₀ (GOLD+i45 weighted) | 3545 | (km/s)²/kpc |
| a₀ (GOLD+i45 hierarchical) | 3374 | (km/s)²/kpc |
| a₀ (GOLD+i45 hier. median) | 4074 | (km/s)²/kpc |
| a₀ uncertainty (v3.0) | ±0.295 | dex |
| a₀ 1σ range | [2454, 9537] | (km/s)²/kpc |
| Hier. tau (between-galaxy) | 0.245 | dex |
| Hier. I² | 89.4 | % |
| c | 299,792,458 | m/s |
| H₀ | 67.4 | km/s/Mpc |
| cH₀ | 6.548×10⁻¹⁰ | m/s² |
| cH₀/2π | 1.042×10⁻¹⁰ | m/s² |
| a₀/(cH₀/2π) | 1.504 (UW), 1.102 (W), 1.049 (Hier) | — |
| 1 (km/s)²/kpc → m/s² | 3.241×10⁻¹⁴ | conversion |
| Υ_disk (fixed) | 0.5 | M☉/L☉ |
| Υ_bulge (fixed) | 0.7 | M☉/L☉ |
| Estimator | McGaugh RAR ONLY (v3.0) | — |

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

### Pipeline v3.0.0 Results (April 2026)

Script: `scripts/final-measurement-v3.cjs`
Results: `public/gold-standard-results.json` (v3.0.0)

**McGaugh-only unweighted fits:**

| Sample | n_gal | n_pts | a₀ (km/s)²/kpc | RMS | a₀/(cH₀/2π) |
|--------|-------|-------|----------------|-----|-------------|
| ALL | 175 | 1698 | 2670 | 0.200 | 0.830 |
| CLEAN | 78 | 1141 | 4463 | 0.144 | 1.388 |
| GOLD | 52 | 872 | 3934 | 0.143 | 1.223 |
| GOLD+i45 | 47 | 767 | 4837 | 0.147 | 1.504 |
| HIGH-MASS | 56 | 820 | 3720 | 0.128 | 1.157 |

**Error-weighted fits:**

| Sample | a₀ (km/s)²/kpc | wRMS | a₀/(cH₀/2π) | Delta vs UW |
|--------|----------------|------|-------------|-------------|
| GOLD+i45 | 3545 | 0.118 | 1.102 | -26.7% |
| HIGH-MASS | 3732 | 0.133 | 1.161 | +0.3% |

**Hierarchical (DerSimonian-Laird) results:**

| Sample | n_gal | Hier. mean | tau (dex) | I² | a₀/(cH₀/2π) |
|--------|-------|-----------|-----------|-----|-------------|
| GOLD+i45 | 39 | 3374 | 0.245 | 89.4% | 1.049 |
| GOLD | 44 | 3370 | 0.242 | 89.5% | 1.048 |

**Three-estimator convergence (GOLD+i45):**

| Estimator | a₀ | log(a₀) | a₀/(cH₀/2π) |
|-----------|-----|---------|-------------|
| Unweighted | 4837 | 3.685 | 1.504 |
| Weighted | 3545 | 3.550 | 1.102 |
| Hierarchical | 3374 | 3.528 | 1.049 |

**Uncertainty Budget (v3.0, 6 components):**

| Source | ±σ (dex) |
|--------|----------|
| Estimator method (UW/W/Hier) | ±0.078 |
| Sample selection | ±0.058 |
| Υ★ (±15%) | ±0.066 |
| Inclination cut | ±0.106 |
| Dynamic range cut | ±0.044 |
| Hier. tau (intrinsic scatter) | ±0.245 |
| **TOTAL** | **±0.295** |

## Key Insight (v3.0)

The transition scale EXISTS and is robust. GOLD+i45 is the recommended sample.
HEADLINE number = hierarchical a₀ = 3374 (most defensible estimate).

Number hierarchy:
- HEADLINE: hierarchical a₀ = 3374 (km/s)²/kpc = 1.09e-10 m/s²
- SUPPORTING: weighted a₀ = 3545 (km/s)²/kpc = 1.15e-10 m/s²
- DIAGNOSTIC/RAW: unweighted a₀ = 4837 (biased high, not headline)

Key findings:
- ALL 6 residual audit variables clean (inclination resolved with inc>=45)
- Total uncertainty: ±0.295 dex (dominated by between-galaxy tau = 0.245)
- Literature (3703 = 1.20e-10 m/s²) within 1σ band
- Hierarchical a₀/(cH₀/2π) = 1.049 — suggestive, not compelling
- I² = 89.4% — substantial real heterogeneity between galaxies
- Low-mass a₀ divergence is a fitting artifact (limited dynamic range)

NOT CLAIMED: a₀ universal exact constant, dark matter solved,
MOND proved, cosmological origin established.

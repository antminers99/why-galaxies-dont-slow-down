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
| a₀ (v4.0 marg. hierarchical) | 3633 | (km/s)²/kpc |
| a₀ (v4.0 marg. hier. median) | 3848 | (km/s)²/kpc |
| a₀ (v3.0 hier. baseline) | 3374 | (km/s)²/kpc |
| a₀ (v3.0 weighted baseline) | 3545 | (km/s)²/kpc |
| a₀ (v3.0 UW diagnostic) | 4837 | (km/s)²/kpc |
| a₀ defensible range (Phase 6) | 2800-4300 | (km/s)²/kpc |
| a₀ uncertainty (v4.0) | ±0.304 | dex |
| a₀ 68% CI (v4.0) | [1846, 7149] | (km/s)²/kpc |
| Hier. tau (v4.0 marginalized) | 0.291 | dex |
| Hier. I² (v4.0) | 92.4 | % |
| c | 299,792,458 | m/s |
| H₀ | 67.4 | km/s/Mpc |
| cH₀ | 6.548×10⁻¹⁰ | m/s² |
| cH₀/2π | 1.042×10⁻¹⁰ | m/s² |
| a₀/(cH₀/2π) | 1.130 (v4.0), 1.049 (v3.0 Hier), 1.504 (v3.0 UW) | — |
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

### Pipeline v4.0.0 Results (April 2026) — CURRENT DEFINITIVE

Script: `scripts/definitive-v4.cjs`
Results: `public/definitive-v4-results.json`

**v4.0 upgrades:** Full sample (3755 pts, not subsampled), per-galaxy Y* marginalization (13 Y* x 5 dD grid = 65 evaluations/galaxy), distance marginalization.

**v4.0 Marginalized Hierarchical (DerSimonian-Laird):**

| Sample | n_gal | a₀(marg) | tau (dex) | I² | a₀/(cH₀/2π) |
|--------|-------|----------|-----------|-----|-------------|
| ALL | 165 | 2630 | 0.449 | 94.9% | 0.818 |
| CLEAN | 84 | 3616 | 0.280 | 91.1% | 1.125 |
| GOLD | 64 | 3626 | 0.270 | 92.0% | 1.128 |
| GOLD+i45 | 59 | 3633 | 0.291 | 92.4% | 1.130 |
| HIGH-MASS | 56 | 3405 | 0.280 | 92.5% | 1.059 |

**Ratio stability test (v4.0, GOLD+i45 splits):**
- 11 splits tested, ratio range: 0.875 - 1.469
- Mean ratio: 1.149, Std: 0.162, CV: 14.1%
- Verdict: MODERATE variation — suggestive but not compelling

**Uncertainty Budget (v4.0, 7 components):**

| Source | v3.0 (dex) | v4.0 (dex) |
|--------|----------|----------|
| Estimator method | ±0.078 | ±0.078 |
| Sample selection | ±0.058 | ±0.014 |
| Y_star (fixed/marg) | ±0.066 | ±0.020 |
| Inclination cut | ±0.106 | ±0.021 |
| Distance (marg) | (n/a) | ±0.015 |
| Dynamic range cut | ±0.044 | ±0.044 |
| Hier. tau (intrinsic) | ±0.245 | ±0.291 |
| **TOTAL** | **±0.295** | **±0.304** |

### Pipeline v3.0.0 Results (archived baseline)

Script: `scripts/final-measurement-v3.cjs`
Results: `public/gold-standard-results.json` (v3.0.0)

v3.0 GOLD+i45: UW=4837, W=3545, Hier=3374, tau=0.245, I²=89.4%
(2062 subsampled points, fixed Y*=0.5, 47 galaxies)

## Phase 5 v5.1 Results (Kinematic Contamination Audit)

Script: `scripts/phase5-kinematic-audit.cjs` (v5.1.0)
Results: `public/phase5-kinematic-results.json`
SPARC master table: `public/sparc-table.json` (175 galaxies with Q, T, fD)

Uses CORRECT v4 GOLD+i45 sample (59 galaxies). Reproduces v4 baseline exactly.

| Test | Split (dex) | Verdict |
|------|-------------|---------|
| Inner vs outer | 0.006 | NEGLIGIBLE |
| SPARC Q=1 vs Q=2 | 0.016 | NEGLIGIBLE |
| Early vs late morph | 0.048 | SMALL |
| Q_kin split | 0.002 | NEGLIGIBLE |
| Distance method | 0.105 | MODERATE (27%) |
| Pressure correction | +0.033 | SMALL (+7.8%) |
| Combined outer+pressure | a₀=3274, ratio=1.018 | tau=0.312 |

Residual correlations: r(a₀, distance) = -0.31 (significant), r(a₀, inc) = 0.16 (moderate).
Red-flag splits: NOT resolved by pressure correction (they worsen).

## Phase 6 v6.0 Results (Matched-Sample Split Resolution)

Script: `scripts/phase6-matched-splits.cjs` (v6.0.0)
Results: `public/phase6-matched-results.json`

| Split | Raw (dex) | Matched (dex) | Paired t | JK t | Verdict |
|-------|-----------|---------------|----------|------|---------|
| Distance | 0.105 | 0.125 | 1.07 | 1.03 | Underpowered |
| Inclination | 0.153 | 0.172 | 1.24 | 1.53 | Underpowered |
| Velocity | 0.135 | 0.061 | 0.77 | 0.77 | Reduced by matching |

Power analysis (galaxies needed to reach t=2.0 / 80% power):
- Inclination: ~101 / ~198 (closest to detection, 1.7x current)
- Distance: ~227 / ~444 (3.8x current)
- Velocity: ~397 / ~777 (6.7x current)

Phase 6 verdict: UNDERDETERMINED by current sample size.
Splits not refuted, not confirmed. Problem is sample size, not pipeline.

Distance binning (continuous gradient):
- D < 10 Mpc: n=18, a₀=4116
- 10-20 Mpc: n=19, a₀=4178
- 20-40 Mpc: n=8, a₀=3402
- D >= 40 Mpc: n=14, a₀=2411

Inclination binning:
- 45-54: n=10, a₀=2271
- 55-64: n=18, a₀=3655
- 65-74: n=15, a₀=4800 (peak)
- >=75: n=16, a₀=3725

Velocity binning:
- 50-80: n=6, a₀=2781
- 80-120: n=5, a₀=6721 (outlier, few galaxies)
- 120-180: n=12, a₀=3415
- 180-250: n=22, a₀=3111
- >=250: n=14, a₀=4319

## Phase 7 v7.0 Results (Anchor-Sample Refit)

Script: `scripts/phase7-anchor-refit.cjs` (v7.0.0)
Results: `public/phase7-anchor-results.json`

Anchor = GOLD+i45 + precise distance (fD=2,3,5) + Q=1 + inc>=60 → 12 galaxies

| Sample | n | a₀ | tau | I² | ratio |
|--------|---|-----|------|-----|-------|
| GOLD+i45 (baseline) | 59 | 3633 | 0.291 | 92.4% | 1.130 |
| Q=1 only | 40 | 3632 | 0.273 | 92.1% | 1.130 |
| Precise dist only | 17 | 4243 | 0.257 | 93.0% | 1.320 |
| Precise+Q1+inc≥45 | 14 | 4124 | 0.248 | 93.3% | 1.283 |
| ANCHOR (strictest) | 12 | 4357 | 0.265 | 93.9% | 1.355 |

Residual correlations:
- GOLD+i45: r(a₀, logD) = -0.235 (marginal)
- ANCHOR: r(a₀, logD) = -0.036 (negligible) → distance correlation DISAPPEARS

Key diagnostics:
- tau: -9% → STABLE (heterogeneity is INTRINSIC)
- a₀: +20% → INCREASES (Hubble-flow biases a₀ downward)
- r(a₀, D): -0.24 → -0.04 → distance split = Hubble-flow systematic

## Key Insight (v4.0 + Phase 5-7)

The transition scale EXISTS and is robust. GOLD+i45 is the recommended sample.
HEADLINE number = v4.0 marginalized hierarchical a₀ = 3633 (most defensible).

Number hierarchy:
- HEADLINE: v4.0 marg. hier. a₀ = 3633 (km/s)²/kpc = 1.18e-10 m/s²
- BASELINE: v3.0 hier. a₀ = 3374 (km/s)²/kpc = 1.09e-10 m/s²
- COMBINED BEST: outer+pressure a₀ = 3274, ratio = 1.018 to cH₀/2π
- DEFENSIBLE RANGE: 2800-4300 (km/s)²/kpc = (0.9-1.4)e-10 m/s²
- DIAGNOSTIC/RAW: v3.0 UW a₀ = 4837 (biased high, not headline)

Key findings:
- v4.0 a₀ is within 2% of literature (3703 = 1.20e-10 m/s²)
- Y* marginalization shifts a₀ +7.7% but does NOT reduce tau
- Population heterogeneity (tau=0.291 dex) is real, not Y*-driven
- Ratio stability: moderate (CV=14.1%), range [0.875, 1.469]
- Total uncertainty: ±0.304 dex (dominated by tau=0.291)
- I² = 92.4% — substantial real heterogeneity
- Low-mass a₀ divergence is a fitting artifact (limited dynamic range)
- Phase 5 v5.1: all kinematic tests PASS, distance method split (27%) remains
- Phase 6: velocity split reduced by confounders, distance/inc PERSIST but underpowered
- Phase 6: ALL splits NOT significant (t < 2.0) — tau swamps them
- Phase 7: tau is INTRINSIC (-9% in anchor) — not measurement error
- Phase 7: distance correlation DISAPPEARS with precise distances (Hubble-flow bias)
- Phase 7: a₀ ~ 4100-4400 with precise distances (headline 3633 biased low ~20%)
- Phase 7: Q flag does NOT drive a₀ (Q=1 alone: 3632 = same as 3633)

NOT CLAIMED: a₀ universal exact constant, dark matter solved,
MOND proved, cosmological origin established.

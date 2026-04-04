---
name: galaxy-summary
description: Quick executive summary of the Galaxy Rotation Curve Analyzer's current definitive result. Load this skill whenever you need a concise recap of the a₀ measurement, its status, and key caveats.
---

# Galaxy Rotation Curve Analysis — Executive Summary

## Current Definitive Result (v4.0 + Phase 5-7)

**Pipeline**: Full sample, per-galaxy Y★ & distance marginalization, DerSimonian-Laird hierarchical model, kinematic contamination audit, matched-sample split resolution, anchor-sample refit.

| Quantity | Value |
|----------|-------|
| a₀ (headline) | 3633 (km/s)²/kpc = 1.18×10⁻¹⁰ m/s² |
| log₁₀(a₀) | 3.560 ± 0.041 (stat) ± 0.304 (total) dex |
| τ (intrinsic scatter) | 0.291 dex |
| I² (heterogeneity) | 92.4% |
| a₀/(cH₀/2π) | 1.130 |
| Sample | GOLD+i≥45°, 59 galaxies, 1789 points |
| Literature (McGaugh+2016) | 3703 (km/s)²/kpc = 1.20×10⁻¹⁰ m/s² |
| Agreement | Within 1σ (Δ = −1.9%) |
| 68% CI | [1846, 7149] (km/s)²/kpc |
| 95% CI | [938, 14067] (km/s)²/kpc |
| Defensible range (Phase 6) | 2800-4300 (km/s)²/kpc = (0.9-1.4)×10⁻¹⁰ m/s² |
| a₀ anchor (precise D, Q=1, i≥60) | 4357 (km/s)²/kpc = 1.41×10⁻¹⁰ m/s² (n=12) |
| tau anchor | 0.265 dex (vs 0.291 baseline, -9%) |
| Anchor ratio | 1.355 |

## One-Paragraph Verdict

We independently recover a₀ ≈ 1.18×10⁻¹⁰ m/s² from 175 SPARC + 22 LITTLE THINGS galaxies, within 2% of the McGaugh+2016 value. The Radial Acceleration Relation is confirmed as a tight empirical law. Phase 5 v5.1 kinematic contamination audit showed a₀ survives non-circular motion tests (0.006 dex) and pressure support corrections (+0.033 dex). Phase 6 matched-sample analysis found: the three remaining red-flag splits have non-trivial effect sizes (distance 0.125 dex, inclination 0.172 dex, velocity 0.061 dex), but NONE reach statistical significance (all t < 2.0). Phase 6 did NOT refute these splits, and did NOT confirm them — it showed the current 59-galaxy sample is UNDERPOWERED to decide. The EXISTENCE of the transition scale is robust. Whether its PRECISE VALUE is strictly universal cannot be determined with current data. The limiting factor is sample size, not methodology. Power analysis: need ~101-397 GOLD+i45 galaxies (1.7-6.7x current) to resolve.

## Phase 6 Key Results (Matched-Sample Split Resolution)

| Split | Raw (dex) | Matched (dex) | JK t-stat | Verdict |
|-------|-----------|---------------|-----------|---------|
| Distance (precise vs HF) | 0.105 (25%) | 0.125 (29%) | 1.03 | Underpowered (need ~227 gal) |
| Inclination (≥60 vs 45-59) | 0.153 (33%) | 0.172 (38%) | 1.53 | Underpowered (need ~101 gal) |
| Velocity (>150 vs 80-150) | 0.135 (35%) | 0.061 (17%) | 0.77 | Reduced by matching (need ~397 gal) |

Key interpretation:
- Phase 6 did NOT kill the splits, and did NOT confirm them
- The splits have non-trivial effect sizes but are NOT statistically significant (all t<2.0)
- This is a SAMPLE SIZE problem: 59 galaxies is underpowered to resolve effects of this size
- Velocity split reduced by matching (confounders), distance/inclination persist
- Power analysis: need 101-397 GOLD+i45 galaxies (currently 59)

## Phase 5 v5.1 Key Results (Kinematic Audit)

| Test | Result | Verdict |
|------|--------|---------|
| Inner vs outer a₀ | 0.006 dex split | Non-circular motions NEGLIGIBLE |
| SPARC Q=1 vs Q=2 | 0.016 dex split | Data quality NO effect |
| Hubble type early vs late | 0.048 dex split | Morphology SMALL effect |
| Distance method split | 0.105 dex (27%) | LARGEST remaining systematic |
| Pressure correction | +7.8% shift, tau unchanged | Minor effect only |

## Decision Framework — ALL 4 BOXES CHECKED

1. Non-circular motions ✓ (Phase 5 — inner-outer split 0.006 dex)
2. Pressure support ✓ (Phase 5 — +0.033 dex shift, tau unchanged)
3. Y★ + distance + inclination ✓ (v4.0 — per-galaxy marginalization)
4. Weighting + covariance ✓ (v4.0 — hierarchical DL model)

## Phase 7 Key Results (Anchor-Sample Refit)

| Diagnostic | GOLD+i45 (n=59) | Anchor (n=12) | Interpretation |
|------------|----------------|---------------|----------------|
| a₀ | 3633 | 4357 (+20%) | Higher with precise distances |
| tau | 0.291 | 0.265 (-9%) | STABLE — heterogeneity is INTRINSIC |
| r(a₀, distance) | -0.235 (t=-1.83) | -0.036 (t=-0.11) | Distance correlation DISAPPEARS |
| r(a₀, inc) | +0.159 | -0.467 | Sign REVERSES (sample-dependent) |

Key findings:
- tau is intrinsic: even cleanest galaxies scatter ~0.27 dex (not measurement error)
- Distance correlation disappears in anchor (r: -0.24→-0.04)
- Q flag does NOT drive a₀: Q=1 alone (n=40) gives a₀=3632 (unchanged)

## Phase 8 v8.0 Results (Distance-Method Meta-Regression)

Phase 8 meta-regression OVERRIDES the Phase 7 "Hubble-flow bias" interpretation:
- b(isPrecise) = +0.061 dex (t=0.56, permutation p=0.55) → NOT SIGNIFICANT
- Raw split (0.105 dex) shrinks to 0.061 after controlling for covariates
- PRIMARY CONFOUNDER: inclination (b=+0.465 per inc/90, t=1.26)
- Phase 7 anchor a₀ increase (3633→4357) was driven by inc≥60 cut, NOT by precise distances
- Propensity-matched paired diff: +0.174 dex (t=1.79) — moderate effect, NOT significant
- LOO: b(isPrecise) never reaches significance (range [0.016, 0.110])
- Headline a₀=3633 is NOT demonstrably biased by any measured systematic

## Phase 9 v9.0 Results (Benchmark-to-Literature Replication)

NO TENSION with literature. Different analyses measure different scatter components:
- Point-level (within-galaxy) scatter: 0.122 dex ≈ McGaugh+2016 (0.13 dex) ✓
- Galaxy-level (between-galaxy) tau: 0.29 dex — DIFFERENT QUANTITY from literature's scatter
- Variance decomposition (ANOVA): total²=0.0239 = within²=0.0148 + between²=0.0091 (point-level)
- SD of per-galaxy logA0 = 0.275 dex (DIFFERENT from between-galaxy point-level component 0.096 dex)
- Ablation: tau not driven by Y* marg (+4%), quality cuts, or dynamic range
- F-test: 39/56 (70%) galaxies individually reject global a₀ (consistent with Rodrigues+2018)
- Position: complements literature — RAR is tight at point level, heterogeneous at galaxy level

## Phase 10 v10.0 Results (Second-Parameter Search)

10 candidates tested. 3 significant but MODEST:
- eta_rot (RC shape): r=0.334, p=0.010, dAIC=-4.6, tau drops 5.7% — STRONGEST
- T (Hubble type): r=0.278, p=0.043, tau drops 4.4%
- S_out (outer slope): r=-0.418, p=0.045, but n=23 only
- Surface brightness, inclination, distance method: NEGLIGIBLE
- 94% of heterogeneity remains UNEXPLAINED

## Phase 13 v13.0 Results (Formal Model Decision)

M3 (per-galaxy a0) wins ALL criteria decisively:
- 46.5% CV RMS improvement over M0 (universal a0)
- ΔAIC = 1686 (overwhelming evidence)
- ΔBIC = 1395
- M0 and M1 fail posterior predictive checks
- eta_rot (M2) worsens out-of-sample prediction
- VERDICT: Universal a0 is INSUFFICIENT — M3 decisively wins

## Phase 14 v14.0 Results (Decompose M3)

M3 decomposed into offset (a0 shift) vs shape (gamma variation):
- OFFSET: 86% of M3's advantage — galaxies differ in SCALE
- SHAPE: 14% of M3's advantage — transition curvature varies
- gamma and a0 are INDEPENDENT (r = -0.14)
- Mean gamma = 0.96 ± 0.30 (standard = 0.5)
- Best predictor of a0 variation: log(MHI) (r = -0.388, LOO-CV R² = 0.108)
- Equation: a0,i = a0 - 0.174 × (log MHI - mean)
- Top-2 model (MHI + n_points): R² = 0.28, tau-reduction = 28%

## Phase 15 v15.0 Results (Model Ladder Decision Map)

Tested parametric models a0,i = a0 + f(Xi) against M3 (free per-galaxy):
- M0 (universal): gap closed = 0%
- M1 (MHI only): gap closed = 7.4%
- M5 (kitchen sink, 8 predictors): gap closed = 15.7%
- M6 (per-galaxy free): gap closed = 100%
- VERDICT: a0 GENUINELY VARIES — no f(Xi) with measured properties can replace M3
- 84.3% of variation remains UNEXPLAINED
- Either an unmeasured variable drives variation, or a0 is physically galaxy-dependent

## Key Files

- v4.0 pipeline: `artifacts/galaxy-analyzer/scripts/definitive-v4.cjs`
- Phase 5 audit: `artifacts/galaxy-analyzer/scripts/phase5-kinematic-audit.cjs` (v5.1.0)
- Phase 6 splits: `artifacts/galaxy-analyzer/scripts/phase6-matched-splits.cjs` (v6.0.0)
- Phase 7 anchor: `artifacts/galaxy-analyzer/scripts/phase7-anchor-refit.cjs` (v7.0.0)
- Phase 13 model decision: `artifacts/galaxy-analyzer/scripts/phase13-model-decision.cjs`
- Phase 14 decompose: `artifacts/galaxy-analyzer/scripts/phase14-decompose-m3.cjs`
- Phase 15 ladder: `artifacts/galaxy-analyzer/scripts/phase15-model-ladder.cjs`
- v4.0 results: `artifacts/galaxy-analyzer/public/definitive-v4-results.json`
- Phase 13 results: `artifacts/galaxy-analyzer/public/phase13-model-decision.json`
- Phase 14 results: `artifacts/galaxy-analyzer/public/phase14-decompose-m3.json`
- Phase 15 results: `artifacts/galaxy-analyzer/public/phase15-model-ladder.json`
- SPARC master table: `artifacts/galaxy-analyzer/public/sparc-table.json` (175 galaxies)
- Full report: `artifacts/galaxy-analyzer/diagnostic-report.txt` (33 sections)
- Experiment log: `artifacts/galaxy-analyzer/experiment-log.txt` (44 sections)

## Phase 16-52 Results (Sequential Door Investigation — COMPLETE)

ALL SIX DOORS CLOSED. 37 phases, ~138 proxies tested across 56 galaxies.

Five variables survive as independently informative:
1. log(MHI) — gas mass (ZERO circularity, Door 1: Galaxy History)
2. rcWiggliness — RC bumpiness (MED circularity, Door 5: Gas/Kinematics)
3. envCode — environment threshold (ZERO circularity, Door 2: Environment)
4. Sigma0_bar — baryon surface density (ZERO circularity, Door 4: Dark Matter Halo)
5. meanRun — RAR residual coherence (MED circularity, Door 6: Internal Structure)

Model ladder (LOO gap-closed):
- Conservative BL (4 vars): 27.8%
- Working BL (5 vars): 38.1%
- Per-galaxy free: 100%
- Unexplained: 61.9%

Key negative results:
- NOT a resolution artifact (distance/beam/Q all fail)
- No cosmic-scale correlate (all cosmic context fails)
- Classical morphology irrelevant (Hubble type, bulge ratio fail)
- Stellar populations absorbed by MHI
- DM halo parameters fail independently

## Status

PROJECT COMPLETE. All six investigation doors closed (Phases 16-52). Five variables confirmed, explaining 38.1% of the M0→M6 gap (LOO). tau=0.22 dex robust, never broken. 61.9% of variation remains unexplained. ~138 proxies tested, 95% rejection rate.

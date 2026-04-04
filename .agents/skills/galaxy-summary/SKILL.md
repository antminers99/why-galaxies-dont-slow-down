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
- Variance decomposition: total² = within² + between² → compatible
- Ablation: tau not driven by Y* marg (+4%), quality cuts, or dynamic range
- F-test: 39/56 (70%) galaxies individually reject global a₀ (consistent with Rodrigues+2018)
- Position: complements literature — RAR is tight at point level, heterogeneous at galaxy level

## Key Files

- v4.0 pipeline: `artifacts/galaxy-analyzer/scripts/definitive-v4.cjs`
- Phase 5 audit: `artifacts/galaxy-analyzer/scripts/phase5-kinematic-audit.cjs` (v5.1.0)
- Phase 6 splits: `artifacts/galaxy-analyzer/scripts/phase6-matched-splits.cjs` (v6.0.0)
- Phase 7 anchor: `artifacts/galaxy-analyzer/scripts/phase7-anchor-refit.cjs` (v7.0.0)
- v4.0 results: `artifacts/galaxy-analyzer/public/definitive-v4-results.json`
- Phase 6 results: `artifacts/galaxy-analyzer/public/phase6-matched-results.json`
- Phase 7 results: `artifacts/galaxy-analyzer/public/phase7-anchor-results.json`
- Phase 8 metareg: `artifacts/galaxy-analyzer/scripts/phase8-distance-metareg.cjs` (v8.0.0)
- Phase 8 results: `artifacts/galaxy-analyzer/public/phase8-metareg-results.json`
- Phase 9 benchmark: `artifacts/galaxy-analyzer/scripts/phase9-benchmark-replication.cjs` (v9.0.0)
- Phase 9 results: `artifacts/galaxy-analyzer/public/phase9-benchmark-results.json`
- SPARC master table: `artifacts/galaxy-analyzer/public/sparc-table.json` (175 galaxies)
- Full report: `artifacts/galaxy-analyzer/diagnostic-report.txt` (25 sections)

## Status

All four systematic checkboxes passed. The EXISTENCE of a transition scale is observationally robust. Phase 9 benchmark replication resolved the apparent tension with literature: our tau (~0.29 dex) and the literature's scatter (~0.04-0.13 dex) are DIFFERENT QUANTITIES (galaxy-level vs point-level). Both are correct. Within-galaxy scatter matches McGaugh+2016 (0.12 vs 0.13 dex). Between-galaxy heterogeneity is intrinsic, not a pipeline artifact. 70% of galaxies individually reject the global a₀. Headline a₀=3633 is not demonstrably biased. The RAR is tight at the data-point level but heterogeneous at the galaxy level.

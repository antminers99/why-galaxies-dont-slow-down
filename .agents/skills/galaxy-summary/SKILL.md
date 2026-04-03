---
name: galaxy-summary
description: Quick executive summary of the Galaxy Rotation Curve Analyzer's current definitive result. Load this skill whenever you need a concise recap of the a₀ measurement, its status, and key caveats.
---

# Galaxy Rotation Curve Analysis — Executive Summary

## Current Definitive Result (v4.0 + Phase 5 v5.1)

**Pipeline**: Full sample, per-galaxy Y★ & distance marginalization, DerSimonian-Laird hierarchical model, kinematic contamination audit with SPARC master table metadata.

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
| Combined best (outer+pressure) | 3274, ratio = 1.018 |

## One-Paragraph Verdict

We independently recover a₀ ≈ 1.18×10⁻¹⁰ m/s² from 175 SPARC + 22 LITTLE THINGS galaxies, within 2% of the McGaugh+2016 value. The Radial Acceleration Relation is confirmed as a tight empirical law. Phase 5 v5.1 kinematic contamination audit (corrected to use proper GOLD+i45 sample of 59 galaxies, matching v4) showed that a₀ survives non-circular motion tests (inner-outer split = 0.006 dex) and pressure support corrections (shift = +0.033 dex). SPARC quality flags (Q=1 vs Q=2: 0.016 dex) and Hubble morphology (early vs late: 0.048 dex) produce negligible effects. The largest remaining systematic is the distance method split (0.105 dex): precise TRGB/Cepheid distances give a₀ ~ 4243 vs Hubble flow a₀ ~ 3331. Red-flag splits (distance, inclination, velocity) persist and are NOT resolved by pressure correction — they reflect intrinsic sample systematics. All four decision framework checkboxes are passed.

## Phase 5 v5.1 Key Results (Corrected Kinematic Audit)

| Test | Result | Verdict |
|------|--------|---------|
| Inner vs outer a₀ | 0.006 dex split | Non-circular motions NEGLIGIBLE |
| SPARC Q=1 vs Q=2 | 0.016 dex split | Data quality NO effect |
| Hubble type early vs late | 0.048 dex split | Morphology SMALL effect |
| Q_kin split (within GOLD+i45) | 0.002 dex split | Kinematic quality NO effect |
| Distance method split | 0.105 dex (27%) | LARGEST remaining systematic |
| Pressure correction | +7.8% shift, tau unchanged | Minor effect only |
| Outer + pressure combined | a₀=3274, ratio=1.018 | Closest to cH₀/2π |
| Residual r(a₀, distance) | -0.31 (significant) | Distance systematic |

## v5.0 → v5.1 Corrections

The v5.0 Phase 5 report used an INCORRECT sample of 126 galaxies (Q_kin>=0.7, inc>=45, n>=5) that did NOT match the v4 GOLD+i45 definition (Vmax>=50, inc>=45, n>=5, gRange>=1.0). Key v5.0 claims that were WRONG:
- "Distance split ELIMINATED (35% → 3.6%)" — WRONG, it persists at ~31-47%
- "Tau drops 75% with combined cleaning" — WRONG, tau is stable at ~0.29-0.31
- "Gas-rich tau=0.259" — was an artifact of the wrong sample

## Decision Framework — ALL 4 BOXES CHECKED

1. Non-circular motions ✓ (Phase 5 — inner-outer split 0.006 dex)
2. Pressure support ✓ (Phase 5 — +0.033 dex shift, tau unchanged)
3. Y★ + distance + inclination ✓ (v4.0 — per-galaxy marginalization)
4. Weighting + covariance ✓ (v4.0 — hierarchical DL model)

## Uncertainty Budget (7 components)

| Component | Contribution (dex) |
|-----------|-------------------|
| Estimator choice | 0.078 |
| Sample selection | 0.014 |
| Y★ (residual after marginalization) | 0.020 |
| Inclination (residual) | 0.021 |
| Distance (residual) | 0.015 |
| Bin range | 0.044 |
| Intrinsic scatter τ | 0.291 |
| **Total (quadrature)** | **0.304** |

## Key Files

- v4.0 pipeline: `artifacts/galaxy-analyzer/scripts/definitive-v4.cjs`
- Phase 5 audit: `artifacts/galaxy-analyzer/scripts/phase5-kinematic-audit.cjs` (v5.1.0)
- v4.0 results: `artifacts/galaxy-analyzer/public/definitive-v4-results.json`
- Phase 5 results: `artifacts/galaxy-analyzer/public/phase5-kinematic-results.json` (v5.1.0)
- SPARC master table: `artifacts/galaxy-analyzer/public/sparc-table.json` (175 galaxies)
- Full report: `artifacts/galaxy-analyzer/diagnostic-report.txt` (21 sections)

## Status

All four systematic checkboxes are passed. The transition scale a₀ ≈ 1.2×10⁻¹⁰ m/s² is observationally robust. Remaining open: distance method split (27%), red-flag splits persist, cosmological ratio interpretation, and whether τ ≈ 0.29 can be further reduced.

---
name: galaxy-summary
description: Quick executive summary of the Galaxy Rotation Curve Analyzer's current definitive result. Load this skill whenever you need a concise recap of the a₀ measurement, its status, and key caveats.
---

# Galaxy Rotation Curve Analysis — Executive Summary

## Current Definitive Result (v4.0 + Phase 5)

**Pipeline**: Full sample, per-galaxy Y★ & distance marginalization, DerSimonian-Laird hierarchical model, kinematic contamination audit.

| Quantity | Value |
|----------|-------|
| a₀ (headline) | 3633 (km/s)²/kpc = 1.18×10⁻¹⁰ m/s² |
| log₁₀(a₀) | 3.560 ± 0.041 (stat) ± 0.304 (total) dex |
| τ (intrinsic scatter) | 0.291 dex (full); 0.178 dex (outer+corrected) |
| I² (heterogeneity) | 92.4% |
| a₀/(cH₀/2π) | 1.130 |
| Sample | GOLD+i≥45°, 59 galaxies, 1789 points |
| Literature (McGaugh+2016) | 3703 (km/s)²/kpc = 1.20×10⁻¹⁰ m/s² |
| Agreement | Within 1σ (Δ = −1.9%) |
| 68% CI | [1846, 7149] (km/s)²/kpc |
| 95% CI | [938, 14067] (km/s)²/kpc |

## One-Paragraph Verdict

We independently recover a₀ ≈ 1.18×10⁻¹⁰ m/s² from 175 SPARC + 22 LITTLE THINGS galaxies, within 2% of the McGaugh+2016 value. The Radial Acceleration Relation is confirmed as a tight empirical law. Phase 5 kinematic contamination audit showed that a₀ survives both non-circular motion tests (inner-outer split = 0.033 dex) and pressure support corrections (shift = +0.053 dex). When both contaminants are cleaned (outer radii + pressure correction), τ drops 75% (0.705 → 0.178) while a₀ shifts only 2.9% — the transition scale is robust. All four systematic checkboxes are now passed. We do NOT claim a₀ is a universal exact constant, that dark matter is ruled out, that MOND is proved, or that the cosmological coincidence is established.

## Phase 5 Key Results (Kinematic Contamination Audit)

| Test | Result | Verdict |
|------|--------|---------|
| Inner vs outer a₀ | 0.033 dex split | Non-circular motions NOT dominant |
| Pressure support correction | +12.9% shift, τ drops 0.068 | Correction helps, does not destroy signal |
| Outer + pressure combined | a₀ shifts 2.9%, τ drops 75% | **Transition scale is robust** |
| Distance red flag | 35% → 3.6% | ELIMINATED |
| Inclination red flag | 42% → 16.4% | REDUCED 61% |
| Velocity red flag | 36% → 29.8% | SLIGHTLY IMPROVED |

## Decision Framework — ALL 4 BOXES CHECKED

1. Non-circular motions ✓ (Phase 5 — inner-outer split tiny)
2. Pressure support ✓ (Phase 5 — correction helps, signal survives)
3. Y★ + distance + inclination ✓ (v4.0 — per-galaxy marginalization)
4. Weighting + covariance ✓ (v4.0 — hierarchical DL model)

## What Changed from v3 → v4 → Phase 5

| | v3.0 | v4.0 | Phase 5 (outer+corr) |
|---|---|---|---|
| Data | 2062 pts (subsampled) | 3755 pts (full) | 3755 (outer subset) |
| Y★ | Fixed 0.5 | Marginalized | Marginalized |
| Pressure corr | No | No | Yes (σ=10 km/s) |
| a₀ | 3374 | 3633 | ~2961 (wider sample) |
| τ | 0.245 | 0.291 | 0.178 |

## Uncertainty Budget (7 components)

| Component | Contribution (dex) |
|-----------|-------------------|
| Estimator choice | 0.078 |
| Sample selection | 0.014 |
| Y★ (residual after marginalization) | 0.020 |
| Inclination (residual) | 0.021 |
| Distance (residual) | 0.015 |
| Bin range | 0.044 |
| Intrinsic scatter τ | 0.291 (0.178 cleaned) |
| **Total (quadrature)** | **0.304** |

## Key Files

- v4.0 pipeline: `artifacts/galaxy-analyzer/scripts/definitive-v4.cjs`
- Phase 5 audit: `artifacts/galaxy-analyzer/scripts/phase5-kinematic-audit.cjs`
- v4.0 results: `artifacts/galaxy-analyzer/public/definitive-v4-results.json`
- Phase 5 results: `artifacts/galaxy-analyzer/public/phase5-kinematic-results.json`
- Full report: `artifacts/galaxy-analyzer/diagnostic-report.txt` (21 sections)

## Status

All four systematic checkboxes are passed. The transition scale a₀ ≈ 1.2×10⁻¹⁰ m/s² is observationally robust. Remaining open: velocity split (30%), cosmological ratio interpretation, and whether τ ≈ 0.18 is irreducible.

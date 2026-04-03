---
name: galaxy-summary
description: Quick executive summary of the Galaxy Rotation Curve Analyzer's current definitive result. Load this skill whenever you need a concise recap of the a₀ measurement, its status, and key caveats.
---

# Galaxy Rotation Curve Analysis — Executive Summary

## Current Definitive Result (v4.0)

**Pipeline**: Full sample, per-galaxy Y★ & distance marginalization, DerSimonian-Laird hierarchical model.

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

## One-Paragraph Verdict

We independently recover a₀ ≈ 1.18×10⁻¹⁰ m/s² from 175 SPARC + 22 LITTLE THINGS galaxies, within 2% of the McGaugh+2016 value. The Radial Acceleration Relation is confirmed as a tight empirical law. However, population heterogeneity remains large (τ ≈ 0.29 dex, I² ≈ 92%), the ratio a₀/(cH₀/2π) ≈ 1.13 is suggestive but not compelling (CV = 14% across subsample splits), and the total uncertainty budget is ±0.30 dex — dominated by intrinsic scatter, not systematics. We do NOT claim a₀ is a universal exact constant, that dark matter is ruled out, that MOND is proved, or that the cosmological coincidence is established.

## What Changed from v3 → v4

| | v3.0 (baseline) | v4.0 (definitive) |
|---|---|---|
| Data | 2062 pts (subsampled) | 3755 pts (full) |
| Y★ | Fixed 0.5 | Marginalized [0.2, 0.8] |
| Distance | Fixed | Marginalized ±0.1 dex |
| a₀ | 3374 | 3633 (+7.7%) |
| τ | 0.245 | 0.291 |
| Key finding | — | Y★ marginalization shifts a₀ toward literature but does NOT reduce τ |

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

- Pipeline script: `artifacts/galaxy-analyzer/scripts/definitive-v4.cjs`
- Results JSON: `artifacts/galaxy-analyzer/public/definitive-v4-results.json`
- Full report: `artifacts/galaxy-analyzer/diagnostic-report.txt` (19 sections, ~1488 lines)

## Status

Phase 1 (diagnostics) and Phase 2 (measurement) are complete. The primary open questions are whether the intrinsic scatter is physical or driven by uncontrolled systematics (non-circular motions, pressure support, gas flaring), and whether the cosmological ratio has deeper significance.

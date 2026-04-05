---
name: galaxy-summary
description: Quick executive summary of the Galaxy Rotation Curve Analyzer's current definitive result. Load this skill whenever you need a concise recap of the a₀ measurement, its status, and key caveats.
---

# Galaxy Rotation Curve Analysis — Executive Summary

## Current Definitive Result (Phase 60 — Final Death Match)

**Pipeline**: Full sample, per-galaxy Y★ & distance marginalization, DerSimonian-Laird hierarchical model, 60-phase investigation program.

| Quantity | Value |
|----------|-------|
| a₀ (headline) | 3633 (km/s)²/kpc = 1.18×10⁻¹⁰ m/s² |
| log₁₀(a₀) mean | 3.563 ± 0.264 dex |
| τ (intrinsic scatter) | 0.264 dex (robust, never broken) |
| Sample | 56 galaxies, 1483 points |
| Literature (McGaugh+2016) | 3703 (km/s)²/kpc = 1.20×10⁻¹⁰ m/s² |
| Agreement | Within 1σ (Δ = −1.9%) |

## Final Model Comparison (Phase 60)

| Model | Parameters | LOO Gap-Closed | Verdict |
|-------|-----------|----------------|---------|
| M0: Universal a₀ | 0 | 0.0% | **REJECTED** |
| M1a: Conservative BL | 4 | 27.8% | Good (publication-grade) |
| M1b: Extended BL | 5 | 38.1% | Better (internal analysis) |
| M1c: Maximal | 6 | 42.9% | **WINNER** |
| M2: Per-galaxy a₀ | 56 | 100.0% | Overfitting |

## Confirmed Variables (6 in maximal model)

| Variable | Domain | Circularity | Baseline |
|----------|--------|-------------|----------|
| log(MHI) | Gas content | ZERO | A+B+C |
| rcWiggliness | RC texture | MED | A+B+C |
| envCode | Environment | ZERO | A+B+C |
| Sigma0_bar | Baryon density | ZERO | A+B+C |
| meanRun | RAR residual coherence | MED | B+C |
| innBarF | Baryon dominance fraction | MED | C only |

## One-Paragraph Verdict

The a₀ correlation structure splits at Vflat≈70 km/s. **High-Vflat (N=112)**: structured 5-predictor law (R²=0.253, LOO gap=8.1%, p<0.001). **Low-Vflat (N=63)**: appears chaotic (sd=0.77 dex) — Phase 89 shows the SCATTER AMPLITUDE is quality-dominated (Q=1 dwarfs sd=0.325 ≈ high regime). BUT Phase 90 shows the CORRELATION STRUCTURE is genuinely different: even among Q=1 dwarfs (N=22), no multi-axis law emerges (best physics LOO gap=0.5%), 2/4 key signs reverse vs high regime, cross-regime transfer fails (−37% to −44%). The Vflat≈70 boundary is BOTH a data-quality transition (inflating scatter) AND a physical transition (changing which properties correlate with a₀). Dwarfs and spirals share comparable scatter when quality-controlled but do NOT share the same predictive law.

## Decisive Negative Results (Phases 55-59)

| Phase | Test | Result |
|-------|------|--------|
| 55 | 10 pairwise interactions | ALL FAIL — variables act additively |
| 57 | 14 kinematic proxies | ALL FAIL — rcWig+meanRun capture all 1D info |
| 58 | 12 environmental processing proxies | ALL FAIL — envCode irreducible |
| 59 | 11 stellar M/L proxies | 2 survive (barGravD, innBarF) but don't absorb Σ₀ |

## Key Numbers (Locked)

- **HEADLINE**: a₀ = 3633 (km/s)²/kpc
- **tau**: 0.264 dex (robust across all subsamples, 95% CI [0.223, 0.299])
- **Conservative BL**: 27.8% gap closed (4 vars, mostly ZERO circ)
- **Extended BL**: 38.1% gap closed (5 vars)
- **Maximal BL**: 42.9% gap closed (6 vars)
- **Unexplained**: 57.1% (maximal) to 72.2% (conservative)
- **~138 proxies tested**, 95% rejection rate

## Equation (Conservative Baseline A)

log(a₀) = 3.13 − 0.23·logMHI + 2.29·wig − 0.12·env + 0.18·logΣ₀

## Key Files

- Phase 56 frozen baselines: `public/phase56-frozen-baselines.json`
- Phase 60 death match: `public/phase60-death-match.json`
- Diagnostic report: `diagnostic-report.txt` (82 sections, 5814 lines)
- Experiment log: `experiment-log.txt` (93 sections, 2285 lines)

## Status

PROGRAM AT PHASE 90. Two-regime structure fully characterized. Key finding: the Vflat≈70 boundary is BOTH physical AND data-quality. Scatter amplitude is quality-dominated (Phase 89: Q=1 sd=0.325 ≈ high regime). But correlation structure is genuinely different (Phase 90: no law in Q=1 dwarfs, 2/4 signs reverse, cross-regime transfer fails −37% to −44% even Q=1-only). Inclination signal (r=−0.473, p=0.025) likely systematic artifact.

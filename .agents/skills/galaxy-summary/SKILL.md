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

The a₀ correlation structure splits at Vflat≈70 km/s. **High-Vflat (N=112)**: structured 5-predictor law (R²=0.253, LOO gap=8.1%, p<0.001). **Low-Vflat (N=63)**: no law even in Q=1 dwarfs (Phase 90). Phase 91 explains WHY: four convergent transitions at Vflat≈70 — (1) fgas crosses 50% (disk loses dominance), (2) gbar range halves (a₀ degenerate), (3) deep MOND regime (RAR insensitive to a₀), (4) 73% rising RCs (Vflat undefined). Per-galaxy a₀ is OPERATIONALLY ILL-DEFINED in dwarfs. The boundary is both data-quality (scatter inflated, Phase 89) and physical (correlation structure different, Phase 90), driven by the fundamental transition from disk-dominated spirals to gas-dominated dwarfs (Phase 91).

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

PROGRAM AT PHASE 98. External replication: PARTIAL. Signal replicates qualitatively within N45 (r=0.625, p<0.0001) but calibration differs 45.7% between samples. Train non-N45→test N45: R²=0.23 (weak but positive). Reverse: R²=0.607 (strong). Within-N45 LOO R²=0.327. Attenuation is consistent with restricted-range effect, not law failure. The inside-to-outside coupling is ROBUST; quantification needs distance/resolution attention. Phase chain: 93→94→95→96→97→98.

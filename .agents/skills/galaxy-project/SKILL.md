---
name: galaxy-project
description: Project context, goals, current status, and strategic direction for the Galaxy Rotation Curve Analyzer. Use when planning work, understanding project history, or making decisions about what to do next.
---

# Galaxy Rotation Curve Analyzer — Project Context

## What This Project Is

A data analysis tool studying the Radial Acceleration Relation (RAR) across
197 galaxies (175 SPARC + 22 LITTLE THINGS). Built with React + Vite.

## Single Defensible Claim

> "There is an acceleration transition scale, robust across 197 galaxies
> and 4 interpolation functions, whose physical interpretation has not
> been settled."

This is the ONLY claim that survives all diagnostic tests. Do not overclaim.

## What We Know For Sure

- A transition acceleration scale exists (all interpolation functions find it)
- The RAR is tight (~0.15 dex in GOLD+i45, ~0.24 dex full sample)
- The BTFR holds at 4.4σ
- Signal survives 5 independent stress tests
- ALL 6 residual audit variables clean in GOLD+i45 sample
- Formal uncertainty: a₀ = 5275 ± 0.185 dex [3443, 8082] (km/s)²/kpc

## What Is Model-Dependent

- The precise value of a₀ (GOLD+i45: 5275, GOLD: 4285, HIGH-MASS: 4035, ALL: 1929)
- The cosmological ratio a₀/(cH₀/2π) = 1.641 (GOLD+i45) or 1.333 (GOLD) or 0.600 (ALL)
- a₀ varies across sub-populations (low-mass: fitting artifact — limited g_bar range)
- Dominant uncertainties: global-vs-median (0.112 dex) and inclination cut (0.106 dex)

## What We Do NOT Know

- Is a₀ ≈ cH₀/2π physics or coincidence? Ratio range: 1.07–2.52 (1σ)
- Does a₀ evolve with redshift?
- What is the physical mechanism behind the transition?
- Can error-weighted / hierarchical fitting reduce the ±0.185 dex uncertainty?

## Project History (Condensed)

1. Built analyzer app with SPARC + LITTLE THINGS data
2. Detected RAR, BTFR, transition scale across 197 galaxies
3. Ran 5 stress tests (null, scramble, M/L, distance MC, quality)
4. Explored theoretical interpretations (Parts XXVI-XXVIII)
5. Self-corrected: "floor" was wrong, model is MOND variant (Part XXIX)
6. Ran novelty audit (Part XXX) — honest about what's new vs known
7. Designed mass-matched BTFR forecast pipeline (Part XXXI)
8. Diagnosed hidden factor sources (Part XXXII)
9. **PIVOT**: Created diagnostic report with 3 actual computational tests
10. Diagnostic tests revealed: residual correlations are the #1 open issue
11. **Phase 2 diagnostics**: Found Phase 1 residual bug (gbar*a₀ error).
    Corrected analysis shows gas/SB correlations were ARTIFACTS.
    Real picture: residuals clean for well-constrained galaxies.
    Low-mass a₀ divergence = fitting artifact (limited g_bar range).
12. **Gold-standard pipeline v1.0.0**: Locked parameters, 4 sample layers,
    blind residual audit. GOLD sample a₀=4285, 5/6 audit variables clean.
    Only inclination marginal (r=-0.31). Pipeline fully reproducible.
13. **Gold-standard pipeline v2.0.0**: Added GOLD+i45 sample (inc≥45°),
    6-component uncertainty budget, sensitivity tables.
    GOLD+i45 a₀=5275 ± 0.185 dex. ALL 6 audit variables CLEAN.
    Inclination signal RESOLVED. Literature value within 1σ band.

## Current Phase: DIAGNOSTIC (not theoretical)

Do NOT pursue:
- Quantum gravity / Unruh temperature derivations
- "Why does a₀ = cH₀/2π" explanations
- New theoretical frameworks

DO pursue:
- Error-weighted fitting and nuisance marginalization
- Hierarchical/Bayesian a₀ estimation
- Full-sample refit (all 4123 points)
- Clean, reproducible analysis (gold-standard pipeline v2.0.0 built)

## Two Separate Questions (Never Mix Them)

**A. Descriptive**: What is the best phenomenological law?
**B. Causal**: Why does this law exist?

Solve A first. B becomes much cleaner after.

## Three Target Products

1. **Reproducible analysis package**: data + code + exact assumptions
2. **Systematics paper**: interpolation, nuisance params, baryon calibration, residual structure
3. **Forecast paper**: mass-matched BTFR at z ~ 0.5-1.0, power analysis

## Language Rules

- Say "framework" not "model" or "theory"
- Say "transition scale" not "floor" or "minimum acceleration"
- Say "cosmologically fixed MOND-like interpolation" not "Cosmic Scale model"
- Keep claims conditional: "if confirmed..." / "pending literature check..."
- Never say "we proved" or "we discovered" for known results (RAR, BTFR, a₀)

## File Reference

| File | Purpose |
|------|---------|
| `diagnostic-report.txt` | PRIMARY report (edit this one) |
| `galaxy-rotation-curve-analysis-full-report.txt` | ARCHIVED (read-only) |
| `scripts/diagnostic-tests.cjs` | Main diagnostic pipeline |
| `scripts/transition-scale-v3.cjs` | Core RAR analysis |
| `public/diagnostic-results.json` | Latest computed results |
| `public/rar-analysis-real.json` | Per-galaxy metadata |
| `public/transition-scale.json` | Full RAR + fitting output |
| `replit.md` | Project documentation |

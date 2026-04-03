---
name: galaxy-project
description: Project context, goals, current status, and strategic direction for the Galaxy Rotation Curve Analyzer. Use when planning work, understanding project history, or making decisions about what to do next.
---

# Galaxy Rotation Curve Analyzer — Project Context

## What This Project Is

A data analysis tool studying the Radial Acceleration Relation (RAR) across
197 galaxies (175 SPARC + 22 LITTLE THINGS). Built with React + Vite.

## Single Defensible Claim (v4.0)

> "There is a real acceleration transition scale in galaxy rotation curves.
> After per-galaxy Y* marginalization and hierarchical modeling on the full
> (not subsampled) dataset, our most defensible estimate is
> a₀ = 3633 (km/s)²/kpc = 1.18e-10 m/s², with substantial between-galaxy
> heterogeneity (tau = 0.291 dex, I²=92.4%) and total systematic
> uncertainty ±0.304 dex. The literature value (1.20e-10 m/s²) falls
> within our 1σ band. The cosmological ratio a₀/(cH₀/2π) = 1.13 is
> suggestive but shows moderate variation across subsamples (CV = 14.1%)."

NOT claimed: a₀ universal constant, dark matter solved, MOND proved,
cosmological origin established. Do not overclaim.

## What We Know For Sure

- A transition acceleration scale exists (McGaugh RAR finds it in all samples)
- HEADLINE: a₀ = 3633 (km/s)²/kpc = 1.18e-10 m/s² (v4.0 marg. hier.)
- v3.0 baseline: hier. a₀ = 3374; diagnostic/raw: UW a₀ = 4837
- v4.0 is within 2% of literature (3703 = 1.20e-10 m/s²)
- The RAR is tight (~0.147 dex in GOLD+i45, ~0.200 dex full sample)
- The BTFR holds at 4.4σ
- Signal survives 5 independent stress tests
- ALL 6 residual audit variables clean in GOLD+i45 sample
- Y* marginalization does NOT reduce tau — heterogeneity is real
- Formal uncertainty: ±0.304 dex (dominated by between-galaxy tau = 0.291)

## What Is Model-Dependent

- Precise a₀ depends on method: v4.0 marg=3633, v3.0 hier=3374, v3.0 UW=4837
- Cosmological ratio: 1.130 (v4.0), 1.049 (v3.0 Hier), 1.504 (v3.0 UW)
- Between-galaxy heterogeneity: I²=92.4%, tau=0.291 dex
- Ratio stability: CV=14.1% across 11 subsample splits

## What We Do NOT Know

- Is a₀ ≈ cH₀/2π physics or coincidence? Ratio=1.13 but CV=14.1%
- What drives the between-galaxy scatter (tau=0.291 dex)?
  (Y* marginalization only accounts for +0.046 dex of tau increase)
- Does a₀ evolve with redshift?

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
13. **Gold-standard pipeline v2.0.0**: Added GOLD+i45 sample (inc>=45),
    6-component uncertainty budget, sensitivity tables.
14. **Pipeline v3.0.0**: Single McGaugh-only estimator (dropped mean-of-3),
    error-weighted fit, DerSimonian-Laird hierarchical random-effects model.
    Three estimators converge: UW=4837, W=3545, Hier=3374.
    tau=0.245 dex, I²=89.4%, total budget ±0.295 dex.
15. **Pipeline v4.0.0**: Full sample (3755 pts, not subsampled), per-galaxy
    Y* + distance marginalization (65 grid evals/galaxy), hierarchical DL.
    a₀=3633, tau=0.291, I²=92.4%, total budget ±0.304 dex.
    Within 2% of literature. Y* marg does NOT reduce tau.
    Ratio stability test: 11 splits, CV=14.1%. Report updated to 19 sections.

## Current Phase: DIAGNOSTIC (not theoretical)

Do NOT pursue:
- Quantum gravity / Unruh temperature derivations
- "Why does a₀ = cH₀/2π" explanations
- New theoretical frameworks

DO pursue:
- Understanding tau sources — see Section 20 (Suspects Table) in diagnostic-report.txt
- Three red-flag splits from v4.0: distance (35%), inclination (42%), velocity (36%)
- Non-circular motions and pressure support testing (Suspects 1 & 2)
- Within-galaxy covariance modeling
- Systematics paper write-up
- High-z test with JWST/ALMA kinematic data

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
| `scripts/definitive-v4.cjs` | v4.0 pipeline (CURRENT DEFINITIVE) |
| `scripts/final-measurement-v3.cjs` | v3.0 pipeline (archived baseline) |
| `scripts/gold-standard-pipeline.cjs` | v2.0 pipeline (archived) |
| `scripts/diagnostic-tests.cjs` | Phase 1 diagnostics |
| `scripts/transition-scale-v3.cjs` | Core RAR analysis |
| `public/definitive-v4-results.json` | v4.0 results (current) |
| `public/gold-standard-results.json` | v3.0 results (archived) |
| `public/diagnostic-results.json` | Phase 1 results |
| `public/rar-analysis-real.json` | Per-galaxy metadata (3391 SPARC pts) |
| `public/transition-scale.json` | Full RAR + fitting output |
| `replit.md` | Project documentation |

---
name: galaxy-project
description: Project context, goals, current status, and strategic direction for the Galaxy Rotation Curve Analyzer. Use when planning work, understanding project history, or making decisions about what to do next.
---

# Galaxy Rotation Curve Analyzer — Project Context

## What This Project Is

A data analysis tool studying the Radial Acceleration Relation (RAR) across
197 galaxies (175 SPARC + 22 LITTLE THINGS). Built with React + Vite.

## Single Defensible Claim (Phase 60 Final)

> "a₀ is NOT universal — it shows structured, galaxy-dependent variation
> (τ=0.264 dex). About 43% of the variation can be predicted from 6 observable
> properties (gas mass, RC texture, environment, baryon density, RAR coherence,
> baryon dominance). The variables act additively and independently. The remaining
> ~57% is irreducible with SPARC data alone, representing either true physical
> variation or an unmeasured deep state variable. Resolution requires external data."

NOT claimed: a₀ universal constant, dark matter solved, MOND proved.

## Final Model Hierarchy (Phase 60)

| Model | k | LOO Gap-Closed | AIC rank | Verdict |
|-------|---|----------------|----------|---------|
| M0: Universal a₀ | 0 | 0.0% | worst | REJECTED |
| M1a: Conservative | 4 | 27.8% | 4th | Publication-grade |
| M1b: Extended | 5 | 38.1% | 3rd | Internal analysis |
| M1c: Maximal | 6 | 42.9% | 2nd | WINNER |
| M2: Per-galaxy | 56 | 100.0% | 1st | Overfitting |

## Confirmed Variables

| Variable | b | beta | partial r | Circ | Domain |
|----------|---|------|-----------|------|--------|
| log(MHI) | -0.231 | -0.515 | -0.541 | ZERO | Gas content |
| rcWiggliness | +2.289 | +0.313 | +0.345 | MED | RC texture |
| envCode | -0.121 | -0.359 | -0.385 | ZERO | Environment |
| Sigma0_bar | +0.183 | +0.322 | +0.351 | ZERO | Baryon density |
| meanRun | +0.367 | +0.301 | +0.368 | MED | RAR coherence |
| innBarF | (in maximal) | | | MED | Baryon dominance |

## Project History (Condensed)

1-15. Built analyzer, detected RAR, ran stress tests, explored models,
      established M3 (per-galaxy a₀) wins decisively over M0 (universal)
16-52. Sequential Door Investigation: 6 doors, 37 phases, ~138 proxies
       5 confirmed variables, 95% rejection rate
53. Final comprehensive review — discovery map
54. Freeze baselines + filter partials (tBst/Rd REVOKED)
55. Interaction terms — ALL 10 FAIL (additive model)
56. Frozen baselines — definitive reference point
57. 2D kinematics — ALL 14 proxies FAIL (rcWig+meanRun optimal)
58. Environmental processing — ALL 12 proxies FAIL (envCode irreducible)
59. Stellar Y* — barGravD + innBarF SURVIVE (don't absorb Σ₀)
60. FINAL DEATH MATCH — M1c wins, a₀ structured but incompletely explained

## Current Phase: PROGRAM COMPLETE

All investigation doors closed. Final model comparison done.
a₀ is structured but incompletely explained.

POSSIBLE FUTURE WORK:
- Enlarge sample (BIG-SPARC when available)
- External data: IFU spectroscopy, multi-band imaging, 2D velocity fields
- True SPS Y* from multi-band SED fitting
- Non-linear models / machine learning approaches
- Systematics paper write-up

## Key Files

| File | Purpose |
|------|---------|
| `diagnostic-report.txt` | PRIMARY report (82 sections, 5814 lines) |
| `experiment-log.txt` | Experiment log (93 sections, 2285 lines) |
| `scripts/phase56-frozen-baselines.cjs` | Definitive reference |
| `scripts/phase60-death-match.cjs` | Final model comparison |
| `public/phase56-frozen-baselines.json` | Frozen baseline data |
| `public/phase60-death-match.json` | Death match results |
| `public/sparc-table.json` | SPARC master table (175 galaxies) |
| `public/phase11-sensitivity-lab.json` | Per-galaxy RAR profiles |
| `public/phase25-group-membership.json` | Environment assignments |

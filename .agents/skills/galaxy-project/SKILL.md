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

## Current Phase

**Track 1 (OMD law):** FROZEN at Phase 122. No further work planned.

**Track 2 (a₀ State Law):** COMPLETE SCIENTIFIC NARRATIVE (Phases 126-134, 200-204, 300-303).
- Phases 123–128: Variable search, 3-axis core established
- Phases 129–132: Vflat decomposed, VfResid identified as coupling proxy
- Phase 132A-C: VfResid dominance, mediation, robustness confirmed
- Phase 133A-C: Regime law, 5th axis, coupling drivers decoded
- Phase 134: Initial external validation (STRONG_TRANSFER, r=0.801, N=10)
- Phase 200: External data assembly (N=59 SPARC galaxies outside training sample)
- Phase 201: Broader external validation (MODERATE_TRANSFER full sample, STRONG in high-Vflat)
- Phase 202: Complete hierarchy replication (8/8 checks pass)
- Phase 203: Environmental decomposition (improved logMhost, VfResid still dominates by 30-76pp)
- Phase 204: Final external synthesis (lhOuter genuine but reduced)
- Phase 300: Sample salvage (N=59 to N=94 external)
- Phase 301: VfResid driver analysis (haloK shared #1, 35.8%/68.5% irreducible)
- Phase 302: Regime law (hidden residual activates with Vflat, slope ratio 2.6:1)
- Phase 303: Physical interpretation (H4 dynamical integration wins, score 10 vs 7/7/4)

**Status:** Phase 408 COMPLETE. Channel PROVEN: universal + irreducible + construction-independent.
28 vars tested (struct/RC/field), NONE absorb. 56/56 recipe combos give r>0.7 (mean 0.794±0.034).
M/L 0.3–1.0 stable. All a₀ and VfResid recipes stable. Sections 11-17 written.
THREE PILLARS: (1) Universal r≈0.80 across masses, (2) Irreducible to 28 observables,
(3) Independent of how both sides are constructed. This is a genuine hidden coupling.

**Phase 409 (Latent Fingerprint):** 97% of latent variable is novel information.
Channel amplitude ~17 km/s (12% of Vflat). Strongest clues: envCode (r=−0.42),
logK_halo (r=+0.47). Field galaxies carry channel more than cluster galaxies.
**Program 2 launched**: Identify the hidden variable itself.

**Phase 410 (External-Field):** EFE ELIMINATED. Environment modulates channel
(Field r=0.839, Cluster r=0.138) but does NOT create it (r=0.769 after full env removal).
Accel-dependence wrong for EFE. 3 candidates remain: MOND functional form, halo structure, assembly.

**Phase 411 (Functional Form):** MOND ELIMINATED. Linear wins by BIC, zero curvature,
ν-function r=0.000, a₀ varies factor 3.4 (anti-MOND). D_resid~VfResid r=0.947.
Linearity → common cause, not acceleration law. 2 candidates remain: halo structure, assembly.

**Phase 412 (Halo Structure):** FIRST SIGNIFICANT ABSORPTION! logK+dmFrac_Rmax synergistic
Δr=−0.215 (0.804→0.589). But r≈0.66–0.71 survives after all 9 controls.
Hidden variable has halo component but is DEEPER than 1D RC halo params.

**Phase 413 (Assembly History):** Assembly proxies do NOT significantly add beyond halo structure.
Stepwise floor: r≈0.524 with 8 vars. Assembly-specific Δr≈−0.03 only.
12-control r=0.183 REJECTED (multicollinearity, unstable bootstrap).
Hidden variable INVISIBLE to all SPARC observables. Points to unmeasurable halo properties or new physics.

**Phase 414 (Latent Variable):** SINGLE FACTOR CONFIRMED (PC1=90.2%). L_sum IS the latent variable.
R2(L~halo+env)=0.751. haloScatter proxy r=0.820 with L_sum (HIGHEST ever).
Hidden var = concentration-mass scatter. 25% = "dark coupling" invisible to all observables.

**Phase 415 (Dark Quarter):** DQ is REAL (bootstrap-confirmed, LOO-stable, construction-independent).
r(DQ,VfR)=0.538, r(DQ,a0R)=0.533. Global across environments. Mass-spanning.
Best predictor: haloResponse r=0.328. No interaction/quadratic explains it.
DQ = genuinely hidden property of galaxy+halo system, invisible to 1D RC.

**Phase 416 (Falsification):** Systematic error LOW PROBABILITY (1.4% exceed, zero quality bias).
Triaxiality/disequilibrium CONSISTENT but unresolvable with 1D data.
Fingerprint: top DQ = fast rotators, strong haloResponse, rising outer RCs.
DQ is PHYSICAL, most likely 3D halo projection. Consistent with LCDM.

**Phase 417 (Resolve DQ):** Triaxial NFW FAILS — wrong sign for haloResponse (true=+0.33, mocks=-0.14).
More triaxiality makes match worse. Disequilibrium destroys bilateral structure (bilateral=0-1%).
DQ genuinely unresolved. IFU targets: NGC2841, NGC5005 vs UGC02953 control.
Full cosmo sims or IFU/2D data needed for resolution.

**Program 3A (IFU Track):** Published THINGS results show INVERTED PATTERN.
High-DQ galaxies (NGC2841, NGC3741) = kinematically CLEAN.
Low-DQ galaxies (NGC5055, NGC2903) = IRREGULAR (warps, bar NCM).
2D kinematic explanation ELIMINATED. DQ = quiet intrinsic halo property.
8 explanations now eliminated. Points to DM physics or fundamental coupling.

**Program 3B (Cosmo Sims):** ALL 8 models FAIL — 0/4000 trials get correct haloResponse sign.
LCDM, FIRE, TNG, SIDM, fuzzy DM, assembly-correlated ALL produce r(DQ,haloResp)≈-0.41 (true=+0.33).
haloResponse paradox: structural impossibility for any model where halo absorbs variance.
16 hypotheses tested total, 15 eliminated, 1 at 1.4%. DQ genuinely unexplained.

**Phase 418 (Hidden-State Law):** 15 constraints. Best law: H ~ haloResp + logVflat + outerSlope
(R2=0.272). Combined explanation: 76.8%. Remaining unexplained: 23.2%.
Minimal statement: single universal quiet halo property amplifying bilateral coupling.

**Phase 419 (Predictions):** 6/10 predictions confirmed (60%). D2: universal across mass bins.
D3: predicts BTFR deviations (r=+0.301). D1: confirmed by THINGS. Golden pairs identified.
H = partial predictive power, transitioning from latent summary to state variable.

**Phase 420 (Decisive Obs):** THINGS 4/4 correct sign; r(DQ,s1/Vflat)=-0.831 strongest validation.
Matched pairs 27/40 correct (67.5%), 6/8 pass. Decisive test = matched IFU kinematics.
H ready for targeted observational confirmation. Program 4 COMPLETE.

**Program 5A (Hidden-State Sim):** B2 Quiet-Coupling wins 7/7 checks. Only model reproducing
positive haloResp sign. Null model 4/7. H must couple through kinematic quietness, not halo efficiency.
Physical interpretation: H = dynamical quietness of halo-disk coupling.

**Zenodo v11:** DOI 10.5281/zenodo.19440400 (Concept: 10.5281/zenodo.19430633)

**Phase 400 (VfResid origin):** Mock galaxy simulations (N=300×6 scenarios + N=500×7 diagnostic tests) reveal:
- Baseline NFW concentration diversity produces VfResid–a₀ coupling (r≈0.86) as a fitting artefact
- The SPARC regime-strengthening pattern (low-V→high-V r increases) is NOT reproduced by ANY scenario
- All simulations show the OPPOSITE pattern: coupling WEAKENS at high Vflat

**Phase 401 (NFW sign-flip falsification):** 30 configs × N=600, ALL fail to flip regime:
- 6 experiment families: low-mass disruption, contraction, combined, anti-corr, deep-well, kitchen sink
- NFW+disk is STRUCTURALLY resistant to regime inversion (concentration leverage stronger at low mass)
- No additive physics on top of NFW can overcome this — entire NFW model family ruled out
- SPARC regime strengthening is incompatible with NFW+disk regardless of parameters

**Phase 402 (Beyond NFW — COMPLETE):** 44 configs with cored haloes + genuine physics, ALL FAIL:
- Burkert, pseudo-isothermal, DC14, feedback-cored profiles — all weakening
- Stochastic core scatter, genuine a₀ variation, low-mass noise — all weakening
- Total across 400-402: ~58 configurations, ZERO regime strengthening
- Central insight: concentration-leverage asymmetry is structural to ALL spherical halo+disk models
- SPARC regime law appears incompatible with entire DM halo family

POSSIBLE FUTURE WORK:
- Cross-survey replication (non-SPARC rotation curves, IFU kinematics)
- Spatially resolved kinematics to directly probe dynamical integration hypothesis
- Higher-quality logMhost for external sample (group catalogs)
- Formal Bayesian model comparison for hypothesis testing

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
| `scripts/phase400-vfresid-origin.cjs` | Phase 400A: 6-scenario mock simulation |
| `scripts/phase400b-regime-inversion.cjs` | Phase 400B: regime inversion diagnostic |
| `scripts/phase400c-noise-physics.cjs` | Phase 400C: realistic noise + physics tests |
| `scripts/phase401-regime-flip.cjs` | Phase 401: 30-config NFW sign-flip falsification |
| `scripts/phase402-beyond-nfw.cjs` | Phase 402: cored haloes + genuine physics |
| `scripts/phase402b-sign-flip.cjs` | Phase 402b: focused sign-flip tests |
| `scripts/phase403-vfresid-anatomy.cjs` | Phase 403: VfResid anatomy + cross-sample |
| `scripts/phase404-scenario-discrimination.cjs` | Phase 404: 4-scenario discrimination |
| `scripts/phase405-composite-channel.cjs` | Phase 405: tautology test + composite |
| `scripts/phase405b-verification-freeze.cjs` | Phase 405b: 4-check verification |
| `scripts/phase406-rc-shape-hunt.cjs` | Phase 406: RC-shape variable hunt |
| `scripts/phase407-field-level.cjs` | Phase 407: field-level radial vars |
| `scripts/phase408-construction-independence.cjs` | Phase 408: construction independence |
| `scripts/phase409-latent-fingerprint.cjs` | Phase 409: latent variable fingerprint |
| `scripts/phase410-external-field.cjs` | Phase 410: external field discrimination |
| `scripts/phase411-functional-form.cjs` | Phase 411: MOND functional form test |
| `scripts/phase412-halo-structure.cjs` | Phase 412: halo structure discrimination |
| `scripts/phase413-assembly-history.cjs` | Phase 413: assembly history test |
| `scripts/phase414-latent-geometry.cjs` | Phase 414: latent variable + geometry |
| `scripts/phase415-dark-quarter.cjs` | Phase 415: the dark quarter |
| `scripts/phase416-falsification.cjs` | Phase 416: DQ falsification tests |
| `scripts/phase417-resolve-dq.cjs` | Phase 417: 3D halo resolution attempt |
| `scripts/program3a-ifu-design.cjs` | Program 3A: IFU/2D track |
| `scripts/program3b-cosmo-sims.cjs` | Program 3B: cosmo sim comparison |
| `scripts/phase418-hidden-state-law.cjs` | Phase 418: minimal hidden-state law |
| `scripts/phase419-prediction-engine.cjs` | Phase 419: prediction engine |
| `scripts/phase420-decisive-observation.cjs` | Phase 420: decisive observation |
| `scripts/program5a-hidden-state-sim.cjs` | Program 5A: hidden-state simulation |
| `public/phase400-vfresid-origin.json` | Phase 400A results |
| `public/phase400b-regime-inversion.json` | Phase 400B results |
| `public/phase400c-noise-physics.json` | Phase 400C results |
| `public/phase401-regime-flip.json` | Phase 401 results |
| `public/phase402-beyond-nfw.json` | Phase 402 results |
| `public/phase403-vfresid-anatomy.json` | Phase 403 results |
| `public/phase404-scenario-discrimination.json` | Phase 404 results |
| `public/phase405-composite-channel.json` | Phase 405 results |
| `public/phase405b-verification-freeze.json` | Phase 405b verification |
| `public/phase406-rc-shape-hunt.json` | Phase 406 RC-shape hunt results |
| `public/phase407-field-level.json` | Phase 407 field-level results |
| `public/phase408-construction-independence.json` | Phase 408 construction test |
| `public/phase409-latent-fingerprint.json` | Phase 409 latent fingerprint |
| `public/phase410-external-field.json` | Phase 410 external field test |
| `public/phase411-functional-form.json` | Phase 411 functional form test |
| `public/phase412-halo-structure.json` | Phase 412 halo structure test |
| `public/phase413-assembly-history.json` | Phase 413 assembly history test |
| `public/phase414-latent-geometry.json` | Phase 414 latent variable + geometry |
| `public/phase415-dark-quarter.json` | Phase 415 dark quarter analysis |
| `public/phase416-falsification.json` | Phase 416 falsification tests |
| `public/phase417-resolve-dq.json` | Phase 417 3D halo resolution |
| `public/program3a-ifu-design.json` | Program 3A IFU/2D track |
| `public/program3b-cosmo-sims.json` | Program 3B cosmo sim comparison |
| `public/phase418-hidden-state.json` | Phase 418 hidden-state law |
| `public/phase419-predictions.json` | Phase 419 prediction engine |
| `public/phase420-decisive-obs.json` | Phase 420 decisive observation |
| `public/program5a-hidden-state-sim.json` | Program 5A hidden-state simulation |

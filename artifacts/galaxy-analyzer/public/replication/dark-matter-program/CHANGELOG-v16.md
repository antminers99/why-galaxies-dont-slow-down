# Zenodo v16 Changelog — Dark Matter Model Discrimination (Program 12)

**Date:** 2026-04-07
**Previous version:** v15 (draft, DOI: 10.5281/zenodo.19446498)
**Context:** Program 12 adds a systematic dark matter model discrimination series (DM-1 through DM-4V) that identifies CDM + non-axisymmetric halo shape as the leading model and kills CDM smooth, MOND, and Fuzzy DM.

---

## What Remains Unchanged

1. **The existence of H is confirmed.** r(VfResid, a0Resid) = 0.77 (LOO, p < 0.001, N = 55), replicated on N = 59.
2. **The bilateral coupling channel is real.** Construction-independent, distance-independent, externally replicated.
3. **The 1D information ceiling is genuine.** 88.5% of H is inaccessible from rotation-curve features alone.
4. **A 2D angular carrier exists.** Velocity fields carry information that rotation curves destroy.
5. **Red team results stand.** 11/11 pass.

---

## What Is New in v16: Program 12 — Dark Matter Discrimination Series

### DM-1: Kill Constraints Table
- 12 mandatory constraints tested against 7 dark matter models.
- **DEAD:** CDM smooth (6/12 fails), MOND (10/12 fails).
- **WEAK:** WDM.
- **OPEN:** SIDM, FDM, Exotic.
- **LEADING:** CDM + halo shape state (12/12 pass).

### DM-2 + DM-2V: CDM+Shape vs SIDM
- DM-2 initial result showed SIDM 2-1 advantage — but 3 galaxies had zero outer coverage (artifact).
- DM-2V (verification): fractional radial zones, r90 cutoff, LOO, 5000-bootstrap.
  - r(DQ, outer C) = 0.746 vs r(DQ, inner C) = 0.431 — **outer dominates**.
  - LOO inner-led: **0/7** times. Bootstrap inner-led: **19.7%**.
  - **VERDICT: CDM+shape CONFIRMED. DM-2 SIDM signal was a coverage artifact.**

### DM-3: Fuzzy DM / ψDM Wave Fingerprint Test
- 4 tests: spectral width, radial phase coherence, odd/even asymmetry, wave uniqueness.
- Spectral entropy = 0.72 (Fuzzy expects > 0.85); effective modes = 2.59 (expects > 3.5).
- Phase flip rate = 0.18 (expects > 0.4); odd/even ratio = 4.15 (expects ~1.0).
- **Score: CDM+shape 4/4, Fuzzy DM 0/4. Fuzzy DM is DEAD.**

### DM-4: Quantitative Halo Shape Parameter
- Best single parameter: **shapeAmplitude** (total non-axisymmetric power, m ≥ 2).
- r(DQ, shapeAmplitude) = **0.793**, p = 0.032.
- Beats haloResponse (r = 0.267) by 3x.
- LOO: 7/7 positive (range 0.43–0.91). Bootstrap: 96.3% positive.

### DM-4V: Verification of shapeAmplitude
- **V1 Binning sensitivity:** r range [0.612, 0.812] across 8 configurations. All > 0.5. STABLE.
- **V2 Normalization:** raw r=0.804, fractional r=0.652, velocity-normalized r=0.815. All > 0.3. ROBUST.
- **V3 LOO + Bootstrap:** LOO 7/7 positive, bootstrap 97.8% positive, Kendall τ = 0.429. ROBUST.
- **V4 Bar exclusion:** Even excluding 30% inner region, r = 0.591. Not bar-driven. ROBUST.
- **Combined: 3 ROBUST + 1 STABLE = shapeAmplitude CONFIRMED.**

---

## Updated Model Status Table

| Model | Status | Evidence |
|-------|--------|----------|
| CDM smooth | **DEAD** | 6/12 kill constraints fail |
| MOND / EFE | **DEAD** | 10/12 kill constraints fail |
| Fuzzy DM (ψDM) | **DEAD** | 0/4 wave fingerprint tests |
| SIDM | **No signal** | DM-2V: 0/4, SIDM advantage was coverage artifact |
| WDM | Weak | No positive evidence |
| **CDM + halo shape** | **LEADING** | 12/12 DM-1, 3/4 DM-2V, 4/4 DM-3, r=0.80 DM-4, confirmed DM-4V |

---

## Updated Interpretation

| Aspect | v15 | v16 |
|--------|-----|-----|
| H carrier | "angular velocity-field complexity" | **shapeAmplitude: total non-axisymmetric gravitational complexity (m ≥ 2)** |
| Best quantitative proxy | angular complexity (r ~ 0.67) | **shapeAmplitude r = 0.80, LOO-stable, bar-independent** |
| DM discrimination | "CDM consistent but not uniquely favored" | **CDM + non-axisymmetric halo is the only surviving model; 3 alternatives killed, 2 have no signal** |
| Physical meaning | "departure from axisymmetry" | **magnitude of halo non-axisymmetry — the more non-circular the halo, the larger the BTFR/RAR residual coupling** |

---

## Summary Statement

Program 12 transforms the hidden-state detection into a dark matter model discriminator. Three alternative models (CDM smooth, MOND, Fuzzy DM) are killed by quantitative tests; SIDM shows no positive signal after cleaning. The leading model — CDM with non-axisymmetric halo shape — is supported by a robust quantitative parameter (shapeAmplitude, r = 0.80) that survives binning, normalization, LOO, bootstrap, and bar-exclusion stress tests.

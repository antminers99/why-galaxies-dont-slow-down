# Decisive IFU Test of the Hidden State H

**A Follow-up Proposal to: "Per-Galaxy a0 Variation and the Hidden Common-Cause State in SPARC Galaxies"**

---

## 1. Scientific Motivation

Our analysis of 55 SPARC galaxies reveals a hidden physical variable H that drives a bilateral coupling between kinematic excess (VfResid) and acceleration discrepancy (a0Resid) at r = 0.77 (p < 0.001). This coupling survives all statistical challenges: leave-one-out cross-validation, construction independence (55/55), permutation testing, and a dedicated distance-error kill test (V+).

However, ~70-80% of H is structurally inaccessible from 1D rotation curves. This is not a methodological limitation but a fundamental information barrier: azimuthal averaging destroys the angular structure where H resides. No rotation-curve feature -- from simple scalars to 71-dimensional state vectors -- recovers more of H than the scalar haloResponse alone.

**The decisive test requires 2D kinematic data.**

## 2. The Inaccessibility-Strength Paradox

Program 8B tested 12 generative models across three physics families (halo redistribution, formation-history common cause, exotic DM). Key finding:

- Models producing strong coupling (r >= 0.65) make H fully predictable from RC features (C6 FAIL)
- Models achieving high inaccessibility (>50% hidden) produce channels too weak to match the data (C1 FAIL)
- No single-layer model achieves both simultaneously

The empirical H achieves both. This implies H operates through channels that affect integral quantities (Vflat, a0) without proportionally affecting the radial profile shape -- a signature of multi-scale physics requiring angular resolution.

## 3. Proposed IFU Test

### 3.1 Target Selection

From our N=55 sample, select matched pairs with:
- Similar structural properties (logMbar, logRdisk, morphT within 0.1 dex / 1 unit)
- Divergent bilateral excess (DQ difference > 1.5 sigma)

Priority targets (existing pairs from our matched-pair analysis):

| High-DQ Galaxy | Low-DQ Galaxy | Structural Match | DQ Difference |
|---------------|---------------|-----------------|---------------|
| NGC 2841 | UGC 02953 | Excellent | 2.95 sigma |
| NGC 3198 | NGC 3521 | Good | ~2 sigma |
| NGC 7331 | NGC 5055 | Good | ~1.5 sigma |

### 3.2 Observable Predictions

If H is real and angular:

1. **Velocity field asymmetry**: High-DQ galaxies should show systematically different velocity field residuals (after tilted-ring subtraction) than their matched low-DQ counterparts
2. **Non-circular motion patterns**: High-H galaxies should exhibit coherent radial or vertical motions not present in low-H galaxies
3. **Angular halo structure**: The ratio V_obs(R,theta)/V_circ(R) should show theta-dependent structure correlated with DQ

If H is NOT angular but rather a global scaling:

4. **No velocity field differences**: Matched pairs should show identical velocity field structure despite different DQ
5. **Global potential offset**: The difference should appear only in the total enclosed mass, not in its angular distribution

### 3.3 Required Data

- **Instrument**: MUSE/VLT, SITELLE/CFHT, or MANGA/SDSS (for galaxies with existing coverage)
- **Spatial resolution**: < 1 kpc at galaxy distance (for radial profile distinction)
- **Spectral resolution**: R > 5000 (for velocity precision < 10 km/s)
- **Coverage**: Full disk to > 2 Rdisk (to capture outer halo region where H effects are strongest)
- **Sample size**: Minimum 6 matched pairs (12 galaxies); ideal 15 pairs (30 galaxies)

### 3.4 Analysis Pipeline

1. **Tilted-ring decomposition**: Extract V_circ(R), V_rad(R,theta), V_vert(R,theta)
2. **Residual velocity field**: V_obs(R,theta) - V_model(R,theta) for each galaxy
3. **Angular power spectrum**: Decompose residuals into Fourier modes m=0,1,2,...
4. **Matched-pair comparison**: Test whether m >= 1 power differs between high-DQ and low-DQ galaxies
5. **H recovery**: Attempt to reconstruct H from 2D features and test whether R-squared exceeds the 1D ceiling (~30%)

## 4. Expected Outcomes

### Scenario A: IFU breaks the ceiling (most informative)

- 2D features recover > 50% of H variance (vs ~30% from 1D)
- Angular Fourier modes m=1,2 correlate with DQ
- H is identified as angular halo asymmetry or non-circular flow pattern
- **Impact**: Identifies the physical carrier of H; enables direct measurement in IFU surveys

### Scenario B: IFU confirms the ceiling (still valuable)

- 2D features recover < 40% of H variance (no improvement over 1D)
- Velocity field residuals are identical between matched pairs
- **Impact**: H is not angular but operates at scales or in quantities not captured by kinematics (e.g., assembly history, DM micro-physics); redirects search to lensing or SFH proxies

### Scenario C: Mixed signal

- Some angular power in matched pairs but insufficient to fully recover H
- **Impact**: H has both angular and non-angular components; partial identification possible

## 5. Existing Survey Coverage

Several SPARC galaxies already have IFU or high-resolution HI data:

- **THINGS** (The HI Nearby Galaxy Survey): 2D HI velocity fields for several SPARC galaxies including NGC 2841, NGC 3198, NGC 7331
- **PHANGS-MUSE**: Optical IFU for NGC 4826, NGC 5055
- **MANGA**: May include some SPARC galaxies; cross-match needed
- **CALIFA**: Several SPARC targets in the catalog

**Immediate first step**: Cross-match our N=55 sample against THINGS, PHANGS, and MANGA to identify galaxies with existing 2D kinematic data. If >= 3 matched pairs exist in archival data, the test can begin immediately without new observations.

## 6. Summary

The SPARC 1D analysis has reached a provable information ceiling. The hidden state H is real (r = 0.77, p < 0.001, distance-independent) but ~70-80% inaccessible from azimuthally-averaged rotation curves. IFU velocity fields provide the most direct path to breaking this ceiling. The test is well-posed, the predictions are falsifiable, and archival data may already contain sufficient coverage for a pilot study.

---

*This memo accompanies the program closure of "Per-Galaxy a0 Variation and the Hidden Common-Cause State in SPARC Galaxies" (Zenodo Concept DOI: 10.5281/zenodo.19430633).*

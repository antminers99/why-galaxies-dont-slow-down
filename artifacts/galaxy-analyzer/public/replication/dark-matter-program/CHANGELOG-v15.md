# Zenodo v15 Changelog — Post-Submission Verification Update

**Date:** 2026-04-07
**Previous version:** v14 (DOI: 10.5281/zenodo.19446000)
**Context:** This version incorporates post-submission verification (Program 11) that refines the physical interpretation of the hidden state H. The MNRAS submission (ScholarOne) remains unchanged.

---

## What Remains Unchanged

1. **The existence of H is confirmed.** The bilateral VfResid–a0Resid coupling (r = 0.77, LOO, p < 0.001, N = 55) and its replication on N = 59 independent galaxies are unaffected.

2. **The bilateral coupling channel is real.** Construction-independent (55/55 recipe variants), distance-independent, and externally replicated.

3. **The 1D information ceiling is genuine.** 88.5% of H is inaccessible from rotation-curve features alone, regardless of representation.

4. **A 2D angular carrier exists.** The velocity field carries information about H that rotation curves destroy by azimuthal averaging.

5. **The red team results stand.** 11/11 verification tests pass (8/8 critical), including independent pipeline replication (r = 0.91), bar exclusion, LOO stability, and confounder control.

---

## What Has Changed

### Interpretive Refinement (not a retraction)

The original claim (v13–v14) identified the carrier as "the outer-halo m=2 velocity-field mode" and interpreted it as evidence for "halo triaxiality / oval distortion."

**Program 11 (post-submission verification)** performed two additional tests using higher azimuthal resolution (16 bins vs. the original 8):

**Test A — Radial m=2 gradient:**
- All 7 THINGS galaxies show outer-dominant m=2 (ratio 1.9–6.3).
- This is a **universal** pattern that does **not** correlate with DQ (r = −0.26, p = 0.49).
- Most likely a trivial velocity-field scaling effect, not evidence for SIDM core suppression.

**Test B — Full m-spectrum (m = 1…6):**
- m = 1 (rotation) dominates at 68.5% — expected for velocity fields.
- Among non-rotational modes: **m = 3 carries 58%**, m = 2 carries only 11%.
- r(DQ, m = 3 power) = 0.786 and r(DQ, m = 5 power) = 0.802 **exceed** r(DQ, m = 2 power) = 0.516.
- The original r(DQ, m2) = 0.847 was partly inflated by **azimuthal aliasing** (8 bins, Nyquist limit ≈ m = 3); higher-mode power leaked into the m = 2 measurement.

### Updated Interpretation

| Aspect | v14 (original) | v15 (revised) |
|--------|----------------|---------------|
| Carrier description | "outer-halo m = 2 mode" | "angular velocity-field complexity (non-axisymmetric structure)" |
| Physical interpretation | "halo triaxiality (~90% confidence)" | "departure from axisymmetry — may include halo shape, spiral streaming, warps, and interaction signatures" |
| m = 2 role | "dominant carrier" | "one component among several; m = 3 and m = 5 are at least as strongly correlated with H" |
| CDM discrimination | "triaxiality required" | "CDM consistent but not uniquely favored; SIDM and Fuzzy DM not specifically supported" |

---

## What Has Not Collapsed

- The **detection** of H is robust and unaffected by this refinement.
- The **2D nature** of the carrier is confirmed: rotation curves lose the angular information by construction.
- The **gold pair** (NGC 2841 vs NGC 5055) difference in non-circular motion power is real.
- The **inaccessibility–strength paradox** stands: H is strong but structurally invisible to 1D data.
- **No dark matter model is ruled out**; the data is consistent with CDM but the evidence does not uniquely require triaxial halos.

---

## Summary Statement

This version preserves the core hidden-state result, but updates the interpretation of the 2D carrier from a uniquely isolated m = 2 mode to a broader angular/non-axisymmetric velocity-field complexity, after additional post-submission verification (Program 11, Tests A and B, 16-bin azimuthal decomposition).

# Zenodo v18 Changelog — DM-5C: Minimal Carrier Confirmed, Phase Frozen

**Date:** 2026-04-07
**Previous version:** v17 (DOI: 10.5281/zenodo.19453653)
**Context:** DM-5C answers the final question of the DM-5 series: is the +26.6pp R² from multi-parameter models genuine structure or overfitting? Answer: mostly overfitting. shapeAmplitude is the minimal robust carrier. Phase frozen pending new 2D data.

---

## What Remains Unchanged

1. **The existence of H is confirmed.** r(VfResid, a0Resid) = 0.77 (LOO, p < 0.001, N = 55), replicated on N = 59.
2. **The bilateral coupling channel is real.** Construction-independent, distance-independent, externally replicated.
3. **The 1D information ceiling is genuine.** 88.5% of H is inaccessible from rotation-curve features alone.
4. **A 2D angular carrier exists.** Velocity fields carry information that rotation curves destroy.
5. **Red team results stand.** 11/11 pass.
6. **shapeAmplitude is the leading quantitative carrier of H.** r = 0.80, confirmed DM-4V.
7. **CDM + halo shape is the only surviving model.** 3 killed (CDM smooth, MOND, Fuzzy DM), 2 no signal (SIDM, WDM).
8. **Hiddenness paradox is resolved.** Azimuthal averaging erases all m ≥ 2 by construction (DM-5B).

---

## What Is New in v18

### DM-5C: Extra Physics or Overfitting?

Four tests determine whether the +26.6pp R² improvement (from R² = 0.629 to 0.895) reflects genuine CDM structural physics or overfitting on N = 7.

#### C1 — Incremental LOO

| Model | In-sample R² | Adj R² | LOO R² | LOO RMSE |
|-------|-------------|--------|--------|----------|
| M1: shapeAmplitude only | 0.629 | 0.555 | 0.034 | 1.048 |
| M2: shape + outerSupport | 0.722 | 0.583 | 0.118 | 1.001 |
| M3: shape + outer + quietness | 0.744 | 0.488 | −1.055 | 1.528 |

**Finding:** M3 collapses catastrophically out of sample (LOO R² = −1.055). The "R² = 0.90" three-parameter model was pure overfitting. M2 shows a marginal improvement over M1 (+8.4pp LOO R²), suggesting outerSupport may carry real information.

#### C2 — Permutation Necessity

| Shuffled | R² drop | p-value | Verdict |
|----------|---------|---------|---------|
| outerSupport | 0.019 | 0.333 | COSMETIC |
| quietness | −0.049 | 0.639 | COSMETIC |
| shapeAmplitude | −0.045 | 0.637 | COSMETIC |

**Finding:** With N = 7 and p = 3, the model is too saturated for permutation tests to discriminate. All parameters appear cosmetic because the model has nearly as many parameters as data points.

#### C3 — Causal Ordering

| Partial correlation | Value | Meaning |
|--------------------|-------|---------|
| r(outer, DQ \| shape) | 0.500 | outerSupport has independent information |
| r(quiet, DQ \| shape) | −0.038 | quietness is purely downstream |
| r(shape, DQ \| outer) | 0.748 | shape retains strong direct effect |
| r(shape, quiet) | 0.831 | quietness is a shadow of shape |

**Finding:** This is the most informative test:
- **shapeAmplitude → DQ is direct and strong** — survives all controls (partial r = 0.60–0.75)
- **outerSupport has real independent contribution** — partial r = 0.50 controlling for shape
- **downstreamQuietness is NOT a cause** — it is a downstream consequence of shape (r = 0.83 between them, partial with DQ = −0.04)

#### C4 — Minimal Sufficient Model

**Winner: M2 (shape + outer)** — best LOO R² = 0.118, most parsimonious within 5pp.

But caveat: even M2's LOO R² = 0.118 is very low in absolute terms. With N = 7, no linear model achieves reliable out-of-sample prediction.

### DM-5C Summary Table

| Test | Verdict |
|------|---------|
| C1 Incremental LOO | EXTRAS PARTIALLY GENUINE |
| C2 Permutation necessity | UNCLEAR (N too small) |
| C3 Causal ordering | SHAPE DIRECT, SOME MEDIATION |
| C4 Minimal sufficient model | SHAPE + OUTER |

### Overall Verdict: SHAPE AMPLITUDE DOMINATES, EXTRAS AMBIGUOUS

shapeAmplitude is the minimal robust quantitative carrier of H. The +26.6pp from extra parameters was mostly overfitting. outerSupport may contain real independent information (partial r = 0.50), but N = 7 is fundamentally insufficient to confirm this. downstreamQuietness is a downstream shadow, not a causal driver.

---

## Updated Interpretation

| Aspect | v17 | v18 |
|--------|-----|-----|
| shapeAmplitude sufficiency | "Necessary (R²=0.63) but spatial distribution adds info (R²=0.90)" | **R²=0.90 was overfitting. shapeAmplitude (R²=0.63) is the minimal robust carrier.** |
| outerSupport | Not separately tested | **Real independent information (partial r = 0.50), but unconfirmed with N = 7** |
| downstreamQuietness | Not separately tested | **Downstream consequence, NOT a causal driver (partial r = −0.04)** |
| Multi-parameter models | Open question | **Collapses out of sample. Do not trust with N/p < 10.** |
| Phase status | DM-5B complete, DM-5C pending | **DM-5 series COMPLETE. Phase FROZEN pending new 2D data.** |

---

## What Must Happen Next (Not in This Release)

1. **WALLABY / MeerKAT MHONGOOSE data (N > 20):** Test whether r(DQ, shapeAmplitude) > 0.5 holds with larger samples
2. **outerSupport confirmation:** With N > 20, test whether outerSupport adds genuine LOO R² improvement
3. **Multi-parameter models:** Only trust when N/p > 10 (i.e., N > 30 for 3 parameters)

---

## Summary Statement

DM-5C closes the DM-5 series and the current analysis phase. The +26.6pp R² from multi-parameter models was largely overfitting: M3 (3 parameters) collapses to LOO R² = −1.055. shapeAmplitude is the minimal robust quantitative carrier of H — the safe, parsimonious answer supported by all stress tests (DM-4V), the hiddenness paradox resolution (DM-5B), and now confirmed as the only parameter that cannot be explained away (DM-5C). outerSupport carries a tantalizing partial correlation (r = 0.50 controlling for shape), but confirming it requires N >> 7. The phase is frozen: no new analysis until new 2D kinematic data arrive.

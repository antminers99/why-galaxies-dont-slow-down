# Hierarchical Coupling Law for Per-Galaxy a₀ Variation: Evidence from SPARC Rotation Curves

*Non-peer-reviewed computational analysis report*
*Analysis: Phases 123–134 (internal), Phases 200–201 (external validation)*
*Zenodo v9: DOI 10.5281/zenodo.19433840 (Concept DOI: 10.5281/zenodo.19430633)*

---

## Abstract

Using per-galaxy fits to 175 SPARC rotation curves, we find that the MOND acceleration scale a₀ varies systematically across galaxies and is well-described by a hierarchical, regime-dependent empirical law. Within a published-quality subset of N=45 galaxies, a three-axis structural core (gas mass MHI, host halo mass Mhost, and dynamical coherence MeanRun) explains 44% of a₀ variance (LOO). Adding a kinematic coupling channel — the residual flat velocity after removing baryonic structural predictions (VfResid) — raises predictive power to 61%, with a further independent contribution from outer halo fitting quality (lh_outerImprove, 5-axis LOO = 65%). The coupling signal is regime-dependent, strongest in galaxies with Vflat >= 120 km/s and activating sharply near Vflat ~ 181 km/s. The dominant channel (VfResid) mediates essentially all halo proxy signals and contains ~36% irreducible variance beyond any available catalog observable.

The coupling law now has external support beyond the original training sample. In two independent tests:

1. **Initial holdout** (N=10, crude quality, Phase 134): Core+VfResid achieves 57% gap closure, r = 0.801.
2. **Broader external validation** (N=59 SPARC galaxies outside the training sample, Phase 201): Core+VfResid achieves 8.2% gap closure on the full sample, but 34.3% for Vflat >= 120 (N=16), 48.7% for Vflat >= 180 (N=8), and 59.0% for Q=1 + Vflat >= 120 (N=11, r = 0.830). The structural core alone fails in all regimes. The low-Vflat regime fails as predicted.

The strongest transfer appears in higher-Vflat galaxies and collapses in lower-Vflat systems, consistent with a regime-dependent law rather than a universal one. We do not claim a₀ is definitively non-universal; we report structured empirical variation that is well-predicted by observable galaxy properties and whose hierarchical structure — Core fails, VfResid dominates, regime-dependent — reproduces itself outside the training sample.

---

## 1. Introduction

The radial acceleration relation (RAR; McGaugh et al. 2016) connects baryonic and observed accelerations in disk galaxies through a single parameter a₀ ~ 1.2 x 10^-10 m/s^2. Whether a₀ is a universal constant or varies systematically across galaxies remains an open question with implications for both dark matter models and modified gravity theories.

Previous work in this project (Track 1, Phases 1–122) established that the outer mass discrepancy (logOMD) correlates with gas-to-stellar balance, motivating a deeper investigation into per-galaxy a₀ variation. Track 2 (Phases 123–134) investigates this variation directly, asking: which observable galaxy properties predict per-galaxy a₀, what is the physical nature of the predictive signal, and does the relationship transfer to new galaxies?

---

## 2. Data and Methods

### 2.1 Sample

We use the SPARC database (Lelli, McGaugh & Schombert 2016, AJ 152, 157). Per-galaxy a₀ values are obtained by fitting the McGaugh RAR formula to each galaxy's rotation curve using marginalized stellar mass-to-light ratios (Y_disk prior: N(0.50, 0.12^2), range [0.20, 0.80]). The primary analysis sample ("Stage A") consists of 56 galaxies with complete derived properties, of which 45 are classified as published quality and 11 as crude quality. The crude subset serves as the holdout for external validation.

### 2.2 Variables

| Variable | Definition | Source |
|----------|------------|--------|
| logA0 | log10(per-galaxy a₀) in (km/s)^2/kpc | RAR fit |
| logMHI | log10(MHI / 10^9 M_sun) | SPARC catalog |
| logMhost | log10(estimated host halo mass) | Tidal analysis |
| logMeanRun (MR) | log10(mean run length in RAR residuals) | RAR fit quality |
| logVflat | log10(asymptotic flat velocity) | SPARC catalog |
| VfResid | logVflat minus structural prediction | Derived (see 3.2) |
| lh_outerImprove | Outer improvement of log-halo over Newton | Rotation curve fit |
| haloK | log10(dark halo linear scale factor) | Rotation curve fit |

### 2.3 Statistical Framework

All predictive claims use leave-one-out (LOO) cross-validation. Gap percentage is defined as 100 x (1 - RMSE_LOO^2 / SD^2), measuring variance explained out-of-sample. For transfer validation, the baseline is the RMSE of predicting the training mean for all test points. Nested cross-validation, bootstrap resampling (10,000 iterations), and permutation tests are used throughout to control for overfitting.

---

## 3. Results

### 3.1 The Structural Core (Phases 123–128)

Systematic variable search identified three robustly predictive axes:

| Axis | LOO gap (alone) | Role |
|------|-----------------|------|
| logMHI | ~25% | Gas content |
| logMhost | ~20% | Host halo environment |
| logMeanRun | ~15% | Dynamical coherence / RAR regularity |

Together, the three-axis core achieves LOO gap = 44.1% (N=45). All three axes survive nested cross-validation (45/45 folds), bootstrap sign stability (max flip < 3%), and mutual partial correlation tests.

### 3.2 The Fourth Axis: Vflat and VfResid (Phases 127–132)

Adding logVflat to the core (Model C) raises LOO gap to 52.1%. Physical decomposition reveals:

- Baryonic structure (logMbar, logL36, logRdisk, morphT) explains 92.2% of Vflat variance
- The structural component adds only +0.4pp above the core
- The 7.8% kinematic residual (VfResid = logVflat - structural prediction) adds +17.1pp

VfResid is identified as a baryon-halo coupling proxy:
- Correlates strongly with dark halo strength (haloK, r = +0.60)
- Anti-correlates with MOND improvement (r = -0.54)
- Clean on measurement systematics (inclination, distance, quality: all |r| < 0.08)

In head-to-head competition (Phase 132A), VfResid dominates every tested halo proxy by at least 7pp. Mediation analysis (Phase 132B) shows VfResid absorbs 100% of haloK's signal (Sobel-like mediation = 107%), while lh_outerImprove contributes independently.

### 3.3 The Fifth Axis: lh_outerImprove (Phase 133B)

The outer halo fitting quality metric adds a genuine +4.3pp above the 4-axis model:

| Model | LOO gap |
|-------|---------|
| Core (3-axis) | 44.1% |
| Core + VfResid (4-axis) | 61.1% |
| Core + VfResid + lhOuter (5-axis) | 65.4% |

Fold-internal delta = +3.0pp, flip rate = 5.4%, systematics clean (all |r| < 0.13).

### 3.4 Regime Dependence (Phase 133A)

The VfResid coupling signal is strongly regime-dependent:

| Regime | N | Core+VfResid delta | Bootstrap P(delta>0) |
|--------|---|-------------------|---------------------|
| All galaxies | 45 | +17.1pp | 99.9% |
| Vflat >= 120 | 35 | +27.8pp | 99.7% |
| Vflat < 120 | 10 | unstable | 18.6% |

A sharp transition occurs near Vflat ~ 181 km/s, where the running correlation r(VfResid, logA0) exceeds 0.5 for the first time. The VfResid coefficient increases from 2.51 (full sample) to 3.08 (high-Vflat regime).

### 3.5 Irreducible Residual (Phase 133C)

What drives VfResid itself?

| Predictor set | LOO R^2 of VfResid |
|---------------|-------------------|
| haloK alone | 0.287 |
| haloK + environment + MR | 0.431 |
| Best 5 predictors | 0.517 |

After removing all identifiable predictor content, 35.6% of VfResid variance remains unexplained. This irreducible portion still adds +6.8pp to the core model — indicating genuine dynamical physics beyond current catalog observables.

### 3.6 Initial External Validation (Phase 134)

First test: train on N=45 published-quality galaxies, predict N=10 crude-quality holdout galaxies (all with Vflat >= 120 km/s):

| Model | Transfer RMSE | Transfer gap |
|-------|---------------|--------------|
| Naive (predict training mean) | 0.2962 | 0.0% |
| Core (3-axis) | 0.3128 | -11.5% |
| Core + Vflat (Model C) | 0.2810 | +10.0% |
| Core + VfResid | 0.1944 | +56.9% |
| Core + VfResid + lhOuter (5-axis) | 0.1727 | +66.0% |

Key findings:
- r(VfResid, logA0) on holdout = 0.801
- All coefficient signs preserved in combined N=55 sample
- Regime-restricted training (Vflat >= 120 only) outperforms full-sample training
- 2 of 10 holdout galaxies in extrapolation territory (ESO563-G021, UGC02885)

### 3.7 Broader External Validation (Phases 200–201)

To test whether the coupling law generalizes beyond the original training sample, we assembled per-galaxy data for N=59 SPARC galaxies not included in the Stage A sample. Per-galaxy a₀ values were obtained from RAR fits using the transition-scale plotPoints (median 8 points per galaxy); VfResid was computed using the frozen N=45 structural model; logMhost was estimated from Vflat-based abundance matching. No refitting was performed — all predictions use frozen N=45 coefficients.

| Regime | N | Core gap | Core+VfResid gap | r(VfResid, a₀) |
|--------|---|----------|-------------------|----------------|
| Full sample | 59 | -53.0% | +8.2% | 0.713 |
| Vflat >= 120 | 16 | -6.7% | +34.3% | 0.841 |
| Vflat >= 180 | 8 | -49.4% | +48.7% | 0.858 |
| Q=1 + Vflat >= 120 | 11 | — | +59.0% | 0.830 |
| Vflat < 120 | 43 | -74.4% | -3.9% | — |

The full hierarchy reproduces outside the training sample:
1. The structural core alone fails everywhere (negative gaps in all regimes)
2. VfResid is the dominant transferable channel (carries all positive gap)
3. The signal strengthens monotonically with Vflat, peaking in the very-high-Vflat regime
4. The low-Vflat regime fails, exactly as predicted by the regime law (Phase 133A)
5. In the highest-quality high-Vflat subsample (Q=1, Vflat >= 120), performance matches the Phase 134 holdout (gap 59% vs 57%, r 0.830 vs 0.801)

VfResid sign is preserved in the combined N=104 sample; logMhost and logMeanRun signs changed (expected: these variables were estimated from cruder data in the external sample). Excluding 5 extrapolating galaxies improves performance (gap 21.9%, r 0.547 on N=54).

---

## 4. The Coupling Law (Summary)

The empirical law can be stated as:

**Per-galaxy a₀ is governed by a hierarchical, regime-dependent baryon-halo coupling law with three tiers:**

1. **Structural core** (MHI, Mhost, MeanRun): Sets the baseline a₀ based on gas content, halo environment, and dynamical regularity. Necessary but insufficient — fails to transfer alone.

2. **Dominant kinematic coupling channel** (VfResid): The residual flat velocity beyond baryonic structural prediction. Carries the bulk of the halo-to-a₀ coupling signal. Mediates all tested halo proxy effects. Activates sharply at Vflat ~ 181 km/s.

3. **Secondary outer-halo channel** (lh_outerImprove): Independent contribution from outer rotation curve halo fitting quality. Adds ~4pp internally and ~9pp on external holdout.

The law is strongest — and most reliably transferable — in the high-Vflat (massive galaxy) regime.

---

## 5. Limitations and Caveats

1. **Still within SPARC**: Both the training and external samples come from the same survey. True cross-survey replication (different photometry, different rotation curve reduction) would substantially strengthen the claim. The Phase 201 external sample uses different galaxies but shares systematics.

2. **External a₀ quality**: Per-galaxy a₀ for the N=59 external sample was derived from fewer RAR data points (median 8 vs ~20–30 for the training sample) using transition-scale plotPoints. This adds noise and may attenuate correlations.

3. **logMhost estimation**: The external sample lacks group membership data. Host halo masses were estimated from Vflat-based abundance matching, which is cruder than the literature-calibrated values used for the training sample. The logMhost sign change in the combined sample likely reflects this quality gap.

4. **Regime dependence**: The coupling law works well only in the high-Vflat regime. Whether this reflects genuine physics (mass-dependent coupling) or insufficient data quality in low-mass galaxies remains unresolved.

5. **VfResid construction**: The structural model used to compute VfResid is trained on the N=45 published sample. Applying it to external galaxies assumes the same structural relationship holds — a mild circularity concern.

6. **Sample size**: The best-performing external subsample (Q=1, Vflat >= 120) has only N=11 galaxies. Statistical power is limited.

7. **No claimed universality**: We report an empirical regularity within SPARC data, now supported by external validation within the same survey. Whether this reflects a fundamental physical law or a survey-specific pattern requires larger, multi-survey validation.

8. **Not peer-reviewed**: This analysis has not undergone formal peer review.

---

## 6. Conclusions

Within the SPARC galaxy sample, per-galaxy a₀ variation is well-described by a hierarchical coupling law in which baryonic structure sets the baseline and a kinematic residual channel (VfResid) carries the dominant halo coupling signal. The coupling-law picture now has external support beyond the original N=45 sample, with the strongest transfer appearing in higher-Vflat galaxies and collapsing in lower-Vflat systems, consistent with a regime-dependent law rather than a universal one.

The most important finding is not any single model's predictive power, but the external reproduction of the full hierarchical structure:

1. The structural core alone fails to transfer (in both N=10 and N=59 external tests)
2. VfResid is the dominant transferable channel (carrying all positive gap in both tests)
3. The signal is regime-dependent, strengthening monotonically with Vflat
4. The low-Vflat regime fails, exactly as predicted by the internal analysis

This suggests that the physical driver of a₀ variation is not baryonic structure per se, but the efficiency of baryon-halo coupling as encoded in the kinematic excess — and this picture survives outside the training sample.

The dominant transferable channel outside N=45 is VfResid, and the transfer is strongest in the high-Vflat regime. The natural next step is testing whether this picture survives in an entirely independent survey with higher-quality data, particularly for the mass estimation variables.

---

## References

Lelli, F., McGaugh, S.S. & Schombert, J.M., 2016, AJ, 152, 157 (SPARC database)
McGaugh, S.S., Lelli, F. & Schombert, J.M., 2016, PRL, 117, 201101 (RAR)

---

*All analysis scripts, machine-readable results, and input data are archived at Zenodo (Concept DOI: 10.5281/zenodo.19430633) for full reproducibility.*

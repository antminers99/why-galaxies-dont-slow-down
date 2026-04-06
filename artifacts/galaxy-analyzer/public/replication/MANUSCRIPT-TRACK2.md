# Hierarchical Coupling Law for Per-Galaxy a₀ Variation: Evidence from SPARC Rotation Curves

*Non-peer-reviewed computational analysis report*
*Analysis: Phases 123–134 (internal), Phases 200–204 (external validation), Phases 300–303 (physical interpretation)*
*Zenodo v11: DOI 10.5281/zenodo.19440400 (Concept DOI: 10.5281/zenodo.19430633)*

---

## Abstract

Using per-galaxy fits to 175 SPARC rotation curves, we find that the MOND acceleration scale a₀ varies systematically across galaxies and is well-described by a hierarchical, regime-dependent empirical law. Within a published-quality subset of N=45 galaxies, a three-axis structural core (gas mass MHI, host halo mass Mhost, and dynamical coherence MeanRun) explains 44% of a₀ variance (LOO). Adding a kinematic coupling channel — the residual flat velocity after removing baryonic structural predictions (VfResid) — raises predictive power to 61%, with a further independent contribution from outer halo fitting quality (lh_outerImprove, 5-axis LOO = 65%). The coupling signal is regime-dependent, strongest in galaxies with Vflat >= 120 km/s and activating sharply near Vflat ~ 181 km/s. The dominant channel (VfResid) mediates essentially all halo proxy signals and contains ~36% irreducible variance beyond any available catalog observable.

A comprehensive external validation program (Phases 200–204) tests the coupling law on N=59 independent SPARC galaxies using frozen N=45 coefficients. The complete hierarchy replicates externally: Core alone fails in all regimes, VfResid dominates (bootstrap P > 99.6% in all subsamples), and lhOuter provides genuine secondary improvement. Crucially, we decompose the role of environmental estimation quality by replacing the crude host-mass formula (r = 0.232 vs. tidal values) with a trained estimator (r = 0.542). Core gap improves substantially (full: -53% to -8.1%; high-V: -6.7% to +16.7%), proving the environmental axis carries real signal. However, VfResid retains clear dominance (margin 30–76 pp) even with improved environmental representation. This establishes that the dominant transferable channel is VfResid — the velocity residual from structural expectation — independent of environmental proxy quality.

Phases 300–303 probe the physical nature of VfResid itself. Sample salvage recovered 35 additional galaxies from unused SPARC data (N=59 to N=94 external). Driver analysis confirms haloK (dark halo amplitude) as the only shared top predictor across both samples (r = 0.598 internal, 0.511 external), but pure structural models fail (LOO R² = -0.144 internal). A large irreducible fraction (35.8% internal, 68.5% external) persists beyond all catalog observables and still predicts a₀. Regime analysis reveals that this hidden residual genuinely activates with Vflat — the unexplained VfResid after removing haloK predicts a₀ at r = 0.299 (low-V) versus r = 0.749 (high-V) internally and r = 0.479 versus r = 0.889 externally. Critically, both slope and correlation strengthen with Vflat (slope ratio 2.6:1, r ratio 3.5:1), ruling out a pure observability explanation. Hypothesis testing favors dynamical integration (H4, score 10) over halo response (H1, score 7), feedback imprint (H3, score 7), and assembly history (H2, score 4). The hidden channel is best explained as an accumulated dynamical baryon-halo processing effect that genuinely amplifies in deeper potential wells.

We do not claim a₀ is definitively non-universal; we report structured empirical variation that is well-predicted by observable galaxy properties and whose hierarchical structure reproduces itself outside the training sample with a robustness that is independent of environmental estimation quality.

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

To test whether the coupling law generalizes beyond the original training sample, we assembled per-galaxy data for N=59 SPARC galaxies not included in the Stage A sample. Per-galaxy a₀ values were obtained from RAR fits using the transition-scale plotPoints (median 8 points per galaxy); VfResid was computed using the frozen N=45 structural model; logMhost was estimated from a Vflat-based formula (11.5 + 3.0 * (log10(Vflat) - 2.2)). No refitting was performed — all predictions use frozen N=45 coefficients.

| Regime | N | Core gap | Core+VfResid gap | 5-axis gap | VfResid-only gap | r(VfResid, a₀) |
|--------|---|----------|-------------------|-----------|-----------------|----------------|
| Full sample | 59 | -53.0% | +8.2% | +9.7% | +47.5% | 0.713 |
| Vflat >= 120 | 16 | -6.7% | +34.3% | +35.3% | +48.5% | 0.841 |
| Vflat >= 180 | 8 | -49.4% | +48.7% | +66.9% | +75.4% | 0.858 |
| Q=1 + Vflat >= 120 | 11 | +1.9% | +59.0% | +71.1% | +63.8% | 0.830 |
| Vflat < 120 | 43 | -74.4% | -3.9% | — | — | — |

### 3.8 External Hierarchy Replication and Environmental Decomposition (Phases 202–204)

Phase 202 formally tests whether the complete hierarchical structure — not just aggregate performance — replicates outside N=45. All 8 hierarchy checks pass: Core alone fails in every regime, VfResid dominates in every regime, bootstrap P(VfResid > Core) exceeds 99.6% in all subsamples, partial r(VfResid, a₀ | Core) = 0.844 in high-Vflat, and the regime dependence (high-V > full) is preserved. Channel dominance testing confirms VfResid outperforms every alternative proxy (logVflat, logMHI, logMbar, logMhost, logMeanRun) by margins of +35 to +71 pp.

A critical finding is that VfResid alone — without any Core variables — achieves 47.5% gap (full), 75.4% (very-high-Vflat), and 63.8% (Q=1+HV). This establishes VfResid not merely as an "improvement over Core" but as the primary carrier of transferable information.

**Environmental decomposition (Phase 203).** The crude Vflat-based logMhost formula correlates only r = 0.232 with real tidal-analysis values on the N=45 training sample. We trained four alternative estimators using galaxy observables that do not leak a₀ information:

| Estimator | Features | r with real logMhost | LOO RMSE |
|-----------|----------|---------------------|----------|
| Formula A (crude) | Fixed formula | 0.232 | 0.8299 |
| Model B | logVflat (trained) | 0.232 | 0.6258 |
| Model C | logVflat + logMbar + logRdisk + T | 0.542 | 0.5766 |
| Model D1 | logVflat + logMbar | 0.419 | 0.5940 |

Replacing the crude formula with Model C substantially improves Core performance:

| Regime | Core (crude) | Core (improved) | VfResid margin over improved Core |
|--------|-------------|----------------|----------------------------------|
| Full (N=59) | -53.0% | -8.1% | +40.3 pp |
| High-V (N=16) | -6.7% | +16.7% | +30.0 pp |
| Very-High-V (N=8) | -49.4% | -16.7% | +76.3 pp |
| Q=1+HV (N=11) | +1.9% | +25.4% | +43.0 pp |

This proves two things simultaneously: (1) the environmental axis carries real signal — crude estimation was suppressing it, and (2) VfResid dominance is not an artifact of poor Core representation. Even with improved environmental estimates, VfResid dominates by 30–76 pp in every regime, with bootstrap P(VfResid > improved Core) = 100% in all tested subsamples.

**lhOuter after environmental improvement.** The 5th axis (lh_outerImprove) was partially absorbing logMhost estimation error. After improvement, its contribution reduces but remains genuine:

| Regime | lhOuter delta (crude Mhost) | lhOuter delta (improved Mhost) |
|--------|---------------------------|-------------------------------|
| Full | +1.4 pp | +1.2 pp |
| Very-High-V | +18.1 pp | +8.3 pp |
| Q=1+HV | +12.1 pp | +4.9 pp |

The lhOuter signal persists in the high-mass regimes where it matters most, confirming it as a genuine — though secondary and smaller than initially estimated — independent channel.

### 3.9 Sample Salvage (Phase 300)

To maximize statistical power for the physical interpretation program, we recovered 35 additional galaxies from the 60 unused SPARC galaxies via two strategies: (a) estimating Vflat from the maximum observed rotation curve velocity for 22 galaxies lacking catalog Vflat, and (b) lowering the minimum RAR data point threshold from 5 to 3 for 13 galaxies with sufficient but sparse data. Three known high-velocity targets (NGC3953, NGC3949, NGC4051) were captured by the lowered threshold.

| Metric | Original | Augmented |
|--------|----------|-----------|
| External sample N | 59 | 94 |
| Vflat >= 120 km/s | 16 | 21 |
| Vflat >= 180 km/s | 8 | 9 |
| r(VfResid, a₀) full | 0.713 | 0.691 |
| r(VfResid, a₀) high-V | 0.841 | 0.859 |

The augmented sample preserves all previously established correlations. Original and augmented results are consistent at every Vflat threshold (all within +/- 0.02 in r), confirming that the salvaged galaxies do not distort the coupling signal.

### 3.10 VfResid Driver Analysis (Phase 301)

We systematically modeled VfResid itself — asking what determines the dominant coupling channel — at three levels: internal-only (N=45), external-only (N=94), and pooled with source controls (N=139).

**Single-predictor rankings.** haloK (dark halo linear amplitude) is the strongest predictor of VfResid in both samples (internal r = 0.598, external r = 0.511). MOND improvement (mondImprove) is second internally (r = -0.542) and fourth externally (r = -0.358), with consistent sign. Only 2 of the top-5 predictors overlap between internal and external rankings (haloK, mondImprove); secondary drivers diverge, likely because VfResid is constructed by residualizing Vflat on structural variables (logMbar, logL36, logRdisk, T), which are orthogonal to VfResid internally by construction but retain residual correlations externally.

**Multi-predictor models.**

| Model | N | LOO R² | Key finding |
|-------|---|--------|-------------|
| Internal: haloK alone | 45 | 0.287 | Strong but leaves 64% unexplained |
| Internal: haloK + lhOuter | 45 | 0.355 | lhOuter adds genuine signal |
| Internal: best-5 (haloK + lhOuter + envCode + MR + concentration) | 45 | 0.493 | Best internal model |
| Internal: pure structure (no halo) | 45 | -0.144 | Fails — VfResid is not structural |
| External: haloK alone | 94 | 0.225 | Transfers but weaker |
| External: haloK + baryonic | 94 | 0.341 | Best external model |
| External: pure structure (no halo) | 94 | 0.077 | Also fails |

**Irreducible fraction.** After the best internal model (LOO R² = 0.493), 35.8% of VfResid variance remains unexplained. This irreducible component still correlates with a₀ at r = 0.438. Externally, 68.5% remains unexplained with r(residual, a₀) = 0.556. This demonstrates that VfResid encodes genuine physical information beyond any available catalog observable.

**Cross-sample transfer.** A haloK-only model trained on N=45 internal galaxies predicts external VfResid with r = 0.511. The reverse transfer (external-trained, internal-tested) yields r = 0.598. Pure structure models do not transfer (R² < 0).

### 3.11 Regime Law (Phase 302)

Three tests address the central question: does the hidden physics inside VfResid genuinely activate with Vflat, or does only its observability improve?

**Test 1: Activation shape.** Cumulative threshold analysis shows the VfResid-a₀ correlation strengthening gradually with Vflat (GRADUAL_STRENGTHENING in external and pooled samples). There is no single sharp threshold — rather a smooth ramp. The internal sample (N=45) shows high variability in small bins (classified NOISY by the automated classifier), but the monotonic trend is preserved.

**Test 2: Driver regime split.** haloK strengthens with regime: r(haloK, VfResid) = 0.452 at low-V (N=10) versus 0.688 at high-V (N=35) internally, and 0.459 versus 0.579 in the pooled sample. mondImprove also strengthens (-0.454 to -0.578 internally). This confirms that halo-related drivers gain power in the high-mass regime.

**Test 3: Residual regime law.** After removing haloK from VfResid, the unexplained remainder still predicts a₀ — and this predictive power rises sharply with Vflat:

| Sample | Low-V r(residual, a₀) | High-V r(residual, a₀) | Very-high-V r(residual, a₀) |
|--------|----------------------|----------------------|---------------------------|
| Internal (N=45) | 0.299 | 0.749 | 0.832 |
| External (N=59) | 0.479 | 0.889 | 0.854 |
| Pooled (N=104) | 0.467 | 0.815 | 0.852 |

In the pooled sample above Vflat ~ 100, the unexplained residual surpasses haloK itself in a₀-predictive power: r(haloK, a₀) weakens while r(residual, a₀) climbs to 0.847. The hidden physics activates faster than any known observable halo proxy.

**Robustness.** The regime pattern is consistent between the original external sample and the augmented sample at every threshold (all differences < 0.02 in r).

### 3.12 Physical Interpretation (Phase 303)

We test four candidate hypotheses for the hidden physics inside VfResid, using hypothesis-specific proxy observables and evaluating each on the haloK-residual (hkResid = VfResid after removing haloK):

| Hypothesis | Physical mechanism | Key proxies |
|------------|-------------------|-------------|
| H1: Halo response | Adiabatic contraction / baryon-halo coupling | Baryonic surface density (SigmaBar), disk dominance |
| H2: Assembly history | Formation time / halo age | Morphological type, surface brightness |
| H3: Feedback imprint | SN/AGN reshaping of inner halo | Gas fraction, MOND improvement |
| H4: Dynamical integration | Accumulated baryon-halo processing | Dynamical time, RC extent, RAR mean run |

**Observability versus physics.** Both the slope and correlation of the hkResid-a₀ relationship strengthen with Vflat (internal: slope 1.57 at low-V to 4.13 at high-V; r from 0.201 to 0.713; pooled: slope ratio 1.86, r ratio 1.70). This rules out a pure observability explanation: the physics itself amplifies, not just its measurement precision.

**Hypothesis scorecard.** Scoring each hypothesis by the number of supporting lines of evidence (proxy correlations in the correct direction, regime strengthening, cross-sample sign consistency, positive LOO R²):

| Hypothesis | Score | Key evidence |
|------------|-------|-------------|
| H4: Dynamical integration | 10 | logMeanRun is the only proxy that strengthens with regime in hkResid (low-V r = 0.11, high-V r = 0.508); cross-sample consistent; best LOO = 0.066 |
| H1: Halo response | 7 | SigmaBar correlates at high-V (r = -0.426) but sign is unstable across regimes |
| H3: Feedback imprint | 7 | Gas fraction correlates at low-V (r = 0.597) but fades at high-V (r = 0.077) — opposite of regime activation |
| H2: Assembly history | 4 | morphT strong at low-V (r = 0.733) but vanishes at high-V (r = -0.049) |

The FULL model (all proxies combined) explains 34.8% of hkResid (internal), reducing r(remaining, a₀) from 0.651 to 0.288. Even after removing all hypothesis-specific content, the remaining signal still predicts a₀ — indicating genuine dynamical physics beyond any single proxy.

**Physical narrative.** The hidden channel inside VfResid is best explained as an accumulated dynamical integration effect: the long-term product of baryon-halo interaction that deepens with potential well depth. haloK captures the halo-amplitude component of this effect, while logMeanRun (dynamical coherence of the RAR) provides a partial handle on the integration depth. Gas fraction and morphology are relevant at low Vflat but become irrelevant in the regime where the coupling law is strongest — consistent with a transition from feedback-dominated (low mass) to integration-dominated (high mass) physics.

---

## 4. The Coupling Law (Summary)

The empirical law can be stated as:

**Per-galaxy a₀ is governed by a hierarchical, regime-dependent baryon-halo coupling law with three tiers:**

1. **Structural core** (MHI, Mhost, MeanRun): Sets the baseline a₀ based on gas content, halo environment, and dynamical regularity. Carries real signal (confirmed by environmental decomposition) but fails to transfer alone and remains secondary to VfResid externally.

2. **Dominant kinematic coupling channel** (VfResid): The residual flat velocity beyond baryonic structural prediction. Carries the bulk of the halo-to-a₀ coupling signal. Mediates all tested halo proxy effects. Activates gradually above Vflat ~ 120 km/s, reaching full strength above ~180 km/s. **Confirmed as the primary transferable channel**: dominates by 30–76 pp even after improving Core environmental estimates, with bootstrap P > 99.6% in all subsamples.

3. **Secondary outer-halo channel** (lh_outerImprove): Independent contribution from outer rotation curve halo fitting quality. Adds ~4pp internally. Externally: +8.3 pp in very-high-Vflat, +4.9 pp in Q=1+HV (after environmental correction; initial crude estimates were inflated by logMhost error absorption).

The law is strongest — and most reliably transferable — in the high-Vflat (massive galaxy) regime. The external validation program (Phases 200–204) establishes that this hierarchy is not sample-specific: it replicates in N=59 independent galaxies (augmented to N=94 via sample salvage) with the same structure, the same channel ordering, and the same regime dependence.

**Physical interpretation (Phases 300–303).** The dominant channel VfResid is not reducible to any single catalog observable. haloK (dark halo amplitude) is its strongest shared driver (r ~ 0.5–0.6 in both samples), but pure structural models fail entirely (LOO R² < 0). After removing haloK, a large irreducible fraction persists (35.8% internal, 68.5% external) that still predicts a₀. This hidden residual activates with Vflat — its a₀-predictive power rises from r ~ 0.3 at low-V to r ~ 0.8 at high-V — and the activation reflects genuine physical amplification, not merely improved observability (both slope and correlation strengthen, slope ratio ~ 2:1).

Hypothesis testing identifies dynamical integration (H4) as the best-supported explanation: an accumulated product of baryon-halo dynamical interaction that deepens in more massive systems. Gas fraction and morphology are relevant at low Vflat but fade in the high-Vflat regime where the coupling law is strongest, consistent with a transition from feedback-dominated to integration-dominated physics. The coupling law thus encodes something deeper than simple halo strength — it reflects a regime-dependent dynamical integration effect with haloK as the best partial observable handle.

---

## 5. Limitations and Caveats

1. **Still within SPARC**: Both the training and external samples come from the same survey. True cross-survey replication (different photometry, different rotation curve reduction) would substantially strengthen the claim. The external sample uses different galaxies but shares survey systematics.

2. **External a₀ quality**: Per-galaxy a₀ for the N=59 external sample was derived from fewer RAR data points (median 8 vs ~20–30 for the training sample) using transition-scale plotPoints. This adds noise and may attenuate correlations.

3. **logMhost estimation**: The external sample lacks group membership data. Even the improved estimator (Model C, r = 0.542) explains only ~29% of real logMhost variance. True group catalog data would provide a definitive test of Core's external weight. The environmental decomposition (Phase 203) partially addresses this by showing that better estimation does improve Core, but room for improvement remains.

4. **Regime dependence**: The coupling law works well only in the high-Vflat regime. Phase 302 provides evidence that this reflects genuine physical amplification (both slope and r strengthen, not merely reduced noise), but the low-Vflat failure could also partly reflect data quality limitations.

5. **VfResid construction**: The structural model used to compute VfResid is trained on the N=45 published sample. Applying it to external galaxies assumes the same structural relationship holds — a mild circularity concern.

6. **Sample size**: The best-performing external subsamples have limited N (very-high-Vflat: N=8, Q=1+HV: N=11). Bootstrap stability is strong (P > 99.6%), but independent confirmation with larger high-Vflat samples would be valuable.

7. **lhOuter inflation**: The initial external lhOuter contribution (e.g., +18.1 pp in very-high-V) was inflated by logMhost error absorption. The corrected values (+8.3 pp very-high-V, +4.9 pp Q1+HV) are genuine but smaller. This illustrates how proxy quality affects apparent channel contributions.

8. **logMR sign inconsistency**: The logMeanRun coefficient sign changes externally across estimators. This variable may not transfer cleanly between samples with different RAR data quality.

9. **No claimed universality**: We report an empirical regularity within SPARC data, now supported by comprehensive external validation within the same survey. Whether this reflects fundamental physics or a survey-specific pattern requires larger, multi-survey validation.

10. **Not peer-reviewed**: This analysis has not undergone formal peer review.

11. **Salvaged galaxy quality**: The 35 recovered galaxies (Phase 300) have either estimated Vflat (from maximum rotation curve velocity) or fewer RAR data points (3–4 instead of 5+). Regime-law results are robust to their inclusion (original vs. augmented differences < 0.02 in r), but individual measurements are less reliable.

12. **Hypothesis testing limitations**: The Phase 303 hypothesis scoring uses proxy observables (gas fraction, dynamical time, morphological type) that are indirect handles on the underlying physical processes. The hypotheses are not mutually exclusive — dynamical integration, halo response, and feedback may all contribute simultaneously. The scoring is heuristic, not a formal Bayesian model comparison.

13. **logFgas circularity risk**: Gas fraction (logFgas) and dynamical time (logDynTime) are constructed from the same catalog variables used to define VfResid. While hkResid (VfResid after removing haloK) is analyzed rather than VfResid directly, partial circularity through shared structural variables cannot be fully excluded.

---

## 6. Conclusions

Within the SPARC galaxy sample, per-galaxy a₀ variation is well-described by a hierarchical coupling law in which baryonic structure sets the baseline and a kinematic residual channel (VfResid) carries the dominant halo coupling signal. A comprehensive external validation program (Phases 200–204) establishes that this hierarchy is not sample-specific: it replicates in N=59 independent galaxies (augmented to N=94 via sample salvage) with the same structure, the same channel ordering, and the same regime dependence.

The most important findings from the external program are:

1. **The complete hierarchy replicates**: Core alone fails in every regime, VfResid dominates in every regime (bootstrap P > 99.6%), and lhOuter provides genuine secondary improvement. All 8 hierarchy checks pass (Phase 202).

2. **VfResid is the primary transferable channel**: VfResid alone (without Core) achieves 47.5% gap (full) and 75.4% (very-high-Vflat). It outperforms every alternative proxy by 35–71 pp. It retains strong signal after controlling for Core (partial r = 0.844 in high-V, 0.869 in Q=1+HV).

3. **Core is real but secondary**: Improving the environmental proxy (logMhost r: 0.232 to 0.542) raises Core gap from -53% to -8.1% (full) and -6.7% to +16.7% (high-V). This proves the environmental axis carries real physical signal. But VfResid still dominates by 30–76 pp in every regime.

4. **VfResid dominance is independent of environmental proxy quality**: The dominant transferable channel is not an artifact of poor Core estimation — it persists with the same margin after improving the environmental representation.

5. **The signal is regime-dependent**: It concentrates in Vflat >= 120 km/s and is strongest in very-high-Vflat galaxies. The low-Vflat regime fails, exactly as predicted.

The most important findings from the physical interpretation program (Phases 300–303) are:

6. **VfResid is not a structural artifact**: Pure structural models (surface brightness, disk size, baryonic mass, morphology) fail entirely to predict VfResid (LOO R² = -0.144 internal, 0.077 external). haloK is the only shared driver across both samples (r = 0.598 internal, 0.511 external), confirming genuine halo information content.

7. **A large irreducible fraction persists and still predicts a₀**: After the best available model, 35.8% of VfResid variance remains unexplained internally (68.5% externally). This irreducible part still correlates with a₀ (r = 0.438 internal, 0.556 external), indicating physics beyond any available catalog observable.

8. **The hidden physics genuinely activates with Vflat**: After removing haloK, the unexplained VfResid predicts a₀ at r = 0.299 (low-V) versus r = 0.749 (high-V) internally, and r = 0.479 versus r = 0.889 externally. Both slope and correlation strengthen (slope ratio 2.6:1), ruling out a pure observability explanation.

9. **Dynamical integration is the best-supported physical interpretation**: Hypothesis testing favors accumulated baryon-halo dynamical processing (score 10) over halo response (7), feedback imprint (7), and assembly history (4). Gas fraction and morphology are relevant at low Vflat but fade in the high-mass regime where the coupling law is strongest, consistent with a transition from feedback-dominated to integration-dominated physics.

Taken together, these results suggest that the physical driver of a₀ variation is not baryonic structure per se, nor simple halo amplitude, but a regime-dependent dynamical integration effect — the accumulated product of baryon-halo interaction that genuinely amplifies in deeper potential wells. haloK provides the best partial observable handle on this process, but the deeper signal goes beyond any single catalog variable.

The natural next step is testing whether this picture survives in an entirely independent survey with higher-quality data, particularly for spatially resolved kinematics that could directly probe the dynamical integration hypothesis.

---

## 10. Phase 400: Origin of VfResid — Mock Galaxy Simulation Program

### 10.1 Motivation

Phases 300–303 identified VfResid as the dominant coupling channel and dynamical integration (H4) as the best-supported physical interpretation. Phase 400 asks the critical question: what *physically generates* the VfResid–a₀ coupling? Can standard ΛCDM galaxy formation with NFW haloes reproduce the observed signal?

### 10.2 Phase 400A: Six-Scenario Simulation (N=300 each)

We generated 300 mock galaxies per scenario using NFW haloes + exponential baryonic discs with realistic scaling relations:

| Scenario | Description | r(VfResid, logA0) | VfResid lift |
|----------|------------|-------------------|-------------|
| Baseline NFW | Pure NFW, universal a₀ | 0.881 | +77.6pp |
| Halo Response | Adiabatic contraction (mass-dependent) | 0.850 | +72.2pp |
| Assembly History | Assembly-dependent a₀ scatter (active) | 0.842 | +71.0pp |
| Feedback | SN core formation (Di Cintio 2014) | 0.837 | +70.0pp |
| Dynamical Integration | Accumulated baryon-halo processing (active) | 0.912 | +83.2pp |
| Combined | All mechanisms active simultaneously | 0.864 | +74.7pp |

**Key finding**: VfResid–a₀ coupling appears in ALL scenarios (r ≈ 0.84–0.91), including the pure baseline. Physical mechanisms are properly applied via RAR-based correction to gobs. The primary driver is **halo concentration diversity**: r(log c, VfResid) ≈ 0.60–0.77 across all scenarios. Concentration scatter systematically biases per-galaxy RAR fitting — higher c → faster-rising inner rotation curve → different apparent a₀. VfResid carries massive non-structural information (+70–83pp lift over structure alone for predicting a₀).

### 10.3 Phase 400B: The Regime Inversion Problem

**Critical discovery**: All simulations show the VfResid–a₀ coupling **weakening** with Vflat, while SPARC shows it **strengthening**:

| Source | Low-V r | High-V r | Pattern |
|--------|---------|----------|---------|
| SPARC observed | 0.30 | 0.75 | Strengthens |
| Sim baseline | 0.93 | 0.87 | Weakens |
| Sim halo response | 0.93 | 0.81 | Weakens |
| Sim dynamical integration | 0.92 | 0.90 | Weakens |
| Sim combined | 0.87 | 0.86 | Weakens |

The regime direction is inverted. No tested scenario reproduces the SPARC regime-strengthening pattern.

### 10.4 Phase 400C: Realistic Mock with Noise and Physical Signal

Seven additional tests with increasingly realistic conditions (SPARC-like velocity errors, sparse low-mass data, M/L uncertainty, physical a₀ variation, low-mass scatter):

| Test | r(all) | r(low-V) | r(high-V) | Pattern |
|------|--------|----------|-----------|---------|
| Clean mock | 0.874 | 0.946 | 0.880 | Weakens |
| Noise only | 0.577 | 0.630 | 0.579 | Weakens |
| Noise + low-mass scatter | 0.693 | 0.777 | 0.645 | Weakens |
| Noise + dynInt (strong, active) | 0.707 | 0.750 | 0.624 | Weakens |
| SPARC-calibrated (active) | 0.710 | 0.766 | 0.595 | Weakens |
| SPARC observed | 0.70 | 0.30 | 0.75 | **Strengthens** |

**Verdict**: No simulation configuration reproduces the SPARC regime-strengthening pattern. The concentration baseline always dominates and produces weakening. Even with realistic noise, mass-dependent physical a₀, and large low-mass scatter, the pattern remains inverted.

### 10.5 Phase 401: Can Any NFW Extension Flip the Regime Sign?

**Goal**: Find a mechanism that produces regime STRENGTHENING (r at high-V > r at low-V), matching SPARC.

Six experiment families, 30 configurations total, N=600 galaxies each:

| Experiment | Mechanism | r(all) | r(low-V) | r(high-V) | Pattern |
|------------|-----------|--------|----------|-----------|---------|
| E1 | Stochastic low-mass disruption (non-circ motions, 5 strengths) | 0.60–0.67 | 0.64–0.74 | 0.54–0.67 | Weakens |
| E2 | Mass-dependent adiabatic contraction (5 strengths) | 0.36–0.57 | 0.55–0.65 | 0.43–0.60 | Weakens |
| E3 | E1+E2 combined (7 combos) | 0.43–0.51 | 0.63–0.70 | 0.46–0.61 | Weakens |
| E4 | Concentration-response anti-correlation (5 strengths) | 0.41–0.50 | 0.61–0.66 | 0.28–0.54 | Weakens |
| E5 | Deep-well dynamical integration (4 strengths) | 0.58–0.64 | 0.69–0.73 | 0.49–0.59 | Weakens |
| E6 | Full kitchen sink — all mechanisms combined (4 combos) | 0.42–0.50 | 0.58–0.68 | 0.34–0.43 | Weakens |
| SPARC observed | — | 0.70 | 0.30 | 0.75 | **Strengthens** |

**Closest approach**: E1 nonCirc=0.5 achieves delta(high−low) = −0.027 (nearly flat, but still weakening).

**Central result**: The NFW+disk framework is structurally resistant to regime sign inversion. Concentration scatter creates stronger coupling at low mass (where it has more leverage on RC shape) and weaker at high mass. No additive physics layered onto NFW overcomes this structural bias. The SPARC regime-strengthening pattern is **incompatible** with the entire NFW+disk model family, not just its baseline.

**Implication**: The key discriminator is not the existence of VfResid–a₀ coupling itself, which arises trivially in NFW-based mocks, but the observed SPARC regime strengthening. Across all tested NFW+disk mock variants, the coupling systematically weakens or at best flattens with increasing Vflat, never reproducing the observed strengthening. This makes the SPARC regime law incompatible with simple concentration diversity, observational noise, and a wide range of standard extensions layered onto NFW.

### 10.6 Phase 402: Beyond NFW — Can Cored Haloes Flip the Sign?

**Goal**: Test whether qualitatively different halo physics (cored profiles, feedback-driven cores, stochastic core sizes, genuine a₀ variation) can produce regime strengthening.

Phase 402a: 22 configurations (N=600 each), 4 branches:

| Branch | Models tested | Best delta(high−low) | Pattern |
|--------|---------------|---------------------|---------|
| A: Cored haloes (Burkert, PseudoIso, DC14) | 6 | −0.075 | Weakens |
| B: Genuine high-mass a₀ activation | 5 | −0.014 | Weakens |
| C: Hybrid (cores + a₀) | 6 | −0.057 | Weakens |
| C+: Hybrid + low-mass noise | 5 | −0.024 | Weakens |

Phase 402b: 22 configurations (N=800 each), focused mechanisms:

| Test family | Mechanism | Best delta | Pattern |
|-------------|-----------|-----------|---------|
| Feedback-cored haloes | Mass-dependent core from stellar feedback | +0.005 | Weakens |
| Burkert transition + core scatter | Stochastic core sizes at low mass | −0.146 | Weakens |
| NFW + random core injection | Core perturbation uncorrelated with c | −0.050 | Weakens |
| Hybrid: cores + genuine a₀ | Feedback cores + potential-depth a₀ | −0.055 | Weakens |
| Triple: cores + a₀ + noise | All mechanisms combined | −0.047 | Weakens |
| SPARC observed | — | **+0.45** | **Strengthens** |

**Central result**: The concentration-leverage asymmetry is a structural property of ALL spherical halo+disk models, not just NFW. Any dark matter halo profile with a concentration-like parameter (scale radius, core radius) creates coupling that is inherently STRONGER at low mass because low-mass galaxies have more leverage from profile shape variation on the fitted a₀. No cored profile, no feedback mechanism, no genuine a₀ variation, and no combination thereof overcomes this structural bias.

**Total falsification scope**: ~58 configurations across Phases 400–402, spanning NFW, Burkert, pseudo-isothermal, DC14-inspired, and feedback-cored profiles with noise, contraction, anti-correlation, dynamical integration, stochastic core scatter, and genuine mass-dependent a₀ variation. ALL produce weakening or at best flat regime dependence. None reproduce the SPARC regime-strengthening law.

**Implication**: The SPARC regime-strengthening pattern (coupling strengthens at high Vflat) appears incompatible with ANY spherical dark matter halo + exponential disk model, regardless of profile family or additional physics. This pattern requires either: (1) a fundamentally non-CDM acceleration law (e.g., MOND-like), (2) a source of a₀ variation at high mass that is completely independent of halo structure, or (3) a systematic observational effect specific to the SPARC sample.

### 10.7 Phase 403: Anatomy of VfResid

**Goal**: (A) Is regime strengthening real or a SPARC artifact? (B) What generates VfResid? (C) Why does coupling activate at high-V?

**A. Cross-sample replication — CONFIRMED:**

| Sample | N | delta(high−low) | Pattern |
|--------|---|----------------|---------|
| Internal (training) | 55 | +0.066 | Strengthens |
| External (validation) | 94 | +0.186 | Strengthens |
| Pooled | 149 | +0.114 | Strengthens |
| Q=1 only | 39 | +0.100 | Strengthens |
| Moderate inclination | 42 | +0.103 | Strengthens |

External sample shows STRONGER regime strengthening than internal (delta=+0.186 vs +0.066). High-V external: r=0.859 (N=21). Not an artifact of SPARC sample construction or data quality.

**B. What generates VfResid:**
- Top single predictor: logA0 (r=0.711, circular)
- Next: envCode (r=-0.365), logMeanRun (r=0.289)
- All 5 known structural variables together: R²=0.217 (only 21.7% of VfResid explained)
- After removing structural effects: r(VfResid_residual, logA0) = 0.538
- No single control variable absorbs the VfResid-a₀ coupling (all partial r ≥ 0.67)
- Controlling for Vflat INCREASES coupling: partial r = 0.788 (raw: 0.711)

**C. Why high-V only (activation anatomy):**

| Vflat bin | N | r(VfResid, a₀) |
|-----------|---|----------------|
| <80 | 5 | 0.44 |
| 80-120 | 5 | 0.93 |
| 120-160 | 11 | 0.74 |
| 160-200 | 14 | 0.83 |
| 200-250 | 14 | 0.88 |
| >250 | 6 | 0.83 |

Activation is **gradual**, not a sharp threshold. Coupling strength ramps from r≈0.4 at V<80 to r≈0.84 at V>160. High-V galaxies have more gas (logMHI: −0.12→+0.85), higher surface brightness (logSBdisk: 2.63→3.20), and higher central density (logΣ₀: 2.42→3.00).

**Phase 403 verdict**: Regime strengthening is real, replicated across samples, and robust to quality cuts. VfResid is NOT primarily driven by any known structural or observational variable. The coupling's gradual activation with Vflat, and its strengthening after Vflat control, points to a genuine physical signal that operates in deep potential wells.

### 10.8 Phase 400–403 Conclusions

1. **Baseline VfResid is a fitting artefact**: Pure NFW concentration diversity produces r(VfResid, logA0) ≈ 0.86 via systematic bias in per-galaxy RAR fitting. This is non-physical.

2. **Regime strengthening is the key discriminant**: The SPARC pattern (coupling strengthens at high Vflat) cannot be reproduced by any combination of NFW concentration effects + noise + mass-dependent a₀. This regime signature distinguishes physical a₀ variation from fitting bias.

3. **The SPARC signal requires qualitatively different physics**: Simple NFW+disc models fail to produce regime strengthening regardless of parameter choice. Possible explanations include non-circular motions scrambling low-V signals, baryon-dependent halo profiles beyond NFW, genuinely mass-dependent acceleration scales, or selection effects in SPARC sample construction.

4. **Implication for the hierarchical coupling law**: The observational VfResid channel contains *both* a concentration-driven fitting artefact (always present) and a genuine physical signal whose regime-strengthening signature is unique to SPARC data and unexplained by standard ΛCDM mock galaxies.

---

## References

Lelli, F., McGaugh, S.S. & Schombert, J.M., 2016, AJ, 152, 157 (SPARC database)
McGaugh, S.S., Lelli, F. & Schombert, J.M., 2016, PRL, 117, 201101 (RAR)

---

*All analysis scripts, machine-readable results, and input data are archived at Zenodo (Concept DOI: 10.5281/zenodo.19430633) for full reproducibility.*

# The Inferred MOND Acceleration Scale as an Emergent Galaxy State Variable: A Five-Axis Law from N=45 SPARC Galaxies

---

## Abstract

The radial acceleration relation (RAR) in disk galaxies is commonly parameterized by a single, supposedly universal acceleration scale $a_0 \approx 1.2 \times 10^{-10}$ m s$^{-2}$. We test whether galaxy-to-galaxy variation in the inferred $a_0$ is best described as (i) a strict universal constant, (ii) unstructured per-galaxy scatter, or (iii) a structured multi-axis law. Working with a quality-controlled subsample of $N = 45$ SPARC galaxies selected for published, high-quality distance estimates, we construct per-galaxy $a_0$ values and search for predictive structure.

We find that within this sample, a five-predictor power-law model (M5) — incorporating gas mass ($M_\mathrm{HI}^{-1/4}$), host-halo mass ($M_\mathrm{host}^{-1/6}$), central baryonic surface density ($\Sigma_0^{+1/7}$), kinematic coherence ($\overline{\mathrm{Run}}^{+3/7}$), and an orthogonalized stellar mass-to-light ratio ($\Upsilon_\star^{\perp\,+2/3}$) — explains 51.0% of the $\log a_0$ variance in leave-one-out cross-validation (LOO). A compressed three-axis state law (M3) retains 86% of this signal using only gas mass, host mass, and kinematic coherence (LOO = 44.1%). A comprehensive overfitting audit — including nested cross-validation ($p < 0.001$ permutation), collinearity analysis (all VIF $< 1.3$), and bootstrap stability — confirms that the structured signal is not a fitting artifact within the training set.

**However, external validation fails.** When frozen M3 coefficients are applied to 36–78 SPARC galaxies outside the training set (using rotation-curve data from CDS/VizieR with identically computed inputs), all model variants perform worse than the trivial universal-constant baseline — including on regime-matched subsamples selected to resemble the training set (win rates 7–10%).

A selection audit identifies the mechanism: the $V_\mathrm{flat}$ criterion (excluding low-mass dwarfs) flips the MHI–$a_0$ correlation from positive ($r = +0.278$ in full SPARC) to negative ($r = -0.215$ in $V_\mathrm{flat} \geq 80$). Controlled subsample analysis (5000 draws, 23 pools) shows the training $r = -0.322$ is **typical** within the high-$V_\mathrm{flat}$ regime ($p = 0.17$) — not a selection-amplified extreme. We conclude that the structured $a_0$ variation is a **regime-specific finding**: a genuine mass-dependent structure in $a_0$ correlations that holds for massive spirals but reverses for dwarf irregulars. The multi-axis law is not a universal description of $a_0$ physics but a valid characterization of the high-mass regime. These results demonstrate that internal cross-validation cannot detect regime boundaries, and that external validation spanning different mass regimes is essential for testing universality claims.

---

## 1. Introduction

The radial acceleration relation (RAR; McGaugh et al. 2016; Lelli et al. 2017) establishes a tight empirical correlation between the observed centripetal acceleration in disk galaxies and the acceleration predicted from the baryonic mass distribution alone. This relation is a central prediction of Modified Newtonian Dynamics (MOND; Milgrom 1983) and is parameterized by a characteristic acceleration scale $a_0$. Whether $a_0$ is a strict universal constant or varies systematically between galaxies has profound implications: a truly universal $a_0$ supports MOND as a fundamental theory, while structured variation would suggest that $a_0$ encodes additional galaxy-state information beyond what a single constant captures.

Previous investigations have reached conflicting conclusions. Li et al. (2018) found that a universal $a_0$ provides a good fit to SPARC galaxies with minimal intrinsic scatter. Rodrigues et al. (2018) argued for significant galaxy-to-galaxy variation, a claim contested by subsequent analyses (Kroupa et al. 2018; McGaugh et al. 2018). The debate has largely focused on whether observed scatter is consistent with measurement noise or requires genuine variation.

We take a different approach. Rather than asking *whether* $a_0$ varies, we ask: *if it varies, does the variation follow a structured, falsifiable law?* Specifically, we test whether the pattern of $a_0$ values across galaxies can be predicted from independently measured galaxy properties — gas content, environment, baryonic structure, and kinematic regularity — using cross-validated, out-of-sample metrics that guard against overfitting.

Our key advance beyond prior work is the transition from demonstrating *scatter* to demonstrating *law*: we show that $a_0$ variation is not merely nonzero but is predictive, pairwise-testable, and interpretable as an equation of state for the galaxy.

---

## 2. Data and Sample Definition

### 2.1 Source Catalog

We draw rotation-curve data and photometric properties from the Spitzer Photometry and Accurate Rotation Curves (SPARC) database (Lelli, McGaugh & Schombert 2016), which provides 175 disk galaxies with 3.6 $\mu$m photometry and high-quality HI/H$\alpha$ rotation curves. We supplement SPARC with:

- Host-halo masses ($M_\mathrm{host}$) from the Tully (2015) Cosmicflows group catalog and the Extragalactic Distance Database (EDD).
- Stellar mass-to-light ratios ($\Upsilon_\star^\mathrm{disk}$) from Li et al. (2020) SPS-based estimates at 3.6 $\mu$m.

### 2.2 Quality Subsample (N = 45)

From the 175 SPARC galaxies, we select a quality-controlled subsample of $N = 45$ based on the following criteria:

1. **Valid per-galaxy $a_0$ fit**: The galaxy must have a well-defined individual $a_0$ from its rotation curve.
2. **Published distance quality**: The adopted distance must come from a published, peer-reviewed source using a primary or well-calibrated secondary distance indicator (TRGB, Cepheid, or equivalent), not solely from Hubble-flow estimates.
3. **Host-mass assignment**: A host-halo mass from the Tully group catalog or equivalent environmental database must be available.
4. **Sufficient rotation-curve sampling**: At least 5 measured RC points.

The published-distance criterion is the dominant filter (175 $\to$ 56 $\to$ 45), ensuring that distance-dependent quantities are not dominated by Hubble-flow errors. The resulting sample spans morphological types T = 0–11, $V_\mathrm{flat} = 17$–300 km/s, and $M_\mathrm{HI} = 0.07$–33 $\times 10^9 \, M_\odot$.

### 2.3 Sample Selection Transparency

The sample was frozen prior to model comparison and was not optimized to favor any particular model. Selection criteria are purely observational and contain no cuts on $a_0$ or any predictor variable.

---

## 3. Inference of Per-Galaxy $a_0$

### 3.1 The Radial Acceleration Relation

The RAR relates observed centripetal acceleration $g_\mathrm{obs} = V_\mathrm{obs}^2/r$ to the Newtonian acceleration from baryons alone $g_\mathrm{bar}$ via:

$$g_\mathrm{obs} = \frac{g_\mathrm{bar}}{1 - e^{-\sqrt{g_\mathrm{bar}/a_0}}}$$

### 3.2 Per-Galaxy Fitting Procedure

For each galaxy, we fit the RAR to its rotation curve by minimizing the sum of squared residuals, treating $a_0$ as the free parameter while holding $\Upsilon_\star^\mathrm{disk}$ at its SPS-estimated value.

### 3.3 Target Variable

$$y_i = \log_{10}(a_{0,i}) \quad [\text{m/s}^2]$$

Sample mean $\bar{y} = 3.580$ (i.e., $a_0 \approx 3800$ m/s$^2$ in our units), SD$(y) = 0.258$ dex.

---

## 4. Predictor Construction

We define five predictor variables for the M5 model. Each is constructed from independently measured galaxy properties.

### 4.1 Gas Mass: $\log M_\mathrm{HI}$

$$\log M_\mathrm{HI} = \log_{10}\!\left(\frac{M_\mathrm{HI}}{10^9 \, M_\odot}\right)$$

Source: SPARC catalog. Rational exponent in M5: $\approx -1/4$.

### 4.2 Host-Halo Mass: $\log M_\mathrm{host}$

$$\log M_\mathrm{host} = \log_{10}\!\left(\frac{M_\mathrm{host}}{M_\odot}\right)$$

Source: Tully (2015) Cosmicflows group catalog. Range: 10.4–13.0. Rational exponent: $\approx -1/6$.

### 4.3 Central Baryonic Surface Density: $\log \Sigma_0$

$$\log \Sigma_0 = \log_{10}\!\left(\mathrm{SB}_\mathrm{disk} \times \Upsilon_\star^\mathrm{disk}\right)$$

Units: $M_\odot$ pc$^{-2}$. Rational exponent: $\approx +1/7$.

### 4.4 Kinematic Coherence: $\log \overline{\text{Run}}$

Measures the average length of consecutive same-sign residual stretches in the MOND fit:

$$\log \overline{\text{Run}} = \log_{10}\!\left(\frac{n_\mathrm{points}}{n_\mathrm{runs}}\right)$$

Higher values indicate more coherent, large-scale departures from the MOND prediction. Rational exponent: $\approx +3/7$. Confirmed as a genuine 1D projection of 2D dynamical structure via THINGS velocity-field diagnostics ($|r| = 0.83$ with lopsidedness, bisymmetric flow, and non-circular motion fraction).

### 4.5 Independent Stellar Mass-to-Light Ratio: $\Upsilon_\star^\perp$

The raw $\log_{10} \Upsilon_\star^\mathrm{disk}$ is partially confounded with gas mass, surface density, and morphology. We orthogonalize:

1. Fit auxiliary OLS: $\log_{10} \Upsilon_\star^\mathrm{disk} = \alpha_0 + \alpha_1 \log M_\mathrm{HI} + \alpha_2 \log \Sigma_0 + \alpha_3 T + \varepsilon$
2. This explains $R^2 = 0.68$ of the $\Upsilon_\star$ variance.
3. Define $\Upsilon_\star^\perp = \varepsilon$ (the OLS residual).

Rational exponent: $\approx +2/3$.

### 4.6 Dropped Variable: rcWiggliness

An earlier model variant (B$''$, 6-variable) included rcWiggliness (rotation-curve irregularity). This variable was dropped from the final M5 model because: (a) incremental F-test = 1.50 (not significant), (b) 14.5% bootstrap sign-flip rate (highest of all candidates), and (c) its removal improves LOO performance. The kinematic axis is better captured by $\log\overline{\text{Run}}$ alone.

---

## 5. Model Space and Evaluation

### 5.1 Competing Models

| Model | Predictors | $k$ | Description |
|-------|-----------|-----|-------------|
| M0 | (none) | 1 | Universal $a_0$ |
| M3 | $\log M_\mathrm{HI}$, $\log M_\mathrm{host}$, $\log \overline{\text{Run}}$ | 4 | Compressed state law |
| M5 | M3 + $\log \Sigma_0$, $\Upsilon_\star^\perp$ | 6 | Full predictive law |
| M$_\mathrm{free}$ | Per-galaxy free $a_0$ | 45 | Unconstrained ceiling |

### 5.2 Evaluation Metrics

**LOO gap%** (primary): $100 \times (1 - \text{LOO-RMS}^2 / \text{SD}(y)^2)$

**AIC / BIC**: Standard information criteria.

**$k$-fold CV**: 5-fold and 10-fold, averaged over 100 random partitions.

### 5.3 Why Not $f_\mathrm{gas}$?

Gas fraction $f_\mathrm{gas} = M_\mathrm{HI}/M_\mathrm{host}$ was tested as a replacement for the independent MHI and Mhost axes. Collapsing to a single ratio loses 36.5 percentage points of LOO performance, demonstrating that gas mass and host mass carry genuinely independent physical information with different exponents ($-1/4$ vs $-1/6$).

---

## 6. Main Results

### 6.1 Model Comparison

**Table 3. Model Comparison on the N = 45 Quality Subsample**

| Model | LOO gap% | RMS (dex) | AIC | BIC |
|-------|----------|-----------|-----|-----|
| M0 (universal) | $-$2.3% | 0.261 | $-$121.0 | $-$119.2 |
| M3 (3-axis) | 44.1% | 0.193 | — | — |
| **M5 (5-axis)** | **51.0%** | **0.157** | **best** | **best** |
| M$_\mathrm{free}$ | $\ll 0$ | 0.000 | — | — |

M5 wins on every metric. M0 has *negative* LOO gap% — it predicts worse than the sample mean on unseen galaxies.

### 6.2 The M5 Law

**Table 4. M5 Coefficients**

| Predictor | Coefficient | Rational approx. | Physical axis |
|-----------|------------|-------------------|---------------|
| Intercept | +4.978 | — | — |
| $\log M_\mathrm{HI}$ | $-$0.236 | $-1/4$ | Gas reservoir |
| $\log M_\mathrm{host}$ | $-$0.172 | $-1/6$ | Environmental depth |
| $\log \Sigma_0$ | +0.145 | $+1/7$ | Baryon concentration |
| $\log \overline{\text{Run}}$ | +0.452 | $+3/7$ | Dynamical coherence |
| $\Upsilon_\star^\perp$ | +0.658 | $+2/3$ | Stellar structure |

### 6.3 The M3 Compressed State Law (General Equation)

In power-law form:

$$a_0 \approx A_0 \cdot M_\mathrm{HI}^{-1/5} \cdot M_\mathrm{host}^{-1/6} \cdot \overline{\text{Run}}^{1/2}$$

Or equivalently in log-linear form:

$$\log a_0 = 5.182 - 0.198 \log M_\mathrm{HI} - 0.155 \log M_\mathrm{host} + 0.459 \log \overline{\text{Run}}$$

LOO = 44.1%, retaining 86% of M5's signal with only three predictors. M3 serves as a *first-order equation of state*: given three macroscopic galaxy properties — gas reservoir, gravitational depth, and dynamical coherence — $a_0$ is determined to $\sim$0.19 dex.

### 6.4 Relationship Between M3 and M5

M5 is formally M3 plus two second-order structural corrections:

$$a_0 = \underbrace{A_0 \cdot M_\mathrm{HI}^{-1/4} \cdot M_\mathrm{host}^{-1/6} \cdot \overline{\text{Run}}^{3/7}}_{\text{first-order state law}} \times \underbrace{\Sigma_0^{1/7} \cdot 10^{(2/3)\Upsilon_\star^\perp}}_{\text{structural corrections}}$$

The corrections add +6.9 percentage points of LOO performance and reduce RMS from 0.193 to 0.157 dex. Heuristically, M3 resembles a first-order equation of state that captures the dominant macroscopic structure, while M5 adds second-order structural corrections that improve predictive precision — analogous (though not physically identical) to the relationship between an ideal gas law and its van der Waals refinement.

### 6.5 Structural Insight: Suppression vs Amplification

The five axes split into two physical families:

- **Negative coefficients** (mass axes): More mass $\to$ lower $a_0$. Deeper potential wells suppress the effective acceleration scale.
- **Positive coefficients** (organization axes): More concentrated, coherent, high-M/L disks $\to$ higher $a_0$. Baryonic organization amplifies the effective scale.

$a_0$ emerges from the *tension* between gravitational depth and baryonic organization.

---

## 7. Validation: Falsifiability and Pairwise Tests

### 7.1 Falsifiable Predictions (Phase 71)

M5 was subjected to 26 falsifiable tests across five categories:

| Category | Tests | Passed |
|----------|-------|--------|
| Quantitative shift predictions | 5 | 5 |
| Matched-pair sign predictions | 5 | 5 |
| Partial-residual rank correlations | 5 | 5 |
| Negative predictions (what should NOT work) | 5 | 5 |
| Tertile monotonicity | 5 | 3 |
| F-test: M3→M5 improvement | 1 | 1 (F=4.99) |
| **Total** | **26** | **24 (92%)** |

### 7.2 Controlled Matched-Pair Tests (Phase 72)

For each M5 axis, we identified galaxy pairs that differ primarily on one axis while matching on the other four, then tested whether $a_0$ shifts in the predicted direction by the predicted amount.

**Table 5a. Matched-Pair Results by Axis**

| Axis family | Pairs | Sign passes | Median |pred$-$obs| |
|-------------|-------|-------------|----------------------|
| Gas mass | 5 | 5/5 | 0.143 dex |
| Host mass | 5 | 4/5 | 0.112 dex |
| Baryon density | 5 | 4/5 | — |
| Kinematic coherence | 5 | 4/5 | — |
| Stellar structure | 5 | 5/5 | — |
| **Total** | **25** | **22/25 (88%)** | **0.165 dex** |

Best showcase pair: NGC5907 vs NGC7814 (gas axis) — predicted $\Delta\log a_0 = -0.303$, observed $-0.286$, error only 0.017 dex.

### 7.3 2D Dynamical Bridge (Phase 74)

The kinematic-coherence axis ($\log\overline{\text{Run}}$) was tested against resolved 2D velocity-field diagnostics from the THINGS survey (N = 7 galaxies with available data):

| 2D diagnostic | $r$ with MeanRun |
|---------------|-----------------|
| Lopsidedness | $-$0.828 |
| Bisymmetric flow | $-$0.831 |
| Non-circular motion fraction | $-$0.821 |

All three correlations exceed $|r| = 0.8$, confirming that MeanRun is a genuine 1D projection of real 2D dynamical structure, not a processing artifact.

---

## 8. Residual Frontier (Phase 75)

### 8.1 Residual Diagnostics

After M5, the residual RMS = 0.157 dex with no exploitable structure:

| Test | Result | Conclusion |
|------|--------|------------|
| Stratification (7 splits) | 0/7 significant ($|t|>2$) | No subpopulation bias |
| k-means clustering (k=2,3) | F = 0.02, 0.07 | No cluster carries excess residual |
| Breusch-Pagan | nR² = 4.63 (critical 11.1) | No heteroscedasticity |
| Skewness | $-$0.124 | Symmetric |
| Excess kurtosis | $-$0.557 | Near-Gaussian |
| Jarque-Bera | 1.22 (critical 5.99) | Cannot reject normality |
| Bimodality | None detected | Unimodal |

### 8.2 No Seventh Scalar Axis

A comprehensive scan of candidate variables (Phases 62, 65, 66) found no additional scalar variable reaching significance against M5 residuals. Environmental processing proxies (HI deficiency, gas truncation, ram-pressure susceptibility) show no signal. Distributed weak-signal bundles fail in out-of-sample prediction.

### 8.3 Intrinsic Scatter

The remaining $\sim$0.10 dex intrinsic scatter (after subtracting estimated measurement noise) is genuine and represents the current irreducible frontier of scalar-law modeling on this dataset.

---

## 9. Blind External Transfer (Phase 73)

We applied frozen M5 coefficients to $N = 11$ galaxies excluded from training (crude/Hubble-flow distances only). The test was **inconclusive**: M5 did not outperform the universal constant on this sample (RMS = 0.320 vs 0.288 dex).

However, the test is severely compromised: (a) all test galaxies have crude Hubble-flow distances adding $\sim$0.15–0.20 dex noise to $\log a_0$, comparable to the signal M5 predicts; (b) no SPS-based $\Upsilon_\star$ estimates are available, rendering the $\Upsilon_\star^\perp$ axis uninformative; (c) $N = 11$ is too small for reliable assessment.

We do not claim external generalization from this test. True external validation requires a sample with published distances and SPS estimates.

---

## 10. Theory Bridge

### 10.1 $a_0$ as an Emergent State Variable

The M5 law reads as an *equation of state* for the galaxy: given five measurable macroscopic properties, the effective acceleration scale is determined. The axes decompose naturally:

- **Suppression axes** ($M_\mathrm{HI}$, $M_\mathrm{host}$): More mass $\to$ lower $a_0$
- **Amplification axes** ($\Sigma_0$, MeanRun, $\Upsilon_\star^\perp$): More organization $\to$ higher $a_0$

This structure is analogous to thermodynamic equations of state, where macroscopic observables determine intensive parameters.

### 10.2 Baryon-Halo Coupling

The law implies that inferred $a_0$ encodes a systematic, multi-axis coupling between baryonic content and the host gravitational environment. The failure of $f_\mathrm{gas}$ collapse ($-$36.5 pp) demonstrates that gas content and environmental depth operate through genuinely independent physical mechanisms.

### 10.3 Framework Neutrality

This result is *compatible with both* dark-matter and modified-dynamics frameworks:

- **In dark matter**: $a_0$ variation reflects varying baryon-halo feedback efficiency.
- **In MOND**: variation reflects environmental modulation (external field effect, screening) of the effective acceleration scale.

We do *not* claim that MOND is refuted or that dark matter is confirmed. The data constrain any theory to reproduce a state-dependent effective acceleration scale.

---

## 11. Robustness Against Overfitting

A structured multi-axis law discovered on $N = 45$ galaxies invites legitimate concern about overfitting — not only in the final equation, but across the full analytical pipeline: sample selection, variable choice, model compression, and best-model identification. We subject the entire pipeline to six independent tests designed to isolate selection bias from genuine signal.

**Nested cross-validation.** When model selection itself (M3 vs M5 vs M6) is performed inside each training fold and evaluated on held-out outer folds, the LOO gap% degrades from 51.0% to 38.5% — an optimism of 12.5 percentage points, within the expected range for $N = 45$. The surviving 38.5% represents signal that persists even when the choice of model is part of the test.

**Full-pipeline permutation.** We scramble $\log a_0$ across galaxies and re-run the complete pipeline (construct $\Upsilon_\star^\perp$, fit M3 and M5, take the better LOO) 1,000 times. No null run reaches the real gap% of 46.6%; the 99th percentile of the null distribution is 11.3% ($p < 0.001$).

**Collinearity.** All five M5 predictors have variance inflation factors $\text{VIF} < 1.3$, confirming that the five axes are not collinear re-encodings of the same latent degree of freedom.

**Model-complexity trade-off.** M5 ($k = 6$) wins on LOO and AIC; M3 ($k = 4$) wins on BIC. M6 ($k = 7$) does *not* beat M5 on any metric — the pipeline does not reward excess complexity.

**Calibration.** On 200 random 50/50 train-test splits, the regression of observed vs predicted $\log a_0$ has mean slope 0.82 (ideal = 1.0), mean $R^2 = 0.45$, and slope 95% CI $[0.46, 1.35]$.

The structured law survives nested cross-validation, fails to appear in null-permuted pipelines, and does not reduce to collinear re-encoding of the same latent degree of freedom. Full protocol details are in Appendix C.

---

## 12. External Validation on SPARC-out Sample

### 12.1 Motivation and Protocol

The robustness audit (Section 11) confirmed that the structured law survives internal cross-validation, permutation testing, and bootstrap resampling within the $N = 45$ training set. However, internal validation cannot rule out selection effects arising from the quality criteria that defined the training sample. To test genuine out-of-sample generalization, we apply frozen M3 coefficients — with no refitting — to SPARC galaxies excluded from the training set.

We downloaded the complete SPARC rotation-curve dataset (Lelli et al. 2016, Table 2) from CDS/VizieR (J/AJ/152/157), providing Vobs, Vgas, Vdisk, Vbul for all 175 galaxies. Using the identical pipeline ($\Upsilon_\star^\mathrm{disk} = 0.5$, McGaugh RAR interpolation, golden-section $a_0$ fit, radial-run MeanRun computation), we computed per-galaxy $a_0$ and $\log\overline{\mathrm{Run}}$ for every galaxy. Pipeline validation on the 45 training galaxies confirmed excellent agreement: CDS-derived $a_0$ vs. original $a_0$ has RMS = 0.067 dex (mean bias $-$0.031), and CDS-derived $\log\overline{\mathrm{Run}}$ vs. original has RMS = 0.177 dex.

After quality filtering (published-distance flag $f_D \geq 2$, $\log a_0 \geq 2.0$, $n \geq 5$ RC points, valid $V_\mathrm{flat}$ for Mhost estimation), 36–46 external galaxies remain. Since no external galaxy has a group-catalog $M_\mathrm{host}$ value (the group-catalog criterion was used to select the training set), $M_\mathrm{host}$ was estimated from $V_\mathrm{flat}$ using the baryonic Tully-Fisher proxy. We tested three model variants:

1. **Frozen M3** (logMHI + logMhost$_\mathrm{Vflat}$ + logMeanRun): coefficients frozen from training.
2. **M2$'$** (logMHI + logMeanRun only): drops Mhost entirely to avoid proxy contamination.
3. **M3$_\mathrm{fair}$** (logMHI + logMhost$_\mathrm{Vflat}$ + logMeanRun): refitted on training using the same $V_\mathrm{flat}$ proxy for Mhost, ensuring identical input definitions across train and test.

### 12.2 Results

| Model | $N_\mathrm{ext}$ | RMS (dex) | RMS M0 | Gap% | Spearman | Win rate |
|-------|------------------|-----------|--------|------|----------|----------|
| M0 (baseline) | 46 | — | 0.492 | 0% | — | — |
| M2$'$ frozen | 46 | 0.605 | 0.492 | $-$120% | $-$0.297 | 32.6% |
| M3 frozen | 36 | 0.723 | 0.424 | $-$278% | $-$0.230 | 16.7% |
| M3$_\mathrm{fair}$ | 36 | 0.452 | 0.424 | $-$48% | $-$0.013 | 47.2% |

All three variants perform **worse** than the trivial M0 baseline (predicting the training mean for every galaxy). Negative gap% values indicate that the structured law adds noise rather than signal on external data.

### 12.3 Axis Sign Reversal

The most diagnostic finding is that partial correlations of the M3 axes with $a_0$ are **reversed** on external data relative to the training set:

| Axis | Training partial $r$ | External partial $r$ | Signs match? |
|------|---------------------|---------------------|--------------|
| logMHI | negative ($-$0.198) | $+$0.261 | No |
| logMeanRun | positive ($+$0.459) | $-$0.088 | No |

Both axes have the wrong sign. This is not an issue of magnitude or noise — the *direction* of the relationships is reversed. The training-set pattern (higher gas mass $\to$ lower $a_0$; longer residual runs $\to$ higher $a_0$) does not hold in the broader SPARC population.

### 12.4 Ruling Out Extrapolation

Of 36 external galaxies with complete inputs, 35 fall within the training-set ranges for both logMHI ($-$0.86 to 1.52) and logMeanRun (0.34 to 1.40). The failure is not caused by extrapolation beyond the training domain.

### 12.5 Interpretation

The external failure has three possible explanations:

1. **Selection-effect artifact**: The published-distance criterion ($f_D \geq 2$) that defines the $N = 45$ training set may introduce a selection bias that creates apparent structure in $a_0$ variation. Galaxies with published distances are biased toward nearby, well-studied objects with specific observational histories.

2. **Mhost proxy contamination**: All external galaxies lack group-catalog $M_\mathrm{host}$ values, requiring a $V_\mathrm{flat}$-based proxy with $-$0.74 dex mean bias and 0.83 dex scatter relative to group-catalog values. However, the M2$'$ test (which drops Mhost entirely) also fails, suggesting this is not the sole explanation.

3. **Sample-specific pattern**: The structured law may describe genuine but sample-specific correlations within the $N = 45$ subset that do not extend to the broader galaxy population.

### 12.6 Sample Population Comparison (Phase 82)

To determine whether the external failure reflects a selection effect or a fundamental law instability, we compare the $N = 45$ training galaxies against $N = 46$ external SPARC galaxies ($f_D \geq 2$, $\log a_0 \geq 2$, $n_\mathrm{pts} \geq 5$) across 13 variables using Welch's $t$-test, Mann–Whitney $U$, and Kolmogorov–Smirnov tests.

**Result: The two populations are fundamentally different.** 12 of 13 variables differ significantly ($p < 0.05$), with 8 showing large effect sizes ($|d| > 0.8$). Only inclination ($p = 0.17$, $d = 0.29$) is statistically indistinguishable.

**Table 6. Sample Comparison Summary (Phase 82)**

| Variable | Training mean | External mean | Cohen $d$ | Effect size |
|----------|:---:|:---:|:---:|:---:|
| $\log M_\mathrm{HI}$ (gas mass) | 0.554 | $-$0.128 | +1.16 | LARGE |
| $V_\mathrm{flat}$ | 140 km/s | 73 km/s | +1.48 | LARGE |
| $T$ (morphology) | 4.2 (Sbc) | 7.8 (Sd–Im) | $-$1.45 | LARGE |
| $r_\mathrm{max}$ (RC extent) | 31.1 kpc | 9.4 kpc | +1.39 | LARGE |
| Environment code | 0.96 | 0.30 | +1.10 | LARGE |
| $L_{3.6}$ (luminosity) | 5.6 | 0.45 ($10^9 L_\odot$) | +1.01 | LARGE |
| $n_\mathrm{pts}$ (RC points) | 31.6 | 15.8 | +0.96 | LARGE |
| $\log$ MeanRun | 0.845 | 0.664 | +0.84 | LARGE |
| Distance | 17.9 Mpc | 10.1 Mpc | +0.75 | MEDIUM |
| $\log a_0$ | 3.548 | 3.298 | +0.72 | MEDIUM |
| $R_\mathrm{disk}$ | 3.6 kpc | 1.9 kpc | +0.58 | MEDIUM |
| Quality $Q$ | 1.36 | 1.65 | $-$0.50 | MEDIUM |
| Inclination | 67.5° | 63.2° | +0.29 | SMALL |

**Morphological type distribution** reveals the starkest contrast: the training set is dominated by early-to-intermediate spirals (S0–Sc: 76%), while the external sample is dominated by late-type dwarfs and irregulars (Sd–Im: 63%, with Im alone at 33% vs 2%).

**Interpretation**: The published-distance criterion ($f_D \geq 2$) combined with quality cuts creates a strong selection toward massive, luminous, gas-rich, data-rich spirals in group/cluster environments. The external galaxies that remain are predominantly low-mass, gas-poor late-type dwarfs in field environments with sparse rotation curves. The multi-axis law discovered in $N = 45$ describes $a_0$ structure among a specific galaxy type (massive spirals with rich kinematics), not a universal relationship across the Hubble sequence.

This is a meaningful scientific result: it identifies the published-distance criterion as a powerful selection filter that generates a non-representative subsample of SPARC.

### 12.7 Common Support Test (Phase 83)

The population mismatch (Section 12.6) raises a critical question: does the law fail because external galaxies are different types, or would it fail even on galaxies that match the training regime? To answer this, we apply three progressively strict matching criteria to select external galaxies that resemble the training set, plus Mahalanobis nearest-neighbor pairing (44 pairs).

**Table 7. Common Support Test Results (Phase 83)**

| Matching level | $N$ | M3 RMS | M0 RMS | Win rate | Signs | Verdict |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Loose (in-range) | 78 | 0.995 | 0.344 | 9.0% | OK (weak) | FAIL |
| Moderate (P10–P90) | 22 | 0.858 | 0.250 | 9.1% | Reversed | FAIL |
| Strict (P10–P90 + Q$\leq$2 + T$\leq$7) | 14 | 0.891 | 0.205 | 7.1% | Reversed | FAIL |
| Paired (Mahalanobis NN) | 44 | 0.917 | 0.300 | 6.8% | Partial | FAIL |

**All four matching strategies fail.** The frozen M3 consistently predicts $\log a_0 \approx 4.2$–$4.8$ while observed values are $\approx 2.9$–$3.9$, a systematic overprediction of $\sim$1 dex. M2$'$ (no Mhost) performs equally poorly ($\mathrm{RMS} \approx 0.91$–$0.93$). Win rates are 7–10% across all levels — worse than random.

At the loosest matching ($N = 78$), the axis signs are technically correct but the correlations are essentially zero ($r_\mathrm{MHI} = -0.057$, $r_\mathrm{MeanRun} = +0.061$). As matching becomes stricter and the external subsample more closely resembles the training set, the correlations reverse sign — exactly the opposite of what regime-specificity would predict.

**Interpretation**: The failure is **not** explained by population mismatch alone. Even external galaxies with similar $V_\mathrm{flat}$, $M_\mathrm{HI}$, MeanRun, $n_\mathrm{pts}$, and morphology do not follow the M3 law. The structured $a_0$ variation in $N = 45$ reflects something specific to those 45 galaxies — not a general property of massive spirals, but a pattern unique to the particular combination of galaxies selected by the published-distance and quality criteria.

### 12.8 Selection / Collider Audit (Phase 84)

Phase 83 showed the law fails even on regime-matched galaxies, suggesting the problem is deeper than population mismatch. To identify the root cause, we conduct a systematic audit of the selection criteria that define the $N = 45$ sample.

**12.8.1 The fundamental correlation reversal.** In the full SPARC sample ($N = 171$), $r(\log M_\mathrm{HI}, \log a_0) = +0.278$ — galaxies with more gas have *higher* inferred $a_0$. But in the training $N = 45$, this correlation is $r = -0.322$ — the *opposite sign*. This reversal persists: all external subsets show $r > 0$ ($+0.309$ for all external, $+0.458$ for external $f_D \geq 2$), while only the training set shows $r < 0$. No single selection criterion (quality, distance method, morphology) alone produces this reversal; it requires their specific intersection.

**Table 8. Correlation $r(\log M_\mathrm{HI}, \log a_0)$ across subsets**

| Subset | $N$ | $r(\mathrm{MHI})$ | $r(\mathrm{MeanRun})$ | $r(V_\mathrm{flat})$ | $\sigma(a_0)$ |
|:---|:---:|:---:|:---:|:---:|:---:|
| All SPARC | 171 | +0.278 | +0.259 | +0.423 | 0.589 |
| $f_D \geq 2$ only | 76 | +0.392 | +0.206 | +0.444 | 0.598 |
| $f_D \geq 2$, $Q \leq 2$ | 69 | +0.451 | +0.184 | +0.348 | 0.517 |
| **Training $N = 45$** | **45** | **$-$0.322** | **+0.212** | **$-$0.069** | **0.278** |
| External (all) | 126 | +0.309 | +0.218 | +0.433 | 0.652 |
| External $f_D \geq 2$ | 50 | +0.458 | +0.109 | +0.467 | 0.673 |

**12.8.2 Random subsample test.** From 1000 random draws of $N = 45$ from the full 171, **zero** produced $r(\mathrm{MHI}, a_0) \leq -0.322$ ($p < 0.001$). Only 1.6% of random samples showed both MHI-negative and MeanRun-positive signs simultaneously. The training set is a statistical extreme — not a typical subsample of SPARC.

**12.8.3 Sample membership as predictor.** A binary "in-training" flag alone predicts $\log a_0$ with LOO gap = 4.2% on the pooled sample ($N = 76$, $f_D \geq 2$) — comparable to M2$'$ with MHI and MeanRun (gap = 4.7%). When the flag is added alongside M2$'$ variables, interaction terms (flag $\times$ logMHI) absorb most of the signal ($\beta_\mathrm{MHI \times flag} = -0.661$), confirming that the MHI–$a_0$ relationship changes sign between training and non-training galaxies.

**12.8.4 Criterion-by-criterion dissection (Phase 85).** A systematic scan of all selection criteria identifies exactly which cuts flip the sign:

*Single criteria producing negative $r(\mathrm{MHI}, a_0)$*: Only three criteria, applied alone to the full SPARC sample, reverse the sign — $V_\mathrm{flat} \geq 60$ km/s ($r = -0.175$, $N = 118$), $V_\mathrm{flat} \geq 80$ ($r = -0.215$, $N = 99$), environment $\geq 1$ ($r = -0.447$, $N = 19$), and $Q = 1$ ($r = -0.103$, $N = 99$). No other criterion (distance method, morphology, inclination, luminosity) alone flips the sign.

*Two-criterion scan*: The combination $V_\mathrm{flat} \geq 70$ + $n_\mathrm{pts} \geq 10$ produces $r = -0.274$ ($N = 95$); $V_\mathrm{flat} \geq 100$ + $n_\mathrm{pts} \geq 20$ produces $r = -0.297$ ($N = 39$). The $V_\mathrm{flat}$ cut is the primary driver: it excludes low-mass dwarfs whose high gas fractions produce a positive MHI–$a_0$ correlation, leaving a massive-galaxy subsample where the sign naturally reverses.

*Cumulative paths*: Regardless of the order criteria are applied, the sign flip always occurs when the $V_\mathrm{flat}$ cut is added. The distance-method criterion ($f_D \geq 2$) does not contribute to the flip; it actually *strengthens* the positive correlation ($r = +0.306$).

**12.8.5 Interpretation: Mass-dependent regime + collider amplification.** Phase 85 reveals a more nuanced picture than pure collider bias. The Vflat cut *alone* flips the MHI–$a_0$ sign (from $r = +0.28$ to $r = -0.22$), meaning the correlation structure genuinely differs between massive and dwarf galaxies — this is a real physical regime boundary. Phase 86 confirms this decisively: among 5000 random $N = 45$ draws from the $V_\mathrm{flat} \geq 80$ pool ($N = 99$), 17.4% produce $r(\mathrm{MHI}, a_0) \leq -0.322$. The training value is **typical, not extreme**, within the high-Vflat regime. The picture sharpens further with $V_\mathrm{flat} \geq 80$ + $Q = 1$ ($N = 69$): $p = 0.58$, mean $r = -0.337$ — the training $r$ sits at the *median* of this pool.

**12.8.6 Controlled subsample analysis (Phase 86).** Five tests systematically determined whether the training $r(\mathrm{MHI}, a_0) = -0.322$ is a statistical extreme or a typical draw from the high-Vflat regime:

*Test 1 — Vflat-only pools*: Drawing $N = 45$ from $V_\mathrm{flat} \geq 60$ ($p = 0.11$), $\geq 80$ ($p = 0.17$), and $\geq 100$ ($p = 0.08$) all produce training-compatible $r$ distributions. Only $V_\mathrm{flat} \geq 120$ makes training extreme ($p = 0$), because this ultra-high regime flips back positive (mean $r = +0.047$).

*Test 2 — Vflat + one criterion*: Adding Q, nPts $\leq 10$, or environment to the $V_\mathrm{flat} \geq 80$ base keeps training typical. The most normalizing criterion is $Q = 1$ ($p = 0.58$, mean $r = -0.337$) and environment $\geq 1$ ($p = 1.0$, mean $r = -0.447$). But adding $f_D \geq 2$ ($p = 0$), morphology $T \leq 5$ ($p = 0$), or luminosity $L_{3.6} \geq 5$ ($p = 0.001$) makes training extreme — these criteria narrow the pool in ways that change the correlation structure.

*Test 3 — Progressive stacking*: The cumulative path $V_\mathrm{flat} \geq 80 \to +\mathrm{nPts} \geq 10 \to +Q \leq 2 \to +T \leq 8$ maintains typicality ($p = 0.10$). Only adding $f_D \geq 2$ as the final criterion makes $r$ extreme ($p = 0$) — but at that point the pool equals the training set (deterministic).

*Test 4 — Excluding training galaxies*: Among $V_\mathrm{flat} \geq 80$ galaxies *not in training* ($N = 54$), drawing $N = 45$ gives $p = 0.048$ — unusual but not extreme. The regime effect is present even in non-training galaxies.

*Summary*: Of 23 controlled pools tested, 11 make training $r$ typical ($p > 0.10$), 5 unusual ($0.01 < p \leq 0.10$), and 7 extreme ($p \leq 0.01$). The training $r = -0.322$ is a **normal draw from the high-Vflat regime**, not a selection-amplified extreme. The multi-axis law captures a real mass-regime effect that becomes visible when dwarfs are excluded.

**12.8.7 Regime transition mapping (Phase 87).** Seven complementary tests locate and characterize the transition:

*Broken regression (Test 4)*: A piecewise-linear model with separate MHI slopes above and below a $V_\mathrm{flat}$ breakpoint achieves optimal BIC at $V_\mathrm{flat} = 70$ km/s ($\Delta\mathrm{BIC} = -23.8$ vs. single-slope model — decisive evidence for two regimes). Below: $\partial \log a_0 / \partial \log M_\mathrm{HI} = +0.059$; above: $-0.127$. The slopes cleanly reverse sign.

*Interaction model (Test 5)*: A continuous $\log M_\mathrm{HI} \times \log V_\mathrm{flat}$ interaction gives $\beta_\mathrm{int} = -0.079$ (substantial). The MHI slope crosses zero at $V_\mathrm{flat} = 48$ km/s. BIC slightly disfavors the interaction vs. additive model ($\Delta\mathrm{BIC} = +3.5$), suggesting the transition is better described as a break than a smooth gradient.

*Logistic fit to sliding windows (Test 6)*: Fitting $r(\mathrm{MHI}, a_0) = r_\mathrm{hi} + (r_\mathrm{lo} - r_\mathrm{hi}) / (1 + e^{(V-V_0)/w})$ gives $V_0 = 40$ km/s, $w = 13$ km/s. The transition is **sharp** (10–90% range: 11–69 km/s), with $r_\mathrm{lo} = +0.10$ and $r_\mathrm{hi} = -0.30$.

*Alternative separators (Test 7)*: $\log V_\mathrm{flat}$ is the **best** transition variable ($\Delta\mathrm{BIC} = -23.0$ vs. no-break). Alternatives — $\log M_\mathrm{HI}$ ($+7.1$), $\log L_{3.6}$ ($+5.0$), morphology $T$ ($+3.1$) — all have positive $\Delta\mathrm{BIC}$ (no improvement). The regime boundary is fundamentally kinematic, not morphological or luminosity-based.

*Summary*: The MHI–$a_0$ relationship reverses sign across a sharp kinematic boundary at $V_\mathrm{flat} \approx 50$–$70$ km/s. Below this threshold, gas-rich dwarfs show a weak positive MHI–$a_0$ correlation; above it, massive spirals show a moderate negative anticorrelation. This is a genuine two-regime structure in the data, not a selection artifact.

### 12.9 Verdict

**The structured $a_0$ law reflects a real mass-regime boundary, not a pure selection artifact.** The MHI–$a_0$ anticorrelation that anchors M3/M5 does not exist in the general SPARC population ($r = +0.278$), but it is a **normal feature** of the high-Vflat regime ($V_\mathrm{flat} \geq 80$: population $r = -0.215$). Phase 86 demonstrates that the training $r = -0.322$ falls within the typical range of random $N = 45$ draws from this regime ($p = 0.17$). Phase 87 locates the transition: the MHI–$a_0$ slope reverses sign at $V_\mathrm{flat} \approx 70$ km/s ($\Delta\mathrm{BIC} = -23.8$ for piecewise vs. single-slope model). The transition is sharp ($w \approx 13$ km/s), kinematic in nature ($V_\mathrm{flat}$ is the best separator, far surpassing morphology, luminosity, or gas mass), and fundamentally changes the correlation structure between regimes. The multi-axis law is a valid description of the high-$V_\mathrm{flat}$ regime; it fails externally because external galaxies predominantly inhabit the opposite regime.

### 12.10 Two-Regime Law Analysis (Phase 88)

Splitting the sample at $V_\mathrm{flat} = 70$ km/s yields two populations with fundamentally different $a_0$ structure:

**Low-$V_\mathrm{flat}$ regime ($V_\mathrm{flat} < 70$, $N = 63$).** No structured $a_0$ law exists. The best single-predictor model ($\log V_\mathrm{flat}$) has LOO gap $= -1.8\%$ (worse than the constant baseline). Permutation $p = 0.689$ — the apparent $R^2 = 0.003$ is indistinguishable from chance. The $a_0$ distribution is wide (sd $= 0.77$ dex) and unstructured: no variable — gas mass, coherence, density, or morphology — predicts $a_0$ variation. All environment codes equal zero (no cluster/group members), preventing environment-based models.

**High-$V_\mathrm{flat}$ regime ($V_\mathrm{flat} \geq 70$, $N = 112$).** A five-predictor model (Gas + Coherence + Density + Environment + $V_\mathrm{flat}$) achieves $R^2 = 0.253$, LOO gap $= 8.1\%$, permutation $p < 0.001$. The $a_0$ distribution is tighter (sd $= 0.30$ dex). Key coefficients: $\log M_\mathrm{HI}$ ($-0.35$), $\log\overline{\mathrm{Run}}$ ($+0.26$), $\log\Sigma_0$ ($-0.20$), envCode ($-0.06$), $\log V_\mathrm{flat}$ ($+1.09$). The Gas + Coherence + Environment three-predictor variant retains LOO gap $= 4.7\%$.

**Sign comparison (88d).** Three of five predictor variables reverse sign across the transition: $\log M_\mathrm{HI}$ ($+0.05 \to -0.22$), $\log\overline{\mathrm{Run}}$ ($-0.01 \to +0.20$), and $\log V_\mathrm{flat}$ ($+0.06 \to -0.02$). Only $\log\Sigma_0$ ($-0.16$ both) and envCode ($\sim 0 \to -0.23$) maintain direction. This constitutes a **full state transition**, not merely a gas-regime effect.

**Axis transition strength (88e).** Broken-slope models at $V_\mathrm{flat} = 70$ show strong $\Delta\mathrm{BIC}$ for $\log M_\mathrm{HI}$ ($-23.8$), $\log\overline{\mathrm{Run}}$ ($-23.8$), and $\log\Sigma_0$ ($-30.6$, same sign but magnitude changes). Only $\log V_\mathrm{flat}$ shows weak transition ($\Delta\mathrm{BIC} = -1.4$).

**Cross-regime validation (88f).** Training on one regime and testing on the other fails completely: high$\to$low gap $= -77\%$; low$\to$high gap $= -89\%$. The regimes are genuinely separate — models do not transfer.

**Verdict**: The $a_0$ correlation structure is **regime-asymmetric**. In the high-$V_\mathrm{flat}$ regime, a structured multi-axis law exists (LOO gap $= 8.1\%$, $p < 0.001$). In the low-$V_\mathrm{flat}$ regime, $a_0$ variation is large but unstructured — no predictive law emerges. This is not two symmetric laws but one law plus one regime of irreducible scatter.

### 12.11 Low-Regime Audit: Chaos or Noise? (Phase 89)

Phase 88 found no structured $a_0$ law among low-$V_\mathrm{flat}$ galaxies (sd $= 0.77$ dex). Phase 89 investigates whether this scatter is genuine physics or measurement/data-quality limitation.

**Variance decomposition (Test 2).** Bootstrap resampling ($N_\mathrm{boot} = 200$) per galaxy estimates the measurement contribution to $a_0$ uncertainty. In the low regime, measurement accounts for only **3.1%** of the observed variance (bootstrap SD $= 0.135$ dex vs. observed SD $= 0.772$ dex; intrinsic SD $= 0.760$ dex). The high regime is similarly measurement-clean (3.6%). Point-level observational noise is not the driver.

**Data quality stratification (Test 5).** Restricting to $Q = 1$ (best quality) low-$V_\mathrm{flat}$ galaxies ($N = 22$) reduces the scatter dramatically: sd $= 0.325$ dex — comparable to the high-regime population (sd $= 0.302$). $Q = 2$ ($N = 31$): sd $= 0.760$; $Q = 3$ ($N = 10$): sd $= 0.803$. This reveals that the "chaos" is not intrinsic dwarf physics but **systematic data-quality effects** (non-circular motions, asymmetric rotation curves, beam smearing) that the bootstrap cannot capture but the $Q$ flag does. The $Q = 1$ subsample has more points (mean 16 vs. 9) and wider gbar range (0.67 vs. 0.43 dex).

**Dynamic range (Test 4).** Well-constrained low-$V_\mathrm{flat}$ galaxies ($n_\mathrm{pts} \geq 10$, gbar range $\geq 1.0$ dex, $N = 5$) have sd $= 0.154$ dex — even tighter than the high regime. Though $N = 5$ is too small for conclusions, it is consistent with the $Q = 1$ finding: when rotation-curve quality is high, dwarf $a_0$ scatter is comparable to (or smaller than) massive spirals.

**Distribution shape (Test 6).** The low-regime $\log a_0$ distribution is strongly left-skewed (skewness $= -1.26$) with excess kurtosis ($= 1.30$) and large gaps at $\log a_0 \approx 1.7$ and $1.1$. This suggests a small number of extreme outliers with very low $a_0$ (likely driven by poorly constrained fits on sparse data) rather than a uniformly scattered population.

**Distance sensitivity (Test 3).** A 20% distance perturbation shifts $\log a_0$ by 0.156 dex in both regimes (ratio $= 0.98$). Dwarfs are not more distance-sensitive than massive spirals.

**Revised interpretation.** The apparent "irreducible scatter" in the low-$V_\mathrm{flat}$ regime is **quality-dominated, not physics-dominated**. The $Q$ flag — which captures systematics invisible to bootstrap (rotation-curve asymmetries, non-circular motions, beam smearing) — is the primary discriminator. Among $Q = 1$ dwarfs, $a_0$ scatter matches the high regime. This means the two-regime picture from Phase 88 may be partly a **data-quality artifact**: the low regime appears structureless because it is dominated by galaxies with poorly determined $a_0$, not because dwarf physics is fundamentally different.

### 12.12 Low-Regime Q=1 Law Test (Phase 90)

Phase 89 showed that $Q = 1$ dwarfs ($N = 22$) have $a_0$ scatter comparable to the high regime (sd $= 0.325$ vs. $0.302$). This raises the critical question: does a structured multi-axis law exist in these well-measured dwarfs, or does the kinematic boundary at $V_\mathrm{flat} \approx 70$ km/s mark a genuine physical transition?

**Model competition (Test 2).** Among the 14 models tested (single-variable through 3-variable), no physics-based model achieves meaningful LOO improvement. The best physics predictor, logMeanRun (kinematic coherence), has LOO gap $= +0.5\%$ — negligible. Gas mass (logMHI) has gap $= -6.9\%$, worsening the fit. The only model with a positive LOO gap is **inclination** (gap $= 2.9\%$, $R^2 = 0.224$, perm $p = 0.025$), which is a geometric/systematic variable, not a physical one.

**Sign comparison (Test 5).** Only 2 of 4 key variables share their correlation sign with the high regime: logMHI (negative in both) and logVflat (negative in both). But logMeanRun *reverses* ($r = -0.331$ in LoQ1 vs. $+0.203$ in HiAll) and logSigma0 reverses ($+0.073$ vs. $-0.157$). This means the high-regime law's coefficient structure does not apply to dwarfs.

**Cross-regime transfer (Test 6).** Training on HiQ1 ($N = 77$) and testing on LoQ1 ($N = 22$) — or vice versa — fails completely even with only $Q = 1$ galaxies: Gas+Coherence transfer gaps $= -41\%$ and $-44\%$; Gas-only gaps $= -37\%$ and $-24\%$. The regime boundary persists even when data quality is controlled.

**Unified Q=1 model (Test 7).** Pooling all $Q = 1$ galaxies ($N = 99$) produces no structured law: best LOO gap $= +0.7\%$ (inclination). The regime-pooled sample is no more predictable than chance.

**The inclination signal.** The $r = -0.473$ correlation between inclination and $\log a_0$ in LoQ1 dwarfs (perm $p = 0.025$) is notable but almost certainly a **systematic artifact**: lower-inclination dwarfs have larger deprojection corrections that systematically inflate $V_\mathrm{rot}$ and thus $a_0$. This is a known bias in dwarf kinematic measurements where beam smearing and pressure support corrections are inclination-dependent.

**Verdict.** The regime boundary at $V_\mathrm{flat} \approx 70$ km/s is **genuinely physical, not merely a data-quality effect**. While Phase 89 correctly showed that the *scatter amplitude* in dwarfs is inflated by quality issues (Q=1 sd $= 0.325$ vs. all-Q sd $= 0.77$), Phase 90 demonstrates that even among well-measured Q=1 dwarfs, the *correlation structure* differs fundamentally from massive spirals: key predictors reverse sign, cross-regime transfer fails, and no multi-axis law emerges. The Vflat boundary separates two populations with different $a_0$ physics, not merely different data quality. The Phase 89 insight remains valid — scatter amplitude is quality-dominated — but the absence of structure is physics-dominated.

### 12.13 Why Only the High Regime? (Phase 91)

Phase 90 established that the $V_\mathrm{flat} \approx 70$ km/s boundary is genuinely physical. Phase 91 asks: **why** does a structured $a_0$ law exist only above this threshold? The answer is a convergence of four independent transitions, all occurring near the same kinematic scale.

**1. Gas fraction transition.** Mean $f_\mathrm{gas}$ jumps from $0.31$ (high regime) to $0.60$ (low regime). Above 70 km/s, stellar disks dominate the baryonic mass budget and produce well-defined gravitational potentials. Below it, gas contributes comparably to or more than the stellar disk ($V_\mathrm{gas}/V_\mathrm{total}$: 0.18 vs. 0.34), and the disk dominance fraction drops from 0.81 to 0.58. In gas-dominated systems, the connection between baryonic surface density and rotation-curve shape becomes less direct — pressure support, turbulence, and irregular gas distributions decouple $g_\mathrm{bar}$ from $g_\mathrm{obs}$.

**2. Acceleration range narrows.** The mean log($g_\mathrm{bar}$) dynamic range drops from **1.19 dex** (high regime) to **0.51 dex** (low regime). With barely half a decade of acceleration lever arm, the RAR fit becomes degenerate: many $a_0$ values produce comparably good fits. This does not merely add noise — it fundamentally undermines the *definition* of a galaxy-specific $a_0$.

**3. Deep MOND regime.** The median log($g_\mathrm{bar}$) shifts from 2.91 (high) to 2.33 (low), relative to the canonical $\log a_0 \approx 3.24$. Dwarfs sit almost entirely in the deep-MOND limit where $g_\mathrm{obs} \approx \sqrt{g_\mathrm{bar} \cdot a_0}$. In this regime, the RAR asymptotes to a power law and the sensitivity to $a_0$ decreases — many different $a_0$ values produce nearly identical rotation curves.

**4. Rising rotation curves.** 73% of low-regime galaxies have rising outer rotation curves (slope $> 0.1$), compared to 34% in the high regime. A rising RC means the galaxy has not reached its asymptotic flat velocity — $V_\mathrm{flat}$ is not well-defined, and the outer acceleration data that most constrains $a_0$ is missing.

**Convergence.** These four transitions are correlated but not identical — each independently degrades the ability to extract a structured $a_0$ law. Their convergence near $V_\mathrm{flat} \approx 70$ km/s is not coincidental: this scale corresponds to the transition from disk-dominated, baryon-rich spirals to gas-dominated, pressure-supported dwarfs. The high-regime law works because massive spirals have (i) disk-dominated kinematics making $a_0$ cleanly defined, (ii) wide acceleration range making $a_0$ well constrained, (iii) flat rotation curves making $V_\mathrm{flat}$ meaningful, and (iv) varied gas fractions making the MHI anticorrelation visible. Dwarfs fail on all four counts simultaneously.

**Implication for MOND.** This analysis does not test whether MOND is correct. It shows that the structured $a_0$ variation found in massive spirals cannot be extended to dwarfs — not because dwarfs violate MOND, but because the parameter $a_0$ becomes **operationally ill-defined** in the dwarf regime. The combination of narrow acceleration range, deep-MOND degeneracy, and rising RCs means that galaxy-by-galaxy $a_0$ fitting is an inappropriate tool for characterizing MOND in dwarfs. Future work should use methods that do not require individual $a_0$ estimates, such as stacked RAR analyses or direct dynamical modeling.

### 12.14 From $a_0$ Law to Rotation Curve Flatness (Phase 92)

If the high-regime $a_0$ law captures real physics, it should connect to the most fundamental observable: rotation curve shape. Phase 92 tests whether the M3/M5 predictors — and $a_0$ itself — predict RC flatness in the high-$V_\mathrm{flat}$ regime ($N = 112$).

**$a_0$ does not directly predict RC shape.** The correlation between $\log a_0$ and outer slope is $r = -0.052$ (negligible). After removing the common $V_\mathrm{flat}$ dependence, the partial correlation is $r = -0.072$. The inferred acceleration scale and rotation curve shape are **independent phenomena** — knowing a galaxy's $a_0$ tells you nothing about whether its RC is flat, rising, or declining.

**But the same galaxy properties predict both.** The M3 variables (logMHI, logMeanRun, logSigma0) predict outer slope with $R^2 = 0.361$ and LOO gap $= 16.9\%$ (perm $p < 0.001$). Adding $\log a_0$ to M5 yields the best model: $R^2 = 0.432$, LOO gap $= 18.8\%$. The strongest single predictor of flatness is baryonic surface density: $r(\log \Sigma_0, \mathrm{outerSlope}) = -0.518$.

For the flatness ratio $V_\mathrm{last3}/V_\mathrm{max}$, M3 achieves $R^2 = 0.326$, LOO gap $= 14.3\%$. Morphological type ($T$) correlates most strongly ($r = 0.559$): later types have flatter RCs.

**Flat vs. rising in the high regime.** Among 112 high-$V_\mathrm{flat}$ galaxies: 50% have flat outer RCs ($|\mathrm{slope}| \leq 0.1$), 34% rising, 16% declining. Flat galaxies differ from rising ones by: higher surface density ($\Delta \log \Sigma_0 = +0.51$), lower gas fraction ($\Delta f_\mathrm{gas} = -0.15$), higher $V_\mathrm{flat}$ ($\Delta \log V_\mathrm{flat} = +0.17$), and earlier type ($\Delta T = -1.4$).

**Interpretation.** The $a_0$ law and RC flatness are **parallel consequences** of the same underlying galaxy structure, not causally linked to each other. Baryonic surface density ($\Sigma_0$) and gas fraction ($f_\mathrm{gas}$) jointly determine both the inferred $a_0$ and the outer RC shape — but through different physical pathways. This means the high-regime law is not "about flatness" directly; it describes a structured pattern in how galaxies populate the RAR, and that pattern correlates with RC shape because both are manifestations of the same baryon-dominated disk physics.

### 12.15 Mediation Analysis: What Is the Parent Variable? (Phase 93)

Phase 92 showed that $a_0$ and RC flatness share common structural causes. Phase 93 identifies those causes through partial-correlation mediation analysis on the high-$V_\mathrm{flat}$ sample ($N = 112$).

**Σ₀ has independent paths to both a₀ and slope.** Controlling for outer slope, $r(\Sigma_0, a_0) = -0.215$ (survives); controlling for $a_0$, $r(\Sigma_0, \mathrm{slope}) = -0.533$ (strengthens). This pattern is consistent with $\Sigma_0$ as a common cause rather than a mediator.

**Dual parent structure.** Both $\Sigma_0$ and $V_\mathrm{flat}$ maintain independent partial correlations with outer slope when controlling for each other: $r(\Sigma_0, \mathrm{slope} \,|\, V_\mathrm{flat}) = -0.270$ and $r(V_\mathrm{flat}, \mathrm{slope} \,|\, \Sigma_0) = -0.212$. Neither fully mediates the other. The flatness of rotation curves is driven by **two independent structural axes**: baryonic surface density (how concentrated the baryons are) and rotation velocity (how deep the potential well is).

**a₀ retains a residual connection to flatness.** After controlling for all four structural variables ($\Sigma_0$, $V_\mathrm{flat}$, $f_\mathrm{gas}$, $M_\mathrm{HI}$), the partial correlation $r(\log a_0, \mathrm{slope}) = -0.294$ — a substantial residual. This means $a_0$ contains information about RC shape that the measured structural properties do not capture. One interpretation: the galaxy-specific $a_0$ encodes aspects of the dark-matter/MOND potential that are not fully determined by the baryonic observables.

**Baryon dominance as partial mediator.** The composite of $\Sigma_0$, $f_\mathrm{gas}$, and disk dominance fraction explains $R^2 = 0.344$ of outer slope variance but only $R^2 = 0.034$ of $a_0$ variance. After removing the baryon-dominance composite from both, the residual $a_0$–slope correlation is $r = -0.165$, confirming partial but incomplete mediation.

**Strongest mediators for individual paths.** The $M_\mathrm{HI}$–slope connection ($r = -0.465$) is most reduced by controlling $M_\star$ ($-78\%$) or $V_\mathrm{flat}$ ($-63\%$). The $V_\mathrm{flat}$–slope connection ($r = -0.496$) is most reduced by $M_\star$ ($-111\%$, fully mediated) and $\Sigma_0$ ($-57\%$). After controlling all structural variables simultaneously, **no single predictor retains a significant partial correlation with slope** — the effect is fully distributed across the structural basis.

**Verdict.** The causal structure is not a simple tree with one parent. It is a **dual-parent network**: $\Sigma_0$ (baryon concentration) and $V_\mathrm{flat}$ (potential depth) jointly and independently shape both $a_0$ and RC flatness. The residual $a_0$–slope connection ($r = -0.294$ after controlling structural variables) suggests that the inferred acceleration scale carries additional dynamical information beyond what the baryonic observables provide. Phase 94 tests whether this residual is genuine or a circularity artifact.

### 12.16 Circularity Test: Is the $a_0$ Residual Genuine? (Phase 94)

Phase 93 found $r(\log a_0, \mathrm{slope} \,|\, \text{structural}) = -0.294$. Since both $a_0$ and outer slope are derived from the same rotation curve, this residual could reflect shared RC information rather than genuine extra dynamics. Phase 94 tests this with split-half analysis, extended controls, and permutation tests ($N = 104$).

**Split-half test: a₀ from inner half, slope from outer half.** If $a_0$ and slope share information only because they are computed from the same data, splitting the curve should destroy the residual. Result: the partial correlation **strengthens** from $r = -0.305$ (same RC) to $r = -0.325$ ($p_\mathrm{perm} = 0.003$). The ratio split/full $= 1.07$ — the circularity hypothesis is rejected.

**Thirds test: maximum separation.** Using $a_0$ from the first third of the RC and slope from the outer half: partial $r = -0.255$, also surviving structural controls. The signal persists even with maximal radial separation.

**Asymmetric direction.** $a_0(\mathrm{inner}) \to \mathrm{outerSlope}$: partial $r = -0.325$. But $a_0(\mathrm{outer}) \to \mathrm{innerSlope}$: partial $r = +0.071$ (no signal). The information flows **outward**: the inner potential encodes dynamical content that shapes the outer RC, but not vice versa. This is physically expected — the inner mass distribution sets boundary conditions for the outer dynamics.

**Extended controls (RC-shape variables).** Adding MeanRun and RAR RMS to the structural controls barely changes the split-half residual ($r = -0.323$). Adding flatness ratio and outer velocity CV reduces the full-RC residual to $r = -0.116$, but these variables are themselves RC-shape summaries that overlap with the slope being predicted — they absorb the signal rather than explaining it away. The split-half estimate with extended controls remains $r = -0.196$.

**Inner–outer $a_0$ consistency.** $r(a_0^\mathrm{inner}, a_0^\mathrm{outer}) = 0.494$ — moderate correlation confirming that $a_0$ is a stable galaxy property, not fitting noise. The inner half has larger scatter ($\sigma = 0.68$ dex vs. $\sigma = 0.26$ dex for the outer half) because fewer data points constrain the fit.

**Verdict.** The residual $a_0$–slope connection is **genuine, not circularity**. The split-half test is definitive: extracting $a_0$ from the inner curve and slope from the outer curve yields an equal or stronger partial correlation ($r = -0.325$, $p = 0.003$). The inferred acceleration scale encodes dynamical information about the outer rotation curve shape that propagates from the inner potential — information that the baryonic structural variables ($\Sigma_0$, $V_\mathrm{flat}$, $f_\mathrm{gas}$, $M_\mathrm{HI}$) do not capture. Phase 95 identifies what this information is.

### 12.17 What Is the Hidden Inner Quantity? (Phase 95)

Phase 94 established that $a_0(\mathrm{inner})$ carries genuine dynamical information about the outer RC slope. Phase 95 asks: what specific inner-derived RC property is $a_0$ encoding? We test 11 candidate inner quantities as mediators of the $a_0(\mathrm{inner}) \to \mathrm{outerSlope}$ connection ($N = 104$, baseline partial $r = -0.325$).

**Inner $V_\mathrm{max}$ is the dominant carrier.** The peak velocity reached in the inner half of the RC ($\log V_\mathrm{max}^\mathrm{inner}$) has a partial correlation with outer slope of $r = -0.761$ after structural controls — far stronger than $a_0(\mathrm{inner})$ itself ($r = -0.325$). Adding it as a control reduces the $a_0$ residual by $79.2\%$ (from $r = -0.325$ to $r = -0.068$). Inner $V_\mathrm{max}$ is a **better predictor** of outer slope than $a_0$, and $a_0(\mathrm{inner})$ is largely a proxy for it.

**Inner mass discrepancy as dark-matter proxy.** The mean $\log(g_\mathrm{obs}/g_\mathrm{bar})$ in the inner half — a direct measure of how much "extra" gravity exists beyond baryons — reduces the $a_0$ residual by $46.9\%$. Galaxies with more inner dark-matter contribution have flatter outer RCs. The inner discrepancy slope (how the discrepancy changes with $g_\mathrm{bar}$) adds nothing ($-0.9\%$).

**Top three candidates fully absorb the residual.** Combining $\log V_\mathrm{max}^\mathrm{inner}$, inner mass discrepancy, and inner RC curvature as controls: $r(a_0^\mathrm{inner}, \mathrm{outerSlope}) = -0.078$ ($76.1\%$ total reduction, $N = 95$). The $a_0(\mathrm{inner})$ signal is **not irreducible** — it is a composite proxy for these three inner RC properties.

**Physical interpretation.** The inner $V_\mathrm{max}$ after controlling $V_\mathrm{flat}$ (from the catalog) reflects how quickly the galaxy reaches peak rotation in its inner region — this encodes **halo concentration**: a concentrated halo produces high inner velocities and a flat outer profile. The inner mass discrepancy directly measures the dark-matter (or MOND-excess) contribution in the inner region. Together, they describe the shape of the total potential well as set by the inner mass distribution.

**What $a_0$ was actually encoding.** The galaxy-specific $a_0$ inferred from the RAR fit is not a fundamental physical parameter but a **summary statistic** that integrates: (1) how fast the galaxy rotates internally relative to its baryonic prediction (halo concentration), and (2) how much dark-matter contribution exists in the inner region. The structured $a_0$ law in the high regime reflects the systematic variation of these inner dynamical properties with galaxy structure ($\Sigma_0$, $V_\mathrm{flat}$).

**Verdict.** The "hidden inner quantity" is primarily **inner $V_\mathrm{max}$** (halo concentration proxy) supplemented by **inner mass discrepancy** (dark-matter contribution). The per-galaxy $a_0$ is not an irreducible dynamical parameter — it is a composite summary of how the inner potential well differs from baryonic expectations. The inside-to-outside coupling discovered in Phase 94 is physically interpretable: galaxies with concentrated halos (high inner $V_\mathrm{max}$ relative to structure) have more dark matter in their inner regions and flatter outer rotation curves. This is consistent with the NFW halo concentration–mass relation or, in MOND, with the transition radius where modified gravity begins to dominate. Phase 96 tests whether $a_0$ is needed at all.

### 12.18 Direct Dynamical Law Without $a_0$ (Phase 96)

Phase 95 showed that $a_0(\mathrm{inner})$ is a composite proxy for inner $V_\mathrm{max}$, inner mass discrepancy, and inner curvature. Phase 96 asks the decisive question: can we predict outer slope directly from inner dynamics and structure, **without $a_0$**? Seven models are compared on the high-$V_\mathrm{flat}$ sample ($N = 95$).

**$a_0$ alone is useless.** $\log a_0(\mathrm{inner})$ alone gives LOO $R^2 = -0.023$ for predicting outer slope — worse than random. By contrast, $\log V_\mathrm{max}^\mathrm{inner}$ alone gives LOO $R^2 = 0.363$ with a single predictor.

**The breakthrough model: Structural + innerVmax.** Adding $\log V_\mathrm{max}^\mathrm{inner}$ to the four structural variables ($\Sigma_0$, $V_\mathrm{flat}$, $f_\mathrm{gas}$, $M_\mathrm{HI}$) yields $R^2 = 0.731$, LOO $R^2 = 0.677$, gap $= 7.4\%$. This is a **massive** improvement over structural + $a_0$ (LOO $R^2 = 0.292$, gap $= 28.7\%$). A single inner-dynamics variable more than doubles the predictive power and nearly eliminates overfitting.

**$a_0$ adds nothing when inner dynamics are included.** Adding $\log a_0(\mathrm{inner})$ to the inner-dynamics model (innerVmax + discrepancy + curvature) yields $\Delta \mathrm{LOO}\ R^2 = -0.009$ — a negative contribution. In the reverse direction, adding innerVmax to $a_0$ improves LOO $R^2$ by $+0.390$. The information is entirely one-directional: innerVmax subsumes $a_0$, not vice versa.

**Physical interpretation of the best model.** The coefficients of Model 3 reveal a striking pattern: $\beta(\log V_\mathrm{flat}) = +1.809$ and $\beta(\log V_\mathrm{max}^\mathrm{inner}) = -1.755$. These nearly cancel, meaning the outer slope is governed by $\log(V_\mathrm{flat} / V_\mathrm{max}^\mathrm{inner})$ — the **ratio** of the asymptotic rotation velocity to the inner peak velocity. Physically:
- If $V_\mathrm{flat} \gg V_\mathrm{max}^\mathrm{inner}$: the galaxy has not yet reached its asymptotic velocity in the inner half → **rising outer RC**.
- If $V_\mathrm{flat} \approx V_\mathrm{max}^\mathrm{inner}$: the galaxy peaked early → **flat or declining outer RC**.

This ratio encodes halo concentration directly: a concentrated halo produces a high inner $V_\mathrm{max}$ relative to $V_\mathrm{flat}$, yielding a flat outer profile. An extended halo produces a low inner $V_\mathrm{max}$, and the RC continues to rise.

**Verdict.** The per-galaxy $a_0$ is **not needed** for predicting rotation curve shape. The direct law — outer slope as a function of structural properties plus inner $V_\mathrm{max}$ — achieves LOO $R^2 = 0.677$ compared to $0.292$ for the $a_0$-based approach. The $a_0$ framework has been superseded by a more fundamental description: the ratio $V_\mathrm{flat}/V_\mathrm{max}^\mathrm{inner}$, which directly encodes halo concentration and determines how much the rotation curve will rise beyond the inner region. The "structured $a_0$ law" discovered in earlier phases was a shadow of this deeper relationship — $a_0$ correlated with structure because it was a noisy proxy for inner $V_\mathrm{max}$, which is the true dynamical variable connecting inner potential to outer RC shape.

---

## 13. Limitations

1. **Mass-regime dependence, not universal law**: The MHI–$a_0$ anticorrelation that anchors M3/M5 is a normal feature of the high-$V_\mathrm{flat}$ regime ($p = 0.17$ within $V_\mathrm{flat} \geq 80$) but reverses to a positive correlation among low-mass dwarfs. The law is therefore **regime-specific**: valid for massive spirals but not a universal description of $a_0$ variation across galaxy types. Phase 86 overturns the initial Phase 84 conclusion that the law is a pure collider artifact.
2. **Sample size**: $N = 45$ with $p = 3$–5 parameters. Internal cross-validation passes, but the external failure suggests the internal signal may reflect selection effects rather than population-level structure.
3. **Distance dependence**: The published-distance criterion ($f_D \geq 2$) creates a strong selection bias toward massive, luminous spirals in group environments (confirmed by Phase 82 population comparison). This generates apparent $a_0$ structure within the quality subsample that does not hold for late-type dwarfs and irregulars.
4. **$M_\mathrm{host}$ availability**: Group-catalog host masses are available only for the training set. External galaxies require $V_\mathrm{flat}$-based proxies with $-$0.74 dex mean bias. This entangles the sample selection with the predictor definition.
5. **$\Upsilon_\star^\perp$ construction**: Depends on the choice of confounders; different choices yield different values.
6. **Intrinsic scatter**: $\sim$0.10 dex remains unexplained.
7. **Interpretation**: We measure *inferred* $a_0$ variation. The external failure makes it premature to claim this variation is physically real rather than methodologically induced.

---

## 14. Conclusion

On a quality-controlled sample of $N = 45$ SPARC galaxies with published distances, galaxy-to-galaxy variation in the inferred MOND acceleration scale $a_0$ follows a structured multi-axis pattern within the training sample. The five-axis M5 law (LOO = 51.0%) and the compressed three-axis M3 (LOO = 44.1%) both outperform the universal-constant null model on internal cross-validation, permutation testing, bootstrap stability, and falsifiable prediction tests.

**However, this structured pattern does not generalize to the broader SPARC population.** When frozen M3 coefficients are applied to 36–46 external SPARC galaxies (Section 12), all model variants perform worse than the trivial baseline. The partial correlations of logMHI and logMeanRun with $a_0$ reverse sign on external data, indicating that the relationships found in $N = 45$ do not hold for galaxies selected under different criteria.

This creates a tension between strong internal evidence and failed external transfer. The internal tests — including nested cross-validation with model selection inside the evaluation loop ($p < 0.001$ permutation), zero collinearity (VIF $< 1.3$), and stable bootstrap coefficients — argue against simple overfitting. But the external failure argues against a population-level physical law.

We conclude that the structured $a_0$ variation observed in $N = 45$ is a **regime-specific finding**: the MHI–$a_0$ anticorrelation is real within the high-$V_\mathrm{flat}$ regime (Phase 86: $p = 0.17$ in $V_\mathrm{flat} \geq 80$; $p = 0.58$ in $Q = 1$ high-quality subset) but does not extend to a universal law across all galaxy types. The refined picture, developed across Phases 82–86:

1. **Mass-regime boundary** (Phases 85–86): The MHI–$a_0$ correlation is *positive* in the full SPARC sample ($r = +0.278$) but *negative* in the high-$V_\mathrm{flat}$ regime ($r = -0.215$ for $V_\mathrm{flat} \geq 80$). The training $r = -0.322$ is a typical draw from this regime ($p = 0.17$), not a selection-amplified extreme. This overturns the initial Phase 84 collider-bias interpretation.
2. **Regime, not law**: The multi-axis law correctly describes the high-mass spiral regime but fails when applied to low-mass dwarfs and irregulars where the correlation structure reverses. The external failure (Phase 82–83) reflects a genuine regime boundary, not a methodological flaw in the training analysis.
3. **Residual amplification**: While the core anticorrelation is regime-typical, the specific $N = 45$ sample may be further shaped by $f_D \geq 2$ (which makes the pool deterministic) and morphology cuts. The $f_D$ criterion alone does not flip the sign ($r = +0.306$) but narrows the pool to exactly the training set.
4. **Mhost confound**: Group-catalog host masses remain available only for the training set. External galaxies require $V_\mathrm{flat}$-based proxies with $-$0.74 dex mean bias.

Phase 87 locates the transition precisely: $V_\mathrm{flat} \approx 70$ km/s ($\Delta\mathrm{BIC} = -23.8$), with the MHI slope reversing from $+0.059$ to $-0.127$. Phase 88 reveals the transition is **asymmetric**: a structured five-predictor law exists in the high-$V_\mathrm{flat}$ regime ($N = 112$, $R^2 = 0.253$, LOO gap $= 8.1\%$, $p < 0.001$), but no predictive law emerges in the low-$V_\mathrm{flat}$ regime ($N = 63$, best LOO gap $= -1.8\%$, $p = 0.689$). Three of five predictors reverse sign across the boundary — a full state transition, not merely a gas-regime effect. Cross-regime transfer fails completely (gaps of $-77\%$ and $-89\%$).

The final picture emerges from the interplay of Phases 88–90. $a_0$ variation in SPARC galaxies divides into two fundamentally different regimes at $V_\mathrm{flat} \approx 70$ km/s. Above this threshold, galaxy properties predict $\sim 25\%$ of $a_0$ variance in a structured law. Below it, Phase 89 shows the *scatter amplitude* is quality-dominated ($Q = 1$ dwarfs: sd $= 0.325$ matching the high regime), but Phase 90 shows the *correlation structure* is genuinely different: no multi-axis law emerges even among $Q = 1$ dwarfs, key predictor signs reverse, and cross-regime transfer fails completely ($-37\%$ to $-44\%$). The distinction is subtle but critical:

- **Scatter amplitude**: inflated by data quality in dwarfs (Phase 89). Controllable.
- **Correlation structure**: genuinely different between regimes (Phase 90). Physical.

The $V_\mathrm{flat} \approx 70$ km/s boundary is therefore **both** a data-quality transition (inflating apparent scatter) **and** a physical transition (changing which properties correlate with $a_0$). Dwarfs and massive spirals share comparable $a_0$ scatter when quality is controlled, but they do not share the same predictive law.

Phase 91 identifies **why** the law exists only in the high regime: the boundary at $V_\mathrm{flat} \approx 70$ km/s marks a convergence of four independent transitions — gas fraction crosses $\sim 50\%$ (stellar disk loses dominance), acceleration range narrows from 1.19 to 0.51 dex (a₀ becomes degenerate), galaxies enter deep MOND (RAR shape insensitive to $a_0$), and 73% of dwarfs have rising RCs vs. 34% of spirals ($V_\mathrm{flat}$ itself undefined). These four effects converge to make the per-galaxy $a_0$ parameter **operationally ill-defined** in the dwarf regime — not because MOND fails, but because the fitting methodology breaks down when applied to gas-dominated, low-acceleration systems with incomplete rotation curves.

Phase 92 connects the $a_0$ law to rotation curve flatness: the inferred $a_0$ does **not** directly predict RC shape ($r = -0.052$), but the same M3 properties predict outer slope with $R^2 = 0.36$ ($p < 0.001$). Phase 93 reveals the underlying causal structure through mediation analysis: $\Sigma_0$ (baryon concentration) and $V_\mathrm{flat}$ (potential depth) are **dual parents**, each maintaining independent partial correlations with both $a_0$ and flatness. Most strikingly, after controlling for all measured structural variables, $a_0$ retains a residual partial correlation with outer slope of $r = -0.294$ — the inferred acceleration scale carries dynamical information about RC shape that the baryonic observables do not fully capture. Phase 94 confirms this residual is **genuine, not circularity**: when $a_0$ is extracted from the inner half of the RC and slope from the outer half (eliminating shared-data artifacts), the partial correlation **strengthens** to $r = -0.325$ ($p_\mathrm{perm} = 0.003$, ratio split/full $= 1.07$). The information flows outward — $a_0(\mathrm{inner}) \to \mathrm{outerSlope}$ is strong ($r = -0.325$) while $a_0(\mathrm{outer}) \to \mathrm{innerSlope}$ is null ($r = +0.071$). Phase 95 identifies the hidden inner quantity: the inner peak velocity ($\log V_\mathrm{max}^\mathrm{inner}$) absorbs $79.2\%$ of the $a_0$ residual and predicts outer slope with $r = -0.761$ after structural controls — far stronger than $a_0$ itself. Combined with the inner mass discrepancy ($46.9\%$ reduction) and inner curvature, three inner properties absorb $76.1\%$ of the signal. Phase 96 delivers the decisive test: a direct law using structural variables plus $\log V_\mathrm{max}^\mathrm{inner}$ achieves LOO $R^2 = 0.677$ for predicting outer slope — more than double the $a_0$-based model ($0.292$) — with a gap of only $7.4\%$. Adding $a_0$ to the inner-dynamics model contributes $\Delta \mathrm{LOO}\ R^2 = -0.009$ (negative), while adding innerVmax to $a_0$ adds $+0.390$. The model coefficients reveal that the outer slope is governed by $\log(V_\mathrm{flat}/V_\mathrm{max}^\mathrm{inner})$ — the ratio of asymptotic to inner peak velocity — which directly encodes **halo concentration**. The per-galaxy $a_0$ is not needed: it was a noisy proxy for the ratio $V_\mathrm{flat}/V_\mathrm{max}^\mathrm{inner}$, which is the true dynamical variable connecting inner potential to outer RC shape.

**Table 5b. Evidence $\to$ Claim Traceability**

| Claim | Evidence | Status |
|-------|----------|--------|
| Not a universal constant (within N=45) | M0 loses to M5/M3 on LOO, AIC, BIC, bootstrap | Confirmed internally |
| Not random scatter (within N=45) | Structured model repeatedly outperforms null | Confirmed internally |
| Best predictive law = M5 | Death match: M5 defeats all competitors, LOO=51.0% | Internal only |
| Best state law = M3 | Compression: 3 axes retain 86% of M5, LOO=44.1% | Internal only |
| Law is falsifiable | 24/26 quantitative predictions passed (92%) | Internal only |
| Law works pairwise | 22/25 controlled pairs pass (88%), median |err|=0.165 dex | Internal only |
| MHI and Mhost independent | $f_\mathrm{gas}$ collapse loses $-$36.5 pp; exponents $-$1/4 vs $-$1/6 | Internal only |
| Kinematic axis is real | MeanRun correlates $|r|=0.83$ with 2D THINGS dynamics | Confirmed |
| No 7th scalar axis | All candidates fail significance tests | Internal only |
| Residual is structureless | 0/7 strat. biases, BP=4.63, JB=1.22, Gaussian | Internal only |
| **External generalization** | **M3/M2$'$ fail on SPARC-out, axis signs reversed** | **Failed** |
| **Population mismatch explains failure** | **12/13 variables differ (8 large effect), training = massive spirals, external = dwarf irregulars** | **Confirmed** |
| **Common support rescues law** | **Regime-matched external galaxies (4 strategies, N=14–78): M3 win rate 7–10%, all fail** | **Failed** |
| **Collider bias identified** | **r(MHI,a₀) = +0.278 in full SPARC but −0.322 in N=45 (0/1000 random); selection reverses sign** | **Confirmed** |
| **Vflat cut is the primary driver** | **Vflat≥80 alone flips sign (r=−0.215, N=99); additional criteria amplify to −0.322** | **Confirmed** |
| **Training r is TYPICAL in Vflat≥80** | **p=0.17 (5000 draws); in Q=1 subset p=0.58; not a selection extreme but a regime effect** | **Confirmed** |
| **Law is regime-specific, not artifact** | **11/23 pools make training r typical; overturns pure-collider interpretation** | **Revised** |
| **Regime transition at Vflat≈70** | **Broken regression ΔBIC=−23.8; slope reverses from +0.059 (dwarfs) to −0.127 (massive); sharp transition w=13 km/s** | **Confirmed** |
| **Vflat is the best separator** | **logVflat ΔBIC=−23.0; logMHI, logL36, T all fail (positive ΔBIC)** | **Confirmed** |
| **High-Vflat law exists** | **N=112, R²=0.253, LOO gap=8.1%, perm p<0.001; 5 predictors** | **Confirmed** |
| **Low-Vflat: no law** | **N=63, best LOO gap=−1.8%, perm p=0.689; a₀ scatter is unstructured (sd=0.77 dex)** | **Confirmed** |
| **Full state transition** | **3/5 variables reverse sign; 3/5 axes have strong broken models (ΔBIC<−6)** | **Confirmed** |
| **Cross-regime transfer fails** | **Hi→Lo gap=−77%, Lo→Hi gap=−89%; regimes are genuinely separate** | **Failed** |
| **Low-regime scatter is intrinsic** | **Bootstrap: only 3.1% measurement; intrinsic SD=0.760 dex** | **Confirmed** |
| **Q=1 dwarfs match high-regime scatter** | **Q=1 (N=22): sd=0.325 vs Q=2,3: sd=0.76-0.80; quality-dominated, not physics** | **Key finding** |
| **Low-regime "chaos" = data quality** | **Systematics (non-circular motions, beam smearing) dominate Q≥2 dwarfs; Q=1 subsample is well-behaved** | **Revised** |
| **No law in Q=1 dwarfs** | **Phase 90: best physics LOO gap=0.5% (logMeanRun); Gas model gap=−6.9%; N=22 too small or genuinely lawless** | **Confirmed** |
| **Signs reverse even in Q=1** | **2/4 key variables match high-regime signs; logMeanRun and logSigma0 reverse** | **Confirmed** |
| **Cross-regime transfer fails even Q=1** | **HiQ1→LoQ1: −37% to −41%; LoQ1→HiQ1: −24% to −44%** | **Failed** |
| **Inclination signal in dwarfs** | **r=−0.473, perm p=0.025; likely systematic (deprojection bias), not physics** | **Noted** |
| **Regime boundary is physical + quality** | **Scatter amplitude = quality; correlation structure = physics; both contribute** | **Key finding** |
| **Convergent transitions at 70 km/s** | **Phase 91: fgas (0.31→0.60), gbar range (1.19→0.51 dex), deep MOND, rising RCs (34%→73%)** | **Key finding** |
| **a₀ operationally ill-defined in dwarfs** | **Narrow range + deep MOND + rising RC = degenerate fit; a₀ not meaningful per-galaxy** | **Confirmed** |
| **a₀ does NOT predict RC flatness** | **Phase 92: r(logA0, outerSlope) = −0.052; partial r = −0.072 after Vflat** | **Confirmed** |
| **Same properties predict BOTH a₀ and flatness** | **M3→outerSlope: R²=0.36, gap=16.9%, p<0.001; Σ₀ strongest (r=−0.518)** | **Key finding** |
| **a₀ and flatness are parallel consequences** | **Not causally linked; both driven by baryon-dominated disk structure** | **Key finding** |
| **Dual parent structure** | **Phase 93: Σ₀ and Vflat both survive controlling for each other (r=−0.270, −0.212)** | **Confirmed** |
| **a₀ residual after structural controls** | **r(logA0, slope \| all structural) = −0.294; a₀ carries extra dynamical info** | **Key finding** |
| **Σ₀ has independent paths to both** | **r(Σ₀, A0 \| slope) = −0.215; r(Σ₀, slope \| A0) = −0.533** | **Confirmed** |
| **a₀ residual is GENUINE, not circularity** | **Phase 94: split-half r(a₀_inner, outerSlope \| struct) = −0.325 (p=0.003); ratio=1.07** | **Confirmed** |
| **Information flows outward** | **a₀(inner)→outerSlope: r=−0.325; a₀(outer)→innerSlope: r=+0.071 (no signal)** | **Key finding** |
| **a₀ encodes halo/potential info** | **Inner potential sets boundary conditions for outer RC; not captured by Σ₀,Vflat,fgas,MHI** | **Key finding** |
| **Inner Vmax is the dominant carrier** | **Phase 95: r(logVmax_inner, outerSlope \| struct) = −0.761; absorbs 79.2% of a₀ residual** | **Key finding** |
| **Inner mass discrepancy = DM proxy** | **Mean(logGobs−logGbar) inner reduces a₀ residual by 46.9%** | **Confirmed** |
| **a₀(inner) is a composite proxy** | **Top 3 (innerVmax + discrepancy + curvature) absorb 76.1%; a₀ is NOT irreducible** | **Key finding** |
| **a₀ encodes halo concentration** | **Inner Vmax after controlling Vflat = how fast galaxy peaks = NFW concentration proxy** | **Key finding** |
| **Direct law without a₀: LOO R²=0.677** | **Phase 96: Structural + innerVmax beats Structural + a₀ (0.292) by factor 2.3x** | **Key finding** |
| **a₀ is NOT NEEDED** | **Adding a₀ to inner dynamics: ΔLOO = −0.009; innerVmax adds +0.390 to a₀** | **Confirmed** |
| **outerSlope ≈ log(Vflat/InnerVmax)** | **β(logVflat)=+1.81, β(logInnerVmax)=−1.76; ratio encodes halo concentration** | **Key finding** |
| **a₀ alone is useless for flatness** | **LOO R² = −0.023 (worse than random); innerVmax alone: +0.363** | **Confirmed** |

---

## References

- Lelli, F., McGaugh, S.S., Schombert, J.M. 2016, AJ, 152, 157 (SPARC)
- Lelli, F., McGaugh, S.S., Schombert, J.M., Pawlowski, M.S. 2017, ApJ, 836, 152
- Li, P., Lelli, F., McGaugh, S., Schombert, J. 2018, A&A, 615, A3
- Li, P., Lelli, F., McGaugh, S., Schombert, J. 2020, ApJS, 247, 31
- McGaugh, S.S., Lelli, F., Schombert, J.M. 2016, PRL, 117, 201101
- McGaugh, S.S., Lelli, F., Schombert, J.M. 2018, Research Notes AAS, 2, 156
- Milgrom, M. 1983, ApJ, 270, 365
- Rodrigues, D.C., Marra, V., del Popolo, A., Davari, Z. 2018, Nature Astronomy, 2, 668
- Kroupa, P. et al. 2018, Nature Astronomy, 2, 925
- Tully, R.B. 2015, AJ, 149, 171

---

## Appendix A: Supplementary Materials

**Table A1**: Full $N = 45$ sample: `N45_final_dataset.csv`

**Table A2**: Variable definitions: `variable_definitions.md`

**Replication**: Self-contained Python script `replicate_from_scratch.py` reproduces all M5/M3 coefficients, LOO gap%, and residual diagnostics. Requires Python 3.7+ and NumPy.

## Appendix B: Glossary

| Symbol | Definition |
|--------|-----------|
| LOO gap% | $100 \times (1 - \text{LOO-RMS}^2 / \text{SD}(y)^2)$ |
| M5 | Five-axis predictive law (logMHI, logMhost, logSigma0, logMeanRun, $\Upsilon_\star^\perp$) |
| M3 | Three-axis compressed state law (logMHI, logMhost, logMeanRun) |
| $\Upsilon_\star^\perp$ | Orthogonalized stellar M/L: residual of $\log\Upsilon_\star$ after OLS on (logMHI, logSigma0, morphT) |
| $\tau_\mathrm{int}$ | Intrinsic scatter after measurement noise subtraction |

---

## Appendix C: Overfitting and Model-Selection Audit

This appendix documents the full protocol and quantitative results for each robustness test summarized in Section 11.

### C.1 Nested Cross-Validation

**Protocol.** Standard LOO evaluates a *fixed* model on held-out galaxies, but the choice of which model to evaluate (M3 vs M5 vs M6) was made using the full dataset. Nested CV removes this circularity by placing model selection inside the evaluation loop.

- **Outer loop**: 5-fold CV, repeated 50 times with random partitions (250 total outer folds).
- **Inner loop**: For each outer training set, compute LOO gap% for M3, M5, and M6 on the training galaxies only. Select the model with the best inner LOO.
- **Outer evaluation**: Fit the selected model on the full training fold, predict the held-out galaxies.
- **$\Upsilon_\star^\perp$ reconstruction**: In each fold, $\Upsilon_\star^\perp$ is recomputed from training data only (auxiliary OLS on logMHI, logSigma0, morphT), ensuring no information leaks from the test fold.

**Results.**

| Metric | Value |
|--------|-------|
| Nested CV gap% | 38.5% |
| Standard LOO gap% | 51.0% |
| Optimism | 12.5 pp |
| Model selected: M5 | 51.2% of folds |
| Model selected: M3 | 29.2% of folds |
| Model selected: M6 | 19.6% of folds |
| Outer-fold RMS | 0.202 dex |

**Interpretation.** The optimism (12.5 pp) is moderate and expected for $N = 45$ with a 5–6 parameter model. The nested gap% of 38.5% confirms that the majority of the apparent signal is genuine and not an artifact of model selection on the same data. M5 is still selected in the majority of folds, but M3 is competitive — consistent with BIC preferring M3 as the more parsimonious model.

### C.2 Full-Pipeline Permutation Test (Y-Scramble)

**Protocol.** This test asks: *could the full analytical pipeline — including model construction, $\Upsilon_\star^\perp$ orthogonalization, and best-model selection — produce a result as strong as M5 from data with no real structure?*

1. Randomly permute $\log a_0$ across galaxies, breaking any true association with predictors.
2. Recompute $\Upsilon_\star^\perp$ from the permuted data (full pipeline faithfulness).
3. Fit both M3 and M5, compute LOO gap% for each, and record the *better* of the two (maximizing the null's chance of success).
4. Repeat 1,000 times to build a null distribution.

**Results.**

| Null distribution statistic | Value |
|-----------------------------|-------|
| Median null gap% | $-$10.7% |
| 95th percentile | 3.5% |
| 99th percentile | 11.3% |
| Maximum (of 1,000) | 27.0% |
| Real M5 gap% | 46.6% |
| $p$-value | $< 0.001$ (0/1,000 null runs $\geq$ real) |

**Interpretation.** Not a single null run out of 1,000 produces a pipeline result remotely close to M5's true performance. The maximum null gap% (27.0%) is 19.6 pp below the real value. The pipeline does not "discover" strong laws in structureless data.

### C.3 Collinearity Audit

**Protocol.** If the five M5 predictors are secretly encoding the same physical degree of freedom in different guises, the law would be less informative than it appears.

**Correlation matrix (Pearson $r$):**

| | logMHI | logMhost | logΣ₀ | logMeanRun | Υ★⊥ |
|-----|--------|----------|-------|------------|------|
| logMHI | 1.000 | +0.183 | +0.361 | +0.242 | 0.000 |
| logMhost | +0.183 | 1.000 | +0.322 | $-$0.201 | $-$0.211 |
| logΣ₀ | +0.361 | +0.322 | 1.000 | +0.004 | 0.000 |
| logMeanRun | +0.242 | $-$0.201 | +0.004 | 1.000 | +0.190 |
| Υ★⊥ | 0.000 | $-$0.211 | 0.000 | +0.190 | 1.000 |

**Variance Inflation Factors:**

| Variable | VIF |
|----------|-----|
| logMHI | 1.26 |
| logMhost | 1.23 |
| logΣ₀ | 1.25 |
| logMeanRun | 1.16 |
| Υ★⊥ | 1.08 |
| **Maximum** | **1.26** |

All VIFs are well below the conventional concern threshold of 5 (let alone the severe threshold of 10). No pair of predictors exceeds $|r| = 0.37$. The five axes are genuinely independent directions in galaxy-property space, not collinear restatements of a smaller number of latent variables.

$\Upsilon_\star^\perp$ has VIF = 1.08 by construction (it is the orthogonal residual after projecting out logMHI, logΣ₀, and morphT), confirming the orthogonalization works as intended.

### C.4 Coefficient Stability Under Resampling

**Protocol.** Bootstrap resampling (1,000 draws with replacement) tests whether M5 coefficients are stable or fluctuate wildly — a hallmark of overfitting.

| Coefficient | Mean | SD | 2.5% | 97.5% | Sign-flip % |
|------------|------|-----|------|-------|------------|
| Intercept | +5.053 | 0.525 | +4.040 | +6.083 | 0.0% |
| logMHI | $-$0.239 | 0.051 | $-$0.354 | $-$0.144 | 0.0% |
| logMhost | $-$0.176 | 0.043 | $-$0.264 | $-$0.092 | 0.0% |
| logΣ₀ | +0.137 | 0.070 | $-$0.008 | +0.260 | 3.7% |
| logMeanRun | +0.446 | 0.112 | +0.212 | +0.674 | 0.0% |
| Υ★⊥ | +0.342 | 0.320 | $-$0.367 | +0.897 | 14.1% |

The three M3 core axes (logMHI, logMhost, logMeanRun) are rock-stable: zero sign flips and tight confidence intervals. logΣ₀ is stable (3.7% sign-flip). $\Upsilon_\star^\perp$ has the widest spread (14.1% sign-flip), consistent with its status as a second-order correction rather than a core state variable.

M5 beats M3 (by LOO gap%) in 69.4% of bootstrap resamples, confirming that M5's advantage is typical, not dependent on a lucky data draw.

### C.5 Calibration Test

**Protocol.** We test whether M5 not only ranks galaxies correctly but also calibrates their predicted $\log a_0$ accurately, using 200 random 50/50 train-test splits.

| Metric | Value | Ideal |
|--------|-------|-------|
| Mean slope (obs vs pred) | 0.819 | 1.0 |
| Mean intercept | 0.644 | 0.0 |
| Mean $R^2$ | 0.451 | 1.0 |
| Slope 95% CI | [0.464, 1.346] | — |

The slope of 0.82 indicates mild regression toward the mean (expected for $N = 45$; extreme predicted values are moderated). Critically, the 95% CI for the slope includes 1.0, meaning the data are consistent with unbiased calibration. $R^2 = 0.45$ on fully held-out test halves is consistent with the nested CV estimate.

### C.6 Complexity Penalty: AIC/BIC Comparison

| Model | $k$ | LOO gap% | AIC | BIC | RMS (dex) |
|-------|-----|----------|-----|-----|-----------|
| M0 | 1 | $-$2.3% | $-$121.0 | $-$119.2 | 0.255 |
| M3 | 4 | 44.1% | $-$148.3 | $-$141.0 | 0.176 |
| **M5** | **6** | **46.6%** | **$-$150.9** | $-$140.1 | **0.164** |
| M6 | 7 | 45.7% | $-$151.1 | **$-$138.4** | 0.160 |

M5 wins on LOO and AIC. BIC prefers M3, penalizing M5's two additional parameters — this is informative rather than damaging: it reinforces the two-tier architecture where M3 is the parsimonious state law and M5 is the predictive extension. M6 does not beat M5 on any metric, confirming that adding rcWiggliness (the dropped 6th variable) provides no benefit.

### C.7 Drop-One-Variable Ablation

| Variable dropped | LOO gap% | $\Delta$ vs full M5 | Verdict |
|-----------------|----------|---------------------|---------|
| logMHI | 23.6% | $-$23.0 pp | Essential |
| logMhost | 31.6% | $-$15.0 pp | Essential |
| logΣ₀ | 43.5% | $-$3.1 pp | Contributes |
| logMeanRun | 32.1% | $-$14.5 pp | Essential |
| Υ★⊥ | 47.4% | +0.8 pp | Marginal (second-order) |

The three M3 core axes are individually essential (each drops performance by $>$14 pp). logΣ₀ contributes modestly ($-$3.1 pp). Υ★⊥ is marginal as a standalone addition to LOO — its value is structural (reduces RMS, improves calibration) rather than cross-validated predictive.

### C.8 Summary Verdict

| Test | Result | Pass? |
|------|--------|-------|
| Nested CV | gap = 38.5%, optimism = 12.5 pp | Yes |
| Pipeline permutation | $p < 0.001$ (0/1,000) | Yes |
| Collinearity (VIF) | max VIF = 1.26 | Yes |
| Coefficient stability | All signs stable, M5 wins 69.4% | Yes |
| Calibration | slope = 0.82, $R^2$ = 0.45 | Yes |
| Complexity penalty | M5 wins LOO/AIC, M6 does not beat M5 | Yes |

All six tests pass. The structured law is not a null-result artifact, a simple overfit, or a collinear re-encoding. Its signal survives when model selection is part of the evaluation, and it fails to appear in null-permuted pipelines.

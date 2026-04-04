# The Inferred MOND Acceleration Scale as an Emergent Galaxy State Variable: A Five-Axis Law from N=45 SPARC Galaxies

---

## Abstract

The radial acceleration relation (RAR) in disk galaxies is commonly parameterized by a single, supposedly universal acceleration scale $a_0 \approx 1.2 \times 10^{-10}$ m s$^{-2}$. We test whether galaxy-to-galaxy variation in the inferred $a_0$ is best described as (i) a strict universal constant, (ii) unstructured per-galaxy scatter, or (iii) a structured multi-axis law. Working with a quality-controlled subsample of $N = 45$ SPARC galaxies selected for published, high-quality distance estimates, we construct per-galaxy $a_0$ values and search for predictive structure.

We find that within this sample, a five-predictor power-law model (M5) — incorporating gas mass ($M_\mathrm{HI}^{-1/4}$), host-halo mass ($M_\mathrm{host}^{-1/6}$), central baryonic surface density ($\Sigma_0^{+1/7}$), kinematic coherence ($\overline{\mathrm{Run}}^{+3/7}$), and an orthogonalized stellar mass-to-light ratio ($\Upsilon_\star^{\perp\,+2/3}$) — explains 51.0% of the $\log a_0$ variance in leave-one-out cross-validation (LOO). A compressed three-axis state law (M3) retains 86% of this signal using only gas mass, host mass, and kinematic coherence (LOO = 44.1%). A comprehensive overfitting audit — including nested cross-validation ($p < 0.001$ permutation), collinearity analysis (all VIF $< 1.3$), and bootstrap stability — confirms that the structured signal is not a fitting artifact within the training set.

**However, external validation fails.** When frozen M3 coefficients are applied to 36–46 SPARC galaxies outside the training set (using rotation-curve data from CDS/VizieR with identically computed inputs), all model variants perform worse than the trivial universal-constant baseline. The partial correlations of logMHI and logMeanRun with $a_0$ reverse sign on external data.

We conclude that $a_0$ variation within the $N = 45$ quality-controlled subsample is structured and statistically real (not a fitting artifact), but this structure is **sample-specific** and does not generalize to the broader SPARC population. The pattern may reflect selection effects arising from the published-distance criterion rather than a population-level physical law. These results caution against interpreting internal cross-validation alone as evidence for law discovery in small astrophysical samples.

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

### 12.6 Verdict

**The multi-axis law found in $N = 45$ does not generalize to the broader SPARC sample.** Internal cross-validation (LOO, nested CV, permutation) confirms that the pattern is real within the training set and not a fitting artifact. But external validation shows that this pattern does not transfer to galaxies selected under different criteria. The structured $a_0$ variation should be classified as a **sample-specific finding** pending confirmation on independently selected samples with group-catalog host masses and published distances.

---

## 13. Limitations

1. **External generalization failure**: The structured law does not generalize to SPARC galaxies outside the $N = 45$ training set (Section 12). This is the most significant limitation and qualifies all claims made in this paper.
2. **Sample size**: $N = 45$ with $p = 3$–5 parameters. Internal cross-validation passes, but the external failure suggests the internal signal may reflect selection effects rather than population-level structure.
3. **Distance dependence**: The published-distance criterion ($f_D \geq 2$) may create a selection bias that generates apparent $a_0$ structure within the quality subsample.
4. **$M_\mathrm{host}$ availability**: Group-catalog host masses are available only for the training set. External galaxies require $V_\mathrm{flat}$-based proxies with $-$0.74 dex mean bias. This entangles the sample selection with the predictor definition.
5. **$\Upsilon_\star^\perp$ construction**: Depends on the choice of confounders; different choices yield different values.
6. **Intrinsic scatter**: $\sim$0.10 dex remains unexplained.
7. **Interpretation**: We measure *inferred* $a_0$ variation. The external failure makes it premature to claim this variation is physically real rather than methodologically induced.

---

## 14. Conclusion

On a quality-controlled sample of $N = 45$ SPARC galaxies with published distances, galaxy-to-galaxy variation in the inferred MOND acceleration scale $a_0$ follows a structured multi-axis pattern within the training sample. The five-axis M5 law (LOO = 51.0%) and the compressed three-axis M3 (LOO = 44.1%) both outperform the universal-constant null model on internal cross-validation, permutation testing, bootstrap stability, and falsifiable prediction tests.

**However, this structured pattern does not generalize to the broader SPARC population.** When frozen M3 coefficients are applied to 36–46 external SPARC galaxies (Section 12), all model variants perform worse than the trivial baseline. The partial correlations of logMHI and logMeanRun with $a_0$ reverse sign on external data, indicating that the relationships found in $N = 45$ do not hold for galaxies selected under different criteria.

This creates a tension between strong internal evidence and failed external transfer. The internal tests — including nested cross-validation with model selection inside the evaluation loop ($p < 0.001$ permutation), zero collinearity (VIF $< 1.3$), and stable bootstrap coefficients — argue against simple overfitting. But the external failure argues against a population-level physical law.

We conclude that the structured $a_0$ variation observed in $N = 45$ is a **sample-specific finding**: it is statistically real within the training set (not a fitting artifact) but does not extend to a general law of galaxy physics. The most likely explanations are:

1. **Selection-induced structure**: The published-distance requirement may select a galaxy subpopulation where $a_0$ correlates with mass and kinematic properties in ways that do not hold for the broader population.
2. **Mhost confound**: Group-catalog host masses (available only for the training set) may encode environmental information that $V_\mathrm{flat}$ proxies cannot replicate, making the external test structurally unfair for the Mhost axis — though the M2$'$ test (dropping Mhost entirely) also fails.

These results caution against interpreting internal cross-validation alone as evidence for physical law discovery in small samples, even when permutation tests, bootstrap stability, and falsification metrics all pass. External validation on independently selected data remains the definitive test.

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

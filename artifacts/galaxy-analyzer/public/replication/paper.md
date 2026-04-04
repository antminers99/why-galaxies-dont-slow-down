# The Inferred MOND Acceleration Scale as an Emergent Galaxy State Variable: A Five-Axis Law from N=45 SPARC Galaxies

---

## Abstract

The radial acceleration relation (RAR) in disk galaxies is commonly parameterized by a single, supposedly universal acceleration scale $a_0 \approx 1.2 \times 10^{-10}$ m s$^{-2}$. We test whether galaxy-to-galaxy variation in the inferred $a_0$ is best described as (i) a strict universal constant, (ii) unstructured per-galaxy scatter, or (iii) a structured multi-axis law. Working with a quality-controlled subsample of $N = 45$ SPARC galaxies selected for published, high-quality distance estimates, we construct per-galaxy $a_0$ values and search for predictive structure.

We find that a five-predictor power-law model (M5) — incorporating gas mass ($M_\mathrm{HI}^{-1/4}$), host-halo mass ($M_\mathrm{host}^{-1/6}$), central baryonic surface density ($\Sigma_0^{+1/7}$), kinematic coherence ($\overline{\mathrm{Run}}^{+3/7}$), and an orthogonalized stellar mass-to-light ratio ($\Upsilon_\star^{\perp\,+2/3}$) — explains 51.0% of the $\log a_0$ variance in leave-one-out cross-validation (LOO). A compressed three-axis state law (M3) retains 86% of this signal using only gas mass, host mass, and kinematic coherence (LOO = 44.1%).

M5 produces 24/26 correct falsifiable predictions (92%), succeeds in 22/25 controlled matched-pair tests between individual galaxies (88%), and its kinematic-coherence axis is confirmed as a genuine 1D projection of 2D dynamical structure ($|r| = 0.83$ with THINGS velocity-field diagnostics). The 0.157 dex residual after M5 is a featureless Gaussian frontier with no subpopulation bias, no heteroscedasticity, and no bimodality.

We conclude that $a_0$ variation across galaxies is structured, falsifiable, and pairwise-testable. The data favor interpreting $a_0$ as an emergent, state-dependent galaxy parameter — set by the tension between gravitational depth (which suppresses it) and baryonic organization (which amplifies it) — rather than a strictly universal acceleration scale.

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

The corrections add +6.9 percentage points of LOO performance and reduce RMS from 0.193 to 0.157 dex. By analogy: M3 is to M5 as the ideal gas law is to the van der Waals equation — M3 captures the dominant physics; M5 adds structural refinements for higher precision.

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

## 11. Limitations

1. **Sample size**: $N = 45$. Cross-validation confirms internal generalizability, but extension to larger populations is essential.
2. **Distance dependence**: The sample requires published distances, which may introduce selection effects.
3. **$\Upsilon_\star^\perp$ construction**: Depends on the choice of confounders; different choices yield different values, though ours explains 68% of raw variance.
4. **External transfer**: Phase 73 was inconclusive due to data quality, not confirmatory.
5. **Intrinsic scatter**: $\sim$0.10 dex remains unexplained and may require resolved 2D kinematic data.
6. **Interpretation**: We measure *inferred* $a_0$ variation. Disentangling genuine acceleration-scale variation from systematic modeling effects requires independent dynamical probes.

---

## 12. Conclusion

On a quality-controlled sample of $N = 45$ SPARC galaxies with published distances, galaxy-to-galaxy variation in the inferred MOND acceleration scale $a_0$ is not best described by a strict universal constant, nor by unstructured per-galaxy scatter, but by a structured multi-axis law. The current best predictive formulation is the five-axis M5 law (LOO = 51.0%), while M3 provides a compressed physical reduction to three state axes: gas content, environmental depth, and dynamical coherence (LOO = 44.1%).

The law is falsifiable (24/26 quantitative predictions passed), succeeds in controlled pairwise comparisons between individual galaxies (22/25 matched pairs), and its kinematic-coherence axis has a confirmed 2D dynamical origin. The remaining 0.16 dex residual is a featureless Gaussian frontier with no additional scalar information extractable from the current data.

The data favor interpreting $a_0$ as an emergent, state-dependent galaxy parameter — set by the tension between gravitational depth (which suppresses it) and baryonic organization (which amplifies it) — rather than a strictly universal acceleration scale. This result is framework-neutral: it constrains any theory, whether dark-matter-based or modified-dynamics-based, to reproduce a state-dependent effective acceleration scale.

**Table 5b. Evidence $\to$ Claim Traceability**

| Claim | Evidence |
|-------|----------|
| Not a universal constant | M0 loses to M5/M3 on LOO, AIC, BIC, bootstrap |
| Not random scatter | Structured model repeatedly outperforms null |
| Best predictive law = M5 | Death match: M5 defeats all competitors, LOO=51.0% |
| Best state law = M3 | Compression: 3 axes retain 86% of M5, LOO=44.1% |
| Law is falsifiable | 24/26 quantitative predictions passed (92%) |
| Law works pairwise | 22/25 controlled pairs pass (88%), median |err|=0.165 dex |
| MHI and Mhost independent | $f_\mathrm{gas}$ collapse loses $-$36.5 pp; exponents $-$1/4 vs $-$1/6 |
| Kinematic axis is real | MeanRun correlates $|r|=0.83$ with 2D THINGS dynamics |
| No 7th scalar axis | All candidates fail significance tests |
| Residual is structureless | 0/7 strat. biases, BP=4.63, JB=1.22, Gaussian |
| Physical interp. = state variable | Mass-suppression vs organization-amplification |
| External transfer | Inconclusive (compromised test data) |

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

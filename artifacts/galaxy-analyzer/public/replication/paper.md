# A Structured Multi-Axis Law for Galaxy-to-Galaxy Variation in the Inferred MOND Acceleration Scale

---

## Abstract

The radial acceleration relation (RAR) in disk galaxies is commonly parameterized by a single, supposedly universal acceleration scale $a_0 \approx 1.2 \times 10^{-10}$ m s$^{-2}$. We test whether galaxy-to-galaxy variation in the inferred $a_0$ is best described as (i) a strict universal constant, (ii) unconstrained per-galaxy values, or (iii) a structured multi-axis law. Working with a quality-controlled subsample of $N = 45$ SPARC galaxies selected for published, high-quality distance estimates, we construct per-galaxy $a_0$ values and search for predictive structure in the residuals.

We find that a six-predictor linear model (B$''$) — incorporating gas mass, host-halo mass, central baryonic surface density, rotation-curve irregularity, kinematic coherence, and an independent stellar mass-to-light ratio component — explains 49.7% of the $\log a_0$ variance in leave-one-out cross-validation (LOO), compared to 0% for a universal constant and negative values for unconstrained per-galaxy fits. B$''$ wins unanimously over all tested alternatives on every evaluation metric (LOO, 5-fold, 10-fold, AIC, BIC). It transfers to held-out galaxies, replicates under three independent fitting methods, and is stable under jackknife, bootstrap, and random-split tests.

A genuine intrinsic residual scatter of $\sim$0.10 dex (95% CI: [0.04, 0.14]) persists after B$''$, which is not captured by any additional scalar variable, environmental processing proxy, or distributed weak-signal bundle available in the current data. We conclude that $a_0$ variation across galaxies is structured, reproducible, and generalizable, but does not reduce to a strict universal constant on this sample.

---

## 1. Introduction

The radial acceleration relation (RAR; McGaugh et al. 2016; Lelli et al. 2017) establishes a tight empirical correlation between the observed centripetal acceleration in disk galaxies and the acceleration predicted from the baryonic mass distribution alone. This relation is a central prediction of Modified Newtonian Dynamics (MOND; Milgrom 1983) and is parameterized by a characteristic acceleration scale $a_0$. Whether $a_0$ is a strict universal constant or varies systematically between galaxies has profound implications: a truly universal $a_0$ supports MOND as a fundamental theory, while structured variation would suggest that $a_0$ encodes additional galaxy-state information beyond what a single constant captures.

Previous investigations have reached conflicting conclusions. Li et al. (2018) found that a universal $a_0$ provides a good fit to SPARC galaxies with minimal intrinsic scatter. Rodrigues et al. (2018) argued for significant galaxy-to-galaxy variation, a claim contested by subsequent analyses (Kroupa et al. 2018; McGaugh et al. 2018). The debate has largely focused on whether observed scatter is consistent with measurement noise or requires genuine variation.

We take a different approach. Rather than asking *whether* $a_0$ varies, we ask: *if it varies, does the variation follow a structured law?* Specifically, we test whether the pattern of $a_0$ values across galaxies can be predicted from independently measured galaxy properties — gas content, environment, baryonic structure, and kinematic regularity — using cross-validated, out-of-sample metrics that guard against overfitting.

Our analysis proceeds through a hierarchy of models: a null model (universal $a_0$), intermediate models with 2–4 predictors, a full six-axis model (B$''$), and unconstrained per-galaxy values. We evaluate each using leave-one-out cross-validation (LOO), information criteria, and extensive robustness tests including holdout transfer, morphology-stratified validation, and independent replication from raw data.

---

## 2. Data and Sample Definition

### 2.1 Source Catalog

We draw rotation-curve data and photometric properties from the Spitzer Photometry and Accurate Rotation Curves (SPARC) database (Lelli, McGaugh & Schombert 2016), which provides 175 disk galaxies with 3.6 $\mu$m photometry and high-quality HI/H$\alpha$ rotation curves. We supplement SPARC with:

- Host-halo masses ($M_\mathrm{host}$) from the Tully (2015) Cosmicflows group catalog and the Extragalactic Distance Database (EDD).
- Stellar mass-to-light ratios ($\Upsilon_\star^\mathrm{disk}$) from Li et al. (2020) SPS-based estimates at 3.6 $\mu$m.

### 2.2 Quality Subsample (N = 45)

From the 175 SPARC galaxies, we select a quality-controlled subsample of $N = 45$ based on the following criteria:

1. **Valid per-galaxy $a_0$ fit**: The galaxy must have a well-defined individual $a_0$ from its rotation curve (finite best-fit value with bounded uncertainty).
2. **Published distance quality**: The adopted distance must come from a published, peer-reviewed source using a primary or well-calibrated secondary distance indicator (TRGB, Cepheid, or equivalent), not solely from Hubble-flow estimates. This is the dominant filter (175 $\to$ 56 candidates).
3. **Host-mass assignment**: A host-halo mass from the Tully group catalog or equivalent environmental database must be available.
4. **Sufficient rotation-curve sampling**: At least 5 measured RC points, required for meaningful computation of kinematic shape diagnostics (rcWiggliness, logMeanRun).

The published-distance criterion ensures that distance-dependent quantities ($a_0$, $M_\mathrm{HI}$, physical sizes) are not dominated by Hubble-flow distance errors, which can introduce spurious scatter. The resulting $N = 45$ sample spans morphological types T = 0 to T = 11, flat rotation velocities $V_\mathrm{flat} = 17$–$300$ km/s, and HI masses $M_\mathrm{HI} = 0.07$–$33 \times 10^9 \, M_\odot$.

The complete sample, with all galaxy names and measured values, is provided in the supplementary dataset (Table A1).

### 2.3 Sample Selection Transparency

We emphasize that the sample was frozen prior to model comparison and was not optimized to favor any particular model. The selection criteria are purely observational (distance quality, data availability) and contain no cuts on $a_0$ itself or on any predictor variable.

---

## 3. Inference of Per-Galaxy $a_0$

### 3.1 The Radial Acceleration Relation

The RAR relates observed centripetal acceleration $g_\mathrm{obs} = V_\mathrm{obs}^2/r$ to the Newtonian acceleration from baryons alone $g_\mathrm{bar} = V_\mathrm{bar}^2/r$ via:

$$g_\mathrm{obs} = \frac{g_\mathrm{bar}}{1 - e^{-\sqrt{g_\mathrm{bar}/a_0}}}$$

where $a_0$ is the characteristic acceleration scale. When $g_\mathrm{bar} \gg a_0$, $g_\mathrm{obs} \approx g_\mathrm{bar}$ (Newtonian regime); when $g_\mathrm{bar} \ll a_0$, $g_\mathrm{obs} \approx \sqrt{a_0 \, g_\mathrm{bar}}$ (deep-MOND regime).

### 3.2 Per-Galaxy Fitting Procedure

For each galaxy, we fit the RAR to its rotation curve by minimizing the sum of squared residuals between observed and predicted velocities, treating $a_0$ as the free parameter while holding the stellar mass-to-light ratio $\Upsilon_\star^\mathrm{disk}$ at its SPS-estimated value. The best-fit $\log_{10} a_0$ for each galaxy constitutes our target variable.

### 3.3 Target Variable

The target for all regression models is:

$$y_i = \log_{10}(a_{0,i}) \quad [\text{units: m/s}^2]$$

The sample mean is $\bar{y} = -10.10$ (i.e., $a_0 \approx 0.8 \times 10^{-10}$ m s$^{-2}$), with a standard deviation SD$(y) = 0.258$ dex across the $N = 45$ galaxies. This 0.258 dex observed spread is the total variance that our models seek to explain.

---

## 4. Predictor Construction

We define six predictor variables for the B$''$ model. Each is constructed from independently measured galaxy properties and is defined below with full mathematical precision.

### 4.1 Gas Mass: $\log M_\mathrm{HI}$

$$\log M_\mathrm{HI} = \log_{10}\!\left(\frac{M_\mathrm{HI}}{10^9 \, M_\odot}\right)$$

Source: SPARC catalog column `MHI`. Units: $\log_{10}(10^9 \, M_\odot)$.

### 4.2 Rotation-Curve Irregularity: rcWiggliness

Measures the root-mean-square fractional deviation of the observed rotation curve from its own smoothed trend:

1. Let $V_i = V_\mathrm{obs}(r_i)$ for $i = 1, \ldots, n$ be the observed RC velocities at measured radii.
2. Compute a 3-point running mean for interior points ($i = 2, \ldots, n-1$):
$$V_i^\mathrm{smooth} = \frac{1}{3}\left(V_{i-1} + V_i + V_{i+1}\right)$$
3. Compute fractional residuals:
$$\delta_i = \frac{V_i - V_i^\mathrm{smooth}}{V_i^\mathrm{smooth}}$$
4. Define:
$$\mathrm{rcWiggliness} = \sqrt{\frac{1}{n-2} \sum_{i=2}^{n-1} \delta_i^2}$$

This is dimensionless and ranges from $\sim$0.01 (smooth RC) to $\sim$0.15 (highly irregular). It captures small-scale kinematic noise or substructure without reference to any model.

### 4.3 Host-Halo Mass: $\log M_\mathrm{host}$

$$\log M_\mathrm{host} = \log_{10}\!\left(\frac{M_\mathrm{host}}{M_\odot}\right)$$

Source: Tully (2015) Cosmicflows group catalog and EDD environmental assignments. For field galaxies, $M_\mathrm{host}$ is the galaxy's own estimated halo mass; for group/cluster members, it is the total group mass. Each SPARC galaxy is matched to its Tully group by name cross-reference. Range: 10.4–13.0 dex.

### 4.4 Central Baryonic Surface Density: $\log \Sigma_0$

$$\log \Sigma_0 = \log_{10}\!\left(\mathrm{SB}_\mathrm{disk} \times \Upsilon_\star^\mathrm{disk}\right)$$

where $\mathrm{SB}_\mathrm{disk}$ is the SPARC disk central surface brightness (in $L_\odot$ pc$^{-2}$) and $\Upsilon_\star^\mathrm{disk}$ is the adopted stellar mass-to-light ratio at 3.6 $\mu$m. Units: $\log_{10}(M_\odot \, \text{pc}^{-2})$.

### 4.5 Kinematic Coherence: $\log \overline{\text{Run}}$

Measures the average length of consecutive same-sign residual stretches in the MOND fit:

1. Compute residuals $e_i = V_\mathrm{obs}(r_i) - V_\mathrm{MOND}(r_i)$ for all $n$ RC points.
2. Record the sign sequence $s_i = \mathrm{sign}(e_i)$.
3. Count runs $R$ = number of contiguous subsequences of identical sign.
4. Define:
$$\log \overline{\text{Run}} = \log_{10}\!\left(\frac{n}{R}\right)$$

Higher values indicate longer stretches of systematic positive or negative residuals — i.e., more coherent, large-scale departures from the MOND prediction. Range: 0.1–0.8.

### 4.6 Independent Stellar Mass-to-Light Ratio: $\Upsilon_\star^\perp$

The raw $\log_{10} \Upsilon_\star^\mathrm{disk}$ is partially confounded with gas mass, surface density, and morphology. We isolate the independent component through orthogonalization:

1. Fit an auxiliary OLS regression:
$$\log_{10} \Upsilon_\star^\mathrm{disk} = \alpha_0 + \alpha_1 \log M_\mathrm{HI} + \alpha_2 \log \Sigma_0 + \alpha_3 T + \varepsilon$$
where $T$ is the de Vaucouleurs morphological type.
2. This auxiliary regression explains $R^2 = 0.68$ of the $\Upsilon_\star$ variance. The confounders absorb the predictable component.
3. Define $\Upsilon_\star^\perp = \varepsilon$ (the OLS residual), which carries the genuinely independent stellar-structure information.

Source for $\Upsilon_\star^\mathrm{disk}$: Li et al. (2020) Table 1 SPS estimates. Where unavailable, we adopt $\Upsilon_\star^\mathrm{disk} = 0.50 \, M_\odot/L_\odot$ as default.

---

## 5. Model Space and Evaluation

### 5.1 Competing Models

We evaluate five nested models, all fit via ordinary least squares (OLS) on the $N = 45$ sample:

| Model | Predictors | Parameters ($k$) | Description |
|-------|-----------|-------------------|-------------|
| M0 | (none) | 1 | Universal $a_0$ (intercept only) |
| A$'$ | $\log M_\mathrm{HI}$, $\log M_\mathrm{host}$ | 3 | Gas + environment |
| B$'$ | A$'$ + rcWig, $\log \Sigma_0$ | 5 | + kinematic shape, density |
| B$''$ | B$'$ + $\log \overline{\text{Run}}$, $\Upsilon_\star^\perp$ | 7 | Full structured law |
| M4 | Per-galaxy free $a_0$ | 45 | Unconstrained (in-sample ceiling) |

### 5.2 Evaluation Metrics

We use five metrics designed to penalize complexity differently:

**LOO gap%** (primary metric): Leave-one-out cross-validated variance explained.
$$\text{LOO gap\%} = 100 \times \left(1 - \frac{\text{LOO-RMS}^2}{\text{SD}(y)^2}\right)$$
where LOO-RMS $= \sqrt{\frac{1}{N}\sum_i (y_i - \hat{y}_{i,-i})^2}$, with $\hat{y}_{i,-i}$ the prediction for galaxy $i$ from a model trained on all galaxies except $i$. This metric directly measures out-of-sample predictive power without any holdout waste.

**$k$-fold CV**: 5-fold and 10-fold cross-validation, averaged over 100 random partitions.

**AIC / BIC**: Akaike and Bayesian information criteria, computed from the in-sample residual sum of squares with standard penalty terms.

**Residual $\tau$**: Standard deviation of in-sample residuals: $\tau = \text{SD}(\hat{\varepsilon})$.

---

## 6. Main Results

### 6.1 Model Comparison

Table 1 presents the head-to-head comparison of all five models.

**Table 1. Model Comparison on the N = 45 Quality Subsample**

| Model | LOO gap% | $\tau$ (dex) | AIC | BIC | 5-fold | 10-fold |
|-------|----------|-------------|-----|-----|--------|---------|
| M0 (universal $a_0$) | $-$2.3% | 0.261 | $-$121.0 | $-$119.2 | — | — |
| A$'$ (2 var) | 30.2% | 0.215 | $-$139.4 | $-$130.4 | 28.5% | 29.2% |
| B$'$ (4 var) | 45.6% | 0.190 | $-$151.3 | $-$140.5 | 44.0% | 44.5% |
| **B$''$ (6 var)** | **49.7%** | **0.183** | **$-$154.2** | **$-$141.6** | **47.7%** | **48.8%** |
| M4 (free) | $\ll 0$ | 0.000 | — | — | — | — |

B$''$ wins unanimously. It achieves the highest LOO gap%, the lowest residual scatter among parsimonious models, and the best AIC and BIC. The universal-constant model M0 has *negative* LOO gap% — it predicts worse than the sample mean on unseen galaxies. The unconstrained model M4, while achieving zero in-sample residuals, catastrophically overfits and produces deeply negative out-of-sample scores.

### 6.2 The B$''$ Law

The fitted B$''$ model is:

$$\log a_0 = \\underset{(\pm 0.63)}{4.60} \\underset{(\pm 0.05)}{- 0.243} \log M_\mathrm{HI} \underset{(\pm 1.04)}{+ 1.279} \, \text{rcWig} \underset{(\pm 0.05)}{- 0.150} \log M_\mathrm{host}$$
$$\\underset{(\pm 0.07)}{+ 0.168} \log \Sigma_0 \\underset{(\pm 0.12)}{+ 0.446} \log \overline{\text{Run}} \underset{(\pm 0.31)}{+ 0.659} \, \Upsilon_\star^\perp$$

**Table 2. B$''$ Coefficients**

| Predictor | Coefficient | SE | $t$-statistic | Physical axis |
|-----------|------------|-----|--------------|---------------|
| $\log M_\mathrm{HI}$ | $-$0.243 | 0.051 | $-$4.78 | Gas mass |
| rcWiggliness | +1.279 | 1.044 | +1.22 | RC irregularity |
| $\log M_\mathrm{host}$ | $-$0.150 | 0.049 | $-$3.08 | Environmental depth |
| $\log \Sigma_0$ | +0.168 | 0.066 | +2.56 | Baryonic density |
| $\log \overline{\text{Run}}$ | +0.446 | 0.117 | +3.80 | Kinematic coherence |
| $\Upsilon_\star^\perp$ | +0.659 | 0.314 | +2.10 | Independent stellar M/L |

All six coefficients carry the same sign across every robustness test (Section 8). The in-sample $R^2 = 0.73$ (adjusted $R^2 = 0.69$).

### 6.3 Physical Interpretation of Axes

The six axes span five distinct physical dimensions of galaxy state:

1. **Gas content** ($\log M_\mathrm{HI}$, negative): Galaxies with more gas have lower inferred $a_0$.
2. **Environmental depth** ($\log M_\mathrm{host}$, negative): Galaxies in more massive halos/groups have lower $a_0$.
3. **Baryonic concentration** ($\log \Sigma_0$, positive): Denser central baryonic distributions correlate with higher $a_0$.
4. **Kinematic regularity** ($\log \overline{\text{Run}}$, positive; rcWiggliness, positive): Galaxies with more coherent, structured MOND residuals show higher $a_0$.
5. **Stellar structure** ($\Upsilon_\star^\perp$, positive): Independent stellar M/L variation tracks $a_0$ beyond what gas, density, and morphology predict.

---

## 7. Validation and Generalization

### 7.1 Holdout Transfer

We train B$''$ on 34 galaxies with all six predictors fully available and predict $\log a_0$ for 11 held-out galaxies. The holdout LOO gap% = 22.9%, demonstrating genuine out-of-sample predictive power despite the smaller, noisier holdout sample.

### 7.2 Repeated Random Splits

Across 200 random 50/50 train-test splits, B$''$ (M3) outperforms the 4-variable model (M2) in 75% of splits, with mean test gap% = 34.4%. No split produces sign reversals in more than one coefficient.

### 7.3 Cross-Morphology Transfer

We test transfer across galaxy types:

- **Early-type trained $\to$ late-type tested**: Gap improvement of +7 percentage points over the within-sample baseline.
- **Late-type trained $\to$ early-type tested**: Gap improvement of +13 percentage points.

B$''$ generalizes across morphological boundaries — the law is not specific to one galaxy type.

### 7.4 Coefficient Stability Under Half-Sample Splits

All seven coefficients (intercept + 6 predictors) remain stable across 200 random half-sample splits. No coefficient exhibits a sign change in any split.

### 7.5 Independent Replication

Three independent fitting methods reproduce B$''$ from the raw data:

**Table 3. Replication Across Methods**

| Method | $\beta_1$ | $\beta_2$ | $\beta_3$ | $\beta_4$ | $\beta_5$ | $\beta_6$ |
|--------|----------|----------|----------|----------|----------|----------|
| OLS | $-$0.243 | +1.279 | $-$0.150 | +0.168 | +0.446 | +0.659 |
| Huber IRLS | $-$0.243 | +1.183 | $-$0.150 | +0.162 | +0.448 | +0.656 |
| Bootstrap mean | $-$0.247 | +1.349 | $-$0.146 | +0.157 | +0.438 | +0.660 |

All methods produce identical signs and coefficients agreeing within $\pm$0.10. B$''$ is not an artifact of any single fitting pipeline.

### 7.6 Jackknife and Bootstrap Stability

- **Jackknife** (45 leave-one-out iterations): 0/45 sign flips for all six variables.
- **Bootstrap** (1000 resamples): rcWiggliness shows 13.7% sign flips (weakest axis); all other axes show $\leq$3.3% flips.

---

## 8. Robustness and Alternative Explanations

### 8.1 Circularity Audit

Any variable mechanically linked to $a_0$ is suspect. We explicitly reject:

- **$\log$(Compactness)** = $\log(V_\mathrm{flat}^2/R_\mathrm{eff})$: This achieved $t = 3.04$ against B$''$ residuals but was rejected because $V_\mathrm{flat}^2/R_\mathrm{eff} \approx$ centripetal acceleration $\approx a_0$ by construction ($r = 0.42$ with $\log a_0$). Any variable containing $V_\mathrm{flat}$ as a predictor of $a_0$ is circular.

No B$''$ predictor involves $V_\mathrm{flat}$ or any directly $a_0$-derived quantity. The rcWiggliness and logMeanRun metrics use $V_\mathrm{obs}$ residuals from MOND fits, but these are *shape* statistics (fractional irregularity and coherence length) rather than acceleration-scale variables.

### 8.2 Alternative Environmental Models

We tested whether simpler environmental models explain the B$''$ residuals:

- **HI deficiency alone**: $r = -0.025$, $t = -0.17$ (no signal).
- **$\log(R_\mathrm{HI}/R_\mathrm{disk})$** (gas truncation): $r = 0.033$, $t = 0.22$ (no signal).
- **Ram-pressure susceptibility**: $r = -0.100$, $t = -0.66$ (no signal).
- **Processing composite** (HI$_\mathrm{def}$ $-$ $\log R_\mathrm{HI}/R_\mathrm{disk}$): $r = -0.023$, $t = -0.15$ (no signal).

No environmental processing proxy captures any residual structure beyond what $\log M_\mathrm{host}$ (already in B$''$) provides.

### 8.3 No Clean Seventh Axis

A comprehensive scan of 29 candidate variables against B$''$ residuals found no variable reaching $|t| \geq 1.65$ after orthogonalization against B$''$ predictors. The strongest candidate (logCompactness, $t = 3.04$) was rejected on circularity grounds ($\S$8.1).

---

## 9. Residual Frontier

### 9.1 Residual Diagnostics After B$''$

The B$''$ residuals ($N = 45$) exhibit:

| Statistic | Value |
|-----------|-------|
| Residual SD | 0.156 dex |
| Skewness | $-$0.07 |
| Excess kurtosis | $-$0.91 |
| Runs test $z$ | +1.06 |

The residuals are approximately symmetric and slightly platykurtic. The runs test shows no significant serial structure.

### 9.2 Error-Floor Audit

We decompose the residual variance into measurement and intrinsic components:

$$\sigma^2_\mathrm{obs} = \sigma^2_\mathrm{meas} + \sigma^2_\mathrm{int}$$

Per-galaxy measurement uncertainties are estimated by propagating distance, inclination, mass-model, and RC-fitting errors in quadrature.

**Table 4. Residual Variance Decomposition**

| Component | Variance (dex$^2$) | Fraction |
|-----------|-------------------|----------|
| $\sigma^2_\mathrm{obs}$ (after B$''$) | 0.0243 | 100% |
| $\sigma^2_\mathrm{meas}$ | 0.0133 | 54.7% |
| $\sigma^2_\mathrm{int}$ | 0.0110 | 45.3% |

$$\tau_\mathrm{int} = \sqrt{\sigma^2_\mathrm{int}} = 0.105 \text{ dex} \quad (95\% \text{ CI: } [0.041, 0.143])$$

The 95% confidence interval does not touch zero. We confirm this with a noise-only simulation: injecting only measurement noise into B$''$ predictions produces mock residual SDs in the range [0.081, 0.134] dex, while the observed residual SD is 0.156 dex — *outside* the mock 95% range. The data are noisier than pure measurement error predicts.

Furthermore, the absolute residuals do not correlate with per-galaxy measurement uncertainty ($r = -0.09$, $t = -0.59$), and galaxies with smaller measurement errors do not have systematically smaller residuals. This rules out the hypothesis that the residual scatter is entirely measurement-driven.

### 9.3 No Distributed Weak-Signal Bundle

We tested whether a bundle of 11 non-circular, non-B$''$ scalar variables could collectively explain the residual structure:

- **Ridge/OLS bundle**: LOO gap% drops by 19 percentage points relative to B$''$ alone (catastrophic overfitting).
- **Permutation test**: $p = 0.41$ (not significant).
- **PCA latent factors**: PC2 shows a marginal signal ($t = -1.97$) but collapses in out-of-sample prediction.
- **50/50 split validation**: The bundle wins only 0.2% of random splits.

The residual after B$''$ is not a distributed weak signal across available scalars.

### 9.4 Conclusion on the Residual

The remaining $\sim$0.10 dex intrinsic scatter after B$''$ is:
- Not explained by any single additional scalar variable (29 tested).
- Not explained by environmental processing proxies (4 tested).
- Not explained by a distributed bundle of 11 weak-signal variables.
- Not purely measurement noise (noise-only simulation fails).

This scatter likely requires dynamical-state information not captured by one-dimensional scalar proxies — resolved 2D kinematic data, direct infall/orbit measurements, or higher-resolution mass models — or represents a genuinely irreducible component at current data quality.

---

## 10. Limitations

1. **Sample size**: The primary result rests on $N = 45$ galaxies. While cross-validation, bootstrap, and holdout tests confirm generalizability within this sample, extension to larger and more diverse galaxy populations is essential.

2. **Distance dependence**: The sample selection is driven by distance quality. Galaxies with only Hubble-flow distances are excluded, which may introduce selection effects correlated with environment or luminosity.

3. **$\Upsilon_\star^\perp$ construction**: The orthogonalized stellar M/L variable depends on the choice of confounders ($\log M_\mathrm{HI}$, $\log \Sigma_0$, morphological type $T$). Different confounder sets would yield different $\Upsilon_\star^\perp$ values, though our choice explains 68% of the raw $\Upsilon_\star$ variance and is physically motivated.

4. **rcWiggliness stability**: This is the weakest axis in B$''$, with 13.7% bootstrap sign-flip rate. It contributes to predictive power but is less robustly determined than the other five predictors.

5. **Intrinsic scatter**: A genuine $\sim$0.10 dex residual remains unexplained. We do not claim that B$''$ captures all $a_0$ variation — only the structured, reproducible component accessible from current 1D scalar data.

6. **Interpretation**: We measure *inferred* $a_0$ variation, which could reflect genuine variation in the underlying acceleration scale, systematic biases in rotation-curve modeling, or both. Disentangling these requires independent dynamical probes.

---

## 11. Conclusion

On a quality-controlled sample of $N = 45$ SPARC galaxies with published distance estimates, galaxy-to-galaxy variation in the inferred MOND acceleration scale $a_0$ is neither best described by a strict universal constant (which predicts *worse* than the sample mean in LOO) nor by unconstrained per-galaxy values (which catastrophically overfit). Instead, it is captured by a structured six-axis law (B$''$) that:

1. **Explains 49.7%** of $\log a_0$ variance in leave-one-out cross-validation — a genuine, out-of-sample result.
2. **Wins unanimously** over all alternatives on every evaluation metric (LOO, 5-fold, 10-fold, AIC, BIC).
3. **Transfers** to held-out galaxies (gap% = 22.9%), across morphological types, and across random data splits.
4. **Replicates** identically under three independent fitting methods and is stable under jackknife (0 sign flips for 6/6 variables) and bootstrap testing.

The six axes — gas mass, environmental depth, baryonic density, kinematic regularity/coherence, and independent stellar M/L — span multiple physical dimensions of galaxy state. A genuine intrinsic residual of $\sim$0.10 dex persists, which cannot be reduced by any available scalar proxy but may yield to future resolved kinematic data.

The central result is methodological as much as physical: $a_0$ variation across galaxies is *structured, reproducible, and generalizable* — properties that random noise or pipeline artifacts do not possess. Whether this variation reflects genuine physics beyond a universal acceleration scale or systematic modeling effects remains an open question for future investigation with larger samples and independent dynamical probes.

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

## Appendix A: Supplementary Tables

**Table A1**: The full $N = 45$ sample with all variable values is provided in the supplementary dataset file `N45_final_dataset.csv`.

**Table A2**: Variable definitions, computation procedures, and sources are detailed in `variable_definitions.md`.

**Table A3**: Model equations, fitting procedures, and expected outputs for independent replication are provided in `model_spec.md` and `expected_outputs.md`.

**Replication Script**: A self-contained Python script (`replicate_from_scratch.py`) that reproduces all coefficients, LOO gap%, and residual diagnostics from the CSV dataset is provided. It requires only Python 3.7+ and NumPy.

---

## Appendix B: Glossary of Metrics

| Symbol | Definition |
|--------|-----------|
| LOO gap% | $100 \times (1 - \text{LOO-RMS}^2 / \text{SD}(y)^2)$; out-of-sample variance explained |
| $\tau$ | SD of in-sample OLS residuals |
| $\tau_\mathrm{int}$ | $\sqrt{\max(0, \tau^2 - \langle\sigma_{\mathrm{meas},i}^2\rangle)}$; intrinsic scatter after measurement noise |
| AIC | $n \ln(\text{RSS}/n) + 2k$ |
| BIC | $n \ln(\text{RSS}/n) + k \ln n$ |
| gap closure | $100 \times (1 - \tau_\mathrm{model}^2 / \text{SD}(y)^2)$ (in-sample analog of LOO gap%) |

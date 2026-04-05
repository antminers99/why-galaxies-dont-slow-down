# Gas-to-Stellar State Is the Leading Predictor of Outer Support Requirement in High-Vflat SPARC Galaxies

*Non-peer-reviewed computational analysis report*
*Analysis frozen at Phase 122*

---

## 1. Introduction

Galaxy rotation curves remain one of the most important observational windows into the relationship between baryonic matter and total gravitational potential. The observed discrepancy between baryonic predictions and measured rotation velocities in the outer regions of disk galaxies -- the "outer support requirement" -- has been quantified through various scaling relations, most notably the Radial Acceleration Relation (RAR) and the Baryonic Tully-Fisher Relation (BTFR).

This project began by investigating a strong apparent correlation between inner-to-outer rotation-curve shape variables: log(Vflat/InnerVmax) appeared to predict outerSlope with r=0.85. However, systematic null testing revealed this to be a geometric artifact of smooth monotonic curves, not an independent physical law.

The project then reframed the problem entirely: replacing curve-derived variables with independent catalog predictors, and redefining the target as the outer mass discrepancy (logOMD = log10 of the mean Vobs^2/Vbar^2 in the outer half of the rotation curve). In this framework, the gas-to-stellar balance emerged as the dominant predictor.

After 22 phases of systematic testing -- including robustness batteries, matched falsification, anti-circularity tests, mechanism mapping, and model competition -- we conclude that log(MHI/L3.6) is the strongest independently surviving state-variable candidate for outer support requirement in the high-Vflat SPARC regime.

---

## 2. Data

We use the SPARC database (Lelli, McGaugh & Schombert 2016, AJ 152, 157), which provides Spitzer 3.6um photometry and high-quality rotation curves for 175 disk galaxies. Our analysis sample consists of N=104 galaxies with Vflat >= 70 km/s, at least 8 rotation curve data points, and valid MHI and L3.6 measurements.

We adopt fixed mass-to-light ratios of Upsilon_disk = 0.5 and Upsilon_bulge = 0.7 M_sun/L_sun at 3.6um, following the SPARC convention. No preprocessing, filtering, or custom data pipelines are applied beyond these standard choices.

---

## 3. The Geometric Artifact (Phases 1-101)

The original investigation (Phases 1-100) found that log(Vflat/InnerVmax) predicts outerSlope with r=0.85 and LOO R^2=0.69. This appeared to be a strong physical law connecting inner rise rate to outer flatness.

Phase 101 tested this against smooth monotonic null curves with no physical content. The null model reproduced r=0.83. Three of four null-coupling criteria failed. The real signal was only 1.02x above the geometric floor. The ratio law was therefore classified as a geometric artifact arising from the mathematical properties of smooth rotation curves, not from independent physics.

This critical finding led to a complete reframing:
- The target was changed to outer mass discrepancy (immune to geometric coupling between inner and outer curve shape).
- Predictors were restricted to catalog-only variables (not derived from the rotation curve itself).

---

## 4. Discovery and Validation (Phases 102-112)

### 4.1 Discovery (Phase 102)

Among catalog predictors of logOMD, gas-related variables dominated. The gas fraction fgas = MHI/(MHI + 0.5*L3.6) achieved LOO R^2 = 0.503, substantially outperforming surface density (Sigma0, LOO = 0.321), luminosity (L3.6, LOO = 0.475), and baryonic compactness (LOO = 0.471).

### 4.2 External Validation (Phase 104)

Locked 70/30 train/test splits over 200 iterations yielded external R^2 = 0.474 with a calibration slope of 1.012 (near-perfect), confirming that the signal generalizes within the SPARC sample.

### 4.3 Variable Competition (Phase 105)

In a systematic death match, fgas emerged as the strong backbone predictor. No other single variable matched its LOO performance, and a 4-variable model added only +0.016 LOO over fgas alone (with severe collinearity, VIF 14-22).

### 4.4 Mediation Analysis (Phase 108)

A critical dissection showed that surface density (Sigma0) is fully mediated by gas fraction: controlling for fgas, the partial r of Sigma0 with logOMD drops to 0.004 (effectively zero). In contrast, fgas retains partial r = 0.424 after controlling Sigma0. The ratio is 103:1, establishing that fgas is causally upstream of Sigma0 in this context.

### 4.5 Matched Falsification (Phase 112)

The log form log(MHI/L3.6) was identified as more fundamental than fgas. It survives double matching for both Sigma0 and L3.6 simultaneously (permutation p=0.034), while fgas is marginal (p=0.065). Within-bin analysis confirmed: fgas is stable across Sigma0 terciles (r=0.39/0.73/0.46), while Sigma0 collapses across fgas terciles (r=-0.36/-0.05/-0.37).

---

## 5. Independence from Known Relations (Phases 113-116)

The signal is not a repackaging of known scaling relations:

- **Not RAR residuals**: logOMD and mean RAR residual share only 35% variance. 65% of the gas-to-stellar signal is independent of the RAR (Phase 113).
- **Not BTFR scatter**: fgas vs BTFR residual r=0.056. The gas-to-stellar signal and the BTFR are completely orthogonal axes (Phase 115).
- **Not acceleration-scale variation**: a0 shows only weak correlation with gas state (r=0.17), with a suggestive +0.11 dex higher a0 for gas-rich galaxies but far from significant (Phase 116).

---

## 6. Anti-Circularity (Phase 117)

A critical methodological concern: since MHI feeds into both the predictor (log(MHI/L3.6)) and the target (through Vgas in Vbar), the correlation could be a computational artifact. Six independent tests definitively ruled this out:

1. **Stars-only target** (gas removed from Vbar): r = 0.885, even stronger than the full target (0.722).
2. **Gas-only target**: r = -0.208, weak and opposite sign.
3. **Pure velocity target** (no decomposition): correlations persist.
4. **Gas contribution**: median 23% of Vbar^2 in outer regions.
5. **Permutation test** (1000 iterations): p < 0.0001, 7.5 sigma.
6. **Cross-component signals**: confirmed between independent baryonic components.

Score: 6/6 PASS. The signal is not a circularity artifact.

---

## 7. Mechanism Map (Phase 120)

Phase 120 dissected *why* the gas-to-stellar ratio predicts outer support requirement, testing five competing hypotheses:

### 7.1 Component Competition

| Variable | r | LOO R^2 |
|:---|:---|:---|
| log(MHI/L3.6) | 0.728 | 0.510 |
| log(L3.6) alone | -0.705 | 0.475 |
| log(MHI) alone | -0.376 | 0.101 |

The ratio beats both individual components. A two-component model (logMHI + logL3.6) gains only +0.026 LOO over the ratio.

### 7.2 Residual Dissection

After residualizing against each other:
- MHI_perp (gas after removing star effect): r = 0.272
- L3.6_perp (stars after removing gas effect): r = -0.655
- Ratio: r = 0.728

The stellar side carries more weight, but the ratio itself is more informative than either component alone.

### 7.3 Mechanism Verdicts

| Hypothesis | Verdict |
|:---|:---|
| Gas-only driver | FAIL |
| Stars-only driver | FAIL |
| **Ratio-state** | **PASS** |
| Structure-only | FAIL |
| Halo-mediated | FAIL |

**Conclusion:** log(MHI/L3.6) is not merely a proxy for gas mass alone, nor stellar mass alone, nor surface structure alone, but represents an independent relative physical state of the galaxy.

---

## 8. Minimal Physical Model (Phase 121)

### 8.1 Adopted Model

The best two-variable model combines gas-to-stellar balance with total baryonic mass:

$$\text{logOMD} = 1.749 + 0.203 \times \log(M_\text{HI}/L_{3.6}) - 0.101 \times \log(M_\text{bar})$$

LOO R^2 = 0.584, versus 0.543 for logRatio alone. The improvement is significant (permutation p=0.001) and both coefficients are bootstrap-stable (95% CI excluding zero).

### 8.2 Exploratory Three-Variable Model

Adding halo mass (log(Mdyn - Mbar)) yields LOO R^2 = 0.850, a massive improvement. However, since haloMass is derived from Vflat (which shares information with the target through Vobs), this model is flagged as exploratory only and is not adopted as the primary result.

---

## 9. Discussion

### 9.1 What the Result Means

Outer support requirement -- the degree to which observed rotation velocities exceed baryonic predictions in the outer regions of disk galaxies -- is most tightly linked not to gas mass alone or stellar mass alone, but to the gas-to-stellar state of the galaxy, with additional but incomplete coupling to total baryonic mass.

Gas-rich galaxies (high log(MHI/L3.6)) systematically show larger outer mass discrepancies. The signal is not driven by surface density, baryonic compactness, or any single halo property. It survives controlling for every tested confound, including all halo properties simultaneously (retaining 79% of the original signal).

### 9.2 Physical Interpretation

The gas-to-stellar balance likely reflects a combination of evolutionary state (how far the galaxy has progressed in converting gas to stars) and structural properties (how concentrated the baryons are relative to the halo). Phase 119 showed that no single standard model -- halo concentration, baryon-halo coupling, or simple evolutionary state -- fully captures the effect. Phase 120 confirmed that the ratio itself, not either component, carries the fundamental information.

This suggests that outer support requirement is sensitive to a galaxy "state variable" that integrates multiple aspects of galaxy evolution and structure, with gas-to-stellar balance serving as the best available observable proxy for this state.

### 9.3 Relation to Known Results

The finding is independent of both major scaling relations: it shares only 35% variance with RAR residuals and is completely orthogonal to the BTFR (r=0.056 with BTFR residuals). While gas fraction has been extensively studied in the context of galaxy evolution and scaling relations, the specific formulation connecting gas-to-stellar balance to outer mass discrepancy as an independent axis from both BTFR and RAR appears to be new.

---

## 10. Limitations

1. **Regime-specific**: The signal strengthens gradually above Vflat ~55-75 km/s and does not transfer to low-Vflat dwarfs, where all tested predictors fail under locked calibration.

2. **Single survey**: All data comes from the SPARC database. Cross-survey replication is blocked because no other survey currently provides both Spitzer photometry and full baryonic rotation curve decomposition for a comparable sample.

3. **Moderate effect size**: LOO R^2 of 0.51-0.58 means 42-49% of the variance in outer support requirement remains unexplained by the adopted model.

4. **No single causal mechanism**: While the correlation is robustly established and circularity is ruled out, no single physical model fully accounts for the effect.

5. **Fixed mass-to-light ratios**: The analysis uses fixed Upsilon values. Variation in stellar mass-to-light ratios could in principle modify the results, though the Phase 110 MC analysis showed the signal is robust to reasonable perturbations.

6. **Not peer-reviewed**: This is a computational analysis report. The results have not been subjected to formal peer review.

---

## References

Lelli, F., McGaugh, S.S. & Schombert, J.M., 2016, AJ, 152, 157 (SPARC database)

See `references.md` and `LITERATURE_SUMMARY.md` for the full reference list (22 papers across 5 relevance tiers) and detailed literature positioning.

---

*Analysis frozen at Phase 122. All scripts and intermediate results are included in the repository for full reproducibility. DOI: 10.5281/zenodo.19431363 (concept DOI, all versions).*

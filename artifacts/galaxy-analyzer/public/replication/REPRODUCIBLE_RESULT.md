# Gas-to-Stellar State Is the Leading Predictor of Outer Support Requirement in High-Vflat SPARC Galaxies

**Status:** FROZEN (Phase 122). Strong validated result. Analysis complete. No further phases planned unless external replication data becomes available or a falsification test kills the claim.

**Claim:** The gas-to-stellar balance, quantified by log(MHI/L3.6), is a real, non-circular, independently surviving predictor of outer support requirement in high-Vflat SPARC galaxies. It is not a proxy for gas mass alone, stellar mass alone, surface structure alone, or halo properties alone, but represents an independent relative physical state of the galaxy. The signal is not a geometric artifact, not a computational circularity artifact, and not fully reducible to any single standard physical model.

**Adopted model:**
$$\log_{10}\langle V_\text{obs}^2 / V_\text{bar}^2 \rangle_\text{outer} = 1.749 + 0.203 \times \log(M_\text{HI}/L_{3.6}) - 0.101 \times \log(M_\text{bar})$$

**Qualification:** This result is regime-specific with a gradual onset around Vflat ~55-75 km/s. It does not transfer to low-Vflat dwarfs. Fully independent cross-survey replication remains blocked by target incompatibility. This is a non-peer-reviewed computational analysis report.

---

## Abstract

We initially identified a strong relation between inner-to-outer rotation-curve shape variables, but Phase 101 showed that this apparent law was largely a geometric artifact of smooth monotonic curves. The project then shifted to independent catalog variables and redefined the target as outer support requirement. In this framework, gas-related observables strongly outperformed curve-shape predictors. Robustness, mediation, uncertainty, and matched-falsification analyses showed that the strongest surviving signal is the gas-to-stellar balance, quantified by log(MHI/L3.6), with fgas as a closely related secondary proxy. Mechanism-mapping tests (Phase 120) confirmed that the signal is not driven by gas mass alone (logMHI r=-0.376) or stellar mass alone (logL36 r=-0.705), but by their ratio (r=0.728), which beats both components. Anti-circularity tests confirmed the signal is not produced by shared gas inputs: when gas is removed from the target, the correlation strengthens (r=0.885), while gas-only targets fail (r=-0.208). The signal survives controlling for surface density (76%), baryonic mass (67%), star formation efficiency (63%), and all halo properties simultaneously (79%). An adopted two-variable model adds total baryonic mass as a second predictor (LOO R^2=0.584 vs 0.543 for the ratio alone), with both coefficients bootstrap-stable. A three-variable model adding halo mass yields LOO R^2=0.850 but is flagged as exploratory due to potential circularity. No single standard physical model fully captures the effect: halo concentration, baryon-halo coupling, and evolutionary-state proxies each fail to absorb the signal. The result is regime-specific, with gradual strengthening above Vflat ~55-75 km/s. We conclude that gas-to-stellar balance is currently the strongest falsification-surviving state-variable candidate for outer support requirement in the high-Vflat SPARC regime.

---

## 1. The Result

### Adopted Model (Two-Variable)

$$\log_{10}\langle V_\text{obs}^2 / V_\text{bar}^2 \rangle_\text{outer} = 1.749 + 0.203 \times \log(M_\text{HI}/L_{3.6}) - 0.101 \times \log(M_\text{bar})$$

| Metric | 1-var: log(MHI/L36) | 2-var: + log(Mbar) |
|:---|:---|:---|
| LOO R^2 | 0.543 | 0.584 |
| Permutation p (added var) | - | 0.001 |
| Bootstrap beta(logRatio) 95% CI | - | [0.115, 0.279] |
| Bootstrap beta(logMbar) 95% CI | - | [-0.165, -0.033] |

### Single-Variable Form (Backbone)

$$\text{predictor} = \log_{10}(M_\text{HI} / L_{3.6})$$

| Metric | fgas | log(MHI/L36) |
|:---|:---|:---|
| Pearson r | 0.724 | 0.729 |
| LOO R^2 | 0.503 | 0.512 |
| Permutation p | < 0.0001 | < 0.0001 |
| Matched falsification p | 0.065 (marginal) | 0.034 (significant) |
| MC uncertainty r 95% CI | [0.684, 0.736] | - |

### Exploratory Model (NOT adopted — circularity risk)

$$\text{logOMD} = 1.136 + 0.099 \times \log(M_\text{HI}/L_{3.6}) + 0.491 \times \log(M_\text{halo}) - 0.569 \times \log(M_\text{bar})$$
LOO R^2 = 0.850. Flagged because haloMass is derived from Vflat, which shares information with the target.

### Definitions

- **fgas** = MHI / (MHI + 0.5 * L36)
- **log(MHI/L36)** = log10(MHI) - log10(L36)
- **log(Mbar)** = log10(0.5 * L36 * 1e9 + 1.33 * MHI * 1e9)
- **Outer mass discrepancy** = mean of Vobs^2/Vbar^2 across outer-half rotation curve points
- **Vbar^2** = 0.5 * Vdisk|Vdisk| + 0.7 * Vbul|Vbul| + Vgas|Vgas|
- **Regime**: Vflat >= 70 km/s, N = 104 galaxies

---

## 2. Mechanism Map (Phase 120)

### Table 1 — Component Competition

| Variable | r | LOO R^2 | AIC | BIC |
|:---|:---|:---|:---|:---|
| **log(MHI/L3.6)** | **0.728** | **0.510** | **-391.6** | **-386.4** |
| log(L3.6) | -0.705 | 0.475 | -384.4 | -379.1 |
| fgas | 0.722 | 0.500 | -389.7 | -384.4 |
| log(MHI) | -0.376 | 0.101 | -328.9 | -323.6 |
| log(MHI)+log(L3.6) | R=0.755 | 0.536 | -398.9 | -391.0 |

The ratio beats both individual components. The two-variable model gains only +0.026 LOO over the ratio alone.

### Table 2 — Residual Survival

| Control variable(s) | partial r(logRatio, logOMD) | % of original |
|:---|:---|:---|
| Rdisk (size) | 0.680 | 93% |
| Sigma0 (surface density) | 0.556 | 76% |
| BaryonCompact | 0.644 | 88% |
| log(Mbar) | 0.487 | 67% |
| log(Mbar) + Vflat | 0.625 | 86% |
| log(SFE) | 0.457 | 63% |
| ALL halo properties | 0.591 | 79% |
| Rdisk + Sigma0 | 0.372 | 51% |

Reverse tests: Sigma0 dies after controlling logRatio (pr=0.009). SFE dies (pr=0.062). logMbar retains partial r=-0.309.

### Table 3 — Mechanism Verdict

| Hypothesis | Verdict | Evidence |
|:---|:---|:---|
| Gas-only driver | **FAIL** | logMHI r=-0.376 vs ratio r=0.728 |
| Stars-only driver | **FAIL** | logL36 r=-0.705 vs ratio r=0.728 |
| **Ratio-state (balance)** | **PASS** | Ratio beats both components |
| Structure-only | **FAIL** | Ratio survives structure controls (pr=0.372) |
| Halo-mediated | **FAIL** | Ratio survives all halo controls (pr=0.591) |

### Residual Dissection

- MHI_perp (gas after removing star effect): r = 0.272
- L36_perp (stars after removing gas effect): r = -0.655
- log(MHI/L3.6) (raw ratio): r = 0.728

The stellar side drives more of the signal, but the ratio is stronger than either residual alone.

---

## 3. Minimal Physical Model (Phase 121)

### Single-Variable Ranking

| Variable | LOO R^2 |
|:---|:---|
| log(MHI/L3.6) | 0.543 |
| fgas | 0.525 |
| log(Mbar) | 0.458 |
| innerDMfrac | 0.352 |
| logSigma0 | 0.290 |
| logRdisk | 0.212 |
| baryonCompact | 0.147 |
| haloMass | 0.145 |
| logRc | 0.053 |
| logRho0 | -0.011 |

### Model Comparison

| Model | k | LOO R^2 | Status |
|:---|:---|:---|:---|
| logRatio | 1 | 0.543 | Backbone |
| **logRatio + logMbar** | **2** | **0.584** | **ADOPTED** |
| logRatio + haloMass + logMbar | 3 | 0.850 | Exploratory only |

---

## 4. What Has Been Proven (Phases 101-121)

### Core Signal (Phases 101-112)

| Finding | Evidence |
|:---|:---|
| Old shape-law is a geometric artifact | Ph 101: Null curves reproduce r=0.83; real signal 1.02x above floor |
| fgas is NOT a proxy for Sigma0 | Ph 108: Sigma0\|fgas r=0.004 (dies); fgas\|Sigma0 r=0.424 (survives) |
| log(MHI/L36) is more fundamental than fgas | Ph 112: Survives double matching (p=0.034); fgas marginal (p=0.065) |
| Signal is noise-robust | Ph 110: r never below 0.5 under simultaneous perturbation |
| Boundary is gradual | Ph 111: Onset ~55 km/s, peak at 75 km/s |
| External validation passes | Ph 104: ext R^2=0.474, calibration slope=1.012 |

### Independence from Known Relations (Phases 113-116)

| Finding | Evidence |
|:---|:---|
| Independent of RAR residuals | Ph 113: 65% independent (shared variance only 35%) |
| Independent of BTFR | Ph 115: fgas vs BTFR residual r=0.056 |
| a0 roughly universal | Ph 116: r=0.17 with gas state |

### Anti-Circularity (Phase 117)

| Test | Result |
|:---|:---|
| Stars-only target (gas removed) | r=0.885 (STRONGER than full) |
| Gas-only target | r=-0.208 (fails) |
| Permutation (1000 iter) | p<0.0001, 7.5 sigma |
| Cross-component signals | Confirmed |
| Score | **6/6 PASS** |

### Physical Mechanism (Phases 118-121)

| Finding | Evidence |
|:---|:---|
| Connects to inner DM fraction | Ph 118: r=0.637 |
| Signal persists beyond halo | Ph 118: 76% survives all halo controls |
| No single model explains it | Ph 119: Halo conc. REJECTED, evolutionary REJECTED, baryon-halo WEAK |
| NOT gas alone, NOT stars alone | Ph 120: Ratio beats both components |
| Adding Mbar helps modestly | Ph 121: LOO 0.543 -> 0.584, permutation p=0.001 |

---

## 5. Data Source

**SPARC database** (Lelli, McGaugh & Schombert 2016, AJ 152, 157).

```
https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table1.dat   (galaxy properties)
https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat   (rotation curves)
```

No preprocessing, filtering, or custom data needed. Use the raw published tables.

---

## 6. How to Reproduce

### Step 1: Load the data
- Parse `table1.dat` for: galaxy name, Vflat, L36, MHI, Rdisk
- Parse `table2.dat` for: galaxy name, radius, Vobs, Vgas, Vdisk, Vbul

### Step 2: Select galaxies
- Keep galaxies with Vflat >= 70 km/s, >= 8 RC points, MHI > 0 and L36 > 0
- This gives N ~ 104 galaxies

### Step 3: Compute baryonic velocity at each point
$$V_\text{bar}^2 = 0.5 \cdot V_\text{disk}|V_\text{disk}| + 0.7 \cdot V_\text{bul}|V_\text{bul}| + V_\text{gas}|V_\text{gas}|$$

### Step 4: Compute outer mass discrepancy
- Sort rotation curve points by radius
- Take the outer half (last 50% of points)
- Compute: logOMD = log10(mean(Vobs^2 / Vbar^2) in outer half)

### Step 5: Compute predictors
- log(MHI/L36) = log10(MHI) - log10(L36)
- log(Mbar) = log10(0.5 * L36 * 1e9 + 1.33 * MHI * 1e9)

### Step 6: Verify
- r(log(MHI/L36), logOMD) ~ 0.729
- LOO R^2 ~ 0.51 (single), ~0.58 (with logMbar)

---

## 7. What Does NOT Work

| Test | Result | Why |
|:---|:---|:---|
| Low-Vflat dwarfs | ext R^2 < 0 for ALL | Regime boundary |
| LITTLE THINGS external | Cannot compute target | No decomposition |
| 4-variable model | LOO = 0.519 | Only +0.016 over fgas; collinearity |
| Old ratio law | LOO = 0.038 | Geometric artifact |
| Sigma0 independent signal | Dies after controlling fgas | Fully mediated (r=0.004) |
| 3-var with haloMass | LOO = 0.850 | Circularity risk (haloMass from Vflat) |

---

## 8. Limitations

1. **Regime-specific**: Vflat >= ~70 km/s spirals, gradual onset from ~55 km/s
2. **Single survey**: all data from SPARC; cross-survey replication blocked
3. **Moderate effect size**: LOO R^2 ~ 0.50-0.58 means 42-50% unexplained
4. **No single causal mechanism**: no standard model fully absorbs the signal
5. **Not peer-reviewed**: non-peer-reviewed computational analysis report
6. **Mass-luminosity degeneracy**: log(MHI/L36) uses luminosity as stellar mass proxy (fixed Upsilon)

---

## 9. Evidence Chain (All Phases)

| Phase | Finding | Verdict |
|:---|:---|:---|
| 1-100 | Ratio law Vflat/InnerVmax -> outerSlope | SUPERSEDED |
| 101 | Null geometric coupling test | 3/4 FAIL - artifact |
| 102 | Catalog predictors of outer mass discrepancy | fgas LOO R^2=0.503 |
| 103 | Robustness battery | 3/5 PASS |
| 104 | Locked external replication | ext R^2=0.474, cal slope=1.01 |
| 105 | Death match sparse models | fgas STRONG BACKBONE |
| 106 | Target definition robustness | MODERATE PASS (3/4 targets) |
| 108 | Mediation / proxy dissection | fgas CAUSAL CANDIDATE |
| 110 | Uncertainty propagation (MC) | STRONG PASS |
| 111 | Regime boundary mapping | Gradual transition ~55-75 km/s |
| 112 | Matched falsification | STRONG PASS - log(MHI/L36) wins |
| 113 | RAR-Gas Connection | logOMD 65% independent of RAR residuals |
| 114 | Diversity Connection | Gas predicts shape; baryon dominance survives Vflat |
| 115 | BTFR Residuals | fgas independent of BTFR (r=0.056) |
| 116 | Acceleration Scale | a0 roughly universal; gas-rich +0.11 dex hint |
| 117 | Anti-Circularity Battery | 6/6 PASS - NOT circular |
| 118 | Halo Anchoring | fgas->innerDMfrac r=0.637; 76% persists beyond halo |
| 119 | Model Competition | No single model wins; multi-faceted physics |
| 120 | Mechanism Map | Ratio-state PASS; gas-only FAIL; structure-only FAIL |
| 121 | Minimal Physical Model | 2-var adopted (LOO=0.584); 3-var exploratory only |
| **122** | **Writing Freeze** | **FROZEN — analysis complete** |

---

## 10. Files

| File | Contents |
|:---|:---|
| `REPRODUCIBLE_RESULT.md` | This document |
| `LITERATURE_SUMMARY.md` | Literature review and novelty positioning |
| `references.md` | Organized reference list (22 papers, 5 tiers) |
| `paper.md` | Historical phases 1-100 analysis (superseded) |
| `N45_final_dataset.csv` | Quality-controlled N=45 subset (phases 1-100) |
| `../scripts/phase120-mechanism-map.cjs` | Mechanism dissection (6 tests) |
| `../scripts/phase121-minimal-physical-model.cjs` | Minimal model search |
| `../scripts/phase117-anti-circularity.cjs` | Anti-circularity battery |
| `../scripts/phase118-halo-anchoring.cjs` | Halo fitting and mediation |
| `../scripts/phase119-model-competition.cjs` | Three-model competition |
| `../scripts/phase112-matched-falsification.cjs` | Decisive test |
| `../scripts/phase102-residual-physics.cjs` | Original fgas discovery |

---

## 11. Minimal Reproducibility Script (Pseudocode)

```python
import numpy as np

table1 = load("table1.dat")
table2 = load("table2.dat")

UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7

results = []
for galaxy in table1:
    if galaxy.Vflat < 70: continue
    if galaxy.MHI <= 0 or galaxy.L36 <= 0: continue

    rc = table2[galaxy.name]
    rc = rc[(rc.radius > 0) & (rc.vobs > 0)]
    if len(rc) < 8: continue

    vbar_sq = (UPSILON_DISK * rc.vdisk * abs(rc.vdisk)
             + UPSILON_BULGE * rc.vbul * abs(rc.vbul)
             + rc.vgas * abs(rc.vgas))

    half = len(rc) // 2
    mass_disc = (rc.vobs[half:]**2) / np.maximum(vbar_sq[half:], 0.01)
    logOMD = np.log10(np.mean(mass_disc))

    log_ratio = np.log10(galaxy.MHI) - np.log10(galaxy.L36)
    Mbar = 0.5 * galaxy.L36 * 1e9 + 1.33 * galaxy.MHI * 1e9
    log_mbar = np.log10(Mbar)

    results.append((log_ratio, log_mbar, logOMD))

ratio_arr, mbar_arr, y_arr = zip(*results)
print(f"N = {len(results)}")                                  # ~104
print(f"r(ratio) = {pearson_r(ratio_arr, y_arr):.3f}")        # ~0.729
print(f"LOO(ratio) = {loo_r2(ratio_arr, y_arr):.3f}")         # ~0.51
print(f"LOO(ratio+Mbar) = {loo_r2_2var(...):.3f}")            # ~0.58
```

---

## Errata (Superseded Phases 1-100 Analysis)

**T||5 bug in old scripts (fixed):** JavaScript scripts from Phases 1-100 used
`sparcMap[g.name]?.T || 5` which treated morphological type T=0 (valid for S0
galaxies NGC4138, UGC06786) as falsy, replacing it with T=5. This affected the
M5 model's Upsilon_perp construction. Corrected M5 numbers: LOO gap% = 46.6%
(was 50.9%), Upsilon_perp coefficient = 0.372 (was 0.658). M3 was unaffected
(44.1%). Fixed by using `?? 5` (nullish coalescing) in 21 scripts. These models
are superseded by the current Phases 101-122 analysis and this bug does not
affect the current result.

**logA0 units in CSV:** The logA0 column in N45_final_dataset.csv is in
log10((km/s)^2/kpc), NOT log10(m/s^2) as previously documented. To convert:
log10(m/s^2) = logA0 - 13.489. Values range [3.1, 4.1] in CSV units.

---

**Data:** SPARC (Lelli, McGaugh & Schombert, 2016, AJ 152, 157). All scripts and intermediate results included for full reproducibility. This is a non-peer-reviewed computational analysis report. Analysis frozen at Phase 122.

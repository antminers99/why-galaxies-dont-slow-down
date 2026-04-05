# Gas-to-Stellar Balance as the Strongest Independent Predictor of Outer Support Requirement in High-Vflat SPARC Galaxies

**Status:** Current backbone result — log(MHI/L36) has superseded fgas as the strongest falsification-surviving independent predictor.

**Claim:** The dominant independent signal driving outer support requirement in high-Vflat SPARC galaxies is gas-to-stellar balance, especially log(MHI/L36). This is not a reflection of surface density or baryonic compactness.

**Qualification:** This result is regime-specific with a gradual onset around Vflat ~55-75 km/s. It does not transfer to low-Vflat dwarfs, where all tested predictors fail under locked calibration. Fully independent cross-survey replication remains blocked by target incompatibility.

---

## 1. The Result

### Primary Model (fgas form)

$$\log_{10}\langle V_\text{obs}^2 / V_\text{bar}^2 \rangle_\text{outer} \approx 0.352 + 0.706 \times f_\text{gas}$$

### More Fundamental Form (ratio)

$$\text{predictor} = \log_{10}(M_\text{HI} / L_{3.6})$$

| Metric | fgas | log(MHI/L36) |
|:---|:---|:---|
| Pearson r | 0.724 | 0.729 |
| LOO R^2 | 0.503 | 0.512 |
| Permutation p | < 0.0001 | < 0.0001 |
| Matched falsification p (Sig0+L36) | 0.065 (marginal) | 0.034 (significant) |
| MC uncertainty r 95% CI | [0.684, 0.736] | - |

### Definitions

- **fgas** = MHI / (MHI + 0.5 * L36)
- **log(MHI/L36)** = log10(MHI) - log10(L36)
- **Outer mass discrepancy** = mean of Vobs^2/Vbar^2 across outer-half rotation curve points
- **Vbar^2** = 0.5 * Vdisk|Vdisk| + 0.7 * Vbul|Vbul| + Vgas|Vgas|
- **Regime**: Vflat >= 70 km/s, N = 104 galaxies

---

## 2. What Has Been Proven (Phases 105-112)

| Finding | Evidence |
|:---|:---|
| fgas is NOT a proxy for Sigma0 | Ph 108: Sigma0\|fgas r=0.004 (dies); fgas\|Sigma0 r=0.424 (survives). 103x ratio. |
| log(MHI/L36) is more fundamental | Ph 112: Survives double matching (p=0.034); fgas marginal (p=0.065) |
| Sigma0 has NO independent signal | Ph 112: Fails shuffle within L36 bins (p=0.109) and compact bins (p=0.348) |
| Signal is noise-robust | Ph 110: r never below 0.5 under simultaneous perturbation of all measurements |
| Boundary is gradual | Ph 111: Onset ~55 km/s, peak at 75 km/s, no sharp cut |
| fgas wins all 1-var comparisons | Ph 105: LOO=0.503 vs Sigma0=0.321, L36=0.475, compact=0.471 |
| External validation passes | Ph 104: ext R^2=0.474, calibration slope=1.012 (near-perfect) |
| Signal stable within Sigma0 bins | Ph 112: fgas r=0.39/0.73/0.46 across terciles (stable); Sigma0 r=-0.36/-0.05/-0.37 (collapses) |

---

## 3. Data Source

**SPARC database** (Lelli, McGaugh & Schombert 2016, AJ 152, 157).

Download from CDS/VizieR:
```
https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table1.dat   (galaxy properties)
https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat   (rotation curves)
```

No preprocessing, filtering, or custom data needed. Use the raw published tables.

---

## 4. How to Reproduce (Step by Step)

### Step 1: Load the data
- Parse `table1.dat` for: galaxy name, Vflat, L36, MHI, Rdisk
- Parse `table2.dat` for: galaxy name, radius, Vobs, Vgas, Vdisk, Vbul

### Step 2: Select galaxies
- Keep galaxies with Vflat >= 70 km/s
- Keep galaxies with >= 8 rotation curve data points
- Require MHI > 0 and L36 > 0
- This gives N ~ 104 galaxies

### Step 3: Compute baryonic velocity at each point
$$V_\text{bar}^2 = 0.5 \cdot V_\text{disk}|V_\text{disk}| + 0.7 \cdot V_\text{bul}|V_\text{bul}| + V_\text{gas}|V_\text{gas}|$$

### Step 4: Compute outer mass discrepancy
- Sort rotation curve points by radius
- Take the outer half (last 50% of points)
- Compute: logOMD = log10(mean(Vobs^2 / Vbar^2) in outer half)

### Step 5: Compute predictors
- fgas = MHI / (MHI + 0.5 * L36)
- log(MHI/L36) = log10(MHI) - log10(L36)

### Step 6: Correlate
- r(fgas, logOMD) ~ 0.724
- r(log(MHI/L36), logOMD) ~ 0.729

### Step 7: Verify with LOO cross-validation
- LOO R^2 ~ 0.50 (fgas), ~0.51 (log(MHI/L36))

---

## 5. Expected Results

| Metric | fgas | log(MHI/L36) |
|:---|:---|:---|
| Pearson r | 0.72-0.73 | 0.72-0.73 |
| LOO R^2 | 0.49-0.51 | 0.50-0.52 |
| Slope (fgas model) | ~0.71 | - |
| Intercept (fgas model) | ~0.35 | - |
| Permutation p | < 0.0001 | < 0.0001 |

---

## 6. What Does NOT Work

| Test | Result | Why |
|:---|:---|:---|
| Low-Vflat dwarfs (locked) | ext R^2 < 0 for ALL predictors | Regime boundary - different physics |
| LITTLE THINGS external | Cannot compute target | No baryonic decomposition available |
| 4-variable model | LOO R^2 = 0.519 | Only +0.016 over fgas alone; severe collinearity |
| Old Vflat/InnerVmax ratio law | LOO R^2 = 0.038 on new target | Geometric artifact (Phase 101) |
| At 2xRdisk target | logSigma0 wins (LOO=0.697 vs fgas 0.474) | Inner boundary is surface-density-sensitive |
| Sigma0 independent signal | Fails after controlling fgas | Fully mediated (Phase 108: r=0.004) |

---

## 7. Historical Context: The Geometric Artifact

### Phases 1-100: The Ratio Law (Superseded)

The original investigation found that log(Vflat/InnerVmax) predicts outerSlope with r=0.85, LOO R^2=0.69.

### Phase 101: Null Geometric Coupling Test

A null test using smooth monotonic curves with no physics reproduced r=0.83. Three of four null-coupling criteria failed. The real signal was only 1.02x above the geometric floor. The ratio law is a geometric artifact.

### Phase 102: Reframing

The target was changed to outer mass discrepancy (immune to geometric coupling) and predictors were restricted to catalog-only variables.

---

## 8. Physical Interpretation

The outer support requirement correlates fundamentally with the galaxy's gas-to-stellar state:

- **High gas-to-stellar ratio** (gas-dominant): large outer mass discrepancy, more "extra" support needed
- **Low gas-to-stellar ratio** (stellar-dominant): smaller outer mass discrepancy

The signal is NOT driven by surface density (Sigma0 dies after controlling for gas fraction) or baryonic compactness. The ratio form log(MHI/L36) outperforms fgas, suggesting the balance between gas reservoir and stellar content is the physically meaningful variable. This points toward evolutionary state / baryon conversion efficiency as the underlying driver.

---

## 9. Caveats

1. **Regime-specific**: works for Vflat >= ~70 km/s spirals, with gradual onset from ~55 km/s
2. **Single survey**: all data from SPARC; cross-survey replication blocked by data incompatibility
3. **Moderate effect size**: LOO R^2 ~ 0.50 means 50% of variance remains unexplained
4. **No causal mechanism**: the correlation is established; the physical mechanism is not
5. **fgas vs log(MHI/L36)**: these are near-equivalent forms; log(MHI/L36) slightly outperforms but both encode gas-to-stellar balance

---

## 10. Evidence Chain (All Phases)

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

---

## 11. Files in This Repository

| File | Contents |
|:---|:---|
| `REPRODUCIBLE_RESULT.md` | This document - current claim and replication guide |
| `paper.md` | Full phases 1-100 analysis (historical, superseded) |
| `N45_final_dataset.csv` | Quality-controlled N=45 subset (phases 1-100) |
| `../scripts/phase112-matched-falsification.cjs` | Decisive test script |
| `../scripts/phase110-uncertainty-propagation.cjs` | Monte Carlo uncertainty |
| `../scripts/phase111-boundary-mapping.cjs` | Regime boundary analysis |
| `../scripts/phase108-mediation-dissection.cjs` | Mediation / proxy test |
| `../scripts/phase105-death-match-sparse.cjs` | Model competition |
| `../scripts/phase106-target-robustness.cjs` | Target definition test |
| `../scripts/phase104-external-replication.cjs` | Locked external validation |
| `../scripts/phase102-residual-physics.cjs` | fgas discovery |

---

## 12. Minimal Reproducibility Script (Pseudocode)

```python
import numpy as np

table1 = load("table1.dat")  # galaxy properties
table2 = load("table2.dat")  # rotation curves

UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7

results = []
for galaxy in table1:
    if galaxy.Vflat < 70: continue
    if galaxy.MHI <= 0 or galaxy.L36 <= 0: continue

    rc = table2[galaxy.name]
    rc = rc[(rc.radius > 0) & (rc.vobs > 0)]
    if len(rc) < 8: continue

    # Baryonic velocity
    vbar_sq = (UPSILON_DISK * rc.vdisk * abs(rc.vdisk)
             + UPSILON_BULGE * rc.vbul * abs(rc.vbul)
             + rc.vgas * abs(rc.vgas))
    vbar = np.sqrt(np.maximum(vbar_sq, 0.01))

    # Outer half
    half = len(rc) // 2
    outer_vobs = rc.vobs[half:]
    outer_vbar = vbar[half:]

    # Outer mass discrepancy
    mass_disc = (outer_vobs**2) / (outer_vbar**2)
    logOMD = np.log10(np.mean(mass_disc))

    # Gas fraction
    fgas = galaxy.MHI / (galaxy.MHI + UPSILON_DISK * galaxy.L36)

    # Gas-to-stellar ratio (more fundamental)
    log_mhi_l36 = np.log10(galaxy.MHI) - np.log10(galaxy.L36)

    results.append((fgas, log_mhi_l36, logOMD))

fgas_arr, ratio_arr, y_arr = zip(*results)
print(f"r(fgas) = {pearson_r(fgas_arr, y_arr):.3f}")      # expect ~0.724
print(f"r(ratio) = {pearson_r(ratio_arr, y_arr):.3f}")     # expect ~0.729
print(f"N = {len(results)}")                                # expect ~104
```

---

**Data:** SPARC (Lelli, McGaugh & Schombert, 2016, AJ 152, 157). All scripts and intermediate results are included for full reproducibility.

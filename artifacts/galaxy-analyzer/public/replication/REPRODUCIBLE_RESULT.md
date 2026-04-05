# Gas Fraction as a Robust High-Regime Predictor of Outer Support Requirement in SPARC Galaxies

**Status:** Robust high-regime result

**Claim:** $f_\text{gas}$ is the strongest transferable single catalog predictor of outer support requirement in high-$V_\text{flat}$ SPARC galaxies.

**Qualification:** This result is regime-specific and does not transfer to low-$V_\text{flat}$ dwarfs, where all tested predictors fail under locked calibration. Fully independent cross-survey replication remains blocked by target incompatibility (LITTLE THINGS lacks baryonic decomposition).

---

## 1. The Result (One Equation)

$$\log_{10}\langle V_\text{obs}^2 / V_\text{bar}^2 \rangle_\text{outer} \approx 0.352 + 0.706 \times f_\text{gas}$$

Where:
- **$f_\text{gas}$** = $M_\text{HI} / (M_\text{HI} + \Upsilon_\star L_{3.6})$ with $\Upsilon_\star = 0.5$
- **Outer mass discrepancy** = mean of $V_\text{obs}^2 / V_\text{bar}^2$ across outer-half rotation curve points
- **$V_\text{bar}^2$** = $0.5 \cdot V_\text{disk}|V_\text{disk}| + 0.7 \cdot V_\text{bul}|V_\text{bul}| + V_\text{gas}|V_\text{gas}|$

**Interpretation:** Galaxies with higher gas fractions (more gas-rich, less stellar-converted) show larger outer mass discrepancies. The need for additional outer support is tied to the diffuse, gas-rich baryonic state.

---

## 2. Data Source

**SPARC database** (Lelli, McGaugh & Schombert 2016, AJ 152, 157).

Download from CDS/VizieR:
```
https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table1.dat   (galaxy properties)
https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat   (rotation curves)
```

No preprocessing, filtering, or custom data needed. Use the raw published tables.

---

## 3. How to Reproduce (Step by Step)

### Step 1: Load the data
- Parse `table1.dat` for: galaxy name, $V_\text{flat}$, $L_{3.6}$, $M_\text{HI}$, $R_\text{disk}$
- Parse `table2.dat` for: galaxy name, radius, $V_\text{obs}$, $V_\text{gas}$, $V_\text{disk}$, $V_\text{bul}$

### Step 2: Select galaxies
- Keep galaxies with $V_\text{flat} \geq 70$ km/s
- Keep galaxies with $\geq 8$ rotation curve data points
- Require $M_\text{HI} > 0$ and $L_{3.6} > 0$
- This gives $N \approx 104$ galaxies

### Step 3: Compute baryonic velocity at each point
$$V_\text{bar}^2 = 0.5 \cdot V_\text{disk}|V_\text{disk}| + 0.7 \cdot V_\text{bul}|V_\text{bul}| + V_\text{gas}|V_\text{gas}|$$

### Step 4: Compute outer mass discrepancy
- Sort rotation curve points by radius
- Take the outer half (last 50% of points)
- Compute: $\text{logOMD} = \log_{10}\left(\text{mean}\left(\frac{V_\text{obs}^2}{V_\text{bar}^2}\right)_\text{outer}\right)$

### Step 5: Compute gas fraction
$$f_\text{gas} = \frac{M_\text{HI}}{M_\text{HI} + 0.5 \times L_{3.6}}$$

### Step 6: Correlate
$$r(f_\text{gas}, \text{logOMD}) \approx 0.72, \quad p < 0.0001$$

### Step 7: Verify with LOO cross-validation
- Fit $\text{logOMD} = a + b \cdot f_\text{gas}$ leaving one galaxy out each time
- LOO $R^2 \approx 0.50$

---

## 4. Expected Results

| Metric | Expected Value |
|:---|:---|
| Pearson $r$ | $0.72$–$0.73$ |
| LOO $R^2$ | $0.49$–$0.51$ |
| Slope $b$ | $\sim 0.71$ |
| Intercept $a$ | $\sim 0.35$ |
| Permutation $p$ (10000 trials) | $< 0.0001$ |

---

## 5. Locked External Validation (Must Also Pass)

### 5a. 70/30 random split (200 iterations, no refit)

| Metric | Expected |
|:---|:---|
| External $r$ | 0.72 [0.58, 0.82] |
| External $R^2$ | 0.47 [0.11, 0.64] |
| Calibration slope | $\sim 1.01$ |
| Calibration offset | $\sim -0.004$ |

Near-perfect calibration (slope $\approx 1$, offset $\approx 0$) confirms model is not overfit.

### 5b. Leave-P%-out cross-validation

| Holdout % | Expected ext $R^2$ |
|:---|:---|
| 10% | $\sim 0.38$ |
| 20% | $\sim 0.45$ |
| 30% | $\sim 0.48$ |
| 40% | $\sim 0.48$ |
| 50% | $\sim 0.48$ |

Stable across all holdout sizes.

### 5c. Competitor knockout

$f_\text{gas}$ wins best external $R^2$ in $\sim 55\%$ of 70/30 splits, beating logSigma0, logBaryonCompact, and logL36.

---

## 6. What Does NOT Work

| Test | Result | Why |
|:---|:---|:---|
| Low-$V_\text{flat}$ dwarfs (locked) | ext $R^2 < 0$ for ALL predictors | Regime boundary — dwarfs are a different system |
| LITTLE THINGS external | Cannot compute target | No baryonic decomposition available |
| 4-variable model ($\Sigma_0$ + $M_\text{HI}$ + $R_\text{disk}$ + $f_\text{gas}$) | LOO $R^2 = 0.519$ | Only +0.016 over $f_\text{gas}$ alone; severe collinearity (VIF 14–22) |
| Old $V_\text{flat}/V_\text{max}^\text{inner}$ ratio law | LOO $R^2 = 0.038$ on new target | Geometric artifact (Phase 101) |

---

## 7. Historical Context: The Geometric Artifact

### Phases 1–100: The Ratio Law (Superseded)

The original investigation found that $\log(V_\text{flat}/V_\text{max}^\text{inner})$ predicts outerSlope with $r = 0.85$, LOO $R^2 = 0.69$. This appeared to be a strong physical law.

### Phase 101: Null Geometric Coupling Test

A null test using smooth monotonic curves with no physics reproduced $r = 0.83$. Three of four null-coupling criteria failed. The real signal was only 1.02$\times$ above the geometric floor. The ratio law is a **geometric artifact**: any monotonically rising curve will show this correlation by construction.

### Phase 102: Reframing

The target was changed to outer mass discrepancy (immune to geometric coupling) and predictors were restricted to catalog-only variables (not derived from rotation curve shape). $f_\text{gas}$ emerged as the backbone predictor.

---

## 8. Physical Interpretation

$f_\text{gas}$ encodes the baryonic state of the galaxy:

- **High $f_\text{gas}$** (gas-dominant): large outer mass discrepancy, more "extra" support needed
- **Low $f_\text{gas}$** (stellar-dominant): smaller outer mass discrepancy, baryons account for more of the observed rotation

This suggests the outer support requirement correlates with the degree of baryonic conversion — galaxies that have converted less gas into stars require proportionally more dark matter (or modified gravity) support in their outer regions.

**Open question:** Is $f_\text{gas}$ the causal driver, or a proxy for something deeper — low surface density state, recent evolutionary history, weak stellar feedback, or baryon-halo coupling efficiency?

---

## 9. Caveats

1. **Regime-specific**: works for $V_\text{flat} \geq 70$ km/s spirals only
2. **Single survey**: all data from SPARC; cross-survey replication blocked by data incompatibility
3. **Collinearity with $\Sigma_0$**: $r = -0.81$ between $f_\text{gas}$ and $\log\Sigma_0$, so causal attribution is ambiguous
4. **Moderate effect size**: LOO $R^2 = 0.50$ means 50% of variance remains unexplained
5. **No causal mechanism**: the *why* behind the correlation is not yet established

---

## 10. Files in This Repository

| File | Contents |
|:---|:---|
| `REPRODUCIBLE_RESULT.md` | This document — current claim and replication guide |
| `paper.md` | Full phases 1–100 analysis (historical, superseded by 101+) |
| `N45_final_dataset.csv` | Quality-controlled N=45 subset (phases 1–100) |
| `../phase104-external-replication.json` | Phase 104 raw results |
| `../phase103-robustness-battery.json` | Phase 103 raw results |
| `../phase102-residual-physics.json` | Phase 102 raw results |
| `../scripts/phase104-external-replication.cjs` | Locked external validation script |
| `../scripts/phase102-residual-physics.cjs` | fgas discovery script |

---

## 11. Minimal Reproducibility Script (Pseudocode)

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
    outer = rc[half:]
    outer_vbar = vbar[half:]

    # Outer mass discrepancy
    mass_disc = (outer.vobs**2) / (outer_vbar**2)
    logOMD = np.log10(np.mean(mass_disc))

    # Gas fraction
    fgas = galaxy.MHI / (galaxy.MHI + UPSILON_DISK * galaxy.L36)

    results.append((fgas, logOMD))

x, y = zip(*results)
print(f"r = {pearson_r(x, y):.3f}")     # expect ~0.72
print(f"N = {len(results)}")             # expect ~104
# LOO R^2 ~ 0.50
```

---

**Data:** SPARC (Lelli, McGaugh & Schombert 2016, AJ 152, 157). All scripts and intermediate results are included in this repository for full reproducibility.

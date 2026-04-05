# A Single-Variable Law Governing Rotation Curve Shape in Disk Galaxies

**Claim:** The ratio $V_\text{flat} / V_\text{max}^\text{inner}$ predicts the outer rotation curve slope of massive spiral galaxies with $r = 0.85$, $p < 0.0001$, surviving all measurement controls.

**Status:** Confirmed physical (not measurement artifact). Score: 8/8 physical, 0/8 measurement.

---

## 1. The Law (One Equation)

$$\text{outerSlope} \approx 1.92 \times \log_{10}\!\left(\frac{V_\text{flat}}{V_\text{max}^\text{inner}}\right) + 0.018$$

Where:
- **$V_\text{flat}$**: asymptotic flat rotation velocity (Table 1 of Lelli+ 2016)
- **$V_\text{max}^\text{inner}$**: maximum observed $V_\text{obs}$ in the inner half of the rotation curve
- **outerSlope**: slope of $\log V_\text{obs}$ vs $\log R$ in the outer half of the rotation curve

**Interpretation:** Galaxies that reach their peak velocity early (high $V_\text{max}^\text{inner}$ relative to $V_\text{flat}$) have flatter or declining outer curves. Galaxies still rising internally have steeper outer curves. This encodes **halo concentration** — the inside-to-outside mass distribution.

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
- Parse `table1.dat` for: galaxy name, $V_\text{flat}$, $R_\text{disk}$, distance $D$, quality flag $Q$
- Parse `table2.dat` for: galaxy name, radius, $V_\text{obs}$

### Step 2: Select galaxies
- Keep galaxies with $V_\text{flat} \geq 70$ km/s (high-regime spirals)
- Keep galaxies with $\geq 8$ rotation curve data points
- This gives $N \approx 104$ galaxies

### Step 3: Split each rotation curve in half
- Sort data points by radius
- **Inner half**: first 50% of points
- **Outer half**: last 50% of points

### Step 4: Compute two quantities
- **$V_\text{max}^\text{inner}$** = maximum $V_\text{obs}$ among inner-half points
- **outerSlope** = OLS slope of $\log_{10}(V_\text{obs})$ vs $\log_{10}(R)$ on outer-half points

### Step 5: Compute the ratio
$$x = \log_{10}\!\left(\frac{V_\text{flat}}{V_\text{max}^\text{inner}}\right)$$

### Step 6: Correlate
$$r(x, \text{outerSlope}) \approx 0.85, \quad p < 0.0001$$

### Step 7: Verify with LOO cross-validation
- Fit $\text{outerSlope} = a + b \cdot x$ leaving one galaxy out each time
- LOO $R^2 \approx 0.69$

---

## 4. Expected Results

| Metric | Expected Value |
|:---|:---|
| Pearson $r$ | $0.84$–$0.85$ |
| LOO $R^2$ (ratio only) | $0.68$–$0.70$ |
| Slope $b$ | $1.8$–$2.0$ |
| Intercept $a$ | $\approx 0.02$ |
| Permutation $p$ (2000 trials) | $< 0.001$ |

---

## 5. Controls That Must Also Pass

### 5a. Not a distance artifact
Split into 3 distance bins (tertiles). Expected:

| Bin | $r$ |
|:---|:---|
| Near ($D < 14$ Mpc) | $\sim 0.93$ |
| Mid ($14$–$21$ Mpc) | $\sim 0.78$ |
| Far ($> 21$ Mpc) | $\sim 0.83$ |

All $p < 0.0001$. If the law were a measurement artifact, it would vanish or weaken dramatically in at least one bin.

### 5b. Partial correlation after controlling D, nPts, radial coverage
$$r_\text{partial}(\text{ratio}, \text{outerSlope} \mid D, N_\text{pts}, R_\text{max}/R_d) \approx 0.81$$

Raw $r = 0.85$ drops by only $\sim 5\%$. Measurement properties explain almost nothing.

### 5c. Distance-matched pairs
Match galaxies within 20% distance. Among $\sim 51$ pairs:
$$r(\Delta\text{ratio}, \Delta\text{slope}) \approx 0.81, \quad p < 0.0001$$

### 5d. Alternative inner proxies
Replace $V_\text{max}^\text{inner}$ with less resolution-sensitive proxies:

| Proxy | Partial $r$ (after all controls) |
|:---|:---|
| Mean inner $V$ | $\sim 0.72$ |
| Deep inner Vmax (inner 25% of inner half) | $\sim 0.66$ |
| Inner baryon $V_\text{max}$ (from mass model) | $\sim 0.53$ |

All survive. The baryon-only proxy is independent of beam smearing.

---

## 6. What Does NOT Work

| Test | Result | Why |
|:---|:---|:---|
| Dwarf galaxies ($V_\text{flat} < 70$ km/s) | $r \approx 0.13$ | $V_\text{flat}$ is ill-defined for rising RCs |
| Per-galaxy $a_0$ as predictor | LOO $R^2 = -0.02$ | $a_0$ is a noisy proxy for this ratio |
| Linear fraction instead of log | Slope instability 64% | Log is the correct functional form |

---

## 7. Relationship to MOND $a_0$

The original investigation asked whether the MOND acceleration scale $a_0$ varies systematically between galaxies. The answer:

1. Per-galaxy $a_0$ does correlate with galaxy structure ($R^2 \approx 0.5$)
2. But $a_0$ is a **noisy proxy** for $V_\text{flat}/V_\text{max}^\text{inner}$
3. Adding $a_0$ to the ratio model contributes $\Delta R^2 = -0.009$ (negative)
4. Adding the ratio to the $a_0$ model contributes $\Delta R^2 = +0.39$
5. **$a_0$ is not needed.** The ratio is the fundamental variable.

---

## 8. Physical Interpretation

The ratio $V_\text{flat}/V_\text{max}^\text{inner}$ encodes **halo concentration**: how the total (baryonic + dark) mass is distributed between the inner and outer galaxy.

- **Ratio $\approx 1$** (high $V_\text{max}^\text{inner}$): concentrated mass, flat/declining outer RC
- **Ratio $> 1$** (low inner peak): extended halo dominates, rising outer RC

This is a **one-parameter family** of rotation curve shapes, governed by a single observable ratio. It predicts 69% of outer slope variance with zero overfitting (LOO), using one variable measured from the rotation curve itself.

---

## 9. Caveats

1. **Regime-specific**: works for $V_\text{flat} \geq 70$ km/s spirals only
2. **Calibration is distance-dependent**: the slope $b$ varies from $\sim 2.2$ (nearby) to $\sim 0.9$ (distant) due to resolution effects on $V_\text{max}^\text{inner}$
3. **Cross-sample transfer**: qualitative relationship is universal, but quantitative calibration depends on sample composition
4. **Not tested outside SPARC**: replication on independent datasets (e.g., THINGS, LITTLE THINGS) would strengthen the claim

---

## 10. Files in This Repository

| File | Contents |
|:---|:---|
| `paper.md` | Full 100-phase analysis (detailed) |
| `N45_final_dataset.csv` | Quality-controlled N=45 subset |
| `../phase100-physical-vs-measurement.json` | Phase 100 raw results |
| `../phase99-calibration-law.json` | Calibration analysis |
| `../phase98-external-replication.json` | Cross-sample replication |
| `../phase97-geometric-test.json` | Geometric artifact tests |
| `../scripts/phase100-physical-vs-measurement.cjs` | Phase 100 script (Node.js) |

---

## 11. Minimal Reproducibility Script (Pseudocode)

```python
# Load SPARC tables
table1 = load("table1.dat")  # galaxy properties
table2 = load("table2.dat")  # rotation curves

results = []
for galaxy in table1:
    if galaxy.Vflat < 70: continue
    
    rc = table2[galaxy.name]  # sorted by radius
    rc = rc[rc.radius > 0 & rc.vobs > 0]
    if len(rc) < 8: continue
    
    half = len(rc) // 2
    inner = rc[:half]
    outer = rc[half:]
    
    inner_vmax = max(inner.vobs)
    
    log_r = log10(outer.radius)
    log_v = log10(outer.vobs)
    outer_slope = OLS_slope(log_r, log_v)
    
    log_ratio = log10(galaxy.Vflat / inner_vmax)
    
    results.append((log_ratio, outer_slope))

x, y = zip(*results)
print(f"r = {pearson_r(x, y):.3f}")   # expect ~0.85
print(f"N = {len(results)}")           # expect ~104
```

---

**Contact / Citation:** This analysis uses publicly available SPARC data (Lelli, McGaugh & Schombert 2016). All scripts and intermediate results are included in this repository for full reproducibility.

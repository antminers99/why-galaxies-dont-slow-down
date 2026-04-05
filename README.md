# Galaxy Rotation Curve Analyzer

A data-driven investigation into the physical drivers of outer support requirement in galaxy rotation curves, using the SPARC database (Lelli, McGaugh & Schombert 2016).

## Status: Robust High-Regime Result

**Gas fraction is the strongest transferable catalog predictor of outer support requirement in high-$V_\text{flat}$ SPARC galaxies.**

### The Finding

Within the high-$V_\text{flat}$ ($\geq 70$ km/s) SPARC population, a single catalog variable — gas fraction $f_\text{gas} = M_\text{HI} / (M_\text{HI} + \Upsilon_\star L_{3.6})$ — predicts the outer mass discrepancy $\log_{10}\langle V_\text{obs}^2 / V_\text{bar}^2 \rangle_\text{outer}$ with:

- Pearson $|r| = 0.724$
- Leave-one-out $R^2 = 0.503$
- Permutation $p < 0.0001$

### Locked External Validation (Phase 104)

The model was locked on 70% training data and applied without refit to 30% held-out data across 200 random splits:

| Metric | Value |
|:---|:---|
| External $r$ | 0.723 [0.584, 0.820] |
| External $R^2$ | 0.474 [0.108, 0.640] |
| Calibration slope | 1.012 (ideal = 1.0) |
| Calibration offset | −0.004 (ideal = 0.0) |

Near-perfect calibration confirms this is not overfit.

### What This Means

Galaxies with higher gas fractions (more gas-rich, less stellar-converted) show larger outer mass discrepancies — they require more "extra" support beyond their baryonic content in the outer disk. This suggests that the need for additional outer support is tied to the diffuse, gas-rich baryonic state, not to the abstract inner shape of the rotation curve.

### What This Does NOT Claim

- This is **not** a universal law for all disk galaxies
- It **does not** transfer to low-$V_\text{flat}$ dwarfs (where all tested predictors fail equally under locked calibration)
- It has **not** been fully replicated on an independent external survey (LITTLE THINGS lacks baryonic decomposition for target variable matching)

## Evidence Chain

| Phase | Finding | Result |
|:---|:---|:---|
| 1–100 | Original $V_\text{flat}/V_\text{max}^\text{inner}$ ratio law | $r = 0.85$, LOO $R^2 = 0.69$ |
| **101** | **Null geometric coupling test** | **3/4 FAIL — ratio law is geometric artifact** |
| 102 | Reframed: catalog-only predictors of outer mass discrepancy | $f_\text{gas}$ LOO $R^2 = 0.503$ vs old ratio LOO $R^2 = 0.038$ |
| 103 | Robustness battery | 3/5 PASS (Collinearity FAIL due to $\Sigma_0$/$f_\text{gas}$ $r = -0.81$) |
| **104** | **Locked external replication** | **ext $R^2 = 0.474$, cal slope = 1.01** |

### The Critical Pivot (Phase 101)

The original ratio law ($V_\text{flat}/V_\text{max}^\text{inner}$ predicts outerSlope) appeared strong ($r = 0.85$) but Phase 101 demonstrated it is a **geometric artifact**: smooth monotonic curves with no physics reproduce $r = 0.83$. Three of four null-coupling tests failed. The real signal was only 1.02× above the geometric floor.

This led to a complete reframing:
- **New target**: $\log_{10}\langle V_\text{obs}^2 / V_\text{bar}^2 \rangle_\text{outer}$ (outer mass discrepancy — immune to geometric coupling)
- **New predictors**: catalog-only variables (no quantities derived from the rotation curve shape itself)
- $f_\text{gas}$ emerged as the backbone predictor

## Repository Structure

```
artifacts/galaxy-analyzer/
├── public/replication/          # Key deliverables
│   ├── REPRODUCIBLE_RESULT.md   # Current result summary and replication guide
│   ├── paper.md                 # Historical phases 1–100 analysis (superseded by 101+)
│   ├── N45_final_dataset.csv    # Quality-controlled N=45 dataset (phases 1–100)
│   └── replicate_from_scratch.py # Python replication script (phases 1–100)
├── scripts/                     # All analysis phase scripts (Node.js)
│   ├── phase104-external-replication.cjs   # Locked external validation
│   ├── phase103-robustness-battery.cjs     # Robustness tests
│   ├── phase102-residual-physics.cjs       # fgas discovery
│   ├── phase101-null-geometric-coupling.cjs # Artifact exposure
│   └── ... (phases 1–100)                  # Earlier investigation
├── public/                      # JSON results for each phase
│   ├── phase104-external-replication.json
│   ├── phase103-robustness-battery.json
│   └── ... (all phase results)
├── src/                         # React + Vite web application
│   ├── pages/                   # Interactive analysis pages
│   └── components/              # UI components
└── experiment-log.txt           # Chronological experiment log
```

## Quick Reproduction

### Current Result (Phases 102–104)

```bash
cd artifacts/galaxy-analyzer/
# Download SPARC data first:
# curl -o /tmp/sparc_table1.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table1.dat
# curl -o /tmp/sparc_cds.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat
node scripts/phase102-residual-physics.cjs     # fgas discovery
node scripts/phase103-robustness-battery.cjs   # robustness tests
node scripts/phase104-external-replication.cjs # locked external validation
```

### Historical Result (Phases 1–100, superseded)

See [`paper.md`](artifacts/galaxy-analyzer/public/replication/paper.md) for the full 100-phase investigation that led to and was ultimately superseded by the geometric artifact discovery.

## Data

All analysis uses the publicly available [SPARC database](http://astroweb.cwru.edu/SPARC/) (175 galaxies, 3600+ rotation curve points).

**Reference:** Lelli, F., McGaugh, S.S. & Schombert, J.M., 2016, AJ, 152, 157

## Technical Stack

- **Analysis**: Node.js scripts (CommonJS)
- **Web app**: React + TypeScript + Vite + Tailwind CSS + shadcn/ui
- **Statistics**: OLS regression, LOO cross-validation, permutation tests, bootstrap stability, partial correlations — all implemented from scratch (no external stats libraries)

## License

This is a research project using publicly available data. The analysis code and results are provided for scientific reproducibility.

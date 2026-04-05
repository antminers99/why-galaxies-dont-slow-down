# Galaxy Rotation Curve Analyzer

A data-driven investigation into the physical drivers of outer support requirement in galaxy rotation curves, using the SPARC database (Lelli, McGaugh & Schombert 2016).

## Main Result

**After ruling out the original curve-shape law as a geometric artifact, the project finds that gas-to-stellar balance, quantified by log(MHI/L36), is the strongest independent predictor of outer support requirement in high-Vflat SPARC galaxies.**

### Status

Current backbone result — log(MHI/L36) has superseded fgas as the strongest falsification-surviving independent predictor of outer support requirement in the high-Vflat SPARC regime.

### The Refined Claim (5 Statements)

1. **The old law fell.** The original Vflat/InnerVmax-to-outerSlope claim is not an independent physical law. It is reproduced almost entirely by geometric null models of smooth rotation curves (Phase 101).

2. **Moving to independent predictors was correct.** When the analysis moved from curve-derived variables to independent catalog variables, the physical signal became stronger and methodologically cleaner (Phase 102).

3. **The winner is NOT surface density.** Sigma0 and baryonic compactness do not retain comparable independent content once gas-related variables are controlled. Their apparent signal is largely mediated (Phases 105, 108, 112).

4. **The strongest variable is log(MHI/L36).** The strongest surviving independent predictor is log(MHI/L36), the gas-to-stellar balance. It survives matched falsification more cleanly than fgas alone (Phase 112: p=0.034 vs fgas p=0.065).

5. **The result is not universal.** The signal is regime-specific, with a gradual onset above roughly 55-75 km/s, and it fails to transfer to low-Vflat dwarfs under locked calibration (Phases 104, 111).

### Key Numbers

| Quantity | fgas | log(MHI/L36) |
|:---|:---|:---|
| Pearson r | 0.724 | 0.729 |
| LOO R^2 | 0.503 | 0.512 |
| Permutation p | < 0.0001 | < 0.0001 |
| External R^2 (70/30) | 0.474 | - |
| Matched falsification p | 0.065 | 0.034 |
| MC uncertainty r CI | [0.684, 0.736] | - |
| N | 104 | 104 |

## Evidence Chain

| Phase | Finding | Verdict |
|:---|:---|:---|
| 1-100 | Ratio law Vflat/InnerVmax -> outerSlope | SUPERSEDED |
| **101** | **Null geometric coupling test** | **3/4 FAIL - artifact** |
| 102 | Catalog predictors of outer mass discrepancy | fgas LOO R^2=0.503 |
| 103 | Robustness battery | 3/5 PASS |
| **104** | **Locked external replication** | **ext R^2=0.474, cal slope=1.01** |
| **105** | **Death match sparse models** | **fgas STRONG BACKBONE** |
| 106 | Target definition robustness | MODERATE PASS (3/4 targets) |
| **108** | **Mediation / proxy dissection** | **fgas CAUSAL CANDIDATE, Sigma0 fully mediated** |
| **110** | **Uncertainty propagation (MC)** | **STRONG PASS, r CI [0.684, 0.736]** |
| **111** | **Regime boundary mapping** | **Gradual transition ~55-75 km/s** |
| **112** | **Matched falsification** | **STRONG PASS - log(MHI/L36) is more fundamental** |

### The Critical Pivot (Phase 101)

The original ratio law (Vflat/InnerVmax predicts outerSlope) appeared strong (r=0.85) but Phase 101 demonstrated it is a geometric artifact: smooth monotonic curves with no physics reproduce r=0.83. Three of four null-coupling tests failed. The real signal was only 1.02x above the geometric floor.

This led to a complete reframing:
- **New target**: outer mass discrepancy (immune to geometric coupling)
- **New predictors**: catalog-only variables (not derived from rotation curve shape)
- **Winner**: gas-to-stellar balance (fgas / log(MHI/L36))

### The Decisive Test (Phase 112)

Matched falsification controlled for surface density and luminosity simultaneously:
- log(MHI/L36) survives double matching (p=0.034 significant)
- fgas survives marginally (p=0.065)
- Sigma0 fails shuffle tests within luminosity bins (p=0.109) and compactness bins (p=0.348)
- Within-bin analysis: fgas stable across all Sigma0 terciles (r=0.39/0.73/0.46); Sigma0 collapses (r=-0.36/-0.05/-0.37)

### Physical Interpretation

The outer support requirement (mass discrepancy) correlates fundamentally with the galaxy's gas-to-stellar state. Gas-rich galaxies show larger outer mass discrepancies. The signal is NOT driven by surface density (Sigma0 dies after controlling for gas fraction) or baryonic compactness. The ratio form log(MHI/L36) slightly outperforms fgas, suggesting the balance between gas reservoir and stellar content is the physically meaningful variable. This points toward evolutionary state / baryon conversion efficiency as the underlying driver.

### What This Does NOT Claim

- This is **not** a universal law for all disk galaxies
- It **does not** transfer to low-Vflat dwarfs (where all tested predictors fail equally)
- It has **not** been fully replicated on an independent external survey
- The boundary at ~70 km/s is **gradual**, not sharp

## Repository Structure

```
artifacts/galaxy-analyzer/
├── public/replication/          # Key deliverables
│   ├── REPRODUCIBLE_RESULT.md   # Current result summary and replication guide
│   ├── paper.md                 # Historical phases 1-100 analysis (superseded)
│   ├── N45_final_dataset.csv    # Quality-controlled N=45 dataset (phases 1-100)
│   └── replicate_from_scratch.py # Python replication script (phases 1-100)
├── scripts/                     # All analysis phase scripts (Node.js)
│   ├── phase112-matched-falsification.cjs
│   ├── phase110-uncertainty-propagation.cjs
│   ├── phase111-boundary-mapping.cjs
│   ├── phase108-mediation-dissection.cjs
│   ├── phase106-target-robustness.cjs
│   ├── phase105-death-match-sparse.cjs
│   ├── phase104-external-replication.cjs
│   ├── phase103-robustness-battery.cjs
│   ├── phase102-residual-physics.cjs
│   ├── phase101-null-geometric-coupling.cjs
│   └── ... (phases 1-100, earlier investigation)
├── public/                      # JSON results for each phase
│   ├── phase112-matched-falsification.json
│   ├── phase110-uncertainty-propagation.json
│   └── ... (all phase results)
├── src/                         # React + Vite web application
└── experiment-log.txt           # Chronological experiment log
```

## Quick Reproduction

```bash
cd artifacts/galaxy-analyzer/
# Download SPARC data first:
# curl -o /tmp/sparc_table1.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table1.dat
# curl -o /tmp/sparc_cds.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat

# Core finding
node scripts/phase102-residual-physics.cjs

# Robustness + external validation
node scripts/phase103-robustness-battery.cjs
node scripts/phase104-external-replication.cjs

# Competition and mediation
node scripts/phase105-death-match-sparse.cjs
node scripts/phase106-target-robustness.cjs
node scripts/phase108-mediation-dissection.cjs

# Uncertainty and boundary
node scripts/phase110-uncertainty-propagation.cjs
node scripts/phase111-boundary-mapping.cjs

# Decisive test
node scripts/phase112-matched-falsification.cjs
```

## Data

All analysis uses the publicly available [SPARC database](http://astroweb.cwru.edu/SPARC/) (175 galaxies, 3600+ rotation curve points).

**Reference:** Lelli, F., McGaugh, S.S. & Schombert, J.M., 2016, AJ, 152, 157

## Technical Stack

- **Analysis**: Node.js scripts (CommonJS)
- **Web app**: React + TypeScript + Vite + Tailwind CSS + shadcn/ui
- **Statistics**: OLS regression, LOO cross-validation, permutation tests, bootstrap, partial correlations, nearest-neighbor matching, within-bin shuffle — all implemented from scratch

## Zenodo

- v1: DOI 10.5281/zenodo.19430634 (phases 96-100)
- v2: DOI 10.5281/zenodo.19431186 (phases 1-104)
- v3: Pending (phases 1-112, current)

## License

This is a research project using publicly available data. The analysis code and results are provided for scientific reproducibility.

# Galaxy Rotation Curve Analyzer

A data-driven investigation into the physical drivers of outer support requirement in galaxy rotation curves, using the SPARC database (Lelli, McGaugh & Schombert 2016).

## Main Result

**Gas-to-stellar state is the leading predictor of outer support requirement in high-Vflat SPARC galaxies.**

After ruling out the original curve-shape law as a geometric artifact, the project finds that gas-to-stellar balance, quantified by log(MHI/L3.6), is the strongest independently surviving, non-circular, falsification-resistant predictor. It is not a proxy for gas mass alone, stellar mass alone, surface structure alone, or halo properties alone, but represents an independent relative physical state of the galaxy.

### Status

FROZEN (Phase 122). Analysis complete. No further phases planned unless external replication data becomes available or a falsification test kills the claim.

### The Claim (5 Statements)

1. **The old law fell.** The original Vflat/InnerVmax-to-outerSlope relation is a geometric artifact, reproduced by smooth null curves (Phase 101).

2. **The winner is log(MHI/L3.6).** Gas-to-stellar balance is the strongest surviving independent predictor (r=0.729, LOO=0.512, matched falsification p=0.034).

3. **It is not circular.** Anti-circularity battery passes 6/6 tests. Stars-only target strengthens to r=0.885. Permutation: 7.5 sigma (Phase 117).

4. **It is a true ratio state, not one component.** The ratio beats gas alone (r=-0.376), stars alone (r=-0.705), and both residuals. The signal is the balance itself (Phase 120).

5. **The result is regime-specific.** Gradual onset above ~55-75 km/s. Does not transfer to low-Vflat dwarfs (Phases 104, 111).

### Adopted Model

```
logOMD = 1.749 + 0.203 * log(MHI/L3.6) - 0.101 * log(Mbar)
LOO R^2 = 0.584, both coefficients bootstrap-stable
```

### Key Numbers

| Quantity | 1-var: log(MHI/L36) | 2-var: + log(Mbar) |
|:---|:---|:---|
| Pearson r / R | 0.729 | 0.755 |
| LOO R^2 | 0.512 | 0.584 |
| Permutation p | < 0.0001 | 0.001 (added var) |
| External R^2 (70/30) | 0.474 | - |
| N | 104 | 97 (complete data) |

## Evidence Chain

| Phase | Finding | Verdict |
|:---|:---|:---|
| 1-100 | Ratio law Vflat/InnerVmax -> outerSlope | SUPERSEDED |
| **101** | **Null geometric coupling test** | **3/4 FAIL - artifact** |
| 102 | Catalog predictors of outer mass discrepancy | fgas LOO R^2=0.503 |
| 103 | Robustness battery | 3/5 PASS |
| **104** | **Locked external replication** | **ext R^2=0.474** |
| **105** | **Death match sparse models** | **fgas STRONG BACKBONE** |
| 106 | Target definition robustness | MODERATE PASS |
| **108** | **Mediation / proxy dissection** | **fgas CAUSAL, Sigma0 mediated** |
| **110** | **Uncertainty propagation (MC)** | **STRONG PASS** |
| **111** | **Regime boundary mapping** | **Gradual ~55-75 km/s** |
| **112** | **Matched falsification** | **log(MHI/L36) wins** |
| 113 | RAR independence | 65% independent |
| 115 | BTFR independence | r=0.056, orthogonal |
| **117** | **Anti-circularity battery** | **6/6 PASS** |
| 118 | Halo anchoring | 76% persists beyond halo |
| 119 | Model competition | No single model wins |
| **120** | **Mechanism map** | **Ratio-state PASS** |
| **121** | **Minimal physical model** | **2-var adopted** |
| **122** | **Writing freeze** | **FROZEN** |

## What This Does NOT Claim

- Not a universal law for all disk galaxies
- Does not transfer to low-Vflat dwarfs
- Not fully replicated on an independent external survey
- No single causal mechanism identified
- Not peer-reviewed

## Repository Structure

```
artifacts/galaxy-analyzer/
├── public/replication/
│   ├── REPRODUCIBLE_RESULT.md   # Full result + replication guide (FROZEN)
│   ├── MANUSCRIPT.md            # Manuscript draft
│   ├── LITERATURE_SUMMARY.md    # Literature review
│   ├── references.md            # 22 papers, 5 tiers
│   ├── paper.md                 # Historical phases 1-100 (superseded)
│   └── N45_final_dataset.csv    # Historical dataset
├── scripts/                     # All analysis scripts (Node.js)
│   ├── phase120-mechanism-map.cjs
│   ├── phase121-minimal-physical-model.cjs
│   ├── phase117-anti-circularity.cjs
│   ├── phase118-halo-anchoring.cjs
│   ├── phase119-model-competition.cjs
│   ├── phase112-matched-falsification.cjs
│   └── ... (all phases)
├── public/                      # JSON results for each phase
└── src/                         # React + Vite web application
```

## Quick Reproduction

```bash
# Download SPARC data:
curl -o /tmp/sparc_table1.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table1.dat
curl -o /tmp/sparc_cds.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat

cd artifacts/galaxy-analyzer/

# Core result
node scripts/phase102-residual-physics.cjs

# Anti-circularity
node scripts/phase117-anti-circularity.cjs

# Mechanism map
node scripts/phase120-mechanism-map.cjs

# Minimal model
node scripts/phase121-minimal-physical-model.cjs
```

## Data

All analysis uses the publicly available [SPARC database](http://astroweb.cwru.edu/SPARC/) (175 galaxies, 3600+ rotation curve points).

**Reference:** Lelli, F., McGaugh, S.S. & Schombert, J.M., 2016, AJ, 152, 157

## Zenodo

- v1: DOI 10.5281/zenodo.19430634 (phases 96-100)
- v2: DOI 10.5281/zenodo.19431186 (phases 1-104)
- v3: DOI 10.5281/zenodo.19431363 (phases 1-116 + literature)
- v4: DOI 10.5281/zenodo.19431868 (phases 1-119 + anti-circularity)
- v5: DOI 10.5281/zenodo.19432087 (phases 1-122, FROZEN)
- **v6: DOI 10.5281/zenodo.19432600 (corrected post-morphT fix, 91 files)**

## License

This is a non-peer-reviewed research project using publicly available data. The analysis code and results are provided for scientific reproducibility.

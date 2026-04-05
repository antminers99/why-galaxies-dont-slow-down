# Galaxy Rotation Curve Analyzer

A data-driven investigation into whether the MOND acceleration scale $a_0$ varies systematically between galaxies, using the SPARC database (Lelli, McGaugh & Schombert 2016).

## Main Result

**A single-variable physical law governs rotation curve shape:**

$$\text{outerSlope} \approx 1.92 \times \log_{10}\!\left(\frac{V_\text{flat}}{V_\text{max}^\text{inner}}\right) + 0.018$$

The ratio of asymptotic velocity to inner peak velocity predicts the outer rotation curve slope with $r = 0.85$, $p < 0.0001$, across 104 massive spiral galaxies. This survives all measurement controls (8/8 physical, 0/8 measurement artifact).

**What it means:** Galaxies that build their rotation velocity quickly in the center have flat or declining outer curves. Those still rising internally have rising outer curves. This encodes halo concentration — the inside-to-outside mass distribution.

**What happened to $a_0$:** The per-galaxy MOND acceleration scale is a noisy proxy for this ratio. Adding $a_0$ to the ratio model gives $\Delta R^2 = -0.009$ (negative). $a_0$ is not needed.

## Repository Structure

```
artifacts/galaxy-analyzer/
├── public/replication/          # Key deliverables
│   ├── REPRODUCIBLE_RESULT.md   # Concise replication guide (start here)
│   ├── paper.md                 # Full 100-phase analysis paper
│   ├── N45_final_dataset.csv    # Quality-controlled N=45 dataset
│   ├── replicate_from_scratch.py # Self-contained Python replication script
│   ├── variable_definitions.md  # Variable definitions
│   ├── model_spec.md            # Model equations
│   └── expected_outputs.md      # Expected numerical results
├── scripts/                     # All 100 analysis phase scripts (Node.js)
│   ├── phase96-direct-law.cjs         # Discovery of the ratio law
│   ├── phase97-geometric-test.cjs     # Geometric artifact tests
│   ├── phase98-external-replication.cjs # Cross-sample replication
│   ├── phase99-calibration-law.cjs    # Distance calibration
│   ├── phase100-physical-vs-measurement.cjs # Final physical vs artifact test
│   └── ... (phases 1-95)              # Earlier investigation phases
├── public/                      # JSON results for each phase
│   ├── phase100-physical-vs-measurement.json
│   ├── phase99-calibration-law.json
│   └── ... (all phase results)
├── src/                         # React + Vite web application
│   ├── pages/                   # Interactive analysis pages
│   └── components/              # UI components
└── experiment-log.txt           # Chronological experiment log
```

## Quick Reproduction

### Option 1: From SPARC raw data (recommended)

See [`REPRODUCIBLE_RESULT.md`](artifacts/galaxy-analyzer/public/replication/REPRODUCIBLE_RESULT.md) for step-by-step instructions.

Download SPARC data from CDS/VizieR:
- [table1.dat](https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table1.dat) (galaxy properties)
- [table2.dat](https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat) (rotation curves)

### Option 2: Python script

```bash
cd artifacts/galaxy-analyzer/public/replication/
python replicate_from_scratch.py
```

### Option 3: Node.js scripts

```bash
# Download SPARC data first, then:
cd artifacts/galaxy-analyzer/
node scripts/phase100-physical-vs-measurement.cjs
```

## Key Evidence Chain

| Phase | Finding | Result |
|:---|:---|:---|
| 96 | $a_0$ superseded by $V_\text{flat}/V_\text{max}^\text{inner}$ | LOO $R^2 = 0.69$ vs $a_0$'s $-0.02$ |
| 97 | Law is genuine, not geometric artifact | 4/4 genuine, 0/4 geometric |
| 98 | Partially replicates on N=45 subset | $r = 0.63$, $p < 0.0001$ |
| 99 | Calibration depends on distance | LOO $R^2 = 0.75$ with correction |
| 100 | **Physical, not measurement artifact** | **8/8 physical, 0/8 measurement** |

## Phase 100 Controls

The law survives all of these:

- **Every distance bin**: Near $r = 0.93$, Mid $r = 0.78$, Far $r = 0.83$ (all $p < 0.0001$)
- **51 distance-matched pairs**: $r = 0.81$, $p < 0.0001$
- **Partial correlation** after D, nPts, radial coverage: $r = 0.81$ (only 5% drop)
- **Alternative proxies**: mean inner V ($r = 0.72$), baryon-only $V_\text{max}$ ($r = 0.53$)

## Data

All analysis uses the publicly available [SPARC database](http://astroweb.cwru.edu/SPARC/) (175 galaxies, 3600+ rotation curve points).

**Reference:** Lelli, F., McGaugh, S.S. & Schombert, J.M., 2016, AJ, 152, 157

## Technical Stack

- **Analysis**: Node.js scripts (CommonJS)
- **Web app**: React + TypeScript + Vite + Tailwind CSS + shadcn/ui
- **Statistics**: OLS regression, LOO cross-validation, permutation tests, Theil-Sen robust regression, partial correlations — all implemented from scratch (no external stats libraries)

## License

This is a research project using publicly available data. The analysis code and results are provided for scientific reproducibility.

# Galaxy Rotation Curve Analyzer

Interactive web application and analysis pipeline for investigating physical drivers of outer mass discrepancy in SPARC galaxies.

## Main Finding

**Gas-to-stellar state log(MHI/L3.6) is the strongest independent predictor of outer mass discrepancy in high-Vflat SPARC galaxies (N=104, Vflat >= 70 km/s).**

The original ratio law (Vflat/InnerVmax -> outerSlope) was shown to be a geometric artifact (Phase 101). The current result uses logOMD = log10(mean Vobs^2/Vbar^2 in outer disk) as target.

See [`public/replication/REPRODUCIBLE_RESULT.md`](public/replication/REPRODUCIBLE_RESULT.md) for the full reproducible result.

## Adopted Model

```
logOMD = 1.749 + 0.203 * log(MHI/L3.6) - 0.101 * log(Mbar)
LOO R^2 = 0.584, both coefficients bootstrap-stable
```

## Status

FROZEN (Phase 122). No further analysis phases unless external replication data or falsification test.

## Structure

- `scripts/` — Analysis phase scripts (Node.js, run independently)
- `public/` — JSON results from each phase
- `public/replication/` — Reproducible result, manuscript, datasets, Python replication script
- `src/` — React + Vite web application for interactive exploration

## Running the Analysis

```bash
curl -o /tmp/sparc_table1.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table1.dat
curl -o /tmp/sparc_cds.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat

node scripts/phase102-residual-physics.cjs
node scripts/phase117-anti-circularity.cjs
node scripts/phase120-mechanism-map.cjs
node scripts/phase121-minimal-physical-model.cjs
```

## Running the Web App

```bash
pnpm install
PORT=3000 BASE_PATH=/ pnpm run dev
```

## Data

Uses the publicly available [SPARC database](http://astroweb.cwru.edu/SPARC/) (Lelli, McGaugh & Schombert 2016, AJ 152, 157).

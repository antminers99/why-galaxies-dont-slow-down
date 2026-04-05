# Galaxy Rotation Curve Analyzer

Interactive web application and analysis pipeline for investigating the relationship between inner rotation curve structure and outer slope in SPARC galaxies.

## Main Finding

The ratio `Vflat / InnerVmax` predicts outer rotation curve slope with r = 0.85 (p < 0.0001) across 104 massive spirals. This single observable encodes halo concentration and supersedes the per-galaxy MOND acceleration scale a0.

See [`public/replication/REPRODUCIBLE_RESULT.md`](public/replication/REPRODUCIBLE_RESULT.md) for the full reproducible result.

## Structure

- `scripts/` — 100 analysis phase scripts (Node.js, run independently)
- `public/` — JSON results from each phase
- `public/replication/` — Reproducible result, full paper, datasets, Python replication script
- `src/` — React + Vite web application for interactive exploration

## Running the Analysis

```bash
# Download SPARC data to /tmp/
curl -o /tmp/sparc_table1.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table1.dat
curl -o /tmp/sparc_cds.dat https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/table2.dat

# Run any phase script
node scripts/phase100-physical-vs-measurement.cjs
```

## Running the Web App

```bash
npm install
PORT=3000 BASE_PATH=/ npm run dev
```

## Data

Uses the publicly available [SPARC database](http://astroweb.cwru.edu/SPARC/) (Lelli, McGaugh & Schombert 2016, AJ 152, 157).

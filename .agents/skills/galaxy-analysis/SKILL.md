---
name: galaxy-analysis
description: How to run analysis scripts, modify reports, and work with the Galaxy Rotation Curve Analyzer codebase. Use when writing new analysis scripts, editing reports, or modifying the app.
---

# Galaxy Rotation Curve Analysis — Complete Experiment Reference

## Quick Access

- **Full experiment log (all equations/calculations):** `artifacts/galaxy-analyzer/experiment-log.txt`
- **Diagnostic report (26 sections):** `artifacts/galaxy-analyzer/diagnostic-report.txt`
- **Canonical summary page:** `/canonical` route
- **Archived full report:** `artifacts/galaxy-analyzer/galaxy-rotation-curve-analysis-full-report.txt` (reference only)

## Project Structure

```
artifacts/galaxy-analyzer/
├── experiment-log.txt                 ← ALL equations + calculations (36 sections, 15 parts)
├── diagnostic-report.txt              ← PRIMARY report (26 sections)
├── galaxy-rotation-curve-analysis-full-report.txt  ← ARCHIVED (32 parts)
├── scripts/
│   ├── definitive-v4.cjs             ← v4.0 DEFINITIVE PIPELINE (current)
│   ├── phase10-second-parameter.cjs  ← Phase 10 second-parameter search
│   ├── phase9-benchmark-replication.cjs
│   ├── phase8-distance-metareg.cjs
│   ├── phase7-anchor-refit.cjs
│   ├── phase6-matched-splits.cjs
│   ├── phase5-kinematic-audit.cjs
│   ├── diagnostic-tests.cjs          ← Phase 1 diagnostics
│   ├── diagnostic-phase2.cjs         ← Phase 2 corrected residuals
│   ├── final-measurement-v3.cjs      ← v3.0 pipeline (archived)
│   ├── transition-scale-v3.cjs       ← Core RAR analysis (SPARC + LT data)
│   ├── measurement-pipeline.cjs      ← BTFR MC simulation
│   ├── break-the-idea.cjs            ← 5 stress tests
│   ├── redshift-feasibility.cjs      ← z-test forecast
│   ├── theoretical-framework.cjs     ← Theoretical audit
│   └── [other analysis scripts]
├── public/
│   ├── definitive-v4-results.json    ← v4.0 results (CURRENT)
│   ├── phase10-second-param-results.json
│   ├── phase9-benchmark-results.json
│   ├── phase8-metareg-results.json
│   ├── phase7-anchor-results.json
│   ├── phase6-matched-results.json
│   ├── phase5-kinematic-results.json
│   ├── diagnostic-results.json
│   ├── diagnostic-phase2-results.json
│   ├── figure1-data.json (v1.1.0)
│   ├── sparc-table.json (175 galaxies)
│   └── [other result files]
└── src/
    ├── App.tsx
    ├── pages/
    │   ├── canonical-summary.tsx      ← Canonical conclusions page
    │   ├── why-a0.tsx
    │   ├── pipeline.tsx
    │   └── break-test.tsx
    └── components/
        └── layout.tsx                 ← GlassCard layout
```

## Equations Reference

### RAR (McGaugh) — sole estimator

```
g_obs = g_bar / (1 - exp(-sqrt(g_bar / a0)))
```

### Baryonic acceleration

```
g_bar = (Y_disk * V_disk^2 + Y_bulge * V_bulge^2 + V_gas^2) / r
```

Y_disk: 0.5 fixed (v3.0), marginalized [0.20, 0.80] with Gaussian prior N(0.50, 0.12^2) (v4.0)
Y_bulge: 0.7 fixed

### DerSimonian-Laird hierarchical model

```
Level 1: r_ij = log10(g_obs,ij) - log10(RAR(g_bar,ij, a0_i)) + e_ij
Level 2: log(a0_i) = mu + u_i,  u_i ~ N(0, tau^2)
Level 3: mu = population mean log(a0)

Q = sum w_i * (theta_i - theta_FE)^2
tau^2 = max(0, (Q-(k-1)) / (S1 - S2/S1))
W_i = 1/(se_i^2 + tau^2)
mu = sum(W_i * theta_i) / sum(W_i)
SE = 1 / sqrt(sum(W_i))
I^2 = max(0, (Q-(k-1))/Q) * 100%
```

### ANOVA decomposition

```
total^2 = within^2 + between^2
0.0239  = 0.0148   + 0.0091   (variances in dex^2)
0.155^2 = 0.122^2  + 0.096^2  (RMS in dex)
```

### Phase 10 second-parameter model

```
logA0_i = mu + beta * (X_i - mean(X)) / sd(X)
dAIC = n*ln(RSS1/n) + 4 - n*ln(RSS0/n) - 2
```

### Unit conversion

```
1 (km/s)^2/kpc = 3.241e-14 m/s^2
c*H0/(2*pi) = 1.042e-10 m/s^2
```

## Headline Numbers (v4.0 GOLD+i45)

| Quantity | Value |
|----------|-------|
| a0 | 3633 (km/s)^2/kpc = 1.18e-10 m/s^2 |
| log(a0) | 3.560 +/- 0.041 (stat) +/- 0.304 (total) dex |
| tau | 0.291 dex |
| I^2 | 92.4% |
| a0/(cH0/2pi) | 1.130 |
| 68% CI | [1846, 7149] (km/s)^2/kpc |
| Within-galaxy scatter | 0.122 dex (matches McGaugh 0.13) |
| Best 2nd param | eta_rot: r=0.334, p=0.010, dAIC=-4.6, -5.7% tau |

## NUMBER HIERARCHY (v4.0 locked)

- **HEADLINE** = v4.0 marginalized hierarchical **3633** — always use this
- **BASELINE** = v3.0 hierarchical **3374** — for comparison only
- **DIAGNOSTIC/RAW** = v3.0 unweighted **4837** — never quote as "the" result
- **Literature** = McGaugh 2016 = **3703** (1.20e-10 m/s^2)

## Sample Definitions

| Sample | Criteria | n_gal (v4.0) |
|--------|----------|------------|
| ALL | no cuts | 165 |
| CLEAN | Vmax>=80, inc>=30, n>=5, gRange>=0.5 | 84 |
| GOLD | Vmax>=50, inc>=30, n>=5, gRange>=1.0 | 64 |
| **GOLD+i45** | Vmax>=50, inc>=45, n>=5, gRange>=1.0 | **59 (PRIMARY)** |
| ANCHOR | GOLD+i45 + precise D + Q=1 + inc>=60 | 12 |

## Running Analysis Scripts

All scripts are CommonJS (.cjs) and run with Node:
```bash
cd artifacts/galaxy-analyzer && node scripts/definitive-v4.cjs
cd artifacts/galaxy-analyzer && node scripts/phase10-second-parameter.cjs
```

Scripts read data from:
- `/tmp/rotmod/` — SPARC rotation curve files (175 galaxies)
- `/tmp/little_things/` — LITTLE THINGS data files
- `public/*.json` — Previously computed results

If `/tmp/rotmod/` is empty, the SPARC data needs to be re-downloaded.

## Data Format

### Plot Points (transition-scale.json -> plotPoints)
```
x = log10(g_bar) in (km/s)^2/kpc
y = log10(g_obs/g_bar)
g = galaxy name
```

To recover physical values:
- g_bar = 10^x (km/s)^2/kpc
- g_obs = g_bar * 10^y = 10^(x+y) (km/s)^2/kpc

## JSX/Template Rules

- NEVER use backticks inside JSX (template literals break the parser)
- GlassCard `glow` prop: only `'cyan'` | `'purple'` | `'amber'` | `'none'`
- Use string concatenation instead of template literals in JSX

## Diagnostic Report Structure

The `diagnostic-report.txt` has 26 sections:
- Sections 1-11: Empirical diagnostics and methodology
- Sections 12-14: v3.0 baseline (archived)
- Sections 15-19: v4.0 definitive pipeline
- Section 20: Systematic suspects audit
- Section 21: Phase 5 kinematic audit
- Section 22: Phase 6 matched splits
- Section 23: Phase 7 anchor refit
- Section 24: Phase 8 distance meta-regression
- Section 25: Phase 9 benchmark replication
- Section 26: Phase 10 second-parameter search

## Experiment Log Structure

The `experiment-log.txt` has 36 sections across 15 parts:
- Part I (1-4): Physical setup and equations
- Part II (5-7): Fitting methodology
- Part III (8-9): Hierarchical model
- Part IV (10): Sample definitions
- Part V (11-15): Complete numerical results
- Part VI (16-20): Diagnostic tests
- Part VII (21-22): Phase 5 kinematic audit
- Part VIII (23-24): Phase 6 matched splits + power analysis
- Part IX (25): Phase 7 anchor refit
- Part X (26-27): Phase 8 distance meta-regression
- Part XI (28-30): Phase 9 benchmark replication
- Part XII (31-33): Phase 10 second-parameter search
- Part XIII (34): Systematic suspects catalog
- Part XIV (35): All equations summary
- Part XV (36): Paper-ready statement

## CRITICAL RULES

- HEADLINE number = v4.0 marg. hier. 3633 (NEVER use raw UW 4837)
- BASELINE for comparison = v3.0 hier. 3374
- NEVER mix estimators (McGaugh RAR ONLY)
- NEVER claim: universality proven, MOND confirmed, dark matter solved
- Four scatter measures are DIFFERENT quantities — never conflate them
- tau DOMINATES the uncertainty budget (91.6% of total variance)
- 94% of heterogeneity remains UNEXPLAINED after Phase 10

## Resolved Questions

- Gas fraction / SB residual correlations: ARTIFACTS (Phase 1 code bug, corrected)
- Low-mass a0 divergence: FITTING ARTIFACT (73% have <0.5 dex g_bar range)
- Quality-dependent a0: NOT an issue (high/low quality give same a0)
- Distance method bias: NOT significant (meta-regression p=0.55)
- Y* as scatter driver: NOT the driver (marginalization changes tau by only 4%)

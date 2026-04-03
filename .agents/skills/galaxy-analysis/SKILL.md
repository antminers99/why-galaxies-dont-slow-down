---
name: galaxy-analysis
description: How to run analysis scripts, modify reports, and work with the Galaxy Rotation Curve Analyzer codebase. Use when writing new analysis scripts, editing reports, or modifying the app.
---

# Galaxy Analysis Workflow

## Project Structure

```
artifacts/galaxy-analyzer/
├── diagnostic-report.txt              ← PRIMARY report (use for all future edits)
├── galaxy-rotation-curve-analysis-full-report.txt  ← ARCHIVED (32 parts, ~3745 lines)
├── scripts/
│   ├── diagnostic-tests.cjs           ← Main diagnostic pipeline (3 tests)
│   ├── transition-scale-v3.cjs        ← Core RAR analysis (SPARC + LT data)
│   ├── measurement-pipeline.cjs       ← BTFR MC simulation
│   ├── break-the-idea.cjs             ← 5 stress tests
│   ├── redshift-feasibility.cjs       ← z-test forecast
│   ├── theoretical-framework.cjs      ← Theoretical audit
│   └── [20+ other analysis scripts]
├── public/
│   ├── diagnostic-results.json        ← Latest test results
│   ├── transition-scale.json          ← Full RAR analysis output
│   ├── rar-analysis-real.json         ← Per-galaxy metadata (175 SPARC)
│   ├── break-test-results.json        ← Stress test results
│   └── [15+ other result files]
└── src/
    ├── App.tsx                         ← Main app with routes
    ├── pages/
    │   ├── why-a0.tsx                  ← Theoretical investigation page
    │   ├── pipeline.tsx                ← Measurement pipeline page
    │   └── break-test.tsx              ← Breaking tests page
    └── components/
        └── layout.tsx                  ← GlassCard layout
```

## Running Analysis Scripts

All scripts are CommonJS (.cjs) and run with Node:
```bash
cd artifacts/galaxy-analyzer && node scripts/diagnostic-tests.cjs
```

Scripts read data from:
- `/tmp/rotmod/` — SPARC rotation curve files (175 galaxies)
- `/tmp/little_things/` — LITTLE THINGS data files
- `public/*.json` — Previously computed results

If `/tmp/rotmod/` is empty, the SPARC data needs to be re-downloaded.
Check `scripts/transition-scale-v3.cjs` lines 47-68 for the data loading logic.

## Data Format

### Plot Points (transition-scale.json → plotPoints)
```
x = log10(g_bar) in (km/s)²/kpc
y = log10(g_obs/g_bar)
g = galaxy name (12 chars max)
```
To recover physical values:
- g_bar = 10^x (km/s)²/kpc
- g_obs = g_bar × 10^y = 10^(x+y) (km/s)²/kpc

### g_bar Computation (SPARC)
```
Vbar² = Υ_disk × Vdisk × |Vdisk| + Υ_bulge × Vbul × |Vbul| + Vgas × |Vgas|
g_bar = |Vbar²| / r
```
This is a proper disk mass model, NOT spherical GM/r².

### Interpolation Functions
```javascript
// McGaugh RAR (BEST FIT):
g_obs = g_bar / (1 - Math.exp(-Math.sqrt(g_bar / a0)))

// Simple MOND:
g_obs = g_bar * (0.5 + Math.sqrt(0.25 + a0 / g_bar))

// Standard MOND:
g_obs = g_bar * (1 + Math.sqrt(1 + 4 * a0 / g_bar)) / 2

// Our framework (WORST FIT — avoid for new work):
g_obs = Math.sqrt(g_bar² + a₀² × g_bar / (g_bar + a₀))
```

## Editing Reports

### Active Report: `diagnostic-report.txt`
- 8 sections, ~350 lines
- Focus: what is real vs model-dependent vs open
- Contains actual computed test results

### Archived Report: `galaxy-rotation-curve-analysis-full-report.txt`
- 32 parts, ~3745 lines
- Contains theoretical investigation, derivation audits, novelty audit
- Do NOT modify — reference only

## JSX/Template Rules

When editing React components:
- NEVER use backticks inside JSX (template literals break the parser)
- GlassCard `glow` prop: only `'cyan'` | `'purple'` | `'amber'` | `'none'`
- Use string concatenation instead of template literals in JSX

## Key Open Questions (priority order)

1. What is the correct a₀ value? (range: 2670 global – 3720 high-mass)
2. Is the weak SB trend (r=-0.35, n=46) in highest-quality galaxies real?
3. Is the cosmological ratio a₀ ≈ cH₀/2π real or method-dependent?
4. Does a₀ evolve with redshift?

## Resolved Questions

- Gas fraction / SB residual correlations: ARTIFACTS (Phase 1 code bug, corrected)
- Low-mass a₀ divergence: FITTING ARTIFACT (73% have <0.5 dex g_bar range)
- Quality-dependent a₀: NOT an issue (high/low quality give same a₀)

## Next Tests Needed

1. Settle a₀ value across methods (global, per-galaxy, high-mass, range>1dex)
2. Investigate weak SB trend in 46 well-constrained galaxies
3. Joint fit of a₀ + f_bar to resolve degeneracy
4. Report a₀ from all methods for method-independent comparison

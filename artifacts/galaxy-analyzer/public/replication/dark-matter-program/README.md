# Program 11: From Observational Carrier to Dark Matter Identity

## Project Status: INITIATED
## Parent Project: Galaxy Rotation Curve Analyzer (Programs 1–10)
## Zenodo Parent: DOI 10.5281/zenodo.19430633 (v14)

---

## Motivation

Programs 1–10 established that:

1. A hidden variable **H** governs per-galaxy deviations from the mean Radial Acceleration Relation (RAR)
2. The outer-halo **m=2 azimuthal mode** of HI velocity fields is the strongest observational carrier of H identified to date (r = 0.847, p = 0.005, N = 7 THINGS galaxies)
3. H survives 11/11 red-team tests including partial-correlation control for inclination, Vflat, and Mbar
4. The signal is consistent with triaxial or oval distortions in the dark matter halo potential

**This closes the detection program.** The next question is fundamentally different:

> **If H traces halo shape, what does halo shape tell us about the nature of dark matter?**

---

## Research Question

**What constraints does the H–m2 coupling place on dark matter particle physics and halo formation models?**

Sub-questions:
1. Which dark matter models predict triaxial halos with the observed m=2 power distribution?
2. Can CDM N-body simulations reproduce the observed r(DQ, m2) = 0.85 correlation?
3. Do alternative dark matter models (SIDM, fuzzy DM, MOND + EFE) predict different m=2 signatures?
4. What observational tests would distinguish between these models?

---

## Program Structure

### Phase 11.1: Theoretical Landscape Survey
- Map which DM models predict triaxial halos
- CDM: expected from hierarchical merging, tidal fields
- SIDM: cores may reduce triaxiality — testable prediction
- Fuzzy DM (ψDM): soliton cores + interference patterns → different m spectrum
- MOND + External Field Effect: could mimic m=2 via EFE direction

### Phase 11.2: CDM Halo Shape Predictions
- Literature survey: axis ratios from Illustris, EAGLE, FIRE simulations
- Compare predicted m=2 power distribution with observed THINGS values
- Key question: does CDM predict the RANGE of m=2 we observe?

### Phase 11.3: SIDM Discriminant
- SIDM with σ/m ~ 1 cm²/g → rounder cores → LOWER m=2 in inner regions
- Testable: inner vs outer m=2 gradient as SIDM diagnostic
- Cross-check with existing THINGS data (inner/outer m=2 already computed in Program 9V)

### Phase 11.4: Fuzzy DM Signature
- ψDM predicts specific interference patterns in density/velocity fields
- m spectrum should differ: power at higher m modes (m=3, m=4) relative to m=2
- Testable with existing THINGS data: check m=1, m=2, m=3, m=4 power ratios

### Phase 11.5: Observational Decision Tree
- Design a flowchart: given m=2 power + m-spectrum + radial gradient → which DM model?
- Identify which new observations (WALLABY, MeerKAT, SKA) would be decisive

### Phase 11.6: Quantitative Model Comparison
- For each model, predict r(DQ, m2) and compare with observed 0.847
- Bayesian model comparison if predictions are quantitative enough

---

## Connection to Parent Project

| Parent (Programs 1–10) | This Program (11+) |
|---|---|
| What is H? | What does H tell us about DM? |
| Empirical detection | Theoretical interpretation |
| r(DQ,m2) = 0.847 is the result | r(DQ,m2) = 0.847 is the input |
| SPARC + THINGS data | Same data + simulations |
| MNRAS paper (submitted) | New paper target: ApJ or MNRAS Letters |

---

## Key Input Files (from Parent)

- `data/phase56-frozen-baselines.json` — frozen SPARC baselines
- `data/program9-phase903.json` — decisive test results
- `data/program9v-red-team.json` — red team validation
- `data/program8a-2d-state.json` — 2D velocity field analysis
- `data/things-2d/*.FITS` — raw THINGS velocity fields

---

## Author
- Fnd89, Independent Researcher
- Contact: antminers99@gmail.com
- CC: Super.centauri@yandex.com

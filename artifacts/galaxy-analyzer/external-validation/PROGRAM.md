# External Validation of the Coupling Law

## Program Goal

Test whether the hierarchical coupling law discovered in N=45 (Track 2, Phases 123-134)
transfers to an independent sample of SPARC galaxies not used in model development.

## Central Question

**Does VfResid remain the dominant coupling channel outside N=45, especially in high-Vflat galaxies?**

Sub-questions:
1. Does Core alone fail (as expected)?
2. Does Core+VfResid dominate?
3. Does the 5th axis (lh_outerImprove) remain secondary but real?
4. Is the law still strongest in the high-Vflat regime?

## Frozen Internal Reference (from Track 2)

| Model | LOO Gap% (N=45) | Transfer Gap% (N=10 crude) |
|-------|-----------------|---------------------------|
| Core alone | 36.8% | -11.5% (fails) |
| Core+VfResid | 61.1% | 56.9% |
| 5-axis | 65.4% | 66.0% |

- r(VfResid, a₀) = 0.801 on N=10 holdout
- VfResid activates at Vflat ~ 181 km/s
- Regime: Vflat >= 120 km/s

## Data Strategy

### Source
SPARC galaxies NOT in the original stageA sample (N=56).
RAR data points from transition-scale.json plotPoints.
Galaxy properties from sparc-table.json and sparc-results.json.

### Variables to Compute
- **logA0** (target): Per-galaxy RAR fit using McGaugh interpolation
- **logMHI**: log10(MHI / 1e9 M☉) from SPARC table
- **logMbar**: log10((L36×0.5 + MHI×1.33) × 1e9 M☉) — same formula as Track 2
- **VfResid**: Vflat_obs - Vflat_expected(logMbar) via BTFR from N=45 training
- **lh_outerImprove**: From sparcResults log_halo model
- **logMeanRun**: From RAR residual run analysis
- **logMhost**: Proxy from dark_halo_linear k parameter (exploratory)

### Quality Filters
- Minimum 5 RAR data points per galaxy (for reliable a₀ fit)
- Valid Vflat > 0
- Valid MHI > 0 and L36 > 0
- Q <= 2 preferred (Q=3 excluded or flagged)

### Expected Sample Sizes
- Total external with >= 5 RAR points: ~75
- With valid Vflat: ~55
- High-Vflat (>= 120): ~15
- Very high-Vflat (>= 180): ~6

## Analysis Phases

### Phase 200: Data Assembly
Build the external validation dataset with all computed variables.
Output: phase200-external-dataset.json

### Phase 201: Blind Prediction
Apply frozen Track 2 coefficients to external galaxies.
NO refitting. Pure prediction using training coefficients.
Test: Core, Core+VfResid, 5-axis models.

### Phase 202: Hierarchy Test
Does the same hierarchy hold?
- Core fails alone?
- VfResid is the dominant channel?
- 5th axis adds incremental value?

### Phase 203: Regime Verification
Is the law strongest in Vflat >= 120?
Does VfResid activate near Vflat ~ 181?

### Phase 204: Robustness
Bootstrap stability, outlier sensitivity, point-count effects.

## Rules
1. NO refitting on external data — prediction only
2. Coefficients frozen from Phase 134
3. Compare to proper baseline (RMSE of training mean)
4. Flag extrapolating galaxies (outside training logMbar/logVflat range)
5. Report honestly: if it fails, it fails

---

## Results (Phases 200-201)

### Phase 200: Data Assembly
- 59 external galaxies assembled (from 119 candidates)
- 39 rejected for missing Vflat, 20 for too few RAR points
- 16 high-Vflat (>=120), 8 very-high-Vflat (>=180)
- logA0 from RAR fit, logMhost estimated from Vflat formula

### Phase 201: Blind Prediction — MODERATE_TRANSFER

| Regime | N | Core Gap | Core+VfResid Gap | r(VfResid,a₀) |
|--------|---|----------|-------------------|----------------|
| Full sample | 59 | -53.0% | +8.2% | 0.713 |
| Vflat >= 120 | 16 | -6.7% | +34.3% | 0.841 |
| Vflat >= 180 | 8 | -49.4% | +48.7% | 0.858 |
| Q=1 + Vflat>=120 | 11 | — | +59.0% | 0.830 |
| Vflat < 120 | 43 | -74.4% | -3.9% | — |

**Key findings:**
1. Hierarchy preserved: Core fails, VfResid dominates
2. Regime-dependence confirmed: strong in high-Vflat, fails in low-Vflat
3. In best subset (Q=1, Vflat>=120, N=11): gap=59%, r=0.830 — matches Phase 134
4. Full sample diluted by low-Vflat regime and crude logA0 estimates
5. VfResid coefficient sign preserved; logMhost and logMR signs changed (expected: crude estimates)

### Phase 202: External Hierarchy Replication — STRONG_HIERARCHY_REPLICATION

**8/8 hierarchy checks passed** across all regimes:

| Check | Result |
|-------|--------|
| Core fails (full sample) | YES (gap=-53%) |
| Core fails (high-Vflat) | YES (gap=-6.7%) |
| VfResid dominates (full) | YES |
| VfResid dominates (high-V) | YES |
| Bootstrap P(VfResid>Core) full | 100% |
| Bootstrap P(VfResid>Core) high-V | 100% |
| Partial r(VfResid,a0|Core) > 0.3 | YES (0.844) |
| Regime preserved (high-V > full) | YES |

**Channel dominance:**
- Full sample: VfResid gap=47.5% (alone), next best=logVflat at -1.6% — margin +49pp
- High-Vflat: VfResid gap=34.3%, next best=logVflat at -1.0% — margin +35pp
- Very-high-Vflat: VfResid gap=48.7%, next best=logVflat at -21.8% — margin +71pp

**Bootstrap stability (10000 resamples):**
- Full: P(VfResid>Core)=100%, mean delta=61.1pp, 90% CI=[43.4, 80.0]pp
- High-V: P(VfResid>Core)=100%, mean delta=45.4pp, 90% CI=[18.8, 85.9]pp
- Very-High-V: P(VfResid>Core)=99.9%, mean delta=101.4pp
- Q=1+HV: P(VfResid>Core)=99.6%, mean delta=61.4pp

**lhOuter (5th axis) external status:**
- Full: adds +1.4pp (gap 8.2% → 9.7%) — modest but positive
- High-V: adds +1.0pp (gap 34.3% → 35.3%) — small
- Very-High-V: adds +18.1pp (gap 48.7% → 66.9%) — substantial
- Q=1+HV: adds +12.1pp (gap 59.0% → 71.1%) — strong
- Bootstrap P(lhOuter adds): 71.5% (full), 90.4% (very-high-V)

**Partial correlations (VfResid after Core removed):**
- Full: r=0.806 | High-V: r=0.844 | Very-High-V: r=0.760 | Q=1+HV: r=0.869

**Key discovery:** VfResid-only (no Core at all) outperforms every model in every regime except Q=1+HV where Core+VfResid is better. This means VfResid alone contains more transferable information than the 3-axis core — consistent with the internal finding that VfResid mediates Core content.

### Phase 203: logMhost Improvement — CORE_IMPROVES_VfResid_STILL_DOMINATES

**Central question answered:** Does better logMhost restore Core's external weight?
**Answer:** Core improves dramatically, but VfResid **still dominates in every regime.**

**Estimator training on N=45 (real logMhost from tidal analysis):**

| Model | Features | RMSE (train) | LOO RMSE | r |
|-------|----------|-------------|----------|---|
| A (current formula) | fixed: 11.5+3×(logV-2.2) | 0.8299 | — | 0.232 |
| B (logVflat trained) | logVflat | 0.5874 | 0.6258 | 0.232 |
| C (4-variable) | logVflat+logMbar+logRdisk+T | 0.5073 | 0.5766 | 0.542 |
| D1 (logVflat+logMbar) | logVflat+logMbar | 0.5482 | 0.5940 | 0.419 |
| D2 (logVflat+T) | logVflat+T | 0.5607 | 0.6133 | 0.371 |

**Core gap improvement with Model C (best):**

| Regime | Core (Model A) | Core (Model C) | Delta | VfResid still dominates? |
|--------|---------------|----------------|-------|------------------------|
| Full (N=59) | -53.0% | -8.1% | +44.9pp | YES (margin 40.3pp) |
| High-V (N=16) | -6.7% | +16.7% | +23.4pp | YES (margin 30.0pp) |
| Very-High-V (N=8) | -49.4% | -16.7% | +32.7pp | YES (margin 76.3pp) |
| Q=1+HV (N=11) | +1.9% | +25.4% | +23.5pp | YES (margin 43.0pp) |

**VfResid-only gap unchanged** at 47.5% (full), 75.4% (very-high-V), 63.8% (Q1+HV)

**Interpretation:** The current formula (r=0.232 vs real) was severely degrading Core. Model C (r=0.542) recovers ~45pp of Core gap in the full sample. Yet even with this improvement, VfResid dominates by 30-76pp in every regime. This confirms: **VfResid is the primary transfer channel, not an artifact of bad Core estimates.**

**Hierarchy scores:** Model A: 8/8 | Model B: 7/8 | Model C: 6/8 | Model D1: 6/8 | Model D2: 7/8
Note: Better logMhost actually *reduces* hierarchy score because Core no longer fully "fails" in high-V — it just becomes less dominant. This is expected and healthy.

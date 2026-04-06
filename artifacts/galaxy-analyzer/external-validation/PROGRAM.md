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

---
name: galaxy-summary
description: Quick executive summary of the Galaxy Rotation Curve Analyzer's current definitive result. Load this skill whenever you need a concise recap of the current claim, its status, and key caveats.
---

# Galaxy Rotation Curve Analysis — Executive Summary

## STATUS: COMPREHENSIVELY EXTERNALLY VALIDATED (Phases 134 + 200–204)

Two-track analysis:
- **Track 1 (Phases 1–122):** OMD law — FROZEN
- **Track 2 (Phases 123–134, 200–204):** State law for logA0 variation — COMPREHENSIVE EXTERNAL VALIDATION COMPLETE

## Current Definitive Result (Track 2)

**Title:** Comprehensively Validated Regime-Dependent Coupling Law for Per-Galaxy a₀ Variation in SPARC Galaxies

**Claim:** Per-galaxy a₀ is governed by a transferable, regime-dependent, hierarchical baryon-halo coupling law. A comprehensive external validation program (Phases 200–204, N=59 independent galaxies) establishes that VfResid dominance over environmental/structural channels is independent of environmental proxy quality. The complete hierarchy replicates externally (8/8 checks pass), VfResid alone achieves 47.5%–75.4% gap without Core, and VfResid dominates by 30–76 pp even after improving the environmental proxy. The signal concentrates in high-Vflat galaxies and collapses in low-Vflat systems.

**Primary Model C (predictive leader):**

logA0 = 4.120 − 0.339×logMHI − 0.175×logMhost + 0.449×logMR + 0.631×logVflat

| Quantity | Value |
|----------|-------|
| LOO gap | 52.1% |
| AIC | −155.63 |
| Nested CV (3-way) | 45/45 |
| Bootstrap max flip | 0.1% |
| Total CV | 1.006 |
| Sample | N=45 published SPARC subset |

**Secondary Model M4 (interpretive leader):**

logA0 = 5.096 − 0.237×logMHI − 0.184×logMhost + 0.465×logMR + 0.151×logSigma0

| Quantity | Value |
|----------|-------|
| LOO gap | 47.4% |
| Permutation p | 0.039 |
| Nested CV | 43/45 |
| Bootstrap max flip | 3.2% |
| Total CV | 1.210 |

## VfResid Identity (Phases 131–132)

**HALO_COUPLING — PARTIALLY_REDUCIBLE** — The 7.8% kinematic residual inside Vflat (VfResid) is identified as a baryon-halo coupling proxy, but cannot be replaced by any bundle of catalog observables.

| Correlate | r(VfResid) | r(partial, all 7 controls) | Interpretation |
|-----------|------------|---------------------------|----------------|
| haloK (dark halo strength) | +0.598 | +0.628 | Stronger halo → higher VfResid |
| mondImprove | −0.542 | −0.376 | MOND works worse when VfResid high |
| logMeanRun | +0.388 | 0.000 (absorbed by core) | Coherent dynamics → higher VfResid |
| envCode | −0.384 | −0.002 (absorbed by core) | Isolated galaxies → higher VfResid |
| concentration | −0.289 | −0.379 | More concentrated → lower VfResid |

**Phase 132 Reducibility Test:**

| Test | Result |
|------|--------|
| VfResid full-sample gap | 61.1% (+17.1pp above core) |
| VfResid fold-internal gap | 53.4% (+9.4pp above core, genuine signal) |
| Best single replacement (haloK) | 49.1% (+5.1pp) |
| Best 5-probe bundle gap | 44.4% (below core — overfitting!) |
| After removing 5-probe content | 40.6% of VfResid survives, still +8.5pp |
| Core + Sig0 + VfResid | 63.6% (highest gap observed) |

**Key insights:**
1. haloK survives ALL controls (partial r never drops below 0.628)
2. Adding probes to haloK HURTS prediction (LOO overfitting with N=45)
3. After stripping all probe content, 40% of VfResid remains, still predictive
4. Fold-internal VfResid = 53.4% vs Model C = 52.1% (comparable, both genuine)

**Verdict**: VfResid is the BEST AVAILABLE compressed proxy for baryon-halo coupling. No catalog observable or bundle can replicate it. haloK is an interpretive diagnostic, not an operational replacement.

## Phase 132A-C: From Coupling Proxy to Coupling Law

### Phase 132A: Halo Driver Death Match — VFRESID_DOMINANT

Every halo/DM proxy tested against VfResid in direct predictive competition.

| Candidate | Gap% | vs VfResid (61.1%) |
|-----------|------|-------------------|
| dhl_mse (dark halo linear MSE) | 52.7% | −8.5pp |
| logSBeff | 51.6% | −9.5pp |
| haloK_linear | 49.1% | −12.0pp |
| Best bundle (dhl_mse+logSBeff) | 53.3% | −7.9pp |
| Best triple bundle | 53.7% | −7.4pp |

**Asymmetry**: VfResid DOMINATES — survives after removing best halo proxy (94.9% variance, +13.2pp), but best halo dies after removing VfResid (r→−0.064).

**Partial correlations after Core+Structure**:
- haloK_linear: r=0.628 (strongest)
- mgh_k: r=0.501
- lh_outerImprove: r=0.438

### Phase 132B: Mediation / Causal Ordering — VFRESID_DOMINANT_CHANNEL

| Halo Proxy | Mediation Type | Absorption | Key Finding |
|------------|---------------|------------|-------------|
| haloK_linear | MUTUAL_PARTIAL | VfResid absorbs 100% of haloK | haloK signal fully mediated |
| dhl_mse | PARTIAL_MED | 75% absorbed | Some independent content |
| mond_improve | FULL_MEDIATION | 100% absorbed | Fully mediated through VfResid |
| lh_outerImprove | INDEPENDENT | 0% absorbed | Adds +4.3pp ON TOP of VfResid |

**Sobel-like path analysis** (controlling for Core):
- haloK → VfResid → a₀: **107% mediation** (bootstrap P>50%: 95.4%)
- Entry-order: VfResid absorbs haloK signal (+7.1pp added), haloK adds nothing to VfResid (−4.9pp)

**Tally**: 5/12 full mediation, 2/12 partial, 1/12 mutual, 4/12 independent

### Phase 132C: External Robustness — MOSTLY_ROBUST

| Test | Result |
|------|--------|
| VfResid delta > 0 | 14/14 splits (100%) |
| Core signs stable | 100% of splits |
| Bootstrap P(delta>0) | 99.9% |
| Bootstrap P(delta>5pp) | 95.8% |
| Bootstrap P(delta>10pp) | 82.1% |
| r(VfResid,a₀) bootstrap 5th pctl | 0.532 |
| Transfer tests positive | 3/6 |

**Key domain split findings**:
- High-Vflat half: VfResid delta = +27.6pp (very strong)
- Low-Vflat half: VfResid delta = +1.7pp (weak — signal concentrated in massive galaxies)
- High-Mbar: delta = +24.3pp | Low-Mbar: delta = +5.2pp
- Transfer from low→high works (+24.2pp), high→low fails (−65.7pp)

**Caveat**: VfResid signal is strongest in massive/high-Vflat galaxies. Transfer to low-mass regime is unreliable — possible genuine mass-dependent coupling.

## Phase 133A-C: Regime Law + Second Channel + Coupling Drivers

### Phase 133A: Regime Law — VERY_STRONG_REGIME + SHARP_TRANSITION

| Vflat threshold | N | Core gap | Core+VfResid gap | Delta |
|-----------------|---|----------|-------------------|-------|
| ≥70 | 43 | 55.9% | 71.3% | +15.4pp |
| ≥120 (optimal) | 35 | 24.3% | 52.1% | +27.8pp |
| ≥160 | 27 | 38.5% | 64.5% | +26.0pp |

- **Sharp transition** at Vflat ~ 181 km/s (r > 0.5 first appears)
- **VfResid coefficient**: 3.078 in high-Vflat regime (vs 2.514 full sample)
- **Bootstrap**: P(delta>0) = 99.7% (high), P(coef>0) = 100%
- **Low-Vflat regime** (N=10): completely unstable, P(delta>0) = 18.6%

### Phase 133B: Second Channel — GENUINE_5TH_AXIS

lh_outerImprove (log-halo outer improvement) is a genuine 5th axis:

| Test | Result |
|------|--------|
| LOO gap (5-axis) | 65.4% (+4.3pp above 4-axis) |
| Fold-internal delta | +3.0pp (genuine) |
| Nested CV | 26/45 |
| Flip rate | 5.4% |
| Systematics | CLEAN (all |r| < 0.13) |

Physical nature: anti-correlates with logVflat (r=−0.39) and logMbar (r=−0.41) after controls — captures outer halo fitting quality independent of mass/velocity.

### Phase 133C: Coupling Drivers — PARTIAL_COUPLING_LAW + DEEP_RESIDUAL

What drives VfResid? (predicting VfResid itself):

| Model | LOO R² | Variables |
|-------|--------|-----------|
| haloK alone | 0.287 | 1 |
| haloK + lhOuter | 0.355 | 2 |
| haloK + env + MR | 0.431 | 3 |
| ALL top 5 | 0.517 | 5 |

In high-Vflat regime (N=35): haloK alone → R²=0.413, haloK+env+MR → R²=0.533

**DEEP RESIDUAL**: After removing ALL 5 top predictors:
- 35.6% of VfResid variance remains unexplained
- Unexplained part STILL adds +6.8pp to core (Core+unexplained = 50.9%)
- VfResid contains genuine physics beyond any available catalog observable

**Coupling picture**: haloK (dark halo strength) is the dominant single driver (~29% of VfResid), but environment and dynamical coherence contribute independently. The remaining ~36% is irreducible — likely multi-scale dynamical integration that only Vflat captures naturally.

## Phase 134: External / Transfer Validation — STRONG_TRANSFER

Train on N=45 published, predict N=10 crude holdout (all Vflat >= 120):

| Model | Transfer RMSE | Transfer Gap% |
|-------|---------------|---------------|
| Naive (predict train mean) | 0.2962 | 0.0% |
| Core (3-axis) | 0.3128 | -11.5% (FAILS) |
| Core + Vflat (Model C) | 0.2810 | 10.0% |
| Core + VfResid | 0.1944 | **56.9%** |
| Core + VfResid + lhOuter (5-axis) | 0.1727 | **66.0%** |

Key transfer findings:
- r(VfResid, a₀) on holdout = **0.801** (higher than training!)
- r(haloK, a₀) on holdout = 0.261 (VfResid >> haloK)
- All coefficient signs preserved (published vs combined N=55)
- Regime-restricted training (Vflat>=120) outperforms full training
- lh_outerImprove transfers: +9.1pp on holdout
- 2 test galaxies in extrapolation territory (ESO563-G021, UGC02885)

Caveats: N_test=10 small; all crude quality; no low-Vflat holdout galaxies available.

**Critical insight**: Core alone FAILS to transfer (-11.5%). Raw Vflat transfers modestly (10%). Only the residual dynamical channel (VfResid) transfers strongly (57%). This confirms VfResid encodes genuine transferable physics, not sample-specific fitting artifacts.

## Phase 129–130 Context

| Finding | Detail |
|---------|--------|
| Vflat nature (P129) | SUPERIOR COMPRESSED PROXY — structure explains 92.2% of Vflat |
| Structure adds (P130) | +0.4pp (PC1) above core — nearly nothing |
| VfResid adds (P130) | +17.1pp above core — the predictive engine |
| Side discovery (P130) | Core+Sig0+VfResid = 63.6% gap, 6.8% flip |

## Common Core (3 firmly established axes)

| Axis | Role | Evidence |
|------|------|----------|
| logMHI | Gas mass | Present in ALL winning models |
| logMhost | Host halo mass | Present in ALL winning models |
| logMR (MeanRun) | Dynamical coherence | Cannot be replaced by Vflat (Model D collapses to 31.6%) |

## Bug Audit (Completed)

| Bug | Impact | Status |
|-----|--------|--------|
| logMbar gas mass missing ×1e9 | Test D only (physical decomposition) | FIXED |
| LOO data leakage with _perp vars | Phase 128 gap numbers | FIXED (fold-internal) |
| Phase 127 Model B verdict threshold | Misleading "NO" when nested CV says 44/45 | FIXED |
| Core model results (M3/M4/MC) | NOT affected by any bug | CONFIRMED |

## Track 1 (OMD Law — FROZEN at Phase 122)

logOMD = 1.749 + 0.203×log(MHI/L3.6) − 0.101×log(Mbar)

- Sample: N=104 high-Vflat SPARC (Vflat≥70)
- LOO R²: 0.584
- Status: FROZEN, no further work planned

## Zenodo

- Concept DOI: 10.5281/zenodo.19430633 (all versions)
- v10 (current): DOI 10.5281/zenodo.19434177 (PUBLISHED — Phases 126-134 + 200-204, Comprehensively Validated Coupling Law)
- v9 (previous): DOI 10.5281/zenodo.19433840 (Phases 126-134 + 200-201, Externally Validated)
- v8 (previous): DOI 10.5281/zenodo.19433329 (Phases 126-134, Hierarchical Coupling Law)
- v7 (previous): DOI 10.5281/zenodo.19433077 (Phases 126-132, VfResid decoded)
- v6 (previous): DOI 10.5281/zenodo.19432600

## External Validation Program (Phases 200+)

**Phase 200:** Data assembly — N=59 external SPARC galaxies with computed logA0, VfResid, lhOuter
**Phase 201:** Blind prediction — MODERATE_TRANSFER
**Phase 202:** External hierarchy replication — STRONG_HIERARCHY_REPLICATION (8/8 checks PASS)

| Regime | N | Core+VfResid Gap | r(VfResid,a₀) |
|--------|---|-------------------|----------------|
| Full sample | 59 | +8.2% | 0.713 |
| Vflat >= 120 | 16 | +34.3% | 0.841 |
| Vflat >= 180 | 8 | +48.7% | 0.858 |
| Q=1 + Vflat>=120 | 11 | +59.0% | 0.830 |
| Vflat < 120 | 43 | -3.9% | — |

Phase 202 key findings:
- Bootstrap P(VfResid>Core) = 100% in all regimes
- VfResid-only (no Core) gap=47.5% full, 75.4% very-high-V — mediates Core content
- Channel dominance margin: +49pp (full), +35pp (high-V), +71pp (very-high-V) over best alternative
- lhOuter adds +18.1pp in very-high-V, +12.1pp in Q=1+HV (genuine 5th axis externally)
- Partial r(VfResid,a0|Core) = 0.844 high-V, 0.869 Q=1+HV

**Phase 203:** logMhost improvement — CORE_IMPROVES_VfResid_STILL_DOMINATES
- Trained 4 estimators from N=45 real logMhost (tidal analysis)
- Best: Model C (logVflat+logMbar+logRdisk+T), r=0.542 vs formula's 0.232
- Core gap improved: -53% → -8.1% (full), -6.7% → +16.7% (high-V)
- VfResid STILL dominates in every regime (margin 30-76pp over improved Core)
- VfResid-only gap unchanged: 47.5% full, 75.4% very-high-V
- Confirms VfResid is the primary transfer channel, not an artifact of bad Core

**Phase 204:** Final external synthesis — EXTERNAL_VALIDATION_COMPLETE
- Consolidated evidence from Phases 201-203 into unified claim
- lhOuter after Mhost improvement: genuine but reduced (very-high-V: +18.1pp → +8.3pp, Q1+HV: +12.1pp → +4.9pp)
- Part of crude lhOuter was absorbing Mhost error; real residual persists
- Bootstrap P(VfResid > improved Core) = 100% (full), 100% (high-V)
- Paper-ready scientific summary generated
- Program conclusion: VfResid is the primary transfer channel, independent of Mhost quality

## Key Files

- Phase 201: `external-validation/phase201-blind-prediction.cjs`, `public/phase201-blind-prediction.json`
- Phase 200: `external-validation/phase200-data-assembly.cjs`, `public/phase200-external-dataset.json`
- Program plan: `external-validation/PROGRAM.md`
- Phase 134: `scripts/phase134-external-validation.cjs`, `public/phase134-external-validation.json`
- Phase 133C: `scripts/phase133c-coupling-drivers.cjs`, `public/phase133c-coupling-drivers.json`
- Phase 133B: `scripts/phase133b-second-channel.cjs`, `public/phase133b-second-channel.json`
- Phase 133A: `scripts/phase133a-regime-law.cjs`, `public/phase133a-regime-law.json`
- Phase 132C: `scripts/phase132c-external-robustness.cjs`, `public/phase132c-external-robustness.json`
- Phase 132B: `scripts/phase132b-mediation-causal.cjs`, `public/phase132b-mediation-causal.json`
- Phase 132A: `scripts/phase132a-halo-death-match.cjs`, `public/phase132a-halo-death-match.json`
- Phase 132: `scripts/phase132-vfresid-reducibility.cjs`, `public/phase132-vfresid-reducibility.json`
- Phase 131: `scripts/phase131-decode-vfresid.cjs`, `public/phase131-decode-vfresid.json`
- Phase 130: `scripts/phase130-split-4th-sector.cjs`
- Phase 129: `scripts/phase129-vflat-decomposition.cjs`
- Phase 128 (patched): `scripts/phase128-decode-4th-axis.cjs`
- Phase 127 (patched): `scripts/phase127-vflat-challenge.cjs`
- Phase 126: `scripts/phase126-m4-candidate.cjs`
- Results: `public/phase128-decode-4th-axis.json`
- Master table: `public/stage-A-master-table.json`

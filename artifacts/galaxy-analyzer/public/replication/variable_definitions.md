# Variable Definitions — Replication Package
# Galaxy Rotation Curve a0 Variation Analysis
# Final M5 Law (N=45)

**STATUS: SUPERSEDED.** The M5/M3 models below are from Phases 1-100.
Phase 101 showed the underlying Vflat/InnerVmax relation is a geometric artifact.
The current active result (Phases 101-122) uses a completely different framework:
target = logOMD (outer mass discrepancy), predictor = log(MHI/L3.6).
See REPRODUCIBLE_RESULT.md for the current claim.

---

## Target Variable

### logA0
- **Definition**: log10 of the MOND acceleration scale a0 for each galaxy
- **Units**: log10((km/s)^2/kpc)
- **WARNING**: Previous documentation incorrectly stated units as log10(m/s^2).
  To convert: log10(m/s^2) = logA0 + log10(3.241e-14) = logA0 - 13.489
- **Source**: Per-galaxy MOND fits to rotation curves from SPARC database
  (Lelli, McGaugh & Schombert 2016), using the RAR (Radial Acceleration
  Relation) methodology.
- **Typical range**: [3.1, 4.1] in log10((km/s)^2/kpc),
  equivalently [-10.4, -9.4] in log10(m/s^2)

## M5 Predictor Variables (5 axes)

### 1. logMHI
- **Definition**: log10 of total HI gas mass
- **Units**: log10(10^9 M_sun)
- **Source**: SPARC database column `MHI`
- **Computation**: logMHI = log10(MHI_sparc)
- **Corrected M5 coefficient**: -0.235
- **Physical role**: Gas reservoir; more gas suppresses a0

### 2. logMhost
- **Definition**: log10 of the host halo/group mass
- **Units**: log10(M_sun)
- **Source**: Tully (2015) Cosmicflows group catalog; EDD
- **Typical range**: [10.4, 13.0]
- **Corrected M5 coefficient**: -0.175
- **Physical role**: Environmental depth; deeper well suppresses a0

### 3. logSigma0
- **Definition**: log10 of central disk surface mass density
- **Units**: log10(M_sun/pc^2)
- **Computation**: logSigma0 = log10(SBdisk * Upsilon_disk)
- **Corrected M5 coefficient**: +0.146
- **Physical role**: Baryon concentration; denser centers amplify a0

### 4. logMeanRun
- **Definition**: log10 of the mean run length in the rotation curve
  residual sign sequence, measuring kinematic coherence
- **Units**: log10(dimensionless count)
- **Computation**:
  1. Compute residuals: e_i = V_obs(r_i) - V_MOND(r_i)
  2. Record sign sequence: s_i = sign(e_i)
  3. Count runs = contiguous subsequences of same sign
  4. MeanRun = n_points / n_runs
  5. logMeanRun = log10(MeanRun)
- **Typical range**: [0.1, 0.8]
- **Corrected M5 coefficient**: +0.447
- **Physical role**: Dynamical coherence; longer coherent departures
  from MOND amplify a0

### 5. Upsilon_perp (Upsilon-star-perpendicular)
- **Definition**: The component of log10(Upsilon_disk) independent of
  (logMHI, logSigma0, morphT)
- **Units**: dex (residual in log space)
- **Computation**:
  1. Let Y = log10(Upsilon_disk) for each galaxy
  2. Let X = [logMHI, logSigma0, morphT] (3 confounders)
  3. Fit OLS: Y ~ beta_0 + beta_1*logMHI + beta_2*logSigma0 + beta_3*morphT
  4. Upsilon_perp_i = Y_i - Y_hat_i (the OLS residual)
- **CRITICAL NOTE**: morphT=0 is a valid value (S0 galaxies). Use nullish
  coalescing (?? 5) NOT logical OR (|| 5) when providing defaults, otherwise
  T=0 galaxies (NGC4138, UGC06786) are incorrectly assigned T=5.
- **Upsilon_disk source**: Li et al. (2020) Table 1 SPS estimates.
  Default: 0.50 M_sun/L_sun where unavailable.
- **Corrected M5 coefficient**: +0.372
- **Physical role**: Independent stellar structure; higher M/L amplifies a0

## Dropped Variable: rcWiggliness
- **Definition**: RMS fractional deviation of RC from its smoothed trend
- **Status**: DROPPED from M5 (retained in dataset for replication)
- **Reason**: F-test = 1.50 (not significant), 14.5% bootstrap sign
  flips, removing it improves LOO. Kinematic information is better
  captured by logMeanRun alone.

## Auxiliary Variables (not in M5, provided for context)

### morphT
- **Definition**: de Vaucouleurs numerical morphological type
- **Source**: SPARC database column `T`
- **CRITICAL**: T=0 is valid (S0 galaxies). Do NOT use || 5 default.
- **Role**: Confounder in Upsilon_perp construction

### Vflat_km_s
- **Definition**: Flat rotation velocity (km/s)
- **WARNING**: Do NOT use as a predictor of a0 (circular variable)

### D_Mpc, inc_deg, RHI_kpc, Rdisk_kpc, L36_1e9Lsun, Q_SPARC
- Context/metadata only

## Corrected M5 Full Equation

log(a0) = 5.016 - 0.235*logMHI - 0.175*logMhost + 0.146*logSigma0
         + 0.447*logMeanRun + 0.372*Upsilon_perp

LOO gap% = 46.6% | LOO RMS = 0.189 dex | N = 45

## M3 Compressed State Law (unaffected by T bug)

log(a0) = 5.182 - 0.198*logMHI - 0.155*logMhost + 0.459*logMeanRun

LOO gap% = 44.1% | LOO RMS = 0.193 dex | N = 45
Retains 95% of corrected M5 signal

## Errata

**T||5 bug (fixed)**: Prior to correction, scripts used `sparcMap[g.name]?.T || 5`
which treated morphological type T=0 (valid for S0 galaxies NGC4138 and UGC06786)
as falsy, replacing it with T=5. This inflated the Upsilon_perp coefficient from
0.372 to 0.658 and the LOO gap% from 46.6% to 50.9%. Fixed by using `?? 5`
(nullish coalescing). M3 was unaffected (does not use morphT).

## Sample Selection Criteria (N=45)

1. Galaxy is in SPARC database (175 galaxies)
2. Has a valid per-galaxy a0 fit
3. Has a published, peer-reviewed distance estimate (not Hubble-flow only)
4. Has a host mass assignment from Tully group catalog
5. Has >= 5 RC points for kinematic diagnostics

# Variable Definitions — Replication Package
# Galaxy Rotation Curve a0 Variation Analysis
# Final M5 Law (N=45)

## Target Variable

### logA0
- **Definition**: log10 of the MOND acceleration scale a0 for each galaxy
- **Units**: log10(m/s^2)
- **Source**: Per-galaxy MOND fits to rotation curves from SPARC database
  (Lelli, McGaugh & Schombert 2016), using the RAR (Radial Acceleration
  Relation) methodology.
- **Typical range**: [-10.8, -9.4] (i.e., a0 ~ 0.4-4.0 x 10^-10 m/s^2)

## M5 Predictor Variables (5 axes)

### 1. logMHI
- **Definition**: log10 of total HI gas mass
- **Units**: log10(10^9 M_sun)
- **Source**: SPARC database column `MHI`
- **Computation**: logMHI = log10(MHI_sparc)
- **M5 coefficient**: -0.236 (rational approx: -1/4)
- **Physical role**: Gas reservoir; more gas suppresses a0

### 2. logMhost
- **Definition**: log10 of the host halo/group mass
- **Units**: log10(M_sun)
- **Source**: Tully (2015) Cosmicflows group catalog; EDD
- **Typical range**: [10.4, 13.0]
- **M5 coefficient**: -0.172 (rational approx: -1/6)
- **Physical role**: Environmental depth; deeper well suppresses a0

### 3. logSigma0
- **Definition**: log10 of central disk surface mass density
- **Units**: log10(M_sun/pc^2)
- **Computation**: logSigma0 = log10(SBdisk * Upsilon_disk)
- **M5 coefficient**: +0.145 (rational approx: +1/7)
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
- **M5 coefficient**: +0.452 (rational approx: +3/7)
- **Physical role**: Dynamical coherence; longer coherent departures
  from MOND amplify a0
- **2D validation**: Correlates |r|=0.83 with THINGS velocity-field
  diagnostics (lopsidedness, bisymmetric flow, non-circular motion)

### 5. Upsilon_perp (Upsilon-star-perpendicular)
- **Definition**: The component of log10(Upsilon_disk) independent of
  (logMHI, logSigma0, morphT)
- **Units**: dex (residual in log space)
- **Computation**:
  1. Let Y = log10(Upsilon_disk) for each galaxy
  2. Let X = [logMHI, logSigma0, morphT] (3 confounders)
  3. Fit OLS: Y ~ beta_0 + beta_1*logMHI + beta_2*logSigma0 + beta_3*morphT
  4. Upsilon_perp_i = Y_i - Y_hat_i (the OLS residual)
- **Upsilon_disk source**: Li et al. (2020) Table 1 SPS estimates.
  Default: 0.50 M_sun/L_sun where unavailable.
- **Orthogonalization**: Removes 68% of Upsilon variance
- **M5 coefficient**: +0.658 (rational approx: +2/3)
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
- **Role**: Confounder in Upsilon_perp construction

### Vflat_km_s
- **Definition**: Flat rotation velocity (km/s)
- **WARNING**: Do NOT use as a predictor of a0 (circular variable)

### D_Mpc, inc_deg, RHI_kpc, Rdisk_kpc, L36_1e9Lsun, Q_SPARC
- Context/metadata only

## M5 Full Equation

log(a0) = 4.978 - 0.236*logMHI - 0.172*logMhost + 0.145*logSigma0
         + 0.452*logMeanRun + 0.658*Upsilon_perp

LOO gap% = 51.0% | RMS = 0.157 dex | N = 45

## M3 Compressed State Law

log(a0) = 5.182 - 0.198*logMHI - 0.155*logMhost + 0.459*logMeanRun

LOO gap% = 44.1% | RMS = 0.193 dex | N = 45
Retains 86% of M5 signal with 3 axes

## Sample Selection Criteria (N=45)

1. Galaxy is in SPARC database (175 galaxies)
2. Has a valid per-galaxy a0 fit
3. Has a published, peer-reviewed distance estimate (not Hubble-flow only)
4. Has a host mass assignment from Tully group catalog
5. Has >= 5 RC points for kinematic diagnostics

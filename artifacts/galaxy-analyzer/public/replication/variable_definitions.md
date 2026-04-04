# Variable Definitions — Replication Package
# Galaxy Rotation Curve a0 Variation Analysis
# Frozen Baseline B″ (N=45)

## Target Variable

### logA0
- **Definition**: log10 of the MOND acceleration scale a0 for each galaxy
- **Units**: log10(m/s^2)
- **Source**: Per-galaxy MOND fits to rotation curves from SPARC database
  (Lelli, McGaugh & Schombert 2016), using the RAR (Radial Acceleration
  Relation) methodology. Each galaxy's a0 is the best-fit acceleration
  parameter that minimizes residuals between observed and MOND-predicted
  rotation curves.
- **Typical range**: [-10.8, -9.4] (i.e., a0 ~ 0.4–4.0 x 10^-10 m/s^2)

## B″ Predictor Variables (6 axes)

### 1. logMHI
- **Definition**: log10 of total HI gas mass
- **Units**: log10(10^9 M_sun)
- **Source**: SPARC database column `MHI` (already in units of 10^9 M_sun)
- **Computation**: logMHI = log10(MHI_sparc)
- **log base**: 10

### 2. rcWiggliness
- **Definition**: Root-mean-square fractional deviation of the observed
  rotation curve from its own smoothed trend, measuring small-scale
  irregularity in the RC shape.
- **Units**: dimensionless (fractional)
- **Computation**:
  1. Take observed RC velocity points V_obs(r) at radii r_1, ..., r_n
  2. Compute a 3-point running mean: V_smooth(r_i) = mean(V_obs(r_{i-1}),
     V_obs(r_i), V_obs(r_{i+1})) for i = 2..n-1
  3. Compute fractional residuals: delta_i = (V_obs(r_i) - V_smooth(r_i))
     / V_smooth(r_i)
  4. rcWiggliness = sqrt(mean(delta_i^2)) over interior points
- **Source**: Computed from SPARC rotation curve data files
- **Typical range**: [0.01, 0.15]
- **Note**: Higher values indicate more irregular/noisy RC shapes

### 3. logMhost
- **Definition**: log10 of the host halo/group mass from the galaxy's
  large-scale environment
- **Units**: log10(M_sun)
- **Source**: Tully (2015) Cosmicflows group catalog; EDD (Extragalactic
  Distance Database) environmental assignments. For field galaxies, this
  is the galaxy's own estimated halo mass. For group members, this is
  the group total mass.
- **Matching**: Each SPARC galaxy matched to its Tully group by name
  cross-reference with tolerance for naming conventions (NGC/UGC/etc.)
- **Typical range**: [10.4, 13.0]
- **Note**: Already in B″; captures environmental depth

### 4. logSigma0
- **Definition**: log10 of central disk surface mass density
- **Units**: log10(M_sun/pc^2)
- **Source**: SPARC database column `SBdisk` (disk central surface
  brightness in L_sun/pc^2), converted to mass density using the
  galaxy's stellar mass-to-light ratio
- **Computation**: logSigma0 = log10(SBdisk * Upsilon_disk)
  where Upsilon_disk is the 3.6um mass-to-light ratio
- **log base**: 10

### 5. logMeanRun
- **Definition**: log10 of the mean run length in the rotation curve
  residual sign sequence, measuring kinematic coherence
- **Units**: log10(dimensionless count)
- **Computation**:
  1. Compute residuals of observed RC from MOND fit: e_i = V_obs(r_i)
     - V_MOND(r_i)
  2. Record the sign sequence: s_i = sign(e_i)
  3. Count "runs" = contiguous subsequences of same sign
  4. MeanRun = n_points / n_runs
  5. logMeanRun = log10(MeanRun)
- **Source**: Computed from SPARC rotation curve data + MOND fits
- **Typical range**: [0.1, 0.8]
- **Interpretation**: Higher values = longer coherent stretches of
  positive or negative residuals = more systematic departure from MOND

### 6. Upsilon_perp (Upsilon-star-perpendicular)
- **Definition**: The component of log10(Upsilon_disk) that is
  statistically independent of (logMHI, logSigma0, morphT).
  This isolates the "unexplained" part of the stellar mass-to-light
  ratio after removing what gas mass, surface density, and morphology
  already predict.
- **Units**: dex (residual in log space)
- **Computation**:
  1. Let Y = log10(Upsilon_disk) for each galaxy
  2. Let X = [logMHI, logSigma0, morphT] (3 confounders)
  3. Fit OLS: Y ~ beta_0 + beta_1*logMHI + beta_2*logSigma0 +
     beta_3*morphT
  4. Upsilon_perp_i = Y_i - Y_hat_i (the OLS residual)
- **Upsilon_disk source**: Li et al. (2020) Table 1 for SPARC galaxies
  with published 3.6um SPS-based mass-to-light ratios. Where not
  available, default Upsilon_disk = 0.50 M_sun/L_sun.
- **Note**: This orthogonalization removes 68% of Upsilon variance,
  ensuring Upsilon_perp carries genuinely independent M/L information

## Auxiliary Variables (not in B″, provided for context)

### morphT
- **Definition**: de Vaucouleurs numerical morphological type
- **Source**: SPARC database column `T`
- **Range**: -3 (E) to 11 (Irr)
- **Role**: Confounder in Upsilon_perp construction (not a direct B″ predictor)

### Vflat_km_s
- **Definition**: Flat rotation velocity
- **Units**: km/s
- **Source**: SPARC database
- **WARNING**: Do NOT use as a predictor of a0. Vflat is mechanically
  linked to a0 (circular variable).

### D_Mpc, inc_deg, RHI_kpc, Rdisk_kpc, L36_1e9Lsun, Q_SPARC
- **Definition**: Distance, inclination, HI radius, disk scale length,
  3.6um luminosity, SPARC quality flag
- **Source**: SPARC database
- **Role**: Context/metadata only

## Sample Selection Criteria (N=45)

### Inclusion requires ALL of:
1. Galaxy is in SPARC database (175 galaxies)
2. Has a valid per-galaxy a0 fit (delta_a0 defined)
3. Has a published, peer-reviewed distance estimate (not just SPARC
   default) — specifically from TRGB, Cepheid, or other primary
   distance indicator, OR from a well-calibrated secondary method
   with published uncertainty
4. Has a host mass assignment from Tully group catalog or equivalent
5. Has sufficient rotation curve quality for rcWiggliness and
   logMeanRun computation (>= 5 RC points)

### The "published distance quality" criterion is the primary filter
   that reduces 175 → 56 → 45 galaxies. It ensures distance-dependent
   quantities (a0, MHI, sizes) are not dominated by Hubble-flow
   distance errors.

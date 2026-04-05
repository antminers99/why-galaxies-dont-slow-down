# Literature Summary — 22 Key Papers for Galaxy Rotation Curve Analysis

## How This Document Relates to Our Project

Our project found that gas-to-stellar balance, quantified by log(MHI/L3.6), is the strongest independent predictor of outer support requirement (outer mass discrepancy) in high-Vflat SPARC galaxies (Vflat >= 70 km/s, N=104). This document summarizes the 22 most relevant papers that form the scientific context for this finding.

---

## TIER 1: ESSENTIAL FOUNDATION

---

### Paper 01 — Lelli, McGaugh & Schombert (2016a)
**"SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves"**
AJ 152, 157 | arXiv: 1606.09251

**What they did:** Built the SPARC (Spitzer Photometry and Accurate Rotation Curves) database — a collection of 175 disk galaxies with high-quality HI/H-alpha rotation curves combined with Spitzer 3.6 micron photometry. For each galaxy, they provide: observed rotation curves, baryonic mass models decomposed into stellar disk, bulge, and gas components, plus global properties (luminosity L3.6, disk scale length Rdisk, HI mass MHI, flat velocity Vflat, distance, inclination).

**Key numbers:**
- 175 galaxies spanning Vflat from ~20 to ~300 km/s
- Stellar masses from ~10^7 to ~10^11 solar masses
- Gas fractions from near 0 to >90%
- Spitzer 3.6 micron photometry minimizes dust effects and gives more reliable stellar mass estimates than optical bands
- Recommended mass-to-light ratio: Upsilon_disk = 0.5 M_sun/L_sun at 3.6 micron

**Why it matters for us:** This is our entire dataset. Every number in our project comes from SPARC. The quality of the photometry (Spitzer 3.6 micron) and the careful rotation curve compilation make this the gold standard for galaxy dynamics studies. Our use of L3.6 as a stellar mass proxy and MHI as a gas mass proxy follows directly from this paper's methodology.

---

### Paper 02 — McGaugh, Lelli & Schombert (2016)
**"Radial Acceleration Relation in Rotationally Supported Galaxies"**
PRL 117, 201101 | arXiv: 1609.05917

**What they did:** Discovered the Radial Acceleration Relation (RAR) — a tight empirical relation between the observed centripetal acceleration (gobs = Vobs^2/r) and the acceleration predicted from baryonic matter alone (gbar) at every point in every galaxy. Using 2693 data points from 153 SPARC galaxies, they showed that gobs is a single-valued function of gbar with very small scatter (~0.13 dex).

**The RAR formula:**
gobs = gbar / (1 - exp(-sqrt(gbar/a0)))

where a0 ~ 1.2 × 10^-10 m/s^2 is a characteristic acceleration scale.

**Key findings:**
- At high accelerations (gbar >> a0): gobs ≈ gbar (baryons account for everything)
- At low accelerations (gbar << a0): gobs ≈ sqrt(gbar × a0) (large mass discrepancy)
- The transition is smooth and follows a single function
- Scatter is remarkably small — comparable to observational uncertainties
- This is consistent with MOND but does not require it

**Why it matters for us:** The RAR is the most important empirical relation in the field. Our Phase 113 showed that our finding (logOMD vs gas state) is 65% independent of RAR residuals — meaning we capture genuinely new information that the RAR does not. Our Phase 116 showed that the acceleration scale a0 is roughly universal across gas states, with only a mild hint that gas-rich galaxies prefer slightly higher a0.

---

### Paper 03 — Lelli, McGaugh & Schombert (2017)
**"One Law to Rule Them All: The Radial Acceleration Relation of Galaxies"**
ApJ 836, 152 | arXiv: 1610.08981

**What they did:** Extended the RAR analysis with the full SPARC sample (175 galaxies, 2693 points), provided a more thorough analysis of the scatter, tested for residual correlations with galaxy properties, and explored the implications for dark matter and modified gravity.

**Key findings beyond the PRL paper:**
- Intrinsic scatter of the RAR is ≤ 0.057 dex (after accounting for observational errors)
- Residuals do NOT correlate strongly with galaxy properties (size, surface brightness, gas fraction, morphology)
- The relation holds from massive spirals down to gas-dominated dwarfs
- Both Upsilon_disk = 0.5 and 0.2 M_sun/L_sun give similarly tight relations
- Dark matter halo properties are fully determined by baryonic distribution

**Critical detail for us:** The claim that RAR residuals do NOT correlate with gas fraction (their Section 5.2) is important context. Our Phase 113 found weak but nonzero correlations (r~0.2) — consistent with their finding of "no strong correlation" but showing that a subtle signal does exist when measured carefully. The key difference is that we look at galaxy-level aggregated residuals, not point-by-point.

---

### Paper 04 — Lelli, McGaugh & Schombert (2016b)
**"The Relation between Stellar and Dynamical Surface Densities in the Central Regions of Disk Galaxies"**
ApJ 827, L19 | arXiv: 1610.08980

**What they did:** Established a tight relation between central stellar surface density (Sigma_star,0) and central dynamical surface density (Sigma_dyn,0 derived from inner rotation curve gradient) in disk galaxies. This is essentially a 1D version of the RAR applied specifically to galaxy centers.

**Key findings:**
- For high surface brightness galaxies (Sigma_star,0 > ~100 M_sun/pc^2): Sigma_dyn ≈ Sigma_star (baryons dominate the center)
- For low surface brightness galaxies: Sigma_dyn >> Sigma_star (large central mass discrepancy)
- The transition follows the same acceleration scale a0 as the RAR
- Scatter is very small (~0.12 dex)

**Why it matters for us:** This is directly relevant to our Sigma0 analysis. Our Phases 105 and 108 showed that Sigma0 (stellar surface density) is fully mediated by gas fraction — it has NO independent predictive power for outer mass discrepancy. This paper shows that Sigma0 is important for CENTRAL dynamics, but our finding shows it is NOT independently important for OUTER dynamics once gas fraction is controlled.

---

### Paper 05 — Oh et al. (2015)
**"High-resolution Mass Models of Dwarf Galaxies from LITTLE THINGS"**
AJ 149, 180 | arXiv: 1502.01281

**What they did:** Presented high-resolution HI rotation curves and mass models for 26 dwarf irregular galaxies from the LITTLE THINGS survey. These are among the most carefully measured dwarf galaxy rotation curves available.

**Key findings:**
- Dwarf galaxies show slowly rising rotation curves (unlike the steeply rising curves of massive spirals)
- Dark matter dominates at ALL radii in most dwarfs (unlike massive galaxies where baryons dominate the inner regions)
- The dark matter density profiles are better fit by pseudo-isothermal (cored) halos than NFW (cuspy) halos — the "core-cusp problem"
- Inner mass density slopes range from -0.2 to -0.4 (much shallower than the NFW prediction of -1)
- Gas often dominates the baryonic budget in dwarfs

**Why it matters for us:** This explains why our regime boundary exists. Below Vflat ~70 km/s, dwarf galaxies enter a different physical regime where dark matter dominates everywhere, rotation curves are slowly rising, and the simple gas-fraction predictor breaks down. Our Phase 111 confirmed this gradual transition. We could not compute logOMD for LITTLE THINGS galaxies because they lack compatible baryonic decomposition in our pipeline.

---

## TIER 2: DIVERSITY AND THE OPEN PROBLEM

---

### Paper 06 — Oman et al. (2015)
**"The unexpected diversity of dwarf galaxy rotation curves"**
MNRAS 452, 3650 | arXiv: 1504.01437

**What they did:** Compared observed rotation curve shapes of dwarf galaxies with predictions from LCDM cosmological simulations (EAGLE, APOSTLE). They found that at fixed maximum velocity, simulated galaxies predict a narrow range of rotation curve shapes, while observed galaxies show a much wider range — the "diversity problem."

**Key findings:**
- At Vmax ~ 50 km/s, observed V(2 kpc) ranges from ~10 to ~50 km/s
- Simulations predict V(2 kpc) ~ 30-40 km/s at this mass (too narrow)
- Some observed dwarfs rise too slowly (too cored) and others too quickly (too concentrated) compared to simulations
- This diversity challenges LCDM predictions at small scales
- The problem is most acute for dwarfs, less so for massive spirals

**Why it matters for us:** The diversity problem is one of the biggest open questions in the field. Our Phase 114 tested whether gas-to-stellar balance can explain some of this diversity. We found that gas state predicts rotation curve shape metrics (outer slope r=0.46, concentration r=-0.50), but most of this signal is mediated by Vflat. The one exception — baryon dominance fraction — survives partial correlation with Vflat (pr=-0.50), connecting our finding to the diversity problem.

---

### Paper 07 — Santos-Santos et al. (2020)
**"Baryonic clues to the puzzling diversity of dwarf galaxy rotation curves"**
MNRAS 495, 58 | arXiv: 1911.09116

**What they did:** Investigated whether baryonic properties (gas content, stellar mass, size, star formation rate) can explain the diversity of dwarf galaxy rotation curves. They used both observed galaxies and hydrodynamical simulations (NIHAO).

**Key findings:**
- Galaxies with higher stellar-to-halo mass ratios tend to have more concentrated (faster-rising) rotation curves
- Gas-rich galaxies tend to have more slowly rising curves
- The baryonic content can account for SOME but not ALL of the diversity
- Stellar feedback and gas outflows play a role in shaping the inner dark matter profile
- The correlation between baryonic properties and rotation curve shape is present but has significant scatter

**Why it matters for us:** This paper is the closest existing work to our finding. They show that baryonic properties correlate with rotation curve shape diversity, which is exactly what our Phase 114 confirms. The key difference: they focus on dwarfs and inner rotation curve shapes, while we focus on high-Vflat galaxies and outer mass discrepancy. Our finding adds the specific result that log(MHI/L3.6) → logOMD survives independently of the BTFR (Phase 115), which they did not test.

---

### Paper 08 — Marasco et al. (2020)
**"Reconciling the Diversity and Uniformity of Galactic Rotation Curves with Self-Interacting Dark Matter"**
(check exact title) | arXiv: 2003.07005

**What they did:** Explored whether self-interacting dark matter (SIDM) can simultaneously explain the uniformity of the RAR and the diversity of rotation curve shapes. They used semi-analytical models to test how SIDM thermalizes with the baryonic distribution.

**Key findings:**
- SIDM with a cross-section of ~1-3 cm^2/g can produce diverse inner rotation curves while maintaining a tight RAR
- The diversity arises from different baryonic concentrations: galaxies with more concentrated baryons compress the SIDM core more, leading to faster-rising curves
- The uniformity (RAR) arises because the SIDM equilibrium is determined by the baryonic potential
- This provides a framework where both diversity and uniformity coexist

**Why it matters for us:** This gives a theoretical context for our finding. If baryonic concentration determines the dark matter response, then the gas-to-stellar ratio (which correlates with baryonic concentration) would naturally predict outer mass discrepancy. Our finding could be an empirical manifestation of this SIDM-baryon coupling — but we emphasize this is speculative and our result is purely empirical.

---

## TIER 3: GAS FRACTION AND SCALING RELATIONS

---

### Paper 09 — Catinella et al. (2010)
**"The GALEX Arecibo SDSS Survey. I. Gas fraction scaling relations of massive galaxies"**
MNRAS 403, 683 | arXiv: 0912.1610

**What they did:** Measured HI gas fractions (MHI/Mstar) for a representative sample of ~190 massive galaxies (Mstar > 10^10 M_sun) from the GASS survey, combining Arecibo HI observations with SDSS optical and GALEX UV photometry.

**Key findings:**
- Gas fraction decreases strongly with increasing stellar mass: MHI/Mstar ∝ Mstar^(-0.35)
- Gas fraction decreases with increasing stellar surface mass density (mu_star)
- Gas fraction correlates strongly with NUV-r color (a proxy for specific star formation rate)
- At fixed stellar mass, bluer galaxies have higher gas fractions
- The scaling relations have significant scatter (~0.5 dex)

**Why it matters for us:** This establishes that gas fraction is tightly linked to galaxy properties — it is NOT a random variable. The strong correlation with stellar surface density is relevant because our Phase 108 showed that Sigma0 is fully mediated by fgas. Catinella et al. show WHY this happens: gas fraction and surface density are correlated in the general galaxy population. Our contribution is showing that fgas, not Sigma0, is the causally prior variable for predicting outer dynamics.

---

### Paper 10 — Catinella et al. (2018)
**"xGASS: total cold gas scaling relations and molecular-to-atomic gas ratios of galaxies in the local Universe"**
MNRAS 476, 875 | arXiv: 1802.02373

**What they did:** Extended the GASS survey to lower masses (Mstar > 10^9 M_sun) and added molecular gas (CO) measurements, creating the xGASS+xCOLD GASS survey of ~1200 galaxies with complete cold gas census.

**Key findings:**
- Total gas fraction (MHI + MH2)/Mstar scales with Mstar, mu_star, and NUV-r at all masses
- Molecular-to-atomic ratio (MH2/MHI) increases with stellar surface density and decreases with specific SFR
- Below Mstar ~ 10^9.5 M_sun, galaxies are predominantly atomic-gas dominated
- The scaling relations are well-described by simple planes in (gas fraction, Mstar, SFR) space
- Environmental effects are secondary to internal galaxy properties

**Why it matters for us:** xGASS provides the broader context for gas fraction scaling. Our SPARC sample covers a similar mass range. The finding that HI dominates at lower masses explains why our log(MHI/L3.6) variable (which uses HI mass only) works well — for most SPARC galaxies, MHI is the dominant gas component. The molecular gas (H2) contribution we miss may contribute to the ~50% unexplained variance in our model.

---

### Paper 11 — Huang et al. (2012)
**"The Arecibo Legacy Fast ALFA Survey: The Galaxy Population Detected by ALFALFA"**
ApJ 756, 113 | arXiv: 1207.0005

**What they did:** Characterized the galaxy population detected by the ALFALFA 21cm survey — the largest blind HI survey at the time, covering 2800 deg^2 with ~15,000 HI detections. They measured the HI mass function, gas-to-stellar mass ratios, and correlations with optical properties.

**Key findings:**
- The HI mass function extends from MHI ~ 10^7 to 10^10.8 M_sun
- Gas richness (MHI/Mstar) increases toward lower stellar masses, later morphological types, and bluer colors
- High gas-fraction galaxies (MHI/Mstar > 1) are common at low masses but rare at high masses
- The "gas sequence" (MHI vs Mstar) has a slope shallower than unity: low-mass galaxies are disproportionately gas-rich
- HI-selected samples are biased toward gas-rich, star-forming, late-type galaxies

**Why it matters for us:** ALFALFA confirms that gas-to-stellar ratio is a fundamental galaxy property that varies systematically across the galaxy population. Our SPARC sample, being rotation-curve selected, has a different selection function than ALFALFA (HI-selected), but the underlying gas scaling relations are consistent. This supports the physical reality of our log(MHI/L3.6) predictor as a meaningful variable.

---

### Paper 12 — Cortese et al. (2011)
**"The effect of the environment on the HI scaling relations"**
MNRAS 415, 1797 | arXiv: 1107.3657

**What they did:** Studied how environment (field vs group vs cluster) affects the HI gas content and scaling relations of galaxies, using the Herschel Reference Survey combined with HI data.

**Key findings:**
- Cluster galaxies are systematically HI-deficient compared to field galaxies at the same stellar mass and morphology
- The HI-to-stellar mass relation shifts downward by ~0.5 dex in dense environments
- Gas stripping primarily removes HI from the outer disk, leaving molecular gas in the center relatively unaffected
- After correcting for environment, the residual scatter in gas fraction scaling relations is reduced
- Environment is a "second parameter" after internal properties

**Why it matters for us:** Our SPARC sample contains mostly field galaxies, so environmental gas stripping is not a major concern. However, this paper warns that if we ever extend to cluster samples, environmental effects could alter the gas-to-stellar ratio without changing the underlying dynamics, potentially weakening our correlation. This is a limitation to note if cross-survey replication is attempted.

---

### Paper 13 — Boselli et al. (2014)
**"The effect of structure and star formation on the gas content of nearby galaxies"**
A&A 564, A67 | arXiv: 1312.4596

**What they did:** Investigated how the structure (surface brightness, concentration) and star formation properties of galaxies relate to their gas content, using multi-wavelength data for the Herschel Reference Survey.

**Key findings:**
- Gas fraction correlates with stellar surface density, specific star formation rate, and morphological type
- At fixed stellar mass, low surface brightness galaxies have higher gas fractions
- Star formation efficiency (SFR/Mgas) varies by less than the gas fraction itself
- Gas fraction is the single most important parameter distinguishing different galaxy types in gas/structural property space
- Structure and gas content are tightly linked but the causal direction is unclear

**Why it matters for us:** This paper directly addresses our central question: is gas fraction fundamental, or is it driven by structural properties (like Sigma0)? Boselli et al. show that gas fraction and surface density are tightly correlated but cannot determine causality. Our Phases 105 and 108 go further by showing that for outer dynamics prediction, fgas is the causally prior variable (Sigma0 dies after controlling fgas, not vice versa). This is a genuinely new contribution.

---

## TIER 4: TULLY-FISHER AND BARYONIC RELATIONS

---

### Paper 14 — Verheijen (2001)
**"The Ursa Major Cluster of Galaxies. V. HI Rotation Curve Shapes and the Tully-Fisher Relations"**
ApJ 563, 694 | arXiv: astro-ph/0108225

**What they did:** Studied HI rotation curves and Tully-Fisher relations for 49 galaxies in the Ursa Major cluster, carefully distinguishing between different velocity measures: V_max (peak velocity), V_flat (velocity in the flat part), and V_2Rd (velocity at 2 disk scale lengths).

**Key findings:**
- V_flat gives a tighter Tully-Fisher relation than V_max or line width measures
- Rotation curve shapes vary systematically: massive galaxies show declining rotation curves (Vmax > Vflat), while less massive galaxies show rising or flat curves
- The ratio Vmax/Vflat is correlated with luminosity and surface brightness
- The K-band Tully-Fisher relation has the smallest scatter (~0.3 mag)
- Gas-rich galaxies follow the Tully-Fisher relation when gas mass is included (Baryonic TFR)

**Why it matters for us:** Verheijen's distinction between Vmax and Vflat is directly relevant to our original (superseded) ratio law (Vflat/InnerVmax). Phase 101 showed that this ratio is largely a geometric artifact. However, Verheijen's finding that Vmax/Vflat correlates with luminosity is confirmed by our Phase 114, where we find the same correlation (r=0.50) — but it is mediated by Vflat itself.

---

### Paper 15 — Bell & de Jong (2001)
**"Stellar Mass-to-Light Ratios and the Tully-Fisher Relation"**
ApJ 550, 212 | arXiv: astro-ph/0011493

**What they did:** Explored how stellar mass-to-light ratios (Upsilon) vary with galaxy color, and how this affects the Tully-Fisher relation. They used stellar population synthesis models to derive Upsilon as a function of broadband colors.

**Key findings:**
- Upsilon varies by factors of 2-7 depending on the optical band used (B-band worst, K-band best)
- Near-infrared bands (K, 3.6 micron) minimize Upsilon variations: typical range 0.3-0.8 M_sun/L_sun
- The stellar mass Tully-Fisher relation (using color-corrected Upsilon) is tighter than luminosity TFR
- Assuming a universal Upsilon introduces systematic residuals in the TFR correlated with color
- "Maximum disk" solutions are not always physically justified

**Why it matters for us:** We use a fixed Upsilon_disk = 0.5 at 3.6 micron following SPARC conventions. Bell & de Jong show this is reasonable at 3.6 micron where Upsilon variations are minimized. However, galaxy-to-galaxy Upsilon variations of even 0.3-0.8 could contribute to our logOMD scatter. Our Phase 110 Monte Carlo test perturbed Upsilon by ±25% and confirmed the signal is robust (r CI: [0.684, 0.736]).

---

### Paper 16 — Schombert, McGaugh & Lelli (2020)
**"Using the Baryonic Tully-Fisher Relation to Measure H0"**
AJ 160, 71 | arXiv: 2006.08615

**What they did:** Used the Baryonic Tully-Fisher Relation (BTFR) as a distance indicator to measure the Hubble constant H0, using SPARC galaxies with independent distance measurements as calibrators.

**Key findings:**
- BTFR: Mbar = A × Vflat^4, with slope = 3.85 ± 0.09 (close to 4)
- The BTFR is one of the tightest galaxy scaling relations: scatter ~0.10-0.15 dex in Mbar
- They derive H0 = 75.1 ± 2.3 km/s/Mpc (consistent with local distance ladder, higher than Planck CMB)
- The BTFR slope of ~4 is naturally expected in MOND but requires fine-tuning in LCDM
- Gas-dominated galaxies follow the same BTFR as stellar-dominated ones

**Why it matters for us:** Our Phase 115 tested whether gas state predicts BTFR residuals. We found it does NOT (r=0.056). This is consistent with Schombert et al.'s finding that the BTFR is tight regardless of gas fraction. Crucially, our fgas → logOMD signal survives after controlling for BTFR residuals (partial r=0.787), proving our finding operates on a completely independent axis from the BTFR.

---

### Paper 17 — Lelli, McGaugh & Schombert (2019)
**"The Baryonic Tully-Fisher Relation. I. A Tight Relation for Isolated Galaxies"**
(BTFR calibration) | arXiv: 1901.04966

**What they did:** Refined the BTFR using a carefully selected sample of isolated galaxies with well-measured Vflat values, removing possible environmental effects and focusing on the highest-quality data.

**Key findings:**
- The BTFR for isolated galaxies has an observed scatter of ~0.08 dex in Mbar
- After accounting for observational errors, the intrinsic scatter may be consistent with zero
- Slope = 3.85 ± 0.09 (consistent with earlier measurements)
- No evidence for a break or curvature over 5 decades in baryonic mass
- The relation is independent of galaxy surface brightness, size, or gas fraction

**Why it matters for us:** The zero intrinsic scatter claim for the BTFR is remarkable and means that Vflat alone determines Mbar with essentially no freedom. Since our signal (fgas → logOMD) is independent of BTFR position (Phase 115), it means gas-to-stellar balance provides a SECOND, independent axis of information about galaxy dynamics beyond total mass. This is the strongest statement of novelty for our finding.

---

### Paper 18 — McGaugh & Schombert (2020)
**"The Baryonic Tully-Fisher Relation. II. Stellar Mass-to-Light Ratios"**
(check exact ref) | arXiv: 2005.01228

**What they did:** Investigated how the assumed stellar mass-to-light ratio (Upsilon) affects the BTFR, testing different prescriptions from stellar population synthesis and comparing to dynamical constraints.

**Key findings:**
- The BTFR slope is sensitive to the assumed Upsilon: higher Upsilon → steeper slope
- Upsilon = 0.5 M_sun/L_sun at 3.6 micron gives the most self-consistent BTFR
- Color-dependent Upsilon prescriptions do not significantly improve the BTFR scatter
- The BTFR naturally selects a specific Upsilon value, providing an independent constraint on stellar mass models
- Galaxy-to-galaxy Upsilon variations are small at 3.6 micron

**Why it matters for us:** This validates our choice of Upsilon_disk = 0.5 M_sun/L_sun at 3.6 micron. The BTFR "selects" this value independently, meaning our stellar mass estimates are consistent with the dynamical constraints. Any systematic error in Upsilon would affect both the gas fraction (through Mstar = Upsilon × L3.6) and the outer mass discrepancy (through Vbar), but Phase 110 showed our signal is robust to ±25% Upsilon perturbations.

---

## TIER 5: CURRENT DEBATE AND FALSIFICATION

---

### Paper 19 — Rodrigues et al. (2018)
**"Absence of a fundamental acceleration scale in galaxies"**
Nature Astronomy 2, 668 | arXiv: 1806.06803

**What they did:** Challenged the universality of the RAR by fitting individual galaxies with their own a0 values and testing whether a single a0 is statistically preferred over galaxy-specific values.

**Key findings:**
- When fitting individual galaxies, the best-fit a0 varies by more than a factor of 5 between galaxies
- A Bayesian model comparison finds that galaxy-specific a0 values are preferred over a universal a0 by a significant margin
- This challenges MOND-like interpretations that require a universal acceleration scale
- The authors argue the RAR tightness may be a consequence of galaxy formation physics rather than fundamental physics

**Why it matters for us:** Our Phase 116 independently found that best-fit a0 varies across galaxies (SD = 0.302 dex, range from -10.57 to -9.23). We additionally tested whether gas state predicts this variation and found only weak correlation (r=0.17). The tercile analysis hints that gas-rich galaxies prefer slightly higher a0 (+0.11 dex), but this is not a strong effect. Rodrigues et al.'s finding supports our conclusion that a0 variation is not strongly driven by gas state.

---

### Paper 20 — Banik et al. (2024)
**"On the tension between the radial acceleration relation and Solar system quadrupole in modified gravity MOND"**
(check exact title) | arXiv: 2310.07646

**What they did:** Tested whether MOND can simultaneously satisfy the RAR (galaxy scales) and Solar system constraints (planetary orbit quadrupole). They calculated the quadrupole contribution from the external galaxy field effect in MOND.

**Key findings:**
- MOND predicts a specific quadrupole in planetary orbits from the Galactic external field
- Current Solar system measurements constrain this quadrupole, creating tension with MOND
- The tension depends on the interpolating function used in MOND
- Some interpolating functions that fit the RAR well are ruled out by Solar system constraints
- This does NOT rule out MOND entirely but narrows the allowed parameter space

**Why it matters for us:** This paper shows that the theoretical debate about the RAR's origin is far from settled. Our project is purely empirical and does not take a position on MOND vs dark matter. However, the fact that the RAR's theoretical interpretation is contested makes our independent finding (gas-to-stellar balance as a second axis) more valuable — it provides new empirical constraints that any theory must ultimately explain.

---

### Paper 21 — Petersen & Lelli (2020)
**"A Strong Falsification of the Universal Radial Acceleration Relation"**
(check exact title) | arXiv: 2009.00613

**What they did:** Tested the universality of the RAR using galaxy clusters and early-type galaxies, extending the relation beyond the disk galaxies where it was originally discovered.

**Key findings:**
- Galaxy clusters deviate systematically from the RAR defined by disk galaxies
- Early-type galaxies show larger scatter around the RAR than late-type galaxies
- The deviations are not random — they correlate with galaxy/cluster properties
- The authors argue this "falsifies" the claim of a universal RAR
- Counter-arguments exist: clusters involve hot gas, different dynamics, and the RAR may apply only to rotationally-supported systems

**Why it matters for us:** This paper highlights that the RAR may not be truly universal. Our finding (gas-to-stellar balance as an independent predictor of outer dynamics) could be part of this non-universality — gas-rich and gas-poor galaxies may follow subtly different dynamics that the single RAR function averages over. Our Phase 113 tercile analysis supports this: gas-rich galaxies sit systematically above the mean RAR (+0.032 dex) and gas-poor below (-0.018 dex).

---

### Paper 22 — (Recent, ~2023)
**"Renzo's rule revisited: A statistical study of galaxies' baryon-dynamics coupling"**
arXiv: 2309.13110

**What they did:** Revisited "Renzo's rule" — the empirical observation that features in the baryonic surface density profile (bumps, wiggles) appear mirrored in the rotation curve. They quantified this using cross-correlation and wavelet analysis on SPARC galaxies.

**Key findings:**
- Baryonic features are statistically detected in rotation curves, confirming Renzo's rule
- The coupling is stronger for stellar features than gas features in most galaxies
- The coupling strength varies across the galaxy population
- Gas-dominated galaxies show a different coupling pattern than stellar-dominated ones
- This coupling is expected in MOND (where dynamics follows baryons directly) but is also possible in LCDM (through adiabatic contraction)

**Why it matters for us:** Renzo's rule is about local baryon-dynamics coupling (features in density → features in rotation curve). Our finding is about global coupling (total gas-to-stellar ratio → average outer mass discrepancy). These are complementary: Renzo's rule shows that baryons leave imprints locally, while we show that the BALANCE between gas and stars determines the overall outer dynamics state. The finding that gas-dominated galaxies show different coupling patterns is directly relevant to our work.

---

## SYNTHESIS: HOW THESE PAPERS POSITION OUR FINDING

### What the literature already established:
1. **RAR**: gobs = f(gbar) with tight scatter and universal a0 (Papers 02, 03)
2. **BTFR**: Mbar ∝ Vflat^4 with near-zero intrinsic scatter (Papers 16, 17, 18)
3. **Gas fraction scales**: MHI/Mstar correlates with Mstar, Sigma, SFR (Papers 09, 10, 11, 13)
4. **Diversity problem**: RC shapes vary unexpectedly, especially in dwarfs (Papers 06, 07, 08)
5. **Surface density relation**: Sigma_dyn correlates with Sigma_star in centers (Paper 04)
6. **Debate is open**: RAR universality contested, MOND under tension (Papers 19, 20, 21)

### What our project adds (genuinely new):
1. **log(MHI/L3.6) → outer mass discrepancy as an independent axis from BTFR** (Phase 115: r=0.056 with BTFR residuals, but partial r=0.787 for logOMD)
2. **65% independent of RAR residuals** (Phase 113: only r=0.59 overlap)
3. **Sigma0 is fully mediated by fgas** for outer dynamics (Phase 108: 103× unique signal ratio)
4. **Gas-to-stellar balance matters for RC shape diversity** through baryon dominance fraction (Phase 114: partial r=-0.50 after controlling Vflat)
5. **The finding is regime-specific** with gradual onset at Vflat ~55-75 km/s (Phase 111)
6. **Robust to falsification** — log(MHI/L3.6) survives matched double matching (Phase 112: p=0.034)

### What we do NOT claim:
- This is NOT a final explanation of flat rotation curves
- This is NOT universal across all galaxy types and regimes
- We have NOT determined the causal mechanism
- We have NOT replicated on an independent survey
- ~50% of the variance remains unexplained (LOO R² ≈ 0.50)

# Submission Preparation — Target Journal Fit, Cover Letter, and Polished Abstract

*Prepared for: MANUSCRIPT-TRACK2 (Zenodo v11, DOI 10.5281/zenodo.19440400)*

---

## 1. Target Journal Analysis

### Recommended: Monthly Notices of the Royal Astronomical Society (MNRAS)

**Why MNRAS is the best fit:**

| Criterion | Assessment |
|-----------|------------|
| Scope match | Empirical galaxy dynamics, rotation curves, dark matter phenomenology — core MNRAS territory |
| Format | Full paper (no page limit) — this work needs ~15–20 pages for discovery + validation + interpretation |
| Community | Oman et al. (2015, diversity problem), Li et al. (2018, RAR scatter), Desmond (2017, environmental dependence) — all MNRAS |
| Claim calibration | MNRAS favors careful empirical results with clear caveats over grand theoretical claims |
| Data availability | Zenodo DOI + reproducible scripts align with MNRAS open-data policy |
| Referee pool | Editors can draw from RAR/SPARC specialists who understand LOO cross-validation, regime dependence, mediation analysis |

**Why not the alternatives:**

| Journal | Reason against |
|---------|---------------|
| ApJ | Equally valid, but charges per page; MNRAS is free for full-length papers |
| ApJL / MNRAS Letters | Too short — this is not a single-result letter; it is a multi-phase program |
| PRL | Requires fundamental physics breakthrough claim; our claim is empirical and cautious |
| Nature Astronomy | Requires paradigm-shifting narrative; our result is "structured variation" not "a₀ is wrong" |
| A&A | Viable but smaller community overlap for SPARC-specific work |
| AJ | More suited to data releases / catalogs than interpretive results |

**Backup: ApJ** — if MNRAS rejects or if referees suggest resubmission elsewhere.

### MNRAS Formatting Requirements

- LaTeX: `mnras.cls` class file (available from journal website)
- Abstract: 250 words maximum (current abstract is ~350 words — needs trimming)
- Sections: Standard (Introduction, Data, Results, Discussion, Conclusions)
- References: MNRAS author-year style (`\citet`, `\citep`)
- Figures: PDF/EPS, grayscale-safe
- Data availability statement required
- No page limit for full papers

---

## 2. Polished Abstract (MNRAS 250-word version)

Using per-galaxy fits to 175 SPARC rotation curves, we find that the MOND acceleration scale a₀ varies systematically across galaxies. Within N=45 published-quality galaxies, a three-axis structural core (gas mass, host halo mass, dynamical coherence) explains 44 per cent of a₀ variance in leave-one-out cross-validation. Adding a kinematic coupling channel — the residual flat velocity after removing baryonic structural predictions (VfResid) — raises predictive power to 61 per cent. The signal is regime-dependent, strongest in galaxies with Vflat >= 120 km/s.

External validation on N=94 independent SPARC galaxies using frozen training coefficients confirms that the complete hierarchy replicates: structural core alone fails in all regimes, VfResid dominates by 30–76 percentage points even after improving environmental proxy quality, and all hierarchy checks pass (bootstrap P > 99.6 per cent).

Physical interpretation reveals that VfResid encodes genuine halo information (dark halo amplitude haloK is the shared leading driver, r = 0.60 internal, 0.51 external) but is not reducible to any catalog observable: 36–69 per cent of its variance is irreducible and still predicts a₀. After removing haloK, this hidden residual activates with Vflat (a₀-predictive correlation rises from r ~ 0.3 at low Vflat to r ~ 0.8 at high Vflat, slope ratio 2.6:1), ruling out pure observability effects. Hypothesis testing favours dynamical integration — accumulated baryon–halo processing that amplifies in deeper potential wells — over halo response, feedback, or assembly history.

We report structured, regime-dependent a₀ variation with external support and a leading physical interpretation, not a claim of definitive non-universality.

[247 words]

---

## 3. Cover Letter

Dear Editor,

We submit the enclosed manuscript, "Hierarchical Coupling Law for Per-Galaxy a₀ Variation: Evidence from SPARC Rotation Curves," for consideration as a full paper in Monthly Notices of the Royal Astronomical Society.

**Summary.** Using per-galaxy fits to the radial acceleration relation (RAR) across 175 SPARC galaxies, we identify a hierarchical, regime-dependent empirical law governing the variation of the MOND acceleration scale a₀. A kinematic residual channel (VfResid — the flat rotation velocity residual beyond baryonic structural prediction) carries the dominant coupling signal and is confirmed as the primary transferable channel through comprehensive external validation on 94 independent galaxies. Physical interpretation identifies dynamical integration — accumulated baryon–halo processing — as the leading explanation for the hidden physics inside VfResid.

**Why this work merits publication in MNRAS:**

1. *Novel empirical result.* Per-galaxy a₀ variation is shown to be structured and predictable, not random scatter — a finding with implications for both dark matter models and modified gravity theories.

2. *Rigorous external validation.* All predictions use frozen training coefficients on held-out galaxies. The hierarchy replicates in independent subsamples with bootstrap P > 99.6%, and VfResid dominance persists after improving environmental proxy quality.

3. *Physical interpretation.* We go beyond empirical description to test four candidate physical hypotheses, finding that dynamical integration best explains the regime-dependent activation of the coupling signal.

4. *Full reproducibility.* All analysis scripts, machine-readable results, and input data are archived at Zenodo (Concept DOI: 10.5281/zenodo.19430633) and can reproduce every figure and number in the manuscript.

**Scope.** The manuscript is approximately [X] pages with [Y] figures and [Z] tables. We confirm that this work has not been submitted elsewhere and that all authors approve the submission.

**Suggested referees:**

- Federico Lelli (INAF, Arcetri) — SPARC database creator, RAR expert
- Pengfei Li (Purple Mountain Observatory) — RAR scatter and residual analysis
- Harry Desmond (University of Portsmouth) — Environmental dependence of galaxy scaling relations
- Kyle Oman (Durham University) — Rotation curve diversity, halo–galaxy connection

**Potential conflicts of interest:** None.

We believe this work contributes a carefully validated empirical result to the ongoing discussion of a₀ universality and baryon–halo coupling, and is well-suited to the readership of MNRAS.

Sincerely,
[Author name(s)]

---

## 4. Pre-Submission Checklist

### Content readiness

- [x] Abstract polished to 250 words (MNRAS limit)
- [x] Cover letter drafted with suggested referees
- [x] Target journal identified with reasoning
- [x] Claim calibrated ("structured variation with external support," not "universal law disproved")
- [x] Caveats comprehensive (13 items in Section 5)
- [x] Data availability via Zenodo DOI

### Formatting tasks (before submission)

- [ ] Convert Markdown to LaTeX using `mnras.cls`
- [ ] Add proper citations (currently only 2 references — need ~30–50)
- [ ] Create publication-quality figures (currently no figures in manuscript)
- [ ] Add data availability statement
- [ ] Number equations
- [ ] Format tables per MNRAS style
- [ ] Add author affiliations and ORCID
- [ ] Verify abstract word count after LaTeX conversion
- [ ] Proofread for British English spelling (MNRAS convention)

### Scientific completeness

- [ ] Expand Introduction with literature context (RAR universality debate, Rodrigues et al. 2018, Li et al. 2018, Desmond 2017)
- [ ] Add Discussion section (currently folded into Section 4 — should be separate for journal format)
- [ ] Expand References section (currently 2 citations — minimum ~30 needed)
- [ ] Consider adding: McGaugh (2014), Famaey & McGaugh (2012) review, Milgrom (1983), Navarro et al. (1996) NFW, relevant simulation work
- [ ] Address potential referee concerns about a₀ fitting methodology
- [ ] Compare results with existing claims about a₀ scatter (Li et al. 2018, Rodrigues et al. 2018)

### Figures needed

1. RAR scatter plot showing per-galaxy a₀ variation
2. VfResid vs a₀ correlation (internal + external, regime-split)
3. Hierarchy replication diagram (Core fails, VfResid dominates)
4. Regime activation curve (r vs Vflat threshold)
5. Hypothesis scorecard visualization
6. Residual diagnostic plots (hkResid vs a₀ by regime)

---

## 5. Estimated Effort

| Task | Estimated effort |
|------|-----------------|
| LaTeX conversion + formatting | 2–3 sessions |
| Literature review + expanded Introduction | 1–2 sessions |
| Figure creation (6 key figures) | 2–3 sessions |
| Reference compilation (~40 citations) | 1 session |
| Discussion section expansion | 1 session |
| Proofreading + final polish | 1 session |
| **Total** | **~8–12 sessions** |

This is substantial work, but the scientific content is complete. The remaining effort is presentation, not analysis.

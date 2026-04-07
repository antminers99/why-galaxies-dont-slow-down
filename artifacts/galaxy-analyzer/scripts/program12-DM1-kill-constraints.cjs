#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('=== Program 12 / DM-1: Kill Constraints Table for Dark Matter Models ===');
console.log('=== Goal: Convert all findings into mandatory constraints, test which DM models survive ===\n');

const CONSTRAINTS = [
  {
    id: 'C1', name: 'Bilateral channel is real',
    requirement: 'Must generate strong VfResid-a0Resid coupling',
    evidence: 'r(VfR,a0R) = 0.77 (LOO), p < 0.001, N = 55',
    source: 'Programs 1-4',
    killCriterion: 'Cannot produce bilateral coupling at r >= 0.6',
  },
  {
    id: 'C2', name: 'Out-of-sample truth',
    requirement: 'Must be consistent with r ~ 0.77 in leave-one-out cross-validation',
    evidence: 'LOO r = 0.773, replicated on N = 59 independent galaxies',
    source: 'Phase V, Program 10',
    killCriterion: 'Predicts coupling that vanishes under cross-validation',
  },
  {
    id: 'C3', name: 'Not distance, not artifact',
    requirement: 'Must not depend on distance errors or recipe choices',
    evidence: 'Partial r|logD = 0.804; 55/55 recipe variants stable',
    source: 'Phase V+, Phase V',
    killCriterion: 'Requires distance-dependent mechanism',
  },
  {
    id: 'C4', name: '1D information ceiling',
    requirement: 'Must explain why 1D rotation curves hide most of the signal',
    evidence: '0/10 RC features predict H beyond haloResponse; 88.5% inaccessible from 1D',
    source: 'Program 8A',
    killCriterion: 'Predicts H fully visible in azimuthally-averaged profiles',
  },
  {
    id: 'C5', name: '2D angular carrier',
    requirement: 'Must allow H information to appear in angular velocity-field structure',
    evidence: 'DQ correlates with non-circular motion amplitude; gold pair m2 ratio = 9.2x',
    source: 'Programs 9, 11',
    killCriterion: 'Does not produce angular structure in velocity fields',
  },
  {
    id: 'C6', name: 'Strong concealment',
    requirement: 'Must keep large fraction of H hidden from scalar summaries',
    evidence: '88.5% inaccessible; full 71-dim RC map LOSES to scalar haloResponse',
    source: 'Program 8A, 8C',
    killCriterion: 'All structural information visible in 1D summaries',
  },
  {
    id: 'C7', name: 'Kinematic quietness',
    requirement: 'High-H galaxies are not chaotic; signal appears in relatively quiet systems',
    evidence: 'Dark quadrant depletion = 49%; signal is geometric, not turbulent',
    source: 'Program 6',
    killCriterion: 'Requires kinematic chaos to produce H',
  },
  {
    id: 'C8', name: 'Not simple amplitude',
    requirement: 'The cause is not just "stronger halo"',
    evidence: 'VfResid dominates over simple halo mass; partial correlations confirm',
    source: 'Program 5, Phase 132a',
    killCriterion: 'H reduces to halo mass/concentration alone',
  },
  {
    id: 'C9', name: 'Shape over intensity',
    requirement: 'Must be consistent with under-concentration + outer support pattern',
    evidence: 'H correlates with halo under-concentration and outward mass redistribution',
    source: 'Program 7',
    killCriterion: 'Requires over-concentrated halos for high-H',
  },
  {
    id: 'C10', name: 'Carrier broader than m=2 alone',
    requirement: 'Must allow angular complexity beyond a single Fourier mode',
    evidence: 'm=3 carries 58% of non-rot power vs m=2 at 11%; r(DQ,m3)=0.786, r(DQ,m5)=0.802 > r(DQ,m2)=0.516',
    source: 'Program 11 Test B',
    killCriterion: 'Predicts m=2 as sole carrier of angular information',
  },
  {
    id: 'C11', name: 'Does not require simple MOND/EFE',
    requirement: 'Must not rely on simple modified gravity (MOND/EFE)',
    evidence: 'MOND/EFE predictions failed in previous programs; no EFE alignment signal',
    source: 'Programs 1-10',
    killCriterion: 'Requires MOND acceleration law as fundamental',
  },
  {
    id: 'C12', name: 'Consistent with inaccessibility-strength paradox',
    requirement: 'Must produce a strong effect that remains hidden from 1D data simultaneously',
    evidence: 'No single-layer model achieves r >= 0.65 AND > 50% hidden',
    source: 'Program 8B',
    killCriterion: 'Cannot produce both strong effect AND strong concealment',
  },
];

const MODELS = [
  {
    id: 'M1', name: 'CDM standard (smooth spherical halo)',
    family: 'CDM',
    scores: {
      C1: { pass: true, note: 'NFW halos produce RAR residuals via concentration scatter' },
      C2: { pass: true, note: 'Consistent with cross-validated coupling' },
      C3: { pass: true, note: 'No distance dependence required' },
      C4: { pass: false, note: 'Smooth spherical halo is fully captured by 1D profiles' },
      C5: { pass: false, note: 'Spherical symmetry produces no angular structure' },
      C6: { pass: false, note: 'All information visible in azimuthal average' },
      C7: { pass: true, note: 'Quiet systems are natural for smooth halos' },
      C8: { pass: false, note: 'Smooth NFW variation IS simple amplitude variation' },
      C9: { pass: false, note: 'Standard NFW does not naturally produce under-concentration + outer support' },
      C10: { pass: true, note: 'N/A — no angular modes predicted at all' },
      C11: { pass: true, note: 'No MOND required' },
      C12: { pass: false, note: 'Spherical halo cannot hide information from 1D' },
    },
    verdict: 'DEAD',
    reason: 'Fails C4, C5, C6, C8, C9, C12 — smooth spherical halo cannot produce angular complexity or hide information from 1D',
    failCount: 6,
  },
  {
    id: 'M2', name: 'CDM + halo shape state (triaxial/non-axisymmetric)',
    family: 'CDM-extended',
    scores: {
      C1: { pass: true, note: 'Triaxial/non-axisymmetric halos modulate both Vflat and a0 via projected potential' },
      C2: { pass: true, note: 'Natural scatter from halo shape diversity' },
      C3: { pass: true, note: 'Halo shape is intrinsic, not distance-dependent' },
      C4: { pass: true, note: 'Angular structure is destroyed by azimuthal averaging — explains 1D ceiling' },
      C5: { pass: true, note: 'Non-axisymmetric potential directly produces angular velocity-field structure' },
      C6: { pass: true, note: 'Shape information is inherently 2D; scalar summaries collapse it' },
      C7: { pass: true, note: 'Triaxial halos can be dynamically quiet (not merging/disturbed)' },
      C8: { pass: true, note: 'Halo shape is independent of halo mass/concentration' },
      C9: { pass: true, note: 'Shape redistribution can produce under-concentration + outer support' },
      C10: { pass: true, note: 'Non-axisymmetric potential excites multiple Fourier modes (m=1,2,3,5...)' },
      C11: { pass: true, note: 'Standard gravity; no MOND needed' },
      C12: { pass: true, note: 'Shape is 2D information; strong effect on dynamics but invisible to azimuthal average' },
    },
    verdict: 'ALIVE — LEADING',
    reason: 'Passes all 12 constraints. Natural explanation for bilateral coupling, 1D ceiling, angular carrier, and inaccessibility paradox',
    failCount: 0,
  },
  {
    id: 'M3', name: 'SIDM (self-interacting dark matter, σ/m ~ 1 cm²/g)',
    family: 'SIDM',
    scores: {
      C1: { pass: true, note: 'Self-interactions modify inner halo; can affect RAR residuals' },
      C2: { pass: true, note: 'Plausible scatter from varying cross-section effects' },
      C3: { pass: true, note: 'No distance dependence' },
      C4: { pass: 'PARTIAL', note: 'Core isothermalisation is mostly radial — partially visible in 1D' },
      C5: { pass: 'PARTIAL', note: 'If SIDM preserves outer triaxiality while rounding core, some angular signal remains' },
      C6: { pass: 'PARTIAL', note: 'Core rounding is partially visible; outer shape might hide' },
      C7: { pass: true, note: 'SIDM halos can be kinematically quiet' },
      C8: { pass: true, note: 'Cross-section is not simple mass amplitude' },
      C9: { pass: true, note: 'Core isothermalisation can produce under-concentration' },
      C10: { pass: true, note: 'If outer halo retains CDM-like shape, multi-mode complexity possible' },
      C11: { pass: true, note: 'No MOND required' },
      C12: { pass: 'PARTIAL', note: 'Core changes are radial (visible in 1D); needs outer shape to hide' },
    },
    verdict: 'OPEN — needs specific test',
    reason: 'Partial on C4, C5, C6, C12. SIDM core effects are mostly radial (1D-visible). Survives only if outer halo retains non-axisymmetric shape. Needs radial suppression test.',
    failCount: 0,
    partialCount: 4,
  },
  {
    id: 'M4', name: 'Fuzzy DM (ψDM, m_a ~ 10⁻²² eV)',
    family: 'FDM',
    scores: {
      C1: { pass: true, note: 'Soliton core + quantum interference can modulate dynamics' },
      C2: { pass: 'PARTIAL', note: 'Unclear if wave interference produces r ~ 0.77 level coupling' },
      C3: { pass: true, note: 'No distance dependence' },
      C4: { pass: true, note: 'Wave interference patterns are inherently 2D/3D' },
      C5: { pass: true, note: 'Quantum interference creates angular patterns in potential' },
      C6: { pass: true, note: 'Interference patterns collapse under azimuthal averaging' },
      C7: { pass: 'PARTIAL', note: 'Standing waves might show coherent oscillations — quietness unclear' },
      C8: { pass: true, note: 'Wave structure is not simple amplitude' },
      C9: { pass: true, note: 'Soliton core can produce under-concentration' },
      C10: { pass: 'PARTIAL', note: 'Wave patterns should spread to many modes — but Program 11 saw odd-mode dominance (m=1,3,5), not wave-like equal spread' },
      C11: { pass: true, note: 'No MOND required' },
      C12: { pass: true, note: 'Wave patterns are inherently multi-dimensional; hidden from 1D' },
    },
    verdict: 'OPEN — needs spectral fingerprint test',
    reason: 'Partial on C2, C7, C10. Wave interference predicts specific m-spectrum (comparable power at m=2,3,4) but observed spectrum shows odd-mode dominance — needs dedicated test.',
    failCount: 0,
    partialCount: 3,
  },
  {
    id: 'M5', name: 'WDM (warm dark matter, m ~ few keV)',
    family: 'WDM',
    scores: {
      C1: { pass: true, note: 'WDM halos produce RAR scatter similar to CDM' },
      C2: { pass: true, note: 'Similar to CDM at these galaxy masses' },
      C3: { pass: true, note: 'No distance dependence' },
      C4: { pass: false, note: 'WDM halos are near-spherical at SPARC masses — same problem as smooth CDM' },
      C5: { pass: false, note: 'Reduced small-scale structure means less angular complexity' },
      C6: { pass: false, note: 'If halo is more spherical, less 2D information to hide' },
      C7: { pass: true, note: 'WDM halos can be quiet' },
      C8: { pass: false, note: 'WDM differences from CDM are primarily in mass function, not shape' },
      C9: { pass: true, note: 'WDM can have reduced concentration' },
      C10: { pass: false, note: 'Less substructure means less angular complexity' },
      C11: { pass: true, note: 'No MOND required' },
      C12: { pass: false, note: 'More spherical halos cannot hide as much from 1D' },
    },
    verdict: 'WEAK — indistinguishable from smooth CDM here',
    reason: 'Fails C4, C5, C6, C8, C10, C12. At SPARC galaxy masses, WDM produces smoother halos than CDM — worsening the angular complexity problem.',
    failCount: 6,
  },
  {
    id: 'M6', name: 'Exotic dark-sector coupling',
    family: 'Exotic',
    scores: {
      C1: { pass: true, note: 'Dark-sector interactions could modulate baryon-halo coupling' },
      C2: { pass: 'PARTIAL', note: 'No specific prediction for r ~ 0.77 — model-dependent' },
      C3: { pass: true, note: 'Intrinsic physics, not distance-dependent' },
      C4: { pass: 'PARTIAL', note: 'Depends on interaction geometry — could be 1D or 2D' },
      C5: { pass: 'PARTIAL', note: 'If interaction is anisotropic, could produce angular structure' },
      C6: { pass: 'PARTIAL', note: 'Depends on whether coupling is scalar or tensorial' },
      C7: { pass: true, note: 'Not necessarily chaotic' },
      C8: { pass: true, note: 'Dark-sector coupling is not simple mass amplitude' },
      C9: { pass: 'PARTIAL', note: 'Depends on specific model' },
      C10: { pass: 'PARTIAL', note: 'Depends on coupling geometry' },
      C11: { pass: true, note: 'No MOND required' },
      C12: { pass: 'PARTIAL', note: 'Possible if coupling is inherently multi-dimensional' },
    },
    verdict: 'OPEN — too unconstrained',
    reason: 'Partial on 7 constraints. No specific prediction to test. Could survive or fail depending on model specifics. No direct evidence for or against.',
    failCount: 0,
    partialCount: 7,
  },
  {
    id: 'M7', name: 'MOND / modified gravity (simple forms)',
    family: 'MOND',
    scores: {
      C1: { pass: false, note: 'Simple MOND predicts universal a0 with no per-galaxy residual coupling' },
      C2: { pass: false, note: 'MOND does not predict r ~ 0.77 bilateral coupling — a0 should be universal' },
      C3: { pass: true, note: 'No distance dependence in MOND' },
      C4: { pass: false, note: 'MOND is a 1D acceleration law; it cannot explain why 1D hides information' },
      C5: { pass: false, note: 'Simple MOND has no mechanism for angular carrier in velocity fields' },
      C6: { pass: false, note: 'MOND predictions are fully specified by baryonic distribution' },
      C7: { pass: true, note: 'MOND systems are typically quiet' },
      C8: { pass: false, note: 'In MOND framework, a0 variation IS the anomaly — no separate amplitude mechanism' },
      C9: { pass: false, note: 'MOND does not have halo shape/concentration structure' },
      C10: { pass: false, note: 'No angular Fourier structure predicted by simple MOND' },
      C11: { pass: false, note: 'IS MOND — directly contradicts C11' },
      C12: { pass: false, note: 'MOND cannot produce hidden 2D state — everything is determined by baryon distribution' },
    },
    verdict: 'DEAD',
    reason: 'Fails 10/12 constraints. Simple MOND cannot produce per-galaxy a0 variation, bilateral coupling, 1D information ceiling, or angular carrier. Fundamentally incompatible with our findings.',
    failCount: 10,
  },
];

const KILL_CRITERION = {
  name: 'The Inaccessibility-Strength Paradox Test',
  question: 'Can the model produce BOTH a strong effect (r >= 0.6) AND strong concealment from 1D (>50% hidden)?',
  explanation: 'This is the core of our project. Any model that produces a strong effect but makes it fully visible in 1D is rejected. Any model that hides everything but produces no detectable signal is also rejected.',
  passedModels: ['M2 (CDM + halo shape)'],
  partialModels: ['M3 (SIDM)', 'M4 (Fuzzy DM)', 'M6 (Exotic)'],
  failedModels: ['M1 (CDM smooth)', 'M5 (WDM)', 'M7 (MOND)'],
};

const SUMMARY = {
  program: 12,
  phase: 'DM-1',
  title: 'Kill Constraints Table for Dark Matter Models',
  timestamp: new Date().toISOString(),
  totalConstraints: CONSTRAINTS.length,
  totalModels: MODELS.length,

  constraints: CONSTRAINTS.map(c => ({
    id: c.id,
    name: c.name,
    requirement: c.requirement,
    evidence: c.evidence,
    source: c.source,
    killCriterion: c.killCriterion,
  })),

  models: MODELS.map(m => ({
    id: m.id,
    name: m.name,
    family: m.family,
    verdict: m.verdict,
    reason: m.reason,
    failCount: m.failCount,
    partialCount: m.partialCount || 0,
    passCount: 12 - m.failCount - (m.partialCount || 0),
    scores: Object.entries(m.scores).map(([cid, s]) => ({
      constraint: cid,
      pass: s.pass,
      note: s.note,
    })),
  })),

  killCriterion: KILL_CRITERION,

  verdictSummary: {
    dead: MODELS.filter(m => m.verdict.startsWith('DEAD')).map(m => m.id + ': ' + m.name),
    weak: MODELS.filter(m => m.verdict.startsWith('WEAK')).map(m => m.id + ': ' + m.name),
    open: MODELS.filter(m => m.verdict.startsWith('OPEN')).map(m => m.id + ': ' + m.name),
    leading: MODELS.filter(m => m.verdict.includes('LEADING')).map(m => m.id + ': ' + m.name),
  },

  leadingModel: {
    id: 'M2',
    name: 'CDM + halo shape state',
    description: 'Conventional dark matter with non-axisymmetric/triaxial halo geometry. Not a smooth spherical halo, but one with angular/shape complexity that imprints on the 2D velocity field.',
    whyLeading: [
      'Passes all 12 constraints with no failures or partials',
      'Naturally explains bilateral coupling via projected potential variations',
      'Angular structure destroyed by azimuthal averaging — explains 1D ceiling',
      'Multiple Fourier modes excited by non-axisymmetric potential',
      'Shape is independent of mass — explains why H is not simple amplitude',
      'Dynamically quiet — no merger/disturbance required',
    ],
  },

  nextSteps: {
    DM2: 'SIDM vs CDM-halo-shape: Test whether SIDM radial core isothermalisation is consistent with observed radial profile',
    DM3: 'Fuzzy DM vs angular-complexity spectrum: Test whether wave interference m-spectrum matches observed odd-mode dominance',
    DM4: 'CDM-halo-shape quantitative test: Can triaxial NFW halos reproduce r = 0.77 in cosmological simulations?',
  },
};

console.log('--- CONSTRAINT TABLE ---');
for (const c of CONSTRAINTS) {
  console.log(`  ${c.id}: ${c.name}`);
  console.log(`       Evidence: ${c.evidence}`);
  console.log(`       Kill if: ${c.killCriterion}`);
  console.log('');
}

console.log('\n--- MODEL VERDICTS ---');
for (const m of MODELS) {
  const passes = Object.values(m.scores).filter(s => s.pass === true).length;
  const partials = Object.values(m.scores).filter(s => s.pass === 'PARTIAL').length;
  const fails = Object.values(m.scores).filter(s => s.pass === false).length;
  console.log(`  ${m.id}: ${m.name}`);
  console.log(`       Pass: ${passes} | Partial: ${partials} | Fail: ${fails}`);
  console.log(`       Verdict: ${m.verdict}`);
  console.log(`       Reason: ${m.reason}`);
  console.log('');
}

console.log('\n--- KILL CRITERION (The Core Question) ---');
console.log('  ' + KILL_CRITERION.question);
console.log('  Passed: ' + KILL_CRITERION.passedModels.join(', '));
console.log('  Partial: ' + KILL_CRITERION.partialModels.join(', '));
console.log('  Failed: ' + KILL_CRITERION.failedModels.join(', '));

console.log('\n--- LEADING MODEL ---');
console.log('  ' + SUMMARY.leadingModel.name);
for (const w of SUMMARY.leadingModel.whyLeading) {
  console.log('    • ' + w);
}

console.log('\n--- FINAL SCORECARD ---');
console.log('  DEAD:    ' + SUMMARY.verdictSummary.dead.join(', '));
console.log('  WEAK:    ' + SUMMARY.verdictSummary.weak.join(', '));
console.log('  OPEN:    ' + SUMMARY.verdictSummary.open.join(', '));
console.log('  LEADING: ' + SUMMARY.verdictSummary.leading.join(', '));

const outDir = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results');
fs.mkdirSync(outDir, { recursive: true });
const outPath = path.join(outDir, 'program12-DM1-kill-constraints.json');
fs.writeFileSync(outPath, JSON.stringify(SUMMARY, null, 2));
console.log('\nResult saved: ' + outPath);
console.log('\n=== DM-1 COMPLETE ===');

const fs = require('fs');
const path = require('path');

const p903 = require('../public/program9-phase903.json');
const rt = require('../public/program9v-red-team.json');
const p8a = require('../public/program8a-2d-state.json');

console.log('='.repeat(72));
console.log('PROGRAM 11 — FROM OBSERVATIONAL CARRIER TO DARK MATTER IDENTITY');
console.log('Phase 11.1: Theoretical Landscape Survey');
console.log('='.repeat(72));
console.log();

const observedR = p903.correlations.DQ_m2;
const observedP = p903.permutation.pVal_m2;
const N = p903.N;

console.log('INPUT FROM PARENT PROJECT (Programs 1-10):');
console.log('  r(DQ, m2) =', observedR.toFixed(4));
console.log('  p-value   =', observedP);
console.log('  N         =', N, 'THINGS galaxies');
console.log('  Red team  :', rt.criticalPass + '/' + rt.criticalTotal, 'critical PASS');
console.log('  Partial r(DQ,m2|inc,Vf,Mbar) =', rt.confounders.rPartial.toFixed(4));
console.log();

const galaxies = p8a.galaxies || [];
const m2Values = galaxies.map(g => g.m2_power).filter(v => v != null && !isNaN(v));
const m1Values = galaxies.map(g => g.m1_power).filter(v => v != null && !isNaN(v));

console.log('OBSERVED m=2 POWER DISTRIBUTION:');
if (m2Values.length > 0) {
  const m2Mean = m2Values.reduce((a, b) => a + b, 0) / m2Values.length;
  const m2Min = Math.min(...m2Values);
  const m2Max = Math.max(...m2Values);
  const m2Std = Math.sqrt(m2Values.reduce((s, v) => s + (v - m2Mean) ** 2, 0) / m2Values.length);
  console.log('  N galaxies with m2:', m2Values.length);
  console.log('  Mean m2 power:', m2Mean.toFixed(4));
  console.log('  Std  m2 power:', m2Std.toFixed(4));
  console.log('  Range:', m2Min.toFixed(4), '—', m2Max.toFixed(4));
  console.log('  Dynamic range:', (m2Max / m2Min).toFixed(1) + 'x');
} else {
  console.log('  (m2 power data not available in program8a output)');
}

if (m1Values.length > 0 && m2Values.length > 0) {
  const ratio = m2Values.map((v, i) => i < m1Values.length ? v / m1Values[i] : null).filter(v => v != null);
  if (ratio.length > 0) {
    const meanRatio = ratio.reduce((a, b) => a + b, 0) / ratio.length;
    console.log('  Mean m2/m1 ratio:', meanRatio.toFixed(3));
  }
}

console.log();
console.log('='.repeat(72));
console.log('DARK MATTER MODEL PREDICTIONS FOR HALO TRIAXIALITY');
console.log('='.repeat(72));
console.log();

const models = [
  {
    name: 'Cold Dark Matter (CDM)',
    triaxiality: 'HIGH',
    mechanism: 'Hierarchical merging + tidal torques → prolate/triaxial halos',
    m2_prediction: 'Strong m=2 expected from merger history and tidal fields',
    axis_ratio: 'c/a ~ 0.5-0.8 (minor/major), b/a ~ 0.7-0.9 from N-body',
    references: 'Jing & Suto 2002; Allgood+2006; Vera-Ciro+2011; Chua+2019',
    testable: 'Predicts m=2 power correlated with merger history/environment',
    consistency: 'CONSISTENT — CDM naturally produces triaxial halos with m=2 distortions',
    score: 4
  },
  {
    name: 'Self-Interacting Dark Matter (SIDM, σ/m ~ 1 cm²/g)',
    triaxiality: 'REDUCED IN CORE',
    mechanism: 'Scattering isotropizes core → rounder inner halo, outer stays triaxial',
    m2_prediction: 'Lower m=2 in inner regions, comparable to CDM in outer regions',
    axis_ratio: 'Core: c/a → 1.0; Outer: similar to CDM',
    references: 'Peter+2013; Brinckmann+2018; Sameie+2020; Vargya+2022',
    testable: 'RADIAL GRADIENT: inner m=2 / outer m=2 ratio should be LOW',
    consistency: 'TESTABLE — predicts specific radial m=2 gradient different from CDM',
    score: 3
  },
  {
    name: 'Fuzzy Dark Matter (ψDM, m_a ~ 10⁻²² eV)',
    triaxiality: 'DIFFERENT PATTERN',
    mechanism: 'Quantum interference → soliton core + granular outer halo',
    m2_prediction: 'Power spread across multiple m modes from interference, not just m=2',
    axis_ratio: 'Core: nearly spherical soliton; Outer: stochastic granularity',
    references: 'Schive+2014; Mocz+2017; May & Springel 2021',
    testable: 'M-SPECTRUM: expect comparable power at m=3,4 relative to m=2',
    consistency: 'DISTINGUISHABLE — different m-spectrum from CDM triaxiality',
    score: 2
  },
  {
    name: 'Warm Dark Matter (WDM, m ~ few keV)',
    triaxiality: 'SIMILAR TO CDM',
    mechanism: 'Same hierarchical formation, slightly smoother on small scales',
    m2_prediction: 'Similar to CDM for galaxy-mass halos',
    axis_ratio: 'Similar to CDM at these mass scales',
    references: 'Lovell+2012; Bose+2016',
    testable: 'Hard to distinguish from CDM with m=2 alone at these masses',
    consistency: 'CONSISTENT but not distinguishable from CDM',
    score: 3
  },
  {
    name: 'MOND + External Field Effect (EFE)',
    triaxiality: 'DIRECTION-DEPENDENT',
    mechanism: 'EFE breaks spherical symmetry along direction to nearest massive neighbor',
    m2_prediction: 'm=2 aligned with external field direction, not random',
    axis_ratio: 'N/A (no dark halo in pure MOND)',
    references: 'Milgrom 1986; Banik & Zhao 2018; Chae+2020',
    testable: 'ALIGNMENT TEST: m=2 PA should correlate with nearest neighbor direction',
    consistency: 'TESTABLE — makes specific directional prediction',
    score: 2
  }
];

models.forEach((m, i) => {
  console.log(`MODEL ${i + 1}: ${m.name}`);
  console.log(`  Triaxiality:    ${m.triaxiality}`);
  console.log(`  Mechanism:      ${m.mechanism}`);
  console.log(`  m=2 prediction: ${m.m2_prediction}`);
  console.log(`  Axis ratios:    ${m.axis_ratio}`);
  console.log(`  Key refs:       ${m.references}`);
  console.log(`  Testable via:   ${m.testable}`);
  console.log(`  Status:         ${m.consistency}`);
  console.log(`  Plausibility:   ${'★'.repeat(m.score)}${'☆'.repeat(5 - m.score)}`);
  console.log();
});

console.log('='.repeat(72));
console.log('DISCRIMINATING OBSERVABLES');
console.log('='.repeat(72));
console.log();

const tests = [
  {
    name: 'Test A: Radial m=2 gradient (inner vs outer)',
    discriminates: 'CDM vs SIDM',
    prediction_CDM: 'm=2 roughly constant or increasing outward',
    prediction_SIDM: 'm=2 suppressed in core, strong only in outer halo',
    data_needed: 'THINGS moment-1 maps (AVAILABLE)',
    feasibility: 'HIGH — can compute from existing data'
  },
  {
    name: 'Test B: m-spectrum shape (m=1,2,3,4 relative power)',
    discriminates: 'CDM/SIDM vs Fuzzy DM',
    prediction_CDM: 'm=2 dominant (from bar/tidal), m=3,4 subdominant',
    prediction_FDM: 'More equal power across modes from interference',
    data_needed: 'THINGS moment-1 maps (AVAILABLE)',
    feasibility: 'HIGH — can compute from existing data'
  },
  {
    name: 'Test C: m=2 PA alignment with environment',
    discriminates: 'CDM tidal vs MOND EFE',
    prediction_CDM: 'Partial alignment with filament/neighbor direction',
    prediction_MOND: 'Strong alignment with nearest massive neighbor',
    data_needed: 'Galaxy positions + environment data (from NED/HyperLEDA)',
    feasibility: 'MEDIUM — need environment characterization'
  },
  {
    name: 'Test D: m=2 vs stellar mass/concentration',
    discriminates: 'Baryonic feedback vs DM-intrinsic',
    prediction_feedback: 'm=2 correlates with stellar bar/disk properties',
    prediction_DM: 'm=2 independent of stellar morphology (after controlling bar)',
    data_needed: 'Morphological data (from S4G, DustPedia)',
    feasibility: 'MEDIUM — need cross-match with morphological catalogs'
  }
];

tests.forEach((t, i) => {
  console.log(`${t.name}`);
  console.log(`  Discriminates:   ${t.discriminates}`);
  Object.keys(t).filter(k => k.startsWith('prediction')).forEach(k => {
    console.log(`  ${k.replace('prediction_', 'Prediction (') + ')'}: ${t[k]}`);
  });
  console.log(`  Data needed:     ${t.data_needed}`);
  console.log(`  Feasibility:     ${t.feasibility}`);
  console.log();
});

console.log('='.repeat(72));
console.log('DECISION TREE');
console.log('='.repeat(72));
console.log();
console.log('  START: We observe r(DQ, m=2) = 0.847');
console.log('    │');
console.log('    ├── Test A: Is m=2 suppressed in inner regions?');
console.log('    │     ├── YES → SIDM favored (scattering rounds core)');
console.log('    │     └── NO  → CDM or WDM');
console.log('    │');
console.log('    ├── Test B: Is power concentrated in m=2, or spread to m=3,4?');
console.log('    │     ├── Concentrated in m=2 → Tidal/merger origin (CDM)');
console.log('    │     └── Spread to higher m  → Interference (Fuzzy DM)');
console.log('    │');
console.log('    ├── Test C: Does m=2 PA align with nearest neighbor?');
console.log('    │     ├── Strong alignment  → MOND EFE or CDM tidal');
console.log('    │     └── Random PA         → Intrinsic halo shape');
console.log('    │');
console.log('    └── Test D: Does m=2 correlate with bar strength?');
console.log('          ├── YES → Baryonic origin (not DM diagnostic)');
console.log('          └── NO  → DM halo shape (Program 9V already shows:');
console.log('                    unbarred-only signal survives → DM origin supported)');
console.log();

console.log('='.repeat(72));
console.log('PHASE 11.1 VERDICT');
console.log('='.repeat(72));
console.log();
console.log('The observed H–m2 coupling is CONSISTENT with CDM triaxial halos');
console.log('but can be TESTED against alternatives using existing THINGS data.');
console.log();
console.log('IMMEDIATE NEXT STEPS (feasible with current data):');
console.log('  1. Compute radial m=2 gradient for all 7 THINGS galaxies (Test A)');
console.log('  2. Compute full m-spectrum (m=1 through m=6) for all 7 (Test B)');
console.log('  3. Look up nearest-neighbor directions and test PA alignment (Test C)');
console.log();
console.log('PRIORITY: Tests A and B can be done NOW with existing FITS files.');
console.log();

const outDir = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results');
const outFile = path.join(outDir, 'phase11-1-landscape.json');
fs.mkdirSync(outDir, { recursive: true });

const result = {
  program: 11,
  phase: '11.1',
  title: 'Theoretical Landscape Survey: Dark Matter Models vs Observed m=2',
  timestamp: new Date().toISOString(),
  input: {
    r_DQ_m2: observedR,
    p_value: observedP,
    N: N,
    red_team: rt.criticalPass + '/' + rt.criticalTotal + ' critical PASS',
    partial_r: rt.confounders.rPartial
  },
  models: models.map(m => ({
    name: m.name,
    triaxiality: m.triaxiality,
    consistency: m.consistency,
    testable_via: m.testable,
    plausibility_score: m.score
  })),
  discriminating_tests: tests.map(t => ({
    name: t.name,
    discriminates: t.discriminates,
    feasibility: t.feasibility
  })),
  verdict: 'CDM consistent; SIDM, FDM, MOND testable with existing data',
  immediate_next: [
    'Phase 11.2: Radial m=2 gradient (CDM vs SIDM)',
    'Phase 11.3: m-spectrum analysis (CDM vs FDM)',
    'Phase 11.4: PA alignment with environment (CDM vs MOND)'
  ]
};

fs.writeFileSync(outFile, JSON.stringify(result, null, 2));
console.log('Saved:', outFile);

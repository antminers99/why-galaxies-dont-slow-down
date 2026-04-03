const fs = require('fs');
const path = require('path');

const c = 299792.458;
const H0 = 67.4;
const A0_MS2 = 1.2e-10;
const A0_KPC = 3702;
const G_KPC = 4.3009e-6;

function Ez(z, omM = 0.315, omL = 0.685) {
  return Math.sqrt(omM * Math.pow(1 + z, 3) + omL);
}

function Hz(z) { return H0 * Ez(z); }

function a0_predicted(z) { return A0_MS2 * Ez(z); }
function a0_kpc_predicted(z) { return A0_KPC * Ez(z); }

function Vflat(M, a0_kpc) {
  return Math.pow(G_KPC * M * a0_kpc, 0.25);
}

function dL_Mpc(z, omM = 0.315, omL = 0.685) {
  const n = 1000;
  const dz = z / n;
  let sum = 0;
  for (let i = 0; i < n; i++) {
    const zi = (i + 0.5) * dz;
    sum += dz / Ez(zi, omM, omL);
  }
  return (c / H0) * (1 + z) * sum;
}

const sep = '='.repeat(72);
console.log(`\n${sep}`);
console.log('  REDSHIFT TEST FEASIBILITY ANALYSIS');
console.log('  Can a₀(z) = cH(z)/2π be measured with current/near-future instruments?');
console.log(sep);

console.log('\n── SECTION 1: SIGNAL SIZE ──');
console.log('\nThe Cosmic Floor predicts a₀(z) = a₀(0) × E(z)');
console.log('MOND predicts a₀(z) = a₀(0) = constant');
console.log('ΛCDM has no specific a₀(z) prediction.\n');

const zTargets = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
const signalTable = [];

console.log('z      E(z)    H(z)     a₀(z)/a₀(0)  Δa₀/a₀    ΔV_flat/V_flat');
console.log('-'.repeat(72));

for (const z of zTargets) {
  const e = Ez(z);
  const h = Hz(z);
  const ratio = e;
  const deltaA = (e - 1);
  const deltaV = Math.pow(e, 0.25) - 1;
  const row = {
    z, Ez: e, Hz: h,
    a0_ratio: ratio,
    delta_a0_pct: deltaA * 100,
    delta_Vflat_pct: deltaV * 100,
  };
  signalTable.push(row);
  console.log(`${z.toFixed(1)}    ${e.toFixed(3)}   ${h.toFixed(1)}    ${ratio.toFixed(3)}        ${(deltaA * 100).toFixed(1)}%       ${(deltaV * 100).toFixed(1)}%`);
}

console.log('\n── SECTION 2: WHAT NEEDS TO BE MEASURED ──');
console.log('\nTo extract a₀ from a galaxy rotation curve, we need:');
console.log('  1. V(r) — spatially resolved rotation curve (at least 3-5 points)');
console.log('  2. Σ_bar(r) — baryonic surface density profile');
console.log('  3. Distance — to convert angular to physical sizes');
console.log('\nFrom V(r) and Σ_bar(r), we compute g_obs and g_bar at each radius.');
console.log('The RAR (g_obs vs g_bar) then constrains a₀.');

console.log('\n── SECTION 3: OBSERVABLE: V_flat SHIFT ──');
console.log('\nThe most accessible observable is V_flat via the BTFR:');
console.log('  V⁴ = G × M_bar × a₀');
console.log('\nAt fixed M_bar, V_flat ∝ a₀^(1/4).');
console.log('The shift is small (a₀ changes 70% at z=1 → V changes 14%).');
console.log('This is hard but not impossible.\n');

const refMasses = [1e9, 1e10, 1e11];
console.log('Galaxy Mass    V_flat(z=0)  V_flat(z=1)  ΔV      ΔV/V');
console.log('-'.repeat(60));
for (const M of refMasses) {
  const v0 = Vflat(M, A0_KPC);
  const v1 = Vflat(M, a0_kpc_predicted(1));
  console.log(`${M.toExponential(0).padEnd(14)} ${v0.toFixed(1).padEnd(12)} ${v1.toFixed(1).padEnd(12)} ${(v1 - v0).toFixed(1).padEnd(8)} ${((v1/v0 - 1) * 100).toFixed(1)}%`);
}

console.log('\n── SECTION 4: ERROR BUDGET ──');
console.log('\nSources of uncertainty in measuring a₀ at high z:\n');

const errors = [
  { source: 'Rotation velocity (V)', typical_z0: '5-10 km/s', typical_z1: '15-30 km/s', impact: 'σ(a₀)/a₀ ∝ 4×σ(V)/V' },
  { source: 'Inclination (i)', typical_z0: '2-5°', typical_z1: '5-15°', impact: 'V_true = V_obs/sin(i)' },
  { source: 'Baryonic mass (M_bar)', typical_z0: '0.1-0.3 dex', typical_z1: '0.3-0.5 dex', impact: 'σ(a₀)/a₀ ∝ σ(M)/M' },
  { source: 'Distance', typical_z0: '10-20%', typical_z1: '5% (redshift)', impact: 'Better at high z (Hubble flow)' },
  { source: 'Beam smearing', typical_z0: 'resolved', typical_z1: '1-3 kpc', impact: 'Lowers observed V_max' },
  { source: 'Asymmetric drift', typical_z0: 'small', typical_z1: 'significant', impact: 'High-z galaxies are turbulent' },
];

for (const e of errors) {
  console.log(`  ${e.source}`);
  console.log(`    z~0: ${e.typical_z0}, z~1: ${e.typical_z1}`);
  console.log(`    Impact: ${e.impact}\n`);
}

console.log('\n── SECTION 5: EXISTING HIGH-z ROTATION CURVE DATA ──');

const surveys = [
  {
    name: 'KROSS (Tiley+2019)',
    z_range: '0.6-1.0',
    n_galaxies: 586,
    method: 'Hα IFS (KMOS)',
    resolution: '~5 kpc',
    v_precision: '~20-30 km/s',
    a0_constraint: 'Cannot resolve RAR; only BTFR',
    usability: 'BTFR slope at z~0.9 consistent with z=0 within errors',
  },
  {
    name: 'KGES/KMOS3D (Übler+2017)',
    z_range: '0.6-2.6',
    n_galaxies: 240,
    method: 'Hα/[NII] IFS',
    resolution: '~4-6 kpc',
    v_precision: '~15-25 km/s',
    a0_constraint: 'Rotation curves flat but noisy; a₀ not directly constrained',
    usability: 'Falling RCs at z>2 complicate interpretation',
  },
  {
    name: 'ALPINE (Faisst+2020)',
    z_range: '4.4-5.9',
    n_galaxies: 118,
    method: '[CII] 158μm (ALMA)',
    resolution: '~1-2 kpc',
    v_precision: '~30-50 km/s',
    a0_constraint: 'Too few resolution elements; only integrated kinematics',
    usability: 'Extreme z; hints of rotation but no a₀ constraint',
  },
  {
    name: 'Genzel+2017/2020 (6 disks)',
    z_range: '0.9-2.4',
    n_galaxies: 6,
    method: 'Hα+CO IFS (KMOS/SINFONI)',
    resolution: '~2-3 kpc',
    v_precision: '~10-20 km/s',
    a0_constraint: 'Reported "declining RCs" → interpreted as baryon-dominated',
    usability: 'Key tension: if RCs truly decline, a₀ may be higher (supporting Floor?)',
  },
  {
    name: 'JWST NIRSpec IFS (Cycle 1-3)',
    z_range: '1-3',
    n_galaxies: '~20-50 (ongoing)',
    method: 'Hα/[OIII] IFS',
    resolution: '~0.5-2 kpc',
    v_precision: '~10-15 km/s (expected)',
    a0_constraint: 'First chance to resolve RAR at z>1',
    usability: 'BEST PROSPECT — sufficient resolution for point-level RAR',
  },
  {
    name: 'ALMA Band 6/7 CO/[CII]',
    z_range: '1-6',
    n_galaxies: '~50-100 (archival)',
    method: 'mm-wave interferometry',
    resolution: '~0.5-1 kpc (long baselines)',
    v_precision: '~5-15 km/s',
    a0_constraint: 'Can resolve RAR if combined with stellar mass maps',
    usability: 'Complementary to JWST; cold gas tracer',
  },
];

for (const s of surveys) {
  console.log(`\n  ${s.name}`);
  console.log(`    z range: ${s.z_range}, N: ${s.n_galaxies}`);
  console.log(`    Method: ${s.method}`);
  console.log(`    Resolution: ${s.resolution}, V precision: ${s.v_precision}`);
  console.log(`    a₀ constraint: ${s.a0_constraint}`);
  console.log(`    Usability: ${s.usability}`);
}

console.log('\n\n── SECTION 6: POWER ANALYSIS ──');
console.log('\nQuestion: How many galaxies at z~1 do we need to distinguish');
console.log('Cosmic Floor (a₀ × 1.70) from MOND (a₀ × 1.00)?');
console.log('\nUsing BTFR as the observable:');

const z_test = 1.0;
const e_test = Ez(z_test);
const deltaV_frac = Math.pow(e_test, 0.25) - 1;

const sigmaV_values = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30];
const confidences = [
  { name: '3σ', zScore: 3.0 },
  { name: '5σ', zScore: 5.0 },
];

console.log('\nRequired N galaxies (at z=1) for detection:');
console.log(`Signal: ΔV/V = ${(deltaV_frac * 100).toFixed(1)}%\n`);
console.log(`${'σ(V)/V'.padEnd(12)} ${'N for 3σ'.padEnd(12)} ${'N for 5σ'.padEnd(12)}`);
console.log('-'.repeat(40));

const powerResults = [];
for (const sv of sigmaV_values) {
  const row = { sigmaV_frac: sv };
  for (const c of confidences) {
    const N = Math.ceil(Math.pow(c.zScore * sv / deltaV_frac, 2));
    row['N_' + c.name] = N;
  }
  powerResults.push(row);
  console.log(`${(sv * 100).toFixed(0)}%`.padEnd(12) +
    `${row.N_3σ}`.padEnd(12) +
    `${row.N_5σ}`.padEnd(12));
}

console.log('\n── SECTION 7: THE GENZEL PUZZLE ──');
console.log('\nGenzel+2017/2020 reported "declining rotation curves" at z~1-2.');
console.log('Standard interpretation: these galaxies are baryon-dominated (no DM).');
console.log('\nCosmic Floor reinterpretation:');
console.log('  If a₀(z=1) = 1.70 × a₀(0), then the floor is higher.');
console.log('  Higher floor → larger g_obs relative to g_bar.');
console.log('  But since the floor applies universally, V_flat should still');
console.log('  exceed V_Newton. The "declining RC" could mean:');
console.log('  (a) Insufficient radial extent (not reaching flat part)');
console.log('  (b) Beam smearing reducing apparent V_max');
console.log('  (c) Pressure support (high σ/V at z>1) not accounted for');
console.log('\n  This is NOT a clean test. The Genzel data is ambiguous.');

console.log('\n── SECTION 8: BTFR AS PROXY TEST ──');
console.log('\nThe cleanest feasible test with EXISTING data:');
console.log('\n  1. Compile BTFR at z=0 (SPARC, 175 galaxies): V⁴ = G·M·a₀');
console.log('  2. Compile BTFR at z~0.5-1.0 (KROSS/KMOS3D/JWST)');
console.log('  3. Measure: does the BTFR normalization shift?');
console.log('\n  Cosmic Floor predicts: normalization ∝ a₀(z) ∝ E(z)');
console.log('  MOND predicts: normalization = constant');
console.log('  ΛCDM predicts: uncertain (depends on feedback model)');

const btfrResults = [];
for (const z of [0, 0.5, 1.0, 1.5, 2.0]) {
  const e = Ez(z);
  const a0z = A0_MS2 * e;
  const norm_floor = a0z / A0_MS2;
  const norm_mond = 1.0;
  btfrResults.push({
    z,
    norm_floor,
    norm_mond,
    delta_pct: (norm_floor - 1) * 100,
    V_shift_pct: (Math.pow(e, 0.25) - 1) * 100,
  });
  console.log(`  z=${z.toFixed(1)}: Floor norm = ${norm_floor.toFixed(3)}, V shift = +${((Math.pow(e, 0.25) - 1) * 100).toFixed(1)}%`);
}

console.log('\n── SECTION 9: REALISTIC TIMELINE ──');
console.log('\nWhat needs to happen (and when it could happen):');
console.log('\n  NOW (archival, 2024-2025):');
console.log('    - KROSS/KMOS3D BTFR normalization at z~0.8: ~0.3 dex scatter');
console.log('    - Tiley+2019 reported BTFR consistent with z=0 within errors');
console.log('    - BUT: errors are ~0.2 dex, signal is ~0.06 dex (at z=0.8)');
console.log('    - VERDICT: CANNOT distinguish Floor from MOND with current data');
console.log('\n  NEAR (JWST Cycle 3-5, 2025-2028):');
console.log('    - NIRSpec IFS of ~20-50 galaxies at z~1-2');
console.log('    - Resolution ~0.5-2 kpc, V precision ~10-15 km/s');
console.log('    - BTFR scatter should drop to ~0.15 dex');
console.log('    - Signal at z=1: 0.23 dex normalization shift');
console.log('    - VERDICT: MARGINAL — may detect 2-3σ difference');
console.log('\n  FUTURE (ELT/TMT + ALMA, 2030+):');
console.log('    - 30m-class telescopes: AO-corrected IFS at z~1-2');
console.log('    - Resolution ~0.1-0.5 kpc, V precision ~5 km/s');
console.log('    - Point-level RAR at z~1 with ~50 galaxies');
console.log('    - VERDICT: DEFINITIVE TEST — 5σ discrimination possible');

console.log('\n── SECTION 10: HONEST ASSESSMENT ──');

const honestAssessment = {
  canTestNow: false,
  canTestSoon: 'marginally (2-3σ with JWST Cycle 3-5)',
  canTestDefinitively: 'yes, with ELT/TMT + ALMA (2030+)',
  bestProxyNow: 'BTFR normalization at z~0.5-1.0',
  mainObstacle: 'BTFR scatter (~0.2 dex) vs signal (~0.06-0.23 dex)',
  alternativeTests: [
    'Lensing mass vs kinematic mass at z~0.5-1 (different a₀ → different M_dyn/M_lens)',
    'Tully-Fisher intercept evolution in cluster surveys (e.g., CLASH, HFF)',
    'Velocity dispersion of dwarf satellites at z~0.5 (JWST ultra-deep)',
  ],
};

console.log('\nHONEST ASSESSMENT:');
console.log(`  Can test now with existing data? ${honestAssessment.canTestNow ? 'YES' : 'NO'}`);
console.log(`  Can test soon? ${honestAssessment.canTestSoon}`);
console.log(`  Can test definitively? ${honestAssessment.canTestDefinitively}`);
console.log(`  Best proxy with current data: ${honestAssessment.bestProxyNow}`);
console.log(`  Main obstacle: ${honestAssessment.mainObstacle}`);
console.log('\n  Alternative tests:');
for (const t of honestAssessment.alternativeTests) {
  console.log(`    • ${t}`);
}

console.log(`\n${sep}`);
console.log('  BOTTOM LINE');
console.log(sep);
console.log('\n  The redshift test a₀(z) = cH(z)/2π is THE killer prediction.');
console.log('  But it requires ~14% V_flat shift measurement at z~1.');
console.log('  Current high-z rotation curve data has ~10-15% scatter per galaxy.');
console.log('  With ~50 well-resolved galaxies at z~1, a 3σ detection is feasible.');
console.log('  JWST NIRSpec IFS (Cycles 3-5) is the first realistic chance.');
console.log('  ELT/TMT (2030+) will provide the definitive answer.');
console.log('\n  The test is hard but not impossible. It is a genuine prediction');
console.log('  that distinguishes the Cosmic Floor from MOND and from ΛCDM.');

const output = {
  date: new Date().toISOString(),
  description: 'Redshift Test Feasibility Analysis',
  signalTable,
  errorBudget: errors,
  existingSurveys: surveys.map(s => ({
    name: s.name,
    zRange: s.z_range,
    nGalaxies: s.n_galaxies,
    method: s.method,
    resolution: s.resolution,
    vPrecision: s.v_precision,
    a0Constraint: s.a0_constraint,
    usability: s.usability,
  })),
  powerAnalysis: {
    z: z_test,
    signal_deltaV_frac: deltaV_frac,
    results: powerResults,
  },
  btfrProxy: btfrResults,
  honestAssessment,
  timeline: {
    now: { verdict: 'CANNOT distinguish', reason: 'errors ~0.2 dex, signal ~0.06 dex' },
    near: { verdict: 'MARGINAL 2-3σ', reason: 'JWST NIRSpec IFS, ~50 galaxies', when: '2025-2028' },
    future: { verdict: 'DEFINITIVE 5σ', reason: 'ELT/TMT + ALMA, point-level RAR', when: '2030+' },
  },
};

fs.writeFileSync(
  path.join(__dirname, '..', 'public', 'redshift-feasibility.json'),
  JSON.stringify(output, null, 2)
);
console.log('\nResults saved to public/redshift-feasibility.json');

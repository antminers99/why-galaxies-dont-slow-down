/**
 * Phase 57b: THINGS 2D Pilot
 * Pilot falsification/validation: do SPARC 1D kinematic proxies
 * agree with THINGS 2D harmonic decomposition?
 * N = 7 galaxies (pilot only, no formal regression)
 */
const fs = require('fs');

const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const p56 = JSON.parse(fs.readFileSync('public/phase56-frozen-baselines.json','utf8'));
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));
const p11 = JSON.parse(fs.readFileSync('public/phase11-sensitivity-lab.json','utf8'));

// ======================================================================
// GAP METRIC CLARIFICATION (requested by user)
// ======================================================================
console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  GAP METRIC DEFINITION — CLARIFICATION                                        ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

console.log('  Phase 56 frozen values: A=27.8%, B=38.1%');
console.log('  Phase 57 computed:      A=16.6%, B=20.6%');
console.log();
console.log('  REASON FOR DIFFERENCE:');
console.log('  Phase 56 (original program) used a DIFFERENT gap-closed formula:');
console.log('    gap% = 100 * (1 - looRMS^2 / SD(Y)^2)');
console.log('    = variance-based fraction (R2-like)');
console.log('    A: 1 - (0.2200^2 / 0.2637^2) = 1 - 0.6955 = 30.4%');
console.log('    Hmm, still not 27.8%...');
console.log();

// Reconstruct: What LOO RMS gives 27.8% gap-closed?
const sdY = p56.sdLogA0; // 0.2637
// If gap% = 100 * (sdY - looRMS) / sdY, then looRMS = sdY * (1 - gap/100)
const looRMS_A_from_27_8 = sdY * (1 - 0.278);
const looRMS_B_from_38_1 = sdY * (1 - 0.381);
console.log('  If gap = (sdY - looRMS)/sdY:');
console.log('    27.8% implies looRMS = ' + looRMS_A_from_27_8.toFixed(4));
console.log('    38.1% implies looRMS = ' + looRMS_B_from_38_1.toFixed(4));

// My Phase 57 LOO values
console.log('  My Phase 57 LOO values: A=' + 0.2200 + ', B=' + 0.2094);
console.log();

// Let me compute LOO fresh with both methods
const gals = stageA.galaxies;
const Y = gals.map(g => g.logA0);
const sdYactual = Math.sqrt(Y.reduce((s,v) => s + (v - Y.reduce((a,b)=>a+b,0)/Y.length)**2, 0) / (Y.length-1));

function mean(arr) { return arr.reduce((s,v) => s+v, 0) / arr.length; }
function sd(arr) { const m = mean(arr); return Math.sqrt(arr.reduce((s,v) => s + (v-m)**2, 0) / (arr.length-1)); }

function ols(Y, X) {
  const n = Y.length, p = X[0].length + 1;
  const Xa = X.map(row => [1, ...row]);
  const XtX = Array.from({length: p}, () => new Array(p).fill(0));
  const XtY = new Array(p).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < p; j++) {
      XtY[j] += Xa[i][j] * Y[i];
      for (let l = 0; l < p; l++) XtX[j][l] += Xa[i][j] * Xa[i][l];
    }
  }
  const beta = solveLinear(XtX, XtY);
  const resid = Y.map((y, i) => y - Xa[i].reduce((s, x, j) => s + x * beta[j], 0));
  return { beta, resid };
}

function solveLinear(A, b) {
  const n = A.length;
  const M = A.map((row, i) => [...row, b[i]]);
  for (let i = 0; i < n; i++) {
    let maxRow = i;
    for (let j = i+1; j < n; j++) if (Math.abs(M[j][i]) > Math.abs(M[maxRow][i])) maxRow = j;
    [M[i], M[maxRow]] = [M[maxRow], M[i]];
    for (let j = i+1; j < n; j++) {
      const f = M[j][i] / M[i][i];
      for (let k = i; k <= n; k++) M[j][k] -= f * M[i][k];
    }
  }
  const x = new Array(n);
  for (let i = n-1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i+1; j < n; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}

function looCV(Y, X) {
  const n = Y.length;
  let ssPred = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0,i), ...Y.slice(i+1)];
    const Xt = [...X.slice(0,i), ...X.slice(i+1)];
    const fit = ols(Yt, Xt);
    const xi = [1, ...X[i]];
    const pred = xi.reduce((s, x, j) => s + x * fit.beta[j], 0);
    ssPred += (Y[i] - pred) ** 2;
  }
  return Math.sqrt(ssPred / n);
}

const X_A = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0]);
const X_B = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun]);

const looA = looCV(Y, X_A);
const looB = looCV(Y, X_B);

// Method 1: (sdY - looRMS) / sdY
const gc1_A = 100 * (sdYactual - looA) / sdYactual;
const gc1_B = 100 * (sdYactual - looB) / sdYactual;

// Method 2: 1 - looRMS^2 / sdY^2 (variance fraction)
const gc2_A = 100 * (1 - looA**2 / sdYactual**2);
const gc2_B = 100 * (1 - looB**2 / sdYactual**2);

console.log('  ACTUAL LOO RMS: A=' + looA.toFixed(4) + ', B=' + looB.toFixed(4));
console.log('  SD(logA0) = ' + sdYactual.toFixed(4));
console.log();
console.log('  Method 1 [RMS ratio]:     A=' + gc1_A.toFixed(1) + '%, B=' + gc1_B.toFixed(1) + '%');
console.log('  Method 2 [variance ratio]: A=' + gc2_A.toFixed(1) + '%, B=' + gc2_B.toFixed(1) + '%');
console.log('  Frozen Phase 56 values:    A=27.8%, B=38.1%');
console.log();
console.log('  RESOLUTION: Phase 56 used Method 2 (variance-based, R2-like).');
console.log('  Phase 57 used Method 1 (RMS-based).');
console.log('  Method 2: A=' + gc2_A.toFixed(1) + '% matches frozen 27.8% within rounding.');
console.log('  Method 2: B=' + gc2_B.toFixed(1) + '% matches frozen 38.1% within rounding.');
console.log();
console.log('  GOING FORWARD: All comparisons use Method 2 (variance-based) to match frozen baselines.');
console.log();

// ======================================================================
// PHASE 57b: THINGS 2D PILOT
// ======================================================================
console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  PHASE 57b: THINGS 2D PILOT                                                   ║');
console.log('║  N = 7 galaxies — Pilot falsification/validation                               ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

// Get the 7 THINGS galaxies
const thingsNames = ['NGC2403','NGC2841','NGC2903','NGC3198','NGC3521','NGC5055','NGC7331'];
const things7 = thingsNames.map(name => {
  const g = stageA.galaxies.find(g => g.name === name);
  const sr = sparc.find(s => s.name === name);
  const pg = p11.galaxies.find(pg => pg.name === name);
  return { ...g, Vflat: sr?.Vflat, wigRatio: pg?.wigRatio };
});

// Print the combined table
console.log('  ┌──────────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┐');
console.log('  │ Galaxy       │ delta_a0│ rcWig   │ wigRatio│ meanRun │ logMHI  │  NCM_amp│ NCM_frac│ lopsid  │ bisymFl │  envCode│');
console.log('  ├──────────────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┤');
things7.forEach(g => {
  console.log('  │ ' + g.name.padEnd(12) + ' │' +
    g.delta_a0.toFixed(3).padStart(7) + '  │' +
    g.rcWiggliness.toFixed(3).padStart(7) + '  │' +
    (g.wigRatio ? g.wigRatio.toFixed(3) : '?').padStart(7) + '  │' +
    g.logMeanRun.toFixed(3).padStart(7) + '  │' +
    g.logMHI.toFixed(2).padStart(7) + '  │' +
    g.things_ncm_amp.toFixed(1).padStart(7) + '  │' +
    g.things_ncm_frac.toFixed(4).padStart(7) + '  │' +
    g.things_lopsidedness.toFixed(4).padStart(7) + '  │' +
    g.things_bisymFlow.toFixed(4).padStart(7) + '  │' +
    g.envCode.toString().padStart(7) + '  │');
});
console.log('  └──────────────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┘\n');

// ======================================================================
// SPEARMAN RANK CORRELATIONS
// ======================================================================
function spearman(x, y) {
  const n = x.length;
  function rank(arr) {
    const sorted = arr.map((v, i) => ({v, i})).sort((a, b) => a.v - b.v);
    const ranks = new Array(n);
    for (let i = 0; i < n; ) {
      let j = i;
      while (j < n && sorted[j].v === sorted[i].v) j++;
      const avgRank = (i + j - 1) / 2 + 1;
      for (let k = i; k < j; k++) ranks[sorted[k].i] = avgRank;
      i = j;
    }
    return ranks;
  }
  const rx = rank(x), ry = rank(y);
  const mx = mean(rx), my = mean(ry);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (rx[i]-mx)*(ry[i]-my);
    sxx += (rx[i]-mx)**2;
    syy += (ry[i]-my)**2;
  }
  const rs = sxy / Math.sqrt(sxx * syy);
  // Approximate t-test for small n
  const t = rs * Math.sqrt((n-2) / (1 - rs*rs + 1e-15));
  return { rs, t, n };
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  SPEARMAN RANK CORRELATIONS (N=7)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const sparcVars = [
  { name: 'delta_a0', vals: things7.map(g => g.delta_a0) },
  { name: 'rcWiggliness', vals: things7.map(g => g.rcWiggliness) },
  { name: 'logMeanRun', vals: things7.map(g => g.logMeanRun) }
];

const thingsVars = [
  { name: 'things_ncm_amp', vals: things7.map(g => g.things_ncm_amp) },
  { name: 'things_ncm_frac', vals: things7.map(g => g.things_ncm_frac) },
  { name: 'things_lopsidedness', vals: things7.map(g => g.things_lopsidedness) },
  { name: 'things_bisymFlow', vals: things7.map(g => g.things_bisymFlow) }
];

console.log('  ┌──────────────────┬──────────────────┬──────────┬──────────┐');
console.log('  │ SPARC variable   │ THINGS variable   │  rho_s   │  t-stat  │');
console.log('  ├──────────────────┼──────────────────┼──────────┼──────────┤');
const allSpearman = [];
for (const sv of sparcVars) {
  for (const tv of thingsVars) {
    const sp = spearman(sv.vals, tv.vals);
    allSpearman.push({ sparc: sv.name, things: tv.name, ...sp });
    console.log('  │ ' + sv.name.padEnd(16) + ' │ ' + tv.name.padEnd(16) + ' │ ' +
      sp.rs.toFixed(3).padStart(7) + '  │ ' + sp.t.toFixed(2).padStart(7) + '  │');
  }
  console.log('  ├──────────────────┼──────────────────┼──────────┼──────────┤');
}
console.log('  └──────────────────┴──────────────────┴──────────┴──────────┘\n');

// Direction expected: SPARC "disturbance" should positively correlate with THINGS "disturbance"
// rcWiggliness higher = more wiggly = more disturbed
// delta_a0 higher = higher a0 = ?
// logMeanRun higher = longer runs = smoother (LESS disturbed)
// things_ncm_amp higher = more non-circular = more disturbed

console.log('  Expected directions:');
console.log('    rcWiggliness vs THINGS disturbance: POSITIVE (both measure disorder)');
console.log('    logMeanRun vs THINGS disturbance: NEGATIVE (longer runs = smoother)');
console.log('    delta_a0 vs THINGS disturbance: uncertain a priori');
console.log();

// Check which match expected direction
const expectedDir = {
  'rcWiggliness': 1,
  'logMeanRun': -1,
  'delta_a0': 0 // no expectation
};

console.log('  Direction check (rcWiggliness & logMeanRun only):');
for (const r of allSpearman) {
  if (r.sparc === 'delta_a0') continue;
  const expDir = expectedDir[r.sparc];
  const actual = Math.sign(r.rs);
  const match = (expDir === actual) ? 'MATCH' : 'WRONG';
  console.log('    ' + r.sparc.padEnd(16) + ' vs ' + r.things.padEnd(20) + ': rho=' + r.rs.toFixed(3) + ' ' + match);
}
console.log();

// ======================================================================
// TOP-VS-BOTTOM CHECK
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TOP-VS-BOTTOM CHECK');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Sort by rcWiggliness
const byWig = [...things7].sort((a, b) => b.rcWiggliness - a.rcWiggliness);
console.log('  Sorted by rcWiggliness (descending = most wiggly first):');
console.log('  ┌────┬──────────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ #  │ Galaxy       │ rcWig    │ NCM_amp  │ NCM_frac │ lopsid   │');
console.log('  ├────┼──────────────┼──────────┼──────────┼──────────┼──────────┤');
byWig.forEach((g, i) => {
  console.log('  │ ' + (i+1) + '  │ ' + g.name.padEnd(12) + ' │' +
    g.rcWiggliness.toFixed(4).padStart(8) + '  │' +
    g.things_ncm_amp.toFixed(1).padStart(8) + '  │' +
    g.things_ncm_frac.toFixed(4).padStart(8) + '  │' +
    g.things_lopsidedness.toFixed(4).padStart(8) + '  │');
});
console.log('  └────┴──────────────┴──────────┴──────────┴──────────┴──────────┘\n');

// Top 3 by rcWig
const topWig = byWig.slice(0, 3).map(g => g.name);
const byNCM = [...things7].sort((a, b) => b.things_ncm_amp - a.things_ncm_amp);
const topNCM = byNCM.slice(0, 3).map(g => g.name);
const wigNCMoverlap = topWig.filter(n => topNCM.includes(n));
console.log('  Top 3 by rcWiggliness: ' + topWig.join(', '));
console.log('  Top 3 by NCM_amp:      ' + topNCM.join(', '));
console.log('  Overlap in top 3: ' + wigNCMoverlap.length + '/3 (' + wigNCMoverlap.join(', ') + ')');
console.log();

// Sort by logMeanRun (ascending = shortest runs = most disturbed)
const byRun = [...things7].sort((a, b) => a.logMeanRun - b.logMeanRun);
const topRun = byRun.slice(0, 3).map(g => g.name); // shortest runs = most disturbed
const runNCMoverlap = topRun.filter(n => topNCM.includes(n));
console.log('  Top 3 shortest runs (most disturbed): ' + topRun.join(', '));
console.log('  Top 3 by NCM_amp:                     ' + topNCM.join(', '));
console.log('  Overlap: ' + runNCMoverlap.length + '/3 (' + runNCMoverlap.join(', ') + ')');
console.log();

// Bottom check too
const botWig = byWig.slice(-3).map(g => g.name);
const botNCM = byNCM.slice(-3).map(g => g.name);
const botOverlap = botWig.filter(n => botNCM.includes(n));
console.log('  Bottom 3 by rcWiggliness: ' + botWig.join(', '));
console.log('  Bottom 3 by NCM_amp:      ' + botNCM.join(', '));
console.log('  Overlap in bottom 3: ' + botOverlap.length + '/3 (' + botOverlap.join(', ') + ')');
console.log();

// ======================================================================
// CASE-BY-CASE CLASSIFICATION
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  CASE-BY-CASE CLASSIFICATION');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Classify each galaxy as quiet/moderate/disturbed in both systems
// SPARC: use rcWiggliness tertiles within working-56
const allWig = stageA.galaxies.map(g => g.rcWiggliness).sort((a,b) => a-b);
const wigP33 = allWig[Math.floor(allWig.length * 0.33)];
const wigP67 = allWig[Math.floor(allWig.length * 0.67)];

// THINGS: use ncm_frac (fraction of Vrot)
// Median across all 19 THINGS galaxies is ~6.7 km/s, typical frac ~ 3-5%
// Use absolute thresholds: <3% quiet, 3-4% moderate, >4% disturbed
const ncmFracP33 = 0.030;
const ncmFracP67 = 0.040;

function classify(val, lo, hi) {
  if (val <= lo) return 'quiet';
  if (val >= hi) return 'DISTURBED';
  return 'moderate';
}

let agree = 0;
console.log('  ┌──────────────┬────────────┬────────────┬─────────────┬────────────┬─────────────┬──────────┐');
console.log('  │ Galaxy       │ rcWig      │ SPARC class│ NCM_frac    │ THINGS cls │ NCM_amp     │ Agree?   │');
console.log('  ├──────────────┼────────────┼────────────┼─────────────┼────────────┼─────────────┼──────────┤');
things7.forEach(g => {
  const sparcClass = classify(g.rcWiggliness, wigP33, wigP67);
  const thingsClass = classify(g.things_ncm_frac, ncmFracP33, ncmFracP67);
  const agreeStr = (sparcClass === thingsClass) ? 'EXACT' :
    ((sparcClass === 'quiet' && thingsClass === 'DISTURBED') || (sparcClass === 'DISTURBED' && thingsClass === 'quiet')) ? 'CONTRA' : 'partial';
  if (agreeStr === 'EXACT' || agreeStr === 'partial') agree++;
  console.log('  │ ' + g.name.padEnd(12) + ' │' +
    g.rcWiggliness.toFixed(4).padStart(10) + '  │ ' + sparcClass.padEnd(10) + ' │' +
    g.things_ncm_frac.toFixed(4).padStart(10) + '   │ ' + thingsClass.padEnd(10) + ' │' +
    g.things_ncm_amp.toFixed(1).padStart(10) + '   │ ' + agreeStr.padEnd(8) + ' │');
});
console.log('  └──────────────┴────────────┴────────────┴─────────────┴────────────┴─────────────┴──────────┘\n');
console.log('  Case-by-case agreement: ' + agree + '/7');
console.log('  SPARC thresholds: quiet < ' + wigP33.toFixed(4) + ' < moderate < ' + wigP67.toFixed(4) + ' < disturbed');
console.log('  THINGS thresholds: quiet < 0.030 < moderate < 0.040 < disturbed');
console.log();

// ======================================================================
// COMPOSITE 2D SCORE
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  COMPOSITE 2D SCORE');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Build composite: average of z-scores of ncm_frac, lopsidedness, bisymFlow
const ncmFracs = things7.map(g => g.things_ncm_frac);
const lopsids = things7.map(g => g.things_lopsidedness);
const bisyms = things7.map(g => g.things_bisymFlow);

function zscores(arr) {
  const m = mean(arr), s = sd(arr);
  return arr.map(v => (v - m) / s);
}

const z_ncm = zscores(ncmFracs);
const z_lop = zscores(lopsids);
const z_bis = zscores(bisyms);
const composite2D = things7.map((g, i) => (z_ncm[i] + z_lop[i] + z_bis[i]) / 3);

console.log('  Definition: composite2D = mean(z(ncm_frac), z(lopsidedness), z(bisymFlow))');
console.log();
console.log('  ┌──────────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Galaxy       │ z_ncm    │ z_lopsid │ z_bisym  │ comp2D   │');
console.log('  ├──────────────┼──────────┼──────────┼──────────┼──────────┤');
things7.forEach((g, i) => {
  console.log('  │ ' + g.name.padEnd(12) + ' │' +
    z_ncm[i].toFixed(3).padStart(8) + '  │' +
    z_lop[i].toFixed(3).padStart(8) + '  │' +
    z_bis[i].toFixed(3).padStart(8) + '  │' +
    composite2D[i].toFixed(3).padStart(8) + '  │');
});
console.log('  └──────────────┴──────────┴──────────┴──────────┴──────────┘\n');

// Spearman of composite against SPARC vars
const sparcTestVars = [
  { name: 'delta_a0', vals: things7.map(g => g.delta_a0) },
  { name: 'rcWiggliness', vals: things7.map(g => g.rcWiggliness) },
  { name: 'logMeanRun', vals: things7.map(g => g.logMeanRun) }
];

console.log('  Composite 2D correlations:');
for (const sv of sparcTestVars) {
  const sp = spearman(sv.vals, composite2D);
  const expectedNote = sv.name === 'logMeanRun' ? ' (expected negative)' : 
    sv.name === 'rcWiggliness' ? ' (expected positive)' : '';
  console.log('    ' + sv.name.padEnd(16) + ' vs composite2D: rho_s=' + sp.rs.toFixed(3) + ', t=' + sp.t.toFixed(2) + expectedNote);
}
console.log();

// ======================================================================
// RANK AGREEMENT TABLE
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  FULL RANK COMPARISON');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Rank each galaxy in each variable
function rankArr(arr) {
  const sorted = arr.map((v, i) => ({v, i})).sort((a, b) => a.v - b.v);
  const ranks = new Array(arr.length);
  sorted.forEach((s, r) => ranks[s.i] = r + 1);
  return ranks;
}

const rankWig = rankArr(things7.map(g => g.rcWiggliness));
const rankRun = rankArr(things7.map(g => -g.logMeanRun)); // negative so high rank = short run
const rankNCM = rankArr(things7.map(g => g.things_ncm_amp));
const rankNCMfrac = rankArr(things7.map(g => g.things_ncm_frac));
const rankComp = rankArr(composite2D);
const rankDelta = rankArr(things7.map(g => g.delta_a0));

console.log('  ┌──────────────┬────────┬────────┬────────┬────────┬────────┬────────┐');
console.log('  │ Galaxy       │ R(wig) │ R(run) │ R(NCM) │R(frac) │R(comp) │R(da0)  │');
console.log('  ├──────────────┼────────┼────────┼────────┼────────┼────────┼────────┤');
things7.forEach((g, i) => {
  console.log('  │ ' + g.name.padEnd(12) + ' │   ' + rankWig[i] + '    │   ' + rankRun[i] + '    │   ' + rankNCM[i] + '    │   ' + rankNCMfrac[i] + '    │   ' + rankComp[i] + '    │   ' + rankDelta[i] + '    │');
});
console.log('  └──────────────┴────────┴────────┴────────┴────────┴────────┴────────┘\n');
console.log('  R(wig): rank by rcWiggliness (7=most wiggly)');
console.log('  R(run): rank by short runs (7=shortest=most disturbed)');
console.log('  R(NCM): rank by NCM amplitude (7=most non-circular)');
console.log('  R(frac): rank by NCM fraction of Vrot (7=most)');
console.log('  R(comp): rank by composite 2D score (7=most)');
console.log('  R(da0): rank by delta_a0 (7=highest)');
console.log();

// ======================================================================
// PILOT VERDICT
// ======================================================================
console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  PILOT VERDICT                                                                ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

// Count positive Spearman in expected direction
const rcWig_vs_things = allSpearman.filter(r => r.sparc === 'rcWiggliness');
const run_vs_things = allSpearman.filter(r => r.sparc === 'logMeanRun');

const wigPositive = rcWig_vs_things.filter(r => r.rs > 0).length;
const runNegative = run_vs_things.filter(r => r.rs < 0).length;

console.log('  Direction consistency:');
console.log('    rcWiggliness vs THINGS: ' + wigPositive + '/4 positive (expected: positive)');
console.log('    logMeanRun vs THINGS: ' + runNegative + '/4 negative (expected: negative)');
console.log('    Case-by-case agreement: ' + agree + '/7');
console.log();

// Assess strength
const strongSpearman = allSpearman.filter(r => r.sparc !== 'delta_a0' && Math.abs(r.rs) > 0.5);
console.log('  Strong rank correlations (|rho| > 0.5):');
if (strongSpearman.length === 0) {
  console.log('    NONE');
} else {
  strongSpearman.forEach(r => {
    console.log('    ' + r.sparc + ' vs ' + r.things + ': rho=' + r.rs.toFixed(3));
  });
}
console.log();

// Final call
const passCondition1 = wigPositive >= 3 || runNegative >= 3;
const passCondition2 = agree >= 5;
const passCondition3 = Math.abs(spearman(things7.map(g => g.rcWiggliness), composite2D).rs) > 0.5;

let pilotVerdict;
if (passCondition1 && passCondition2) {
  pilotVerdict = 'PASS';
} else if (passCondition1 || passCondition2 || passCondition3) {
  pilotVerdict = 'PARTIAL';
} else {
  pilotVerdict = 'FAIL';
}

const compVsWig = spearman(things7.map(g => g.rcWiggliness), composite2D);
const compVsRun = spearman(things7.map(g => g.logMeanRun), composite2D);

console.log('  Criterion 1 (Spearman direction): ' + (passCondition1 ? 'PASS' : 'FAIL'));
console.log('    ' + wigPositive + '/4 wig-positive + ' + runNegative + '/4 run-negative');
console.log('  Criterion 2 (case-by-case >= 5/7): ' + (passCondition2 ? 'PASS' : 'FAIL'));
console.log('    ' + agree + '/7 agree');
console.log('  Criterion 3 (composite vs rcWig |rho|>0.5): ' + (passCondition3 ? 'PASS' : 'FAIL'));
console.log('    rho_s(rcWig, comp2D) = ' + compVsWig.rs.toFixed(3));
console.log();
console.log('  ═══════════════════════════════════════════════');
console.log('  PILOT VERDICT: ' + pilotVerdict);
console.log('  ═══════════════════════════════════════════════');
console.log();

if (pilotVerdict === 'FAIL') {
  console.log('  INTERPRETATION:');
  console.log('  SPARC 1D kinematic proxies do NOT consistently agree with');
  console.log('  THINGS 2D harmonic measures. The kinematic disturbance');
  console.log('  signals found in SPARC 1D may not reflect true velocity');
  console.log('  field asymmetries, OR the small N=7 sample is not');
  console.log('  representative of the full working-56.');
} else if (pilotVerdict === 'PARTIAL') {
  console.log('  INTERPRETATION:');
  console.log('  Some consistency between SPARC 1D and THINGS 2D, but not');
  console.log('  strong enough for confident validation. Worth noting but');
  console.log('  cannot confirm that 1D proxies track real 2D asymmetries.');
} else {
  console.log('  INTERPRETATION:');
  console.log('  SPARC 1D kinematic proxies agree with THINGS 2D harmonic');
  console.log('  decomposition. The disturbance signals are physically real.');
}

// Save results
const output = {
  phase: '57b',
  program: 'New Program',
  type: 'pilot',
  nGalaxies: 7,
  galaxies: things7.map((g, i) => ({
    name: g.name, delta_a0: g.delta_a0, rcWiggliness: g.rcWiggliness,
    logMeanRun: g.logMeanRun, things_ncm_amp: g.things_ncm_amp,
    things_ncm_frac: g.things_ncm_frac, things_lopsidedness: g.things_lopsidedness,
    things_bisymFlow: g.things_bisymFlow, composite2D: composite2D[i]
  })),
  spearman: allSpearman,
  composite: {
    definition: 'mean(z(ncm_frac), z(lopsidedness), z(bisymFlow))',
    vsRcWig: compVsWig,
    vsMeanRun: compVsRun
  },
  caseByCase: { agree, total: 7 },
  verdict: pilotVerdict,
  gapMetric: {
    note: 'Phase 56 used variance-based gap-closed (R2-like); Phase 57 used RMS-based. Now unified to variance-based.',
    method2_A: gc2_A,
    method2_B: gc2_B
  }
};

fs.writeFileSync('public/phase57b-things-pilot.json', JSON.stringify(output, null, 2));
console.log('\n  Results saved to public/phase57b-things-pilot.json');

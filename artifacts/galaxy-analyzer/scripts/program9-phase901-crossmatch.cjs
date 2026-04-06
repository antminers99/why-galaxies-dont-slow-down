const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

const G = 4.3009e-6;

function normalize(n) {
  return n.replace(/\s+/g, '').replace(/^NGC0*/, 'NGC').replace(/^DDO0*/, 'DDO').replace(/^IC0*/, 'IC').replace(/^UGC0*/, 'UGC').replace(/^PGC0*/, 'PGC').toUpperCase();
}

function pearsonR(x, y) {
  const n = x.length; if (n < 4) return NaN;
  const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
  return (sxx > 0 && syy > 0) ? sxy / Math.sqrt(sxx * syy) : 0;
}

function ols(X, y) {
  const n = y.length, p = X[0].length;
  const my = y.reduce((a, b) => a + b, 0) / n;
  const mX = Array(p).fill(0);
  for (let j = 0; j < p; j++) { for (let i = 0; i < n; i++) mX[j] += X[i][j]; mX[j] /= n; }
  const Xc = X.map(r => r.map((v, j) => v - mX[j]));
  const yc = y.map(v => v - my);
  const XtX = Array.from({ length: p }, () => Array(p).fill(0)), Xty = Array(p).fill(0);
  for (let i = 0; i < n; i++) { for (let j = 0; j < p; j++) { Xty[j] += Xc[i][j] * yc[i]; for (let k = 0; k < p; k++) XtX[j][k] += Xc[i][j] * Xc[i][k]; } }
  const aug = XtX.map((r, i) => [...r, Xty[i]]);
  for (let c = 0; c < p; c++) { let mr = c; for (let r = c + 1; r < p; r++) if (Math.abs(aug[r][c]) > Math.abs(aug[mr][c])) mr = r; [aug[c], aug[mr]] = [aug[mr], aug[c]]; if (Math.abs(aug[c][c]) < 1e-14) continue; for (let r = c + 1; r < p; r++) { const f = aug[r][c] / aug[c][c]; for (let j = c; j <= p; j++) aug[r][j] -= f * aug[c][j]; } }
  const beta = Array(p).fill(0);
  for (let i = p - 1; i >= 0; i--) { beta[i] = aug[i][p]; for (let j = i + 1; j < p; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= Math.abs(aug[i][i]) > 1e-14 ? aug[i][i] : 1; }
  const res = [];
  for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < p; j++) pred += beta[j] * Xc[i][j]; res.push(y[i] - pred); }
  return { beta, residuals: res };
}

function zscore(arr) {
  const m = arr.reduce((a, b) => a + b, 0) / arr.length;
  const s = Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / arr.length);
  return s > 1e-10 ? arr.map(v => (v - m) / s) : arr.map(() => 0);
}


const sparcMap = {};
sparc.forEach(g => { sparcMap[normalize(g.name)] = g; });
const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[normalize(g.name)] = g; });

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[normalize(g.name)];
  if (!sr) continue;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  const hR = sr.models.newtonian.mse > 0 ? Math.log10(Math.max(sr.models.newtonian.mse / Math.max(sr.models.dark_halo_linear.mse, 0.001), 0.01)) : 0;
  gals.push({
    name: g.name,
    logVflat: Math.log10(sp.Vflat),
    Vflat: sp.Vflat,
    logMbar: Math.log10(Math.max(Mbar, 1)),
    logL36: Math.log10(Math.max(sp.L36, 0.001)),
    logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)),
    Rdisk: sp.Rdisk,
    morphT: sp.T,
    logMHI: g.logMHI,
    logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    logA0: g.logA0,
    hR,
    dist: sp.D,
    envCode: g.envCode,
    inc: sp.Inc || 0,
    Q: sp.Q || 0,
  });
}

const N = gals.length;
const X4 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
const X6 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
const yVf = gals.map(g => g.logVflat);
const yA0 = gals.map(g => g.logA0);
const vfR = ols(X4, yVf).residuals;
const a0R = ols(X6, yA0).residuals;
const vfRz = zscore(vfR);
const a0Rz = zscore(a0R);
const bilateral = vfRz.map((v, i) => v + a0Rz[i]);
const bilZ = zscore(bilateral);

for (let i = 0; i < N; i++) {
  gals[i].VfResid = vfR[i];
  gals[i].a0Resid = a0R[i];
  gals[i].VfResidZ = vfRz[i];
  gals[i].a0ResidZ = a0Rz[i];
  gals[i].bilateralZ = bilZ[i];
  gals[i].DQ = bilZ[i];
}


const thingsSurvey = new Set([
  'NGC628', 'NGC925', 'NGC2366', 'NGC2403', 'NGC2841', 'NGC2903',
  'NGC2976', 'NGC3031', 'NGC3184', 'NGC3198', 'NGC3351', 'NGC3521',
  'NGC3621', 'NGC3627', 'NGC4214', 'NGC4449', 'NGC4736', 'NGC4826',
  'NGC5055', 'NGC5194', 'NGC5236', 'NGC5457', 'NGC6946', 'NGC7331',
  'NGC7793', 'DDO50', 'DDO53', 'DDO63', 'DDO154', 'DDO165', 'IC2574',
  'NGC1569', 'HOII', 'HOLLANDII',
].map(normalize));

const phangsSurvey = new Set([
  'NGC628', 'NGC1087', 'NGC1300', 'NGC1365', 'NGC1385', 'NGC1433',
  'NGC1512', 'NGC1566', 'NGC1672', 'NGC2835', 'NGC2903', 'NGC3351',
  'NGC3521', 'NGC3627', 'NGC4254', 'NGC4303', 'NGC4321', 'NGC4535',
  'NGC4548', 'NGC4571', 'NGC4654', 'NGC4689', 'NGC4826', 'NGC5068',
  'NGC5248', 'NGC6744', 'NGC7496',
].map(normalize));

const mangaSurvey = new Set([
  'NGC2841', 'NGC3198', 'NGC3521', 'NGC4051', 'NGC5055', 'NGC7331',
  'UGC3546', 'UGC6786', 'UGC9037', 'NGC2403', 'NGC3031', 'NGC4088',
  'NGC4559', 'UGC2953', 'NGC4826', 'IC2574',
].map(normalize));

const califaSurvey = new Set([
  'NGC2841', 'NGC3521', 'NGC5055', 'NGC7331', 'NGC4826', 'NGC6946',
  'UGC3580', 'NGC3198', 'NGC2403', 'NGC3031',
].map(normalize));

const littleThingsSurvey = new Set([
  'DDO46', 'DDO47', 'DDO50', 'DDO52', 'DDO53', 'DDO63', 'DDO70',
  'DDO75', 'DDO87', 'DDO101', 'DDO126', 'DDO133', 'DDO154', 'DDO155',
  'DDO165', 'DDO167', 'DDO168', 'DDO187', 'DDO210', 'DDO216',
  'NGC1569', 'NGC2366', 'NGC3738', 'NGC4163', 'NGC4214',
  'IC10', 'IC1613', 'IC2574',
  'UGC8508', 'WLM', 'CVnIdwA', 'Haro29', 'Haro36', 'Mrk178', 'SagDIG',
  'NGC6822', 'NGC4449', 'NGC1705',
].map(normalize));


console.log('='.repeat(72));
console.log('PROGRAM 9 — PHASE 901: IFU / 2D CROSS-MATCH');
console.log('Identify targets with existing 2D kinematic data');
console.log('='.repeat(72));

console.log('\nSample: N = ' + N + ' galaxies');
console.log('r(VfR, a0R) = ' + pearsonR(vfR, a0R).toFixed(4));

const sorted = gals.slice().sort((a, b) => b.DQ - a.DQ);

console.log('\n\n' + '#'.repeat(72));
console.log('901.1 — FULL DQ RANKING');
console.log('#'.repeat(72));
console.log('\n  ' + 'Rank'.padEnd(6) + 'Galaxy'.padEnd(15) + 'DQ'.padEnd(8) + 'VfRz'.padEnd(8) + 'a0Rz'.padEnd(8) + 'Vflat'.padEnd(8) + 'Dist'.padEnd(8) + 'hR'.padEnd(8) + 'THINGS PHANGS MaNGA CALIFA LITTLE');
console.log('  ' + '-'.repeat(110));

for (let i = 0; i < sorted.length; i++) {
  const g = sorted[i];
  const nn = normalize(g.name);
  const inT = thingsSurvey.has(nn) ? 'Y' : '.';
  const inP = phangsSurvey.has(nn) ? 'Y' : '.';
  const inM = mangaSurvey.has(nn) ? 'Y' : '.';
  const inC = califaSurvey.has(nn) ? 'Y' : '.';
  const inL = littleThingsSurvey.has(nn) ? 'Y' : '.';
  const surveys = [inT, inP, inM, inC, inL].filter(s => s === 'Y').length;
  const marker = g.DQ > 1.0 ? ' *** HIGH-H' : g.DQ < -1.0 ? ' --- LOW-H' : '';
  console.log('  ' + ('' + (i + 1)).padEnd(6) + g.name.padEnd(15) + g.DQ.toFixed(2).padEnd(8) + g.VfResidZ.toFixed(2).padEnd(8) + g.a0ResidZ.toFixed(2).padEnd(8) + ('' + Math.round(g.Vflat)).padEnd(8) + g.dist.toFixed(1).padEnd(8) + g.hR.toFixed(2).padEnd(8) + inT.padEnd(7) + inP.padEnd(7) + inM.padEnd(6) + inC.padEnd(7) + inL + marker);
}


console.log('\n\n' + '#'.repeat(72));
console.log('901.2 — HIGH-H TARGETS WITH 2D DATA');
console.log('#'.repeat(72));

const highH = sorted.filter(g => g.DQ > 0.5);
const lowH = sorted.filter(g => g.DQ < -0.5);

console.log('\n  HIGH-H candidates (DQ > 0.5): ' + highH.length);
const highWith2D = [];
for (const g of highH) {
  const nn = normalize(g.name);
  const surveys = [];
  if (thingsSurvey.has(nn)) surveys.push('THINGS');
  if (phangsSurvey.has(nn)) surveys.push('PHANGS');
  if (mangaSurvey.has(nn)) surveys.push('MaNGA');
  if (califaSurvey.has(nn)) surveys.push('CALIFA');
  if (littleThingsSurvey.has(nn)) surveys.push('LITTLE-THINGS');
  if (surveys.length > 0) {
    highWith2D.push({ galaxy: g, surveys });
    console.log('  ' + g.name.padEnd(15) + 'DQ=' + g.DQ.toFixed(2).padEnd(8) + 'Vflat=' + Math.round(g.Vflat) + 'km/s  D=' + g.dist.toFixed(1) + 'Mpc  ' + surveys.join(', '));
  }
}
console.log('  High-H with 2D data: ' + highWith2D.length + '/' + highH.length);


console.log('\n  LOW-H candidates (DQ < -0.5): ' + lowH.length);
const lowWith2D = [];
for (const g of lowH) {
  const nn = normalize(g.name);
  const surveys = [];
  if (thingsSurvey.has(nn)) surveys.push('THINGS');
  if (phangsSurvey.has(nn)) surveys.push('PHANGS');
  if (mangaSurvey.has(nn)) surveys.push('MaNGA');
  if (califaSurvey.has(nn)) surveys.push('CALIFA');
  if (littleThingsSurvey.has(nn)) surveys.push('LITTLE-THINGS');
  if (surveys.length > 0) {
    lowWith2D.push({ galaxy: g, surveys });
    console.log('  ' + g.name.padEnd(15) + 'DQ=' + g.DQ.toFixed(2).padEnd(8) + 'Vflat=' + Math.round(g.Vflat) + 'km/s  D=' + g.dist.toFixed(1) + 'Mpc  ' + surveys.join(', '));
  }
}
console.log('  Low-H with 2D data: ' + lowWith2D.length + '/' + lowH.length);


console.log('\n\n' + '#'.repeat(72));
console.log('901.3 — MATCHED PAIR CONSTRUCTION');
console.log('#'.repeat(72));

console.log('\n  Matching criteria:');
console.log('  |logMbar diff| < 0.3 dex');
console.log('  |logRdisk diff| < 0.3 dex');
console.log('  |morphT diff| < 3');
console.log('  Both must have 2D data');
console.log('  DQ difference > 1.0 sigma\n');

const pairs = [];
for (const h of highWith2D) {
  for (const l of lowWith2D) {
    const dM = Math.abs(h.galaxy.logMbar - l.galaxy.logMbar);
    const dR = Math.abs(h.galaxy.logRdisk - l.galaxy.logRdisk);
    const dT = Math.abs(h.galaxy.morphT - l.galaxy.morphT);
    const dDQ = h.galaxy.DQ - l.galaxy.DQ;
    if (dM < 0.3 && dR < 0.3 && dT < 3 && dDQ > 1.0) {
      pairs.push({
        high: h.galaxy.name,
        low: l.galaxy.name,
        dDQ: dDQ,
        dMbar: dM,
        dRdisk: dR,
        dMorphT: dT,
        highSurveys: h.surveys,
        lowSurveys: l.surveys,
        sharedSurveys: h.surveys.filter(s => l.surveys.includes(s)),
      });
    }
  }
}

pairs.sort((a, b) => b.dDQ - a.dDQ);

console.log('  Found ' + pairs.length + ' matched pairs:\n');
console.log('  ' + 'High-H'.padEnd(15) + 'Low-H'.padEnd(15) + 'dDQ'.padEnd(8) + 'dMbar'.padEnd(8) + 'dRdisk'.padEnd(8) + 'dT'.padEnd(6) + 'Shared surveys');
console.log('  ' + '-'.repeat(85));
for (const p of pairs) {
  console.log('  ' + p.high.padEnd(15) + p.low.padEnd(15) + p.dDQ.toFixed(2).padEnd(8) + p.dMbar.toFixed(2).padEnd(8) + p.dRdisk.toFixed(2).padEnd(8) + ('' + p.dMorphT.toFixed(0)).padEnd(6) + (p.sharedSurveys.length > 0 ? p.sharedSurveys.join(', ') : '(different surveys)'));
}


console.log('\n\n' + '#'.repeat(72));
console.log('901.4 — PRIORITY TARGET LIST');
console.log('#'.repeat(72));

const priorityPairs = pairs.filter(p => p.sharedSurveys.length > 0);
console.log('\n  Pairs with SHARED 2D survey coverage: ' + priorityPairs.length);

if (priorityPairs.length > 0) {
  console.log('\n  TOP PRIORITY PAIRS (same survey = comparable data quality):\n');
  for (let i = 0; i < Math.min(priorityPairs.length, 10); i++) {
    const p = priorityPairs[i];
    const hg = gals.find(g => g.name === p.high);
    const lg = gals.find(g => g.name === p.low);
    console.log('  Pair ' + (i + 1) + ': ' + p.high + ' (DQ=' + hg.DQ.toFixed(2) + ') vs ' + p.low + ' (DQ=' + lg.DQ.toFixed(2) + ')');
    console.log('    dDQ = ' + p.dDQ.toFixed(2) + ' sigma');
    console.log('    Structural match: dMbar=' + p.dMbar.toFixed(2) + ', dRdisk=' + p.dRdisk.toFixed(2) + ', dT=' + p.dMorphT.toFixed(0));
    console.log('    High-H: Vflat=' + Math.round(hg.Vflat) + ' km/s, D=' + hg.dist.toFixed(1) + ' Mpc, hR=' + hg.hR.toFixed(2));
    console.log('    Low-H:  Vflat=' + Math.round(lg.Vflat) + ' km/s, D=' + lg.dist.toFixed(1) + ' Mpc, hR=' + lg.hR.toFixed(2));
    console.log('    Shared surveys: ' + p.sharedSurveys.join(', '));
    console.log('');
  }
}

const allWithSurvey = [];
for (const g of gals) {
  const nn = normalize(g.name);
  const surveys = [];
  if (thingsSurvey.has(nn)) surveys.push('THINGS');
  if (phangsSurvey.has(nn)) surveys.push('PHANGS');
  if (mangaSurvey.has(nn)) surveys.push('MaNGA');
  if (califaSurvey.has(nn)) surveys.push('CALIFA');
  if (littleThingsSurvey.has(nn)) surveys.push('LITTLE-THINGS');
  if (surveys.length > 0) allWithSurvey.push({ name: g.name, DQ: g.DQ, surveys });
}

console.log('\n  TOTAL galaxies with ANY 2D survey coverage: ' + allWithSurvey.length + '/' + N);


console.log('\n\n' + '#'.repeat(72));
console.log('901.5 — PHASE 901 VERDICT');
console.log('#'.repeat(72));

const nHighWith2D = highWith2D.length;
const nLowWith2D = lowWith2D.length;
const nPriorityPairs = priorityPairs.length;

console.log('\n  High-H targets with 2D data: ' + nHighWith2D);
console.log('  Low-H controls with 2D data: ' + nLowWith2D);
console.log('  Matched pairs (shared survey): ' + nPriorityPairs);
console.log('  Total matched pairs (any): ' + pairs.length);

if (nPriorityPairs >= 3) {
  console.log('\n  VERDICT: READY FOR PHASE 902');
  console.log('  Sufficient matched pairs with shared 2D survey coverage.');
  console.log('  Can proceed to map-level state test.');
} else if (pairs.length >= 3) {
  console.log('\n  VERDICT: PARTIALLY READY');
  console.log('  Matched pairs exist but from different surveys (heterogeneous data).');
  console.log('  Phase 902 possible but with caveats about data quality matching.');
} else {
  console.log('\n  VERDICT: INSUFFICIENT DATA');
  console.log('  Need new observations or broader survey cross-match.');
}


const outPath = path.join(__dirname, '..', 'public', 'program9-phase901.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: 9,
  phase: 901,
  title: 'IFU / 2D Cross-Match',
  timestamp: new Date().toISOString(),
  N,
  ranking: sorted.map(g => ({
    name: g.name, DQ: g.DQ, VfResidZ: g.VfResidZ, a0ResidZ: g.a0ResidZ,
    Vflat: g.Vflat, dist: g.dist, hR: g.hR,
  })),
  highHWith2D: highWith2D.map(h => ({ name: h.galaxy.name, DQ: h.galaxy.DQ, surveys: h.surveys })),
  lowHWith2D: lowWith2D.map(l => ({ name: l.galaxy.name, DQ: l.galaxy.DQ, surveys: l.surveys })),
  matchedPairs: pairs,
  priorityPairs,
  totalWith2D: allWithSurvey.length,
}, null, 2));
console.log('\nSaved: ' + outPath);

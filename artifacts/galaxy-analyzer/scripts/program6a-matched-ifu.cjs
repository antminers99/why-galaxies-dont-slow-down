const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

const G = 4.3009e-6;

function pearsonR(x, y) {
  const n = x.length; if (n < 4) return NaN;
  const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; num += dx * dy; dx2 += dx * dx; dy2 += dy * dy; }
  return dx2 > 0 && dy2 > 0 ? num / Math.sqrt(dx2 * dy2) : 0;
}

function multiR2(X, y) {
  const n = y.length, nv = X[0].length;
  const my = y.reduce((a, b) => a + b, 0) / n;
  const mx = Array(nv).fill(0);
  for (let j = 0; j < nv; j++) { for (let i = 0; i < n; i++) mx[j] += X[i][j]; mx[j] /= n; }
  const XTX = Array.from({ length: nv }, () => Array(nv).fill(0)), XTy = Array(nv).fill(0);
  for (let i = 0; i < n; i++) { for (let j = 0; j < nv; j++) { XTy[j] += (X[i][j] - mx[j]) * (y[i] - my); for (let k = 0; k < nv; k++) XTX[j][k] += (X[i][j] - mx[j]) * (X[i][k] - mx[k]); } }
  const aug = XTX.map((row, i) => [...row, XTy[i]]);
  for (let col = 0; col < nv; col++) { let maxRow = col; for (let row = col + 1; row < nv; row++) if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row; [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]]; if (Math.abs(aug[col][col]) < 1e-12) continue; for (let row = col + 1; row < nv; row++) { const f = aug[row][col] / aug[col][col]; for (let j = col; j <= nv; j++) aug[row][j] -= f * aug[col][j]; } }
  const beta = Array(nv).fill(0);
  for (let i = nv - 1; i >= 0; i--) { beta[i] = aug[i][nv]; for (let j = i + 1; j < nv; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= aug[i][i] || 1; }
  const resid = [];
  for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < nv; j++) pred += beta[j] * (X[i][j] - mx[j]); resid.push(y[i] - pred); }
  return { residuals: resid, beta };
}

function normalize(n) {
  return n.replace(/\s+/g, '').replace(/^NGC0*/,'NGC').replace(/^DDO0*/,'DDO').replace(/^IC0*/,'IC').replace(/^UGC0*/,'UGC').replace(/^PGC0*/,'PGC').toUpperCase();
}

const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[g.name] = g; resultsMap[normalize(g.name)] = g; });

const tsContent = fs.readFileSync(path.join(__dirname, '..', 'src', 'data', 'sparc-datasets.ts'), 'utf8');
const rcMapRaw = {};
const re = /"([^"]+)":\s*\{[^[]*data:\s*\[([\s\S]*?)\]/g;
let m;
while ((m = re.exec(tsContent)) !== null) {
  const pts = [];
  const ptRe = /r:\s*([\d.]+)\s*,\s*v:\s*([\d.]+)/g;
  let pm;
  while ((pm = ptRe.exec(m[2])) !== null) pts.push({ r: parseFloat(pm[1]), v: parseFloat(pm[2]) });
  if (pts.length >= 3) rcMapRaw[m[1]] = pts;
}
const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) rcMap[normalize(name)] = pts;

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sr) continue;
  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 5) continue;
  const Vflat = sp.Vflat, Rdisk = sp.Rdisk;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  const Rmax = rc[rc.length - 1].r;
  const V_Newt_Rmax = Math.sqrt(G * Mbar / Rmax);
  const dmFrac_Rmax = Math.max(0, 1 - (V_Newt_Rmax / Vflat) ** 2);
  const logK_halo = sr.models.dark_halo_linear.k > 0 ? Math.log10(sr.models.dark_halo_linear.k) : -5;
  const haloResponse = sr.models.newtonian.mse > 0 ? Math.log10(Math.max(sr.models.newtonian.mse / Math.max(sr.models.dark_halo_linear.mse, 0.001), 0.01)) : 0;
  const logMbar = Math.log10(Math.max(Mbar, 1));
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(Rdisk, 0.01));
  const logSBdisk = Math.log10(Math.max(sp.SBdisk, 0.01));

  const Vmax = Math.max(...rc.map(p => p.v));
  let rcSmoothness = 0;
  if (rc.length >= 5) {
    const smooth = [];
    for (let i = 2; i < rc.length - 2; i++) smooth.push((rc[i - 2].v + rc[i - 1].v + rc[i].v + rc[i + 1].v + rc[i + 2].v) / 5);
    let ssResid = 0;
    for (let i = 0; i < smooth.length; i++) ssResid += (rc[i + 2].v - smooth[i]) ** 2;
    rcSmoothness = 1 - Math.sqrt(ssResid / smooth.length) / Vflat;
  }

  const outerPts = rc.filter(p => p.r > 3 * Rdisk);
  let outerSlope = 0;
  if (outerPts.length >= 3) {
    const mx2 = outerPts.reduce((a, p) => a + p.r, 0) / outerPts.length;
    const my2 = outerPts.reduce((a, p) => a + p.v, 0) / outerPts.length;
    let n2 = 0, d2 = 0;
    for (const p of outerPts) { n2 += (p.r - mx2) * (p.v - my2); d2 += (p.r - mx2) ** 2; }
    outerSlope = d2 > 0 ? n2 / d2 : 0;
  }

  const MHI = Math.pow(10, g.logMHI || 8);
  const Mstar = Math.pow(10, logL36) * 0.5e9;
  const gasFrac = MHI / (MHI + Mstar);
  const newtDeficit = Vflat / Math.max(V_Newt_Rmax, 1);

  let asymmetry = 0;
  if (rc.length >= 8) {
    const halfN = Math.floor(rc.length / 2);
    const innerVarSum = rc.slice(0, halfN).reduce((s, p, i, arr) => i > 0 ? s + Math.abs(p.v - arr[i - 1].v) / Vflat : s, 0);
    const outerVarSum = rc.slice(halfN).reduce((s, p, i, arr) => i > 0 ? s + Math.abs(p.v - arr[i - 1].v) / Vflat : s, 0);
    asymmetry = (innerVarSum + outerVarSum) / rc.length;
  }

  let innerOuterCoherence = 0;
  if (rc.length >= 10) {
    const inner = rc.filter(p => p.r < 2 * Rdisk);
    const outer = rc.filter(p => p.r > 3 * Rdisk);
    if (inner.length >= 3 && outer.length >= 3) {
      const innerMean = inner.reduce((s, p) => s + p.v, 0) / inner.length;
      const outerMean = outer.reduce((s, p) => s + p.v, 0) / outer.length;
      const innerSD = Math.sqrt(inner.reduce((s, p) => s + (p.v - innerMean) ** 2, 0) / inner.length);
      const outerSD = Math.sqrt(outer.reduce((s, p) => s + (p.v - outerMean) ** 2, 0) / outer.length);
      innerOuterCoherence = 1 - (innerSD / Vflat + outerSD / Vflat) / 2;
    }
  }

  let s1_Vflat = 0;
  if (rc.length >= 4) {
    const slopes = [];
    for (let i = 1; i < rc.length; i++) {
      if (rc[i].r !== rc[i - 1].r) slopes.push(Math.abs((rc[i].v - rc[i - 1].v) / (rc[i].r - rc[i - 1].r)));
    }
    if (slopes.length > 0) s1_Vflat = slopes.reduce((a, b) => a + b, 0) / slopes.length / Vflat;
  }

  gals.push({
    name: g.name, logA0: g.logA0,
    Vflat, Rdisk, Rmax, Mbar, rc,
    logVflat: Math.log10(Vflat), logL36, logRdisk, logMbar,
    logMHI: g.logMHI, morphT: sp.T,
    logSBdisk, envCode: g.envCode,
    logK_halo, dmFrac_Rmax,
    haloResponse, outerSlope,
    gasFrac, newtDeficit,
    rcSmoothness, asymmetry, innerOuterCoherence, s1_Vflat,
    Npts: rc.length, D: sp.D, L36: sp.L36, MHI: sp.MHI,
  });
}

const N = gals.length;
const struct4 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
const struct6 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
const vfModel = multiR2(struct4, gals.map(g => g.logVflat));
const a0Model = multiR2(struct6, gals.map(g => g.logA0));
for (let i = 0; i < N; i++) { gals[i].VfResid = vfModel.residuals[i]; gals[i].a0Resid = a0Model.residuals[i]; }
const sdVf = Math.sqrt(gals.reduce((a, g) => a + g.VfResid ** 2, 0) / N);
const sdA0 = Math.sqrt(gals.reduce((a, g) => a + g.a0Resid ** 2, 0) / N);
for (let i = 0; i < N; i++) {
  gals[i].VfResid_z = gals[i].VfResid / sdVf;
  gals[i].a0Resid_z = gals[i].a0Resid / sdA0;
  gals[i].L_sum = gals[i].VfResid_z + gals[i].a0Resid_z;
}
const bestCtrl = gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.envCode]);
const Lresid = multiR2(bestCtrl, gals.map(g => g.L_sum));
for (let i = 0; i < N; i++) gals[i].dq = Lresid.residuals[i];

gals.sort((a, b) => b.dq - a.dq);

console.log('='.repeat(70));
console.log('PROGRAM 6A: MATCHED IFU DECISIVE MEASUREMENT');
console.log('Does H leave a distinct 2D/kinematic signature in matched pairs?');
console.log('='.repeat(70));
console.log('\nN = ' + N + ' galaxies with DQ computed');


console.log('\n\n' + '#'.repeat(70));
console.log('6A-1: TARGET SELECTION (High-H galaxies)');
console.log('#'.repeat(70));

const targetNames = ['NGC2841', 'NGC3741', 'ESO563-G021'];
const targets = [];
for (const tn of targetNames) {
  const g = gals.find(x => x.name === tn);
  if (g) { targets.push(g); }
  else { console.log('  WARNING: ' + tn + ' not in DQ sample'); }
}

console.log('\n  HIGH-H TARGETS:');
console.log('  ' + '-'.repeat(110));
console.log('  ' + 'Name'.padEnd(18) + 'DQ'.padEnd(8) + 'Rank'.padEnd(7) + 'Vflat'.padEnd(8) + 'logMbar'.padEnd(9) + 'Rdisk'.padEnd(7) + 'T'.padEnd(5) + 'hR'.padEnd(8) + 'VfR'.padEnd(9) + 'a0R'.padEnd(9) + 'smooth'.padEnd(8) + 's1/Vf'.padEnd(8) + 'D(Mpc)');
console.log('  ' + '-'.repeat(110));
for (const t of targets) {
  const rank = gals.indexOf(t) + 1;
  console.log('  ' + t.name.padEnd(18) + t.dq.toFixed(3).padEnd(8) + (rank + '/' + N).padEnd(7) + t.Vflat.toFixed(1).padEnd(8) + t.logMbar.toFixed(2).padEnd(9) + t.Rdisk.toFixed(2).padEnd(7) + t.morphT.toFixed(0).padEnd(5) + t.haloResponse.toFixed(3).padEnd(8) + t.VfResid.toFixed(4).padEnd(9) + t.a0Resid.toFixed(4).padEnd(9) + t.rcSmoothness.toFixed(3).padEnd(8) + t.s1_Vflat.toFixed(4).padEnd(8) + (t.D || 0).toFixed(1));
}


console.log('\n\n' + '#'.repeat(70));
console.log('6A-2: MATCHED CONTROL SELECTION');
console.log('For each target, find the best-matched LOW-DQ control');
console.log('#'.repeat(70));

function matchScore(target, candidate) {
  const dVflat = Math.abs(Math.log10(target.Vflat) - Math.log10(candidate.Vflat)) / 0.3;
  const dMbar = Math.abs(target.logMbar - candidate.logMbar) / 0.5;
  const dRdisk = Math.abs(Math.log10(target.Rdisk) - Math.log10(candidate.Rdisk)) / 0.3;
  const dT = Math.abs(target.morphT - candidate.morphT) / 3;
  return Math.sqrt(dVflat ** 2 + dMbar ** 2 + dRdisk ** 2 + dT ** 2);
}

const controls = [];
const usedNames = new Set(targetNames);

for (const target of targets) {
  console.log('\n  Matching for ' + target.name + ' (DQ=' + target.dq.toFixed(3) + ', Vflat=' + target.Vflat.toFixed(1) + ', logMbar=' + target.logMbar.toFixed(2) + ', T=' + target.morphT + '):');

  const lowDQ = gals.filter(g => g.dq < 0 && !usedNames.has(g.name));
  const scored = lowDQ.map(g => ({ g, score: matchScore(target, g) }));
  scored.sort((a, b) => a.score - b.score);

  console.log('  Top 5 matched low-DQ candidates:');
  console.log('  ' + 'Name'.padEnd(18) + 'DQ'.padEnd(8) + 'Vflat'.padEnd(8) + 'logMbar'.padEnd(9) + 'Rdisk'.padEnd(7) + 'T'.padEnd(5) + 'Match'.padEnd(8) + 'hR'.padEnd(8) + 'D(Mpc)');
  for (let i = 0; i < Math.min(5, scored.length); i++) {
    const c = scored[i].g;
    console.log('  ' + c.name.padEnd(18) + c.dq.toFixed(3).padEnd(8) + c.Vflat.toFixed(1).padEnd(8) + c.logMbar.toFixed(2).padEnd(9) + c.Rdisk.toFixed(2).padEnd(7) + c.morphT.toFixed(0).padEnd(5) + scored[i].score.toFixed(3).padEnd(8) + c.haloResponse.toFixed(3).padEnd(8) + (c.D || 0).toFixed(1));
  }

  const best = scored[0].g;
  controls.push(best);
  usedNames.add(best.name);
  console.log('  SELECTED: ' + best.name + ' (match score: ' + scored[0].score.toFixed(3) + ')');
}


console.log('\n\n' + '#'.repeat(70));
console.log('6A-3: MATCHED PAIR COMPARISON — 1D ROTATION CURVE DIAGNOSTICS');
console.log('#'.repeat(70));

const pairs = targets.map((t, i) => ({ target: t, control: controls[i] }));

for (const pair of pairs) {
  const t = pair.target, c = pair.control;
  console.log('\n  PAIR: ' + t.name + ' (H+) vs ' + c.name + ' (H-)');
  console.log('  ' + '-'.repeat(80));
  console.log('  Property'.padEnd(30) + t.name.padEnd(18) + c.name.padEnd(18) + 'Diff'.padEnd(10) + 'Expected');
  console.log('  ' + '-'.repeat(80));

  const metrics = [
    { name: 'DQ (H proxy)', tv: t.dq, cv: c.dq, expect: 'H+ >> H-' },
    { name: 'Vflat (km/s)', tv: t.Vflat, cv: c.Vflat, expect: 'matched' },
    { name: 'logMbar', tv: t.logMbar, cv: c.logMbar, expect: 'matched' },
    { name: 'Rdisk (kpc)', tv: t.Rdisk, cv: c.Rdisk, expect: 'matched' },
    { name: 'morphT', tv: t.morphT, cv: c.morphT, expect: 'matched' },
    { name: 'VfResid', tv: t.VfResid, cv: c.VfResid, expect: 'H+ > H-' },
    { name: 'a0Resid', tv: t.a0Resid, cv: c.a0Resid, expect: 'H+ > H-' },
    { name: 'haloResponse', tv: t.haloResponse, cv: c.haloResponse, expect: 'H+ > H-' },
    { name: 'rcSmoothness', tv: t.rcSmoothness, cv: c.rcSmoothness, expect: 'H+ > H-' },
    { name: 's1/Vflat', tv: t.s1_Vflat, cv: c.s1_Vflat, expect: 'H+ < H-' },
    { name: 'outerSlope', tv: t.outerSlope, cv: c.outerSlope, expect: 'H+ flatter' },
    { name: 'asymmetry', tv: t.asymmetry, cv: c.asymmetry, expect: 'H+ < H-' },
    { name: 'innerOuterCoherence', tv: t.innerOuterCoherence, cv: c.innerOuterCoherence, expect: 'H+ > H-' },
    { name: 'gasFrac', tv: t.gasFrac, cv: c.gasFrac, expect: 'variable' },
  ];

  for (const met of metrics) {
    const diff = met.tv - met.cv;
    const diffStr = (diff >= 0 ? '+' : '') + diff.toFixed(4);
    console.log('  ' + met.name.padEnd(30) + met.tv.toFixed(4).padEnd(18) + met.cv.toFixed(4).padEnd(18) + diffStr.padEnd(10) + met.expect);
  }
}


console.log('\n\n' + '#'.repeat(70));
console.log('6A-4: AGGREGATE DECISIVE TESTS');
console.log('Across ALL pairs, does H leave consistent signatures?');
console.log('#'.repeat(70));

const decisiveTests = [
  { name: 'T1: VfResid(H+) > VfResid(H-)', key: 'VfResid', dir: 1 },
  { name: 'T2: a0Resid(H+) > a0Resid(H-)', key: 'a0Resid', dir: 1 },
  { name: 'T3: haloResponse(H+) > haloResponse(H-)', key: 'haloResponse', dir: 1 },
  { name: 'T4: rcSmoothness(H+) > rcSmoothness(H-)', key: 'rcSmoothness', dir: 1 },
  { name: 'T5: s1/Vflat(H+) < s1/Vflat(H-)', key: 's1_Vflat', dir: -1 },
  { name: 'T6: asymmetry(H+) < asymmetry(H-)', key: 'asymmetry', dir: -1 },
  { name: 'T7: innerOuterCoherence(H+) > innerOuterCoherence(H-)', key: 'innerOuterCoherence', dir: 1 },
];

console.log('\n  DECISIVE TEST RESULTS:');
console.log('  ' + '-'.repeat(90));
console.log('  ' + 'Test'.padEnd(50) + 'Pair1'.padEnd(8) + 'Pair2'.padEnd(8) + 'Pair3'.padEnd(8) + 'Total'.padEnd(8) + 'PASS?');
console.log('  ' + '-'.repeat(90));

let totalPass = 0;
const testResults = [];
for (const test of decisiveTests) {
  let passCount = 0;
  const pairResults = [];
  for (const pair of pairs) {
    const tv = pair.target[test.key];
    const cv = pair.control[test.key];
    const diff = tv - cv;
    const pass = test.dir > 0 ? diff > 0 : diff < 0;
    if (pass) passCount++;
    pairResults.push(pass ? 'YES' : 'no');
  }
  const allPass = passCount >= 2;
  if (allPass) totalPass++;
  testResults.push({ name: test.name, passCount, total: 3, allPass });
  console.log('  ' + test.name.padEnd(50) + pairResults[0].padEnd(8) + pairResults[1].padEnd(8) + pairResults[2].padEnd(8) + (passCount + '/3').padEnd(8) + (allPass ? 'PASS' : 'FAIL'));
}

console.log('  ' + '-'.repeat(90));
console.log('  TOTAL DECISIVE TESTS PASSED: ' + totalPass + '/7');


console.log('\n\n' + '#'.repeat(70));
console.log('6A-5: EFFECT SIZES AND SIGNIFICANCE');
console.log('#'.repeat(70));

console.log('\n  MEAN DIFFERENCES (High-H minus Low-H, across pairs):');
const effectKeys = ['dq', 'VfResid', 'a0Resid', 'haloResponse', 'rcSmoothness', 's1_Vflat', 'asymmetry', 'innerOuterCoherence'];
for (const key of effectKeys) {
  const diffs = pairs.map(p => p.target[key] - p.control[key]);
  const mean = diffs.reduce((a, b) => a + b, 0) / diffs.length;
  const sd = Math.sqrt(diffs.reduce((a, d) => a + (d - mean) ** 2, 0) / Math.max(diffs.length - 1, 1));
  const t = sd > 0 ? mean / (sd / Math.sqrt(diffs.length)) : 0;
  console.log('  ' + key.padEnd(25) + 'mean diff = ' + (mean >= 0 ? '+' : '') + mean.toFixed(4) + '  SD = ' + sd.toFixed(4) + '  t = ' + (t >= 0 ? '+' : '') + t.toFixed(2));
}


console.log('\n\n' + '#'.repeat(70));
console.log('6A-6: DQ-RANKED FULL SAMPLE OVERVIEW');
console.log('#'.repeat(70));

console.log('\n  ALL ' + N + ' GALAXIES RANKED BY DQ:');
console.log('  ' + 'Rank'.padEnd(6) + 'Name'.padEnd(18) + 'DQ'.padEnd(8) + 'Vflat'.padEnd(8) + 'logMbar'.padEnd(9) + 'T'.padEnd(5) + 'hR'.padEnd(8) + 'VfR'.padEnd(9) + 'a0R'.padEnd(9) + 'smooth'.padEnd(8) + 'Role');
console.log('  ' + '-'.repeat(100));
for (let i = 0; i < N; i++) {
  const g = gals[i];
  let role = '';
  if (targetNames.includes(g.name)) role = '*** TARGET ***';
  else if (controls.find(c => c.name === g.name)) role = '--- control ---';
  console.log('  ' + (i + 1 + '').padEnd(6) + g.name.padEnd(18) + g.dq.toFixed(3).padEnd(8) + g.Vflat.toFixed(1).padEnd(8) + g.logMbar.toFixed(2).padEnd(9) + g.morphT.toFixed(0).padEnd(5) + g.haloResponse.toFixed(3).padEnd(8) + g.VfResid.toFixed(4).padEnd(9) + g.a0Resid.toFixed(4).padEnd(9) + g.rcSmoothness.toFixed(3).padEnd(8) + role);
}


console.log('\n\n' + '#'.repeat(70));
console.log('6A-7: IFU OBSERVING PROPOSAL SPECIFICATIONS');
console.log('#'.repeat(70));

console.log('\n  OBSERVING TARGETS AND CONTROLS:');
console.log('  ' + '-'.repeat(90));
console.log('  ' + 'Galaxy'.padEnd(18) + 'Role'.padEnd(12) + 'RA/DEC'.padEnd(18) + 'D(Mpc)'.padEnd(10) + 'Vflat'.padEnd(10) + 'DQ'.padEnd(8) + 'logMbar');
console.log('  ' + '-'.repeat(90));
for (const pair of pairs) {
  console.log('  ' + pair.target.name.padEnd(18) + 'TARGET'.padEnd(12) + '(from NED)'.padEnd(18) + (pair.target.D || 0).toFixed(1).padEnd(10) + pair.target.Vflat.toFixed(1).padEnd(10) + pair.target.dq.toFixed(3).padEnd(8) + pair.target.logMbar.toFixed(2));
  console.log('  ' + pair.control.name.padEnd(18) + 'CONTROL'.padEnd(12) + '(from NED)'.padEnd(18) + (pair.control.D || 0).toFixed(1).padEnd(10) + pair.control.Vflat.toFixed(1).padEnd(10) + pair.control.dq.toFixed(3).padEnd(8) + pair.control.logMbar.toFixed(2));
  console.log('');
}

console.log('\n  REQUIRED MEASUREMENTS PER GALAXY:');
console.log('  1. s1/Vflat — velocity dispersion normalised by flat velocity');
console.log('       → from IFU: sigma_los across FOV / Vrot');
console.log('  2. PA twist — position angle variation with radius');
console.log('       → from IFU: kinematic PA profile from tilted-ring fits');
console.log('  3. Lopsidedness / asymmetry map');
console.log('       → from IFU: A1/A0 Fourier decomposition of velocity field');
console.log('  4. 2D residual symmetry');
console.log('       → from IFU: residuals after tilted-ring model subtraction');
console.log('  5. Inner–outer kinematic coherence');
console.log('       → from IFU: compare kinematic PA and inclination inner vs outer');

console.log('\n  DECISIVE PREDICTIONS:');
console.log('  For EACH matched pair (high-H target vs low-H control):');
console.log('  P1: target shows LOWER s1/Vflat (quieter kinematics)');
console.log('  P2: target shows SMALLER PA twist (more coherent rotation axis)');
console.log('  P3: target shows LOWER lopsidedness A1/A0 (more symmetric velocity field)');
console.log('  P4: target shows SMALLER 2D residuals (smoother velocity field)');
console.log('  P5: target shows HIGHER inner–outer coherence');
console.log('  P6: target has HIGHER haloResponse (more DM-supported)');
console.log('  P7: target has HIGHER VfResid AND a0Resid (bilateral channel active)');

console.log('\n  PASS/FAIL CRITERION:');
console.log('  Decisive confirmation: >= 5/7 predictions correct in >= 2/3 pairs');
console.log('  Strong support:        >= 4/7 predictions correct in >= 2/3 pairs');
console.log('  Inconclusive:          3/7 or ambiguous pair results');
console.log('  Disconfirmation:       <= 2/7 correct, or controls outperform targets');


console.log('\n\n' + '='.repeat(70));
console.log('PROGRAM 6A GRAND VERDICT (1D ROTATION-CURVE PREVIEW)');
console.log('='.repeat(70));

console.log('\n  Using available 1D rotation curve data as PREVIEW of what IFU would show:');
console.log('  Decisive tests passed: ' + totalPass + '/7');
const verdict6a = totalPass >= 5 ? 'STRONG SUPPORT' : totalPass >= 4 ? 'MODERATE SUPPORT' : totalPass >= 3 ? 'INCONCLUSIVE' : 'DISCONFIRMATION';
console.log('  PREVIEW VERDICT: ' + verdict6a);

if (totalPass >= 4) {
  console.log('\n  The 1D rotation-curve preview SUPPORTS the H hypothesis.');
  console.log('  An IFU follow-up with 2D kinematic measurements would provide');
  console.log('  the decisive test at higher fidelity.');
}

console.log('\n  MATCH QUALITY:');
for (let i = 0; i < pairs.length; i++) {
  const p = pairs[i];
  const dVflat = Math.abs(p.target.Vflat - p.control.Vflat);
  const dLogMbar = Math.abs(p.target.logMbar - p.control.logMbar);
  const dT = Math.abs(p.target.morphT - p.control.morphT);
  const dDQ = p.target.dq - p.control.dq;
  console.log('  Pair ' + (i + 1) + ': ' + p.target.name + ' vs ' + p.control.name);
  console.log('    DQ gap = ' + dDQ.toFixed(3) + ', dVflat = ' + dVflat.toFixed(1) + ' km/s, dlogMbar = ' + dLogMbar.toFixed(2) + ', dT = ' + dT.toFixed(0));
}


const outPath = path.join(__dirname, '..', 'public', 'program6a-matched-ifu.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: '6A',
  title: 'Matched IFU Decisive Measurement',
  timestamp: new Date().toISOString(),
  N,
  targets: targets.map(t => ({ name: t.name, dq: t.dq, Vflat: t.Vflat, logMbar: t.logMbar, morphT: t.morphT, haloResponse: t.haloResponse, VfResid: t.VfResid, a0Resid: t.a0Resid, rcSmoothness: t.rcSmoothness, s1_Vflat: t.s1_Vflat, D: t.D })),
  controls: controls.map(c => ({ name: c.name, dq: c.dq, Vflat: c.Vflat, logMbar: c.logMbar, morphT: c.morphT, haloResponse: c.haloResponse, VfResid: c.VfResid, a0Resid: c.a0Resid, rcSmoothness: c.rcSmoothness, s1_Vflat: c.s1_Vflat, D: c.D })),
  pairs: pairs.map(p => ({ target: p.target.name, control: p.control.name })),
  decisiveTests: testResults,
  totalPass: totalPass,
  verdict: verdict6a,
  goldenSentence: 'H is not a quietness variable. H is a hidden common-cause state that jointly drives both the rotation-velocity residual and the acceleration residual; kinematic calmness emerges as a consequence, not as the mechanism.',
}, null, 2));
console.log('\nSaved: ' + outPath);

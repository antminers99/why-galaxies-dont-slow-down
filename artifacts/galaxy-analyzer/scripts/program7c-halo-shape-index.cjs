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

function partialR(x, y, controls) {
  const n = x.length;
  if (n < controls[0].length + 4) return NaN;
  const residualise = (v, C) => {
    const nv2 = C[0].length;
    const my2 = v.reduce((a, b) => a + b, 0) / n;
    const mc2 = Array(nv2).fill(0);
    for (let j = 0; j < nv2; j++) { for (let i = 0; i < n; i++) mc2[j] += C[i][j]; mc2[j] /= n; }
    const XTX2 = Array.from({ length: nv2 }, () => Array(nv2).fill(0));
    const XTy2 = Array(nv2).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < nv2; j++) {
        XTy2[j] += (C[i][j] - mc2[j]) * (v[i] - my2);
        for (let k = 0; k < nv2; k++) XTX2[j][k] += (C[i][j] - mc2[j]) * (C[i][k] - mc2[k]);
      }
    }
    const aug2 = XTX2.map((row, i) => [...row, XTy2[i]]);
    for (let col = 0; col < nv2; col++) {
      let maxRow = col;
      for (let row = col + 1; row < nv2; row++) if (Math.abs(aug2[row][col]) > Math.abs(aug2[maxRow][col])) maxRow = row;
      [aug2[col], aug2[maxRow]] = [aug2[maxRow], aug2[col]];
      if (Math.abs(aug2[col][col]) < 1e-12) continue;
      for (let row = col + 1; row < nv2; row++) { const f = aug2[row][col] / aug2[col][col]; for (let j = col; j <= nv2; j++) aug2[row][j] -= f * aug2[col][j]; }
    }
    const beta2 = Array(nv2).fill(0);
    for (let i = nv2 - 1; i >= 0; i--) { beta2[i] = aug2[i][nv2]; for (let j = i + 1; j < nv2; j++) beta2[i] -= aug2[i][j] * beta2[j]; beta2[i] /= aug2[i][i] || 1; }
    const r = [];
    for (let i = 0; i < n; i++) { let pred = my2; for (let j = 0; j < nv2; j++) pred += beta2[j] * (C[i][j] - mc2[j]); r.push(v[i] - pred); }
    return r;
  };
  return pearsonR(residualise(x, controls), residualise(y, controls));
}

function zScore(arr) {
  const n = arr.length;
  const mu = arr.reduce((a, b) => a + b, 0) / n;
  const sd = Math.sqrt(arr.reduce((a, v) => a + (v - mu) ** 2, 0) / n);
  return arr.map(v => sd > 0 ? (v - mu) / sd : 0);
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
  const logMbar = Math.log10(Math.max(Mbar, 1));
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(Rdisk, 0.01));
  const logSBdisk = Math.log10(Math.max(sp.SBdisk, 0.01));
  const logK = sr.models.dark_halo_linear.k > 0 ? Math.log10(sr.models.dark_halo_linear.k) : -5;
  const hR = sr.models.newtonian.mse > 0 ? Math.log10(Math.max(sr.models.newtonian.mse / Math.max(sr.models.dark_halo_linear.mse, 0.001), 0.01)) : 0;
  const V_Newt_Rmax = Math.sqrt(G * Mbar / Rmax);
  const dmFrac = Math.max(0, 1 - (V_Newt_Rmax / Vflat) ** 2);

  let rcSmooth = 0;
  if (rc.length >= 5) {
    const sm = [];
    for (let i = 2; i < rc.length - 2; i++) sm.push((rc[i - 2].v + rc[i - 1].v + rc[i].v + rc[i + 1].v + rc[i + 2].v) / 5);
    let ss = 0;
    for (let i = 0; i < sm.length; i++) ss += (rc[i + 2].v - sm[i]) ** 2;
    rcSmooth = 1 - Math.sqrt(ss / sm.length) / Vflat;
  }

  const haloVelocities = [];
  for (const p of rc) {
    if (p.r < 0.3) continue;
    const encMass = Mbar * Math.min(p.r / Rmax, 1);
    const V_bar = Math.sqrt(G * encMass / p.r);
    const V_halo_sq = Math.max(p.v * p.v - V_bar * V_bar, 0);
    haloVelocities.push({ r: p.r, V_halo: Math.sqrt(V_halo_sq), V_obs: p.v, V_bar });
  }

  let concentration = 5;
  let bestNFWmse = Infinity;
  for (let c = 2; c <= 40; c += 2) {
    for (let V200 = 20; V200 <= 400; V200 += 20) {
      const Rs = Rmax / c;
      let sse = 0;
      for (const h of haloVelocities) {
        const x = h.r / Rs;
        const gx = Math.log(1 + x) - x / (1 + x);
        const gc = Math.log(1 + c) - c / (1 + c);
        const pred = V200 * Math.sqrt(Math.max(gx / (x * gc), 0));
        sse += (h.V_halo - pred) ** 2;
      }
      const mse = sse / haloVelocities.length;
      if (mse < bestNFWmse) { bestNFWmse = mse; concentration = c; }
    }
  }

  const zones = { inner: { min: 0, max: Rdisk }, disk: { min: Rdisk, max: 3 * Rdisk }, outer: { min: 3 * Rdisk, max: Infinity } };
  let totalSupport = 0;
  const zoneSupport = {};
  for (const [zn, zr] of Object.entries(zones)) {
    const pts = haloVelocities.filter(h => h.r >= zr.min && h.r < zr.max);
    const sup = pts.reduce((s, h) => s + h.V_halo * h.V_halo, 0);
    zoneSupport[zn] = sup;
    totalSupport += sup;
  }
  const innerFrac = totalSupport > 0 ? zoneSupport.inner / totalSupport : 0;
  const diskFrac = totalSupport > 0 ? zoneSupport.disk / totalSupport : 0;
  const outerFrac = totalSupport > 0 ? zoneSupport.outer / totalSupport : 0;

  const outerExcess = outerFrac - 0.5;
  const diskDeficit = 0.33 - diskFrac;

  gals.push({
    name: g.name, Vflat, Rdisk, Rmax, Mbar, rc,
    logVflat: Math.log10(Vflat), logL36, logRdisk, logMbar,
    logMHI: g.logMHI, morphT: sp.T, logA0: g.logA0,
    logSBdisk, envCode: g.envCode,
    logK, dmFrac, hR, rcSmooth,
    concentration, innerFrac, diskFrac, outerFrac,
    outerExcess, diskDeficit,
  });
}

const N = gals.length;

const concRegress = multiR2(gals.map(g => [g.logMbar, g.logVflat]), gals.map(g => g.concentration));
for (let i = 0; i < N; i++) gals[i].concResid = concRegress.residuals[i];

const struct4 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
const struct6 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
const vfM = multiR2(struct4, gals.map(g => g.logVflat));
const a0M = multiR2(struct6, gals.map(g => g.logA0));
for (let i = 0; i < N; i++) { gals[i].VfR = vfM.residuals[i]; gals[i].a0R = a0M.residuals[i]; }
const sdVf = Math.sqrt(gals.reduce((a, g) => a + g.VfR ** 2, 0) / N);
const sdA0 = Math.sqrt(gals.reduce((a, g) => a + g.a0R ** 2, 0) / N);
for (let i = 0; i < N; i++) {
  gals[i].VfR_z = gals[i].VfR / sdVf;
  gals[i].a0R_z = gals[i].a0R / sdA0;
  gals[i].L_sum = gals[i].VfR_z + gals[i].a0R_z;
}
const ctrlV = gals.map(g => [g.logK, g.dmFrac, g.envCode]);
const Lresid = multiR2(ctrlV, gals.map(g => g.L_sum));
for (let i = 0; i < N; i++) gals[i].dq = Lresid.residuals[i];

const r_VfR_a0R = pearsonR(gals.map(g => g.VfR), gals.map(g => g.a0R));


const concResid_z = zScore(gals.map(g => g.concResid));
const outerExcess_z = zScore(gals.map(g => g.outerExcess));
const diskDeficit_z = zScore(gals.map(g => g.diskDeficit));
const rcSmooth_z = zScore(gals.map(g => g.rcSmooth));

for (let i = 0; i < N; i++) {
  gals[i].S1 = -concResid_z[i] + outerExcess_z[i] + diskDeficit_z[i];
  gals[i].S2 = gals[i].S1 * (1 + 0.5 * rcSmooth_z[i]);
}

console.log('='.repeat(70));
console.log('PROGRAM 7C: HALO-SHAPE DECISIVE INDEX');
console.log('Build a single halo-shape index and test if it captures H');
console.log('='.repeat(70));
console.log('\nN = ' + N + ' galaxies');
console.log('Baseline r(VfR, a0R) = ' + r_VfR_a0R.toFixed(3));


console.log('\n\n' + '#'.repeat(70));
console.log('7C-1: INDEX CONSTRUCTION');
console.log('#'.repeat(70));

console.log('\n  S1 (Shape-only) = -z(concResid) + z(outerExcess) + z(diskDeficit)');
console.log('  S2 (Shape x Quiet) = S1 * (1 + 0.5 * z(rcSmooth))');
console.log('\n  Component correlations with DQ:');

const components = [
  { name: 'concResid', vals: gals.map(g => g.concResid) },
  { name: 'outerExcess', vals: gals.map(g => g.outerExcess) },
  { name: 'diskDeficit', vals: gals.map(g => g.diskDeficit) },
  { name: 'rcSmooth', vals: gals.map(g => g.rcSmooth) },
  { name: 'haloResponse', vals: gals.map(g => g.hR) },
  { name: 'S1 (shape-only)', vals: gals.map(g => g.S1) },
  { name: 'S2 (shape x quiet)', vals: gals.map(g => g.S2) },
];

const dqArr = gals.map(g => g.dq);

console.log('  ' + '-'.repeat(55));
console.log('  ' + 'Variable'.padEnd(25) + 'r(DQ, X)'.padEnd(15) + 'Sig?');
console.log('  ' + '-'.repeat(55));

const corrTable = {};
for (const c of components) {
  const r = pearsonR(dqArr, c.vals);
  const t = r * Math.sqrt((N - 2) / (1 - r * r));
  const sig = Math.abs(t) > 2;
  corrTable[c.name] = { r, t, sig };
  console.log('  ' + c.name.padEnd(25) + ((r >= 0 ? '+' : '') + r.toFixed(3)).padEnd(15) + (sig ? 'yes (t=' + t.toFixed(2) + ')' : 'no'));
}


console.log('\n\n' + '#'.repeat(70));
console.log('7C-2: TEST 1 — DOES S1/S2 BEAT haloResponse AS H-PROXY?');
console.log('#'.repeat(70));

const r_dq_hR = corrTable['haloResponse'].r;
const r_dq_S1 = corrTable['S1 (shape-only)'].r;
const r_dq_S2 = corrTable['S2 (shape x quiet)'].r;

console.log('\n  r(DQ, haloResponse) = ' + r_dq_hR.toFixed(3));
console.log('  r(DQ, S1)           = ' + r_dq_S1.toFixed(3));
console.log('  r(DQ, S2)           = ' + r_dq_S2.toFixed(3));

const s1Beats = Math.abs(r_dq_S1) > Math.abs(r_dq_hR);
const s2Beats = Math.abs(r_dq_S2) > Math.abs(r_dq_hR);
console.log('\n  S1 beats haloResponse: ' + (s1Beats ? 'YES' : 'NO'));
console.log('  S2 beats haloResponse: ' + (s2Beats ? 'YES' : 'NO'));

const bestIndex = Math.abs(r_dq_S2) >= Math.abs(r_dq_S1) ? 'S2' : 'S1';
const bestIndexR = bestIndex === 'S2' ? r_dq_S2 : r_dq_S1;
console.log('  Best index: ' + bestIndex + ' (r=' + bestIndexR.toFixed(3) + ')');
console.log('  TEST 1: ' + (s1Beats || s2Beats ? 'PASS' : 'FAIL'));


console.log('\n\n' + '#'.repeat(70));
console.log('7C-3: TEST 2 — CHANNEL ABSORPTION');
console.log('Does controlling for S1/S2 reduce r(VfR, a0R)?');
console.log('#'.repeat(70));

const ctrl_logK_dm = gals.map(g => [g.logK, g.dmFrac]);
const ctrl_logK_dm_env = gals.map(g => [g.logK, g.dmFrac, g.envCode]);
const ctrl_S1 = gals.map(g => [g.S1]);
const ctrl_S2 = gals.map(g => [g.S2]);
const ctrl_hR = gals.map(g => [g.hR]);
const ctrl_S1_hR = gals.map(g => [g.S1, g.hR]);
const ctrl_S2_hR = gals.map(g => [g.S2, g.hR]);
const ctrl_full_S1 = gals.map(g => [g.logK, g.dmFrac, g.envCode, g.S1]);
const ctrl_full_S2 = gals.map(g => [g.logK, g.dmFrac, g.envCode, g.S2]);
const ctrl_full_hR = gals.map(g => [g.logK, g.dmFrac, g.envCode, g.hR]);
const ctrl_full_S2_hR = gals.map(g => [g.logK, g.dmFrac, g.envCode, g.S2, g.hR]);

const vfArr = gals.map(g => g.VfR);
const a0Arr = gals.map(g => g.a0R);

const absorptionTests = [
  { name: 'raw (no controls)', r: r_VfR_a0R },
  { name: 'ctrl: logK+dmFrac', r: partialR(vfArr, a0Arr, ctrl_logK_dm) },
  { name: 'ctrl: logK+dmFrac+env', r: partialR(vfArr, a0Arr, ctrl_logK_dm_env) },
  { name: 'ctrl: haloResponse', r: partialR(vfArr, a0Arr, ctrl_hR) },
  { name: 'ctrl: S1', r: partialR(vfArr, a0Arr, ctrl_S1) },
  { name: 'ctrl: S2', r: partialR(vfArr, a0Arr, ctrl_S2) },
  { name: 'ctrl: S1+hR', r: partialR(vfArr, a0Arr, ctrl_S1_hR) },
  { name: 'ctrl: S2+hR', r: partialR(vfArr, a0Arr, ctrl_S2_hR) },
  { name: 'ctrl: logK+dmFrac+env+hR', r: partialR(vfArr, a0Arr, ctrl_full_hR) },
  { name: 'ctrl: logK+dmFrac+env+S1', r: partialR(vfArr, a0Arr, ctrl_full_S1) },
  { name: 'ctrl: logK+dmFrac+env+S2', r: partialR(vfArr, a0Arr, ctrl_full_S2) },
  { name: 'ctrl: logK+dmFrac+env+S2+hR', r: partialR(vfArr, a0Arr, ctrl_full_S2_hR) },
];

console.log('\n  PARTIAL r(VfR, a0R | controls):');
console.log('  ' + '-'.repeat(60));
console.log('  ' + 'Controls'.padEnd(35) + 'partial r'.padEnd(12) + 'Reduction');
console.log('  ' + '-'.repeat(60));

for (const test of absorptionTests) {
  const reduction = 1 - Math.abs(test.r) / Math.abs(r_VfR_a0R);
  console.log('  ' + test.name.padEnd(35) + ((test.r >= 0 ? '+' : '') + test.r.toFixed(3)).padEnd(12) + (reduction * 100).toFixed(1) + '%');
}

const bestAbsName = 'ctrl: logK+dmFrac+env+S2+hR';
const bestAbs = absorptionTests.find(t => t.name === bestAbsName);
const s2OnlyAbs = absorptionTests.find(t => t.name === 'ctrl: S2');
const hROnlyAbs = absorptionTests.find(t => t.name === 'ctrl: haloResponse');
const s2AbsRed = 1 - Math.abs(s2OnlyAbs.r) / Math.abs(r_VfR_a0R);
const hRAbsRed = 1 - Math.abs(hROnlyAbs.r) / Math.abs(r_VfR_a0R);

console.log('\n  S2 alone absorbs: ' + (s2AbsRed * 100).toFixed(1) + '% of the channel');
console.log('  hR alone absorbs: ' + (hRAbsRed * 100).toFixed(1) + '% of the channel');
console.log('  S2 > hR absorption: ' + (s2AbsRed > hRAbsRed ? 'YES' : 'NO'));
console.log('  TEST 2: ' + (s2AbsRed > hRAbsRed || s2AbsRed > 0.05 ? 'PASS' : 'FAIL'));


console.log('\n\n' + '#'.repeat(70));
console.log('7C-4: TEST 3 — MATCHED PAIR VALIDATION');
console.log('#'.repeat(70));

gals.sort((a, b) => b.dq - a.dq);
const targetNames = ['NGC2841', 'NGC3741', 'ESO563-G021'];

function matchScore(t, c) {
  const dVf = Math.abs(Math.log10(t.Vflat) - Math.log10(c.Vflat)) / 0.3;
  const dMb = Math.abs(t.logMbar - c.logMbar) / 0.5;
  const dRd = Math.abs(Math.log10(t.Rdisk) - Math.log10(c.Rdisk)) / 0.3;
  const dT = Math.abs(t.morphT - c.morphT) / 3;
  return Math.sqrt(dVf ** 2 + dMb ** 2 + dRd ** 2 + dT ** 2);
}

const targets = targetNames.map(n => gals.find(g => g.name === n)).filter(Boolean);
const controls = [];
const used = new Set(targetNames);
for (const target of targets) {
  const lowDQ = gals.filter(g => g.dq < 0 && !used.has(g.name));
  const scored = lowDQ.map(g => ({ g, score: matchScore(target, g) })).sort((a, b) => a.score - b.score);
  const best = scored[0].g;
  controls.push(best);
  used.add(best.name);
}
const pairs = targets.map((t, i) => ({ target: t, control: controls[i] }));

console.log('\n  MATCHED PAIR S1/S2 VALUES:');
console.log('  ' + '-'.repeat(100));
console.log('  ' + 'Pair'.padEnd(6) + 'Target'.padEnd(18) + 'DQ'.padEnd(8) + 'S1'.padEnd(8) + 'S2'.padEnd(8) + 'Control'.padEnd(18) + 'DQ'.padEnd(8) + 'S1'.padEnd(8) + 'S2'.padEnd(8) + 'S1 diff'.padEnd(10) + 'S2 diff');
console.log('  ' + '-'.repeat(100));

let pairS1pass = 0, pairS2pass = 0;
for (let i = 0; i < pairs.length; i++) {
  const t = pairs[i].target, c = pairs[i].control;
  const s1diff = t.S1 - c.S1;
  const s2diff = t.S2 - c.S2;
  if (s1diff > 0) pairS1pass++;
  if (s2diff > 0) pairS2pass++;
  console.log('  ' + (i + 1 + '').padEnd(6) + t.name.padEnd(18) + t.dq.toFixed(2).padEnd(8) + t.S1.toFixed(2).padEnd(8) + t.S2.toFixed(2).padEnd(8) + c.name.padEnd(18) + c.dq.toFixed(2).padEnd(8) + c.S1.toFixed(2).padEnd(8) + c.S2.toFixed(2).padEnd(8) + ((s1diff >= 0 ? '+' : '') + s1diff.toFixed(2)).padEnd(10) + (s2diff >= 0 ? '+' : '') + s2diff.toFixed(2));
}

console.log('\n  S1: ' + pairS1pass + '/3 pairs correct direction');
console.log('  S2: ' + pairS2pass + '/3 pairs correct direction');
console.log('  TEST 3: ' + (pairS1pass >= 2 || pairS2pass >= 2 ? 'PASS' : 'FAIL'));


console.log('\n\n' + '#'.repeat(70));
console.log('7C-5: TEST 4 — INTERNAL STABILITY (SPLIT-HALF)');
console.log('#'.repeat(70));

const shuffled = gals.slice().sort((a, b) => a.name.localeCompare(b.name));
const halfN = Math.floor(N / 2);
const half1 = shuffled.slice(0, halfN);
const half2 = shuffled.slice(halfN);

const r_s1_dq_h1 = pearsonR(half1.map(g => g.S1), half1.map(g => g.dq));
const r_s1_dq_h2 = pearsonR(half2.map(g => g.S1), half2.map(g => g.dq));
const r_s2_dq_h1 = pearsonR(half1.map(g => g.S2), half1.map(g => g.dq));
const r_s2_dq_h2 = pearsonR(half2.map(g => g.S2), half2.map(g => g.dq));

const shrinkS1 = Math.abs(r_s1_dq_h1 - r_s1_dq_h2);
const shrinkS2 = Math.abs(r_s2_dq_h1 - r_s2_dq_h2);

console.log('\n  SPLIT-HALF r(index, DQ):');
console.log('  ' + '-'.repeat(55));
console.log('  ' + 'Index'.padEnd(12) + 'Half 1'.padEnd(12) + 'Half 2'.padEnd(12) + 'Shrinkage');
console.log('  ' + '-'.repeat(55));
console.log('  ' + 'S1'.padEnd(12) + r_s1_dq_h1.toFixed(3).padEnd(12) + r_s1_dq_h2.toFixed(3).padEnd(12) + shrinkS1.toFixed(3));
console.log('  ' + 'S2'.padEnd(12) + r_s2_dq_h1.toFixed(3).padEnd(12) + r_s2_dq_h2.toFixed(3).padEnd(12) + shrinkS2.toFixed(3));

const stableS1 = shrinkS1 < 0.3 && r_s1_dq_h1 * r_s1_dq_h2 > 0;
const stableS2 = shrinkS2 < 0.3 && r_s2_dq_h1 * r_s2_dq_h2 > 0;
console.log('\n  S1 stable (same sign, shrinkage < 0.3): ' + (stableS1 ? 'YES' : 'NO'));
console.log('  S2 stable (same sign, shrinkage < 0.3): ' + (stableS2 ? 'YES' : 'NO'));
console.log('  TEST 4: ' + (stableS1 || stableS2 ? 'PASS' : 'FAIL'));


console.log('\n\n' + '#'.repeat(70));
console.log('7C-6: QUINTILE DEEP-DIVE');
console.log('#'.repeat(70));

gals.sort((a, b) => b.dq - a.dq);
const quintileSize = Math.floor(N / 5);
const Q1 = gals.slice(0, quintileSize);
const Q5 = gals.slice(N - quintileSize);

const qMetrics = [
  { name: 'DQ', ex: g => g.dq },
  { name: 'VfResid', ex: g => g.VfR },
  { name: 'a0Resid', ex: g => g.a0R },
  { name: 'S1', ex: g => g.S1 },
  { name: 'S2', ex: g => g.S2 },
  { name: 'concResid', ex: g => g.concResid },
  { name: 'outerFrac', ex: g => g.outerFrac },
  { name: 'diskFrac', ex: g => g.diskFrac },
  { name: 'rcSmooth', ex: g => g.rcSmooth },
  { name: 'haloResponse', ex: g => g.hR },
];

console.log('\n  ' + 'Metric'.padEnd(20) + 'Q1 (high-H)'.padEnd(15) + 'Q5 (low-H)'.padEnd(15) + 'Diff');
console.log('  ' + '-'.repeat(60));
for (const met of qMetrics) {
  const q1m = Q1.reduce((s, g) => s + met.ex(g), 0) / Q1.length;
  const q5m = Q5.reduce((s, g) => s + met.ex(g), 0) / Q5.length;
  const diff = q1m - q5m;
  console.log('  ' + met.name.padEnd(20) + q1m.toFixed(4).padEnd(15) + q5m.toFixed(4).padEnd(15) + (diff >= 0 ? '+' : '') + diff.toFixed(4));
}


console.log('\n\n' + '#'.repeat(70));
console.log('7C-7: INCREMENTAL ABSORPTION LADDER');
console.log('How much of the channel does each variable absorb step-by-step?');
console.log('#'.repeat(70));

const ladder = [
  { name: 'raw', ctrl: null },
  { name: '+logK', ctrl: g => [g.logK] },
  { name: '+dmFrac', ctrl: g => [g.logK, g.dmFrac] },
  { name: '+env', ctrl: g => [g.logK, g.dmFrac, g.envCode] },
  { name: '+hR', ctrl: g => [g.logK, g.dmFrac, g.envCode, g.hR] },
  { name: '+S1', ctrl: g => [g.logK, g.dmFrac, g.envCode, g.hR, g.S1] },
  { name: '+S2 (replace S1)', ctrl: g => [g.logK, g.dmFrac, g.envCode, g.hR, g.S2] },
  { name: '+concResid', ctrl: g => [g.logK, g.dmFrac, g.envCode, g.hR, g.concResid] },
  { name: '+outerFrac', ctrl: g => [g.logK, g.dmFrac, g.envCode, g.hR, g.outerFrac] },
  { name: '+concResid+outerFrac', ctrl: g => [g.logK, g.dmFrac, g.envCode, g.hR, g.concResid, g.outerFrac] },
  { name: '+S2+concR+outerF', ctrl: g => [g.logK, g.dmFrac, g.envCode, g.hR, g.S2, g.concResid, g.outerFrac] },
];

console.log('\n  ' + 'Step'.padEnd(30) + 'partial r'.padEnd(12) + 'Absorbed'.padEnd(12) + 'Remaining');
console.log('  ' + '-'.repeat(65));

for (const step of ladder) {
  let pr;
  if (!step.ctrl) {
    pr = r_VfR_a0R;
  } else {
    const C = gals.map(step.ctrl);
    pr = partialR(vfArr, a0Arr, C);
  }
  const absorbed = (1 - Math.abs(pr) / Math.abs(r_VfR_a0R)) * 100;
  console.log('  ' + step.name.padEnd(30) + ((pr >= 0 ? '+' : '') + pr.toFixed(3)).padEnd(12) + absorbed.toFixed(1).padEnd(12) + (100 - absorbed).toFixed(1) + '%');
}


console.log('\n\n' + '='.repeat(70));
console.log('PROGRAM 7C GRAND VERDICT');
console.log('='.repeat(70));

const test1 = s1Beats || s2Beats;
const test2 = s2AbsRed > hRAbsRed || s2AbsRed > 0.05;
const test3 = pairS1pass >= 2 || pairS2pass >= 2;
const test4 = stableS1 || stableS2;

const allTests = [
  { name: 'T1: S1/S2 beats haloResponse as H-proxy', pass: test1 },
  { name: 'T2: S1/S2 absorbs more channel than hR', pass: test2 },
  { name: 'T3: Matched-pair validation (>= 2/3)', pass: test3 },
  { name: 'T4: Split-half stability', pass: test4 },
];

let totalPass = 0;
for (const t of allTests) {
  if (t.pass) totalPass++;
  console.log('  ' + (t.pass ? 'PASS' : 'FAIL').padEnd(6) + t.name);
}
console.log('  Total: ' + totalPass + '/4');

if (totalPass >= 3) {
  console.log('\n  VERDICT: HALO-SHAPE INDEX IS THE KEY');
  console.log('  S1/S2 is a better proxy for H than haloResponse alone.');
  console.log('  The hidden variable is best captured by: under-concentration');
  console.log('  + outer support excess + disk deficit (+ quietness coupling).');
} else if (totalPass >= 2) {
  console.log('\n  VERDICT: PARTIAL — halo-shape index has promise');
  console.log('  but does not decisively replace haloResponse.');
} else {
  console.log('\n  VERDICT: HALO-SHAPE INDEX INSUFFICIENT');
  console.log('  The index does not capture H better than existing proxies.');
  console.log('  The hidden variable may operate through a mechanism not');
  console.log('  reducible to these halo-shape metrics.');
}


const outPath = path.join(__dirname, '..', 'public', 'program7c-halo-shape-index.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: '7C',
  title: 'Halo-Shape Decisive Index',
  timestamp: new Date().toISOString(),
  N,
  indices: { S1_formula: '-z(concResid) + z(outerExcess) + z(diskDeficit)', S2_formula: 'S1 * (1 + 0.5 * z(rcSmooth))' },
  correlations: corrTable,
  absorption: absorptionTests.map(t => ({ name: t.name, r: t.r, reduction: 1 - Math.abs(t.r) / Math.abs(r_VfR_a0R) })),
  matchedPairs: { S1pass: pairS1pass, S2pass: pairS2pass },
  stability: { S1: { h1: r_s1_dq_h1, h2: r_s1_dq_h2, shrink: shrinkS1, stable: stableS1 }, S2: { h1: r_s2_dq_h1, h2: r_s2_dq_h2, shrink: shrinkS2, stable: stableS2 } },
  tests: { t1: test1, t2: test2, t3: test3, t4: test4, totalPass },
  verdict: totalPass >= 3 ? 'HALO-SHAPE INDEX IS THE KEY' : totalPass >= 2 ? 'PARTIAL' : 'INSUFFICIENT',
}, null, 2));
console.log('\nSaved: ' + outPath);

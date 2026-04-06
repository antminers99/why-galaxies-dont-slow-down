#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 133A: REGIME LAW — High-Vflat Dedicated Model');
console.log('');
console.log('  Does the coupling law become cleaner/stronger when restricted');
console.log('  to the massive galaxy regime where VfResid is strongest?');
console.log('  Is there a sharp regime boundary or a smooth transition?');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'stage-A-master-table.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const sparcTable = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-results.json'), 'utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparcTable.forEach(s => { sparcMap[s.name] = s; });
const resMap = {}; sparcResults.perGalaxy.forEach(g => { if (g.rawName) resMap[g.rawName] = g; });
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
const N = gals45.length;

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function pearsonR(x, y) {
  const n = x.length, mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}
function solveLinear(A, b) {
  const n = A.length;
  const M = A.map((r, i) => [...r, b[i]]);
  for (let i = 0; i < n; i++) {
    let mx = i;
    for (let j = i + 1; j < n; j++) if (Math.abs(M[j][i]) > Math.abs(M[mx][i])) mx = j;
    [M[i], M[mx]] = [M[mx], M[i]];
    if (Math.abs(M[i][i]) < 1e-15) continue;
    for (let j = i + 1; j < n; j++) {
      const f = M[j][i] / M[i][i];
      for (let k = i; k <= n; k++) M[j][k] -= f * M[i][k];
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}
function ols(Y, X) {
  const n = Y.length, p = X[0].length + 1;
  const Xa = X.map(r => [1, ...r]);
  const XtX = Array.from({ length: p }, () => new Array(p).fill(0));
  const XtY = new Array(p).fill(0);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < p; j++) {
      XtY[j] += Xa[i][j] * Y[i];
      for (let l = 0; l < p; l++) XtX[j][l] += Xa[i][j] * Xa[i][l];
    }
  const beta = solveLinear(XtX, XtY);
  const resid = Y.map((y, i) => y - Xa[i].reduce((s, x, j) => s + x * beta[j], 0));
  const rss = resid.reduce((s, r) => s + r * r, 0);
  const tss = Y.reduce((s, y) => s + (y - mean(Y)) ** 2, 0);
  return { beta, resid, rss, tss, r2: tss > 0 ? 1 - rss / tss : 0 };
}
function looCV(Y, X) {
  const n = Y.length; let ss = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xt = [...X.slice(0, i), ...X.slice(i + 1)];
    try {
      const f = ols(Yt, Xt);
      const xi = [1, ...X[i]];
      ss += (Y[i] - xi.reduce((s, x, j) => s + x * f.beta[j], 0)) ** 2;
    } catch (e) { ss += (Y[i] - mean(Y)) ** 2; }
  }
  return Math.sqrt(ss / n);
}
function gapPct(rms, sdy) { return sdy > 0 ? 100 * (1 - rms ** 2 / sdy ** 2) : 0; }

const Y_all = gals45.map(g => g.logA0);
const logMHI_all = gals45.map(g => g.logMHI);
const logMhost_all = gals45.map(g => tdMap[g.name].logMhost);
const logMR_all = gals45.map(g => g.logMeanRun);
const logVflat_all = gals45.map(g => Math.log10(sparcMap[g.name].Vflat));
const Vflat_all = gals45.map(g => sparcMap[g.name].Vflat);
const logMbar_all = gals45.map(g => {
  const s = sparcMap[g.name];
  return Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
});
const logL36_all = gals45.map(g => Math.log10(Math.max(sparcMap[g.name].L36, 0.01)));
const logRdisk_all = gals45.map(g => Math.log10(sparcMap[g.name].Rdisk));
const morphT_all = gals45.map(g => sparcMap[g.name].T);

function computeVfResid(indices) {
  const logVflat = indices.map(i => logVflat_all[i]);
  const structX = indices.map(i => [logMbar_all[i], logL36_all[i], logRdisk_all[i], morphT_all[i]]);
  return ols(logVflat, structX).resid;
}

function evalSubset(indices, label) {
  const n = indices.length;
  if (n < 8) return null;
  const Y = indices.map(i => Y_all[i]);
  const sdY = sd(Y);
  if (sdY < 0.005) return null;
  const VfResid = computeVfResid(indices);
  const core3X = indices.map((_, j) => [logMHI_all[indices[j]], logMhost_all[indices[j]], logMR_all[indices[j]]]);
  const coreVfResidX = indices.map((_, j) => [...core3X[j], VfResid[j]]);
  const coreVflatX = indices.map((_, j) => [...core3X[j], logVflat_all[indices[j]]]);

  const coreGap = gapPct(looCV(Y, core3X), sdY);
  const vfResidGap = gapPct(looCV(Y, coreVfResidX), sdY);
  const vflatGap = gapPct(looCV(Y, coreVflatX), sdY);
  const rVfResid = pearsonR(VfResid, Y);
  const coreFit = ols(Y, core3X);
  const fullFit = ols(Y, coreVfResidX);

  return {
    label, n, sdY: +sdY.toFixed(4),
    coreGap: +coreGap.toFixed(1), vflatGap: +vflatGap.toFixed(1),
    vfResidGap: +vfResidGap.toFixed(1), delta: +(vfResidGap - coreGap).toFixed(1),
    rVfResid: +rVfResid.toFixed(3),
    coreBeta: coreFit.beta.map(b => +b.toFixed(4)),
    fullBeta: fullFit.beta.map(b => +b.toFixed(4)),
    coreR2: +coreFit.r2.toFixed(3), fullR2: +fullFit.r2.toFixed(3)
  };
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: VFLAT THRESHOLD SCAN — Where does VfResid activate?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const thresholds = [70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170];
const scanResults = [];

console.log('  ' + 'Vflat>'.padEnd(10) + '  N   core%  VfRes%  delta   r(VfR)  VfR coef');
console.log('  ' + '─'.repeat(70));

for (const thresh of thresholds) {
  const idx = gals45.map((_, i) => i).filter(i => Vflat_all[i] >= thresh);
  const r = evalSubset(idx, 'Vflat>=' + thresh);
  if (r) {
    scanResults.push({ threshold: thresh, ...r });
    console.log('  ' + ('>=' + thresh).padEnd(10) + r.n.toString().padStart(3) + '   ' +
      r.coreGap.toFixed(1).padStart(5) + '  ' + r.vfResidGap.toFixed(1).padStart(6) + '  ' +
      (r.delta > 0 ? '+' : '') + r.delta.toFixed(1).padStart(5) + '   ' +
      r.rVfResid.toFixed(3).padStart(6) + '  ' + r.fullBeta[4].toFixed(3));
  }
}
console.log('');

const bestScan = scanResults.reduce((a, b) => b.delta > a.delta && b.n >= 15 ? b : a, scanResults[0]);
console.log('  OPTIMAL REGIME: Vflat >= ' + bestScan.threshold + ' (N=' + bestScan.n + ')');
console.log('    Core gap: ' + bestScan.coreGap.toFixed(1) + '%');
console.log('    Core+VfResid gap: ' + bestScan.vfResidGap.toFixed(1) + '%');
console.log('    Delta: +' + bestScan.delta.toFixed(1) + 'pp');
console.log('    r(VfResid, a₀): ' + bestScan.rVfResid.toFixed(3));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: SMOOTH vs SHARP TRANSITION');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const VfResid_full = computeVfResid(gals45.map((_, i) => i));

const sorted = gals45.map((_, i) => ({ i, vf: Vflat_all[i], vfr: VfResid_full[i], y: Y_all[i] }))
  .sort((a, b) => a.vf - b.vf);

const windowSize = 15;
const slideResults = [];

for (let start = 0; start <= sorted.length - windowSize; start++) {
  const window = sorted.slice(start, start + windowSize);
  const idx = window.map(w => w.i);
  const vfMid = window[Math.floor(windowSize / 2)].vf;
  const vfRange = [window[0].vf, window[window.length - 1].vf];

  const Y = idx.map(i => Y_all[i]);
  const sdY_w = sd(Y);
  if (sdY_w < 0.005) continue;

  const VfR = computeVfResid(idx);
  const rVfR = pearsonR(VfR, Y);
  slideResults.push({ vfMid, vfRange, n: windowSize, rVfResid: rVfR });
}

console.log('  Sliding window r(VfResid, a₀) across Vflat (window=' + windowSize + '):');
console.log('  ' + 'Vflat(mid)'.padEnd(12) + '  range'.padEnd(16) + '  r(VfResid)');
console.log('  ' + '─'.repeat(45));
for (const sr of slideResults) {
  const bar = '█'.repeat(Math.max(0, Math.round((sr.rVfResid + 0.5) * 20)));
  console.log('  ' + sr.vfMid.toFixed(0).padStart(8) + '    ' +
    (sr.vfRange[0].toFixed(0) + '-' + sr.vfRange[1].toFixed(0)).padEnd(12) + sr.rVfResid.toFixed(3).padStart(7) + '  ' + bar);
}
console.log('');

const transitionVflat = slideResults.find(s => s.rVfResid > 0.5);
if (transitionVflat) {
  console.log('  TRANSITION POINT: r > 0.5 first appears at Vflat ~ ' + transitionVflat.vfMid.toFixed(0) + ' km/s');
} else {
  console.log('  NO CLEAR TRANSITION: r never consistently exceeds 0.5');
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: REGIME-SPECIFIC MODEL — Coefficients comparison');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const fullIdx = gals45.map((_, i) => i);
const highIdx = gals45.map((_, i) => i).filter(i => Vflat_all[i] >= bestScan.threshold);
const lowIdx = gals45.map((_, i) => i).filter(i => Vflat_all[i] < bestScan.threshold);

const fullEval = evalSubset(fullIdx, 'Full N=45');
const highEval = evalSubset(highIdx, 'High-Vflat (>=' + bestScan.threshold + ')');
const lowEval = lowIdx.length >= 8 ? evalSubset(lowIdx, 'Low-Vflat (<' + bestScan.threshold + ')') : null;

const axes = ['intercept', 'logMHI', 'logMhost', 'logMR', 'VfResid'];

console.log('  Core+VfResid coefficients:');
console.log('  ' + 'Axis'.padEnd(12) + '  Full'.padStart(10) + '  High-Vf'.padStart(10) + (lowEval ? '  Low-Vf'.padStart(10) : ''));
console.log('  ' + '─'.repeat(45));
for (let j = 0; j < axes.length; j++) {
  const line = '  ' + axes[j].padEnd(12) +
    fullEval.fullBeta[j].toFixed(3).padStart(10) +
    highEval.fullBeta[j].toFixed(3).padStart(10) +
    (lowEval ? lowEval.fullBeta[j].toFixed(3).padStart(10) : '');
  console.log(line);
}
console.log('');
console.log('  R²:  Full=' + fullEval.fullR2.toFixed(3) + '  High=' + highEval.fullR2.toFixed(3) + (lowEval ? '  Low=' + lowEval.fullR2.toFixed(3) : ''));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: BOOTSTRAP STABILITY — Regime-specific');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function bootstrapRegime(indices, nBoot) {
  const n = indices.length;
  const deltas = [];
  const vfResidCoefs = [];
  for (let b = 0; b < nBoot; b++) {
    const bIdx = Array.from({ length: n }, () => indices[Math.floor(Math.random() * n)]);
    const Y = bIdx.map(i => Y_all[i]);
    const sdY_b = sd(Y);
    if (sdY_b < 0.005) continue;
    try {
      const VfR = computeVfResid(bIdx);
      const core3X = bIdx.map(i => [logMHI_all[i], logMhost_all[i], logMR_all[i]]);
      const coreVfRX = bIdx.map((_, j) => [...core3X[j], VfR[j]]);
      const coreGap = gapPct(looCV(Y, core3X), sdY_b);
      const vfResidGap = gapPct(looCV(Y, coreVfRX), sdY_b);
      deltas.push(vfResidGap - coreGap);
      vfResidCoefs.push(ols(Y, coreVfRX).beta[4]);
    } catch (e) { }
  }
  deltas.sort((a, b) => a - b);
  vfResidCoefs.sort((a, b) => a - b);
  return {
    n: deltas.length,
    deltaP50: deltas[Math.floor(deltas.length * 0.5)],
    deltaP5: deltas[Math.floor(deltas.length * 0.05)],
    deltaP95: deltas[Math.floor(deltas.length * 0.95)],
    pDeltaGt0: deltas.filter(d => d > 0).length / deltas.length,
    coefP50: vfResidCoefs[Math.floor(vfResidCoefs.length * 0.5)],
    coefP5: vfResidCoefs[Math.floor(vfResidCoefs.length * 0.05)],
    coefP95: vfResidCoefs[Math.floor(vfResidCoefs.length * 0.95)],
    coefPosRate: vfResidCoefs.filter(c => c > 0).length / vfResidCoefs.length
  };
}

const nBoot = 1000;

console.log('  Full sample (N=45):');
const bFull = bootstrapRegime(fullIdx, nBoot);
console.log('    Delta: 5th=' + bFull.deltaP5.toFixed(1) + '  50th=' + bFull.deltaP50.toFixed(1) + '  95th=' + bFull.deltaP95.toFixed(1));
console.log('    P(delta>0): ' + (bFull.pDeltaGt0 * 100).toFixed(1) + '%');
console.log('    VfResid coef: 5th=' + bFull.coefP5.toFixed(3) + '  50th=' + bFull.coefP50.toFixed(3) + '  95th=' + bFull.coefP95.toFixed(3));
console.log('    P(coef>0): ' + (bFull.coefPosRate * 100).toFixed(1) + '%\n');

console.log('  High-Vflat regime (N=' + highIdx.length + '):');
const bHigh = bootstrapRegime(highIdx, nBoot);
console.log('    Delta: 5th=' + bHigh.deltaP5.toFixed(1) + '  50th=' + bHigh.deltaP50.toFixed(1) + '  95th=' + bHigh.deltaP95.toFixed(1));
console.log('    P(delta>0): ' + (bHigh.pDeltaGt0 * 100).toFixed(1) + '%');
console.log('    VfResid coef: 5th=' + bHigh.coefP5.toFixed(3) + '  50th=' + bHigh.coefP50.toFixed(3) + '  95th=' + bHigh.coefP95.toFixed(3));
console.log('    P(coef>0): ' + (bHigh.coefPosRate * 100).toFixed(1) + '%\n');

if (lowIdx.length >= 10) {
  console.log('  Low-Vflat regime (N=' + lowIdx.length + '):');
  const bLow = bootstrapRegime(lowIdx, nBoot);
  console.log('    Delta: 5th=' + bLow.deltaP5.toFixed(1) + '  50th=' + bLow.deltaP50.toFixed(1) + '  95th=' + bLow.deltaP95.toFixed(1));
  console.log('    P(delta>0): ' + (bLow.pDeltaGt0 * 100).toFixed(1) + '%');
  console.log('    VfResid coef: 5th=' + bLow.coefP5.toFixed(3) + '  50th=' + bLow.coefP50.toFixed(3) + '  95th=' + bLow.coefP95.toFixed(3));
  console.log('    P(coef>0): ' + (bLow.coefPosRate * 100).toFixed(1) + '%\n');
}

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 133A: VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const regimeStrength = bestScan.delta > 25 ? 'VERY_STRONG_REGIME' :
  bestScan.delta > 15 ? 'STRONG_REGIME' :
    bestScan.delta > 8 ? 'MODERATE_REGIME' : 'WEAK_REGIME';

const transitionType = slideResults.every(s => s.rVfResid > 0.3) ? 'SMOOTH_ALWAYS_ON' :
  slideResults.some(s => s.rVfResid > 0.6) && slideResults.some(s => s.rVfResid < 0.2) ? 'SHARP_TRANSITION' :
    'GRADUAL_TRANSITION';

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  REGIME: ' + regimeStrength.padEnd(56) + '║');
console.log('  ║  TRANSITION: ' + transitionType.padEnd(52) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════════╝\n');

console.log('  OPTIMAL REGIME: Vflat >= ' + bestScan.threshold + ' km/s');
console.log('    N = ' + bestScan.n + ' galaxies');
console.log('    Core+VfResid gap = ' + bestScan.vfResidGap.toFixed(1) + '% (+' + bestScan.delta.toFixed(1) + 'pp)');
console.log('    r(VfResid, a₀) = ' + bestScan.rVfResid.toFixed(3));
console.log('');

const output = {
  phase: '133A',
  title: 'Regime Law — High-Vflat Dedicated Model',
  regimeStrength, transitionType,
  optimalThreshold: bestScan.threshold,
  optimalN: bestScan.n,
  scanResults: scanResults.map(s => ({
    threshold: s.threshold, n: s.n,
    coreGap: s.coreGap, vfResidGap: s.vfResidGap,
    delta: s.delta, rVfResid: s.rVfResid
  })),
  slidingWindow: slideResults.map(s => ({ vfMid: +s.vfMid.toFixed(0), rVfResid: +s.rVfResid.toFixed(3) })),
  coefficients: { full: fullEval, high: highEval, low: lowEval },
  bootstrap: { full: bFull, high: bHigh }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase133a-regime-law.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase133a-regime-law.json');

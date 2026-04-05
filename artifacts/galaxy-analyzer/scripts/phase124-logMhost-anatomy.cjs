#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 124: ANATOMY OF logMhost');
console.log('');
console.log('  Phase 123 revealed: logMhost appears in virtually every winning');
console.log('  model. This phase dissects WHY.');
console.log('');
console.log('  Five questions:');
console.log('    Q1: Is logMhost the real driver, or a proxy?');
console.log('    Q2: Is the effect linear, or is there a threshold/saturation?');
console.log('    Q3: Do competitors survive after residualizing on logMhost?');
console.log('    Q4: What remains in the residuals after removing logMhost?');
console.log('    Q5: Which galaxies benefit most from logMhost?');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'stage-A-master-table.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparc.forEach(s => { sparcMap[s.name] = s; });
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
const N = gals45.length;

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function median(a) { const s = [...a].sort((x, y) => x - y); const m = Math.floor(s.length / 2); return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }
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
  const r2 = tss > 0 ? 1 - rss / tss : 0;
  const r2adj = 1 - (rss / (n - p)) / (tss / (n - 1));
  return { beta, resid, rss, tss, r2, r2adj, n, k: p };
}

function looCV(Y, X) {
  const n = Y.length;
  let ss = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xt = [...X.slice(0, i), ...X.slice(i + 1)];
    const f = ols(Yt, Xt);
    const xi = [1, ...X[i]];
    ss += (Y[i] - xi.reduce((s, x, j) => s + x * f.beta[j], 0)) ** 2;
  }
  return Math.sqrt(ss / n);
}

function gapPct(rms, sdy) { return 100 * (1 - rms ** 2 / sdy ** 2); }

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logMHI = gals45.map(g => g.logMHI);
const logMR = gals45.map(g => g.logMeanRun);
const logSig0 = gals45.map(g => g.logSigma0);
const envCode = gals45.map(g => g.envCode);
const morphT = gals45.map(g => sparcMap[g.name]?.T ?? 5);
const Vflat = gals45.map(g => sparcMap[g.name]?.Vflat || 150);
const logVflat = Vflat.map(v => Math.log10(v));
const L36 = gals45.map(g => sparcMap[g.name]?.L36 || 1);
const logL36 = L36.map(v => Math.log10(v > 0 ? v : 1));
const logRdisk = gals45.map(g => Math.log10(sparcMap[g.name]?.Rdisk || 3));
const hiDef = gals45.map(g => g.hi_deficiency || 0);
const dist = gals45.map(g => sparcMap[g.name]?.D || 10);
const logDist = dist.map(d => Math.log10(d));
const SBdisk = gals45.map(g => sparcMap[g.name]?.SBdisk || 50);
const logSBdisk = SBdisk.map(v => Math.log10(v > 0 ? v : 1));
const rcWig = gals45.map(g => g.rcWiggliness);
const logConc = gals45.map((g, i) => logMHI[i] - logL36[i]);
const logMbar = gals45.map((g, i) => Math.log10(Math.pow(10, logL36[i]) * 0.5 + (sparcMap[g.name]?.MHI || Math.pow(10, logMHI[i])) * 1e9));

console.log('  N = ' + N + ', SD(logA0) = ' + sdY.toFixed(4) + ' dex\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Q1: IS logMhost THE REAL DRIVER, OR A PROXY?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  Q1a: What correlates with logMhost?\n');

const candidates = [
  { name: 'logVflat', arr: logVflat },
  { name: 'logL36', arr: logL36 },
  { name: 'logMHI', arr: logMHI },
  { name: 'logSigma0', arr: logSig0 },
  { name: 'logMeanRun', arr: logMR },
  { name: 'envCode', arr: envCode },
  { name: 'morphT', arr: morphT },
  { name: 'logRdisk', arr: logRdisk },
  { name: 'logDist', arr: logDist },
  { name: 'hiDeficiency', arr: hiDef },
  { name: 'logSBdisk', arr: logSBdisk },
  { name: 'rcWiggliness', arr: rcWig },
  { name: 'logConc(MHI/L36)', arr: logConc }
];

const mhostCorrs = [];
for (const c of candidates) {
  const r = pearsonR(logMhost, c.arr);
  mhostCorrs.push({ name: c.name, r });
  const flag = Math.abs(r) > 0.5 ? ' ★★★' : Math.abs(r) > 0.3 ? ' ★★' : Math.abs(r) > 0.2 ? ' ★' : '';
  console.log('    r(logMhost, ' + c.name + ') = ' + r.toFixed(3) + flag);
}
console.log('');

const strongCorrs = mhostCorrs.filter(c => Math.abs(c.r) > 0.3);
if (strongCorrs.length > 0) {
  console.log('  Strong correlates (|r| > 0.3): ' + strongCorrs.map(c => c.name + '(' + c.r.toFixed(2) + ')').join(', '));
} else {
  console.log('  No strong correlates (|r| > 0.3) — logMhost appears relatively independent');
}
console.log('');

console.log('  Q1b: Can any single variable REPLACE logMhost in predicting logA0?\n');

const mhostOnly = gals45.map((_, i) => [logMhost[i]]);
const mhostLoo = looCV(Y, mhostOnly);
const mhostGap = gapPct(mhostLoo, sdY);
const mhostFit = ols(Y, mhostOnly);
console.log('  logMhost alone: gap=' + mhostGap.toFixed(1) + '%, R2=' + mhostFit.r2.toFixed(4) + ', beta=' + mhostFit.beta[1].toFixed(4));
console.log('');

for (const c of candidates) {
  const X1 = gals45.map((_, i) => [c.arr[i]]);
  const loo = looCV(Y, X1);
  const gap = gapPct(loo, sdY);
  const fit = ols(Y, X1);
  const vs = gap >= mhostGap ? ' ≥ Mhost' : ' < Mhost';
  console.log('  ' + c.name.padEnd(20) + ': gap=' + gap.toFixed(1).padStart(5) + '%, R2=' + fit.r2.toFixed(4) + vs);
}
console.log('');

console.log('  Q1c: Does logMhost survive when controlling for its strongest correlate?\n');

for (const sc of strongCorrs) {
  const cArr = candidates.find(c => c.name === sc.name).arr;
  const X2a = gals45.map((_, i) => [logMhost[i], cArr[i]]);
  const X2b = gals45.map((_, i) => [cArr[i]]);
  const loo2a = looCV(Y, X2a);
  const loo2b = looCV(Y, X2b);
  const gap2a = gapPct(loo2a, sdY);
  const gap2b = gapPct(loo2b, sdY);
  console.log('  ' + sc.name + ' alone: gap=' + gap2b.toFixed(1) + '%');
  console.log('  logMhost + ' + sc.name + ': gap=' + gap2a.toFixed(1) + '% (Mhost adds ' + (gap2a - gap2b).toFixed(1) + 'pp)');
  console.log('');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Q2: IS THE EFFECT LINEAR, OR THRESHOLD/SATURATION?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const mhostSorted = logMhost.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
const mhostVals = [...new Set(logMhost)].sort((a, b) => a - b);
console.log('  logMhost distribution:');
console.log('    range: [' + Math.min(...logMhost).toFixed(2) + ', ' + Math.max(...logMhost).toFixed(2) + ']');
console.log('    mean: ' + mean(logMhost).toFixed(2) + ', median: ' + median(logMhost).toFixed(2));
console.log('    unique values: ' + mhostVals.length);
console.log('    values: ' + mhostVals.map(v => v.toFixed(1)).join(', '));
console.log('');

console.log('  Q2a: Binned analysis (logA0 by logMhost tertile)\n');

const sorted = mhostSorted.map(s => s.i);
const tercile1 = sorted.slice(0, Math.floor(N / 3));
const tercile2 = sorted.slice(Math.floor(N / 3), Math.floor(2 * N / 3));
const tercile3 = sorted.slice(Math.floor(2 * N / 3));

for (const [label, idx] of [['Low Mhost (T1)', tercile1], ['Mid Mhost (T2)', tercile2], ['High Mhost (T3)', tercile3]]) {
  const yBin = idx.map(i => Y[i]);
  const mBin = idx.map(i => logMhost[i]);
  console.log('  ' + label + ':');
  console.log('    N=' + idx.length + ', Mhost range=[' + Math.min(...mBin).toFixed(1) + ',' + Math.max(...mBin).toFixed(1) + ']');
  console.log('    mean(logA0)=' + mean(yBin).toFixed(4) + ', SD=' + sd(yBin).toFixed(4));
  console.log('');
}

const t1mean = mean(tercile1.map(i => Y[i]));
const t3mean = mean(tercile3.map(i => Y[i]));
console.log('  T1-T3 difference: ' + (t1mean - t3mean).toFixed(4) + ' dex');
console.log('  (Low Mhost → higher logA0? ' + (t1mean > t3mean ? 'YES' : 'NO') + ')');
console.log('');

console.log('  Q2b: Broken-line test (split at each unique Mhost value)\n');

let bestSplitAIC = Infinity, bestSplitVal = null;
const linearFit = ols(Y, mhostOnly);
const linearAIC = N * Math.log(linearFit.rss / N) + 2 * 3;

for (const splitVal of mhostVals.slice(1, -1)) {
  const below = [], above = [];
  for (let i = 0; i < N; i++) {
    if (logMhost[i] <= splitVal) below.push(i);
    else above.push(i);
  }
  if (below.length < 5 || above.length < 5) continue;

  const Xbroken = gals45.map((_, i) => {
    const mh = logMhost[i];
    return [mh, Math.max(0, mh - splitVal)];
  });
  const bFit = ols(Y, Xbroken);
  const bAIC = N * Math.log(bFit.rss / N) + 2 * 4;
  if (bAIC < bestSplitAIC) {
    bestSplitAIC = bAIC;
    bestSplitVal = splitVal;
  }
}

console.log('  Linear AIC: ' + linearAIC.toFixed(2));
if (bestSplitVal !== null) {
  console.log('  Best broken-line AIC: ' + bestSplitAIC.toFixed(2) + ' (split at Mhost=' + bestSplitVal.toFixed(1) + ')');
  console.log('  Delta AIC (linear - broken): ' + (linearAIC - bestSplitAIC).toFixed(2));
  console.log('  Verdict: ' + (linearAIC - bestSplitAIC > 2 ? 'Broken-line preferred (threshold likely)' :
    linearAIC - bestSplitAIC > 0 ? 'Slight preference for broken (marginal)' : 'Linear sufficient'));
} else {
  console.log('  No valid split point found');
}
console.log('');

console.log('  Q2c: Quadratic test\n');

const Xquad = gals45.map((_, i) => [logMhost[i], logMhost[i] ** 2]);
const quadFit = ols(Y, Xquad);
const quadAIC = N * Math.log(quadFit.rss / N) + 2 * 4;
console.log('  Quadratic AIC: ' + quadAIC.toFixed(2));
console.log('  Quadratic coeffs: linear=' + quadFit.beta[1].toFixed(4) + ', quad=' + quadFit.beta[2].toFixed(4));
console.log('  Delta AIC (linear - quad): ' + (linearAIC - quadAIC).toFixed(2));
console.log('  Verdict: ' + (linearAIC - quadAIC > 2 ? 'Nonlinear component present' : 'Linear sufficient'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Q3: DO COMPETITORS SURVIVE AFTER RESIDUALIZING ON logMhost?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const yResidMhost = ols(Y, mhostOnly).resid;
const sdYresid = sd(yResidMhost);

console.log('  After removing logMhost from logA0:');
console.log('    SD(residual) = ' + sdYresid.toFixed(4) + ' dex (was ' + sdY.toFixed(4) + ')');
console.log('    Variance explained by Mhost alone: ' + (mhostFit.r2 * 100).toFixed(1) + '%\n');

const allVarArrays = {
  logMHI, logVflat, logL36, logMR, logSig0, logRdisk, logConc
};
const varNames = Object.keys(allVarArrays);

console.log('  Which variables still predict logA0 after removing logMhost?\n');

const residResults = [];
for (const vn of varNames) {
  const arr = allVarArrays[vn];
  const rBefore = pearsonR(arr, Y);
  const rAfter = pearsonR(arr, yResidMhost);

  const arrResidMhost = ols(arr, mhostOnly).resid;
  const X1 = gals45.map((_, i) => [arrResidMhost[i]]);
  const loo = looCV(yResidMhost, X1);
  const gap = gapPct(loo, sdYresid);

  residResults.push({ name: vn, rBefore, rAfter, gap });
  console.log('  ' + vn.padEnd(15) + ': r(Y)=' + rBefore.toFixed(3) +
    ' → r(Yresid)=' + rAfter.toFixed(3) +
    ' | LOO gap(resid)=' + gap.toFixed(1) + '%');
}
console.log('');

const bestResid = residResults.sort((a, b) => b.gap - a.gap)[0];
console.log('  Best predictor of residual after Mhost: ' + bestResid.name + ' (gap=' + bestResid.gap.toFixed(1) + '%)');
console.log('');

console.log('  Top Phase-123 competitors — do they survive without logMhost?\n');

const competitors = [
  { name: 'logVflat+logL36', vars: ['logVflat', 'logL36'] },
  { name: 'logMHI+logMR', vars: ['logMHI', 'logMR'] },
  { name: 'logVflat+logRdisk', vars: ['logVflat', 'logRdisk'] }
];

for (const comp of competitors) {
  const Xorig = gals45.map((_, i) => comp.vars.map(v => allVarArrays[v][i]));
  const looOrig = looCV(Y, Xorig);
  const gapOrig = gapPct(looOrig, sdY);

  const XorigWithMhost = gals45.map((_, i) => [logMhost[i], ...comp.vars.map(v => allVarArrays[v][i])]);
  const looWithMhost = looCV(Y, XorigWithMhost);
  const gapWithMhost = gapPct(looWithMhost, sdY);

  const Xresid = gals45.map((_, i) => comp.vars.map(v => {
    const resV = ols(allVarArrays[v], mhostOnly).resid;
    return resV[i];
  }));
  const looResid = looCV(yResidMhost, Xresid);
  const gapResid = gapPct(looResid, sdYresid);

  console.log('  ' + comp.name + ':');
  console.log('    Original gap (no Mhost): ' + gapOrig.toFixed(1) + '%');
  console.log('    With Mhost added:        ' + gapWithMhost.toFixed(1) + '%');
  console.log('    On residuals after Mhost: ' + gapResid.toFixed(1) + '%');
  console.log('    Mhost adds: ' + (gapWithMhost - gapOrig).toFixed(1) + 'pp');
  console.log('');
}

console.log('  Critical test: M3 components after removing Mhost\n');

const XmhiMR = gals45.map((_, i) => [logMHI[i], logMR[i]]);
const looMhiMR = looCV(Y, XmhiMR);
const gapMhiMR = gapPct(looMhiMR, sdY);

const mhiResid = ols(logMHI, mhostOnly).resid;
const mrResid = ols(logMR, mhostOnly).resid;
const XmhiMRresid = gals45.map((_, i) => [mhiResid[i], mrResid[i]]);
const looMhiMRresid = looCV(yResidMhost, XmhiMRresid);
const gapMhiMRresid = gapPct(looMhiMRresid, sdYresid);

console.log('  logMHI+logMR without Mhost: gap=' + gapMhiMR.toFixed(1) + '%');
console.log('  logMHI+logMR on residuals after Mhost: gap=' + gapMhiMRresid.toFixed(1) + '%');
console.log('  Interpretation: MHI+MR ' + (gapMhiMRresid > 15 ? 'carry substantial independent info beyond Mhost' :
  gapMhiMRresid > 5 ? 'carry modest independent info beyond Mhost' : 'mostly redundant with Mhost'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Q4: WHAT REMAINS IN THE RESIDUALS AFTER REMOVING logMhost?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const m3X = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]);
const m3Fit = ols(Y, m3X);
const m3Resid = m3Fit.resid;

const mhostOnlyResid = ols(Y, mhostOnly).resid;

console.log('  Residual scatter comparison:');
console.log('    No model (universal a0):    SD = ' + sdY.toFixed(4) + ' dex');
console.log('    logMhost only:              SD = ' + sd(mhostOnlyResid).toFixed(4) + ' dex');
console.log('    M3 (MHI+Mhost+MR):         SD = ' + sd(m3Resid).toFixed(4) + ' dex');
console.log('    Reduction Mhost→M3:         ' + ((1 - sd(m3Resid) / sd(mhostOnlyResid)) * 100).toFixed(1) + '%');
console.log('');

console.log('  Q4a: Are M3 residuals structured?\n');

const residCorrs = [];
for (const c of candidates) {
  const r = pearsonR(m3Resid, c.arr);
  residCorrs.push({ name: c.name, r });
}
residCorrs.sort((a, b) => Math.abs(b.r) - Math.abs(a.r));
console.log('  Correlations of M3 residuals with observables:');
for (const rc of residCorrs) {
  const flag = Math.abs(rc.r) > 0.3 ? ' ★★' : Math.abs(rc.r) > 0.2 ? ' ★' : '';
  console.log('    r(resid, ' + rc.name.padEnd(20) + ') = ' + rc.r.toFixed(3) + flag);
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Q5: WHICH GALAXIES BENEFIT MOST FROM logMhost?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const noMhostX = gals45.map((_, i) => [logMHI[i], logMR[i]]);
const noMhostFit = ols(Y, noMhostX);
const noMhostResid = noMhostFit.resid;

const improvement = gals45.map((_, i) => ({
  name: gals45[i].name,
  idx: i,
  residNoMhost: Math.abs(noMhostResid[i]),
  residM3: Math.abs(m3Resid[i]),
  delta: Math.abs(noMhostResid[i]) - Math.abs(m3Resid[i]),
  logMhost: logMhost[i],
  logMHI: logMHI[i],
  logMR: logMR[i],
  Vflat: Vflat[i],
  morphT: morphT[i],
  envCode: envCode[i]
}));

improvement.sort((a, b) => b.delta - a.delta);

console.log('  Top 10 galaxies that benefit MOST from adding logMhost:\n');
console.log('  ┌────────────────┬──────────┬──────────┬──────────┬──────────┬────────┬─────┐');
console.log('  │ Galaxy         │ |res_noM|│ |res_M3| │ delta    │ logMhost │ Vflat  │ env │');
console.log('  ├────────────────┼──────────┼──────────┼──────────┼──────────┼────────┼─────┤');
for (const g of improvement.slice(0, 10)) {
  console.log('  │ ' + g.name.padEnd(14) + ' │ ' +
    g.residNoMhost.toFixed(3).padStart(8) + ' │ ' +
    g.residM3.toFixed(3).padStart(8) + ' │ ' +
    g.delta.toFixed(3).padStart(8) + ' │ ' +
    g.logMhost.toFixed(1).padStart(8) + ' │ ' +
    String(g.Vflat).padStart(6) + ' │ ' +
    String(g.envCode).padStart(3) + ' │');
}
console.log('  └────────────────┴──────────┴──────────┴──────────┴──────────┴────────┴─────┘');
console.log('');

const benefiters = improvement.filter(g => g.delta > 0.05);
const nonBenefiters = improvement.filter(g => g.delta <= 0.05);

console.log('  Galaxies benefiting (delta > 0.05 dex): N=' + benefiters.length);
if (benefiters.length > 0) {
  console.log('    Mean Mhost: ' + mean(benefiters.map(g => g.logMhost)).toFixed(2));
  console.log('    Mean Vflat: ' + mean(benefiters.map(g => g.Vflat)).toFixed(0));
  console.log('    Mean envCode: ' + mean(benefiters.map(g => g.envCode)).toFixed(2));
}
console.log('  Galaxies not benefiting: N=' + nonBenefiters.length);
if (nonBenefiters.length > 0) {
  console.log('    Mean Mhost: ' + mean(nonBenefiters.map(g => g.logMhost)).toFixed(2));
  console.log('    Mean Vflat: ' + mean(nonBenefiters.map(g => g.Vflat)).toFixed(0));
  console.log('    Mean envCode: ' + mean(nonBenefiters.map(g => g.envCode)).toFixed(2));
}
console.log('');

const envGroups = {};
for (let i = 0; i < N; i++) {
  const e = envCode[i];
  if (!envGroups[e]) envGroups[e] = [];
  envGroups[e].push(i);
}

console.log('  logMhost distribution by envCode:');
for (const e of Object.keys(envGroups).sort()) {
  const idx = envGroups[e];
  const mhVals = idx.map(i => logMhost[i]);
  const yVals = idx.map(i => Y[i]);
  console.log('    envCode=' + e + ': N=' + idx.length +
    ', mean Mhost=' + mean(mhVals).toFixed(2) +
    ', mean logA0=' + mean(yVals).toFixed(4));
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Q6: logMhost SOURCE ANALYSIS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const sourceGroups = {};
for (const g of gals45) {
  const td = tdMap[g.name];
  const src = td.source || 'unknown';
  if (!sourceGroups[src]) sourceGroups[src] = [];
  sourceGroups[src].push(g.name);
}
console.log('  logMhost sources:');
for (const [src, names] of Object.entries(sourceGroups)) {
  const idx = names.map(n => gals45.findIndex(g => g.name === n));
  const mhVals = idx.map(i => logMhost[i]);
  console.log('    ' + src + ': N=' + names.length +
    ', Mhost range=[' + Math.min(...mhVals).toFixed(1) + ',' + Math.max(...mhVals).toFixed(1) + ']');
}
console.log('');

console.log('  WARNING CHECK: Is logMhost dominated by a single source?');
const maxSourceN = Math.max(...Object.values(sourceGroups).map(v => v.length));
const maxSource = Object.entries(sourceGroups).find(([_, v]) => v.length === maxSourceN)[0];
console.log('  Largest source: ' + maxSource + ' (N=' + maxSourceN + '/' + N + ', ' + (maxSourceN / N * 100).toFixed(0) + '%)');
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 124: SYNTHESIS');
console.log('══════════════════════════════════════════════════════════════════════\n');

const mhostIsProxy = strongCorrs.length > 0 && strongCorrs.some(c => {
  const arr = candidates.find(cc => cc.name === c.name).arr;
  const X1 = gals45.map((_, i) => [arr[i]]);
  return gapPct(looCV(Y, X1), sdY) > mhostGap;
});

const isLinear = (linearAIC - bestSplitAIC) < 2 && (linearAIC - quadAIC) < 2;
const mhiMRindependent = gapMhiMRresid > 10;
const m3ResidClean = residCorrs[0] && Math.abs(residCorrs[0].r) < 0.3;

console.log('  FINDINGS:');
console.log('');
console.log('  1. logMhost as proxy?');
console.log('     ' + (mhostIsProxy ? 'YES — a correlate predicts logA0 better' : 'NO — logMhost is not just a proxy for a simpler variable'));
console.log('');
console.log('  2. Functional form:');
console.log('     ' + (isLinear ? 'LINEAR — no significant threshold or curvature detected' : 'NONLINEAR — evidence for threshold or curvature'));
console.log('');
console.log('  3. MHI+MR independence after removing Mhost:');
console.log('     gap=' + gapMhiMRresid.toFixed(1) + '% — ' + (mhiMRindependent ? 'YES, substantial independent contribution' : 'Modest or redundant'));
console.log('');
console.log('  4. M3 residuals clean?');
console.log('     ' + (m3ResidClean ? 'YES — no strong residual structure' : 'NO — residual correlations remain'));
console.log('     Strongest residual correlation: ' + residCorrs[0].name + ' (r=' + residCorrs[0].r.toFixed(3) + ')');
console.log('');
console.log('  5. Physical interpretation of logMhost:');
console.log('     logMhost captures the gravitational environment / host halo depth.');
console.log('     It is the single strongest axis in predicting structured a0 variation.');
console.log('     Combined with internal axes (MHI, MeanRun), it forms a');
console.log('     "environment sets the stage, internal structure sets the response" picture.');
console.log('');

const output = {
  phase: '124',
  title: 'Anatomy of logMhost',
  q1_proxy: {
    correlations: mhostCorrs.map(c => ({ name: c.name, r: +c.r.toFixed(3) })),
    mhost_alone_gap: +mhostGap.toFixed(1),
    is_proxy: mhostIsProxy
  },
  q2_form: {
    linear_aic: +linearAIC.toFixed(2),
    broken_aic: bestSplitAIC !== Infinity ? +bestSplitAIC.toFixed(2) : null,
    broken_split: bestSplitVal,
    quad_aic: +quadAIC.toFixed(2),
    quad_coeff: +quadFit.beta[2].toFixed(4),
    is_linear: isLinear
  },
  q3_residualization: {
    mhost_r2: +mhostFit.r2.toFixed(4),
    residual_sd: +sdYresid.toFixed(4),
    mhi_mr_resid_gap: +gapMhiMRresid.toFixed(1),
    mhi_mr_independent: mhiMRindependent,
    competitors: competitors.map(c => {
      const Xorig = gals45.map((_, i) => c.vars.map(v => allVarArrays[v][i]));
      const XwithMhost = gals45.map((_, i) => [logMhost[i], ...c.vars.map(v => allVarArrays[v][i])]);
      return {
        name: c.name,
        gap_alone: +gapPct(looCV(Y, Xorig), sdY).toFixed(1),
        gap_with_mhost: +gapPct(looCV(Y, XwithMhost), sdY).toFixed(1)
      };
    })
  },
  q4_residuals: {
    universal_sd: +sdY.toFixed(4),
    mhost_only_sd: +sd(mhostOnlyResid).toFixed(4),
    m3_sd: +sd(m3Resid).toFixed(4),
    top_residual_correlations: residCorrs.slice(0, 5).map(c => ({ name: c.name, r: +c.r.toFixed(3) }))
  },
  q5_benefiters: {
    n_benefit: benefiters.length,
    n_total: N,
    top10: improvement.slice(0, 10).map(g => ({
      name: g.name, delta: +g.delta.toFixed(3),
      logMhost: g.logMhost, Vflat: g.Vflat, envCode: g.envCode
    }))
  },
  q6_sources: Object.entries(sourceGroups).map(([src, names]) => ({ source: src, n: names.length })),
  synthesis: {
    mhost_is_proxy: mhostIsProxy,
    form_is_linear: isLinear,
    mhi_mr_independent_of_mhost: mhiMRindependent,
    m3_residuals_clean: m3ResidClean
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase124-logMhost-anatomy.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase124-logMhost-anatomy.json');

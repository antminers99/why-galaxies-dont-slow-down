#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 126: M4 CANDIDATE — M3 + logSigma0');
console.log('');
console.log('  M3: logMHI + logMhost + logMeanRun');
console.log('  M4: logMHI + logMhost + logMeanRun + logSigma0');
console.log('');
console.log('  Question: Does logSigma0 add a genuine, stable fourth axis?');
console.log('');
console.log('  Verdict criteria:');
console.log('    ADOPTED  — nested CV improvement + stable signs');
console.log('    EXTENDED — clear LOO gain but nested selection unstable');
console.log('    REJECTED — fragile or negligible gain');
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
  const aic = n * Math.log(rss / n) + 2 * p;
  const bic = n * Math.log(rss / n) + p * Math.log(n);
  return { beta, resid, rss, tss, r2, r2adj, aic, bic, n, k: p };
}

function looCV(Y, X) {
  const n = Y.length; let ss = 0;
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

function mulberry32(a) {
  return function () {
    a |= 0; a = a + 0x6D2B79F5 | 0;
    var t = Math.imul(a ^ a >>> 15, 1 | a);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

function shuffle(arr, rng) {
  const a = [...arr]; for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1)); [a[i], a[j]] = [a[j], a[i]];
  } return a;
}

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logMHI = gals45.map(g => g.logMHI);
const logMR = gals45.map(g => g.logMeanRun);
const logSig0 = gals45.map(g => g.logSigma0);
const logVflat = gals45.map(g => Math.log10(sparcMap[g.name]?.Vflat || 150));
const logL36 = gals45.map(g => Math.log10(Math.max(sparcMap[g.name]?.L36 || 1, 0.01)));
const logRdisk = gals45.map(g => Math.log10(sparcMap[g.name]?.Rdisk || 3));
const morphT = gals45.map(g => sparcMap[g.name]?.T ?? 5);
const rcWig = gals45.map(g => g.rcWiggliness);
const envCode = gals45.map(g => g.envCode);

const m3X = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]);
const m4X = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], logSig0[i]]);

const m3Names = ['logMHI', 'logMhost', 'logMR'];
const m4Names = ['logMHI', 'logMhost', 'logMR', 'logSigma0'];

console.log('  N = ' + N + ', SD(logA0) = ' + sdY.toFixed(4) + ' dex\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: HEAD-TO-HEAD — M3 vs M4');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const m3Fit = ols(Y, m3X);
const m4Fit = ols(Y, m4X);
const m3Loo = looCV(Y, m3X);
const m4Loo = looCV(Y, m4X);
const m3Gap = gapPct(m3Loo, sdY);
const m4Gap = gapPct(m4Loo, sdY);

function printModel(name, fit, loo, gap, varNames) {
  console.log('  ' + name + ':');
  console.log('    Coefficients:');
  console.log('      intercept = ' + fit.beta[0].toFixed(4));
  for (let j = 0; j < varNames.length; j++) {
    console.log('      ' + varNames[j].padEnd(12) + ' = ' + fit.beta[j + 1].toFixed(4));
  }
  console.log('    LOO gap%:  ' + gap.toFixed(1) + '%');
  console.log('    LOO RMS:   ' + loo.toFixed(4) + ' dex');
  console.log('    R2:        ' + fit.r2.toFixed(4));
  console.log('    R2adj:     ' + fit.r2adj.toFixed(4));
  console.log('    AIC:       ' + fit.aic.toFixed(2));
  console.log('    BIC:       ' + fit.bic.toFixed(2));
  console.log('    Resid SD:  ' + sd(fit.resid).toFixed(4) + ' dex');
  console.log('');
}

printModel('M3 (3-axis)', m3Fit, m3Loo, m3Gap, m3Names);
printModel('M4 (4-axis)', m4Fit, m4Loo, m4Gap, m4Names);

console.log('  Comparison:');
console.log('    Delta gap:    ' + (m4Gap - m3Gap > 0 ? '+' : '') + (m4Gap - m3Gap).toFixed(1) + 'pp');
console.log('    Delta RMS:    ' + (m4Loo - m3Loo > 0 ? '+' : '') + ((m4Loo - m3Loo) * 1000).toFixed(1) + ' millidex');
console.log('    Delta AIC:    ' + (m3Fit.aic - m4Fit.aic).toFixed(2) + ' (positive = M4 preferred)');
console.log('    Delta BIC:    ' + (m3Fit.bic - m4Fit.bic).toFixed(2) + ' (positive = M4 preferred)');
console.log('    Delta R2adj:  ' + (m4Fit.r2adj - m3Fit.r2adj > 0 ? '+' : '') + (m4Fit.r2adj - m3Fit.r2adj).toFixed(4));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: BOOTSTRAP STABILITY');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const NBOOT = 2000;
const rng = mulberry32(42);

function fullBootstrap(Y, X, varNames, nBoot, rngFn) {
  const n = Y.length, p = X[0].length;
  const baseFit = ols(Y, X);
  const baseSigns = baseFit.beta.slice(1).map(b => Math.sign(b));
  const flips = new Array(p).fill(0);
  const coeffSamples = Array.from({ length: p + 1 }, () => []);

  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(rngFn() * n));
    const Yb = idx.map(i => Y[i]); const Xb = idx.map(i => X[i]);
    const fit = ols(Yb, Xb);
    for (let j = 0; j <= p; j++) coeffSamples[j].push(fit.beta[j]);
    for (let j = 0; j < p; j++) {
      if (Math.sign(fit.beta[j + 1]) !== baseSigns[j]) flips[j]++;
    }
  }

  coeffSamples.forEach(s => s.sort((a, b) => a - b));
  return {
    flipRates: flips.map(f => f / nBoot),
    medians: coeffSamples.map(s => s[Math.floor(s.length / 2)]),
    ci025: coeffSamples.map(s => s[Math.floor(s.length * 0.025)]),
    ci975: coeffSamples.map(s => s[Math.floor(s.length * 0.975)]),
    sds: coeffSamples.map(s => sd(s))
  };
}

const m3Boot = fullBootstrap(Y, m3X, m3Names, NBOOT, rng);
const m4Boot = fullBootstrap(Y, m4X, m4Names, NBOOT, mulberry32(42));

console.log('  M3 bootstrap (N=' + NBOOT + '):\n');
for (let j = 0; j < m3Names.length; j++) {
  console.log('    ' + m3Names[j].padEnd(12) + ': ' +
    m3Boot.medians[j + 1].toFixed(4) + ' [' + m3Boot.ci025[j + 1].toFixed(4) + ', ' + m3Boot.ci975[j + 1].toFixed(4) + ']' +
    '  flip=' + (m3Boot.flipRates[j] * 100).toFixed(1) + '%');
}
console.log('');

console.log('  M4 bootstrap (N=' + NBOOT + '):\n');
for (let j = 0; j < m4Names.length; j++) {
  console.log('    ' + m4Names[j].padEnd(12) + ': ' +
    m4Boot.medians[j + 1].toFixed(4) + ' [' + m4Boot.ci025[j + 1].toFixed(4) + ', ' + m4Boot.ci975[j + 1].toFixed(4) + ']' +
    '  flip=' + (m4Boot.flipRates[j] * 100).toFixed(1) + '%');
}
const m4AllStable = m4Boot.flipRates.every(f => f < 0.10);
console.log('    All signs stable (< 10%)? ' + (m4AllStable ? 'YES ✓' : 'NO ✗'));
console.log('');

console.log('  M3 coefficient stability across LOO folds:\n');
const m3LooCoeffs = Array.from({ length: 4 }, () => []);
for (let i = 0; i < N; i++) {
  const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
  const Xt = [...m3X.slice(0, i), ...m3X.slice(i + 1)];
  const f = ols(Yt, Xt);
  for (let j = 0; j <= 3; j++) m3LooCoeffs[j].push(f.beta[j]);
}
console.log('    ' + 'Variable'.padEnd(12) + '  mean       SD         range');
for (let j = 0; j <= 3; j++) {
  const label = j === 0 ? 'intercept' : m3Names[j - 1];
  const vals = m3LooCoeffs[j];
  console.log('    ' + label.padEnd(12) + '  ' + mean(vals).toFixed(4).padStart(8) + '   ' +
    sd(vals).toFixed(4).padStart(8) + '   [' + Math.min(...vals).toFixed(4) + ', ' + Math.max(...vals).toFixed(4) + ']');
}
console.log('');

console.log('  M4 coefficient stability across LOO folds:\n');
const m4LooCoeffs = Array.from({ length: 5 }, () => []);
for (let i = 0; i < N; i++) {
  const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
  const Xt = [...m4X.slice(0, i), ...m4X.slice(i + 1)];
  const f = ols(Yt, Xt);
  for (let j = 0; j <= 4; j++) m4LooCoeffs[j].push(f.beta[j]);
}
console.log('    ' + 'Variable'.padEnd(12) + '  mean       SD         range');
const m4LooLabels = ['intercept', ...m4Names];
for (let j = 0; j <= 4; j++) {
  const vals = m4LooCoeffs[j];
  console.log('    ' + m4LooLabels[j].padEnd(12) + '  ' + mean(vals).toFixed(4).padStart(8) + '   ' +
    sd(vals).toFixed(4).padStart(8) + '   [' + Math.min(...vals).toFixed(4) + ', ' + Math.max(...vals).toFixed(4) + ']');
}
const m4LooSignsStable = m4Names.every((_, j) => {
  const vals = m4LooCoeffs[j + 1];
  const baseSign = Math.sign(m4Fit.beta[j + 1]);
  return vals.every(v => Math.sign(v) === baseSign);
});
console.log('    All LOO folds same sign? ' + (m4LooSignsStable ? 'YES ✓' : 'NO ✗'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: NESTED CV — Does M4 survive model selection?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const modelsNested = {
  M3: { X: m3X },
  M4: { X: m4X }
};
let nestedSS = 0;
const selCount = { M3: 0, M4: 0 };

for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  let bestGap = -Infinity, bestPred = 0, bestMod = '';

  for (const [mk, mod] of Object.entries(modelsNested)) {
    const Xfull = mod.X;
    const Xtrain = [...Xfull.slice(0, i), ...Xfull.slice(i + 1)];
    let innerSS = 0;
    for (let j = 0; j < Ytrain.length; j++) {
      const Yinn = [...Ytrain.slice(0, j), ...Ytrain.slice(j + 1)];
      const Xinn = [...Xtrain.slice(0, j), ...Xtrain.slice(j + 1)];
      const f = ols(Yinn, Xinn);
      const xi = [1, ...Xtrain[j]];
      innerSS += (Ytrain[j] - xi.reduce((s, x, jj) => s + x * f.beta[jj], 0)) ** 2;
    }
    const innerGap = gapPct(Math.sqrt(innerSS / Ytrain.length), sd(Ytrain));
    if (innerGap > bestGap) {
      bestGap = innerGap;
      bestMod = mk;
      const f = ols(Ytrain, Xtrain);
      const xi = [1, ...Xfull[i]];
      bestPred = xi.reduce((s, x, j) => s + x * f.beta[j], 0);
    }
  }
  selCount[bestMod]++;
  nestedSS += (Y[i] - bestPred) ** 2;
}

const nestedRms = Math.sqrt(nestedSS / N);
const nestedGap = gapPct(nestedRms, sdY);

console.log('  Nested CV (outer LOO, inner LOO for M3 vs M4):');
console.log('    Nested gap% = ' + nestedGap.toFixed(1) + '%');
console.log('    M3 selected: ' + selCount.M3 + '/' + N + ' (' + (selCount.M3 / N * 100).toFixed(0) + '%)');
console.log('    M4 selected: ' + selCount.M4 + '/' + N + ' (' + (selCount.M4 / N * 100).toFixed(0) + '%)');
console.log('    M3 standard gap: ' + m3Gap.toFixed(1) + '%');
console.log('    M4 standard gap: ' + m4Gap.toFixed(1) + '%');
console.log('    Nested optimism: ' + (Math.max(m3Gap, m4Gap) - nestedGap).toFixed(1) + 'pp');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: RESIDUAL STRUCTURE — Does Sigma0 clean up the residuals?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const residVars = [
  { name: 'logSigma0', arr: logSig0 },
  { name: 'logVflat', arr: logVflat },
  { name: 'logL36', arr: logL36 },
  { name: 'morphT', arr: morphT },
  { name: 'logRdisk', arr: logRdisk },
  { name: 'rcWig', arr: rcWig },
  { name: 'envCode', arr: envCode }
];

console.log('  Residual correlations:\n');
console.log('    ' + 'Variable'.padEnd(14) + '  r(M3 resid)  r(M4 resid)  change');
for (const rv of residVars) {
  const rM3 = pearsonR(m3Fit.resid, rv.arr);
  const rM4 = pearsonR(m4Fit.resid, rv.arr);
  const change = Math.abs(rM4) < Math.abs(rM3) ? 'reduced' : Math.abs(rM4) > Math.abs(rM3) + 0.02 ? 'INCREASED' : 'same';
  console.log('    ' + rv.name.padEnd(14) + '  ' +
    rM3.toFixed(3).padStart(10) + '  ' +
    rM4.toFixed(3).padStart(10) + '  ' + change);
}
console.log('');

const m3MaxResidR = Math.max(...residVars.map(rv => Math.abs(pearsonR(m3Fit.resid, rv.arr))));
const m4MaxResidR = Math.max(...residVars.map(rv => Math.abs(pearsonR(m4Fit.resid, rv.arr))));
console.log('  Max |residual correlation|:');
console.log('    M3: ' + m3MaxResidR.toFixed(3));
console.log('    M4: ' + m4MaxResidR.toFixed(3));
console.log('    Improvement: ' + ((1 - m4MaxResidR / m3MaxResidR) * 100).toFixed(0) + '%');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 5: COLLINEARITY CHECK — Is Sigma0 independent enough?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const sig0Corrs = [
  { name: 'logMHI', arr: logMHI },
  { name: 'logMhost', arr: logMhost },
  { name: 'logMR', arr: logMR },
  { name: 'logVflat', arr: logVflat },
  { name: 'logL36', arr: logL36 }
];

console.log('  Correlations of logSigma0 with M3/M4 predictors:\n');
for (const sc of sig0Corrs) {
  const r = pearsonR(logSig0, sc.arr);
  const flag = Math.abs(r) > 0.7 ? ' ★★★ DANGER' : Math.abs(r) > 0.5 ? ' ★★ HIGH' : Math.abs(r) > 0.3 ? ' ★' : '';
  console.log('    r(Sigma0, ' + sc.name.padEnd(10) + ') = ' + r.toFixed(3) + flag);
}
console.log('');

const sig0ResidM3pred = ols(logSig0, m3X);
console.log('  Sigma0 variance explained by M3 predictors: R2 = ' + sig0ResidM3pred.r2.toFixed(4));
console.log('  Sigma0 independent fraction: ' + ((1 - sig0ResidM3pred.r2) * 100).toFixed(1) + '%');
const rSig0Indep = pearsonR(sig0ResidM3pred.resid, Y);
console.log('  r(Sigma0_residual, logA0) = ' + rSig0Indep.toFixed(3) + ' — this is the unique info Sigma0 adds');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 6: PERMUTATION TEST — Is Sigma0 addition significant?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const NPERM = 3000;
const realDeltaGap = m4Gap - m3Gap;
const rngPerm = mulberry32(999);
let nBeatPerm = 0;

for (let p = 0; p < NPERM; p++) {
  const sig0Perm = shuffle(logSig0, rngPerm);
  const m4Perm = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], sig0Perm[i]]);
  const looPerm = looCV(Y, m4Perm);
  const gapPerm = gapPct(looPerm, sdY);
  const deltaPerm = gapPerm - m3Gap;
  if (deltaPerm >= realDeltaGap) nBeatPerm++;
}

const permP = (nBeatPerm + 1) / (NPERM + 1);
console.log('  Real delta gap (M4 - M3): ' + (realDeltaGap > 0 ? '+' : '') + realDeltaGap.toFixed(1) + 'pp');
console.log('  Permutations matching or beating: ' + nBeatPerm + '/' + NPERM);
console.log('  Permutation p = ' + permP.toFixed(4));
console.log('  Verdict: ' + (permP < 0.05 ? 'SIGNIFICANT (p < 0.05) ✓' : 'NOT SIGNIFICANT ✗'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 7: ALTERNATIVE 4-VAR MODELS — Is Sigma0 the best 4th axis?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const alt4th = [
  { name: 'logSigma0', arr: logSig0 },
  { name: 'logVflat', arr: logVflat },
  { name: 'logL36', arr: logL36 },
  { name: 'morphT', arr: morphT },
  { name: 'logRdisk', arr: logRdisk },
  { name: 'rcWig', arr: rcWig },
  { name: 'envCode', arr: envCode }
];

const alt4Results = [];
for (const a4 of alt4th) {
  const Xa = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], a4.arr[i]]);
  const fit = ols(Y, Xa);
  const loo = looCV(Y, Xa);
  const gap = gapPct(loo, sdY);
  alt4Results.push({ name: a4.name, gap, aic: fit.aic, bic: fit.bic, r2adj: fit.r2adj });
}

alt4Results.sort((a, b) => b.gap - a.gap);
console.log('  ┌────────────────┬──────────┬──────────┬──────────┐');
console.log('  │ 4th axis       │ LOO gap% │ AIC      │ R2adj    │');
console.log('  ├────────────────┼──────────┼──────────┼──────────┤');
for (const ar of alt4Results) {
  const marker = ar.name === 'logSigma0' ? ' ◄' : '';
  console.log('  │ ' + (ar.name + marker).padEnd(14) + ' │ ' +
    ar.gap.toFixed(1).padStart(6) + '% │ ' +
    ar.aic.toFixed(2).padStart(8) + ' │ ' +
    ar.r2adj.toFixed(4).padStart(8) + ' │');
}
console.log('  └────────────────┴──────────┴──────────┴──────────┘');
console.log('  M3 baseline: gap=' + m3Gap.toFixed(1) + '%, AIC=' + m3Fit.aic.toFixed(2));
console.log('');

const sig0Rank = alt4Results.findIndex(a => a.name === 'logSigma0') + 1;
console.log('  logSigma0 rank as 4th axis: #' + sig0Rank + '/' + alt4Results.length);
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 126: FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const deltaGapSignificant = m4Gap - m3Gap > 1.5;
const aicPrefers = m3Fit.aic - m4Fit.aic > 2;
const bicPrefers = m3Fit.bic - m4Fit.bic > 0;
const nestedM4Wins = selCount.M4 > selCount.M3;
const nestedGapImproved = nestedGap > m3Gap - 5;
const signsStable = m4AllStable;
const looSignsStable = m4LooSignsStable;
const residImproved = m4MaxResidR < m3MaxResidR - 0.03;
const permSignificant = permP < 0.05;
const bestFourthAxis = sig0Rank === 1;

console.log('  SCORECARD:\n');
console.log('    LOO gap improvement:       ' + (m4Gap - m3Gap).toFixed(1) + 'pp ' + (deltaGapSignificant ? '✓' : '—'));
console.log('    AIC prefers M4:            ' + (m3Fit.aic - m4Fit.aic).toFixed(2) + ' ' + (aicPrefers ? '✓' : '✗'));
console.log('    BIC prefers M4:            ' + (m3Fit.bic - m4Fit.bic).toFixed(2) + ' ' + (bicPrefers ? '✓' : '✗'));
console.log('    Nested CV M4 selected:     ' + selCount.M4 + '/' + N + ' ' + (nestedM4Wins ? '✓' : '✗'));
console.log('    Bootstrap signs stable:    ' + (signsStable ? '✓' : '✗'));
console.log('    LOO fold signs stable:     ' + (looSignsStable ? '✓' : '✗'));
console.log('    Residual structure reduced: ' + (residImproved ? '✓' : '✗'));
console.log('    Permutation significant:   p=' + permP.toFixed(4) + ' ' + (permSignificant ? '✓' : '✗'));
console.log('    Best 4th axis:             #' + sig0Rank + ' ' + (bestFourthAxis ? '✓' : '✗'));
console.log('');

const strongCriteria = [nestedM4Wins, signsStable, permSignificant].filter(Boolean).length;
const supportCriteria = [deltaGapSignificant, aicPrefers, residImproved, bestFourthAxis, looSignsStable].filter(Boolean).length;

let verdict;
if (strongCriteria >= 2 && supportCriteria >= 3) {
  verdict = 'ADOPTED';
} else if (strongCriteria >= 1 && supportCriteria >= 2) {
  verdict = 'EXTENDED';
} else {
  verdict = 'REJECTED';
}

console.log('  Strong criteria passed: ' + strongCriteria + '/3');
console.log('  Support criteria passed: ' + supportCriteria + '/5');
console.log('');
console.log('  ╔══════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + verdict.padEnd(49) + '║');
if (verdict === 'ADOPTED') {
  console.log('  ║  M4 is the new baseline state law.                         ║');
  console.log('  ║  logA0 = f(logMHI, logMhost, logMeanRun, logSigma0)        ║');
} else if (verdict === 'EXTENDED') {
  console.log('  ║  M4 improves on M3 but with insufficient nested stability. ║');
  console.log('  ║  M3 remains baseline; M4 is the extended descriptive model.║');
} else {
  console.log('  ║  Sigma0 does not add reliable predictive value.            ║');
  console.log('  ║  M3 remains the definitive state law.                      ║');
}
console.log('  ╚══════════════════════════════════════════════════════════════╝');
console.log('');

if (verdict !== 'REJECTED') {
  console.log('  M4 PHYSICAL INTERPRETATION:');
  console.log('    logMhost  = external gravitational environment (host halo depth)');
  console.log('    logMHI    = baryonic/gas evolutionary state');
  console.log('    logMeanRun = internal dynamical coherence');
  console.log('    logSigma0 = baryon surface concentration/structure');
  console.log('');
  console.log('    Picture: Environment + Fuel + Dynamics + Concentration → a0');
  console.log('');
}

const output = {
  phase: '126',
  title: 'M4 Candidate — M3 + logSigma0',
  m3: {
    formula: 'logA0 = ' + m3Fit.beta.map((b, j) => (j === 0 ? '' : ' + ') + b.toFixed(4) + (j === 0 ? '' : '*' + m3Names[j - 1])).join(''),
    gap: +m3Gap.toFixed(1), loo_rms: +m3Loo.toFixed(4), r2adj: +m3Fit.r2adj.toFixed(4),
    aic: +m3Fit.aic.toFixed(2), bic: +m3Fit.bic.toFixed(2),
    coefficients: Object.fromEntries([['intercept', +m3Fit.beta[0].toFixed(4)], ...m3Names.map((v, j) => [v, +m3Fit.beta[j + 1].toFixed(4)])])
  },
  m4: {
    formula: 'logA0 = ' + m4Fit.beta.map((b, j) => (j === 0 ? '' : ' + ') + b.toFixed(4) + (j === 0 ? '' : '*' + m4Names[j - 1])).join(''),
    gap: +m4Gap.toFixed(1), loo_rms: +m4Loo.toFixed(4), r2adj: +m4Fit.r2adj.toFixed(4),
    aic: +m4Fit.aic.toFixed(2), bic: +m4Fit.bic.toFixed(2),
    coefficients: Object.fromEntries([['intercept', +m4Fit.beta[0].toFixed(4)], ...m4Names.map((v, j) => [v, +m4Fit.beta[j + 1].toFixed(4)])])
  },
  comparison: {
    delta_gap: +(m4Gap - m3Gap).toFixed(1),
    delta_aic: +(m3Fit.aic - m4Fit.aic).toFixed(2),
    delta_bic: +(m3Fit.bic - m4Fit.bic).toFixed(2)
  },
  bootstrap: {
    m4_flip_rates: Object.fromEntries(m4Names.map((v, j) => [v, +(m4Boot.flipRates[j] * 100).toFixed(1)])),
    all_stable: m4AllStable
  },
  loo_stability: {
    all_signs_consistent: m4LooSignsStable
  },
  nested_cv: {
    gap: +nestedGap.toFixed(1),
    m3_selected: selCount.M3, m4_selected: selCount.M4,
    m4_wins: nestedM4Wins
  },
  permutation: { p: +permP.toFixed(4), significant: permSignificant, n_perms: NPERM },
  fourth_axis_ranking: alt4Results.map(a => ({ name: a.name, gap: +a.gap.toFixed(1) })),
  sigma0_independence: {
    r2_explained_by_m3_predictors: +sig0ResidM3pred.r2.toFixed(4),
    independent_fraction: +((1 - sig0ResidM3pred.r2) * 100).toFixed(1),
    r_unique_with_logA0: +rSig0Indep.toFixed(3)
  },
  verdict: {
    result: verdict,
    strong_criteria: strongCriteria,
    support_criteria: supportCriteria
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase126-m4-candidate.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase126-m4-candidate.json');

#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 127: VFLAT CHALLENGE AGAINST M4');
console.log('');
console.log('  M4 (adopted): logMHI + logMhost + logMR + logSigma0');
console.log('  Challenger:   logVflat gave 52.1% gap as 4th axis on M3');
console.log('');
console.log('  Question: Is logVflat independently real, or a proxy?');
console.log('');
console.log('  Models:');
console.log('    A) M4:           MHI + Mhost + MR + Sigma0');
console.log('    B) M4+Vflat:     MHI + Mhost + MR + Sigma0 + Vflat');
console.log('    C) M3+Vflat:     MHI + Mhost + MR + Vflat (replace Sigma0)');
console.log('    D) MHI+Mhost+Sigma0+Vflat (replace MR)');
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

function bootstrapSigns(Y, X, nBoot, seed) {
  const rng = mulberry32(seed);
  const n = Y.length, p = X[0].length;
  const baseFit = ols(Y, X);
  const baseSigns = baseFit.beta.slice(1).map(b => Math.sign(b));
  const flips = new Array(p).fill(0);
  const coeffs = Array.from({ length: p + 1 }, () => []);
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(rng() * n));
    const Yb = idx.map(i => Y[i]); const Xb = idx.map(i => X[i]);
    const fit = ols(Yb, Xb);
    for (let j = 0; j <= p; j++) coeffs[j].push(fit.beta[j]);
    for (let j = 0; j < p; j++) if (Math.sign(fit.beta[j + 1]) !== baseSigns[j]) flips[j]++;
  }
  coeffs.forEach(s => s.sort((a, b) => a - b));
  return {
    flipRates: flips.map(f => f / nBoot),
    ci025: coeffs.map(s => s[Math.floor(s.length * 0.025)]),
    ci975: coeffs.map(s => s[Math.floor(s.length * 0.975)])
  };
}

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logMHI = gals45.map(g => g.logMHI);
const logMR = gals45.map(g => g.logMeanRun);
const logSig0 = gals45.map(g => g.logSigma0);
const Vflat = gals45.map(g => sparcMap[g.name]?.Vflat || 150);
const logVflat = Vflat.map(v => Math.log10(v));

const models = {
  A: { name: 'M4 (MHI+Mhost+MR+Sig0)', vars: ['logMHI', 'logMhost', 'logMR', 'logSig0'],
       X: gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], logSig0[i]]) },
  B: { name: 'M4+Vflat (5-var)', vars: ['logMHI', 'logMhost', 'logMR', 'logSig0', 'logVflat'],
       X: gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], logSig0[i], logVflat[i]]) },
  C: { name: 'M3+Vflat (replace Sig0)', vars: ['logMHI', 'logMhost', 'logMR', 'logVflat'],
       X: gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], logVflat[i]]) },
  D: { name: 'MHI+Mhost+Sig0+Vflat (replace MR)', vars: ['logMHI', 'logMhost', 'logSig0', 'logVflat'],
       X: gals45.map((_, i) => [logMHI[i], logMhost[i], logSig0[i], logVflat[i]]) }
};

console.log('  N = ' + N + ', SD(logA0) = ' + sdY.toFixed(4) + ' dex\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 1: HEAD-TO-HEAD COMPARISON');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const NBOOT = 2000;
const results = {};

for (const [key, mod] of Object.entries(models)) {
  const fit = ols(Y, mod.X);
  const loo = looCV(Y, mod.X);
  const gap = gapPct(loo, sdY);
  const boot = bootstrapSigns(Y, mod.X, NBOOT, 42);

  results[key] = { fit, loo, gap, boot, name: mod.name, vars: mod.vars };

  console.log('  Model ' + key + ': ' + mod.name);
  console.log('    Coefficients:');
  console.log('      intercept = ' + fit.beta[0].toFixed(4));
  for (let j = 0; j < mod.vars.length; j++) {
    console.log('      ' + mod.vars[j].padEnd(12) + ' = ' + fit.beta[j + 1].toFixed(4) +
      '  [' + boot.ci025[j + 1].toFixed(4) + ', ' + boot.ci975[j + 1].toFixed(4) + ']' +
      '  flip=' + (boot.flipRates[j] * 100).toFixed(1) + '%');
  }
  console.log('    LOO gap%:  ' + gap.toFixed(1) + '%');
  console.log('    LOO RMS:   ' + loo.toFixed(4) + ' dex');
  console.log('    R2adj:     ' + fit.r2adj.toFixed(4));
  console.log('    AIC:       ' + fit.aic.toFixed(2));
  console.log('    BIC:       ' + fit.bic.toFixed(2));
  console.log('    Resid SD:  ' + sd(fit.resid).toFixed(4) + ' dex');
  console.log('');
}

console.log('  ┌───────┬──────────────────────────────────────┬────────┬──────────┬──────────┬──────────┐');
console.log('  │ Model │ Description                          │ LOO gap│ AIC      │ BIC      │ k        │');
console.log('  ├───────┼──────────────────────────────────────┼────────┼──────────┼──────────┼──────────┤');
for (const [key, mr] of Object.entries(results)) {
  console.log('  │   ' + key + '   │ ' + mr.name.padEnd(36) + ' │ ' +
    mr.gap.toFixed(1).padStart(5) + '% │ ' +
    mr.fit.aic.toFixed(2).padStart(8) + ' │ ' +
    mr.fit.bic.toFixed(2).padStart(8) + ' │ ' +
    String(mr.vars.length).padStart(8) + ' │');
}
console.log('  └───────┴──────────────────────────────────────┴────────┴──────────┴──────────┴──────────┘');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 2: WHAT DOES VFLAT ABSORB IN MODEL B (5-var)?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const m4Coeffs = results.A.fit.beta;
const m5Coeffs = results.B.fit.beta;
const m4Vars = ['intercept', ...models.A.vars];
const m5Vars = ['intercept', ...models.B.vars];

console.log('  Coefficient comparison (M4 → M4+Vflat):\n');
console.log('    ' + 'Variable'.padEnd(14) + '  M4          M4+Vflat     change');
for (let j = 0; j < m4Vars.length; j++) {
  const c4 = m4Coeffs[j];
  const c5 = m5Coeffs[j];
  const pctChange = c4 !== 0 ? ((c5 - c4) / Math.abs(c4) * 100) : 0;
  const flag = Math.abs(pctChange) > 50 ? ' ★★★ ABSORBED' : Math.abs(pctChange) > 25 ? ' ★★ SHIFTED' : '';
  console.log('    ' + m4Vars[j].padEnd(14) + '  ' + c4.toFixed(4).padStart(8) + '     ' +
    c5.toFixed(4).padStart(8) + '     ' + (pctChange > 0 ? '+' : '') + pctChange.toFixed(0) + '%' + flag);
}
console.log('    logVflat                          ' + m5Coeffs[5].toFixed(4).padStart(8) + '     (new)');
console.log('');

const vflatBoot = results.B.boot;
const vflatFlip = vflatBoot.flipRates[4];
const vflatCI = [vflatBoot.ci025[5], vflatBoot.ci975[5]];
const crossesZero = vflatCI[0] * vflatCI[1] < 0;
console.log('  Vflat in 5-var model:');
console.log('    Coefficient: ' + m5Coeffs[5].toFixed(4));
console.log('    95% CI: [' + vflatCI[0].toFixed(4) + ', ' + vflatCI[1].toFixed(4) + ']');
console.log('    Crosses zero: ' + (crossesZero ? 'YES ✗ — not robustly signed' : 'NO ✓'));
console.log('    Sign flip rate: ' + (vflatFlip * 100).toFixed(1) + '%');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 3: COLLINEARITY AUDIT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const allVars = [
  { name: 'logMHI', arr: logMHI },
  { name: 'logMhost', arr: logMhost },
  { name: 'logMR', arr: logMR },
  { name: 'logSig0', arr: logSig0 },
  { name: 'logVflat', arr: logVflat }
];

console.log('  Pairwise correlations:\n');
for (let i = 0; i < allVars.length; i++) {
  for (let j = i + 1; j < allVars.length; j++) {
    const r = pearsonR(allVars[i].arr, allVars[j].arr);
    const flag = Math.abs(r) > 0.7 ? ' ★★★ DANGER' : Math.abs(r) > 0.5 ? ' ★★ HIGH' : Math.abs(r) > 0.3 ? ' ★' : '';
    console.log('    r(' + allVars[i].name + ', ' + allVars[j].name + ') = ' + r.toFixed(3) + flag);
  }
}
console.log('');

const vflatResidM4 = ols(logVflat, models.A.X);
console.log('  Vflat variance explained by M4 predictors: R2 = ' + vflatResidM4.r2.toFixed(4));
console.log('  Vflat independent fraction: ' + ((1 - vflatResidM4.r2) * 100).toFixed(1) + '%');
const rVflatIndep = pearsonR(vflatResidM4.resid, Y);
console.log('  r(Vflat_residual, logA0) = ' + rVflatIndep.toFixed(3));
console.log('  → This is the UNIQUE info Vflat adds beyond M4');
console.log('');

const sig0ResidM3Vflat = ols(logSig0, models.C.X);
console.log('  Sigma0 variance explained by M3+Vflat: R2 = ' + sig0ResidM3Vflat.r2.toFixed(4));
const mrResidM3Vflat = ols(logMR, gals45.map((_, i) => [logMHI[i], logMhost[i], logVflat[i]]));
console.log('  MeanRun variance explained by MHI+Mhost+Vflat: R2 = ' + mrResidM3Vflat.r2.toFixed(4));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 4: NESTED CV — M4 vs C vs D');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const nestedModels = { A: models.A, C: models.C, D: models.D };
let nestedSS = 0;
const selCount = { A: 0, C: 0, D: 0 };

for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  let bestGap = -Infinity, bestPred = 0, bestMod = '';

  for (const [mk, mod] of Object.entries(nestedModels)) {
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

console.log('  Nested CV (outer LOO, inner LOO for A vs C vs D):');
console.log('    Nested gap% = ' + nestedGap.toFixed(1) + '%');
console.log('    Selection frequency:');
for (const mk of ['A', 'C', 'D']) {
  console.log('      ' + mk + ' (' + nestedModels[mk].name + '): ' + selCount[mk] + '/' + N + ' (' + (selCount[mk] / N * 100).toFixed(0) + '%)');
}
console.log('');

console.log('  Nested CV — M4 vs M4+Vflat (5-var):');
const sel2 = { A: 0, B: 0 };
let nested2SS = 0;
for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  let bestGap = -Infinity, bestPred = 0, bestMod = '';
  for (const mk of ['A', 'B']) {
    const Xfull = models[mk].X;
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
  sel2[bestMod]++;
  nested2SS += (Y[i] - bestPred) ** 2;
}
const nested2Gap = gapPct(Math.sqrt(nested2SS / N), sdY);
console.log('    Nested gap% = ' + nested2Gap.toFixed(1) + '%');
console.log('    M4 selected: ' + sel2.A + '/' + N + ' (' + (sel2.A / N * 100).toFixed(0) + '%)');
console.log('    M4+Vflat selected: ' + sel2.B + '/' + N + ' (' + (sel2.B / N * 100).toFixed(0) + '%)');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 5: RESIDUAL STRUCTURE COMPARISON');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const residVars = [
  { name: 'logVflat', arr: logVflat },
  { name: 'logSig0', arr: logSig0 },
  { name: 'logMR', arr: logMR },
  { name: 'logMHI', arr: logMHI },
  { name: 'logMhost', arr: logMhost },
  { name: 'morphT', arr: gals45.map(g => sparcMap[g.name]?.T ?? 5) },
  { name: 'logRdisk', arr: gals45.map(g => Math.log10(sparcMap[g.name]?.Rdisk || 3)) },
  { name: 'rcWig', arr: gals45.map(g => g.rcWiggliness) },
  { name: 'envCode', arr: gals45.map(g => g.envCode) }
];

console.log('  Residual correlations comparison:\n');
console.log('    ' + 'Variable'.padEnd(12) + '  M4(A)       M3+Vflat(C)  MHI+Mhost+Sig0+Vflat(D)');
for (const rv of residVars) {
  const rA = pearsonR(results.A.fit.resid, rv.arr);
  const rC = pearsonR(results.C.fit.resid, rv.arr);
  const rD = pearsonR(results.D.fit.resid, rv.arr);
  console.log('    ' + rv.name.padEnd(12) + '  ' +
    rA.toFixed(3).padStart(8) + '     ' +
    rC.toFixed(3).padStart(8) + '     ' +
    rD.toFixed(3).padStart(8));
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 6: CIRCULARITY / PROXY CHECK');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  logA0 is derived from rotation curve fitting.');
console.log('  Vflat is the asymptotic flat rotation velocity.');
console.log('  Both are extracted from the same velocity field.\n');

const rDirectVflat = pearsonR(logVflat, Y);
console.log('  Direct r(logVflat, logA0) = ' + rDirectVflat.toFixed(3));
console.log('  (By comparison: r(logSig0, logA0) = ' + pearsonR(logSig0, Y).toFixed(3) + ')');
console.log('');

console.log('  If Vflat were a pure kinematic proxy, we would expect:');
console.log('    1. High direct correlation with logA0 → r = ' + rDirectVflat.toFixed(3) +
  (Math.abs(rDirectVflat) < 0.3 ? ' — actually low, suggesting indirect role' : ' — moderate'));
console.log('    2. It would absorb MeanRun (both kinematic) → check Model D');

const mrCoeffD = results.D.fit.beta;
const mrInC = results.C.fit.beta;
console.log('');
console.log('  In Model C (M3+Vflat), logMR coefficient: ' + mrInC[3].toFixed(4));
console.log('  In Model D (no MR, has Vflat), are remaining coefficients stable?');
const m4Coeffs2 = results.A.fit.beta;
console.log('    logMHI:  M4=' + m4Coeffs2[1].toFixed(4) + ' → D=' + mrCoeffD[1].toFixed(4) + ' (change: ' + ((mrCoeffD[1] - m4Coeffs2[1]) / Math.abs(m4Coeffs2[1]) * 100).toFixed(0) + '%)');
console.log('    logMhost: M4=' + m4Coeffs2[2].toFixed(4) + ' → D=' + mrCoeffD[2].toFixed(4) + ' (change: ' + ((mrCoeffD[2] - m4Coeffs2[2]) / Math.abs(m4Coeffs2[2]) * 100).toFixed(0) + '%)');
console.log('');

const rVflatMR = pearsonR(logVflat, logMR);
console.log('  r(logVflat, logMeanRun) = ' + rVflatMR.toFixed(3));
console.log('  → ' + (Math.abs(rVflatMR) > 0.5 ? 'Vflat and MR share substantial variance — potential redundancy' :
  Math.abs(rVflatMR) > 0.3 ? 'Moderate overlap between Vflat and MR' :
  'Vflat and MR are largely independent'));
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 127: FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const cBeatsA = results.C.gap > results.A.gap;
const cNestedWins = selCount.C > selCount.A;
const cStable = results.C.boot.flipRates.every(f => f < 0.10);
const bAddsOverA = results.B.gap > results.A.gap + 2;
const bNestedWins = sel2.B > sel2.A;
const vflatStableInB = vflatFlip < 0.10 && !crossesZero;

console.log('  QUESTION 1: Does Vflat REPLACE Sigma0 (Model C vs A)?');
console.log('    C gap: ' + results.C.gap.toFixed(1) + '% vs A gap: ' + results.A.gap.toFixed(1) + '%');
console.log('    C nested selection: ' + selCount.C + '/' + N);
console.log('    C signs stable: ' + (cStable ? 'YES' : 'NO'));
console.log('    Verdict: ' + (cBeatsA && cNestedWins ? 'YES — Vflat can replace Sigma0' :
  cBeatsA ? 'PARTIAL — higher gap but nested selection unclear' :
  'NO — Sigma0 remains preferred'));
console.log('');

console.log('  QUESTION 2: Does Vflat ADD to M4 (Model B vs A)?');
console.log('    B gap: ' + results.B.gap.toFixed(1) + '% vs A gap: ' + results.A.gap.toFixed(1) + '%');
console.log('    B nested selection: ' + sel2.B + '/' + N);
console.log('    Vflat CI crosses zero in B: ' + (crossesZero ? 'YES ✗' : 'NO ✓'));
console.log('    Vflat sign flip in B: ' + (vflatFlip * 100).toFixed(1) + '%');
console.log('    Verdict: ' + (bAddsOverA && vflatStableInB ? 'YES — Vflat adds independent value over M4' :
  'NO — Vflat does not robustly add to M4'));
console.log('');

console.log('  QUESTION 3: Is Vflat a kinematic proxy (circularity concern)?');
console.log('    r(Vflat, logA0) = ' + rDirectVflat.toFixed(3));
console.log('    Vflat independent of M4 predictors: ' + ((1 - vflatResidM4.r2) * 100).toFixed(1) + '%');
console.log('    r(Vflat_resid, logA0) = ' + rVflatIndep.toFixed(3));
console.log('    Verdict: ' + (Math.abs(rDirectVflat) > 0.5 ? 'HIGH RISK — strong direct coupling' :
  Math.abs(rVflatIndep) > 0.3 ? 'MODERATE RISK — but carries unique info' :
  'LOW RISK — Vflat contribution is indirect'));
console.log('');

let finalVerdict;
if (cBeatsA && cNestedWins && cStable) {
  finalVerdict = 'RIVAL — Model C (M3+Vflat) is a co-equal alternative to M4';
} else if (cBeatsA && !cNestedWins) {
  finalVerdict = 'BENCHMARK — Vflat version has higher raw gap but M4 wins nested selection';
} else if (!cBeatsA) {
  finalVerdict = 'DISMISSED — Vflat does not outperform M4; M4 is definitive';
} else {
  finalVerdict = 'INCONCLUSIVE';
}

console.log('  ╔══════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + finalVerdict.substring(0, 49).padEnd(49) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════╝');
console.log('');

if (finalVerdict.startsWith('BENCHMARK') || finalVerdict.startsWith('DISMISSED')) {
  console.log('  M4 CONFIRMED as the definitive adopted state law:');
  console.log('  logA0 = ' + results.A.fit.beta[0].toFixed(4) +
    ' + ' + results.A.fit.beta[1].toFixed(4) + '*logMHI' +
    ' + ' + results.A.fit.beta[2].toFixed(4) + '*logMhost' +
    ' + ' + results.A.fit.beta[3].toFixed(4) + '*logMeanRun' +
    ' + ' + results.A.fit.beta[4].toFixed(4) + '*logSigma0');
  console.log('');
}

const output = {
  phase: '127',
  title: 'Vflat Challenge Against M4',
  models: {},
  collinearity: {
    r_vflat_mr: +rVflatMR.toFixed(3),
    r_vflat_logA0: +rDirectVflat.toFixed(3),
    vflat_r2_by_m4: +vflatResidM4.r2.toFixed(4),
    vflat_independent_pct: +((1 - vflatResidM4.r2) * 100).toFixed(1),
    r_vflat_resid_logA0: +rVflatIndep.toFixed(3)
  },
  nested_cv_acd: {
    gap: +nestedGap.toFixed(1),
    selection: selCount
  },
  nested_cv_ab: {
    gap: +nested2Gap.toFixed(1),
    selection: sel2
  },
  vflat_in_5var: {
    coefficient: +m5Coeffs[5].toFixed(4),
    ci95: vflatCI.map(v => +v.toFixed(4)),
    crosses_zero: crossesZero,
    flip_rate: +(vflatFlip * 100).toFixed(1)
  },
  verdict: {
    result: finalVerdict,
    c_beats_a: cBeatsA,
    c_nested_wins: cNestedWins,
    b_adds_over_a: bAddsOverA,
    vflat_stable_in_b: vflatStableInB
  }
};

for (const [key, mr] of Object.entries(results)) {
  output.models[key] = {
    name: mr.name,
    gap: +mr.gap.toFixed(1),
    loo_rms: +mr.loo.toFixed(4),
    r2adj: +mr.fit.r2adj.toFixed(4),
    aic: +mr.fit.aic.toFixed(2),
    bic: +mr.fit.bic.toFixed(2),
    coefficients: Object.fromEntries(
      [['intercept', +mr.fit.beta[0].toFixed(4)],
      ...mr.vars.map((v, j) => [v, +mr.fit.beta[j + 1].toFixed(4)])]
    ),
    bootstrap_flip_rates: Object.fromEntries(mr.vars.map((v, j) => [v, +(mr.boot.flipRates[j] * 100).toFixed(1)]))
  };
}

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase127-vflat-challenge.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase127-vflat-challenge.json');

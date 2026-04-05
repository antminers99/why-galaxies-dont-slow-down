#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 128: DECODE THE 4TH AXIS (PATCHED)');
console.log('');
console.log('  Fixes applied:');
console.log('    - Bug 1: logMbar gas mass units corrected (×1e9)');
console.log('    - Bug 2: LOO with _perp variables now residualizes INSIDE fold');
console.log('    - Verdict thresholds use functional criteria (gap impact)');
console.log('');
console.log('  Core question: Is Vflat an independent physical axis,');
console.log('  or a kinematic super-proxy compressing Sigma0 + other structure?');
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

function vif(X, idx) {
  const n = X.length;
  const Yv = X.map(r => r[idx]);
  const Xv = X.map(r => r.filter((_, j) => j !== idx));
  const fit = ols(Yv, Xv);
  return fit.r2 < 1 ? 1 / (1 - fit.r2) : Infinity;
}

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logMHI = gals45.map(g => g.logMHI);
const logMR = gals45.map(g => g.logMeanRun);
const logSig0 = gals45.map(g => g.logSigma0);
const logVflat = gals45.map(g => Math.log10(sparcMap[g.name].Vflat));
const logL36 = gals45.map(g => Math.log10(Math.max(sparcMap[g.name].L36, 0.01)));
const logRdisk = gals45.map(g => Math.log10(sparcMap[g.name].Rdisk));
const morphT = gals45.map(g => sparcMap[g.name].T);
const logMbar = gals45.map(g => {
  const s = sparcMap[g.name];
  const mStar = s.L36 * 0.5 * 1e9;
  const mGas = Math.pow(10, g.logMHI) * 1.33 * 1e9;
  return Math.log10(mStar + mGas);
});

const m4X = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], logSig0[i]]);
const mCX = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], logVflat[i]]);

console.log('  N = ' + N + ', SD(logA0) = ' + sdY.toFixed(4) + ' dex\n');

const m4Fit = ols(Y, m4X);
const mCFit = ols(Y, mCX);
const m4Loo = looCV(Y, m4X);
const mCLoo = looCV(Y, mCX);
const m4Gap = gapPct(m4Loo, sdY);
const mCGap = gapPct(mCLoo, sdY);

console.log('  Reference gaps:');
console.log('    M4: ' + m4Gap.toFixed(1) + '%    Model C: ' + mCGap.toFixed(1) + '%\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST A: VFLAT_PERP — Orthogonalized Vflat beyond M4');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  Global residualization (diagnostic only — NOT for LOO):');
const vflatPerp = ols(logVflat, m4X).resid;
const rVflatPerpA0 = pearsonR(vflatPerp, Y);
const rVflatPerpM4resid = pearsonR(vflatPerp, m4Fit.resid);
console.log('    r(Vflat_perp, logA0) = ' + rVflatPerpA0.toFixed(3));
console.log('    r(Vflat_perp, M4_residuals) = ' + rVflatPerpM4resid.toFixed(3));
console.log('');

console.log('  Proper LOO with fold-internal residualization:');
console.log('  (Vflat_perp recomputed on training set only in each fold)\n');

let m4PlusVperpSS_proper = 0;
for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  const m4Xtrain = [...m4X.slice(0, i), ...m4X.slice(i + 1)];
  const vflatTrain = [...logVflat.slice(0, i), ...logVflat.slice(i + 1)];

  const vflatFitTrain = ols(vflatTrain, m4Xtrain);
  const vflatPerpTrain = vflatFitTrain.resid;

  const vflatTest_pred = [1, ...m4X[i]].reduce((s, x, j) => s + x * vflatFitTrain.beta[j], 0);
  const vflatPerpTest = logVflat[i] - vflatTest_pred;

  const augXtrain = m4Xtrain.map((r, j) => [...r, vflatPerpTrain[j]]);
  const augFit = ols(Ytrain, augXtrain);
  const xiTest = [1, ...m4X[i], vflatPerpTest];
  const pred = xiTest.reduce((s, x, j) => s + x * augFit.beta[j], 0);
  m4PlusVperpSS_proper += (Y[i] - pred) ** 2;
}
const m4PlusVperpGap_proper = gapPct(Math.sqrt(m4PlusVperpSS_proper / N), sdY);

console.log('  M4:                           gap = ' + m4Gap.toFixed(1) + '%');
console.log('  M4 + Vflat_perp (proper LOO): gap = ' + m4PlusVperpGap_proper.toFixed(1) + '%');
console.log('  Delta:                         ' + (m4PlusVperpGap_proper - m4Gap > 0 ? '+' : '') + (m4PlusVperpGap_proper - m4Gap).toFixed(1) + 'pp');
console.log('');

const vflatPerpAlive = m4PlusVperpGap_proper > m4Gap;
console.log('  VERDICT A: Vflat_perp ' + (vflatPerpAlive ?
  'RETAINS predictive value after removing Sigma0 → independent axis (' + (m4PlusVperpGap_proper - m4Gap).toFixed(1) + 'pp gain)' :
  'LOSES predictive value after removing Sigma0 → proxy'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST B: SIGMA0_PERP — Orthogonalized Sigma0 beyond M3+Vflat');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  Global residualization (diagnostic only):');
const sig0Perp = ols(logSig0, mCX).resid;
const rSig0PerpA0 = pearsonR(sig0Perp, Y);
const rSig0PerpCresid = pearsonR(sig0Perp, mCFit.resid);
console.log('    r(Sigma0_perp, logA0) = ' + rSig0PerpA0.toFixed(3));
console.log('    r(Sigma0_perp, ModelC_residuals) = ' + rSig0PerpCresid.toFixed(3));
console.log('');

console.log('  Proper LOO with fold-internal residualization:\n');

let mCPlusSigPerpSS_proper = 0;
for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  const mCXtrain = [...mCX.slice(0, i), ...mCX.slice(i + 1)];
  const sig0Train = [...logSig0.slice(0, i), ...logSig0.slice(i + 1)];

  const sig0FitTrain = ols(sig0Train, mCXtrain);
  const sig0PerpTrain = sig0FitTrain.resid;

  const sig0Test_pred = [1, ...mCX[i]].reduce((s, x, j) => s + x * sig0FitTrain.beta[j], 0);
  const sig0PerpTest = logSig0[i] - sig0Test_pred;

  const augXtrain = mCXtrain.map((r, j) => [...r, sig0PerpTrain[j]]);
  const augFit = ols(Ytrain, augXtrain);
  const xiTest = [1, ...mCX[i], sig0PerpTest];
  const pred = xiTest.reduce((s, x, j) => s + x * augFit.beta[j], 0);
  mCPlusSigPerpSS_proper += (Y[i] - pred) ** 2;
}
const mCPlusSigPerpGap_proper = gapPct(Math.sqrt(mCPlusSigPerpSS_proper / N), sdY);

console.log('  Model C:                          gap = ' + mCGap.toFixed(1) + '%');
console.log('  Model C + Sig0_perp (proper LOO): gap = ' + mCPlusSigPerpGap_proper.toFixed(1) + '%');
console.log('  Delta:                             ' + (mCPlusSigPerpGap_proper - mCGap > 0 ? '+' : '') + (mCPlusSigPerpGap_proper - mCGap).toFixed(1) + 'pp');
console.log('');

const sig0PerpAlive = mCPlusSigPerpGap_proper > mCGap;
console.log('  VERDICT B: Sigma0_perp ' + (sig0PerpAlive ?
  'SURVIVES after removing Vflat → carries independent info' :
  'DIES after removing Vflat → Vflat subsumes Sigma0 (' + (mCPlusSigPerpGap_proper - mCGap).toFixed(1) + 'pp)'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST C: EXPLANATORY CLEANLINESS — VIF + STABILITY');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const m4VarNames = ['logMHI', 'logMhost', 'logMR', 'logSig0'];
const mCVarNames = ['logMHI', 'logMhost', 'logMR', 'logVflat'];

console.log('  VIF (Variance Inflation Factors):\n');
console.log('  M4 (MHI+Mhost+MR+Sigma0):');
for (let j = 0; j < 4; j++) {
  const v = vif(m4X, j);
  const flag = v > 5 ? ' ★★★ PROBLEMATIC' : v > 2.5 ? ' ★★ ELEVATED' : '';
  console.log('    VIF(' + m4VarNames[j].padEnd(10) + ') = ' + v.toFixed(2) + flag);
}
console.log('');

console.log('  Model C (MHI+Mhost+MR+Vflat):');
for (let j = 0; j < 4; j++) {
  const v = vif(mCX, j);
  const flag = v > 5 ? ' ★★★ PROBLEMATIC' : v > 2.5 ? ' ★★ ELEVATED' : '';
  console.log('    VIF(' + mCVarNames[j].padEnd(10) + ') = ' + v.toFixed(2) + flag);
}
console.log('');

const NBOOT = 2000;
function bootCoeffSpread(Y, X, nBoot, seed) {
  const rng = mulberry32(seed);
  const n = Y.length, p = X[0].length;
  const baseFit = ols(Y, X);
  const baseSigns = baseFit.beta.slice(1).map(b => Math.sign(b));
  const flips = new Array(p).fill(0);
  const coeffs = Array.from({ length: p + 1 }, () => []);
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(rng() * n));
    const fit = ols(idx.map(i => Y[i]), idx.map(i => X[i]));
    for (let j = 0; j <= p; j++) coeffs[j].push(fit.beta[j]);
    for (let j = 0; j < p; j++) if (Math.sign(fit.beta[j + 1]) !== baseSigns[j]) flips[j]++;
  }
  return {
    sds: coeffs.map(c => sd(c)),
    cvs: coeffs.map((c, j) => { const m = mean(c); return m !== 0 ? sd(c) / Math.abs(m) : Infinity; }),
    flipRates: flips.map(f => f / nBoot)
  };
}

const m4Boot = bootCoeffSpread(Y, m4X, NBOOT, 42);
const mCBoot = bootCoeffSpread(Y, mCX, NBOOT, 42);

console.log('  Bootstrap coefficient variation (CV = SD/|mean|):\n');
console.log('    ' + 'Variable'.padEnd(12) + '  M4 CV     M4 flip%   |  ModelC CV   C flip%');
for (let j = 0; j < 3; j++) {
  const m4CV = m4Boot.cvs[j + 1];
  const mcCV = mCBoot.cvs[j + 1];
  const m4Flip = m4Boot.flipRates[j] * 100;
  const mcFlip = mCBoot.flipRates[j] * 100;
  console.log('    ' + m4VarNames[j].padEnd(12) + '  ' + m4CV.toFixed(3).padStart(6) + '   ' + m4Flip.toFixed(1).padStart(6) + '%  |  ' + mcCV.toFixed(3).padStart(6) + '     ' + mcFlip.toFixed(1).padStart(6) + '%');
}
console.log('    ' + 'logSig0'.padEnd(12) + '  ' + m4Boot.cvs[4].toFixed(3).padStart(6) + '   ' + (m4Boot.flipRates[3] * 100).toFixed(1).padStart(6) + '%  |');
console.log('    ' + 'logVflat'.padEnd(12) + '                        |  ' + mCBoot.cvs[4].toFixed(3).padStart(6) + '     ' + (mCBoot.flipRates[3] * 100).toFixed(1).padStart(6) + '%');
console.log('');

const m4SumCV = m4Boot.cvs.slice(1).reduce((s, v) => s + v, 0);
const mcSumCV = mCBoot.cvs.slice(1).reduce((s, v) => s + v, 0);
const m4MaxFlip = Math.max(...m4Boot.flipRates);
const mcMaxFlip = Math.max(...mCBoot.flipRates);

console.log('  Summary:');
console.log('    M4 total CV: ' + m4SumCV.toFixed(3) + ', max flip: ' + (m4MaxFlip * 100).toFixed(1) + '%');
console.log('    MC total CV: ' + mcSumCV.toFixed(3) + ', max flip: ' + (mcMaxFlip * 100).toFixed(1) + '%');
const m4Cleaner = m4SumCV < mcSumCV;
console.log('    Cleaner model: ' + (m4Cleaner ? 'M4' : 'Model C'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST D: PHYSICAL DECOMPOSITION OF VFLAT (corrected Mbar)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  logMbar now computed with correct gas mass units:\n');
gals45.slice(0, 3).forEach(g => {
  const s = sparcMap[g.name];
  const mStar = s.L36 * 0.5 * 1e9;
  const mGas = Math.pow(10, g.logMHI) * 1.33 * 1e9;
  const mbar = mStar + mGas;
  console.log('    ' + g.name + ': Mstar=' + mStar.toExponential(2) + ' Mgas=' + mGas.toExponential(2) + ' logMbar=' + Math.log10(mbar).toFixed(3));
});
console.log('');

console.log('  What is logVflat correlated with?\n');

const physVars = [
  { name: 'logMHI', arr: logMHI, desc: 'gas mass' },
  { name: 'logMhost', arr: logMhost, desc: 'host halo mass' },
  { name: 'logMR', arr: logMR, desc: 'dynamical coherence' },
  { name: 'logSigma0', arr: logSig0, desc: 'surface density' },
  { name: 'logL36', arr: logL36, desc: 'stellar luminosity' },
  { name: 'logMbar', arr: logMbar, desc: 'total baryon mass (CORRECTED)' },
  { name: 'logRdisk', arr: logRdisk, desc: 'disk scale length' },
  { name: 'morphT', arr: morphT, desc: 'morphological type' }
];

for (const pv of physVars) {
  const r = pearsonR(logVflat, pv.arr);
  const flag = Math.abs(r) > 0.7 ? ' ★★★' : Math.abs(r) > 0.5 ? ' ★★' : Math.abs(r) > 0.3 ? ' ★' : '';
  console.log('    r(Vflat, ' + pv.name.padEnd(12) + ') = ' + r.toFixed(3) + flag + '  (' + pv.desc + ')');
}
console.log('');

console.log('  What does Vflat_perp (global) correlate with?\n');
for (const pv of physVars) {
  const r = pearsonR(vflatPerp, pv.arr);
  const flag = Math.abs(r) > 0.3 ? ' ★' : '';
  console.log('    r(Vflat_perp, ' + pv.name.padEnd(12) + ') = ' + r.toFixed(3) + flag);
}
console.log('');

console.log('  Decomposition: What fraction of Vflat is each component?\n');

const decompTargets = [
  { name: 'logMHI only', X: gals45.map((_, i) => [logMHI[i]]) },
  { name: 'logMhost only', X: gals45.map((_, i) => [logMhost[i]]) },
  { name: 'logSigma0 only', X: gals45.map((_, i) => [logSig0[i]]) },
  { name: 'logMR only', X: gals45.map((_, i) => [logMR[i]]) },
  { name: 'logL36 only', X: gals45.map((_, i) => [logL36[i]]) },
  { name: 'logMbar only', X: gals45.map((_, i) => [logMbar[i]]) },
  { name: 'MHI+Mhost', X: gals45.map((_, i) => [logMHI[i], logMhost[i]]) },
  { name: 'MHI+Mhost+Sig0', X: gals45.map((_, i) => [logMHI[i], logMhost[i], logSig0[i]]) },
  { name: 'MHI+Mhost+MR+Sig0', X: m4X },
  { name: 'L36+Rdisk', X: gals45.map((_, i) => [logL36[i], logRdisk[i]]) },
  { name: 'Mbar+Rdisk', X: gals45.map((_, i) => [logMbar[i], logRdisk[i]]) }
];

for (const dt of decompTargets) {
  const fit = ols(logVflat, dt.X);
  console.log('    R2(Vflat ~ ' + dt.name.padEnd(22) + ') = ' + fit.r2.toFixed(3) + ' (' + (fit.r2 * 100).toFixed(1) + '%)');
}
console.log('');

const fitVflatByMbar = ols(logVflat, gals45.map((_, i) => [logMbar[i]]));
const fitVflatByMHI = ols(logVflat, gals45.map((_, i) => [logMHI[i]]));
const fitVflatBySig0 = ols(logVflat, gals45.map((_, i) => [logSig0[i]]));
const fitVflatByL36 = ols(logVflat, gals45.map((_, i) => [logL36[i]]));
const fitVflatByM4 = ols(logVflat, m4X);

console.log('  Summary:');
console.log('    Vflat = mass-scale indicator: R2(Mbar)=' + fitVflatByMbar.r2.toFixed(3) +
  ', R2(MHI)=' + fitVflatByMHI.r2.toFixed(3) + ', R2(L36)=' + fitVflatByL36.r2.toFixed(3));
console.log('    M4 predictors explain ' + (fitVflatByM4.r2 * 100).toFixed(1) + '% of Vflat');
console.log('    Unique kinematic fraction: ' + ((1 - fitVflatByM4.r2) * 100).toFixed(1) + '%');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST E: ASYMMETRY TEST — Who subsumes whom?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  Diagnostic correlations (global):');
console.log('    r(Vflat_perp, M4_resid) = ' + rVflatPerpM4resid.toFixed(3));
console.log('    r(Sig0_perp, MC_resid)  = ' + rSig0PerpCresid.toFixed(3));
console.log('');

console.log('  Functional test (proper LOO gap impact):');
console.log('    Vflat_perp added to M4:  ' + (m4PlusVperpGap_proper - m4Gap > 0 ? '+' : '') + (m4PlusVperpGap_proper - m4Gap).toFixed(1) + 'pp ' + (vflatPerpAlive ? '→ ALIVE' : '→ DEAD'));
console.log('    Sig0_perp added to MC:   ' + (mCPlusSigPerpGap_proper - mCGap > 0 ? '+' : '') + (mCPlusSigPerpGap_proper - mCGap).toFixed(1) + 'pp ' + (sig0PerpAlive ? '→ ALIVE' : '→ DEAD'));
console.log('');

let asymVerdict;
if (vflatPerpAlive && !sig0PerpAlive) {
  asymVerdict = 'VFLAT_SUBSUMES_SIGMA0';
  console.log('  Pattern: Vflat alive after Sigma0, Sigma0 dead after Vflat');
  console.log('  → Vflat strictly subsumes Sigma0 — it contains Sigma0 + extra');
} else if (vflatPerpAlive && sig0PerpAlive) {
  asymVerdict = 'BOTH_INDEPENDENT';
  console.log('  Pattern: Both survive after residualization');
  console.log('  → Two partially overlapping but distinct 4th axes');
} else if (!vflatPerpAlive && sig0PerpAlive) {
  asymVerdict = 'SIGMA0_SUBSUMES_VFLAT';
  console.log('  Pattern: Sigma0 alive after Vflat, Vflat dead after Sigma0');
  console.log('  → Sigma0 is the deeper axis');
} else {
  asymVerdict = 'MUTUAL_PROXY';
  console.log('  Pattern: Both dead after residualization');
  console.log('  → They are the same axis expressed differently');
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST F: LOO GAP DECOMPOSITION — WHERE does Vflat win?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const m4Errors = [];
const mCErrors = [];
for (let i = 0; i < N; i++) {
  const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];

  const Xt4 = [...m4X.slice(0, i), ...m4X.slice(i + 1)];
  const f4 = ols(Yt, Xt4);
  const xi4 = [1, ...m4X[i]];
  m4Errors.push(Y[i] - xi4.reduce((s, x, j) => s + x * f4.beta[j], 0));

  const XtC = [...mCX.slice(0, i), ...mCX.slice(i + 1)];
  const fC = ols(Yt, XtC);
  const xiC = [1, ...mCX[i]];
  mCErrors.push(Y[i] - xiC.reduce((s, x, j) => s + x * fC.beta[j], 0));
}

const diffs = m4Errors.map((e, i) => Math.abs(e) - Math.abs(mCErrors[i]));
const nCwins = diffs.filter(d => d > 0.005).length;
const nM4wins = diffs.filter(d => d < -0.005).length;
const nTied = N - nCwins - nM4wins;

console.log('  Galaxy-by-galaxy LOO error comparison (M4 vs C):');
console.log('    Model C more accurate: ' + nCwins + '/' + N);
console.log('    M4 more accurate:      ' + nM4wins + '/' + N);
console.log('    Effectively tied:      ' + nTied + '/' + N);
console.log('');

const ranked = gals45.map((g, i) => ({
  name: g.name,
  m4err: Math.abs(m4Errors[i]),
  mCerr: Math.abs(mCErrors[i]),
  diff: Math.abs(m4Errors[i]) - Math.abs(mCErrors[i]),
  vflat: Math.pow(10, logVflat[i]),
  sig0: Math.pow(10, logSig0[i])
})).sort((a, b) => b.diff - a.diff);

console.log('  Top 5 galaxies where Model C wins most:');
for (let i = 0; i < 5; i++) {
  const g = ranked[i];
  console.log('    ' + g.name.padEnd(18) + ' M4err=' + g.m4err.toFixed(3) + ' Cerr=' + g.mCerr.toFixed(3) +
    ' diff=' + g.diff.toFixed(3) + ' Vflat=' + g.vflat.toFixed(0) + ' Sig0=' + g.sig0.toFixed(1));
}
console.log('');

console.log('  Top 5 galaxies where M4 wins most:');
const rankedRev = [...ranked].reverse();
for (let i = 0; i < 5; i++) {
  const g = rankedRev[i];
  console.log('    ' + g.name.padEnd(18) + ' M4err=' + g.m4err.toFixed(3) + ' Cerr=' + g.mCerr.toFixed(3) +
    ' diff=' + g.diff.toFixed(3) + ' Vflat=' + g.vflat.toFixed(0) + ' Sig0=' + g.sig0.toFixed(1));
}
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 128: FINAL SYNTHESIS (PATCHED)');
console.log('══════════════════════════════════════════════════════════════════════\n');

console.log('  EVIDENCE SUMMARY:\n');
console.log('    A) Vflat_perp beyond M4 (proper LOO):  ' +
  (m4PlusVperpGap_proper - m4Gap > 0 ? '+' : '') + (m4PlusVperpGap_proper - m4Gap).toFixed(1) + 'pp' +
  (vflatPerpAlive ? ' — ALIVE' : ' — DEAD'));
console.log('    B) Sigma0_perp beyond C (proper LOO):   ' +
  (mCPlusSigPerpGap_proper - mCGap > 0 ? '+' : '') + (mCPlusSigPerpGap_proper - mCGap).toFixed(1) + 'pp' +
  (sig0PerpAlive ? ' — ALIVE' : ' — DEAD'));
console.log('    C) Cleaner model:                       ' + (m4Cleaner ? 'M4' : 'Model C') + ' (CV sum: M4=' + m4SumCV.toFixed(3) + ', MC=' + mcSumCV.toFixed(3) + ')');
console.log('    D) Vflat nature:                        mass integrator (R2_Mbar=' + fitVflatByMbar.r2.toFixed(3) + ', R2_L36=' + fitVflatByL36.r2.toFixed(3) + ')');
console.log('    E) Asymmetry pattern:                   ' + asymVerdict);
console.log('    F) Galaxy-level advantage:               C wins ' + nCwins + '/' + N + ', M4 wins ' + nM4wins + '/' + N);
console.log('');

let finalConclusion;
if (asymVerdict === 'VFLAT_SUBSUMES_SIGMA0') {
  finalConclusion = 'Vflat is a deeper 4th axis that strictly contains Sigma0. Model C (MHI+Mhost+MR+Vflat) is the primary law; M4 (with Sigma0) is a structurally cleaner but informationally weaker projection.';
} else if (asymVerdict === 'BOTH_INDEPENDENT') {
  finalConclusion = 'Vflat and Sigma0 carry partially distinct information. Two co-equal 4-variable laws exist with different physical emphases: C for prediction, M4 for interpretation.';
} else if (asymVerdict === 'MUTUAL_PROXY') {
  finalConclusion = 'Vflat and Sigma0 are essentially the same axis. Choice is a matter of interpretive preference.';
} else {
  finalConclusion = 'Sigma0 is the deeper axis; M4 remains primary.';
}

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
const lines = finalConclusion.match(/.{1,62}/g) || [finalConclusion];
for (const line of lines) {
  console.log('  ║  ' + line.padEnd(62) + '║');
}
console.log('  ╚══════════════════════════════════════════════════════════════════╝');
console.log('');

console.log('  FRAMEWORK STATUS:\n');
if (asymVerdict === 'VFLAT_SUBSUMES_SIGMA0') {
  console.log('    PRIMARY (predictive):    Model C = MHI + Mhost + MR + Vflat  (gap=' + mCGap.toFixed(1) + '%)');
  console.log('    SECONDARY (structural):  M4     = MHI + Mhost + MR + Sigma0 (gap=' + m4Gap.toFixed(1) + '%)');
  console.log('    INTERPRETATION:  Sigma0 is a structural projection of what Vflat');
  console.log('                     captures more fully as a kinematic integrator');
} else if (asymVerdict === 'BOTH_INDEPENDENT') {
  console.log('    PREDICTIVE LAW:   Model C = MHI + Mhost + MR + Vflat  (gap=' + mCGap.toFixed(1) + '%)');
  console.log('    INTERPRETIVE LAW: M4     = MHI + Mhost + MR + Sigma0 (gap=' + m4Gap.toFixed(1) + '%)');
  console.log('    RELATIONSHIP:    Two views of a 4D state law');
} else {
  console.log('    PRIMARY:   M4 = MHI + Mhost + MR + Sigma0 (gap=' + m4Gap.toFixed(1) + '%)');
  console.log('    ALT:       Model C = MHI + Mhost + MR + Vflat  (gap=' + mCGap.toFixed(1) + '%)');
}
console.log('    COMMON CORE:     MHI + Mhost + MR (3 axes firmly established)');
console.log('');

const output = {
  phase: '128',
  title: 'Decode the 4th Axis — Vflat vs Sigma0 (PATCHED)',
  patches: ['Bug1: logMbar gas mass ×1e9 corrected', 'Bug2: LOO residualization done inside each fold', 'Verdict: functional criterion (gap impact) replaces raw correlation threshold'],
  test_A: {
    vflat_perp_r_with_logA0: +rVflatPerpA0.toFixed(3),
    vflat_perp_r_with_m4_resid: +rVflatPerpM4resid.toFixed(3),
    m4_gap: +m4Gap.toFixed(1),
    m4_plus_vflat_perp_gap_proper: +m4PlusVperpGap_proper.toFixed(1),
    delta_pp: +(m4PlusVperpGap_proper - m4Gap).toFixed(1),
    vflat_perp_alive: vflatPerpAlive
  },
  test_B: {
    sig0_perp_r_with_logA0: +rSig0PerpA0.toFixed(3),
    sig0_perp_r_with_mc_resid: +rSig0PerpCresid.toFixed(3),
    mc_gap: +mCGap.toFixed(1),
    mc_plus_sig0_perp_gap_proper: +mCPlusSigPerpGap_proper.toFixed(1),
    delta_pp: +(mCPlusSigPerpGap_proper - mCGap).toFixed(1),
    sig0_perp_alive: sig0PerpAlive
  },
  test_C: {
    m4_total_cv: +m4SumCV.toFixed(3),
    mc_total_cv: +mcSumCV.toFixed(3),
    m4_max_flip: +(m4MaxFlip * 100).toFixed(1),
    mc_max_flip: +(mcMaxFlip * 100).toFixed(1),
    cleaner: m4Cleaner ? 'M4' : 'ModelC'
  },
  test_D: {
    vflat_r2_by_mbar_corrected: +fitVflatByMbar.r2.toFixed(3),
    vflat_r2_by_mhi: +fitVflatByMHI.r2.toFixed(3),
    vflat_r2_by_sigma0: +fitVflatBySig0.r2.toFixed(3),
    vflat_r2_by_l36: +fitVflatByL36.r2.toFixed(3),
    vflat_r2_by_m4: +fitVflatByM4.r2.toFixed(3),
    vflat_unique_pct: +((1 - fitVflatByM4.r2) * 100).toFixed(1)
  },
  test_E: {
    asymmetry_verdict: asymVerdict,
    vflat_perp_signal: +rVflatPerpM4resid.toFixed(3),
    sig0_perp_signal: +rSig0PerpCresid.toFixed(3),
    vflat_perp_gap_delta: +(m4PlusVperpGap_proper - m4Gap).toFixed(1),
    sig0_perp_gap_delta: +(mCPlusSigPerpGap_proper - mCGap).toFixed(1)
  },
  test_F: {
    c_wins: nCwins,
    m4_wins: nM4wins,
    tied: nTied
  },
  conclusion: finalConclusion,
  models: {
    m4: { gap: +m4Gap.toFixed(1), formula: 'MHI + Mhost + MR + Sigma0' },
    modelC: { gap: +mCGap.toFixed(1), formula: 'MHI + Mhost + MR + Vflat' }
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase128-decode-4th-axis.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase128-decode-4th-axis.json');

#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 129: PHYSICAL DECOMPOSITION OF VFLAT');
console.log('');
console.log('  Central question: Is Vflat a fundamental 4th axis,');
console.log('  or a superior kinematic proxy for integrated baryonic structure?');
console.log('');
console.log('  Tests:');
console.log('    1) Substitution — can structure bundles replace Vflat?');
console.log('    2) Residual asymmetry — what survives after cross-removal?');
console.log('    3) Compression — PC1 vs Vflat at equal parameter count');
console.log('    4) Interpretation map — what IS the surviving Vflat fraction?');
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
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(rng() * n));
    const fit = ols(idx.map(i => Y[i]), idx.map(i => X[i]));
    for (let j = 0; j < p; j++) if (Math.sign(fit.beta[j + 1]) !== baseSigns[j]) flips[j]++;
  }
  return { flipRates: flips.map(f => f / nBoot), maxFlip: Math.max(...flips.map(f => f / nBoot)) };
}

function computePC1(vars) {
  const n = vars[0].length;
  const p = vars.length;
  const means = vars.map(v => mean(v));
  const sds2 = vars.map(v => sd(v));
  const Z = Array.from({ length: n }, (_, i) => vars.map((v, j) => (v[i] - means[j]) / sds2[j]));

  const cov = Array.from({ length: p }, () => new Array(p).fill(0));
  for (let i = 0; i < p; i++)
    for (let j = 0; j < p; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += Z[k][i] * Z[k][j];
      cov[i][j] = s / (n - 1);
    }

  let vec = new Array(p).fill(1 / Math.sqrt(p));
  for (let iter = 0; iter < 200; iter++) {
    const newVec = new Array(p).fill(0);
    for (let i = 0; i < p; i++)
      for (let j = 0; j < p; j++) newVec[i] += cov[i][j] * vec[j];
    const norm = Math.sqrt(newVec.reduce((s, v) => s + v * v, 0));
    for (let i = 0; i < p; i++) newVec[i] /= norm;
    vec = newVec;
  }

  let eigenvalue = 0;
  const Av = new Array(p).fill(0);
  for (let i = 0; i < p; i++)
    for (let j = 0; j < p; j++) Av[i] += cov[i][j] * vec[j];
  eigenvalue = vec.reduce((s, v, i) => s + v * Av[i], 0);

  const scores = Z.map(row => row.reduce((s, z, j) => s + z * vec[j], 0));
  return { scores, loadings: vec, eigenvalue, varExplained: eigenvalue / p, means, sds: sds2 };
}

function computePC1fold(vars, trainIdx) {
  const p = vars.length;
  const trainVals = vars.map(v => trainIdx.map(i => v[i]));
  const means = trainVals.map(v => mean(v));
  const sds2 = trainVals.map(v => sd(v));
  const nTr = trainIdx.length;
  const Z = Array.from({ length: nTr }, (_, i) => vars.map((v, j) => (trainVals[j][i] - means[j]) / sds2[j]));

  const cov = Array.from({ length: p }, () => new Array(p).fill(0));
  for (let i = 0; i < p; i++)
    for (let j = 0; j < p; j++) {
      let s = 0;
      for (let k = 0; k < nTr; k++) s += Z[k][i] * Z[k][j];
      cov[i][j] = s / (nTr - 1);
    }

  let vec = new Array(p).fill(1 / Math.sqrt(p));
  for (let iter = 0; iter < 200; iter++) {
    const newVec = new Array(p).fill(0);
    for (let i = 0; i < p; i++)
      for (let j = 0; j < p; j++) newVec[i] += cov[i][j] * vec[j];
    const norm = Math.sqrt(newVec.reduce((s, v) => s + v * v, 0));
    for (let i = 0; i < p; i++) newVec[i] /= norm;
    vec = newVec;
  }

  return { loadings: vec, means, sds: sds2 };
}

function applyPC1(vars, idx, loadings, means, sds2) {
  return vars.map((v, j) => (v[idx] - means[j]) / sds2[j]).reduce((s, z, j) => s + z * loadings[j], 0);
}

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const logMHI = gals45.map(g => g.logMHI);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
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

const structVars = [logMbar, logL36, logRdisk, morphT];
const structNames = ['logMbar', 'logL36', 'logRdisk', 'morphT'];
const pc1Global = computePC1(structVars);
const pc1Scores = pc1Global.scores;

console.log('  N = ' + N + ', SD(logA0) = ' + sdY.toFixed(4) + ' dex\n');
console.log('  PC1_structure loadings (global):');
structNames.forEach((n, j) => console.log('    ' + n.padEnd(10) + ': ' + pc1Global.loadings[j].toFixed(3)));
console.log('    Variance explained: ' + (pc1Global.varExplained * 100).toFixed(1) + '%\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: SUBSTITUTION — Can structure bundles replace Vflat?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const NBOOT = 2000;
const core = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]);

const models = {
  C:  { name: 'MHI+Mhost+MR+Vflat',              X: gals45.map((_, i) => [...core[i], logVflat[i]]) },
  S1: { name: 'MHI+Mhost+MR+Mbar',               X: gals45.map((_, i) => [...core[i], logMbar[i]]) },
  S2: { name: 'MHI+Mhost+MR+L36+Rdisk',           X: gals45.map((_, i) => [...core[i], logL36[i], logRdisk[i]]) },
  S3: { name: 'MHI+Mhost+MR+Mbar+Rdisk',          X: gals45.map((_, i) => [...core[i], logMbar[i], logRdisk[i]]) },
  S4: { name: 'MHI+Mhost+MR+L36+Rdisk+morphT',    X: gals45.map((_, i) => [...core[i], logL36[i], logRdisk[i], morphT[i]]) },
  S5: { name: 'MHI+Mhost+MR+PC1_structure',        X: gals45.map((_, i) => [...core[i], pc1Scores[i]]) }
};

const results = {};
console.log('  ' + 'Model'.padEnd(6) + ' ' + 'Description'.padEnd(34) + '  k   LOO gap%  AIC       BIC       maxFlip%');
console.log('  ' + '─'.repeat(95));

for (const [key, mod] of Object.entries(models)) {
  const fit = ols(Y, mod.X);
  const loo = looCV(Y, mod.X);
  const gap = gapPct(loo, sdY);
  const boot = bootstrapSigns(Y, mod.X, NBOOT, 42);
  results[key] = { fit, loo, gap, boot, name: mod.name, X: mod.X };
  console.log('  ' + key.padEnd(6) + ' ' + mod.name.padEnd(34) + '  ' +
    fit.k + '   ' + gap.toFixed(1).padStart(6) + '%  ' +
    fit.aic.toFixed(2).padStart(8) + '  ' + fit.bic.toFixed(2).padStart(8) + '  ' +
    (boot.maxFlip * 100).toFixed(1).padStart(6) + '%');
}
console.log('');

const cGap = results.C.gap;
console.log('  Substitution summary (vs Model C gap=' + cGap.toFixed(1) + '%):');
for (const key of ['S1', 'S2', 'S3', 'S4', 'S5']) {
  const delta = results[key].gap - cGap;
  const status = delta >= -1 ? 'CATCHES C' : delta >= -5 ? 'CLOSE' : 'FAR BEHIND';
  console.log('    ' + key + ': gap=' + results[key].gap.toFixed(1) + '% (' + (delta > 0 ? '+' : '') + delta.toFixed(1) + 'pp) → ' + status);
}
console.log('');

console.log('  Nested CV: C vs S1 vs S2 vs S3 vs S4 vs S5');
const allKeys = ['C', 'S1', 'S2', 'S3', 'S4', 'S5'];
const nestedSel = {}; allKeys.forEach(k => nestedSel[k] = 0);
let nestedSS = 0;

for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  let bestGap = -Infinity, bestPred = 0, bestMod = '';

  for (const mk of allKeys) {
    let Xfull;
    if (mk === 'S5') {
      const trainIdx = Array.from({ length: N }, (_, j) => j).filter(j => j !== i);
      const pc1f = computePC1fold(structVars, trainIdx);
      Xfull = gals45.map((_, j) => [...core[j], applyPC1(structVars, j, pc1f.loadings, pc1f.means, pc1f.sds)]);
    } else {
      Xfull = models[mk].X;
    }

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
  nestedSel[bestMod]++;
  nestedSS += (Y[i] - bestPred) ** 2;
}
const nestedGap = gapPct(Math.sqrt(nestedSS / N), sdY);

console.log('    Nested gap: ' + nestedGap.toFixed(1) + '%');
console.log('    Selection frequency:');
for (const mk of allKeys) {
  console.log('      ' + mk.padEnd(4) + ': ' + nestedSel[mk] + '/' + N + ' (' + (nestedSel[mk] / N * 100).toFixed(0) + '%)');
}
console.log('');

const bestSub = ['S1', 'S2', 'S3', 'S4', 'S5'].reduce((best, k) => results[k].gap > results[best].gap ? k : best, 'S1');
const bestSubGap = results[bestSub].gap;
const t1verdict = cGap - bestSubGap <= 1 ? 'REPLACEABLE — structure bundle matches Vflat' :
  cGap - bestSubGap <= 3 ? 'NEARLY_REPLACEABLE — close but Vflat still leads' :
  'IRREPLACEABLE — Vflat carries unique information beyond structure';
console.log('  TEST 1 VERDICT: ' + t1verdict);
console.log('    Best challenger: ' + bestSub + ' (' + results[bestSub].name + ') gap=' + bestSubGap.toFixed(1) + '%, deficit=' + (cGap - bestSubGap).toFixed(1) + 'pp');
console.log('    C nested wins: ' + nestedSel.C + '/' + N);
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: RESIDUAL ASYMMETRY — fold-internal orthogonalization');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const structBundleX = gals45.map((_, i) => [logMbar[i], logL36[i], logRdisk[i], morphT[i]]);
const bestBundleKey = bestSub;
const bestBundleX = models[bestBundleKey].X;

console.log('  A) Vflat_perp_structure: Does Vflat add after removing structure bundle?\n');
console.log('     Residualizing Vflat against [logMbar, logL36, logRdisk, morphT] inside each fold\n');

let ssBase_struct = 0;
let ssAug_vfperp = 0;

for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  const coreTrain = [...core.slice(0, i), ...core.slice(i + 1)];
  const structTrain = [...structBundleX.slice(0, i), ...structBundleX.slice(i + 1)];
  const vflatTrain = [...logVflat.slice(0, i), ...logVflat.slice(i + 1)];

  const baseX = coreTrain.map((c, j) => [...c, ...structTrain[j]]);
  const baseFit = ols(Ytrain, baseX);
  const xiBase = [1, ...core[i], ...structBundleX[i]];
  ssBase_struct += (Y[i] - xiBase.reduce((s, x, j) => s + x * baseFit.beta[j], 0)) ** 2;

  const vfResidFit = ols(vflatTrain, structTrain);
  const vfPerpTrain = vfResidFit.resid;
  const vfTestPred = [1, ...structBundleX[i]].reduce((s, x, j) => s + x * vfResidFit.beta[j], 0);
  const vfPerpTest = logVflat[i] - vfTestPred;

  const augX = coreTrain.map((c, j) => [...c, ...structTrain[j], vfPerpTrain[j]]);
  const augFit = ols(Ytrain, augX);
  const xiAug = [1, ...core[i], ...structBundleX[i], vfPerpTest];
  ssAug_vfperp += (Y[i] - xiAug.reduce((s, x, j) => s + x * augFit.beta[j], 0)) ** 2;
}

const gapBase_struct = gapPct(Math.sqrt(ssBase_struct / N), sdY);
const gapAug_vfperp = gapPct(Math.sqrt(ssAug_vfperp / N), sdY);
const deltaVfPerp = gapAug_vfperp - gapBase_struct;
const vfPerpAlive = deltaVfPerp > 0;

console.log('     Core+Structure bundle:              gap = ' + gapBase_struct.toFixed(1) + '%');
console.log('     Core+Structure+Vflat_perp_struct:   gap = ' + gapAug_vfperp.toFixed(1) + '%');
console.log('     Delta:                              ' + (deltaVfPerp > 0 ? '+' : '') + deltaVfPerp.toFixed(1) + 'pp');
console.log('     Vflat_perp_structure: ' + (vfPerpAlive ? 'ALIVE (' + deltaVfPerp.toFixed(1) + 'pp gain)' : 'DEAD'));
console.log('');

console.log('  B) Structure_perp_Vflat: Does structure add after removing Vflat?\n');
console.log('     Residualizing each structure var against [MHI, Mhost, MR, Vflat] inside each fold\n');

let ssBase_C = 0;
let ssAug_structPerp = 0;

for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  const cX = models.C.X;
  const cXtrain = [...cX.slice(0, i), ...cX.slice(i + 1)];

  const fitC = ols(Ytrain, cXtrain);
  const xiC = [1, ...cX[i]];
  ssBase_C += (Y[i] - xiC.reduce((s, x, j) => s + x * fitC.beta[j], 0)) ** 2;

  const structPerpTrain = [];
  const structPerpTest = [];
  for (let sv = 0; sv < 4; sv++) {
    const svArr = structVars[sv];
    const svTrain = [...svArr.slice(0, i), ...svArr.slice(i + 1)];
    const svFit = ols(svTrain, cXtrain);
    structPerpTrain.push(svFit.resid);
    const svTestPred = [1, ...cX[i]].reduce((s, x, j) => s + x * svFit.beta[j], 0);
    structPerpTest.push(svArr[i] - svTestPred);
  }

  const augX = cXtrain.map((c, j) => [...c, ...structPerpTrain.map(sp => sp[j])]);
  const augFit = ols(Ytrain, augX);
  const xiAug = [1, ...cX[i], ...structPerpTest];
  ssAug_structPerp += (Y[i] - xiAug.reduce((s, x, j) => s + x * augFit.beta[j], 0)) ** 2;
}

const gapBase_C = gapPct(Math.sqrt(ssBase_C / N), sdY);
const gapAug_structPerp = gapPct(Math.sqrt(ssAug_structPerp / N), sdY);
const deltaStructPerp = gapAug_structPerp - gapBase_C;
const structPerpAlive = deltaStructPerp > 0;

console.log('     Model C:                          gap = ' + gapBase_C.toFixed(1) + '%');
console.log('     C + Structure_perp (4 vars):       gap = ' + gapAug_structPerp.toFixed(1) + '%');
console.log('     Delta:                             ' + (deltaStructPerp > 0 ? '+' : '') + deltaStructPerp.toFixed(1) + 'pp');
console.log('     Structure_perp_Vflat: ' + (structPerpAlive ? 'ALIVE (' + deltaStructPerp.toFixed(1) + 'pp gain)' : 'DEAD'));
console.log('');

let t2verdict;
if (vfPerpAlive && !structPerpAlive) {
  t2verdict = 'VFLAT_DEEPER — Vflat contains structure + more';
} else if (vfPerpAlive && structPerpAlive) {
  t2verdict = 'PARTIAL_OVERLAP — both carry some unique info';
} else if (!vfPerpAlive && structPerpAlive) {
  t2verdict = 'STRUCTURE_DEEPER — structure bundle is the real axis';
} else {
  t2verdict = 'MUTUAL_PROXY — both express the same underlying axis';
}
console.log('  TEST 2 VERDICT: ' + t2verdict);
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: COMPRESSION — PC1 vs Vflat at equal parameter count');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  Both models have k=4 predictors (+ intercept):\n');
console.log('    Model C:  MHI + Mhost + MR + Vflat         gap = ' + results.C.gap.toFixed(1) + '%');
console.log('    Model S5: MHI + Mhost + MR + PC1_structure gap = ' + results.S5.gap.toFixed(1) + '%');
console.log('');

const comprDelta = results.C.gap - results.S5.gap;
console.log('    Gap difference: ' + comprDelta.toFixed(1) + 'pp in favor of ' + (comprDelta > 0 ? 'Vflat' : 'PC1'));
console.log('    AIC: C=' + results.C.fit.aic.toFixed(2) + ' vs S5=' + results.S5.fit.aic.toFixed(2));
console.log('    BIC: C=' + results.C.fit.bic.toFixed(2) + ' vs S5=' + results.S5.fit.bic.toFixed(2));
console.log('    Max flip: C=' + (results.C.boot.maxFlip * 100).toFixed(1) + '% vs S5=' + (results.S5.boot.maxFlip * 100).toFixed(1) + '%');
console.log('');

console.log('  Nested CV: C vs S5 (head-to-head):');
const cs5Sel = { C: 0, S5: 0 };
let cs5SS = 0;
for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  let bestGap = -Infinity, bestPred = 0, bestMod = '';

  for (const mk of ['C', 'S5']) {
    let Xfull;
    if (mk === 'S5') {
      const trainIdx = Array.from({ length: N }, (_, j) => j).filter(j => j !== i);
      const pc1f = computePC1fold(structVars, trainIdx);
      Xfull = gals45.map((_, j) => [...core[j], applyPC1(structVars, j, pc1f.loadings, pc1f.means, pc1f.sds)]);
    } else {
      Xfull = models[mk].X;
    }
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
  cs5Sel[bestMod]++;
  cs5SS += (Y[i] - bestPred) ** 2;
}
console.log('    C selected:  ' + cs5Sel.C + '/' + N);
console.log('    S5 selected: ' + cs5Sel.S5 + '/' + N);
console.log('');

let t3verdict;
if (comprDelta <= 1 && cs5Sel.S5 >= cs5Sel.C) {
  t3verdict = 'PC1_EQUIVALENT — linear compression matches Vflat';
} else if (comprDelta <= 3) {
  t3verdict = 'VFLAT_SLIGHTLY_BETTER — Vflat edges out but PC1 is close';
} else {
  t3verdict = 'VFLAT_CLEARLY_BETTER — Vflat carries more than linear compression';
}
console.log('  TEST 3 VERDICT: ' + t3verdict);
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: PHYSICAL INTERPRETATION MAP');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  A) What fraction of Vflat is explained by structure?\n');

const decomp = [
  { name: 'logMbar only', X: gals45.map((_, i) => [logMbar[i]]) },
  { name: 'logL36 only', X: gals45.map((_, i) => [logL36[i]]) },
  { name: 'logRdisk only', X: gals45.map((_, i) => [logRdisk[i]]) },
  { name: 'morphT only', X: gals45.map((_, i) => [morphT[i]]) },
  { name: 'logMbar+Rdisk', X: gals45.map((_, i) => [logMbar[i], logRdisk[i]]) },
  { name: 'logL36+Rdisk', X: gals45.map((_, i) => [logL36[i], logRdisk[i]]) },
  { name: 'Mbar+Rdisk+morphT', X: gals45.map((_, i) => [logMbar[i], logRdisk[i], morphT[i]]) },
  { name: 'L36+Rdisk+morphT', X: gals45.map((_, i) => [logL36[i], logRdisk[i], morphT[i]]) },
  { name: 'Mbar+L36+Rdisk+morphT', X: structBundleX },
  { name: 'PC1_structure', X: gals45.map((_, i) => [pc1Scores[i]]) }
];

for (const dt of decomp) {
  const fit = ols(logVflat, dt.X);
  console.log('    R2(Vflat ~ ' + dt.name.padEnd(24) + ') = ' + fit.r2.toFixed(3) + ' (' + (fit.r2 * 100).toFixed(1) + '%)');
}
const fullStructR2 = ols(logVflat, structBundleX).r2;
console.log('\n    Full structure bundle explains ' + (fullStructR2 * 100).toFixed(1) + '% of Vflat');
console.log('    Unique kinematic fraction: ' + ((1 - fullStructR2) * 100).toFixed(1) + '%');
console.log('');

console.log('  B) What does the residual Vflat (after structure) correlate with?\n');

const vfResidStruct = ols(logVflat, structBundleX).resid;

const probeVars = [
  { name: 'logMHI', arr: logMHI },
  { name: 'logMhost', arr: logMhost },
  { name: 'logMR', arr: logMR },
  { name: 'logSigma0', arr: logSig0 },
  { name: 'envCode', arr: gals45.map(g => g.envCode) },
  { name: 'rcWig', arr: gals45.map(g => g.rcWiggliness) },
  { name: 'logA0', arr: Y }
];

for (const pv of probeVars) {
  const r = pearsonR(vfResidStruct, pv.arr);
  const flag = Math.abs(r) > 0.3 ? ' ★★' : Math.abs(r) > 0.15 ? ' ★' : '';
  console.log('    r(Vflat_resid_struct, ' + pv.name.padEnd(10) + ') = ' + r.toFixed(3) + flag);
}
console.log('');

console.log('  C) Regression chain: stepwise structure explanation of Vflat\n');

const chain = [
  { name: 'logMbar', X: gals45.map((_, i) => [logMbar[i]]) },
  { name: '+logRdisk', X: gals45.map((_, i) => [logMbar[i], logRdisk[i]]) },
  { name: '+morphT', X: gals45.map((_, i) => [logMbar[i], logRdisk[i], morphT[i]]) },
  { name: '+logSig0', X: gals45.map((_, i) => [logMbar[i], logRdisk[i], morphT[i], logSig0[i]]) }
];

let prevR2 = 0;
for (const ch of chain) {
  const fit = ols(logVflat, ch.X);
  const delta = fit.r2 - prevR2;
  console.log('    ' + ch.name.padEnd(12) + ': R2 = ' + fit.r2.toFixed(3) + ' (delta = +' + (delta * 100).toFixed(1) + '%)');
  prevR2 = fit.r2;
}
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 129: FINAL SYNTHESIS');
console.log('══════════════════════════════════════════════════════════════════════\n');

console.log('  EVIDENCE SUMMARY:\n');
console.log('    Test 1 (Substitution):     ' + t1verdict.split(' — ')[0]);
console.log('      Best challenger: ' + bestSub + ' gap=' + bestSubGap.toFixed(1) + '% vs C=' + cGap.toFixed(1) + '% (deficit=' + (cGap - bestSubGap).toFixed(1) + 'pp)');
console.log('      C nested wins: ' + nestedSel.C + '/' + N);
console.log('');
console.log('    Test 2 (Asymmetry):        ' + t2verdict.split(' — ')[0]);
console.log('      Vflat_perp_structure:  ' + (vfPerpAlive ? 'ALIVE +' + deltaVfPerp.toFixed(1) + 'pp' : 'DEAD ' + deltaVfPerp.toFixed(1) + 'pp'));
console.log('      Structure_perp_Vflat:  ' + (structPerpAlive ? 'ALIVE +' + deltaStructPerp.toFixed(1) + 'pp' : 'DEAD ' + deltaStructPerp.toFixed(1) + 'pp'));
console.log('');
console.log('    Test 3 (Compression):      ' + t3verdict.split(' — ')[0]);
console.log('      C vs S5 gap diff: ' + comprDelta.toFixed(1) + 'pp, nested C=' + cs5Sel.C + ' S5=' + cs5Sel.S5);
console.log('');
console.log('    Test 4 (Interpretation):');
console.log('      Structure explains ' + (fullStructR2 * 100).toFixed(1) + '% of Vflat');
console.log('      Unique kinematic fraction: ' + ((1 - fullStructR2) * 100).toFixed(1) + '%');
console.log('');

let finalVerdict;
if (vfPerpAlive && !structPerpAlive && comprDelta > 3) {
  finalVerdict = 'A';
} else if (!vfPerpAlive && structPerpAlive) {
  finalVerdict = 'C';
} else if (!vfPerpAlive && !structPerpAlive) {
  finalVerdict = comprDelta <= 1 ? 'C' : 'B';
} else if (vfPerpAlive && comprDelta > 1) {
  finalVerdict = fullStructR2 > 0.90 ? 'B' : 'A';
} else {
  finalVerdict = 'B';
}

const verdictDescriptions = {
  A: 'Vflat is FUNDAMENTAL — it carries unique kinematic information that structure bundles cannot replicate. The 4th axis IS Vflat.',
  B: 'Vflat is a SUPERIOR COMPRESSED PROXY — most of its power comes from baryonic structure (Mbar, L36, Rdisk, morphT), but Vflat remains the best practical summary. The 4th axis is integrated baryonic structure; Vflat is the best single-variable encoding of it.',
  C: 'Vflat is REPLACEABLE — a structure bundle can fully substitute for Vflat. The 4th axis is baryonic structure itself, not any kinematic quantity.'
};

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT ' + finalVerdict + ':');
const vdLines = verdictDescriptions[finalVerdict].match(/.{1,60}/g);
for (const line of vdLines) {
  console.log('  ║  ' + line.padEnd(62) + '║');
}
console.log('  ╚══════════════════════════════════════════════════════════════════╝');
console.log('');

if (finalVerdict === 'A') {
  console.log('  FRAMEWORK CONSEQUENCE:');
  console.log('    Model C confirmed as primary: logA0 = f(MHI, Mhost, MR, Vflat)');
  console.log('    Vflat is an irreducible kinematic axis');
} else if (finalVerdict === 'B') {
  console.log('  FRAMEWORK CONSEQUENCE:');
  console.log('    Model C remains the best practical model');
  console.log('    But the 4th axis should be understood as:');
  console.log('      "integrated baryonic structure" ≈ mass + size + morphology');
  console.log('    Vflat is the kinematic readout of this structure');
  console.log('    ' + (fullStructR2 * 100).toFixed(0) + '% structural, ' + ((1 - fullStructR2) * 100).toFixed(0) + '% unique kinematic');
} else {
  console.log('  FRAMEWORK CONSEQUENCE:');
  console.log('    Model C should be replaced by a structure-based model');
  console.log('    Best candidate: ' + bestSub + ' (' + results[bestSub].name + ')');
}
console.log('');

const output = {
  phase: '129',
  title: 'Physical Decomposition of Vflat',
  N: N,
  baseline: { model: 'C', gap: +cGap.toFixed(1), formula: 'MHI+Mhost+MR+Vflat' },
  test1_substitution: {
    models: {},
    nested_cv: { gap: +nestedGap.toFixed(1), selection: nestedSel },
    best_challenger: bestSub,
    best_challenger_gap: +bestSubGap.toFixed(1),
    deficit_pp: +(cGap - bestSubGap).toFixed(1),
    verdict: t1verdict.split(' — ')[0]
  },
  test2_asymmetry: {
    vflat_perp_structure: {
      base_gap: +gapBase_struct.toFixed(1),
      augmented_gap: +gapAug_vfperp.toFixed(1),
      delta_pp: +deltaVfPerp.toFixed(1),
      alive: vfPerpAlive
    },
    structure_perp_vflat: {
      base_gap: +gapBase_C.toFixed(1),
      augmented_gap: +gapAug_structPerp.toFixed(1),
      delta_pp: +deltaStructPerp.toFixed(1),
      alive: structPerpAlive
    },
    verdict: t2verdict.split(' — ')[0]
  },
  test3_compression: {
    c_gap: +results.C.gap.toFixed(1),
    s5_gap: +results.S5.gap.toFixed(1),
    gap_diff: +comprDelta.toFixed(1),
    nested_c_vs_s5: cs5Sel,
    verdict: t3verdict.split(' — ')[0]
  },
  test4_interpretation: {
    structure_r2_of_vflat: +fullStructR2.toFixed(3),
    unique_kinematic_pct: +((1 - fullStructR2) * 100).toFixed(1),
    pc1_loadings: Object.fromEntries(structNames.map((n, j) => [n, +pc1Global.loadings[j].toFixed(3)])),
    pc1_var_explained: +(pc1Global.varExplained * 100).toFixed(1)
  },
  final_verdict: finalVerdict,
  final_verdict_description: verdictDescriptions[finalVerdict]
};

for (const [key, res] of Object.entries(results)) {
  output.test1_substitution.models[key] = {
    name: res.name,
    gap: +res.gap.toFixed(1),
    aic: +res.fit.aic.toFixed(2),
    bic: +res.fit.bic.toFixed(2),
    max_flip: +(res.boot.maxFlip * 100).toFixed(1)
  };
}

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase129-vflat-decomposition.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase129-vflat-decomposition.json');

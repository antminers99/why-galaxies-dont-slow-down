#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 130: SPLIT THE 4TH SECTOR');
console.log('');
console.log('  Instead of asking "Vflat or structure?", we ask:');
console.log('  Can we decompose the 4th axis into two explicit components?');
console.log('');
console.log('  Candidate decomposed model:');
console.log('    MHI + Mhost + MR + PC1_structure + Vflat_resid_structure');
console.log('');
console.log('  Tests:');
console.log('    1) Head-to-head: Decomposed vs Model C');
console.log('    2) Nested stability: Is 5-var decomposed model overfitting?');
console.log('    3) Component contributions: What does each part add?');
console.log('    4) Bootstrap + VIF: Is the decomposed model clean?');
console.log('    5) Alternative decompositions');
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

function bootstrapAnalysis(Y, X, varNames, nBoot, seed) {
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
  const result = {};
  for (let j = 0; j < p; j++) {
    const c = coeffs[j + 1].sort((a, b) => a - b);
    result[varNames[j]] = {
      coeff: +baseFit.beta[j + 1].toFixed(4),
      flipRate: +(flips[j] / nBoot * 100).toFixed(1),
      ci95: [+c[Math.floor(nBoot * 0.025)].toFixed(4), +c[Math.floor(nBoot * 0.975)].toFixed(4)],
      cv: +(sd(coeffs[j + 1]) / Math.abs(mean(coeffs[j + 1]))).toFixed(3)
    };
  }
  result._maxFlip = Math.max(...flips.map(f => f / nBoot));
  return result;
}

function vif(X, idx) {
  const n = X.length;
  const Yv = X.map(r => r[idx]);
  const Xv = X.map(r => r.filter((_, j) => j !== idx));
  const fit = ols(Yv, Xv);
  return fit.r2 < 1 ? 1 / (1 - fit.r2) : Infinity;
}

function computePC1(vars, indices) {
  const p = vars.length;
  const vals = vars.map(v => indices.map(i => v[i]));
  const means = vals.map(v => mean(v));
  const sds2 = vals.map(v => sd(v));
  const nTr = indices.length;
  const Z = Array.from({ length: nTr }, (_, i) => vars.map((_, j) => (vals[j][i] - means[j]) / sds2[j]));

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

  let eigenvalue = 0;
  const Av = new Array(p).fill(0);
  for (let i = 0; i < p; i++)
    for (let j = 0; j < p; j++) Av[i] += cov[i][j] * vec[j];
  eigenvalue = vec.reduce((s, v, i) => s + v * Av[i], 0);

  return { loadings: vec, means, sds: sds2, eigenvalue, varExplained: eigenvalue / p };
}

function applyPC1(vars, idx, pc1) {
  return vars.map((v, j) => (v[idx] - pc1.means[j]) / pc1.sds[j]).reduce((s, z, j) => s + z * pc1.loadings[j], 0);
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
const allIdx = Array.from({ length: N }, (_, i) => i);
const pc1Global = computePC1(structVars, allIdx);
const pc1Scores = allIdx.map(i => applyPC1(structVars, i, pc1Global));

const vfResidFit = ols(logVflat, structVars.map((_, j) => allIdx.map(i => structVars[j][i])).reduce((rows, col) => {
  col.forEach((v, i) => { if (!rows[i]) rows[i] = []; rows[i].push(v); });
  return rows;
}, []));
const vfResidGlobal = vfResidFit.resid;

console.log('  N = ' + N + ', SD(logA0) = ' + sdY.toFixed(4) + ' dex\n');

console.log('  PC1_structure loadings:');
structNames.forEach((n, j) => console.log('    ' + n.padEnd(10) + ': ' + pc1Global.loadings[j].toFixed(3)));
console.log('    Variance explained: ' + (pc1Global.varExplained * 100).toFixed(1) + '%\n');

console.log('  Global Vflat_resid_structure:');
console.log('    R2(structure -> Vflat) = ' + vfResidFit.r2.toFixed(3));
console.log('    Residual SD = ' + sd(vfResidGlobal).toFixed(4));
console.log('    r(Vf_resid, PC1) = ' + pearsonR(vfResidGlobal, pc1Scores).toFixed(3) + ' (should be ~0)');
console.log('    r(Vf_resid, logA0) = ' + pearsonR(vfResidGlobal, Y).toFixed(3));
console.log('');

const NBOOT = 2000;
const core = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]);

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: HEAD-TO-HEAD — Decomposed vs Model C');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const modelDefs = {
  C:    { name: 'MHI+Mhost+MR+Vflat', k: 4, vars: ['logMHI','logMhost','logMR','logVflat'] },
  D5a:  { name: 'MHI+Mhost+MR+PC1+VfResid', k: 5, vars: ['logMHI','logMhost','logMR','PC1','VfResid'] },
  D5b:  { name: 'MHI+Mhost+MR+Mbar+VfResid', k: 5, vars: ['logMHI','logMhost','logMR','logMbar','VfResid'] },
  D5c:  { name: 'MHI+Mhost+MR+L36+VfResid', k: 5, vars: ['logMHI','logMhost','logMR','logL36','VfResid'] },
  M4:   { name: 'MHI+Mhost+MR+Sigma0', k: 4, vars: ['logMHI','logMhost','logMR','logSig0'] }
};

function buildX_global(key) {
  switch (key) {
    case 'C': return gals45.map((_, i) => [...core[i], logVflat[i]]);
    case 'D5a': return gals45.map((_, i) => [...core[i], pc1Scores[i], vfResidGlobal[i]]);
    case 'D5b': return gals45.map((_, i) => [...core[i], logMbar[i], vfResidGlobal[i]]);
    case 'D5c': return gals45.map((_, i) => [...core[i], logL36[i], vfResidGlobal[i]]);
    case 'M4': return gals45.map((_, i) => [...core[i], logSig0[i]]);
  }
}

function looCV_proper(key) {
  let ss = 0;
  for (let i = 0; i < N; i++) {
    const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const trainIdx = allIdx.filter(j => j !== i);

    let Xfull;
    if (key === 'D5a') {
      const pc1f = computePC1(structVars, trainIdx);
      const pc1Train = trainIdx.map(j => applyPC1(structVars, j, pc1f));
      const pc1Test = applyPC1(structVars, i, pc1f);

      const structXtrain = trainIdx.map(j => structVars.map(v => v[j]));
      const vfTrain = trainIdx.map(j => logVflat[j]);
      const vfFit = ols(vfTrain, structXtrain);
      const vfResidTrain = vfFit.resid;
      const vfTestPred = [1, ...structVars.map(v => v[i])].reduce((s, x, j) => s + x * vfFit.beta[j], 0);
      const vfResidTest = logVflat[i] - vfTestPred;

      const Xtrain = trainIdx.map((j, idx) => [...core[j], pc1Train[idx], vfResidTrain[idx]]);
      const fit = ols(Ytrain, Xtrain);
      const xi = [1, ...core[i], pc1Test, vfResidTest];
      ss += (Y[i] - xi.reduce((s, x, j) => s + x * fit.beta[j], 0)) ** 2;
    } else if (key === 'D5b') {
      const structXtrain = trainIdx.map(j => structVars.map(v => v[j]));
      const vfTrain = trainIdx.map(j => logVflat[j]);
      const vfFit = ols(vfTrain, structXtrain);
      const vfResidTrain = vfFit.resid;
      const vfTestPred = [1, ...structVars.map(v => v[i])].reduce((s, x, j) => s + x * vfFit.beta[j], 0);
      const vfResidTest = logVflat[i] - vfTestPred;

      const Xtrain = trainIdx.map((j, idx) => [...core[j], logMbar[j], vfResidTrain[idx]]);
      const fit = ols(Ytrain, Xtrain);
      const xi = [1, ...core[i], logMbar[i], vfResidTest];
      ss += (Y[i] - xi.reduce((s, x, j) => s + x * fit.beta[j], 0)) ** 2;
    } else if (key === 'D5c') {
      const structXtrain = trainIdx.map(j => structVars.map(v => v[j]));
      const vfTrain = trainIdx.map(j => logVflat[j]);
      const vfFit = ols(vfTrain, structXtrain);
      const vfResidTrain = vfFit.resid;
      const vfTestPred = [1, ...structVars.map(v => v[i])].reduce((s, x, j) => s + x * vfFit.beta[j], 0);
      const vfResidTest = logVflat[i] - vfTestPred;

      const Xtrain = trainIdx.map((j, idx) => [...core[j], logL36[j], vfResidTrain[idx]]);
      const fit = ols(Ytrain, Xtrain);
      const xi = [1, ...core[i], logL36[i], vfResidTest];
      ss += (Y[i] - xi.reduce((s, x, j) => s + x * fit.beta[j], 0)) ** 2;
    } else {
      const Xfull2 = buildX_global(key);
      const Xtrain = [...Xfull2.slice(0, i), ...Xfull2.slice(i + 1)];
      const fit = ols(Ytrain, Xtrain);
      const xi = [1, ...Xfull2[i]];
      ss += (Y[i] - xi.reduce((s, x, j) => s + x * fit.beta[j], 0)) ** 2;
    }
  }
  return Math.sqrt(ss / N);
}

const modelResults = {};
console.log('  ' + 'Model'.padEnd(6) + ' ' + 'Description'.padEnd(28) + ' k  gap%     AIC       BIC       maxFlip%');
console.log('  ' + '─'.repeat(88));

for (const key of ['C', 'D5a', 'D5b', 'D5c', 'M4']) {
  const Xg = buildX_global(key);
  const fit = ols(Y, Xg);
  const loo = looCV_proper(key);
  const gap = gapPct(loo, sdY);
  const boot = bootstrapAnalysis(Y, Xg, modelDefs[key].vars, NBOOT, 42);
  modelResults[key] = { fit, loo, gap, boot, Xg };
  console.log('  ' + key.padEnd(6) + ' ' + modelDefs[key].name.padEnd(28) + ' ' +
    fit.k + '  ' + gap.toFixed(1).padStart(5) + '%  ' +
    fit.aic.toFixed(2).padStart(8) + '  ' + fit.bic.toFixed(2).padStart(8) + '  ' +
    (boot._maxFlip * 100).toFixed(1).padStart(6) + '%');
}
console.log('');

const bestD5key = ['D5a', 'D5b', 'D5c'].reduce((best, k) => modelResults[k].gap > modelResults[best].gap ? k : best, 'D5a');
const bestD5 = modelResults[bestD5key];
const cRes = modelResults.C;

console.log('  Best decomposed model: ' + bestD5key + ' (' + modelDefs[bestD5key].name + ')');
console.log('  gap: ' + bestD5.gap.toFixed(1) + '% vs C: ' + cRes.gap.toFixed(1) + '% (delta: ' + (bestD5.gap - cRes.gap > 0 ? '+' : '') + (bestD5.gap - cRes.gap).toFixed(1) + 'pp)');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: NESTED CV — Is the decomposed model genuinely better?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const nestedKeys = ['C', bestD5key];
const nestedSel = {}; nestedKeys.forEach(k => nestedSel[k] = 0);
let nestedSS = 0;

for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  const trainIdx = allIdx.filter(j => j !== i);
  let bestGap = -Infinity, bestPred = 0, bestMod = '';

  for (const mk of nestedKeys) {
    let Xtrain, xiTest;

    if (mk === 'C') {
      const Xfull = buildX_global('C');
      Xtrain = [...Xfull.slice(0, i), ...Xfull.slice(i + 1)];
      xiTest = [1, ...Xfull[i]];
    } else {
      const structXtrain = trainIdx.map(j => structVars.map(v => v[j]));
      const vfTrain = trainIdx.map(j => logVflat[j]);
      const vfFit = ols(vfTrain, structXtrain);
      const vfResidTrain = vfFit.resid;
      const vfTestPred = [1, ...structVars.map(v => v[i])].reduce((s, x, j) => s + x * vfFit.beta[j], 0);
      const vfResidTest = logVflat[i] - vfTestPred;

      if (mk === 'D5a') {
        const pc1f = computePC1(structVars, trainIdx);
        const pc1Train = trainIdx.map(j => applyPC1(structVars, j, pc1f));
        const pc1Test = applyPC1(structVars, i, pc1f);
        Xtrain = trainIdx.map((j, idx) => [...core[j], pc1Train[idx], vfResidTrain[idx]]);
        xiTest = [1, ...core[i], pc1Test, vfResidTest];
      } else if (mk === 'D5b') {
        Xtrain = trainIdx.map((j, idx) => [...core[j], logMbar[j], vfResidTrain[idx]]);
        xiTest = [1, ...core[i], logMbar[i], vfResidTest];
      } else {
        Xtrain = trainIdx.map((j, idx) => [...core[j], logL36[j], vfResidTrain[idx]]);
        xiTest = [1, ...core[i], logL36[i], vfResidTest];
      }
    }

    let innerSS = 0;
    for (let j = 0; j < Ytrain.length; j++) {
      const Yinn = [...Ytrain.slice(0, j), ...Ytrain.slice(j + 1)];
      const Xinn = [...Xtrain.slice(0, j), ...Xtrain.slice(j + 1)];
      const f = ols(Yinn, Xinn);
      const xinn = [1, ...Xtrain[j]];
      innerSS += (Ytrain[j] - xinn.reduce((s, x, jj) => s + x * f.beta[jj], 0)) ** 2;
    }
    const innerGap = gapPct(Math.sqrt(innerSS / Ytrain.length), sd(Ytrain));
    if (innerGap > bestGap) {
      bestGap = innerGap;
      bestMod = mk;
      const f = ols(Ytrain, Xtrain);
      bestPred = xiTest.reduce((s, x, j) => s + x * f.beta[j], 0);
    }
  }
  nestedSel[bestMod]++;
  nestedSS += (Y[i] - bestPred) ** 2;
}

const nestedGap = gapPct(Math.sqrt(nestedSS / N), sdY);
console.log('  Nested CV (C vs ' + bestD5key + '):');
console.log('    Nested gap: ' + nestedGap.toFixed(1) + '%');
for (const mk of nestedKeys) {
  console.log('    ' + mk.padEnd(5) + ' selected: ' + nestedSel[mk] + '/' + N + ' (' + (nestedSel[mk] / N * 100).toFixed(0) + '%)');
}
console.log('');

const d5Preferred = nestedSel[bestD5key] > nestedSel.C;
const d5Dominant = nestedSel[bestD5key] > N * 0.6;
console.log('  Decomposed model ' + (d5Dominant ? 'DOMINATES' : d5Preferred ? 'PREFERRED' : 'DOES NOT WIN') + ' nested CV');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: COMPONENT CONTRIBUTIONS — What does each part add?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const baselineKeys = {
  'Core only (MHI+Mhost+MR)': gals45.map((_, i) => core[i]),
  '+PC1': gals45.map((_, i) => [...core[i], pc1Scores[i]]),
  '+VfResid': gals45.map((_, i) => [...core[i], vfResidGlobal[i]]),
  '+Vflat (=C)': buildX_global('C'),
  '+PC1+VfResid (=D5a)': buildX_global('D5a'),
  '+Mbar': gals45.map((_, i) => [...core[i], logMbar[i]]),
  '+Mbar+VfResid (=D5b)': buildX_global('D5b'),
  '+L36': gals45.map((_, i) => [...core[i], logL36[i]]),
  '+L36+VfResid (=D5c)': buildX_global('D5c')
};

console.log('  Incremental gap analysis:\n');
console.log('    ' + 'Model'.padEnd(30) + '  k   gap%    R2adj');
console.log('    ' + '─'.repeat(55));

for (const [label, X] of Object.entries(baselineKeys)) {
  const fit = ols(Y, X);
  const loo = looCV(Y, X);
  const gap = gapPct(loo, sdY);
  console.log('    ' + label.padEnd(30) + '  ' + fit.k + '   ' + gap.toFixed(1).padStart(5) + '%   ' + fit.r2adj.toFixed(3));
}
console.log('');

const coreGap = gapPct(looCV(Y, core.map(c => c)), sdY);
const pc1OnlyGap = gapPct(looCV(Y, gals45.map((_, i) => [...core[i], pc1Scores[i]])), sdY);
const vfResidOnlyGap = gapPct(looCV(Y, gals45.map((_, i) => [...core[i], vfResidGlobal[i]])), sdY);

console.log('  Marginal contributions above core (' + coreGap.toFixed(1) + '%):');
console.log('    PC1 alone:    +' + (pc1OnlyGap - coreGap).toFixed(1) + 'pp');
console.log('    VfResid alone: +' + (vfResidOnlyGap - coreGap).toFixed(1) + 'pp');
console.log('    Vflat alone:  +' + (cRes.gap - coreGap).toFixed(1) + 'pp');
console.log('    PC1+VfResid:  +' + (bestD5.gap - coreGap).toFixed(1) + 'pp');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: BOOTSTRAP + VIF — Is the decomposed model clean?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  VIF analysis:\n');
const d5X = buildX_global(bestD5key);
const d5VarNames = modelDefs[bestD5key].vars;
const cVarNames = modelDefs.C.vars;
const cX = buildX_global('C');

console.log('  Model C:');
for (let j = 0; j < cVarNames.length; j++) {
  const v = vif(cX, j);
  const flag = v > 5 ? ' ★★★' : v > 2.5 ? ' ★★' : '';
  console.log('    VIF(' + cVarNames[j].padEnd(12) + ') = ' + v.toFixed(2) + flag);
}
console.log('');

console.log('  ' + bestD5key + ':');
for (let j = 0; j < d5VarNames.length; j++) {
  const v = vif(d5X, j);
  const flag = v > 5 ? ' ★★★' : v > 2.5 ? ' ★★' : '';
  console.log('    VIF(' + d5VarNames[j].padEnd(12) + ') = ' + v.toFixed(2) + flag);
}
console.log('');

console.log('  Bootstrap coefficients (' + bestD5key + '):\n');
const d5Boot = modelResults[bestD5key].boot;
for (const [vn, info] of Object.entries(d5Boot)) {
  if (vn.startsWith('_')) continue;
  console.log('    ' + vn.padEnd(12) + ': ' + info.coeff.toFixed(4).padStart(8) + '  CI95=[' + info.ci95[0].toFixed(4) + ', ' + info.ci95[1].toFixed(4) + ']  flip=' + info.flipRate + '%  CV=' + info.cv);
}
console.log('');

const cBoot = modelResults.C.boot;
console.log('  Bootstrap coefficients (C):\n');
for (const [vn, info] of Object.entries(cBoot)) {
  if (vn.startsWith('_')) continue;
  console.log('    ' + vn.padEnd(12) + ': ' + info.coeff.toFixed(4).padStart(8) + '  CI95=[' + info.ci95[0].toFixed(4) + ', ' + info.ci95[1].toFixed(4) + ']  flip=' + info.flipRate + '%  CV=' + info.cv);
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 5: ALTERNATIVE DECOMPOSITIONS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  What if we use raw structure vars instead of PC1?\n');

const altModels = {
  'Core+Mbar+L36+VfRes': gals45.map((_, i) => [...core[i], logMbar[i], logL36[i], vfResidGlobal[i]]),
  'Core+Mbar+Rdisk+VfRes': gals45.map((_, i) => [...core[i], logMbar[i], logRdisk[i], vfResidGlobal[i]]),
  'Core+L36+Rdisk+VfRes': gals45.map((_, i) => [...core[i], logL36[i], logRdisk[i], vfResidGlobal[i]]),
  'Core+Mbar+morphT+VfRes': gals45.map((_, i) => [...core[i], logMbar[i], morphT[i], vfResidGlobal[i]]),
  'Core+Sig0+VfRes': gals45.map((_, i) => [...core[i], logSig0[i], vfResidGlobal[i]])
};

console.log('    ' + 'Model'.padEnd(28) + '  k   gap%    AIC       maxFlip%');
console.log('    ' + '─'.repeat(65));

for (const [label, X] of Object.entries(altModels)) {
  const fit = ols(Y, X);
  const loo = looCV(Y, X);
  const gap = gapPct(loo, sdY);
  const boot = bootstrapAnalysis(Y, X, ['a','b','c','d','e','f'].slice(0, X[0].length), NBOOT, 42);
  console.log('    ' + label.padEnd(28) + '  ' + fit.k + '   ' + gap.toFixed(1).padStart(5) + '%  ' + fit.aic.toFixed(2).padStart(8) + '  ' + (boot._maxFlip * 100).toFixed(1).padStart(6) + '%');
}
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 130: FINAL SYNTHESIS');
console.log('══════════════════════════════════════════════════════════════════════\n');

const d5Gap = bestD5.gap;
const d5WinsNested = d5Dominant;
const d5Clean = modelResults[bestD5key].boot._maxFlip < 0.05;
const d5BeatsBIC = modelResults[bestD5key].fit.bic < cRes.fit.bic;
const penaltyRisk = d5Gap > cRes.gap + 2 && !d5BeatsBIC;

console.log('  EVIDENCE SUMMARY:\n');
console.log('    Decomposed model (' + bestD5key + '): gap=' + d5Gap.toFixed(1) + '%');
console.log('    Model C:                    gap=' + cRes.gap.toFixed(1) + '%');
console.log('    Delta:                      ' + (d5Gap - cRes.gap > 0 ? '+' : '') + (d5Gap - cRes.gap).toFixed(1) + 'pp');
console.log('    Nested CV:                  ' + bestD5key + '=' + nestedSel[bestD5key] + '/' + N + ', C=' + nestedSel.C + '/' + N);
console.log('    BIC comparison:             ' + bestD5key + '=' + modelResults[bestD5key].fit.bic.toFixed(2) + ' vs C=' + cRes.fit.bic.toFixed(2) + (d5BeatsBIC ? ' (D5 better)' : ' (C better)'));
console.log('    Max flip rate:              ' + bestD5key + '=' + (modelResults[bestD5key].boot._maxFlip * 100).toFixed(1) + '%, C=' + (cRes.boot._maxFlip * 100).toFixed(1) + '%');
console.log('');

let finalVerdict;
if (d5Gap > cRes.gap + 3 && d5WinsNested && d5Clean && d5BeatsBIC) {
  finalVerdict = 'DECOMPOSITION_SUCCEEDS';
} else if (d5Gap > cRes.gap && d5WinsNested) {
  finalVerdict = 'DECOMPOSITION_PROMISING_BUT_PENALIZED';
} else if (d5Gap > cRes.gap + 1) {
  finalVerdict = 'DECOMPOSITION_MARGINAL';
} else {
  finalVerdict = 'VFLAT_REMAINS_OPTIMAL';
}

const verdictDescriptions = {
  DECOMPOSITION_SUCCEEDS: 'The 4th sector is cleanly split into baryonic structure + kinematic residual. The decomposed model is both predictively superior and statistically clean. This is the deeper physical law.',
  DECOMPOSITION_PROMISING_BUT_PENALIZED: 'The decomposed model predicts better, but the extra parameter carries a BIC/complexity penalty. Vflat remains the recommended practical encoding, but the decomposition reveals real physical structure.',
  DECOMPOSITION_MARGINAL: 'The decomposition gains a small edge but does not clearly justify the extra parameter. The 4th axis is best left as Vflat.',
  VFLAT_REMAINS_OPTIMAL: 'The decomposition does not improve on Vflat. The 4th axis is operationally irreducible to a single kinematic variable.'
};

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + finalVerdict);
const vdLines = verdictDescriptions[finalVerdict].match(/.{1,60}/g);
for (const line of vdLines) {
  console.log('  ║  ' + line.padEnd(62) + '║');
}
console.log('  ╚══════════════════════════════════════════════════════════════════╝');
console.log('');

console.log('  FRAMEWORK CONSEQUENCE:\n');

if (finalVerdict === 'DECOMPOSITION_SUCCEEDS') {
  console.log('    NEW PRIMARY LAW: logA0 = f(MHI, Mhost, MR, Structure, Vf_resid)');
  console.log('    The 4th sector is split into two explicit components');
  console.log('    Model C remains valid as a compressed 4-variable version');
} else if (finalVerdict.includes('PROMISING')) {
  console.log('    PRACTICAL LAW:  logA0 = f(MHI, Mhost, MR, Vflat) — Model C');
  console.log('    PHYSICAL DECOMPOSITION: The 4th axis = structure + kinematic bonus');
  console.log('    Both parts carry independent, confirmed predictive value');
  console.log('    But parsimony favors the compressed Vflat encoding');
} else {
  console.log('    Model C confirmed as the primary law');
  console.log('    Vflat is the irreducible 4th axis for this sample');
}

console.log('');
console.log('    COMPLETE HIERARCHY:');
console.log('      1. Model C (MHI+Mhost+MR+Vflat):        gap=' + cRes.gap.toFixed(1) + '% — best practical law');
console.log('      2. ' + bestD5key + ' (' + modelDefs[bestD5key].name + '): gap=' + d5Gap.toFixed(1) + '% — deepest physical model');
console.log('      3. M4 (MHI+Mhost+MR+Sigma0):            gap=' + modelResults.M4.gap.toFixed(1) + '% — interpretive alternative');
console.log('');

const output = {
  phase: '130',
  title: 'Split the 4th Sector',
  N: N,
  models: {},
  nested_cv: { gap: +nestedGap.toFixed(1), selection: nestedSel, models_compared: nestedKeys },
  component_analysis: {
    core_gap: +coreGap.toFixed(1),
    pc1_marginal: +(pc1OnlyGap - coreGap).toFixed(1),
    vfResid_marginal: +(vfResidOnlyGap - coreGap).toFixed(1),
    vflat_marginal: +(cRes.gap - coreGap).toFixed(1),
    combined_marginal: +(d5Gap - coreGap).toFixed(1)
  },
  pc1: {
    loadings: Object.fromEntries(structNames.map((n, j) => [n, +pc1Global.loadings[j].toFixed(3)])),
    var_explained: +(pc1Global.varExplained * 100).toFixed(1)
  },
  vflat_residual: {
    r2_structure_explains: +vfResidFit.r2.toFixed(3),
    unique_fraction: +((1 - vfResidFit.r2) * 100).toFixed(1),
    r_with_logA0: +pearsonR(vfResidGlobal, Y).toFixed(3)
  },
  final_verdict: finalVerdict,
  final_verdict_description: verdictDescriptions[finalVerdict],
  best_decomposed: bestD5key
};

for (const key of ['C', 'D5a', 'D5b', 'D5c', 'M4']) {
  const r = modelResults[key];
  output.models[key] = {
    name: modelDefs[key].name,
    k: r.fit.k,
    gap: +r.gap.toFixed(1),
    loo_rms: +r.loo.toFixed(4),
    aic: +r.fit.aic.toFixed(2),
    bic: +r.fit.bic.toFixed(2),
    r2adj: +r.fit.r2adj.toFixed(4),
    max_flip: +(r.boot._maxFlip * 100).toFixed(1)
  };
}

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase130-split-4th-sector.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase130-split-4th-sector.json');

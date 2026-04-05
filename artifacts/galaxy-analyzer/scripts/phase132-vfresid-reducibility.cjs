#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 132: IS VFRESID REDUCIBLE TO HALO COUPLING OBSERVABLES?');
console.log('');
console.log('  Phase 131 identified VfResid as a baryon-halo coupling proxy.');
console.log('  Now we test: can haloK, mondImprove, envCode, or any bundle');
console.log('  of observables REPLACE VfResid? Or is it irreducibly the best');
console.log('  available compressed representation of coupling?');
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
  const r2 = tss > 0 ? 1 - rss / tss : 0;
  const aic = n * Math.log(rss / n) + 2 * p;
  return { beta, resid, rss, tss, r2, aic, p };
}

function looCV(Y, X) {
  const n = Y.length; let ss = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xt = [...X.slice(0, i), ...X.slice(i + 1)];
    const f = ols(Yt, Xt);
    const xi = [1, ...X[i]];
    const pred = xi.reduce((s, x, j) => s + x * f.beta[j], 0);
    ss += (Y[i] - pred) ** 2;
  }
  return Math.sqrt(ss / n);
}

function nestedCV(Y, Xbase, Xfull) {
  const n = Y.length; let wins = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xbt = [...Xbase.slice(0, i), ...Xbase.slice(i + 1)];
    const Xft = [...Xfull.slice(0, i), ...Xfull.slice(i + 1)];
    const fb = ols(Yt, Xbt);
    const ff = ols(Yt, Xft);
    const xbi = [1, ...Xbase[i]];
    const xfi = [1, ...Xfull[i]];
    const pb = xbi.reduce((s, x, j) => s + x * fb.beta[j], 0);
    const pf = xfi.reduce((s, x, j) => s + x * ff.beta[j], 0);
    if (Math.abs(Y[i] - pf) <= Math.abs(Y[i] - pb)) wins++;
  }
  return wins;
}

function gapPct(rms, sdy) { return 100 * (1 - rms ** 2 / sdy ** 2); }

function bootstrap(Y, X, nBoot) {
  const n = Y.length;
  const refRms = looCV(Y, X);
  let worse = 0;
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(Math.random() * n));
    const Yb = idx.map(i => Y[i]);
    const Xb = idx.map(i => X[i]);
    try {
      const rms = looCV(Yb, Xb);
      if (rms > refRms * 1.5) worse++;
    } catch (e) { worse++; }
  }
  return worse / nBoot;
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
const envCode = gals45.map(g => g.envCode);
const logVmaxToVflat = gals45.map(g => {
  const r = resMap[g.name];
  return Math.log10((r ? r.maxV : sparcMap[g.name].Vflat) / sparcMap[g.name].Vflat);
});
const logExtent = gals45.map(g => {
  const r = resMap[g.name];
  return Math.log10((r ? r.maxR : 10) / sparcMap[g.name].Rdisk);
});
const concentration = gals45.map(g => Math.log10(sparcMap[g.name].Reff / sparcMap[g.name].Rdisk));
const logSBeff = gals45.map(g => Math.log10(Math.max(sparcMap[g.name].SBeff, 0.01)));
const rcWig = gals45.map(g => g.rcWiggliness);

const haloK = gals45.map(g => {
  const r = resMap[g.name];
  return r && r.models && r.models.dark_halo_linear ? Math.log10(Math.max(r.models.dark_halo_linear.k, 1)) : 0;
});
const mondImprove = gals45.map(g => {
  const r = resMap[g.name];
  return r && r.models && r.models.mond ? r.models.mond.improvementVsNewton : 0;
});

const structX = gals45.map((_, i) => [logMbar[i], logL36[i], logRdisk[i], morphT[i]]);
const vfResidFit = ols(logVflat, structX);
const VfResid = vfResidFit.resid;

console.log('  VfResid: SD = ' + sd(VfResid).toFixed(4) + ', r(logA0) = ' + pearsonR(VfResid, Y).toFixed(3));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST A: PARTIAL CORRELATIONS — Does haloK survive controls?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function partialR(x, y, controls) {
  const rx = ols(x, controls).resid;
  const ry = ols(y, controls).resid;
  return pearsonR(rx, ry);
}

const core3X = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]);
const structOnly = gals45.map((_, i) => [logMbar[i], logL36[i], logRdisk[i], morphT[i]]);
const envOnly = gals45.map((_, i) => [envCode[i]]);
const mrOnly = gals45.map((_, i) => [logMR[i]]);
const coreStruct = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], logMbar[i], logL36[i], logRdisk[i], morphT[i]]);

const partialTests = [
  { name: 'raw (no control)', ctrl: null },
  { name: 'control: Core (MHI+Mhost+MR)', ctrl: core3X },
  { name: 'control: Structure (Mbar+L36+Rdisk+T)', ctrl: structOnly },
  { name: 'control: envCode', ctrl: envOnly },
  { name: 'control: logMeanRun', ctrl: mrOnly },
  { name: 'control: Core+Structure (all 7)', ctrl: coreStruct }
];

const partialTargets = [
  { name: 'haloK', arr: haloK },
  { name: 'mondImprove', arr: mondImprove },
  { name: 'envCode', arr: envCode },
  { name: 'logMeanRun', arr: logMR },
  { name: 'concentration', arr: concentration },
  { name: 'logSBeff', arr: logSBeff }
];

console.log('  ' + 'Control'.padEnd(40) + partialTargets.map(t => t.name.padStart(12)).join(''));
console.log('  ' + '─'.repeat(40 + partialTargets.length * 12));

for (const pt of partialTests) {
  let line = '  ' + pt.name.padEnd(40);
  for (const tgt of partialTargets) {
    let r;
    if (!pt.ctrl) {
      r = pearsonR(VfResid, tgt.arr);
    } else {
      r = partialR(VfResid, tgt.arr, pt.ctrl);
    }
    const flag = Math.abs(r) > 0.3 ? '★' : ' ';
    line += (r.toFixed(3) + flag).padStart(12);
  }
  console.log(line);
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST B: REPLACEMENT — Can haloK bundles replace VfResid?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const replacementModels = [
  { name: 'Core only (baseline)', X: core3X },
  { name: 'Core + VfResid', X: gals45.map((_, i) => [...core3X[i], VfResid[i]]) },
  { name: 'Core + Vflat (Model C)', X: gals45.map((_, i) => [...core3X[i], logVflat[i]]) },
  { name: 'Core + haloK', X: gals45.map((_, i) => [...core3X[i], haloK[i]]) },
  { name: 'Core + mondImprove', X: gals45.map((_, i) => [...core3X[i], mondImprove[i]]) },
  { name: 'Core + haloK + envCode', X: gals45.map((_, i) => [...core3X[i], haloK[i], envCode[i]]) },
  { name: 'Core + haloK + mondImprove', X: gals45.map((_, i) => [...core3X[i], haloK[i], mondImprove[i]]) },
  { name: 'Core + haloK + env + mondI', X: gals45.map((_, i) => [...core3X[i], haloK[i], envCode[i], mondImprove[i]]) },
  { name: 'Core + haloK + conc + SBeff', X: gals45.map((_, i) => [...core3X[i], haloK[i], concentration[i], logSBeff[i]]) },
  { name: 'Core + Sig0 (M4)', X: gals45.map((_, i) => [...core3X[i], logSig0[i]]) },
  { name: 'Core + Sig0 + VfResid', X: gals45.map((_, i) => [...core3X[i], logSig0[i], VfResid[i]]) },
  { name: 'Core + logSBeff', X: gals45.map((_, i) => [...core3X[i], logSBeff[i]]) }
];

const coreGap = gapPct(looCV(Y, core3X), sdY);
const vfResidGap = gapPct(looCV(Y, replacementModels[1].X), sdY);

console.log('  ' + 'Model'.padEnd(35) + '  gap%  delta  nested  AIC');
console.log('  ' + '─'.repeat(72));

for (const rm of replacementModels) {
  const rms = looCV(Y, rm.X);
  const gap = gapPct(rms, sdY);
  const delta = gap - coreGap;
  const nested = rm.X[0].length > core3X[0].length ? nestedCV(Y, core3X, rm.X) : N;
  const fit = ols(Y, rm.X);
  console.log('  ' + rm.name.padEnd(35) + gap.toFixed(1).padStart(6) + '%  ' + (delta > 0 ? '+' : '') + delta.toFixed(1).padStart(5) + 'pp  ' + (nested < N ? nested + '/' + N : '  —  ') + '  ' + fit.aic.toFixed(1));
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST C: MULTIVARIATE PROBE BUNDLE — Can probes together match VfResid?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  Strategy: Build a synthetic "coupling proxy" from the best halo');
console.log('  correlates and test whether it can match VfResid predictively.\n');

const bundleTests = [
  { name: 'haloK alone', vars: gals45.map((_, i) => [haloK[i]]) },
  { name: 'haloK + envCode', vars: gals45.map((_, i) => [haloK[i], envCode[i]]) },
  { name: 'haloK + mondImprove', vars: gals45.map((_, i) => [haloK[i], mondImprove[i]]) },
  { name: 'haloK + envCode + mondImprove', vars: gals45.map((_, i) => [haloK[i], envCode[i], mondImprove[i]]) },
  { name: 'haloK + env + mondI + conc', vars: gals45.map((_, i) => [haloK[i], envCode[i], mondImprove[i], concentration[i]]) },
  { name: 'haloK + env + mondI + conc + SBeff', vars: gals45.map((_, i) => [haloK[i], envCode[i], mondImprove[i], concentration[i], logSBeff[i]]) },
  { name: 'all 6 probes + Vmax/Vf + rcWig', vars: gals45.map((_, i) => [haloK[i], envCode[i], mondImprove[i], concentration[i], logSBeff[i], logVmaxToVflat[i], rcWig[i]]) }
];

console.log('  How much of VfResid variance can these bundles explain? (in-sample R2)\n');

for (const bt of bundleTests) {
  const fit = ols(VfResid, bt.vars);
  const rms = looCV(VfResid, bt.vars);
  const looR2 = 1 - rms ** 2 / (sd(VfResid) ** 2);
  console.log('  ' + bt.name.padEnd(42) + '  R2=' + fit.r2.toFixed(3) + '  LOO-R2=' + (looR2 > 0 ? looR2.toFixed(3) : '<0') + '  (' + (fit.r2 * 100).toFixed(1) + '% in-sample)');
}
console.log('');

console.log('  Now: can any bundle, when added to Core, match VfResid predictive gap?\n');
console.log('  ' + 'Bundle (added to Core)'.padEnd(42) + '  gap%    vs VfResid(' + vfResidGap.toFixed(1) + '%)');
console.log('  ' + '─'.repeat(68));

const bundleModels = [
  { name: 'VfResid alone (target)', X: gals45.map((_, i) => [...core3X[i], VfResid[i]]) },
  { name: 'haloK', X: gals45.map((_, i) => [...core3X[i], haloK[i]]) },
  { name: 'haloK + envCode', X: gals45.map((_, i) => [...core3X[i], haloK[i], envCode[i]]) },
  { name: 'haloK + mondImprove', X: gals45.map((_, i) => [...core3X[i], haloK[i], mondImprove[i]]) },
  { name: 'haloK + env + mondI', X: gals45.map((_, i) => [...core3X[i], haloK[i], envCode[i], mondImprove[i]]) },
  { name: 'haloK + env + mondI + conc', X: gals45.map((_, i) => [...core3X[i], haloK[i], envCode[i], mondImprove[i], concentration[i]]) },
  { name: 'haloK + env + mondI + conc + SBeff', X: gals45.map((_, i) => [...core3X[i], haloK[i], envCode[i], mondImprove[i], concentration[i], logSBeff[i]]) },
  { name: 'ALL 7 probes', X: gals45.map((_, i) => [...core3X[i], haloK[i], envCode[i], mondImprove[i], concentration[i], logSBeff[i], logVmaxToVflat[i], rcWig[i]]) }
];

for (const bm of bundleModels) {
  const gap = gapPct(looCV(Y, bm.X), sdY);
  const deficit = gap - vfResidGap;
  console.log('  ' + bm.name.padEnd(42) + gap.toFixed(1).padStart(6) + '%    ' + (deficit > 0 ? '+' : '') + deficit.toFixed(1) + 'pp');
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST D: RESIDUALIZATION — Does VfResid survive after removing');
console.log('  the haloK-explainable portion?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const bestBundle = gals45.map((_, i) => [haloK[i], envCode[i], mondImprove[i], concentration[i], logSBeff[i]]);
const vfResidAfterBundle = ols(VfResid, bestBundle).resid;

console.log('  VfResid after removing best 5-probe bundle:');
console.log('    SD remaining = ' + sd(vfResidAfterBundle).toFixed(4) + ' (from ' + sd(VfResid).toFixed(4) + ')');
console.log('    Fraction remaining = ' + (sd(vfResidAfterBundle) ** 2 / sd(VfResid) ** 2 * 100).toFixed(1) + '%');
console.log('    r(residual, logA0) = ' + pearsonR(vfResidAfterBundle, Y).toFixed(3));
console.log('');

const coreVfResidRemaining = gals45.map((_, i) => [...core3X[i], vfResidAfterBundle[i]]);
const gapRemaining = gapPct(looCV(Y, coreVfResidRemaining), sdY);
console.log('    Core + VfResid_remaining gap = ' + gapRemaining.toFixed(1) + '% (vs VfResid full ' + vfResidGap.toFixed(1) + '%)');
console.log('    Delta = ' + (gapRemaining - coreGap).toFixed(1) + 'pp above core');
console.log('');

if (gapRemaining - coreGap > 3.0) {
  console.log('    ⚠ VfResid SURVIVES even after removing best bundle content!');
  console.log('      There is irreducible information in VfResid beyond haloK+env+mondI+conc+SBeff.');
} else {
  console.log('    ✓ VfResid is largely absorbed by the probe bundle.');
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST E: FOLD-INTERNAL CHECK — Honest VfResid predictive power');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  VfResid is computed from a full-sample OLS on structure vars.');
console.log('  To be rigorous, we recompute VfResid fold-internally in LOO.\n');

let foldInternalSS = 0;
for (let i = 0; i < N; i++) {
  const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
  const VfTrain = [...logVflat.slice(0, i), ...logVflat.slice(i + 1)];
  const StructTrain = [...structX.slice(0, i), ...structX.slice(i + 1)];
  const CoreTrain = [...core3X.slice(0, i), ...core3X.slice(i + 1)];

  const vfFit = ols(VfTrain, StructTrain);
  const vfPred_i = [1, ...structX[i]].reduce((s, x, j) => s + x * vfFit.beta[j], 0);
  const vfResid_i = logVflat[i] - vfPred_i;

  const trainResids = VfTrain.map((v, j) => v - [1, ...StructTrain[j]].reduce((s, x, k) => s + x * vfFit.beta[k], 0));

  const coreResidX_train = CoreTrain.map((c, j) => [...c, trainResids[j]]);
  const coreResidX_test = [...core3X[i], vfResid_i];

  const mainFit = ols(Yt, coreResidX_train);
  const pred = [1, ...coreResidX_test].reduce((s, x, j) => s + x * mainFit.beta[j], 0);
  foldInternalSS += (Y[i] - pred) ** 2;
}

const foldInternalRMS = Math.sqrt(foldInternalSS / N);
const foldInternalGap = gapPct(foldInternalRMS, sdY);
console.log('  Fold-internal Core+VfResid gap = ' + foldInternalGap.toFixed(1) + '%');
console.log('  Standard (full-sample VfResid)  = ' + vfResidGap.toFixed(1) + '%');
console.log('  Shrinkage from fold-internal    = ' + (vfResidGap - foldInternalGap).toFixed(1) + 'pp');
console.log('  Still above core (' + coreGap.toFixed(1) + '%) by ' + (foldInternalGap - coreGap).toFixed(1) + 'pp');
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 132: VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const bestBundleGap = gapPct(looCV(Y, gals45.map((_, i) => [...core3X[i], haloK[i], envCode[i], mondImprove[i], concentration[i], logSBeff[i]])), sdY);
const deficit = vfResidGap - bestBundleGap;
const foldInternalDelta = foldInternalGap - coreGap;

let verdict;
if (deficit > 5 && foldInternalDelta > 10) {
  verdict = 'IRREDUCIBLE';
} else if (deficit > 2 && foldInternalDelta > 5) {
  verdict = 'PARTIALLY_REDUCIBLE';
} else if (deficit <= 2 && bestBundleGap >= vfResidGap - 2) {
  verdict = 'REDUCIBLE';
} else {
  verdict = 'AMBIGUOUS';
}

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + verdict.padEnd(55) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════════╝\n');

if (verdict === 'IRREDUCIBLE') {
  console.log('  VfResid is the BEST AVAILABLE compressed proxy for baryon-halo');
  console.log('  coupling. No bundle of catalog observables (haloK, envCode,');
  console.log('  mondImprove, concentration, SBeff) can replicate its predictive');
  console.log('  power. The gap deficit vs best bundle = ' + deficit.toFixed(1) + 'pp.');
  console.log('');
  console.log('  Fold-internal VfResid still adds +' + foldInternalDelta.toFixed(1) + 'pp above core,');
  console.log('  confirming this is genuine signal, not LOO leakage.');
  console.log('');
  console.log('  PRACTICAL IMPLICATION:');
  console.log('    Vflat remains the operationally best 4th axis because:');
  console.log('    (a) It is a standard observable (no model fitting needed)');
  console.log('    (b) The coupling information it carries is IRREDUCIBLE to');
  console.log('        simpler catalog quantities');
  console.log('    (c) Its physical meaning is now identified: the excess');
  console.log('        dynamical support beyond baryonic structure = baryon-halo');
  console.log('        coupling efficiency');
} else if (verdict === 'PARTIALLY_REDUCIBLE') {
  console.log('  VfResid is partly captured by halo-coupling observables but');
  console.log('  retains irreducible content beyond them. Best bundle reaches');
  console.log('  ' + bestBundleGap.toFixed(1) + '% vs VfResid ' + vfResidGap.toFixed(1) + '% (deficit ' + deficit.toFixed(1) + 'pp).');
  console.log('  Fold-internal VfResid adds +' + foldInternalDelta.toFixed(1) + 'pp above core.');
} else if (verdict === 'REDUCIBLE') {
  console.log('  VfResid can be adequately replaced by a probe bundle.');
  console.log('  Best bundle gap = ' + bestBundleGap.toFixed(1) + '% ~ VfResid ' + vfResidGap.toFixed(1) + '%.');
}

console.log('');
console.log('  HIERARCHY UPDATE:\n');
console.log('    PRIMARY: Model C = Core + Vflat    gap = 52.1%');
console.log('    DECODED: Core + VfResid            gap = ' + vfResidGap.toFixed(1) + '% (fold-internal: ' + foldInternalGap.toFixed(1) + '%)');
console.log('    BEST BUNDLE: Core + 5 probes       gap = ' + bestBundleGap.toFixed(1) + '%');
console.log('    M4: Core + Sig0                    gap = 47.4%');
console.log('    Core baseline                      gap = ' + coreGap.toFixed(1) + '%');
console.log('');

const output = {
  phase: '132',
  title: 'Is VfResid Reducible to Halo Coupling Observables?',
  verdict: verdict,
  vfresid_gap: +vfResidGap.toFixed(1),
  vfresid_fold_internal_gap: +foldInternalGap.toFixed(1),
  core_gap: +coreGap.toFixed(1),
  best_bundle_gap: +bestBundleGap.toFixed(1),
  deficit_vs_best_bundle: +deficit.toFixed(1),
  fold_internal_delta_vs_core: +foldInternalDelta.toFixed(1),
  partial_correlations: {},
  replacement_models: {},
  bundle_r2_for_vfresid: {},
  vfresid_after_bundle_removal: {
    fraction_remaining_pct: +(sd(vfResidAfterBundle) ** 2 / sd(VfResid) ** 2 * 100).toFixed(1),
    r_with_logA0: +pearsonR(vfResidAfterBundle, Y).toFixed(3),
    gap_with_core: +gapRemaining.toFixed(1)
  }
};

for (const pt of partialTests) {
  const row = {};
  for (const tgt of partialTargets) {
    row[tgt.name] = pt.ctrl ? +partialR(VfResid, tgt.arr, pt.ctrl).toFixed(3) : +pearsonR(VfResid, tgt.arr).toFixed(3);
  }
  output.partial_correlations[pt.name] = row;
}

for (const rm of replacementModels) {
  output.replacement_models[rm.name] = +gapPct(looCV(Y, rm.X), sdY).toFixed(1);
}

for (const bt of bundleTests) {
  output.bundle_r2_for_vfresid[bt.name] = +ols(VfResid, bt.vars).r2.toFixed(3);
}

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase132-vfresid-reducibility.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase132-vfresid-reducibility.json');

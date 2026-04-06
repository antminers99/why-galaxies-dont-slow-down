#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 134: EXTERNAL / TRANSFER VALIDATION');
console.log('');
console.log('  Train on N=45 published, predict N=10 crude holdout.');
console.log('  Does the hierarchical coupling law transfer outside the sample?');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'stage-A-master-table.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const sparcTable = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-results.json'), 'utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparcTable.forEach(s => { sparcMap[s.name] = s; });
const resMap = {}; sparcResults.perGalaxy.forEach(g => { if (g.rawName) resMap[g.rawName] = g; });

const allGals = stageA.galaxies;
const pubGals = allGals.filter(g => pubNames.has(g.name));
const crudeGals = allGals.filter(g => !pubNames.has(g.name) && sparcMap[g.name] && sparcMap[g.name].Vflat > 0);
const N_train = pubGals.length;
const N_test = crudeGals.length;

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

function transferGap(rmse, baselineRmse) {
  return baselineRmse > 0 ? 100 * (1 - rmse ** 2 / baselineRmse ** 2) : 0;
}

function getGalData(g) {
  const s = sparcMap[g.name]; const r = resMap[g.name]; const td = tdMap[g.name];
  return {
    name: g.name, logA0: g.logA0,
    logMHI: g.logMHI, logMhost: td.logMhost, logMR: g.logMeanRun,
    logVflat: Math.log10(s.Vflat), Vflat: s.Vflat,
    logMbar: Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9),
    logL36: Math.log10(Math.max(s.L36, 0.01)),
    logRdisk: Math.log10(s.Rdisk), morphT: s.T,
    logSig0: g.logSigma0, envCode: g.envCode,
    lhOuter: r.models.log_halo.outerImprovement,
    haloK: Math.log10(Math.max(r.models.dark_halo_linear.k, 1)),
    mondImprove: r.models.mond.improvementVsNewton
  };
}

const train = pubGals.map(getGalData);
const test = crudeGals.map(getGalData);

console.log('  Training set: N=' + N_train + ' (published)');
console.log('  Test set: N=' + N_test + ' (crude holdout)');
console.log('  Test Vflat range: ' + Math.min(...test.map(g => g.Vflat)).toFixed(0) + ' - ' + Math.max(...test.map(g => g.Vflat)).toFixed(0));
console.log('  Test Vflat >= 120: ' + test.filter(g => g.Vflat >= 120).length);
console.log('  Test Vflat >= 180: ' + test.filter(g => g.Vflat >= 180).length);
console.log('  NOTE: All 10 crude galaxies have Vflat >= 120 (no low-Vflat holdout).\n');

const Y_tr = train.map(g => g.logA0);
const Y_te = test.map(g => g.logA0);
const sdY_te = sd(Y_te);
const meanY_tr = mean(Y_tr);
const meanY_te = mean(Y_te);

const naiveRMSE = Math.sqrt(Y_te.reduce((s, y) => s + (y - meanY_tr) ** 2, 0) / Y_te.length);

console.log('  Baselines:');
console.log('    Training mean logA0: ' + meanY_tr.toFixed(4));
console.log('    Test mean logA0: ' + meanY_te.toFixed(4) + '  (delta=' + (meanY_te - meanY_tr).toFixed(4) + ')');
console.log('    Naive RMSE (predict training mean): ' + naiveRMSE.toFixed(4));
console.log('    Test SD: ' + sdY_te.toFixed(4) + '\n');

const structModel = ols(
  train.map(g => g.logVflat),
  train.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT])
);

function computeVfResid(gal) {
  const pred = [1, gal.logMbar, gal.logL36, gal.logRdisk, gal.morphT]
    .reduce((s, x, j) => s + x * structModel.beta[j], 0);
  return gal.logVflat - pred;
}

const VfResid_tr = train.map(computeVfResid);
const VfResid_te = test.map(computeVfResid);

const trLogMbar = train.map(g => g.logMbar);
const teLogMbar = test.map(g => g.logMbar);
const trLogVf = train.map(g => g.logVflat);
const teLogVf = test.map(g => g.logVflat);

const extrapolating = test.filter(g =>
  g.logMbar > Math.max(...trLogMbar) * 1.01 ||
  g.logVflat > Math.max(...trLogVf) * 1.01 ||
  g.logMbar < Math.min(...trLogMbar) * 0.99 ||
  g.logVflat < Math.min(...trLogVf) * 0.99
);
if (extrapolating.length > 0) {
  console.log('  EXTRAPOLATION WARNING: ' + extrapolating.length + ' test galaxies outside training range:');
  extrapolating.forEach(g => console.log('    ' + g.name + ' (logMbar=' + g.logMbar.toFixed(2) + ', logVflat=' + g.logVflat.toFixed(3) + ')'));
  console.log('  Training ranges: logMbar [' + Math.min(...trLogMbar).toFixed(2) + ', ' + Math.max(...trLogMbar).toFixed(2) + ']' +
    '  logVflat [' + Math.min(...trLogVf).toFixed(3) + ', ' + Math.max(...trLogVf).toFixed(3) + ']\n');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: FULL MODEL TRANSFER — Train on published, predict crude');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const core_tr = train.map(g => [g.logMHI, g.logMhost, g.logMR]);
const core_te = test.map(g => [g.logMHI, g.logMhost, g.logMR]);
const coreVflat_tr = train.map((g, i) => [...core_tr[i], g.logVflat]);
const coreVflat_te = test.map((g, i) => [...core_te[i], g.logVflat]);
const coreVfResid_tr = train.map((g, i) => [...core_tr[i], VfResid_tr[i]]);
const coreVfResid_te = test.map((g, i) => [...core_te[i], VfResid_te[i]]);
const core5_tr = train.map((g, i) => [...core_tr[i], VfResid_tr[i], g.lhOuter]);
const core5_te = test.map((g, i) => [...core_te[i], VfResid_te[i], g.lhOuter]);

function computeTransfer(Y_tr, X_tr, X_te, Y_te) {
  const model = ols(Y_tr, X_tr);
  const preds = X_te.map(x => [1, ...x].reduce((s, v, j) => s + v * model.beta[j], 0));
  const errs = Y_te.map((y, i) => y - preds[i]);
  const rmse = Math.sqrt(errs.reduce((s, e) => s + e * e, 0) / Y_te.length);
  return { rmse, preds, errs, model };
}

const coreResult = computeTransfer(Y_tr, core_tr, core_te, Y_te);
const vflatResult = computeTransfer(Y_tr, coreVflat_tr, coreVflat_te, Y_te);
const vfResidResult = computeTransfer(Y_tr, coreVfResid_tr, coreVfResid_te, Y_te);
const fiveResult = computeTransfer(Y_tr, core5_tr, core5_te, Y_te);

console.log('  Model'.padEnd(35) + 'RMSE     Transfer-Gap%');
console.log('  ' + '─'.repeat(55));
console.log('  ' + 'Naive (predict train mean)'.padEnd(35) + naiveRMSE.toFixed(4) + '      0.0%');
console.log('  ' + 'Core (3-axis)'.padEnd(35) + coreResult.rmse.toFixed(4) + '     ' + transferGap(coreResult.rmse, naiveRMSE).toFixed(1) + '%');
console.log('  ' + 'Core + Vflat (Model C)'.padEnd(35) + vflatResult.rmse.toFixed(4) + '     ' + transferGap(vflatResult.rmse, naiveRMSE).toFixed(1) + '%');
console.log('  ' + 'Core + VfResid'.padEnd(35) + vfResidResult.rmse.toFixed(4) + '     ' + transferGap(vfResidResult.rmse, naiveRMSE).toFixed(1) + '%');
console.log('  ' + 'Core + VfResid + lhOuter (5-ax)'.padEnd(35) + fiveResult.rmse.toFixed(4) + '     ' + transferGap(fiveResult.rmse, naiveRMSE).toFixed(1) + '%');
console.log('');

console.log('  Per-galaxy predictions (Core + VfResid):');
console.log('  ' + 'Galaxy'.padEnd(15) + 'Actual  Pred    Error   Vflat  VfResid');
console.log('  ' + '─'.repeat(65));
for (let i = 0; i < N_test; i++) {
  const err = vfResidResult.errs[i];
  const flag = Math.abs(err) > 0.3 ? '!' : Math.abs(err) > 0.15 ? '*' : ' ';
  console.log('  ' + flag + ' ' + test[i].name.padEnd(14) +
    Y_te[i].toFixed(3).padStart(6) + '  ' + vfResidResult.preds[i].toFixed(3).padStart(6) + '  ' +
    (err > 0 ? '+' : '') + err.toFixed(3).padStart(6) + '   ' + test[i].Vflat.toFixed(0).padStart(3) +
    '   ' + VfResid_te[i].toFixed(4));
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: REGIME-SPECIFIC TRAINING TRANSFER');
console.log('  (Note: all 10 holdout galaxies have Vflat >= 120, so the test is');
console.log('   whether regime-restricted TRAINING improves transfer quality.)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const highTrain = train.filter(g => g.Vflat >= 120);
const Y_htr = highTrain.map(g => g.logA0);
const meanY_htr = mean(Y_htr);
const naiveHRMSE = Math.sqrt(Y_te.reduce((s, y) => s + (y - meanY_htr) ** 2, 0) / Y_te.length);

const core_htr = highTrain.map(g => [g.logMHI, g.logMhost, g.logMR]);
const vflat_htr = highTrain.map(g => [g.logMHI, g.logMhost, g.logMR, g.logVflat]);
const VfR_htr = highTrain.map(computeVfResid);
const vfr_htr = highTrain.map((g, i) => [g.logMHI, g.logMhost, g.logMR, VfR_htr[i]]);

const hCoreR = computeTransfer(Y_htr, core_htr, core_te, Y_te);
const hVflatR = computeTransfer(Y_htr, vflat_htr, coreVflat_te, Y_te);
const hVfResidR = computeTransfer(Y_htr, vfr_htr, coreVfResid_te, Y_te);

console.log('  Regime-trained (Vflat>=120, N=' + highTrain.length + ') vs full-trained (N=45):');
console.log('');
console.log('  ' + 'Model'.padEnd(25) + 'Full(N=45)  Regime(N=' + highTrain.length + ')   Better?');
console.log('  ' + '─'.repeat(60));
console.log('  ' + 'Core RMSE'.padEnd(25) +
  coreResult.rmse.toFixed(4).padStart(8) + '    ' + hCoreR.rmse.toFixed(4).padStart(8) +
  '     ' + (hCoreR.rmse < coreResult.rmse ? 'Regime' : 'Full'));
console.log('  ' + 'Core+Vflat RMSE'.padEnd(25) +
  vflatResult.rmse.toFixed(4).padStart(8) + '    ' + hVflatR.rmse.toFixed(4).padStart(8) +
  '     ' + (hVflatR.rmse < vflatResult.rmse ? 'Regime' : 'Full'));
console.log('  ' + 'Core+VfResid RMSE'.padEnd(25) +
  vfResidResult.rmse.toFixed(4).padStart(8) + '    ' + hVfResidR.rmse.toFixed(4).padStart(8) +
  '     ' + (hVfResidR.rmse < vfResidResult.rmse ? 'Regime' : 'Full'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: SIGN & COEFFICIENT STABILITY');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const trainVfResidModel = ols(Y_tr, coreVfResid_tr);

const allGalsData = [...train, ...test];
const Y_all = allGalsData.map(g => g.logA0);
const VfR_all = allGalsData.map(computeVfResid);
const vfr_all = allGalsData.map((g, i) => [g.logMHI, g.logMhost, g.logMR, VfR_all[i]]);
const fullVfResidModel = ols(Y_all, vfr_all);

const labels = ['intercept', 'MHI', 'Mhost', 'MR', 'VfResid'];
console.log('  Core+VfResid coefficients:');
console.log('  ' + 'Coeff'.padEnd(12) + 'Pub(N=45)'.padStart(10) + 'Comb(N=55)'.padStart(12) + '  Sign');
console.log('  ' + '─'.repeat(40));
let allSignsSame = true;
for (let j = 0; j < labels.length; j++) {
  const a = trainVfResidModel.beta[j], b = fullVfResidModel.beta[j];
  const same = Math.sign(a) === Math.sign(b);
  if (!same) allSignsSame = false;
  console.log('  ' + labels[j].padEnd(12) + a.toFixed(4).padStart(10) + b.toFixed(4).padStart(12) + '  ' + (same ? 'SAME' : 'CHANGED'));
}
console.log('');
console.log('  All signs preserved: ' + (allSignsSame ? 'YES' : 'NO'));

const trainVflatModel = ols(Y_tr, coreVflat_tr);
const vflat_all = allGalsData.map(g => [g.logMHI, g.logMhost, g.logMR, g.logVflat]);
const fullVflatModel = ols(Y_all, vflat_all);
const signsSameC = trainVflatModel.beta.slice(1).every((b, j) =>
  Math.sign(b) === Math.sign(fullVflatModel.beta[j + 1]));
console.log('  Model C signs preserved: ' + (signsSameC ? 'YES' : 'NO'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: COMBINED SAMPLE LOO — Does adding crude improve?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function looRMS(Y, X) {
  const n = Y.length; let ss = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xt = [...X.slice(0, i), ...X.slice(i + 1)];
    try {
      const f = ols(Yt, Xt);
      ss += (Y[i] - [1, ...X[i]].reduce((s, x, j) => s + x * f.beta[j], 0)) ** 2;
    } catch (e) { ss += (Y[i] - mean(Yt)) ** 2; }
  }
  return Math.sqrt(ss / n);
}

function gapPctLOO(rms, sdy) { return sdy > 0 ? 100 * (1 - rms ** 2 / sdy ** 2) : 0; }

const sdY_all = sd(Y_all);
const sdY_p = sd(Y_tr);
const core_all = allGalsData.map(g => [g.logMHI, g.logMhost, g.logMR]);
const core5_all = allGalsData.map((g, i) => [g.logMHI, g.logMhost, g.logMR, VfR_all[i], g.lhOuter]);

const looCoreP = looRMS(Y_tr, core_tr);
const looVflatP = looRMS(Y_tr, coreVflat_tr);
const looVfResidP = looRMS(Y_tr, coreVfResid_tr);

const looCoreAll = looRMS(Y_all, core_all);
const looVflatAll = looRMS(Y_all, vflat_all);
const looVfResidAll = looRMS(Y_all, vfr_all);
const loo5All = looRMS(Y_all, core5_all);

console.log('  NOTE: Combined LOO uses VfResid computed from N=45 structModel.');
console.log('  For published galaxies this is in-sample (minor leakage). Known limitation.\n');
console.log('  ' + 'Model'.padEnd(30) + '  Pub(N=45)  Combined(N=55)');
console.log('  ' + '─'.repeat(55));
console.log('  ' + 'Core gap'.padEnd(30) + gapPctLOO(looCoreP, sdY_p).toFixed(1).padStart(8) + '%   ' + gapPctLOO(looCoreAll, sdY_all).toFixed(1).padStart(8) + '%');
console.log('  ' + 'Core+Vflat gap'.padEnd(30) + gapPctLOO(looVflatP, sdY_p).toFixed(1).padStart(8) + '%   ' + gapPctLOO(looVflatAll, sdY_all).toFixed(1).padStart(8) + '%');
console.log('  ' + 'Core+VfResid gap'.padEnd(30) + gapPctLOO(looVfResidP, sdY_p).toFixed(1).padStart(8) + '%   ' + gapPctLOO(looVfResidAll, sdY_all).toFixed(1).padStart(8) + '%');
console.log('  ' + 'Core+VfResid+lhOuter gap'.padEnd(30) + '     N/A' + '   ' + gapPctLOO(loo5All, sdY_all).toFixed(1).padStart(8) + '%');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 5: VFRESID ON HOLDOUT — Does the coupling signal exist?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  VfResid on holdout galaxies:');
console.log('  ' + 'Galaxy'.padEnd(15) + 'Vflat  VfResid  logA0   haloK');
console.log('  ' + '─'.repeat(55));
for (let i = 0; i < N_test; i++) {
  const g = test[i];
  console.log('  ' + g.name.padEnd(15) + g.Vflat.toFixed(0).padStart(5) + '  ' +
    VfResid_te[i].toFixed(4).padStart(7) + '  ' + g.logA0.toFixed(3).padStart(6) + '  ' + g.haloK.toFixed(2));
}
console.log('');
const rVfResidHoldout = pearsonR(VfResid_te, Y_te);
const rHaloKHoldout = pearsonR(test.map(g => g.haloK), Y_te);
console.log('  r(VfResid, logA0) on holdout: ' + rVfResidHoldout.toFixed(3));
console.log('  r(haloK, logA0) on holdout: ' + rHaloKHoldout.toFixed(3));
console.log('  r(VfResid, haloK) on holdout: ' + pearsonR(VfResid_te, test.map(g => g.haloK)).toFixed(3));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 6: lh_outerImprove ON HOLDOUT — Does the 5th axis transfer?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const lhOuter_te = test.map(g => g.lhOuter);
console.log('  r(lhOuter, logA0) on holdout: ' + pearsonR(lhOuter_te, Y_te).toFixed(3));
const deltaGap5 = transferGap(fiveResult.rmse, naiveRMSE) - transferGap(vfResidResult.rmse, naiveRMSE);
console.log('  Transfer:');
console.log('    Core+VfResid RMSE (holdout): ' + vfResidResult.rmse.toFixed(4));
console.log('    Core+VfResid+lhOuter RMSE:   ' + fiveResult.rmse.toFixed(4));
console.log('    5-axis improves: ' + (fiveResult.rmse < vfResidResult.rmse ? 'YES' : 'NO') +
  ' (delta=' + deltaGap5.toFixed(1) + 'pp)');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 7: CROSS-VALIDATION SANITY (k-fold within published)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function kfoldRMSE(Y, X, k) {
  const n = Y.length;
  const idx = Array.from({ length: n }, (_, i) => i);
  for (let i = n - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [idx[i], idx[j]] = [idx[j], idx[i]];
  }
  const foldSize = Math.floor(n / k);
  let ss = 0, count = 0;
  for (let f = 0; f < k; f++) {
    const testIdx = new Set(idx.slice(f * foldSize, (f + 1) * foldSize));
    const Yt = [], Xt = [], Yv = [], Xv = [];
    for (let i = 0; i < n; i++) {
      if (testIdx.has(i)) { Yv.push(Y[i]); Xv.push(X[i]); }
      else { Yt.push(Y[i]); Xt.push(X[i]); }
    }
    try {
      const model = ols(Yt, Xt);
      for (let i = 0; i < Yv.length; i++) {
        const pred = [1, ...Xv[i]].reduce((s, v, j) => s + v * model.beta[j], 0);
        ss += (Yv[i] - pred) ** 2;
        count++;
      }
    } catch (e) { /* skip fold */ }
  }
  return count > 0 ? Math.sqrt(ss / count) : NaN;
}

const nIter = 200;
let coreSum = 0, vfrSum = 0, five5Sum = 0;
const core5_tr_cv = train.map((g, i) => [g.logMHI, g.logMhost, g.logMR, VfResid_tr[i], g.lhOuter]);
for (let iter = 0; iter < nIter; iter++) {
  coreSum += kfoldRMSE(Y_tr, core_tr, 5);
  vfrSum += kfoldRMSE(Y_tr, coreVfResid_tr, 5);
  five5Sum += kfoldRMSE(Y_tr, core5_tr_cv, 5);
}
const coreCV = coreSum / nIter;
const vfrCV = vfrSum / nIter;
const fiveCV = five5Sum / nIter;
console.log('  5-fold CV (avg over ' + nIter + ' shuffles, within published N=45):');
console.log('  ' + 'Model'.padEnd(30) + 'CV-RMSE   CV-Gap%');
console.log('  ' + '─'.repeat(48));
console.log('  ' + 'Core'.padEnd(30) + coreCV.toFixed(4) + '    ' + gapPctLOO(coreCV, sdY_p).toFixed(1) + '%');
console.log('  ' + 'Core+VfResid'.padEnd(30) + vfrCV.toFixed(4) + '    ' + gapPctLOO(vfrCV, sdY_p).toFixed(1) + '%');
console.log('  ' + 'Core+VfResid+lhOuter'.padEnd(30) + fiveCV.toFixed(4) + '    ' + gapPctLOO(fiveCV, sdY_p).toFixed(1) + '%');
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 134: VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const tgCore = transferGap(coreResult.rmse, naiveRMSE);
const tgVflat = transferGap(vflatResult.rmse, naiveRMSE);
const tgVfResid = transferGap(vfResidResult.rmse, naiveRMSE);
const tg5 = transferGap(fiveResult.rmse, naiveRMSE);
const hierarchy = tgVfResid > tgCore && tgVflat > tgCore;

let verdict;
if (hierarchy && allSignsSame && tgVfResid > 40 && rVfResidHoldout > 0.5) {
  verdict = 'STRONG_TRANSFER';
} else if (hierarchy && allSignsSame && tgVfResid > 20) {
  verdict = 'MODERATE_TRANSFER';
} else if (allSignsSame && tgCore > 0) {
  verdict = 'CORE_TRANSFERS';
} else {
  verdict = 'TRANSFER_FAILS';
}

console.log('  +------------------------------------------------------------------+');
console.log('  |  VERDICT: ' + verdict.padEnd(55) + '|');
console.log('  +------------------------------------------------------------------+\n');

console.log('  TRANSFER SCORECARD (baseline = predict training mean):');
console.log('    Naive RMSE (baseline): ' + naiveRMSE.toFixed(4));
console.log('    Core gap (holdout): ' + tgCore.toFixed(1) + '%');
console.log('    Model C gap (holdout): ' + tgVflat.toFixed(1) + '%');
console.log('    Core+VfResid gap (holdout): ' + tgVfResid.toFixed(1) + '%');
console.log('    5-axis gap (holdout): ' + tg5.toFixed(1) + '%');
console.log('    Signs preserved: ' + (allSignsSame ? 'YES' : 'NO'));
console.log('    Hierarchy preserved: ' + (hierarchy ? 'YES' : 'NO'));
console.log('    r(VfResid, a0) holdout: ' + rVfResidHoldout.toFixed(3));
console.log('');

console.log('  CAVEATS:');
console.log('    1. N_test=10 is small. Results directional, not conclusive.');
console.log('    2. All test galaxies Vflat>=120. Regime split not independently tested.');
console.log('    3. Crude quality — measurement errors larger than published sample.');
if (extrapolating.length > 0) {
  console.log('    4. ' + extrapolating.length + ' test galaxy(s) in extrapolation territory.');
}
console.log('');

const output = {
  phase: '134',
  title: 'External / Transfer Validation',
  verdict,
  trainN: N_train, testN: N_test,
  baselines: {
    naiveRMSE: +naiveRMSE.toFixed(4),
    testSD: +sdY_te.toFixed(4),
    trainMean: +meanY_tr.toFixed(4),
    testMean: +meanY_te.toFixed(4),
    meanDiff: +(meanY_te - meanY_tr).toFixed(4)
  },
  transfer: {
    core: { rmse: +coreResult.rmse.toFixed(4), gap: +tgCore.toFixed(1) },
    modelC: { rmse: +vflatResult.rmse.toFixed(4), gap: +tgVflat.toFixed(1) },
    vfResid: { rmse: +vfResidResult.rmse.toFixed(4), gap: +tgVfResid.toFixed(1) },
    fiveAxis: { rmse: +fiveResult.rmse.toFixed(4), gap: +tg5.toFixed(1) }
  },
  regimeTraining: {
    coreRMSE_full: +coreResult.rmse.toFixed(4),
    coreRMSE_regime: +hCoreR.rmse.toFixed(4),
    vfResidRMSE_full: +vfResidResult.rmse.toFixed(4),
    vfResidRMSE_regime: +hVfResidR.rmse.toFixed(4)
  },
  signsPreserved: allSignsSame,
  hierarchyPreserved: hierarchy,
  rVfResidHoldout: +rVfResidHoldout.toFixed(3),
  rHaloKHoldout: +rHaloKHoldout.toFixed(3),
  perGalaxy: test.map((g, i) => ({
    name: g.name, Vflat: g.Vflat,
    actual: +g.logA0.toFixed(3),
    predVfResid: +vfResidResult.preds[i].toFixed(3),
    errorVfResid: +vfResidResult.errs[i].toFixed(3),
    VfResid: +VfResid_te[i].toFixed(4),
    lhOuter: +g.lhOuter.toFixed(3)
  })),
  combinedLOO: {
    coreGap_pub: +gapPctLOO(looCoreP, sdY_p).toFixed(1),
    coreGap_all: +gapPctLOO(looCoreAll, sdY_all).toFixed(1),
    vfResidGap_pub: +gapPctLOO(looVfResidP, sdY_p).toFixed(1),
    vfResidGap_all: +gapPctLOO(looVfResidAll, sdY_all).toFixed(1),
    fiveAxisGap_all: +gapPctLOO(loo5All, sdY_all).toFixed(1)
  },
  crossValidation: {
    coreCV: +coreCV.toFixed(4),
    vfResidCV: +vfrCV.toFixed(4),
    fiveAxisCV: +fiveCV.toFixed(4)
  },
  caveats: [
    'N_test=10 is small',
    'All test galaxies Vflat>=120, regime split not independently tested',
    'Crude quality, larger measurement errors',
    extrapolating.length > 0 ? extrapolating.length + ' test galaxy(s) in extrapolation territory' : null
  ].filter(Boolean)
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase134-external-validation.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase134-external-validation.json');

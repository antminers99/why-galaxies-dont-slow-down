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
function gapPct(rms, sdy) { return sdy > 0 ? 100 * (1 - rms ** 2 / sdy ** 2) : 0; }

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
console.log('  Test Vflat >= 180: ' + test.filter(g => g.Vflat >= 180).length + '\n');

const Y_tr = train.map(g => g.logA0);
const Y_te = test.map(g => g.logA0);
const sdY_te = sd(Y_te);

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

function transferRMSE(Y_tr, X_tr, X_te, Y_te) {
  const model = ols(Y_tr, X_tr);
  const preds = X_te.map(x => [1, ...x].reduce((s, v, j) => s + v * model.beta[j], 0));
  const errs = Y_te.map((y, i) => y - preds[i]);
  const rmse = Math.sqrt(errs.reduce((s, e) => s + e * e, 0) / Y_te.length);
  return { rmse, preds, errs, model };
}

const baselineRMSE = sdY_te;
const coreResult = transferRMSE(Y_tr, core_tr, core_te, Y_te);
const vflatResult = transferRMSE(Y_tr, coreVflat_tr, coreVflat_te, Y_te);
const vfResidResult = transferRMSE(Y_tr, coreVfResid_tr, coreVfResid_te, Y_te);
const fiveResult = transferRMSE(Y_tr, core5_tr, core5_te, Y_te);

console.log('  Model'.padEnd(30) + 'RMSE    Gap%');
console.log('  ' + '─'.repeat(45));
console.log('  ' + 'Baseline (SD)'.padEnd(30) + baselineRMSE.toFixed(4) + '    0.0%');
console.log('  ' + 'Core (3-axis)'.padEnd(30) + coreResult.rmse.toFixed(4) + '   ' + gapPct(coreResult.rmse, sdY_te).toFixed(1) + '%');
console.log('  ' + 'Core + Vflat (Model C)'.padEnd(30) + vflatResult.rmse.toFixed(4) + '   ' + gapPct(vflatResult.rmse, sdY_te).toFixed(1) + '%');
console.log('  ' + 'Core + VfResid'.padEnd(30) + vfResidResult.rmse.toFixed(4) + '   ' + gapPct(vfResidResult.rmse, sdY_te).toFixed(1) + '%');
console.log('  ' + 'Core + VfResid + lhOuter (5-ax)'.padEnd(30) + fiveResult.rmse.toFixed(4) + '   ' + gapPct(fiveResult.rmse, sdY_te).toFixed(1) + '%');
console.log('');

console.log('  Per-galaxy predictions (Core + Vflat):');
console.log('  ' + 'Galaxy'.padEnd(15) + 'Actual  Pred    Error   Vflat');
console.log('  ' + '─'.repeat(55));
for (let i = 0; i < N_test; i++) {
  const err = vflatResult.errs[i];
  const flag = Math.abs(err) > 0.3 ? '⚠' : Math.abs(err) > 0.15 ? '●' : '✓';
  console.log('  ' + flag + ' ' + test[i].name.padEnd(14) +
    Y_te[i].toFixed(3).padStart(6) + '  ' + vflatResult.preds[i].toFixed(3).padStart(6) + '  ' +
    (err > 0 ? '+' : '') + err.toFixed(3).padStart(6) + '   ' + test[i].Vflat.toFixed(0));
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: REGIME-SPECIFIC TRANSFER');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const highTest = test.filter(g => g.Vflat >= 120);
const lowTest = test.filter(g => g.Vflat < 120);
const nHighTest = highTest.length;

console.log('  High-Vflat test set (Vflat>=120): N=' + nHighTest);
if (nHighTest >= 3) {
  const Y_hte = highTest.map(g => g.logA0);
  const sdY_hte = sd(Y_hte);

  const core_hte = highTest.map(g => [g.logMHI, g.logMhost, g.logMR]);
  const vflat_hte = highTest.map(g => [g.logMHI, g.logMhost, g.logMR, g.logVflat]);
  const vfr_hte = highTest.map(g => [g.logMHI, g.logMhost, g.logMR, computeVfResid(g)]);

  const highTrain = train.filter(g => g.Vflat >= 120);
  const Y_htr = highTrain.map(g => g.logA0);
  const core_htr = highTrain.map(g => [g.logMHI, g.logMhost, g.logMR]);
  const vflat_htr = highTrain.map(g => [g.logMHI, g.logMhost, g.logMR, g.logVflat]);
  const VfR_htr = highTrain.map(computeVfResid);
  const vfr_htr = highTrain.map((g, i) => [g.logMHI, g.logMhost, g.logMR, VfR_htr[i]]);

  const hCoreR = transferRMSE(Y_htr, core_htr, core_hte, Y_hte);
  const hVflatR = transferRMSE(Y_htr, vflat_htr, vflat_hte, Y_hte);
  const hVfResidR = transferRMSE(Y_htr, vfr_htr, vfr_hte, Y_hte);

  console.log('  Trained on high-Vflat published (N=' + highTrain.length + '), predict high-Vflat crude:\n');
  console.log('  ' + 'Model'.padEnd(25) + 'RMSE    Gap%');
  console.log('  ' + '─'.repeat(40));
  console.log('  ' + 'Baseline'.padEnd(25) + sdY_hte.toFixed(4) + '    0.0%');
  console.log('  ' + 'Core'.padEnd(25) + hCoreR.rmse.toFixed(4) + '   ' + gapPct(hCoreR.rmse, sdY_hte).toFixed(1) + '%');
  console.log('  ' + 'Core + Vflat'.padEnd(25) + hVflatR.rmse.toFixed(4) + '   ' + gapPct(hVflatR.rmse, sdY_hte).toFixed(1) + '%');
  console.log('  ' + 'Core + VfResid'.padEnd(25) + hVfResidR.rmse.toFixed(4) + '   ' + gapPct(hVfResidR.rmse, sdY_hte).toFixed(1) + '%');

  const fullTrainCoreR = transferRMSE(Y_tr, core_tr, core_hte, Y_hte);
  const fullTrainVflatR = transferRMSE(Y_tr, coreVflat_tr, vflat_hte, Y_hte);
  const fullTrainVfResidR = transferRMSE(Y_tr, coreVfResid_tr, vfr_hte, Y_hte);

  console.log('\n  Trained on ALL published (N=45), predict high-Vflat crude:\n');
  console.log('  ' + 'Model'.padEnd(25) + 'RMSE    Gap%');
  console.log('  ' + '─'.repeat(40));
  console.log('  ' + 'Baseline'.padEnd(25) + sdY_hte.toFixed(4) + '    0.0%');
  console.log('  ' + 'Core'.padEnd(25) + fullTrainCoreR.rmse.toFixed(4) + '   ' + gapPct(fullTrainCoreR.rmse, sdY_hte).toFixed(1) + '%');
  console.log('  ' + 'Core + Vflat'.padEnd(25) + fullTrainVflatR.rmse.toFixed(4) + '   ' + gapPct(fullTrainVflatR.rmse, sdY_hte).toFixed(1) + '%');
  console.log('  ' + 'Core + VfResid'.padEnd(25) + fullTrainVfResidR.rmse.toFixed(4) + '   ' + gapPct(fullTrainVfResidR.rmse, sdY_hte).toFixed(1) + '%');
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: SIGN & HIERARCHY TRANSFER');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const trainCoreModel = ols(Y_tr, core_tr);
const trainVflatModel = ols(Y_tr, coreVflat_tr);
const trainVfResidModel = ols(Y_tr, coreVfResid_tr);

console.log('  Trained coefficients (published):');
console.log('    Core: intercept=' + trainCoreModel.beta[0].toFixed(3) +
  ' MHI=' + trainCoreModel.beta[1].toFixed(3) +
  ' Mhost=' + trainCoreModel.beta[2].toFixed(3) +
  ' MR=' + trainCoreModel.beta[3].toFixed(3));
console.log('    ModelC: intercept=' + trainVflatModel.beta[0].toFixed(3) +
  ' MHI=' + trainVflatModel.beta[1].toFixed(3) +
  ' Mhost=' + trainVflatModel.beta[2].toFixed(3) +
  ' MR=' + trainVflatModel.beta[3].toFixed(3) +
  ' Vflat=' + trainVflatModel.beta[4].toFixed(3));
console.log('');

const allGalsData = [...train, ...test];
const Y_all = allGalsData.map(g => g.logA0);
const core_all = allGalsData.map(g => [g.logMHI, g.logMhost, g.logMR]);
const vflat_all = allGalsData.map(g => [g.logMHI, g.logMhost, g.logMR, g.logVflat]);
const VfR_all = allGalsData.map(computeVfResid);
const vfr_all = allGalsData.map((g, i) => [g.logMHI, g.logMhost, g.logMR, VfR_all[i]]);

const fullCoreModel = ols(Y_all, core_all);
const fullVflatModel = ols(Y_all, vflat_all);
const fullVfResidModel = ols(Y_all, vfr_all);

console.log('  Combined (N=55) coefficients:');
console.log('    Core: intercept=' + fullCoreModel.beta[0].toFixed(3) +
  ' MHI=' + fullCoreModel.beta[1].toFixed(3) +
  ' Mhost=' + fullCoreModel.beta[2].toFixed(3) +
  ' MR=' + fullCoreModel.beta[3].toFixed(3));
console.log('    ModelC: intercept=' + fullVflatModel.beta[0].toFixed(3) +
  ' MHI=' + fullVflatModel.beta[1].toFixed(3) +
  ' Mhost=' + fullVflatModel.beta[2].toFixed(3) +
  ' MR=' + fullVflatModel.beta[3].toFixed(3) +
  ' Vflat=' + fullVflatModel.beta[4].toFixed(3));
console.log('');

const signsSame = trainVflatModel.beta.slice(1).every((b, j) =>
  Math.sign(b) === Math.sign(fullVflatModel.beta[j + 1]));
console.log('  Sign stability (Model C published vs combined): ' + (signsSame ? 'ALL SAME ✓' : 'CHANGED ✗'));
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
    } catch (e) { ss += (Y[i] - mean(Y)) ** 2; }
  }
  return Math.sqrt(ss / n);
}

const sdY_all = sd(Y_all);
const core5_all = allGalsData.map((g, i) => [g.logMHI, g.logMhost, g.logMR, VfR_all[i], g.lhOuter]);

const looCoreP = looRMS(Y_tr, core_tr);
const looVflatP = looRMS(Y_tr, coreVflat_tr);
const looVfResidP = looRMS(Y_tr, coreVfResid_tr);
const sdY_p = sd(Y_tr);

const looCoreAll = looRMS(Y_all, core_all);
const looVflatAll = looRMS(Y_all, vflat_all);
const looVfResidAll = looRMS(Y_all, vfr_all);
const loo5All = looRMS(Y_all, core5_all);

console.log('  ' + 'Model'.padEnd(30) + '  Pub(N=45)  Combined(N=55)');
console.log('  ' + '─'.repeat(55));
console.log('  ' + 'Core gap'.padEnd(30) + gapPct(looCoreP, sdY_p).toFixed(1).padStart(8) + '%   ' + gapPct(looCoreAll, sdY_all).toFixed(1).padStart(8) + '%');
console.log('  ' + 'Core+Vflat gap'.padEnd(30) + gapPct(looVflatP, sdY_p).toFixed(1).padStart(8) + '%   ' + gapPct(looVflatAll, sdY_all).toFixed(1).padStart(8) + '%');
console.log('  ' + 'Core+VfResid gap'.padEnd(30) + gapPct(looVfResidP, sdY_p).toFixed(1).padStart(8) + '%   ' + gapPct(looVfResidAll, sdY_all).toFixed(1).padStart(8) + '%');
console.log('  ' + 'Core+VfResid+lhOuter gap'.padEnd(30) + '     N/A' + '   ' + gapPct(loo5All, sdY_all).toFixed(1).padStart(8) + '%');
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
console.log('  r(VfResid, logA0) on holdout: ' + pearsonR(VfResid_te, Y_te).toFixed(3));
console.log('  r(haloK, logA0) on holdout: ' + pearsonR(test.map(g => g.haloK), Y_te).toFixed(3));
console.log('  r(VfResid, haloK) on holdout: ' + pearsonR(VfResid_te, test.map(g => g.haloK)).toFixed(3));
console.log('');

if (nHighTest >= 5) {
  const Y_hte = highTest.map(g => g.logA0);
  const VfR_hte = highTest.map(computeVfResid);
  console.log('  High-Vflat holdout only (N=' + nHighTest + '):');
  console.log('    r(VfResid, logA0): ' + pearsonR(VfR_hte, Y_hte).toFixed(3));
  console.log('    r(haloK, logA0): ' + pearsonR(highTest.map(g => g.haloK), Y_hte).toFixed(3));
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 6: lh_outerImprove ON HOLDOUT — Does the 5th axis transfer?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const lhOuter_te = test.map(g => g.lhOuter);
console.log('  r(lhOuter, logA0) on holdout: ' + pearsonR(lhOuter_te, Y_te).toFixed(3));
console.log('  Transfer:');
console.log('    Core+VfResid RMSE (holdout): ' + vfResidResult.rmse.toFixed(4));
console.log('    Core+VfResid+lhOuter RMSE: ' + fiveResult.rmse.toFixed(4));
console.log('    5-axis improves: ' + (fiveResult.rmse < vfResidResult.rmse ? 'YES' : 'NO') +
  ' (delta=' + (gapPct(fiveResult.rmse, sdY_te) - gapPct(vfResidResult.rmse, sdY_te)).toFixed(1) + 'pp)');
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 134: VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const transferGapVflat = gapPct(vflatResult.rmse, sdY_te);
const transferGapVfResid = gapPct(vfResidResult.rmse, sdY_te);
const transferGapCore = gapPct(coreResult.rmse, sdY_te);
const rVfResidHoldout = pearsonR(VfResid_te, Y_te);
const hierarchy = transferGapVflat > transferGapCore && transferGapVfResid > transferGapCore;
const signsOK = signsSame;

let verdict;
if (hierarchy && signsOK && transferGapVflat > 20 && rVfResidHoldout > 0.3) {
  verdict = 'STRONG_TRANSFER';
} else if (hierarchy && signsOK && transferGapVflat > 0) {
  verdict = 'MODERATE_TRANSFER';
} else if (signsOK && transferGapCore > 0) {
  verdict = 'CORE_TRANSFERS';
} else {
  verdict = 'TRANSFER_FAILS';
}

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + verdict.padEnd(55) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════════╝\n');

console.log('  TRANSFER SCORECARD:');
console.log('    Core gap (holdout): ' + transferGapCore.toFixed(1) + '%');
console.log('    Model C gap (holdout): ' + transferGapVflat.toFixed(1) + '%');
console.log('    Core+VfResid gap (holdout): ' + transferGapVfResid.toFixed(1) + '%');
console.log('    Signs preserved: ' + (signsOK ? 'YES' : 'NO'));
console.log('    Hierarchy preserved: ' + (hierarchy ? 'YES' : 'NO'));
console.log('    r(VfResid, a₀) holdout: ' + rVfResidHoldout.toFixed(3));
console.log('');

const output = {
  phase: '134',
  title: 'External / Transfer Validation',
  verdict,
  trainN: N_train, testN: N_test,
  transfer: {
    baseline: +sdY_te.toFixed(4),
    core: { rmse: +coreResult.rmse.toFixed(4), gap: +transferGapCore.toFixed(1) },
    modelC: { rmse: +vflatResult.rmse.toFixed(4), gap: +transferGapVflat.toFixed(1) },
    vfResid: { rmse: +vfResidResult.rmse.toFixed(4), gap: +transferGapVfResid.toFixed(1) },
    fiveAxis: { rmse: +fiveResult.rmse.toFixed(4), gap: +gapPct(fiveResult.rmse, sdY_te).toFixed(1) }
  },
  signsPreserved: signsOK,
  hierarchyPreserved: hierarchy,
  rVfResidHoldout: +rVfResidHoldout.toFixed(3),
  perGalaxy: test.map((g, i) => ({
    name: g.name, Vflat: g.Vflat,
    actual: +g.logA0.toFixed(3),
    predModelC: +vflatResult.preds[i].toFixed(3),
    error: +vflatResult.errs[i].toFixed(3),
    VfResid: +VfResid_te[i].toFixed(4)
  })),
  combinedLOO: {
    coreGap_pub: +gapPct(looCoreP, sdY_p).toFixed(1),
    coreGap_all: +gapPct(looCoreAll, sdY_all).toFixed(1),
    vflatGap_pub: +gapPct(looVflatP, sdY_p).toFixed(1),
    vflatGap_all: +gapPct(looVflatAll, sdY_all).toFixed(1),
    vfResidGap_pub: +gapPct(looVfResidP, sdY_p).toFixed(1),
    vfResidGap_all: +gapPct(looVfResidAll, sdY_all).toFixed(1)
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase134-external-validation.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase134-external-validation.json');

#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const pub = p => path.join(__dirname, '..', 'public', p);

console.log('======================================================================');
console.log('  PHASE 201: BLIND PREDICTION — FROZEN COEFFICIENTS');
console.log('  Apply Track 2 models trained on N=45 to N=59 external galaxies');
console.log('  NO refitting. Pure out-of-sample prediction.');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(pub('stage-A-master-table.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(pub('phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const sparcTable = JSON.parse(fs.readFileSync(pub('sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(pub('sparc-results.json'), 'utf8'));
const extDataset = JSON.parse(fs.readFileSync(pub('phase200-external-dataset.json'), 'utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparcTable.forEach(s => { sparcMap[s.name] = s; });
const resMap = {}; sparcResults.perGalaxy.forEach(g => { if (g.rawName) resMap[g.rawName] = g; });

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : 0; }
function sd(a) { const m = mean(a); return a.length > 1 ? Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)) : 0; }
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
  return { beta, resid };
}

function transferGap(rmse, baselineRmse) {
  return baselineRmse > 0 ? 100 * (1 - rmse ** 2 / baselineRmse ** 2) : 0;
}

const pubGals = stageA.galaxies.filter(g => pubNames.has(g.name));

const structModel = ols(
  pubGals.map(g => Math.log10(sparcMap[g.name].Vflat)),
  pubGals.map(g => {
    const s = sparcMap[g.name];
    return [
      Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9),
      Math.log10(Math.max(s.L36, 0.01)),
      Math.log10(s.Rdisk),
      s.T
    ];
  })
);

function getTrainGal(g) {
  const s = sparcMap[g.name]; const r = resMap[g.name]; const td = tdMap[g.name];
  const logMbar = Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
  const logVflat = Math.log10(s.Vflat);
  const predLV = [1, logMbar, Math.log10(Math.max(s.L36, 0.01)), Math.log10(s.Rdisk), s.T]
    .reduce((sum, x, j) => sum + x * structModel.beta[j], 0);
  return {
    name: g.name, logA0: g.logA0,
    logMHI: g.logMHI, logMhost: td.logMhost, logMR: g.logMeanRun,
    VfResid: logVflat - predLV,
    lhOuter: r.models.log_halo.outerImprovement,
    Vflat: s.Vflat
  };
}

const train = pubGals.map(getTrainGal);
const Y_tr = train.map(g => g.logA0);
const meanY_tr = mean(Y_tr);

const ext = extDataset.galaxies;

const naiveRMSE = Math.sqrt(ext.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / ext.length);

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  BASELINES');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
console.log('  Training N=' + train.length + '  External N=' + ext.length);
console.log('  Training mean logA0: ' + meanY_tr.toFixed(4));
console.log('  External mean logA0: ' + mean(ext.map(g => g.logA0)).toFixed(4));
console.log('  Naive RMSE (predict training mean): ' + naiveRMSE.toFixed(4));
console.log('  External SD logA0: ' + sd(ext.map(g => g.logA0)).toFixed(4));

const coreModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR]));
const vfResidModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR, g.VfResid]));
const fiveAxisModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR, g.VfResid, g.lhOuter]));

console.log('\n  Frozen coefficients (trained on N=45):');
console.log('  Core: [' + coreModel.beta.map(b => b.toFixed(4)).join(', ') + ']');
console.log('  Core+VfResid: [' + vfResidModel.beta.map(b => b.toFixed(4)).join(', ') + ']');
console.log('  5-axis: [' + fiveAxisModel.beta.map(b => b.toFixed(4)).join(', ') + ']');

function predictModel(model, features) {
  return [1, ...features].reduce((s, x, j) => s + x * model.beta[j], 0);
}

function runTransfer(name, model, featureFn, galaxies, baseline) {
  const preds = galaxies.map(g => predictModel(model, featureFn(g)));
  const errs = galaxies.map((g, i) => g.logA0 - preds[i]);
  const rmse = Math.sqrt(errs.reduce((s, e) => s + e * e, 0) / errs.length);
  const gap = transferGap(rmse, baseline);
  const r = pearsonR(galaxies.map(g => g.logA0), preds);
  return { name, rmse, gap, r, preds, errs };
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: FULL SAMPLE TRANSFER (N=' + ext.length + ')');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const coreResult = runTransfer('Core', coreModel, g => [g.logMHI, g.logMhost, g.logMeanRun], ext, naiveRMSE);
const vfrResult = runTransfer('Core+VfResid', vfResidModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid], ext, naiveRMSE);

const fiveFeatures = g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid, g.lhOuter];
const ext5 = ext.filter(g => g.lhOuter !== null);
const naiveRMSE5 = Math.sqrt(ext5.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / ext5.length);
const fiveResult = runTransfer('5-axis', fiveAxisModel, fiveFeatures, ext5, naiveRMSE5);

console.log('  Model'.padEnd(25) + 'RMSE'.padStart(8) + 'Gap%'.padStart(8) + 'r(pred,act)'.padStart(12));
console.log('  ' + '─'.repeat(55));
console.log('  ' + 'Naive (training mean)'.padEnd(25) + naiveRMSE.toFixed(4).padStart(8));
console.log('  ' + coreResult.name.padEnd(25) + coreResult.rmse.toFixed(4).padStart(8) + (coreResult.gap.toFixed(1) + '%').padStart(8) + coreResult.r.toFixed(3).padStart(12));
console.log('  ' + vfrResult.name.padEnd(25) + vfrResult.rmse.toFixed(4).padStart(8) + (vfrResult.gap.toFixed(1) + '%').padStart(8) + vfrResult.r.toFixed(3).padStart(12));
console.log('  ' + fiveResult.name.padEnd(25) + fiveResult.rmse.toFixed(4).padStart(8) + (fiveResult.gap.toFixed(1) + '%').padStart(8) + fiveResult.r.toFixed(3).padStart(12));

console.log('\n  Phase 134 reference (N=10 holdout):');
console.log('    Core gap: -11.5%  |  Core+VfResid gap: 56.9%  |  5-axis gap: 66.0%');

const hierPreserved = vfrResult.gap > coreResult.gap && vfrResult.rmse < coreResult.rmse;
console.log('\n  Hierarchy preserved (VfResid > Core): ' + (hierPreserved ? 'YES' : 'NO'));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: HIGH-VFLAT REGIME (Vflat >= 120)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const highV = ext.filter(g => g.Vflat >= 120);
const naiveHV = Math.sqrt(highV.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / highV.length);
const coreHV = runTransfer('Core', coreModel, g => [g.logMHI, g.logMhost, g.logMeanRun], highV, naiveHV);
const vfrHV = runTransfer('Core+VfResid', vfResidModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid], highV, naiveHV);
const hv5 = highV.filter(g => g.lhOuter !== null);
const naiveHV5 = Math.sqrt(hv5.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / hv5.length);
const fiveHV = runTransfer('5-axis', fiveAxisModel, fiveFeatures, hv5, naiveHV5);

console.log('  N=' + highV.length + ' galaxies with Vflat >= 120');
console.log('  Model'.padEnd(25) + 'RMSE'.padStart(8) + 'Gap%'.padStart(8) + 'r(pred,act)'.padStart(12));
console.log('  ' + '─'.repeat(55));
console.log('  ' + 'Naive'.padEnd(25) + naiveHV.toFixed(4).padStart(8));
console.log('  ' + coreHV.name.padEnd(25) + coreHV.rmse.toFixed(4).padStart(8) + (coreHV.gap.toFixed(1) + '%').padStart(8) + coreHV.r.toFixed(3).padStart(12));
console.log('  ' + vfrHV.name.padEnd(25) + vfrHV.rmse.toFixed(4).padStart(8) + (vfrHV.gap.toFixed(1) + '%').padStart(8) + vfrHV.r.toFixed(3).padStart(12));
console.log('  ' + fiveHV.name.padEnd(25) + fiveHV.rmse.toFixed(4).padStart(8) + (fiveHV.gap.toFixed(1) + '%').padStart(8) + fiveHV.r.toFixed(3).padStart(12));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: VERY HIGH-VFLAT (Vflat >= 180)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const vhV = ext.filter(g => g.Vflat >= 180);
if (vhV.length >= 4) {
  const naiveVH = Math.sqrt(vhV.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / vhV.length);
  const coreVH = runTransfer('Core', coreModel, g => [g.logMHI, g.logMhost, g.logMeanRun], vhV, naiveVH);
  const vfrVH = runTransfer('Core+VfResid', vfResidModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid], vhV, naiveVH);

  console.log('  N=' + vhV.length + ' galaxies with Vflat >= 180');
  console.log('  Model'.padEnd(25) + 'RMSE'.padStart(8) + 'Gap%'.padStart(8) + 'r(pred,act)'.padStart(12));
  console.log('  ' + '─'.repeat(55));
  console.log('  ' + 'Naive'.padEnd(25) + naiveVH.toFixed(4).padStart(8));
  console.log('  ' + coreVH.name.padEnd(25) + coreVH.rmse.toFixed(4).padStart(8) + (coreVH.gap.toFixed(1) + '%').padStart(8) + coreVH.r.toFixed(3).padStart(12));
  console.log('  ' + vfrVH.name.padEnd(25) + vfrVH.rmse.toFixed(4).padStart(8) + (vfrVH.gap.toFixed(1) + '%').padStart(8) + vfrVH.r.toFixed(3).padStart(12));
} else {
  console.log('  Insufficient galaxies (N=' + vhV.length + ')');
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: LOW-VFLAT REGIME (Vflat < 120) — EXPECTED TO FAIL');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const lowV = ext.filter(g => g.Vflat < 120);
const naiveLV = Math.sqrt(lowV.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / lowV.length);
const coreLV = runTransfer('Core', coreModel, g => [g.logMHI, g.logMhost, g.logMeanRun], lowV, naiveLV);
const vfrLV = runTransfer('Core+VfResid', vfResidModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid], lowV, naiveLV);

console.log('  N=' + lowV.length + ' galaxies with Vflat < 120');
console.log('  Model'.padEnd(25) + 'RMSE'.padStart(8) + 'Gap%'.padStart(8) + 'r(pred,act)'.padStart(12));
console.log('  ' + '─'.repeat(55));
console.log('  ' + 'Naive'.padEnd(25) + naiveLV.toFixed(4).padStart(8));
console.log('  ' + coreLV.name.padEnd(25) + coreLV.rmse.toFixed(4).padStart(8) + (coreLV.gap.toFixed(1) + '%').padStart(8) + coreLV.r.toFixed(3).padStart(12));
console.log('  ' + vfrLV.name.padEnd(25) + vfrLV.rmse.toFixed(4).padStart(8) + (vfrLV.gap.toFixed(1) + '%').padStart(8) + vfrLV.r.toFixed(3).padStart(12));
console.log('  (Negative gap = model worse than naive. Expected for low-Vflat regime.)');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 5: PER-GALAXY PREDICTIONS (Core+VfResid)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  ' + 'Galaxy'.padEnd(16) + 'Vflat'.padStart(5) + '  actual  pred   error  VfRes  nPts Q');
console.log('  ' + '─'.repeat(68));
ext.sort((a, b) => b.Vflat - a.Vflat);
ext.forEach((g, i) => {
  const pred = predictModel(vfResidModel, [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid]);
  const err = g.logA0 - pred;
  console.log('  ' + g.name.padEnd(16) +
    g.Vflat.toFixed(0).padStart(5) + '  ' +
    g.logA0.toFixed(3).padStart(6) + ' ' +
    pred.toFixed(3).padStart(6) + '  ' +
    (err >= 0 ? '+' : '') + err.toFixed(3) + '  ' +
    g.VfResid.toFixed(3).padStart(6) + '  ' +
    String(g.nRARpts).padStart(3) + ' ' + g.Q);
});

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 6: SIGN PRESERVATION');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const allData = [...train.map(g => ({
  logA0: g.logA0, logMHI: g.logMHI, logMhost: g.logMhost,
  logMeanRun: g.logMR, VfResid: g.VfResid, lhOuter: g.lhOuter
})), ...ext.filter(g => g.lhOuter !== null)];

const allY = allData.map(g => g.logA0);
const allVfResidModel = ols(allY, allData.map(g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid]));

const labels = ['intercept', 'logMHI', 'logMhost', 'logMR', 'VfResid'];
console.log('  ' + 'Coeff'.padEnd(12) + 'Train(N=45)'.padStart(12) + 'Combined(N=' + allData.length + ')'.padStart(14) + '  Sign');
console.log('  ' + '─'.repeat(45));
let allSigns = true;
for (let j = 0; j < labels.length; j++) {
  const a = vfResidModel.beta[j], b = allVfResidModel.beta[j];
  const same = Math.sign(a) === Math.sign(b);
  if (!same) allSigns = false;
  console.log('  ' + labels[j].padEnd(12) + a.toFixed(4).padStart(12) + b.toFixed(4).padStart(14) + '  ' + (same ? 'SAME' : 'CHANGED'));
}
console.log('\n  All signs preserved: ' + (allSigns ? 'YES' : 'NO'));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 7: EXCLUDING EXTRAPOLATING GALAXIES');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const noExtrap = ext.filter(g => !g.extrapolating);
const naiveNE = Math.sqrt(noExtrap.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / noExtrap.length);
const coreNE = runTransfer('Core', coreModel, g => [g.logMHI, g.logMhost, g.logMeanRun], noExtrap, naiveNE);
const vfrNE = runTransfer('Core+VfResid', vfResidModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid], noExtrap, naiveNE);

console.log('  N=' + noExtrap.length + ' galaxies (excluding ' + (ext.length - noExtrap.length) + ' extrapolating)');
console.log('  Model'.padEnd(25) + 'RMSE'.padStart(8) + 'Gap%'.padStart(8) + 'r(pred,act)'.padStart(12));
console.log('  ' + '─'.repeat(55));
console.log('  ' + 'Naive'.padEnd(25) + naiveNE.toFixed(4).padStart(8));
console.log('  ' + coreNE.name.padEnd(25) + coreNE.rmse.toFixed(4).padStart(8) + (coreNE.gap.toFixed(1) + '%').padStart(8) + coreNE.r.toFixed(3).padStart(12));
console.log('  ' + vfrNE.name.padEnd(25) + vfrNE.rmse.toFixed(4).padStart(8) + (vfrNE.gap.toFixed(1) + '%').padStart(8) + vfrNE.r.toFixed(3).padStart(12));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 8: QUALITY FILTER (Q <= 1 only)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const q1Only = ext.filter(g => g.Q <= 1);
const naiveQ1 = Math.sqrt(q1Only.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / q1Only.length);
const vfrQ1 = runTransfer('Core+VfResid', vfResidModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid], q1Only, naiveQ1);

console.log('  N=' + q1Only.length + ' galaxies with Q=1');
console.log('  Core+VfResid: RMSE=' + vfrQ1.rmse.toFixed(4) + ' Gap=' + vfrQ1.gap.toFixed(1) + '% r=' + vfrQ1.r.toFixed(3));
const hvQ1 = q1Only.filter(g => g.Vflat >= 120);
if (hvQ1.length >= 5) {
  const naiveHVQ1 = Math.sqrt(hvQ1.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / hvQ1.length);
  const vfrHVQ1 = runTransfer('Core+VfResid', vfResidModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid], hvQ1, naiveHVQ1);
  console.log('  Q=1 + Vflat>=120 (N=' + hvQ1.length + '): RMSE=' + vfrHVQ1.rmse.toFixed(4) + ' Gap=' + vfrHVQ1.gap.toFixed(1) + '% r=' + vfrHVQ1.r.toFixed(3));
}

console.log('\n======================================================================');
console.log('  VERDICT');
console.log('======================================================================\n');

const fullGap = vfrResult.gap;
const hvGap = vfrHV.gap;
const hierOK = hierPreserved;
const signOK = allSigns;

let verdict;
if (hvGap > 40 && hierOK && signOK) verdict = 'STRONG_TRANSFER';
else if (hvGap > 20 && hierOK) verdict = 'MODERATE_TRANSFER';
else if (fullGap > 0 && hierOK) verdict = 'WEAK_TRANSFER';
else verdict = 'NO_TRANSFER';

console.log('  Verdict: ' + verdict);
console.log('  Full sample gap: ' + fullGap.toFixed(1) + '%');
console.log('  High-Vflat gap: ' + hvGap.toFixed(1) + '%');
console.log('  Hierarchy preserved: ' + hierOK);
console.log('  Signs preserved: ' + signOK);
console.log('  r(VfResid, logA0) full: ' + pearsonR(ext.map(g => g.VfResid), ext.map(g => g.logA0)).toFixed(3));
console.log('  r(VfResid, logA0) high-V: ' + pearsonR(highV.map(g => g.VfResid), highV.map(g => g.logA0)).toFixed(3));

const output = {
  phase: '201',
  title: 'Blind Prediction — Frozen Coefficients',
  verdict,
  fullSample: {
    N: ext.length,
    naiveRMSE: +naiveRMSE.toFixed(4),
    core: { rmse: +coreResult.rmse.toFixed(4), gap: +coreResult.gap.toFixed(1), r: +coreResult.r.toFixed(3) },
    coreVfResid: { rmse: +vfrResult.rmse.toFixed(4), gap: +vfrResult.gap.toFixed(1), r: +vfrResult.r.toFixed(3) },
    fiveAxis: { rmse: +fiveResult.rmse.toFixed(4), gap: +fiveResult.gap.toFixed(1), r: +fiveResult.r.toFixed(3) }
  },
  highVflat: {
    N: highV.length,
    naiveRMSE: +naiveHV.toFixed(4),
    core: { rmse: +coreHV.rmse.toFixed(4), gap: +coreHV.gap.toFixed(1), r: +coreHV.r.toFixed(3) },
    coreVfResid: { rmse: +vfrHV.rmse.toFixed(4), gap: +vfrHV.gap.toFixed(1), r: +vfrHV.r.toFixed(3) },
    fiveAxis: { rmse: +fiveHV.rmse.toFixed(4), gap: +fiveHV.gap.toFixed(1), r: +fiveHV.r.toFixed(3) }
  },
  lowVflat: {
    N: lowV.length,
    core: { rmse: +coreLV.rmse.toFixed(4), gap: +coreLV.gap.toFixed(1), r: +coreLV.r.toFixed(3) },
    coreVfResid: { rmse: +vfrLV.rmse.toFixed(4), gap: +vfrLV.gap.toFixed(1), r: +vfrLV.r.toFixed(3) }
  },
  hierarchyPreserved: hierOK,
  signsPreserved: signOK,
  rVfResidFull: +pearsonR(ext.map(g => g.VfResid), ext.map(g => g.logA0)).toFixed(3),
  rVfResidHighV: +pearsonR(highV.map(g => g.VfResid), highV.map(g => g.logA0)).toFixed(3),
  phase134Reference: {
    coreGap: -11.5, vfResidGap: 56.9, fiveAxisGap: 66.0, rVfResid: 0.801, N: 10
  }
};

fs.writeFileSync(pub('phase201-blind-prediction.json'), JSON.stringify(output, null, 2));
console.log('\n  Output: public/phase201-blind-prediction.json');

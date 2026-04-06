#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const pub = p => path.join(__dirname, '..', 'public', p);

console.log('======================================================================');
console.log('  PHASE 203: logMhost IMPROVEMENT FOR EXTERNAL SAMPLE');
console.log('  Does better logMhost restore Core weight externally?');
console.log('  Or does VfResid remain the dominant transfer channel?');
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
  const rmse = Math.sqrt(resid.reduce((s, e) => s + e * e, 0) / n);
  return { beta, resid, rmse };
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
    Vflat: s.Vflat, logVflat,
    logMbar, logRdisk: Math.log10(s.Rdisk), morphT: s.T
  };
}

const train = pubGals.map(getTrainGal);
const Y_tr = train.map(g => g.logA0);
const meanY_tr = mean(Y_tr);

const ext = extDataset.galaxies;
const ext5 = ext.filter(g => g.lhOuter !== null);

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 1: BASELINE — Current logMhost formula');
console.log('  Formula: 11.5 + 3.0 * (log10(Vflat) - 2.2)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function formulaA(logVflat) {
  return 11.5 + 3.0 * (logVflat - 2.2);
}

const trainRealMhost = train.map(g => g.logMhost);
const trainFormulaA = train.map(g => formulaA(g.logVflat));
const rFormulaA = pearsonR(trainRealMhost, trainFormulaA);
const rmseFormulaA = Math.sqrt(train.reduce((s, g) => {
  const pred = formulaA(g.logVflat);
  return s + (g.logMhost - pred) ** 2;
}, 0) / train.length);

console.log('  Model A (current formula) vs real logMhost on N=' + train.length + ':');
console.log('    r(formula, real):   ' + rFormulaA.toFixed(3));
console.log('    RMSE:               ' + rmseFormulaA.toFixed(4));
console.log('    Real logMhost range: [' + Math.min(...trainRealMhost).toFixed(2) + ', ' + Math.max(...trainRealMhost).toFixed(2) + ']');
console.log('    Formula range:       [' + Math.min(...trainFormulaA).toFixed(2) + ', ' + Math.max(...trainFormulaA).toFixed(2) + ']');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 2: TRAIN logMhost ESTIMATORS FROM N=45');
console.log('  Target = real logMhost from tidal analysis');
console.log('  Features: ONLY galaxy observables (no a0 leakage)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const mhostReal = train.map(g => g.logMhost);

const estB = ols(mhostReal, train.map(g => [g.logVflat]));
console.log('  Model B: logMhost ~ logVflat (trained)');
console.log('    beta = [' + estB.beta.map(b => b.toFixed(4)).join(', ') + ']');
console.log('    R² = ' + (1 - estB.rmse ** 2 / sd(mhostReal) ** 2 * (train.length - 1) / (train.length - 2)).toFixed(3));
console.log('    RMSE = ' + estB.rmse.toFixed(4));
const rB = pearsonR(mhostReal, train.map(g => [1, g.logVflat].reduce((s, x, j) => s + x * estB.beta[j], 0)));
console.log('    r = ' + rB.toFixed(3));
console.log('');

const estC = ols(mhostReal, train.map(g => [g.logVflat, g.logMbar, g.logRdisk, g.morphT]));
console.log('  Model C: logMhost ~ logVflat + logMbar + logRdisk + T');
console.log('    beta = [' + estC.beta.map(b => b.toFixed(4)).join(', ') + ']');
console.log('    R² = ' + (1 - estC.rmse ** 2 / sd(mhostReal) ** 2 * (train.length - 1) / (train.length - 5)).toFixed(3));
console.log('    RMSE = ' + estC.rmse.toFixed(4));
const rC = pearsonR(mhostReal, train.map(g => [1, g.logVflat, g.logMbar, g.logRdisk, g.morphT].reduce((s, x, j) => s + x * estC.beta[j], 0)));
console.log('    r = ' + rC.toFixed(3));
console.log('');

const estD1 = ols(mhostReal, train.map(g => [g.logVflat, g.logMbar]));
console.log('  Model D1: logMhost ~ logVflat + logMbar');
console.log('    beta = [' + estD1.beta.map(b => b.toFixed(4)).join(', ') + ']');
console.log('    RMSE = ' + estD1.rmse.toFixed(4));
const rD1 = pearsonR(mhostReal, train.map(g => [1, g.logVflat, g.logMbar].reduce((s, x, j) => s + x * estD1.beta[j], 0)));
console.log('    r = ' + rD1.toFixed(3));
console.log('');

const estD2 = ols(mhostReal, train.map(g => [g.logVflat, g.morphT]));
console.log('  Model D2: logMhost ~ logVflat + T');
console.log('    beta = [' + estD2.beta.map(b => b.toFixed(4)).join(', ') + ']');
console.log('    RMSE = ' + estD2.rmse.toFixed(4));
const rD2 = pearsonR(mhostReal, train.map(g => [1, g.logVflat, g.morphT].reduce((s, x, j) => s + x * estD2.beta[j], 0)));
console.log('    r = ' + rD2.toFixed(3));
console.log('');

console.log('  LOO cross-validation on N=45 (to check overfitting):');
function looRMSE(features) {
  let ss = 0;
  for (let i = 0; i < train.length; i++) {
    const trLoo = train.filter((_, j) => j !== i);
    const yLoo = trLoo.map(g => g.logMhost);
    const xLoo = trLoo.map(features);
    const mLoo = ols(yLoo, xLoo);
    const xi = features(train[i]);
    const pred = [1, ...xi].reduce((s, x, j) => s + x * mLoo.beta[j], 0);
    ss += (train[i].logMhost - pred) ** 2;
  }
  return Math.sqrt(ss / train.length);
}
const looBrmse = looRMSE(g => [g.logVflat]);
const looCrmse = looRMSE(g => [g.logVflat, g.logMbar, g.logRdisk, g.morphT]);
const looD1rmse = looRMSE(g => [g.logVflat, g.logMbar]);
const looD2rmse = looRMSE(g => [g.logVflat, g.morphT]);
console.log('    Model B LOO RMSE:  ' + looBrmse.toFixed(4));
console.log('    Model C LOO RMSE:  ' + looCrmse.toFixed(4));
console.log('    Model D1 LOO RMSE: ' + looD1rmse.toFixed(4));
console.log('    Model D2 LOO RMSE: ' + looD2rmse.toFixed(4));
console.log('    Formula A RMSE:    ' + rmseFormulaA.toFixed(4) + ' (no training, fixed formula)');
console.log('');

const estimators = {
  'A_formula': { name: 'Model A (current formula)', fn: g => formulaA(g.logVflat) },
  'B_logVflat': { name: 'Model B (logVflat trained)', fn: g => [1, g.logVflat].reduce((s, x, j) => s + x * estB.beta[j], 0) },
  'C_multivar': { name: 'Model C (4-variable)', fn: g => [1, g.logVflat, g.logMbar, Math.log10(Math.max(g.logRdisk ? Math.pow(10, g.logRdisk) : sparcMap[g.name] ? sparcMap[g.name].Rdisk : 1, 0.01)), g.morphT].reduce((s, x, j) => s + x * estC.beta[j], 0) },
  'D1_VfMbar': { name: 'Model D1 (logVflat+logMbar)', fn: g => [1, g.logVflat, g.logMbar].reduce((s, x, j) => s + x * estD1.beta[j], 0) },
  'D2_VfT': { name: 'Model D2 (logVflat+T)', fn: g => [1, g.logVflat, g.morphT].reduce((s, x, j) => s + x * estD2.beta[j], 0) },
};

function applyExtMhost(estFn, gal) {
  return estFn({
    logVflat: gal.logVflat,
    logMbar: gal.logMbar,
    logRdisk: gal.logRdisk !== undefined ? gal.logRdisk : Math.log10(sparcMap[gal.name] ? sparcMap[gal.name].Rdisk : 1),
    morphT: gal.morphT !== undefined ? gal.morphT : (sparcMap[gal.name] ? sparcMap[gal.name].T : 5),
    name: gal.name
  });
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 3: APPLY EACH ESTIMATOR TO N=59 EXTERNAL');
console.log('  Re-run full hierarchy test for each');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function predictModel(model, features) {
  return [1, ...features].reduce((s, x, j) => s + x * model.beta[j], 0);
}

function runHierarchy(estKey, estObj, regime, gals, gals5) {
  if (gals.length < 4) return null;

  const galsMod = gals.map(g => {
    const newMhost = applyExtMhost(estObj.fn, g);
    return { ...g, logMhost: newMhost };
  });
  const gals5Mod = gals5.map(g => {
    const newMhost = applyExtMhost(estObj.fn, g);
    return { ...g, logMhost: newMhost };
  });

  const coreModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR]));
  const vfResidModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR, g.VfResid]));
  const fiveAxisModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR, g.VfResid, g.lhOuter]));
  const vfOnlyModel = ols(Y_tr, train.map(g => [g.VfResid]));

  const naive = Math.sqrt(galsMod.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / galsMod.length);
  const naive5 = gals5Mod.length >= 4 ? Math.sqrt(gals5Mod.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / gals5Mod.length) : naive;

  function transfer(model, featureFn, gs, bl) {
    const preds = gs.map(g => predictModel(model, featureFn(g)));
    const rmse = Math.sqrt(gs.reduce((s, g, i) => s + (g.logA0 - preds[i]) ** 2, 0) / gs.length);
    const gap = transferGap(rmse, bl);
    const r = pearsonR(gs.map(g => g.logA0), preds);
    return { rmse: +rmse.toFixed(4), gap: +gap.toFixed(1), r: +r.toFixed(3) };
  }

  const core = transfer(coreModel, g => [g.logMHI, g.logMhost, g.logMeanRun], galsMod, naive);
  const vfr = transfer(vfResidModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid], galsMod, naive);
  const five = gals5Mod.length >= 4
    ? transfer(fiveAxisModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid, g.lhOuter], gals5Mod, naive5)
    : null;
  const vfOnly = transfer(vfOnlyModel, g => [g.VfResid], galsMod, naive);

  const coreBeta = coreModel.beta;
  const vfrBeta = vfResidModel.beta;

  return {
    regime,
    N: gals.length,
    core, vfResid: vfr, fiveAxis: five, vfOnly,
    coreFails: core.gap < 5,
    vfrDominates: vfr.gap > core.gap,
    vfOnlyBetter: vfOnly.gap > core.gap,
    lhAdds: five ? five.gap > vfr.gap : false,
    coreCoeffs: { intercept: +coreBeta[0].toFixed(4), logMHI: +coreBeta[1].toFixed(4), logMhost: +coreBeta[2].toFixed(4), logMR: +coreBeta[3].toFixed(4) },
    vfrCoeffs: { intercept: +vfrBeta[0].toFixed(4), logMHI: +vfrBeta[1].toFixed(4), logMhost: +vfrBeta[2].toFixed(4), logMR: +vfrBeta[3].toFixed(4), VfResid: +vfrBeta[4].toFixed(4) }
  };
}

const allResults = {};

for (const [key, est] of Object.entries(estimators)) {
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  ' + est.name.padEnd(58) + '║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');

  const regimes = [
    { label: 'Full', gals: ext, g5: ext5 },
    { label: 'High-V (>=120)', gals: ext.filter(g => g.Vflat >= 120), g5: ext.filter(g => g.Vflat >= 120 && g.lhOuter !== null) },
    { label: 'Very-High-V (>=180)', gals: ext.filter(g => g.Vflat >= 180), g5: ext.filter(g => g.Vflat >= 180 && g.lhOuter !== null) },
    { label: 'Q1+HV', gals: ext.filter(g => g.Q <= 1 && g.Vflat >= 120), g5: ext.filter(g => g.Q <= 1 && g.Vflat >= 120 && g.lhOuter !== null) },
  ];

  const estResults = {};
  console.log('  ' + 'Regime'.padEnd(22) + 'Core'.padStart(8) + 'VfOnly'.padStart(8) + 'C+Vf'.padStart(8) + '5-ax'.padStart(8) + ' CoreFail VfDom lhAdd');
  console.log('  ' + '─'.repeat(75));

  for (const { label, gals, g5 } of regimes) {
    const h = runHierarchy(key, est, label, gals, g5);
    if (!h) continue;
    estResults[label] = h;

    const fiveStr = h.fiveAxis ? (h.fiveAxis.gap.toFixed(1) + '%').padStart(8) : '   n/a  ';
    console.log('  ' + label.padEnd(22)
      + (h.core.gap.toFixed(1) + '%').padStart(8)
      + (h.vfOnly.gap.toFixed(1) + '%').padStart(8)
      + (h.vfResid.gap.toFixed(1) + '%').padStart(8)
      + fiveStr
      + (h.coreFails ? '  YES' : '   NO').padStart(9)
      + (h.vfrDominates ? '  YES' : '   NO').padStart(6)
      + (h.lhAdds ? '  YES' : '   NO').padStart(6));
  }

  allResults[key] = estResults;
  console.log('');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 4: COMPARISON TABLE — Core gap across estimators');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const compareRegimes = ['Full', 'High-V (>=120)', 'Very-High-V (>=180)', 'Q1+HV'];
const estKeys = Object.keys(estimators);

for (const regime of compareRegimes) {
  console.log('  ' + regime + ':');
  console.log('  ' + 'Estimator'.padEnd(32) + 'Core gap'.padStart(10) + 'C+Vf gap'.padStart(10) + '5-ax gap'.padStart(10) + 'VfOnly'.padStart(10) + ' Delta(Vf-Core)');
  console.log('  ' + '─'.repeat(85));

  for (const key of estKeys) {
    const h = allResults[key][regime];
    if (!h) continue;
    const fiveStr = h.fiveAxis ? (h.fiveAxis.gap.toFixed(1) + '%').padStart(10) : '      n/a ';
    const delta = (h.vfResid.gap - h.core.gap).toFixed(1);
    console.log('  ' + estimators[key].name.padEnd(32)
      + (h.core.gap.toFixed(1) + '%').padStart(10)
      + (h.vfResid.gap.toFixed(1) + '%').padStart(10)
      + fiveStr
      + (h.vfOnly.gap.toFixed(1) + '%').padStart(10)
      + ('  +' + delta + 'pp'));
  }
  console.log('');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 5: SIGN CONSISTENCY CHECK');
console.log('  Internal reference signs: logMHI(+), logMhost(-), logMR(-), VfResid(+)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const coreModelRef = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR]));
const vfrModelRef = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR, g.VfResid]));

console.log('  Internal N=45 coefficients (reference):');
console.log('    Core:    intercept=' + coreModelRef.beta[0].toFixed(3) + '  logMHI=' + coreModelRef.beta[1].toFixed(3) + '  logMhost=' + coreModelRef.beta[2].toFixed(3) + '  logMR=' + coreModelRef.beta[3].toFixed(3));
console.log('    C+VfR:   intercept=' + vfrModelRef.beta[0].toFixed(3) + '  logMHI=' + vfrModelRef.beta[1].toFixed(3) + '  logMhost=' + vfrModelRef.beta[2].toFixed(3) + '  logMR=' + vfrModelRef.beta[3].toFixed(3) + '  VfResid=' + vfrModelRef.beta[4].toFixed(3));
console.log('');

const internalSigns = {
  logMHI: Math.sign(coreModelRef.beta[1]),
  logMhost: Math.sign(coreModelRef.beta[2]),
  logMR: Math.sign(coreModelRef.beta[3])
};

console.log('  External coefficient signs by estimator (Core model, Full sample):');
console.log('  ' + 'Estimator'.padEnd(32) + 'logMHI'.padStart(8) + 'logMhost'.padStart(10) + 'logMR'.padStart(8) + '  Match?');
console.log('  ' + '─'.repeat(65));

const signResults = {};
for (const key of estKeys) {
  const h = allResults[key]['Full'];
  if (!h) continue;
  const c = h.coreCoeffs;
  const mhiSign = Math.sign(c.logMHI) === internalSigns.logMHI;
  const mhostSign = Math.sign(c.logMhost) === internalSigns.logMhost;
  const mrSign = Math.sign(c.logMR) === internalSigns.logMR;
  const matchCount = [mhiSign, mhostSign, mrSign].filter(Boolean).length;
  console.log('  ' + estimators[key].name.padEnd(32)
    + (c.logMHI > 0 ? '+' : '-').padStart(8)
    + (c.logMhost > 0 ? '+' : '-').padStart(10)
    + (c.logMR > 0 ? '+' : '-').padStart(8)
    + '  ' + matchCount + '/3');
  signResults[key] = { mhiSign, mhostSign, mrSign, matchCount };
}
console.log('  Internal reference signs:      ' + (internalSigns.logMHI > 0 ? '+' : '-').padStart(8) + (internalSigns.logMhost > 0 ? '+' : '-').padStart(10) + (internalSigns.logMR > 0 ? '+' : '-').padStart(8));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 6: HIERARCHY SCORE (Phase 202 style) PER ESTIMATOR');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const hierScores = {};
for (const key of estKeys) {
  const r = allResults[key];
  const checks = {
    coreFails_full: r['Full'] ? r['Full'].coreFails : false,
    coreFails_hv: r['High-V (>=120)'] ? r['High-V (>=120)'].coreFails : false,
    vfrDominates_full: r['Full'] ? r['Full'].vfrDominates : false,
    vfrDominates_hv: r['High-V (>=120)'] ? r['High-V (>=120)'].vfrDominates : false,
    vfOnlyBeatsCore_full: r['Full'] ? r['Full'].vfOnlyBetter : false,
    vfOnlyBeatsCore_hv: r['High-V (>=120)'] ? r['High-V (>=120)'].vfOnlyBetter : false,
    regimePreserved: r['High-V (>=120)'] && r['Full'] ? r['High-V (>=120)'].vfResid.gap > r['Full'].vfResid.gap : false,
    lhAdds_hv: r['High-V (>=120)'] ? r['High-V (>=120)'].lhAdds : false,
  };
  const passed = Object.values(checks).filter(Boolean).length;
  hierScores[key] = { checks, passed, total: 8 };

  console.log('  ' + estimators[key].name + ': ' + passed + '/8');
  for (const [ck, v] of Object.entries(checks)) {
    console.log('    ' + ck.padEnd(28) + (v ? 'PASS' : 'FAIL'));
  }
  console.log('');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  VERDICT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const bestEstKey = estKeys.reduce((best, k) => {
  const fullCore = allResults[k]['Full'] ? allResults[k]['Full'].core.gap : -999;
  const bestCore = allResults[best]['Full'] ? allResults[best]['Full'].core.gap : -999;
  return fullCore > bestCore ? k : best;
});

const bestHV = estKeys.reduce((best, k) => {
  const hvCore = allResults[k]['High-V (>=120)'] ? allResults[k]['High-V (>=120)'].core.gap : -999;
  const bestCore = allResults[best]['High-V (>=120)'] ? allResults[best]['High-V (>=120)'].core.gap : -999;
  return hvCore > bestCore ? k : best;
});

const baselineCore = allResults['A_formula']['Full'].core.gap;
const bestCore = allResults[bestEstKey]['Full'].core.gap;
const coreImproved = bestCore > baselineCore + 5;
const vfrStillDominates = allResults[bestEstKey]['Full'].vfrDominates;
const vfrStillDominatesHV = allResults[bestEstKey]['High-V (>=120)'] ? allResults[bestEstKey]['High-V (>=120)'].vfrDominates : true;

let verdict;
if (coreImproved && vfrStillDominates) {
  verdict = 'CORE_IMPROVES_VfResid_STILL_DOMINATES';
} else if (coreImproved && !vfrStillDominates) {
  verdict = 'CORE_RESTORED_HIERARCHY_SHIFTS';
} else if (!coreImproved && vfrStillDominates) {
  verdict = 'CORE_UNCHANGED_VfResid_DOMINANT';
} else {
  verdict = 'NO_CLEAR_IMPROVEMENT';
}

console.log('  Best estimator for Core (full): ' + estimators[bestEstKey].name);
console.log('    Core gap improved:           ' + baselineCore.toFixed(1) + '% → ' + bestCore.toFixed(1) + '% (' + (coreImproved ? 'YES' : 'NO') + ')');
console.log('  Best estimator for Core (HV):  ' + estimators[bestHV].name);
const hvBaseCore = allResults['A_formula']['High-V (>=120)'].core.gap;
const hvBestCore = allResults[bestHV]['High-V (>=120)'].core.gap;
console.log('    Core gap improved (HV):      ' + hvBaseCore.toFixed(1) + '% → ' + hvBestCore.toFixed(1) + '% (' + (hvBestCore > hvBaseCore + 5 ? 'YES' : 'NO') + ')');
console.log('  VfResid still dominates (full): ' + (vfrStillDominates ? 'YES' : 'NO'));
console.log('  VfResid still dominates (HV):   ' + (vfrStillDominatesHV ? 'YES' : 'NO'));
console.log('');

const resA_full = allResults['A_formula']['Full'];
const resBest_full = allResults[bestEstKey]['Full'];
console.log('  Core (full) delta from best estimator:    ' + (bestCore - baselineCore).toFixed(1) + 'pp');
console.log('  VfResid (full) gap unchanged at:          ' + resBest_full.vfResid.gap.toFixed(1) + '%');
console.log('  VfResid margin over Core (full):          ' + (resBest_full.vfResid.gap - resBest_full.core.gap).toFixed(1) + 'pp');
console.log('');
console.log('  VERDICT: ' + verdict);
console.log('');

const q1hvA = allResults['A_formula']['Q1+HV'];
const q1hvBest = allResults[bestEstKey]['Q1+HV'];
if (q1hvA && q1hvBest) {
  console.log('  Best regime (Q1+HV):');
  console.log('    Baseline Core gap:  ' + q1hvA.core.gap.toFixed(1) + '% → Best: ' + q1hvBest.core.gap.toFixed(1) + '%');
  console.log('    C+VfResid gap:      ' + q1hvA.vfResid.gap.toFixed(1) + '% → Best: ' + q1hvBest.vfResid.gap.toFixed(1) + '%');
  console.log('    5-axis gap:         ' + (q1hvA.fiveAxis ? q1hvA.fiveAxis.gap.toFixed(1) + '%' : 'n/a') + ' → Best: ' + (q1hvBest.fiveAxis ? q1hvBest.fiveAxis.gap.toFixed(1) + '%' : 'n/a'));
  console.log('');
}

const output = {
  phase: '203',
  title: 'logMhost Improvement for External Sample',
  verdict,
  estimatorTraining: {
    B: { features: ['logVflat'], beta: estB.beta, rmse: +estB.rmse.toFixed(4), looRMSE: +looBrmse.toFixed(4), r: +rB.toFixed(3) },
    C: { features: ['logVflat', 'logMbar', 'logRdisk', 'T'], beta: estC.beta, rmse: +estC.rmse.toFixed(4), looRMSE: +looCrmse.toFixed(4), r: +rC.toFixed(3) },
    D1: { features: ['logVflat', 'logMbar'], beta: estD1.beta, rmse: +estD1.rmse.toFixed(4), looRMSE: +looD1rmse.toFixed(4), r: +rD1.toFixed(3) },
    D2: { features: ['logVflat', 'T'], beta: estD2.beta, rmse: +estD2.rmse.toFixed(4), looRMSE: +looD2rmse.toFixed(4), r: +rD2.toFixed(3) },
    formulaA: { rmse: +rmseFormulaA.toFixed(4), r: +rFormulaA.toFixed(3) }
  },
  hierarchyByEstimator: allResults,
  hierScores,
  signConsistency: signResults,
  bestEstimatorFull: bestEstKey,
  bestEstimatorHV: bestHV,
  coreImproved,
  coreDelta: +(bestCore - baselineCore).toFixed(1),
  vfrStillDominates,
  vfrStillDominatesHV,
  internalCoreCoeffs: {
    intercept: +coreModelRef.beta[0].toFixed(4),
    logMHI: +coreModelRef.beta[1].toFixed(4),
    logMhost: +coreModelRef.beta[2].toFixed(4),
    logMR: +coreModelRef.beta[3].toFixed(4)
  }
};

fs.writeFileSync(pub('phase203-logMhost-improvement.json'), JSON.stringify(output, null, 2));
console.log('  Output: public/phase203-logMhost-improvement.json');
console.log('\n======================================================================');
console.log('  PHASE 203 COMPLETE');
console.log('======================================================================');

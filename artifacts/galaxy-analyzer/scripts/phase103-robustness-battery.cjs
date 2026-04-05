#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 103: ROBUSTNESS BATTERY');
console.log('  Tests 1-4: Bootstrap, Ablation, Collinearity, Threshold Sweep');
console.log('  Target: logOuterMassDiscrepancy (independent vars only)');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;
const N_BOOTSTRAP = 2000;

let rngState = 42;
function seededRandom() {
  rngState = (rngState * 1664525 + 1013904223) & 0x7fffffff;
  return rngState / 0x7fffffff;
}

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function percentile(a, p) { const s = [...a].sort((x, y) => x - y); const i = p * (s.length - 1); const lo = Math.floor(i); return lo === s.length - 1 ? s[lo] : s[lo] + (i - lo) * (s[lo + 1] - s[lo]); }
function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function olsRegress(X, y) {
  const n = y.length, p = X[0].length;
  if (n <= p + 1) return null;
  const Xt = Array.from({ length: p }, (_, i) => X.map(row => row[i]));
  const XtX = Array.from({ length: p }, (_, i) =>
    Array.from({ length: p }, (_, j) =>
      Xt[i].reduce((s, _, k) => s + Xt[i][k] * Xt[j][k], 0)));
  const Xty = Xt.map(col => col.reduce((s, v, k) => s + v * y[k], 0));
  const aug = XtX.map((row, i) => [...row, Xty[i]]);
  for (let i = 0; i < p; i++) {
    let maxRow = i;
    for (let k = i + 1; k < p; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) maxRow = k;
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];
    if (Math.abs(aug[i][i]) < 1e-12) return null;
    for (let k = 0; k < p; k++) {
      if (k === i) continue;
      const f = aug[k][i] / aug[i][i];
      for (let j = i; j <= p; j++) aug[k][j] -= f * aug[i][j];
    }
  }
  const beta = aug.map((row, i) => row[p] / row[i]);
  const yhat = X.map(row => row.reduce((s, v, j) => s + v * beta[j], 0));
  const ss_res = y.reduce((s, v, i) => s + (v - yhat[i]) ** 2, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, residuals: y.map((v, i) => v - yhat[i]), yhat };
}

function looR2(X_fn, y_fn, data) {
  const n = data.length;
  let ss_res = 0, ss_tot = 0;
  const yAll = data.map(y_fn);
  const yMean = mean(yAll);
  for (let i = 0; i < n; i++) {
    const train = data.filter((_, j) => j !== i);
    const fit = olsRegress(train.map(X_fn), train.map(y_fn));
    if (!fit) return NaN;
    const xTest = X_fn(data[i]);
    const yPred = xTest.reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss_res += (yAll[i] - yPred) ** 2;
    ss_tot += (yAll[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const D = parseFloat(line.substring(15, 21).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  const RHI = parseFloat(line.substring(78, 86).trim());
  table1[name] = { T, D, L36, Rdisk, MHI, Vflat, Q, RHI };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul });
}

const allGalaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0) continue;
  if (t1.Vflat < VFLAT_BREAK) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    const vBar = Math.sqrt(Math.max(vBarSq, 0.01));
    pts.push({ r: pt.rad, vobs: pt.vobs, vBar });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);

  const outerMassDisc = outerPts.map(p => (p.vobs * p.vobs) / (p.vBar * p.vBar));
  const logOMD = Math.log10(mean(outerMassDisc.filter(v => isFinite(v) && v > 0)));
  if (!isFinite(logOMD)) continue;

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const baryonCompact = t1.L36 * UPSILON_DISK / (t1.Rdisk > 0 ? t1.Rdisk : 1);
  const gasExtent = t1.RHI / (t1.Rdisk > 0 ? t1.Rdisk : 1);

  if (!isFinite(Sigma0) || Sigma0 <= 0 || !isFinite(fgas) || !isFinite(t1.RHI) || t1.RHI <= 0) continue;

  allGalaxies.push({
    name, nPts: sorted.length,
    logOMD,
    logVflat: Math.log10(t1.Vflat),
    Vflat: t1.Vflat,
    logMHI: Math.log10(t1.MHI),
    logL36: Math.log10(t1.L36),
    logSigma0: Math.log10(Sigma0),
    fgas,
    logBaryonCompact: Math.log10(baryonCompact > 0 ? baryonCompact : 1e-3),
    logGasExtent: Math.log10(gasExtent > 0 ? gasExtent : 0.1),
    logRdisk: Math.log10(t1.Rdisk > 0 ? t1.Rdisk : 0.1),
    logRHI: Math.log10(t1.RHI > 0 ? t1.RHI : 0.1),
    logD: Math.log10(t1.D > 0 ? t1.D : 1),
    T: t1.T
  });
}

const hi = allGalaxies.filter(g =>
  isFinite(g.logSigma0) && isFinite(g.fgas) && isFinite(g.logMHI) &&
  isFinite(g.logRdisk) && isFinite(g.logBaryonCompact) && isFinite(g.logGasExtent));

console.log('  Valid galaxies: N=' + hi.length + '\n');

const targetFn = g => g.logOMD;
const fullModelX = g => [1, g.logSigma0, g.logMHI, g.logRdisk, g.fgas];
const fullModelVars = ['intercept', 'logSigma0', 'logMHI', 'logRdisk', 'fgas'];

const fullFit = olsRegress(hi.map(fullModelX), hi.map(targetFn));
const fullLOO = looR2(fullModelX, targetFn, hi);
console.log('  Full model (Sigma0+MHI+Rdisk+fgas):');
console.log('    R²=' + fullFit.r2.toFixed(4) + ', LOO R²=' + fullLOO.toFixed(4));
console.log('    Coefficients: ' + fullModelVars.map((v, i) => v + '=' + fullFit.beta[i].toFixed(4)).join(', '));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 1: BOOTSTRAP STABILITY (N=' + N_BOOTSTRAP + ')');
console.log('══════════════════════════════════════════════════════════════\n');

const bootBetas = fullModelVars.map(() => []);
const bootR2s = [];
const bootLOOs = [];

for (let b = 0; b < N_BOOTSTRAP; b++) {
  const sample = [];
  for (let i = 0; i < hi.length; i++) {
    const idx = Math.floor(seededRandom() * hi.length);
    sample.push(hi[idx]);
  }
  const fit = olsRegress(sample.map(fullModelX), sample.map(targetFn));
  if (!fit) continue;
  for (let j = 0; j < fit.beta.length; j++) bootBetas[j].push(fit.beta[j]);
  bootR2s.push(fit.r2);

  if (b < 200) {
    const loo = looR2(fullModelX, targetFn, sample);
    bootLOOs.push(loo);
  }
}

console.log('  Coefficient stability:');
for (let j = 0; j < fullModelVars.length; j++) {
  const arr = bootBetas[j];
  const ci025 = percentile(arr, 0.025);
  const ci975 = percentile(arr, 0.975);
  const m = mean(arr);
  const s = sd(arr);
  const crossesZero = (ci025 <= 0 && ci975 >= 0);
  console.log('    ' + fullModelVars[j].padEnd(12) +
    ' mean=' + m.toFixed(4).padStart(8) +
    ' sd=' + s.toFixed(4).padStart(7) +
    ' 95%CI=[' + ci025.toFixed(4) + ', ' + ci975.toFixed(4) + ']' +
    (crossesZero ? ' ⚠ CROSSES ZERO' : ' ✓ stable sign'));
}

console.log('\n  R² stability:');
console.log('    mean=' + mean(bootR2s).toFixed(4) + ', sd=' + sd(bootR2s).toFixed(4) +
  ', 95%CI=[' + percentile(bootR2s, 0.025).toFixed(4) + ', ' + percentile(bootR2s, 0.975).toFixed(4) + ']');

if (bootLOOs.length > 0) {
  console.log('  LOO R² stability (first 200 boots):');
  console.log('    mean=' + mean(bootLOOs).toFixed(4) + ', sd=' + sd(bootLOOs).toFixed(4) +
    ', 95%CI=[' + percentile(bootLOOs, 0.025).toFixed(4) + ', ' + percentile(bootLOOs, 0.975).toFixed(4) + ']');
}

const stableCoeffs = fullModelVars.filter((_, j) => {
  const ci025 = percentile(bootBetas[j], 0.025);
  const ci975 = percentile(bootBetas[j], 0.975);
  return j === 0 || !(ci025 <= 0 && ci975 >= 0);
}).length - 1;

console.log('\n  VERDICT: ' + stableCoeffs + '/' + (fullModelVars.length - 1) + ' coefficients have stable sign');
console.log('  Bootstrap ' + (stableCoeffs >= 3 ? 'PASS' : stableCoeffs >= 2 ? 'MARGINAL' : 'FAIL'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: ABLATION (leave-one-predictor-out)');
console.log('══════════════════════════════════════════════════════════════\n');

const ablationModels = [
  { name: 'Full (Sigma0+MHI+Rdisk+fgas)', X: fullModelX },
  { name: 'Without fgas', X: g => [1, g.logSigma0, g.logMHI, g.logRdisk] },
  { name: 'Without logSigma0', X: g => [1, g.logMHI, g.logRdisk, g.fgas] },
  { name: 'Without logMHI', X: g => [1, g.logSigma0, g.logRdisk, g.fgas] },
  { name: 'Without logRdisk', X: g => [1, g.logSigma0, g.logMHI, g.fgas] },
  { name: 'fgas alone', X: g => [1, g.fgas] },
  { name: 'logSigma0 alone', X: g => [1, g.logSigma0] },
  { name: 'logBaryonCompact alone', X: g => [1, g.logBaryonCompact] },
  { name: 'fgas + logSigma0 only', X: g => [1, g.fgas, g.logSigma0] },
  { name: 'fgas + logBaryonCompact only', X: g => [1, g.fgas, g.logBaryonCompact] },
];

const ablationResults = [];
for (const model of ablationModels) {
  const fit = olsRegress(hi.map(model.X), hi.map(targetFn));
  const loo = looR2(model.X, targetFn, hi);
  ablationResults.push({ name: model.name, r2: fit ? fit.r2 : NaN, loo });
  console.log('  ' + model.name.padEnd(40) +
    ' R²=' + (fit ? fit.r2.toFixed(3) : 'N/A').padStart(6) +
    '  LOO=' + loo.toFixed(3).padStart(7) +
    '  drop=' + ((fullLOO - loo) * 100 / fullLOO).toFixed(1).padStart(6) + '%');
}

const ablationDrops = ablationResults.slice(1, 5).map(a => ({
  name: a.name,
  drop: ((fullLOO - a.loo) / fullLOO * 100)
}));
ablationDrops.sort((a, b) => b.drop - a.drop);
console.log('\n  Most impactful when removed:');
for (const a of ablationDrops) {
  console.log('    ' + a.name + ': ' + a.drop.toFixed(1) + '% LOO drop');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: COLLINEARITY CHECK');
console.log('══════════════════════════════════════════════════════════════\n');

const varFns = [
  { name: 'logSigma0', fn: g => g.logSigma0 },
  { name: 'logMHI', fn: g => g.logMHI },
  { name: 'logRdisk', fn: g => g.logRdisk },
  { name: 'fgas', fn: g => g.fgas },
  { name: 'logL36', fn: g => g.logL36 },
  { name: 'logBaryonCompact', fn: g => g.logBaryonCompact },
  { name: 'logGasExtent', fn: g => g.logGasExtent },
  { name: 'logVflat', fn: g => g.logVflat },
];

console.log('  Pairwise correlations among predictors:');
console.log('  ' + ''.padEnd(18) + varFns.slice(0, 4).map(v => v.name.padStart(12)).join(''));
for (let i = 0; i < 4; i++) {
  let row = '  ' + varFns[i].name.padEnd(18);
  for (let j = 0; j < 4; j++) {
    const r = pearsonR(hi.map(varFns[i].fn), hi.map(varFns[j].fn));
    row += (i === j ? '  1.000' : r.toFixed(3)).padStart(12);
  }
  console.log(row);
}

console.log('\n  VIF (Variance Inflation Factor):');
const modelVarFns = [
  { name: 'logSigma0', fn: g => g.logSigma0 },
  { name: 'logMHI', fn: g => g.logMHI },
  { name: 'logRdisk', fn: g => g.logRdisk },
  { name: 'fgas', fn: g => g.fgas },
];

const vifResults = [];
for (let k = 0; k < modelVarFns.length; k++) {
  const yv = hi.map(modelVarFns[k].fn);
  const Xv = hi.map(g => {
    const row = [1];
    for (let j = 0; j < modelVarFns.length; j++) {
      if (j !== k) row.push(modelVarFns[j].fn(g));
    }
    return row;
  });
  const fit = olsRegress(Xv, yv);
  const vif = fit && fit.r2 < 1 ? 1 / (1 - fit.r2) : Infinity;
  vifResults.push({ name: modelVarFns[k].name, r2: fit ? fit.r2 : NaN, vif });
  console.log('    ' + modelVarFns[k].name.padEnd(16) + ' VIF=' + vif.toFixed(2).padStart(6) +
    ' (R² with others=' + (fit ? fit.r2.toFixed(3) : 'N/A') + ')' +
    (vif > 10 ? ' ⚠ HIGH' : vif > 5 ? ' ⚠ MODERATE' : ' ✓ OK'));
}

const maxVIF = Math.max(...vifResults.map(v => v.vif));
console.log('\n  Max VIF: ' + maxVIF.toFixed(2));
console.log('  Collinearity ' + (maxVIF < 5 ? 'PASS - no severe collinearity' :
  maxVIF < 10 ? 'MARGINAL - some collinearity but manageable' :
  'FAIL - severe collinearity, model interpretation compromised'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: THRESHOLD SWEEP (Vflat cut)');
console.log('══════════════════════════════════════════════════════════════\n');

const thresholds = [0, 40, 50, 60, 70, 80, 90, 100, 120];
const sweepResults = [];
for (const thresh of thresholds) {
  const subset = allGalaxies.filter(g => g.Vflat >= thresh &&
    isFinite(g.logSigma0) && isFinite(g.fgas) && isFinite(g.logMHI) &&
    isFinite(g.logRdisk) && isFinite(g.logBaryonCompact));

  if (subset.length < 15) {
    console.log('  Vflat >= ' + (thresh + ' km/s').padEnd(10) + ' N=' + (subset.length + '').padStart(4) + '  (too few)');
    sweepResults.push({ threshold: thresh, n: subset.length, loo: NaN });
    continue;
  }

  const fit = olsRegress(subset.map(fullModelX), subset.map(targetFn));
  const loo = looR2(fullModelX, targetFn, subset);

  const rFgas = pearsonR(subset.map(g => g.fgas), subset.map(targetFn));

  console.log('  Vflat >= ' + (thresh + ' km/s').padEnd(10) +
    ' N=' + (subset.length + '').padStart(4) +
    '  R²=' + (fit ? fit.r2.toFixed(3) : 'N/A').padStart(6) +
    '  LOO=' + loo.toFixed(3).padStart(7) +
    '  r(fgas)=' + rFgas.toFixed(3).padStart(7));

  sweepResults.push({ threshold: thresh, n: subset.length, r2: fit ? fit.r2 : NaN, loo, rFgas });
}

const mainThreshLOO = sweepResults.find(s => s.threshold === 70);
const stableRange = sweepResults.filter(s => s.loo > 0.3);
console.log('\n  Stable range: Vflat >= ' + (stableRange.length > 0 ? stableRange[0].threshold : 'N/A') + ' km/s');
console.log('  At Vflat >= 70: LOO R² = ' + (mainThreshLOO ? mainThreshLOO.loo.toFixed(3) : 'N/A'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: INFLUENCE / OUTLIER DIAGNOSTICS');
console.log('══════════════════════════════════════════════════════════════\n');

const yAll = hi.map(targetFn);
const yMean = mean(yAll);
const cooksD = [];
const jackR2 = [];
const mse_full = fullFit.residuals.reduce((s, r) => s + r * r, 0) / (hi.length - fullModelVars.length);

for (let i = 0; i < hi.length; i++) {
  const train = hi.filter((_, j) => j !== i);
  const fit_i = olsRegress(train.map(fullModelX), train.map(targetFn));
  if (!fit_i) { cooksD.push(0); jackR2.push(NaN); continue; }

  const yhat_full_i = fullModelX(hi[i]).reduce((s, v, j) => s + v * fullFit.beta[j], 0);
  const yhat_jack_i = fullModelX(hi[i]).reduce((s, v, j) => s + v * fit_i.beta[j], 0);
  const diff = yhat_full_i - yhat_jack_i;
  const cook = (diff * diff) / (fullModelVars.length * mse_full);
  cooksD.push(cook);

  const ss_res_j = train.map(targetFn).reduce((s, v, k) => s + (v - fit_i.yhat[k]) ** 2, 0);
  const ss_tot_j = train.map(targetFn).reduce((s, v) => s + (v - mean(train.map(targetFn))) ** 2, 0);
  jackR2.push(ss_tot_j > 0 ? 1 - ss_res_j / ss_tot_j : 0);
}

const cooksSorted = cooksD.map((c, i) => ({ name: hi[i].name, cook: c, idx: i }))
  .sort((a, b) => b.cook - a.cook);

console.log('  Top 10 influential galaxies (Cook\'s D):');
for (let i = 0; i < Math.min(10, cooksSorted.length); i++) {
  const g = cooksSorted[i];
  console.log('    ' + (i + 1 + '.').padEnd(4) + g.name.padEnd(14) +
    ' Cook\'s D=' + g.cook.toFixed(4).padStart(8) +
    (g.cook > 4 / hi.length ? ' ⚠ influential' : ' ✓ OK'));
}

const highInfluence = cooksSorted.filter(g => g.cook > 4 / hi.length);
console.log('\n  Influential points (Cook\'s D > 4/N=' + (4 / hi.length).toFixed(4) + '): ' + highInfluence.length);

const without5 = hi.filter((_, i) => !cooksSorted.slice(0, 5).some(c => c.idx === i));
const fit_w5 = olsRegress(without5.map(fullModelX), without5.map(targetFn));
const loo_w5 = looR2(fullModelX, targetFn, without5);

const without10 = hi.filter((_, i) => !cooksSorted.slice(0, 10).some(c => c.idx === i));
const fit_w10 = olsRegress(without10.map(fullModelX), without10.map(targetFn));
const loo_w10 = looR2(fullModelX, targetFn, without10);

console.log('\n  After removing top 5 influential:');
console.log('    N=' + without5.length + ', R²=' + fit_w5.r2.toFixed(3) + ', LOO R²=' + loo_w5.toFixed(3));
console.log('  After removing top 10 influential:');
console.log('    N=' + without10.length + ', R²=' + fit_w10.r2.toFixed(3) + ', LOO R²=' + loo_w10.toFixed(3));
console.log('  Change: ' + ((loo_w5 - fullLOO) / Math.abs(fullLOO) * 100).toFixed(1) + '% (rm 5), ' +
  ((loo_w10 - fullLOO) / Math.abs(fullLOO) * 100).toFixed(1) + '% (rm 10)');

const influenceStable = Math.abs(loo_w5 - fullLOO) / Math.abs(fullLOO) < 0.25;
console.log('  Influence ' + (influenceStable ? 'PASS - result survives removal of top influencers' :
  'FAIL - result depends heavily on few galaxies'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  OVERALL SUMMARY');
console.log('══════════════════════════════════════════════════════════════\n');

const bootPass = stableCoeffs >= 3;
const ablBest2 = ablationResults.find(a => a.name === 'fgas + logSigma0 only');
const ablPass = ablBest2 && ablBest2.loo > 0.35;
const collinPass = maxVIF < 10;
const threshPass = sweepResults.filter(s => s.threshold >= 50 && s.threshold <= 100 && s.loo > 0.3).length >= 3;
const influencePass = influenceStable;

const tests = [
  { name: 'Bootstrap stability', pass: bootPass },
  { name: 'Ablation (2-var sufficient)', pass: ablPass },
  { name: 'Collinearity (VIF<10)', pass: collinPass },
  { name: 'Threshold robustness', pass: threshPass },
  { name: 'Influence diagnostics', pass: influencePass },
];

let totalPass = 0;
for (const t of tests) {
  console.log('  ' + (t.pass ? 'PASS' : 'FAIL') + ': ' + t.name);
  if (t.pass) totalPass++;
}

console.log('\n  Score: ' + totalPass + '/5');
console.log('  Status: ' + (totalPass >= 4 ? 'ROBUST — result is solid' :
  totalPass >= 3 ? 'MOSTLY ROBUST — minor concerns remain' :
  'FRAGILE — significant issues found'));

const results = {
  phase: 103,
  title: 'Robustness Battery: Bootstrap, Ablation, Collinearity, Threshold, Influence',
  nGalaxies: hi.length,
  fullModel: {
    vars: fullModelVars.slice(1),
    R2: parseFloat(fullFit.r2.toFixed(4)),
    LOO_R2: parseFloat(fullLOO.toFixed(4)),
    coefficients: Object.fromEntries(fullModelVars.map((v, i) => [v, parseFloat(fullFit.beta[i].toFixed(4))]))
  },
  bootstrap: {
    nTrials: N_BOOTSTRAP,
    coefficients: Object.fromEntries(fullModelVars.map((v, j) => [v, {
      mean: parseFloat(mean(bootBetas[j]).toFixed(4)),
      sd: parseFloat(sd(bootBetas[j]).toFixed(4)),
      ci025: parseFloat(percentile(bootBetas[j], 0.025).toFixed(4)),
      ci975: parseFloat(percentile(bootBetas[j], 0.975).toFixed(4)),
      crossesZero: j > 0 && percentile(bootBetas[j], 0.025) <= 0 && percentile(bootBetas[j], 0.975) >= 0
    }])),
    R2: { mean: parseFloat(mean(bootR2s).toFixed(4)), sd: parseFloat(sd(bootR2s).toFixed(4)),
      ci025: parseFloat(percentile(bootR2s, 0.025).toFixed(4)), ci975: parseFloat(percentile(bootR2s, 0.975).toFixed(4)) },
    stableCoefficients: stableCoeffs + '/' + (fullModelVars.length - 1),
    verdict: bootPass ? 'PASS' : 'FAIL'
  },
  ablation: ablationResults.map(a => ({
    name: a.name, R2: parseFloat((a.r2 || 0).toFixed(3)), LOO_R2: parseFloat(a.loo.toFixed(3)),
    dropPct: parseFloat(((fullLOO - a.loo) / fullLOO * 100).toFixed(1))
  })),
  collinearity: {
    VIF: Object.fromEntries(vifResults.map(v => [v.name, parseFloat(v.vif.toFixed(2))])),
    maxVIF: parseFloat(maxVIF.toFixed(2)),
    verdict: collinPass ? 'PASS' : 'FAIL'
  },
  thresholdSweep: sweepResults.map(s => ({
    threshold: s.threshold, n: s.n,
    R2: s.r2 ? parseFloat(s.r2.toFixed(3)) : null,
    LOO_R2: isFinite(s.loo) ? parseFloat(s.loo.toFixed(3)) : null,
    rFgas: s.rFgas ? parseFloat(s.rFgas.toFixed(3)) : null
  })),
  influence: {
    top5: cooksSorted.slice(0, 5).map(g => ({ name: g.name, cooksD: parseFloat(g.cook.toFixed(4)) })),
    nInfluential: highInfluence.length,
    afterRemove5: { n: without5.length, R2: parseFloat(fit_w5.r2.toFixed(3)), LOO_R2: parseFloat(loo_w5.toFixed(3)) },
    afterRemove10: { n: without10.length, R2: parseFloat(fit_w10.r2.toFixed(3)), LOO_R2: parseFloat(loo_w10.toFixed(3)) },
    verdict: influencePass ? 'PASS' : 'FAIL'
  },
  overallScore: totalPass + '/5',
  overallVerdict: totalPass >= 4 ? 'ROBUST' : totalPass >= 3 ? 'MOSTLY ROBUST' : 'FRAGILE'
};

const outPath = path.join(__dirname, '..', 'public', 'phase103-robustness-battery.json');
fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
console.log('\n  Results saved to: ' + outPath);
console.log('\n  DONE.');

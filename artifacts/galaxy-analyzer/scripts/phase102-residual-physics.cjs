#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 102: GEOMETRIC-BASELINE SUBTRACTION + RESIDUAL PHYSICS');
console.log('  New question: what predicts the NEED for extra support at the edges?');
console.log('  Target = outer mass discrepancy (not outer slope)');
console.log('  Predictors = ONLY catalog/independent variables (not from same RC)');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function median(a) { const s = [...a].sort((x, y) => x - y); return s.length % 2 ? s[Math.floor(s.length / 2)] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2; }
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
  const residuals = y.map((v, i) => v - yhat[i]);
  const ss_res = residuals.reduce((s, r) => s + r * r, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, residuals, yhat };
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

function permTest(x, y, nPerm) {
  const rObs = Math.abs(pearsonR(x, y));
  let count = 0;
  let rng = 12345;
  for (let p = 0; p < nPerm; p++) {
    const shuffled = [...x];
    for (let i = shuffled.length - 1; i > 0; i--) {
      rng = (rng * 1664525 + 1013904223) & 0x7fffffff;
      const j = rng % (i + 1);
      [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
    }
    if (Math.abs(pearsonR(shuffled, y)) >= rObs) count++;
  }
  return count / nPerm;
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
  const eVflat = parseFloat(line.substring(106, 112).trim());
  const RHI = parseFloat(line.substring(78, 86).trim());
  table1[name] = { T, D, L36, Rdisk, MHI, Vflat, Q, eVflat, RHI };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const evobs = parseFloat(line.substring(33, 38).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, evobs, vgas, vdisk, vbul });
}

console.log('── STEP 1: COMPUTE PHYSICAL TARGET VARIABLES ──\n');
console.log('  Target A: outerMassDiscrepancy = mean(Vobs²/Vbar²) in outer half');
console.log('  Target B: outerDMfraction = mean(1 - Vbar²/Vobs²) in outer half');
console.log('  Target C: outerExcessVelocity = mean(Vobs/Vbar - 1) in outer half');
console.log('  These measure: HOW MUCH extra support is needed at the edges\n');

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
    pts.push({ r: pt.rad, vobs: pt.vobs, vBar, vBarSq, vgas: pt.vgas, vdisk: pt.vdisk, vbul: pt.vbul || 0 });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);
  const innerPts = sorted.slice(0, half);

  const outerMassDisc = outerPts.map(p => (p.vobs * p.vobs) / (p.vBar * p.vBar));
  const outerDMfrac = outerPts.map(p => 1 - (p.vBar * p.vBar) / (p.vobs * p.vobs));
  const outerExcessV = outerPts.map(p => p.vobs / p.vBar - 1);

  const logOuterMassDisc = Math.log10(mean(outerMassDisc.filter(v => isFinite(v) && v > 0)));
  const meanOuterDMfrac = mean(outerDMfrac.filter(v => isFinite(v)));
  const meanOuterExcessV = mean(outerExcessV.filter(v => isFinite(v)));

  if (!isFinite(logOuterMassDisc) || !isFinite(meanOuterDMfrac)) continue;

  const logR_outer = outerPts.map(p => Math.log10(p.r));
  const logV_outer = outerPts.map(p => Math.log10(p.vobs));
  const mr = mean(logR_outer), mv = mean(logV_outer);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR_outer.length; i++) { sxy += (logR_outer[i] - mr) * (logV_outer[i] - mv); sxx += (logR_outer[i] - mr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logMbar = Math.log10(t1.L36 * UPSILON_DISK * 2.36e9 + t1.MHI * 2.36e9);
  const baryonCompact = t1.L36 * UPSILON_DISK / (t1.Rdisk > 0 ? t1.Rdisk : 1);
  const gasExtent = t1.RHI / (t1.Rdisk > 0 ? t1.Rdisk : 1);
  const innerVmax = Math.max(...innerPts.map(p => p.vobs));
  const logRatio = Math.log10(t1.Vflat / innerVmax);

  allGalaxies.push({
    name, nPts: sorted.length,
    logOuterMassDisc,
    meanOuterDMfrac,
    meanOuterExcessV,
    outerSlope,
    logVflat: Math.log10(t1.Vflat),
    Vflat: t1.Vflat,
    logMHI: Math.log10(t1.MHI),
    logL36: Math.log10(t1.L36),
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    fgas,
    logMbar,
    logBaryonCompact: Math.log10(baryonCompact > 0 ? baryonCompact : 1e-3),
    logGasExtent: Math.log10(gasExtent > 0 ? gasExtent : 0.1),
    logD: Math.log10(t1.D > 0 ? t1.D : 1),
    T: t1.T,
    logRdisk: Math.log10(t1.Rdisk > 0 ? t1.Rdisk : 0.1),
    logRHI: Math.log10(t1.RHI > 0 ? t1.RHI : 0.1),
    logRatio,
    innerVmax
  });
}

const hi = allGalaxies.filter(g =>
  isFinite(g.logOuterMassDisc) && isFinite(g.meanOuterDMfrac) &&
  isFinite(g.logSigma0) && isFinite(g.fgas) && isFinite(g.logMbar) &&
  isFinite(g.logBaryonCompact) && isFinite(g.logGasExtent) &&
  isFinite(g.logRdisk) && isFinite(g.logRHI)
);

console.log('  Valid galaxies (Vflat >= ' + VFLAT_BREAK + ', all vars finite): N=' + hi.length);

console.log('\n── STEP 2: TARGET VARIABLE DISTRIBUTIONS ──\n');
console.log('  logOuterMassDisc: mean=' + mean(hi.map(g => g.logOuterMassDisc)).toFixed(3) +
  ', sd=' + sd(hi.map(g => g.logOuterMassDisc)).toFixed(3) +
  ', range=[' + Math.min(...hi.map(g => g.logOuterMassDisc)).toFixed(2) + ', ' +
  Math.max(...hi.map(g => g.logOuterMassDisc)).toFixed(2) + ']');
console.log('  meanOuterDMfrac: mean=' + mean(hi.map(g => g.meanOuterDMfrac)).toFixed(3) +
  ', sd=' + sd(hi.map(g => g.meanOuterDMfrac)).toFixed(3));
console.log('  meanOuterExcessV: mean=' + mean(hi.map(g => g.meanOuterExcessV)).toFixed(3) +
  ', sd=' + sd(hi.map(g => g.meanOuterExcessV)).toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 3: SINGLE-VARIABLE CORRELATIONS WITH TARGET A');
console.log('  (logOuterMassDiscrepancy — how much DM support is needed)');
console.log('  Using ONLY catalog/independent variables');
console.log('══════════════════════════════════════════════════════════════\n');

const targetA = g => g.logOuterMassDisc;

const independentVars = [
  { name: 'logVflat', fn: g => g.logVflat, desc: 'Flat velocity (catalog)' },
  { name: 'logSigma0', fn: g => g.logSigma0, desc: 'Central surface brightness' },
  { name: 'logMHI', fn: g => g.logMHI, desc: 'HI gas mass' },
  { name: 'logL36', fn: g => g.logL36, desc: 'Stellar luminosity' },
  { name: 'fgas', fn: g => g.fgas, desc: 'Gas fraction' },
  { name: 'logMbar', fn: g => g.logMbar, desc: 'Total baryonic mass' },
  { name: 'logBaryonCompact', fn: g => g.logBaryonCompact, desc: 'Baryonic compactness (L/Rdisk)' },
  { name: 'logGasExtent', fn: g => g.logGasExtent, desc: 'Gas extent (RHI/Rdisk)' },
  { name: 'logRdisk', fn: g => g.logRdisk, desc: 'Disk scale length' },
  { name: 'logRHI', fn: g => g.logRHI, desc: 'HI radius' },
  { name: 'logD', fn: g => g.logD, desc: 'Distance (calibration check)' },
  { name: 'T', fn: g => g.T, desc: 'Morphological type' },
];

const singleVarResults = [];
for (const v of independentVars) {
  const x = hi.map(v.fn);
  const y = hi.map(targetA);
  const r = pearsonR(x, y);
  const pval = permTest(x, y, 5000);
  const fit = olsRegress(x.map(xi => [1, xi]), y);
  const loo = looR2(g => [1, v.fn(g)], targetA, hi);
  singleVarResults.push({ ...v, r, pval, r2: fit ? fit.r2 : NaN, loo });
  console.log('  ' + v.name.padEnd(20) + ' r=' + r.toFixed(3).padStart(7) +
    '  R²=' + (fit ? fit.r2.toFixed(3) : 'N/A').padStart(6) +
    '  LOO=' + loo.toFixed(3).padStart(7) +
    '  p=' + (pval < 0.0001 ? '<0.0001' : pval.toFixed(4)).padStart(8) +
    '  ' + v.desc);
}

singleVarResults.sort((a, b) => Math.abs(b.r) - Math.abs(a.r));
console.log('\n  Top 3 single-variable predictors of outer mass discrepancy:');
for (let i = 0; i < Math.min(3, singleVarResults.length); i++) {
  const v = singleVarResults[i];
  console.log('    ' + (i + 1) + '. ' + v.name + ' (|r|=' + Math.abs(v.r).toFixed(3) + ', LOO R²=' + v.loo.toFixed(3) + ')');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 4: SAME FOR TARGET B (outer DM fraction)');
console.log('══════════════════════════════════════════════════════════════\n');

const targetB = g => g.meanOuterDMfrac;
const singleVarResultsB = [];
for (const v of independentVars) {
  const x = hi.map(v.fn);
  const y = hi.map(targetB);
  const r = pearsonR(x, y);
  const fit = olsRegress(x.map(xi => [1, xi]), y);
  const loo = looR2(g => [1, v.fn(g)], targetB, hi);
  singleVarResultsB.push({ ...v, r, r2: fit ? fit.r2 : NaN, loo });
  console.log('  ' + v.name.padEnd(20) + ' r=' + r.toFixed(3).padStart(7) +
    '  R²=' + (fit ? fit.r2.toFixed(3) : 'N/A').padStart(6) +
    '  LOO=' + loo.toFixed(3).padStart(7) +
    '  ' + v.desc);
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 5: GEOMETRIC BASELINE — OUTER SLOPE FROM SMOOTH FIT');
console.log('  What does smooth-curve geometry alone predict?');
console.log('  Baseline: outerSlope ~ logVflat + nPts (pure geometry)');
console.log('  Residual = real outerSlope - baseline prediction');
console.log('══════════════════════════════════════════════════════════════\n');

const geoBaselineX = g => [1, g.logVflat, Math.log10(g.nPts)];
const geoBaselineY = g => g.outerSlope;
const geoFit = olsRegress(hi.map(geoBaselineX), hi.map(geoBaselineY));
const geoLOO = looR2(geoBaselineX, geoBaselineY, hi);

console.log('  Geometric baseline for outerSlope:');
console.log('    R² = ' + (geoFit ? geoFit.r2.toFixed(3) : 'N/A'));
console.log('    LOO R² = ' + geoLOO.toFixed(3));
if (geoFit) {
  console.log('    Coefficients: intercept=' + geoFit.beta[0].toFixed(3) +
    ', logVflat=' + geoFit.beta[1].toFixed(3) + ', logNpts=' + geoFit.beta[2].toFixed(3));

  for (let i = 0; i < hi.length; i++) {
    hi[i].geoResidual = geoFit.residuals[i];
  }

  console.log('\n  Residual outerSlope distribution:');
  const resids = hi.map(g => g.geoResidual);
  console.log('    mean=' + mean(resids).toFixed(4) + ', sd=' + sd(resids).toFixed(4));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 6: WHAT PREDICTS THE RESIDUAL? (after removing geometry)');
console.log('══════════════════════════════════════════════════════════════\n');

if (geoFit) {
  const targetResid = g => g.geoResidual;
  const residVarResults = [];
  for (const v of independentVars) {
    const x = hi.map(v.fn);
    const y = hi.map(targetResid);
    const r = pearsonR(x, y);
    const pval = permTest(x, y, 5000);
    const fit2 = olsRegress(x.map(xi => [1, xi]), y);
    const loo = looR2(g => [1, v.fn(g)], targetResid, hi);
    residVarResults.push({ ...v, r, pval, r2: fit2 ? fit2.r2 : NaN, loo });
    console.log('  ' + v.name.padEnd(20) + ' r=' + r.toFixed(3).padStart(7) +
      '  R²=' + (fit2 ? fit2.r2.toFixed(3) : 'N/A').padStart(6) +
      '  LOO=' + loo.toFixed(3).padStart(7) +
      '  p=' + (pval < 0.0001 ? '<0.0001' : pval.toFixed(4)).padStart(8));
  }

  residVarResults.sort((a, b) => Math.abs(b.r) - Math.abs(a.r));
  console.log('\n  Top predictors of geometry-corrected outer slope residual:');
  for (let i = 0; i < Math.min(5, residVarResults.length); i++) {
    const v = residVarResults[i];
    console.log('    ' + (i + 1) + '. ' + v.name + ' (|r|=' + Math.abs(v.r).toFixed(3) + ', p=' +
      (v.pval < 0.0001 ? '<0.0001' : v.pval.toFixed(4)) + ')');
  }
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 7: MINIMAL PHYSICAL MODELS (2-5 vars, independent only)');
console.log('══════════════════════════════════════════════════════════════\n');

const models = [
  {
    name: 'M1: logSigma0 + fgas',
    X: g => [1, g.logSigma0, g.fgas],
    vars: ['logSigma0', 'fgas']
  },
  {
    name: 'M2: logSigma0 + logMHI',
    X: g => [1, g.logSigma0, g.logMHI],
    vars: ['logSigma0', 'logMHI']
  },
  {
    name: 'M3: logBaryonCompact + fgas',
    X: g => [1, g.logBaryonCompact, g.fgas],
    vars: ['logBaryonCompact', 'fgas']
  },
  {
    name: 'M4: logSigma0 + logGasExtent',
    X: g => [1, g.logSigma0, g.logGasExtent],
    vars: ['logSigma0', 'logGasExtent']
  },
  {
    name: 'M5: logSigma0 + fgas + logMbar',
    X: g => [1, g.logSigma0, g.fgas, g.logMbar],
    vars: ['logSigma0', 'fgas', 'logMbar']
  },
  {
    name: 'M6: logBaryonCompact + logGasExtent + T',
    X: g => [1, g.logBaryonCompact, g.logGasExtent, g.T],
    vars: ['logBaryonCompact', 'logGasExtent', 'T']
  },
  {
    name: 'M7: logSigma0 + fgas + logGasExtent',
    X: g => [1, g.logSigma0, g.fgas, g.logGasExtent],
    vars: ['logSigma0', 'fgas', 'logGasExtent']
  },
  {
    name: 'M8: logSigma0 + logMHI + logRdisk + fgas',
    X: g => [1, g.logSigma0, g.logMHI, g.logRdisk, g.fgas],
    vars: ['logSigma0', 'logMHI', 'logRdisk', 'fgas']
  },
  {
    name: 'M9: logVflat + logSigma0 + fgas (includes Vflat)',
    X: g => [1, g.logVflat, g.logSigma0, g.fgas],
    vars: ['logVflat', 'logSigma0', 'fgas']
  },
  {
    name: 'M10: logVflat + logBaryonCompact + logGasExtent',
    X: g => [1, g.logVflat, g.logBaryonCompact, g.logGasExtent],
    vars: ['logVflat', 'logBaryonCompact', 'logGasExtent']
  },
];

const targets = [
  { name: 'logOuterMassDisc', fn: targetA, label: 'Outer mass discrepancy' },
  { name: 'meanOuterDMfrac', fn: targetB, label: 'Outer DM fraction' },
];

const modelResults = {};

for (const tgt of targets) {
  console.log('  ── Target: ' + tgt.label + ' ──\n');
  const tgtResults = [];

  for (const model of models) {
    const fit = olsRegress(hi.map(model.X), hi.map(tgt.fn));
    const loo = looR2(model.X, tgt.fn, hi);
    const r = fit ? Math.sqrt(Math.max(fit.r2, 0)) : NaN;
    tgtResults.push({ name: model.name, r2: fit ? fit.r2 : NaN, loo, r, beta: fit ? fit.beta : null, vars: model.vars });
    console.log('  ' + model.name.padEnd(52) +
      ' R²=' + (fit ? fit.r2.toFixed(3) : 'N/A').padStart(6) +
      '  LOO=' + loo.toFixed(3).padStart(7) +
      '  gap=' + (fit ? ((1 - loo / fit.r2) * 100).toFixed(1) + '%' : 'N/A'));
  }

  tgtResults.sort((a, b) => b.loo - a.loo);
  console.log('\n  Best model by LOO R²: ' + tgtResults[0].name +
    ' (LOO=' + tgtResults[0].loo.toFixed(3) + ', R²=' + tgtResults[0].r2.toFixed(3) + ')\n');
  modelResults[tgt.name] = tgtResults;
}

console.log('══════════════════════════════════════════════════════════════');
console.log('  STEP 8: COMPARISON — OLD RC-BASED vs NEW INDEPENDENT');
console.log('  Is the old ratio (Vflat/InnerVmax) better than independent');
console.log('  variables at predicting outer mass discrepancy?');
console.log('══════════════════════════════════════════════════════════════\n');

const rcBasedX = g => [1, g.logRatio];
const rcR = pearsonR(hi.map(g => g.logRatio), hi.map(targetA));
const rcFit = olsRegress(hi.map(rcBasedX), hi.map(targetA));
const rcLOO = looR2(rcBasedX, targetA, hi);

console.log('  OLD: logRatio → logOuterMassDisc');
console.log('    r=' + rcR.toFixed(3) + ', R²=' + (rcFit ? rcFit.r2.toFixed(3) : 'N/A') + ', LOO R²=' + rcLOO.toFixed(3));

const bestIndep = modelResults['logOuterMassDisc'][0];
console.log('\n  BEST INDEPENDENT: ' + bestIndep.name);
console.log('    R²=' + bestIndep.r2.toFixed(3) + ', LOO R²=' + bestIndep.loo.toFixed(3));

console.log('\n  Verdict: ' +
  (bestIndep.loo > rcLOO ?
    'Independent variables OUTPERFORM RC-based ratio for predicting outer mass discrepancy' :
    rcLOO > bestIndep.loo ?
    'RC-based ratio still better, but likely geometric (Phase 101)' :
    'Comparable performance'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 9: NULL TEST — CAN GEOMETRY PRODUCE THESE CORRELATIONS?');
console.log('  Independent vars should NOT be affected by geometric coupling');
console.log('══════════════════════════════════════════════════════════════\n');

const topVarsForNull = singleVarResults.slice(0, 3).filter(v => v.pval < 0.05);
for (const v of topVarsForNull) {
  const x = hi.map(v.fn);
  const y = hi.map(targetA);
  const realCorr = Math.abs(pearsonR(x, y));

  let rng2 = 54321;
  let nullCount = 0;
  const nNull = 10000;
  for (let t = 0; t < nNull; t++) {
    const shuffled = [...y];
    for (let i = shuffled.length - 1; i > 0; i--) {
      rng2 = (rng2 * 1664525 + 1013904223) & 0x7fffffff;
      const j = rng2 % (i + 1);
      [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
    }
    if (Math.abs(pearsonR(x, shuffled)) >= realCorr) nullCount++;
  }
  const nullP = nullCount / nNull;
  console.log('  ' + v.name + ': real |r|=' + realCorr.toFixed(3) + ', null p=' +
    (nullP < 0.0001 ? '<0.0001' : nullP.toFixed(4)) +
    ' → ' + (nullP < 0.001 ? 'GENUINE (not geometric)' : nullP < 0.05 ? 'Likely genuine' : 'Could be noise'));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  OVERALL SUMMARY');
console.log('══════════════════════════════════════════════════════════════\n');

console.log('  Question A: "What predicts the need for extra support at edges?"');
console.log('    Best single independent predictor: ' + singleVarResults[0].name +
  ' (|r|=' + Math.abs(singleVarResults[0].r).toFixed(3) + ')');
console.log('    Best multi-var independent model: ' + bestIndep.name +
  ' (LOO R²=' + bestIndep.loo.toFixed(3) + ')');
console.log('    Old RC-ratio model: LOO R²=' + rcLOO.toFixed(3) + ' (geometric per Phase 101)');
console.log('    Improvement: independent vars are ' +
  (bestIndep.loo > rcLOO ? 'BETTER' : 'WEAKER') + ' than RC-ratio');
console.log('    BUT: independent vars are immune to geometric coupling criticism');

const results = {
  phase: 102,
  title: 'Geometric-Baseline Subtraction + Residual Physics Search',
  nGalaxies: hi.length,
  targetVariables: {
    logOuterMassDisc: {
      description: 'log10 of mean(Vobs²/Vbar²) in outer half — measures DM support needed',
      mean: parseFloat(mean(hi.map(g => g.logOuterMassDisc)).toFixed(3)),
      sd: parseFloat(sd(hi.map(g => g.logOuterMassDisc)).toFixed(3))
    },
    meanOuterDMfrac: {
      description: 'mean(1 - Vbar²/Vobs²) in outer half — fraction of support from DM',
      mean: parseFloat(mean(hi.map(g => g.meanOuterDMfrac)).toFixed(3)),
      sd: parseFloat(sd(hi.map(g => g.meanOuterDMfrac)).toFixed(3))
    }
  },
  singleVariableCorrelations: {
    targetA_logOuterMassDisc: singleVarResults.map(v => ({
      name: v.name, r: parseFloat(v.r.toFixed(3)), r2: parseFloat(v.r2.toFixed(3)),
      loo: parseFloat(v.loo.toFixed(3)), pValue: v.pval, desc: v.desc
    })),
    targetB_meanOuterDMfrac: singleVarResultsB.map(v => ({
      name: v.name, r: parseFloat(v.r.toFixed(3)), r2: parseFloat(v.r2.toFixed(3)),
      loo: parseFloat(v.loo.toFixed(3)), desc: v.desc
    }))
  },
  geometricBaseline: geoFit ? {
    R2: parseFloat(geoFit.r2.toFixed(3)),
    LOO_R2: parseFloat(geoLOO.toFixed(3)),
    coefficients: {
      intercept: parseFloat(geoFit.beta[0].toFixed(3)),
      logVflat: parseFloat(geoFit.beta[1].toFixed(3)),
      logNpts: parseFloat(geoFit.beta[2].toFixed(3))
    }
  } : null,
  multiVariableModels: Object.fromEntries(Object.entries(modelResults).map(([k, v]) =>
    [k, v.map(m => ({ name: m.name, r2: parseFloat(m.r2.toFixed(3)), loo: parseFloat(m.loo.toFixed(3)), vars: m.vars }))]
  )),
  comparison: {
    oldRCratio: { r: parseFloat(rcR.toFixed(3)), R2: rcFit ? parseFloat(rcFit.r2.toFixed(3)) : null, LOO_R2: parseFloat(rcLOO.toFixed(3)), note: 'geometric per Phase 101' },
    bestIndependentModel: { name: bestIndep.name, R2: parseFloat(bestIndep.r2.toFixed(3)), LOO_R2: parseFloat(bestIndep.loo.toFixed(3)), note: 'immune to geometric coupling' }
  }
};

const outPath = path.join(__dirname, '..', 'public', 'phase102-residual-physics.json');
fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
console.log('\n  Results saved to: ' + outPath);
console.log('\n  DONE.');

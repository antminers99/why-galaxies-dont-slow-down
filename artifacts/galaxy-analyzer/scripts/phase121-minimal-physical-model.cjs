#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 121: MINIMAL PHYSICAL MODEL');
console.log('');
console.log('  Goal: Build the SMALLEST physical model (2-3 vars) that');
console.log('  still captures the gas-to-stellar balance signal.');
console.log('');
console.log('  NOT the best fit — the simplest PHYSICAL PICTURE.');
console.log('');
console.log('  Candidate variables:');
console.log('    1. log(MHI/L3.6)   — gas-to-stellar balance (established)');
console.log('    2. innerDMfrac     — halo presence in inner regions');
console.log('    3. logRho0         — central DM density');
console.log('    4. logRc           — halo core radius');
console.log('    5. haloMassProxy   — total halo mass estimate');
console.log('    6. logSigma0       — baryon surface density');
console.log('    7. logRdisk        — disk scale length');
console.log('    8. baryonCompact   — baryon concentration');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
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
  const r2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
  const aic = n * Math.log(ss_res / n) + 2 * p;
  const bic = n * Math.log(ss_res / n) + p * Math.log(n);
  return { beta, yhat, r2, aic, bic, ss_res, n, p };
}

function looR2_single(xArr, yArr) {
  const n = xArr.length;
  if (n < 8) return NaN;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(yArr);
  for (let i = 0; i < n; i++) {
    const xTr = xArr.filter((_, j) => j !== i);
    const yTr = yArr.filter((_, j) => j !== i);
    const fit = olsRegress(xTr.map(v => [1, v]), yTr);
    if (!fit) return NaN;
    ss_loo += (yArr[i] - (fit.beta[0] + fit.beta[1] * xArr[i])) ** 2;
    ss_tot += (yArr[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_loo / ss_tot : 0;
}

function looR2_multi(features, yArr) {
  const n = yArr.length;
  if (n < 10) return NaN;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(yArr);
  for (let i = 0; i < n; i++) {
    const xTr = features.filter((_, j) => j !== i);
    const yTr = yArr.filter((_, j) => j !== i);
    const fit = olsRegress(xTr.map(r => [1, ...r]), yTr);
    if (!fit) return NaN;
    const pred = [1, ...features[i]].reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss_loo += (yArr[i] - pred) ** 2;
    ss_tot += (yArr[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_loo / ss_tot : 0;
}

function partialR(x, y, controls) {
  const n = x.length;
  const X_ctrl = x.map((_, i) => [1, ...controls.map(c => c[i])]);
  const fitX = olsRegress(X_ctrl, x);
  const fitY = olsRegress(X_ctrl, y);
  if (!fitX || !fitY) return NaN;
  const resX = x.map((v, i) => v - fitX.yhat[i]);
  const resY = y.map((v, i) => v - fitY.yhat[i]);
  return pearsonR(resX, resY);
}

function fitIsoHalo(dmPoints) {
  let bestRho0 = NaN, bestRc = NaN, bestChiSq = Infinity;
  for (let logRho0 = 5; logRho0 <= 9; logRho0 += 0.1) {
    for (let logRc = -1; logRc <= 2; logRc += 0.1) {
      const rho0 = Math.pow(10, logRho0);
      const rc = Math.pow(10, logRc);
      let chiSq = 0;
      for (const p of dmPoints) {
        const x = p.r / rc;
        const vdm_sq = 4 * Math.PI * 4.3009e-6 * rho0 * rc * rc * (1 - Math.atan(x) / x);
        chiSq += (p.vdm_sq - Math.max(vdm_sq, 0)) ** 2 / (Math.max(p.vdm_sq, 1) ** 2);
      }
      if (chiSq < bestChiSq) { bestChiSq = chiSq; bestRho0 = rho0; bestRc = rc; }
    }
  }
  return { rho0: bestRho0, rc: bestRc };
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  table1[name] = { L36, Rdisk, MHI, Vflat };
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
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul: vbul || 0 });
}

const galaxies = [];

for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat <= 0) continue;
  const sorted = rcPoints.filter(p => p.rad > 0 && p.vobs > 0).sort((a, b) => a.rad - b.rad);
  if (sorted.length < 8) continue;

  const logRatio = Math.log10(t1.MHI) - Math.log10(t1.L36);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
  const logRdisk = Math.log10(t1.Rdisk > 0 ? t1.Rdisk : 0.1);
  const Mstar = t1.L36 * UPSILON_DISK * 1e9;
  const Mbar = Mstar + t1.MHI * 1.33 * 1e9;
  const logMbar = Math.log10(Mbar);
  const baryonCompact = t1.Rdisk > 0 ? Math.log10(Mbar / (t1.Rdisk ** 2)) : NaN;

  if (!isFinite(logRatio)) continue;

  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);
  const outerMD = outerPts.map(p => {
    const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                   UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                   p.vgas * Math.abs(p.vgas);
    return (p.vobs ** 2) / Math.max(vBarSq, 0.01);
  }).filter(v => isFinite(v) && v > 0);
  const logOMD = Math.log10(mean(outerMD));
  if (!isFinite(logOMD)) continue;

  const dmPoints = sorted.map(p => {
    const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                   UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                   p.vgas * Math.abs(p.vgas);
    return { r: p.rad, vdm_sq: p.vobs ** 2 - Math.max(vBarSq, 0) };
  }).filter(p => p.vdm_sq > 0);

  let logRho0 = NaN, logRc = NaN;
  if (dmPoints.length >= 4) {
    const fit = fitIsoHalo(dmPoints);
    if (isFinite(fit.rho0) && isFinite(fit.rc)) {
      logRho0 = Math.log10(fit.rho0);
      logRc = Math.log10(fit.rc);
    }
  }

  let innerDMfrac = NaN;
  const innerPts = sorted.filter(p => p.rad <= t1.Rdisk);
  if (innerPts.length > 0) {
    const vals = innerPts.map(p => {
      const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                     UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                     p.vgas * Math.abs(p.vgas);
      return Math.max(p.vobs ** 2 - vBarSq, 0) / (p.vobs ** 2);
    });
    innerDMfrac = mean(vals);
  }

  const Mdyn = 2.326e5 * t1.Vflat ** 2 * sorted[sorted.length - 1].rad;
  const haloMassProxy = Math.log10(Math.max(Mdyn - Mbar, 1));

  galaxies.push({
    name, logRatio, fgas, logSigma0, logRdisk, logMbar,
    baryonCompact, logOMD, logRho0, logRc, innerDMfrac,
    haloMassProxy, Vflat: t1.Vflat,
  });
}

const allValid = galaxies.filter(g =>
  g.Vflat >= 70 &&
  isFinite(g.innerDMfrac) && isFinite(g.logRho0) && isFinite(g.logRc) &&
  isFinite(g.baryonCompact) && isFinite(g.haloMassProxy)
);

const y = allValid.map(g => g.logOMD);
const N = allValid.length;
console.log('  Complete-data sample: N=' + N + '\n');

function fmt(v) { return isFinite(v) ? v.toFixed(3) : 'N/A'; }

const candidateVars = {
  'logRatio': allValid.map(g => g.logRatio),
  'innerDMfrac': allValid.map(g => g.innerDMfrac),
  'logRho0': allValid.map(g => g.logRho0),
  'logRc': allValid.map(g => g.logRc),
  'haloMass': allValid.map(g => g.haloMassProxy),
  'logSigma0': allValid.map(g => g.logSigma0),
  'logRdisk': allValid.map(g => g.logRdisk),
  'baryonCompact': allValid.map(g => g.baryonCompact),
  'logMbar': allValid.map(g => g.logMbar),
  'fgas': allValid.map(g => g.fgas),
};

// ============================================================
// STEP 1: SINGLE-VARIABLE RANKING
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  STEP 1: SINGLE-VARIABLE RANKING');
console.log('══════════════════════════════════════════════════════════════\n');

const singleResults = [];
for (const [name, vals] of Object.entries(candidateVars)) {
  const r = pearsonR(vals, y);
  const loo = looR2_single(vals, y);
  const fit = olsRegress(vals.map(v => [1, v]), y);
  singleResults.push({ name, r, loo, aic: fit?.aic, bic: fit?.bic });
}
singleResults.sort((a, b) => b.loo - a.loo);

console.log('  Variable'.padEnd(20) + 'r'.padStart(8) + 'LOO'.padStart(8) + 'AIC'.padStart(10) + 'BIC'.padStart(10));
console.log('  ' + '-'.repeat(54));
for (const s of singleResults) {
  console.log('  ' + s.name.padEnd(20) + fmt(s.r).padStart(8) + fmt(s.loo).padStart(8) +
    fmt(s.aic).padStart(10) + fmt(s.bic).padStart(10));
}

// ============================================================
// STEP 2: EXHAUSTIVE 2-VARIABLE SEARCH
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 2: EXHAUSTIVE 2-VARIABLE MODEL SEARCH');
console.log('  (logRatio + one other variable)');
console.log('══════════════════════════════════════════════════════════════\n');

const ratio = candidateVars['logRatio'];
const twoVarResults = [];

for (const [name, vals] of Object.entries(candidateVars)) {
  if (name === 'logRatio' || name === 'fgas') continue;
  const features = allValid.map((_, i) => [ratio[i], vals[i]]);
  const fit = olsRegress(features.map(r => [1, ...r]), y);
  const loo = looR2_multi(features, y);
  if (!fit) continue;

  const pr_ratio = partialR(ratio, y, [vals]);
  const pr_other = partialR(vals, y, [ratio]);

  const r_ratio_other = pearsonR(ratio, vals);

  twoVarResults.push({
    name, r2: fit.r2, loo, aic: fit.aic, bic: fit.bic,
    beta_ratio: fit.beta[1], beta_other: fit.beta[2],
    pr_ratio, pr_other, collinearity: r_ratio_other,
  });
}
twoVarResults.sort((a, b) => b.loo - a.loo);

console.log('  Model: logRatio +'.padEnd(35) + 'R2'.padStart(7) + 'LOO'.padStart(8) + 'AIC'.padStart(10) +
  'pr(ratio)'.padStart(11) + 'pr(other)'.padStart(11) + 'collin'.padStart(8));
console.log('  ' + '-'.repeat(88));
for (const t of twoVarResults) {
  console.log('  logRatio + ' + t.name.padEnd(20) +
    fmt(t.r2).padStart(7) + fmt(t.loo).padStart(8) + fmt(t.aic).padStart(10) +
    fmt(t.pr_ratio).padStart(11) + fmt(t.pr_other).padStart(11) +
    fmt(t.collinearity).padStart(8));
}

// ============================================================
// STEP 3: EXHAUSTIVE 3-VARIABLE SEARCH
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 3: 3-VARIABLE MODEL SEARCH');
console.log('  (logRatio + halo var + structure var)');
console.log('══════════════════════════════════════════════════════════════\n');

const haloNames = ['innerDMfrac', 'logRho0', 'logRc', 'haloMass'];
const structNames = ['logSigma0', 'logRdisk', 'baryonCompact', 'logMbar'];
const threeVarResults = [];

for (const hName of haloNames) {
  for (const sName of structNames) {
    const hVals = candidateVars[hName];
    const sVals = candidateVars[sName];
    const features = allValid.map((_, i) => [ratio[i], hVals[i], sVals[i]]);
    const fit = olsRegress(features.map(r => [1, ...r]), y);
    const loo = looR2_multi(features, y);
    if (!fit) continue;

    const pr_ratio = partialR(ratio, y, [hVals, sVals]);
    const pr_halo = partialR(hVals, y, [ratio, sVals]);
    const pr_struct = partialR(sVals, y, [ratio, hVals]);

    threeVarResults.push({
      halo: hName, struct: sName, r2: fit.r2, loo, aic: fit.aic, bic: fit.bic,
      pr_ratio, pr_halo, pr_struct,
      betas: { ratio: fit.beta[1], halo: fit.beta[2], struct: fit.beta[3] },
    });
  }
}
threeVarResults.sort((a, b) => b.loo - a.loo);

console.log('  Halo var'.padEnd(16) + 'Struct var'.padEnd(16) + 'R2'.padStart(7) + 'LOO'.padStart(8) +
  'pr(rat)'.padStart(9) + 'pr(hal)'.padStart(9) + 'pr(str)'.padStart(9));
console.log('  ' + '-'.repeat(72));
for (const t of threeVarResults.slice(0, 12)) {
  console.log('  ' + t.halo.padEnd(16) + t.struct.padEnd(16) +
    fmt(t.r2).padStart(7) + fmt(t.loo).padStart(8) +
    fmt(t.pr_ratio).padStart(9) + fmt(t.pr_halo).padStart(9) + fmt(t.pr_struct).padStart(9));
}

// ============================================================
// STEP 4: THE WINNING MINIMAL MODEL
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 4: WINNING MINIMAL MODEL');
console.log('══════════════════════════════════════════════════════════════\n');

const baseline_loo = looR2_single(ratio, y);
console.log('  Baseline: logRatio alone, LOO = ' + fmt(baseline_loo));

const best2 = twoVarResults[0];
console.log('\n  Best 2-var: logRatio + ' + best2.name);
console.log('    LOO = ' + fmt(best2.loo) + ' (delta = +' + fmt(best2.loo - baseline_loo) + ')');
console.log('    AIC = ' + fmt(best2.aic) + ', BIC = ' + fmt(best2.bic));
console.log('    partial r(logRatio) = ' + fmt(best2.pr_ratio));
console.log('    partial r(' + best2.name + ') = ' + fmt(best2.pr_other));
console.log('    Collinearity r = ' + fmt(best2.collinearity));

const best3 = threeVarResults[0];
console.log('\n  Best 3-var: logRatio + ' + best3.halo + ' + ' + best3.struct);
console.log('    LOO = ' + fmt(best3.loo) + ' (delta = +' + fmt(best3.loo - baseline_loo) + ')');
console.log('    AIC = ' + fmt(best3.aic) + ', BIC = ' + fmt(best3.bic));
console.log('    partial r(logRatio) = ' + fmt(best3.pr_ratio));
console.log('    partial r(' + best3.halo + ') = ' + fmt(best3.pr_halo));
console.log('    partial r(' + best3.struct + ') = ' + fmt(best3.pr_struct));

const gain2 = best2.loo - baseline_loo;
const gain3 = best3.loo - best2.loo;
console.log('\n  Marginal gain: 1->2 var = +' + fmt(gain2) + ', 2->3 var = +' + fmt(gain3));

// ============================================================
// STEP 5: SIMPLICITY vs COMPLEXITY DECISION
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 5: SIMPLICITY vs COMPLEXITY');
console.log('══════════════════════════════════════════════════════════════\n');

const models = [
  { name: '1-var: logRatio', k: 1, loo: baseline_loo },
  { name: '2-var: logRatio + ' + best2.name, k: 2, loo: best2.loo },
  { name: '3-var: logRatio + ' + best3.halo + ' + ' + best3.struct, k: 3, loo: best3.loo },
];

console.log('  Model'.padEnd(50) + 'k'.padStart(4) + 'LOO'.padStart(8) + 'deltaLOO'.padStart(10));
console.log('  ' + '-'.repeat(70));
for (const m of models) {
  const d = m.loo - baseline_loo;
  console.log('  ' + m.name.padEnd(50) + String(m.k).padStart(4) + fmt(m.loo).padStart(8) + (d > 0 ? '+' : '') + fmt(d).padStart(9));
}

const worthIt2 = gain2 > 0.03;
const worthIt3 = gain3 > 0.02;

console.log('\n  Adding 2nd variable worth it (>0.03 LOO gain)? ' + (worthIt2 ? 'YES' : 'NO'));
console.log('  Adding 3rd variable worth it (>0.02 LOO gain)? ' + (worthIt3 ? 'YES' : 'NO'));

let recommended;
if (!worthIt2) {
  recommended = '1-var: logRatio alone';
  console.log('\n  RECOMMENDED MODEL: 1-variable (logRatio alone)');
  console.log('  Reason: No second variable adds substantial LOO improvement.');
} else if (!worthIt3) {
  recommended = '2-var: logRatio + ' + best2.name;
  console.log('\n  RECOMMENDED MODEL: 2-variable (logRatio + ' + best2.name + ')');
  console.log('  Reason: Second variable helps but third does not add enough.');
} else {
  recommended = '3-var: logRatio + ' + best3.halo + ' + ' + best3.struct;
  console.log('\n  RECOMMENDED MODEL: 3-variable');
  console.log('  Reason: Both additional variables contribute meaningfully.');
}

// ============================================================
// STEP 6: DETAILED WINNING MODEL
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 6: DETAILED WINNING MODEL COEFFICIENTS');
console.log('══════════════════════════════════════════════════════════════\n');

if (worthIt2) {
  const bestVals = candidateVars[best2.name];
  const features2 = allValid.map((_, i) => [ratio[i], bestVals[i]]);
  const fit2 = olsRegress(features2.map(r => [1, ...r]), y);
  if (fit2) {
    console.log('  logOMD = ' + fmt(fit2.beta[0]) + ' + ' + fmt(fit2.beta[1]) + ' * logRatio + ' +
      fmt(fit2.beta[2]) + ' * ' + best2.name);
    console.log('  R2 = ' + fmt(fit2.r2) + ', LOO = ' + fmt(best2.loo));

    const residuals2 = y.map((v, i) => v - fit2.yhat[i]);
    const rmse = Math.sqrt(residuals2.reduce((s, v) => s + v * v, 0) / N);
    console.log('  RMSE = ' + fmt(rmse));

    const resid1 = olsRegress(ratio.map(v => [1, v]), y);
    const rmse1 = Math.sqrt(resid1.ss_res / N);
    console.log('  (1-var RMSE = ' + fmt(rmse1) + ')');
  }
}

if (worthIt3) {
  const hVals = candidateVars[best3.halo];
  const sVals = candidateVars[best3.struct];
  const features3 = allValid.map((_, i) => [ratio[i], hVals[i], sVals[i]]);
  const fit3 = olsRegress(features3.map(r => [1, ...r]), y);
  if (fit3) {
    console.log('\n  logOMD = ' + fmt(fit3.beta[0]) + ' + ' + fmt(fit3.beta[1]) + ' * logRatio + ' +
      fmt(fit3.beta[2]) + ' * ' + best3.halo + ' + ' + fmt(fit3.beta[3]) + ' * ' + best3.struct);
    console.log('  R2 = ' + fmt(fit3.r2) + ', LOO = ' + fmt(best3.loo));
  }
}

// ============================================================
// STEP 7: PERMUTATION TEST FOR ADDED VARIABLES
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 7: PERMUTATION TEST — IS THE ADDED VARIABLE REAL?');
console.log('══════════════════════════════════════════════════════════════\n');

if (worthIt2) {
  const bestVals = candidateVars[best2.name];
  const realFeatures = allValid.map((_, i) => [ratio[i], bestVals[i]]);
  const realLOO = best2.loo;

  let countBetter = 0;
  const nPerm = 1000;
  for (let p = 0; p < nPerm; p++) {
    const shuffled = [...bestVals];
    for (let i = shuffled.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
    }
    const permFeatures = allValid.map((_, i) => [ratio[i], shuffled[i]]);
    const permLOO = looR2_multi(permFeatures, y);
    if (permLOO >= realLOO) countBetter++;
  }
  const pVal = countBetter / nPerm;
  console.log('  Added variable: ' + best2.name);
  console.log('  Real LOO = ' + fmt(realLOO) + ' vs baseline = ' + fmt(baseline_loo));
  console.log('  Permutation p-value (1000 iter) = ' + pVal.toFixed(4));
  console.log('  Significant? ' + (pVal < 0.05 ? 'YES' : 'NO'));
}

// ============================================================
// STEP 8: PHYSICAL INTERPRETATION
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 8: PHYSICAL INTERPRETATION OF THE MINIMAL MODEL');
console.log('══════════════════════════════════════════════════════════════\n');

console.log('  The minimal physical model for outer support requirement is:\n');
console.log('  RECOMMENDED: ' + recommended);

if (worthIt2) {
  console.log('\n  Physical picture:');
  console.log('    Component 1: log(MHI/L3.6) — gas-to-stellar balance');
  console.log('      What it captures: the relative state of the galaxy');
  console.log('      in terms of how much gas vs stars it has built.');
  console.log('      Gas-rich galaxies need more outer support.');
  console.log('');
  console.log('    Component 2: ' + best2.name);

  if (best2.name === 'innerDMfrac') {
    console.log('      What it captures: how dark-matter dominated the');
    console.log('      inner galaxy is. Higher inner DM fraction =');
    console.log('      halo extends deeper into galaxy center.');
  } else if (best2.name === 'logSigma0') {
    console.log('      What it captures: how concentrated the stellar');
    console.log('      component is. Higher surface density = more');
    console.log('      concentrated baryons.');
  } else if (best2.name === 'logRdisk') {
    console.log('      What it captures: physical size of the disk.');
    console.log('      Larger disks may have different halo coupling.');
  } else if (best2.name === 'logMbar') {
    console.log('      What it captures: total baryonic mass scale.');
    console.log('      More massive galaxies have different halo response.');
  } else if (best2.name === 'baryonCompact') {
    console.log('      What it captures: how concentrated the baryons');
    console.log('      are relative to disk size.');
  } else if (best2.name.startsWith('log')) {
    console.log('      What it captures: halo structural property');
    console.log('      that modulates the outer support requirement');
    console.log('      independently of gas-to-stellar balance.');
  }
}

console.log('\n  SYNTHESIS SENTENCE:');
if (worthIt2) {
  console.log('  "Outer support requirement is most tightly linked not');
  console.log('  to gas mass alone or stellar mass alone, but to the');
  console.log('  gas-to-stellar state of the galaxy, with additional');
  console.log('  but incomplete coupling to ' + best2.name + '."');
} else {
  console.log('  "Outer support requirement is most tightly linked to');
  console.log('  the gas-to-stellar state of the galaxy alone.');
  console.log('  No second variable adds substantial predictive power."');
}

// ============================================================
// STEP 9: STABILITY CHECK (BOOTSTRAP)
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  STEP 9: BOOTSTRAP STABILITY (500 resamples)');
console.log('══════════════════════════════════════════════════════════════\n');

if (worthIt2) {
  const bestVals = candidateVars[best2.name];
  const nBoot = 500;
  const betas_ratio = [], betas_other = [];

  for (let b = 0; b < nBoot; b++) {
    const indices = Array.from({ length: N }, () => Math.floor(Math.random() * N));
    const xBoot = indices.map(i => [1, ratio[i], bestVals[i]]);
    const yBoot = indices.map(i => y[i]);
    const fit = olsRegress(xBoot, yBoot);
    if (fit) {
      betas_ratio.push(fit.beta[1]);
      betas_other.push(fit.beta[2]);
    }
  }

  betas_ratio.sort((a, b) => a - b);
  betas_other.sort((a, b) => a - b);
  const lo = Math.floor(betas_ratio.length * 0.025);
  const hi = Math.floor(betas_ratio.length * 0.975);

  console.log('  beta(logRatio): ' + fmt(mean(betas_ratio)) +
    ' [' + fmt(betas_ratio[lo]) + ', ' + fmt(betas_ratio[hi]) + ']');
  console.log('  beta(' + best2.name + '): ' + fmt(mean(betas_other)) +
    ' [' + fmt(betas_other[lo]) + ', ' + fmt(betas_other[hi]) + ']');

  const ratioStable = betas_ratio[lo] * betas_ratio[hi] > 0;
  const otherStable = betas_other[lo] * betas_other[hi] > 0;
  console.log('  logRatio sign stable (95% CI excludes 0)? ' + (ratioStable ? 'YES' : 'NO'));
  console.log('  ' + best2.name + ' sign stable? ' + (otherStable ? 'YES' : 'NO'));
}

// Save results
const output = {
  phase: 121,
  title: 'Minimal Physical Model',
  N,
  recommended_model: recommended,
  single_ranking: singleResults,
  best_2var: best2,
  best_3var: { halo: best3.halo, struct: best3.struct, loo: best3.loo, r2: best3.r2 },
  baseline_loo,
  gain_1to2: gain2,
  gain_2to3: gain3,
};

const outPath = path.join(__dirname, '..', 'public', 'phase121-minimal-model.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);

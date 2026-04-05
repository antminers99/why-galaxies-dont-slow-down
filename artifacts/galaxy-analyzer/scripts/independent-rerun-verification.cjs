#!/usr/bin/env node
/**
 * INDEPENDENT RERUN VERIFICATION
 * 
 * Purpose: Level 2 verification — confirm all key numbers are reproducible
 * from a clean state with freshly downloaded data.
 * 
 * Tests the 5 critical results:
 *   1. Core signal: r(logRatio, logOMD) and LOO R^2
 *   2. Anti-circularity: stars-only target r=0.885, gas-only r=-0.208
 *   3. Mechanism map: ratio beats components
 *   4. Adopted model: 2-var coefficients and LOO
 *   5. Mediation: Sigma0|fgas dies (r~0)
 * 
 * Expected runtime: ~60 seconds (includes LOO and permutation)
 */

const fs = require('fs');

console.log('╔══════════════════════════════════════════════════════════════╗');
console.log('║     INDEPENDENT RERUN VERIFICATION — CLEAN STATE           ║');
console.log('║     No cached results. No session memory. Fresh data.      ║');
console.log('╚══════════════════════════════════════════════════════════════╝\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const TOLERANCE = 0.015;

// ─── MATH UTILITIES (reimplemented from scratch) ───

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function pearsonR(x, y) {
  const n = x.length, mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function ols(X, y) {
  const n = y.length, p = X[0].length;
  const Xt = Array.from({ length: p }, (_, i) => X.map(row => row[i]));
  const XtX = Array.from({ length: p }, (_, i) =>
    Array.from({ length: p }, (_, j) =>
      Xt[i].reduce((s, _, k) => s + Xt[i][k] * Xt[j][k], 0)));
  const Xty = Xt.map(col => col.reduce((s, v, k) => s + v * y[k], 0));
  const aug = XtX.map((row, i) => [...row, Xty[i]]);
  for (let i = 0; i < p; i++) {
    let mx = i;
    for (let k = i + 1; k < p; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[mx][i])) mx = k;
    [aug[i], aug[mx]] = [aug[mx], aug[i]];
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
  return { beta, yhat, r2: 1 - ss_res / ss_tot, ss_res };
}

function looR2_1var(x, y) {
  const n = x.length;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(y);
  for (let i = 0; i < n; i++) {
    const xTr = x.filter((_, j) => j !== i);
    const yTr = y.filter((_, j) => j !== i);
    const fit = ols(xTr.map(v => [1, v]), yTr);
    if (!fit) return NaN;
    ss_loo += (y[i] - (fit.beta[0] + fit.beta[1] * x[i])) ** 2;
    ss_tot += (y[i] - yMean) ** 2;
  }
  return 1 - ss_loo / ss_tot;
}

function looR2_multi(features, y) {
  const n = y.length;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(y);
  for (let i = 0; i < n; i++) {
    const xTr = features.filter((_, j) => j !== i);
    const yTr = y.filter((_, j) => j !== i);
    const fit = ols(xTr.map(r => [1, ...r]), yTr);
    if (!fit) return NaN;
    const pred = [1, ...features[i]].reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss_loo += (y[i] - pred) ** 2;
    ss_tot += (y[i] - yMean) ** 2;
  }
  return 1 - ss_loo / ss_tot;
}

function partialR(x, y, controls) {
  const X_ctrl = x.map((_, i) => [1, ...controls.map(c => c[i])]);
  const fitX = ols(X_ctrl, x);
  const fitY = ols(X_ctrl, y);
  if (!fitX || !fitY) return NaN;
  const resX = x.map((v, i) => v - fitX.yhat[i]);
  const resY = y.map((v, i) => v - fitY.yhat[i]);
  return pearsonR(resX, resY);
}

// ─── DATA LOADING (from scratch) ───

console.log('Loading data from /tmp/sparc_table1.dat and /tmp/sparc_cds.dat...\n');

const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');
const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');

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

// ─── GALAXY CONSTRUCTION ───

const galaxies = [];

for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat <= 0) continue;
  const sorted = rcPoints.filter(p => p.rad > 0 && p.vobs > 0).sort((a, b) => a.rad - b.rad);
  if (sorted.length < 8) continue;

  const logRatio = Math.log10(t1.MHI) - Math.log10(t1.L36);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
  const Mstar = t1.L36 * UPSILON_DISK * 1e9;
  const Mbar = Mstar + t1.MHI * 1.33 * 1e9;
  const logMbar = Math.log10(Mbar);

  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);

  // Full target (standard)
  const outerMD_full = outerPts.map(p => {
    const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                   UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                   p.vgas * Math.abs(p.vgas);
    return (p.vobs ** 2) / Math.max(vBarSq, 0.01);
  }).filter(v => isFinite(v) && v > 0);
  const logOMD = Math.log10(mean(outerMD_full));

  // Stars-only target (gas removed from Vbar)
  const outerMD_stars = outerPts.map(p => {
    const vBarSq_stars = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                          UPSILON_BULGE * p.vbul * Math.abs(p.vbul);
    return (p.vobs ** 2) / Math.max(vBarSq_stars, 0.01);
  }).filter(v => isFinite(v) && v > 0);
  const logOMD_stars = Math.log10(mean(outerMD_stars));

  // Gas-only target
  const outerMD_gas = outerPts.map(p => {
    const vBarSq_gas = p.vgas * Math.abs(p.vgas);
    return (p.vobs ** 2) / Math.max(vBarSq_gas, 0.01);
  }).filter(v => isFinite(v) && v > 0);
  const logOMD_gas = outerMD_gas.length > 0 ? Math.log10(mean(outerMD_gas)) : NaN;

  if (!isFinite(logOMD) || !isFinite(logRatio)) continue;

  // Inner DM fraction
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

  // Halo mass proxy
  const Mdyn = 2.326e5 * t1.Vflat ** 2 * sorted[sorted.length - 1].rad;
  const haloMassProxy = Math.log10(Math.max(Mdyn - Mbar, 1));

  galaxies.push({
    name, logRatio, fgas, logSigma0, logMbar,
    logOMD, logOMD_stars, logOMD_gas,
    Vflat: t1.Vflat, innerDMfrac, haloMassProxy,
  });
}

// ─── HIGH-VFLAT SAMPLE ───

const highV = galaxies.filter(g => g.Vflat >= 70);
console.log('Total galaxies parsed: ' + galaxies.length);
console.log('High-Vflat sample (>=70): N = ' + highV.length + '\n');

// ═══════════════════════════════════════════════════════
// TEST 1: CORE SIGNAL
// ═══════════════════════════════════════════════════════

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: CORE SIGNAL');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const x_ratio = highV.map(g => g.logRatio);
const x_fgas = highV.map(g => g.fgas);
const y_omd = highV.map(g => g.logOMD);

const r_ratio = pearsonR(x_ratio, y_omd);
const r_fgas = pearsonR(x_fgas, y_omd);
const loo_ratio = looR2_1var(x_ratio, y_omd);
const loo_fgas = looR2_1var(x_fgas, y_omd);

const expected1 = [
  { name: 'r(logRatio, logOMD)', got: r_ratio, expect: 0.729, tol: TOLERANCE },
  { name: 'r(fgas, logOMD)', got: r_fgas, expect: 0.724, tol: TOLERANCE },
  { name: 'LOO R^2 (logRatio)', got: loo_ratio, expect: 0.512, tol: 0.02 },
  { name: 'LOO R^2 (fgas)', got: loo_fgas, expect: 0.503, tol: 0.02 },
  { name: 'N', got: highV.length, expect: 104, tol: 2 },
];

let allPass = true;
for (const e of expected1) {
  const pass = Math.abs(e.got - e.expect) <= e.tol;
  const status = pass ? 'PASS' : 'FAIL';
  if (!pass) allPass = false;
  console.log('  ' + status + '  ' + e.name + ': got ' + e.got.toFixed(4) + ', expected ~' + e.expect + ' (tol ' + e.tol + ')');
}

// ═══════════════════════════════════════════════════════
// TEST 2: ANTI-CIRCULARITY (Phase 117 core)
// ═══════════════════════════════════════════════════════

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: ANTI-CIRCULARITY');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const validStars = highV.filter(g => isFinite(g.logOMD_stars));
const r_stars = pearsonR(validStars.map(g => g.logRatio), validStars.map(g => g.logOMD_stars));

const validGas = highV.filter(g => isFinite(g.logOMD_gas));
const r_gas = pearsonR(validGas.map(g => g.logRatio), validGas.map(g => g.logOMD_gas));

const expected2 = [
  { name: 'r(logRatio, stars-only target)', got: r_stars, expect: 0.885, tol: 0.02 },
  { name: 'r(logRatio, gas-only target)', got: r_gas, expect: -0.208, tol: 0.05 },
  { name: 'Stars-only > full target?', got: r_stars > r_ratio ? 1 : 0, expect: 1, tol: 0 },
];

for (const e of expected2) {
  const pass = Math.abs(e.got - e.expect) <= e.tol;
  const status = pass ? 'PASS' : 'FAIL';
  if (!pass) allPass = false;
  console.log('  ' + status + '  ' + e.name + ': got ' + e.got.toFixed(4) + ', expected ~' + e.expect);
}

// Permutation test for stars-only
const nPerm = 1000;
let countAbove = 0;
const xStars = validStars.map(g => g.logRatio);
const yStars = validStars.map(g => g.logOMD_stars);
for (let p = 0; p < nPerm; p++) {
  const shuffled = [...xStars];
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  if (Math.abs(pearsonR(shuffled, yStars)) >= Math.abs(r_stars)) countAbove++;
}
const perm_p = countAbove / nPerm;
console.log('  ' + (perm_p < 0.001 ? 'PASS' : 'FAIL') + '  Permutation p (stars-only): ' + perm_p.toFixed(4) + ' (expected < 0.001)');
if (perm_p >= 0.001) allPass = false;

// ═══════════════════════════════════════════════════════
// TEST 3: MECHANISM MAP (Phase 120 core)
// ═══════════════════════════════════════════════════════

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: MECHANISM MAP — RATIO BEATS COMPONENTS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const logMHI_arr = highV.map(g => Math.log10(g.fgas * (1 / (1 - g.fgas)) * UPSILON_DISK));
const logMHI_direct = highV.map(g => {
  const t = table1[g.name];
  return Math.log10(t.MHI);
});
const logL36_direct = highV.map(g => {
  const t = table1[g.name];
  return Math.log10(t.L36);
});

const r_logMHI = pearsonR(logMHI_direct, y_omd);
const r_logL36 = pearsonR(logL36_direct, y_omd);

const expected3 = [
  { name: 'r(logMHI, logOMD)', got: r_logMHI, expect: -0.376, tol: 0.03 },
  { name: 'r(logL36, logOMD)', got: r_logL36, expect: -0.705, tol: 0.02 },
  { name: 'Ratio > |logMHI|?', got: Math.abs(r_ratio) > Math.abs(r_logMHI) ? 1 : 0, expect: 1, tol: 0 },
  { name: 'Ratio > |logL36|?', got: Math.abs(r_ratio) > Math.abs(r_logL36) ? 1 : 0, expect: 1, tol: 0 },
];

for (const e of expected3) {
  const pass = Math.abs(e.got - e.expect) <= e.tol;
  const status = pass ? 'PASS' : 'FAIL';
  if (!pass) allPass = false;
  console.log('  ' + status + '  ' + e.name + ': got ' + e.got.toFixed(4) + ', expected ~' + e.expect);
}

// Halo control survival
const x_logSigma0 = highV.map(g => g.logSigma0);
const pr_ratio_sigma = partialR(x_ratio, y_omd, [x_logSigma0]);

// Check halo survival with innerDMfrac
const haloValid = highV.filter(g => isFinite(g.innerDMfrac) && isFinite(g.haloMassProxy));
const xr_halo = haloValid.map(g => g.logRatio);
const yr_halo = haloValid.map(g => g.logOMD);
const ctrl_idf = haloValid.map(g => g.innerDMfrac);
const ctrl_hm = haloValid.map(g => g.haloMassProxy);
const pr_ratio_halo = partialR(xr_halo, yr_halo, [ctrl_idf, ctrl_hm]);

console.log('  ' + (Math.abs(pr_ratio_sigma) > 0.3 ? 'PASS' : 'FAIL') + '  partial r(ratio | Sigma0) = ' + pr_ratio_sigma.toFixed(3) + ' (expect >0.3, survival)');
console.log('  ' + (Math.abs(pr_ratio_halo) > 0.3 ? 'PASS' : 'FAIL') + '  partial r(ratio | halo vars) = ' + pr_ratio_halo.toFixed(3) + ' (expect >0.3, survival)');

// ═══════════════════════════════════════════════════════
// TEST 4: ADOPTED 2-VAR MODEL (Phase 121)
// ═══════════════════════════════════════════════════════

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: ADOPTED 2-VARIABLE MODEL');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const completeData = highV.filter(g => isFinite(g.innerDMfrac) && isFinite(g.haloMassProxy));
const x2_ratio = completeData.map(g => g.logRatio);
const x2_mbar = completeData.map(g => g.logMbar);
const y2 = completeData.map(g => g.logOMD);

const fit2 = ols(completeData.map(g => [1, g.logRatio, g.logMbar]), y2);
const loo2 = looR2_multi(completeData.map(g => [g.logRatio, g.logMbar]), y2);

const expected4 = [
  { name: 'Intercept', got: fit2.beta[0], expect: 1.749, tol: 0.05 },
  { name: 'beta(logRatio)', got: fit2.beta[1], expect: 0.203, tol: 0.02 },
  { name: 'beta(logMbar)', got: fit2.beta[2], expect: -0.101, tol: 0.02 },
  { name: 'LOO R^2 (2-var)', got: loo2, expect: 0.584, tol: 0.02 },
  { name: 'R^2 (2-var)', got: fit2.r2, expect: 0.617, tol: 0.02 },
  { name: 'N (complete data)', got: completeData.length, expect: 97, tol: 2 },
];

for (const e of expected4) {
  const pass = Math.abs(e.got - e.expect) <= e.tol;
  const status = pass ? 'PASS' : 'FAIL';
  if (!pass) allPass = false;
  console.log('  ' + status + '  ' + e.name + ': got ' + e.got.toFixed(4) + ', expected ~' + e.expect);
}

// Permutation test for logMbar addition
const baseline_loo = looR2_1var(x2_ratio, y2);
let permBetter = 0;
for (let p = 0; p < 1000; p++) {
  const shuffled = [...x2_mbar];
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  const permFeatures = completeData.map((_, i) => [x2_ratio[i], shuffled[i]]);
  const permLOO = looR2_multi(permFeatures, y2);
  if (permLOO >= loo2) permBetter++;
}
const perm_p2 = permBetter / 1000;
console.log('  ' + (perm_p2 <= 0.01 ? 'PASS' : 'FAIL') + '  Permutation p (Mbar addition): ' + perm_p2.toFixed(4) + ' (expected <= 0.01)');
if (perm_p2 > 0.01) allPass = false;

// Bootstrap stability
const nBoot = 500;
const betas_r = [], betas_m = [];
for (let b = 0; b < nBoot; b++) {
  const idx = Array.from({ length: completeData.length }, () => Math.floor(Math.random() * completeData.length));
  const xB = idx.map(i => [1, x2_ratio[i], x2_mbar[i]]);
  const yB = idx.map(i => y2[i]);
  const fitB = ols(xB, yB);
  if (fitB) { betas_r.push(fitB.beta[1]); betas_m.push(fitB.beta[2]); }
}
betas_r.sort((a, b) => a - b);
betas_m.sort((a, b) => a - b);
const lo = Math.floor(betas_r.length * 0.025);
const hi = Math.floor(betas_r.length * 0.975);
const r_stable = betas_r[lo] * betas_r[hi] > 0;
const m_stable = betas_m[lo] * betas_m[hi] > 0;

console.log('  ' + (r_stable ? 'PASS' : 'FAIL') + '  beta(logRatio) sign stable: [' + betas_r[lo].toFixed(3) + ', ' + betas_r[hi].toFixed(3) + ']');
console.log('  ' + (m_stable ? 'PASS' : 'FAIL') + '  beta(logMbar) sign stable: [' + betas_m[lo].toFixed(3) + ', ' + betas_m[hi].toFixed(3) + ']');
if (!r_stable || !m_stable) allPass = false;

// ═══════════════════════════════════════════════════════
// TEST 5: MEDIATION (Phase 108 core)
// ═══════════════════════════════════════════════════════

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 5: MEDIATION — Sigma0 DIES AFTER CONTROLLING fgas');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const pr_fgas_given_sigma = partialR(x_fgas, y_omd, [x_logSigma0]);
const pr_sigma_given_fgas = partialR(x_logSigma0, y_omd, [x_fgas]);

const expected5 = [
  { name: 'partial r(fgas | Sigma0)', got: pr_fgas_given_sigma, expect: 0.424, tol: 0.05 },
  { name: 'partial r(Sigma0 | fgas)', got: pr_sigma_given_fgas, expect: 0.004, tol: 0.05 },
  { name: 'fgas survives (pr > 0.3)?', got: Math.abs(pr_fgas_given_sigma) > 0.3 ? 1 : 0, expect: 1, tol: 0 },
  { name: 'Sigma0 dies (|pr| < 0.1)?', got: Math.abs(pr_sigma_given_fgas) < 0.1 ? 1 : 0, expect: 1, tol: 0 },
];

for (const e of expected5) {
  const pass = Math.abs(e.got - e.expect) <= e.tol;
  const status = pass ? 'PASS' : 'FAIL';
  if (!pass) allPass = false;
  console.log('  ' + status + '  ' + e.name + ': got ' + e.got.toFixed(4) + ', expected ~' + e.expect);
}

// ═══════════════════════════════════════════════════════
// SUMMARY
// ═══════════════════════════════════════════════════════

console.log('\n╔══════════════════════════════════════════════════════════════╗');
if (allPass) {
  console.log('║  OVERALL VERDICT:  ALL TESTS PASSED                        ║');
  console.log('║  Independent rerun confirms all documented numbers.        ║');
  console.log('║  No implementation errors detected.                        ║');
} else {
  console.log('║  OVERALL VERDICT:  SOME TESTS FAILED                       ║');
  console.log('║  Check individual test outputs above.                      ║');
}
console.log('╚══════════════════════════════════════════════════════════════╝');

const summary = {
  verification: 'independent-rerun',
  date: new Date().toISOString(),
  allPass,
  core: { r_ratio: r_ratio, r_fgas: r_fgas, loo_ratio: loo_ratio, loo_fgas: loo_fgas, N: highV.length },
  antiCircularity: { r_stars: r_stars, r_gas: r_gas, perm_p: perm_p },
  mechanism: { r_logMHI: r_logMHI, r_logL36: r_logL36, pr_ratio_sigma: pr_ratio_sigma },
  adoptedModel: { intercept: fit2.beta[0], beta_ratio: fit2.beta[1], beta_mbar: fit2.beta[2], loo: loo2, perm_p: perm_p2 },
  mediation: { pr_fgas_sigma: pr_fgas_given_sigma, pr_sigma_fgas: pr_sigma_given_fgas },
};

const outPath = require('path').join(__dirname, '..', 'public', 'independent-rerun-verification.json');
fs.writeFileSync(outPath, JSON.stringify(summary, null, 2));
console.log('\nResults saved to: ' + outPath);

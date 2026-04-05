#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 120: MECHANISM MAP');
console.log('');
console.log('  WHY does log(MHI/L3.6) predict outer support requirement?');
console.log('  Is it gas alone? Stars alone? The ratio itself?');
console.log('  Or a proxy for structure / halo / evolutionary state?');
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
  const resid = y.map((v, i) => v - yhat[i]);
  const ss_res = resid.reduce((s, v) => s + v * v, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  const r2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
  const k = p - 1;
  const aic = n * Math.log(ss_res / n) + 2 * p;
  const bic = n * Math.log(ss_res / n) + p * Math.log(n);
  return { beta, yhat, resid, r2, aic, bic };
}

function residualize(x, controls) {
  const n = x.length;
  const X = x.map((_, i) => [1, ...controls.map(c => c[i])]);
  const fit = olsRegress(X, x);
  if (!fit) return x;
  return x.map((v, i) => v - fit.yhat[i]);
}

function partialR(x, y, controls) {
  const resX = residualize(x, controls);
  const resY = residualize(y, controls);
  return pearsonR(resX, resY);
}

function looR2(xArr, yArr) {
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

function looR2_multi(XArr, yArr) {
  const n = yArr.length;
  if (n < 10) return NaN;
  const p = XArr[0].length;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(yArr);
  for (let i = 0; i < n; i++) {
    const xTr = XArr.filter((_, j) => j !== i);
    const yTr = yArr.filter((_, j) => j !== i);
    const fit = olsRegress(xTr.map(r => [1, ...r]), yTr);
    if (!fit) return NaN;
    const pred = [1, ...XArr[i]].reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss_loo += (yArr[i] - pred) ** 2;
    ss_tot += (yArr[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_loo / ss_tot : 0;
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

  const logMHI = Math.log10(t1.MHI);
  const logL36 = Math.log10(t1.L36);
  const logRatio = logMHI - logL36;
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
  const logRdisk = Math.log10(t1.Rdisk > 0 ? t1.Rdisk : 0.1);
  const Mstar = t1.L36 * UPSILON_DISK * 1e9;
  const Mbar = Mstar + t1.MHI * 1.33 * 1e9;
  const logMbar = Math.log10(Mbar);
  const SFE = t1.L36 / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logSFE = Math.log10(SFE > 0 ? SFE : 1e-10);
  const baryonCompact = t1.Rdisk > 0 ? Math.log10(Mbar / (t1.Rdisk ** 2)) : NaN;

  if (!isFinite(logRatio) || !isFinite(fgas)) continue;

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
    name, logMHI, logL36, logRatio, fgas, logSigma0, logRdisk,
    logMbar, logSFE, baryonCompact, logOMD, logRho0, logRc,
    innerDMfrac, haloMassProxy, Vflat: t1.Vflat,
  });
}

const highV = galaxies.filter(g => g.Vflat >= 70);
const N = highV.length;
console.log('  High-Vflat sample: N=' + N + '\n');

const y = highV.map(g => g.logOMD);

function fmt(v) { return isFinite(v) ? v.toFixed(3) : 'N/A'; }

// ============================================================
// TEST 1: COMPONENT COMPETITION
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: COMPONENT COMPETITION');
console.log('  log(MHI/L3.6) vs log(MHI) vs log(L3.6) vs both');
console.log('══════════════════════════════════════════════════════════════\n');

const predictors1 = [
  { name: 'log(MHI/L3.6)', vals: highV.map(g => g.logRatio) },
  { name: 'log(MHI)', vals: highV.map(g => g.logMHI) },
  { name: 'log(L3.6)', vals: highV.map(g => g.logL36) },
  { name: 'fgas', vals: highV.map(g => g.fgas) },
];

const table1Results = [];

for (const pred of predictors1) {
  const r = pearsonR(pred.vals, y);
  const loo = looR2(pred.vals, y);
  const fit = olsRegress(pred.vals.map(v => [1, v]), y);
  table1Results.push({ name: pred.name, r, loo, aic: fit?.aic, bic: fit?.bic });
  console.log('  ' + pred.name.padEnd(18) + 'r=' + fmt(r) + '  LOO=' + fmt(loo) +
    '  AIC=' + fmt(fit?.aic) + '  BIC=' + fmt(fit?.bic));
}

const twoVar = highV.map(g => [g.logMHI, g.logL36]);
const fit2 = olsRegress(twoVar.map(r => [1, ...r]), y);
const loo2 = looR2_multi(twoVar, y);
table1Results.push({ name: 'log(MHI)+log(L3.6)', r: Math.sqrt(Math.max(fit2?.r2 || 0, 0)), loo: loo2, aic: fit2?.aic, bic: fit2?.bic });
console.log('  ' + 'log(MHI)+log(L3.6)'.padEnd(18) + 'R=' + fmt(Math.sqrt(Math.max(fit2?.r2 || 0, 0))) +
  '  LOO=' + fmt(loo2) + '  AIC=' + fmt(fit2?.aic) + '  BIC=' + fmt(fit2?.bic));

if (fit2) {
  console.log('\n  Two-var model coefficients:');
  console.log('    intercept = ' + fmt(fit2.beta[0]));
  console.log('    beta(logMHI) = ' + fmt(fit2.beta[1]));
  console.log('    beta(logL36) = ' + fmt(fit2.beta[2]));
  const ratioLike = Math.abs(fit2.beta[1] + fit2.beta[2]) < 0.3 * Math.max(Math.abs(fit2.beta[1]), Math.abs(fit2.beta[2]));
  console.log('    beta(MHI) + beta(L36) = ' + fmt(fit2.beta[1] + fit2.beta[2]));
  console.log('    Ratio-like (betas near equal-opposite)? ' + (ratioLike ? 'YES' : 'NO'));
}

const pr_mhi = partialR(highV.map(g => g.logMHI), y, [highV.map(g => g.logL36)]);
const pr_l36 = partialR(highV.map(g => g.logL36), y, [highV.map(g => g.logMHI)]);
console.log('\n  Partial correlations:');
console.log('    partial r(logMHI, logOMD | logL36) = ' + fmt(pr_mhi));
console.log('    partial r(logL36, logOMD | logMHI) = ' + fmt(pr_l36));

const ratioWins = table1Results[0].loo >= table1Results[1].loo && table1Results[0].loo >= table1Results[2].loo;
console.log('\n  VERDICT: Ratio beats both components? ' + (ratioWins ? 'YES' : 'NO'));
const twoVarBeats = loo2 > table1Results[0].loo + 0.02;
console.log('  Two-var substantially beats ratio (+>0.02 LOO)? ' + (twoVarBeats ? 'YES — ratio is lossy' : 'NO — ratio captures the information'));

console.log('');

// ============================================================
// TEST 2: RESIDUAL DISSECTION — "More gas" or "Less stars"?
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: RESIDUAL DISSECTION');
console.log('  Is it "more gas" or "less stars" or the ratio itself?');
console.log('══════════════════════════════════════════════════════════════\n');

const logMHI_perp = residualize(highV.map(g => g.logMHI), [highV.map(g => g.logL36)]);
const logL36_perp = residualize(highV.map(g => g.logL36), [highV.map(g => g.logMHI)]);
const ratio_perp_mhi = residualize(highV.map(g => g.logRatio), [highV.map(g => g.logMHI)]);
const ratio_perp_l36 = residualize(highV.map(g => g.logRatio), [highV.map(g => g.logL36)]);

const r_mhi_perp = pearsonR(logMHI_perp, y);
const r_l36_perp = pearsonR(logL36_perp, y);
const r_ratio = pearsonR(highV.map(g => g.logRatio), y);
const loo_mhi_perp = looR2(logMHI_perp, y);
const loo_l36_perp = looR2(logL36_perp, y);

console.log('  MHI_perp (gas after removing star effect):');
console.log('    r = ' + fmt(r_mhi_perp) + ', LOO = ' + fmt(loo_mhi_perp));
console.log('  L36_perp (stars after removing gas effect):');
console.log('    r = ' + fmt(r_l36_perp) + ', LOO = ' + fmt(loo_l36_perp));
console.log('  log(MHI/L3.6) (raw ratio):');
console.log('    r = ' + fmt(r_ratio) + ', LOO = ' + fmt(looR2(highV.map(g => g.logRatio), y)));

const pr_ratio_both = partialR(highV.map(g => g.logRatio), y,
  [highV.map(g => g.logMHI), highV.map(g => g.logL36)]);
console.log('\n  partial r(ratio, logOMD | logMHI + logL36) = ' + fmt(pr_ratio_both));

console.log('\n  Which residual is stronger? ' +
  (Math.abs(r_mhi_perp) > Math.abs(r_l36_perp) ? 'MHI_perp (gas side drives more)' : 'L36_perp (stellar side drives more)'));

const ratioBeatsBoth = Math.abs(r_ratio) > Math.abs(r_mhi_perp) && Math.abs(r_ratio) > Math.abs(r_l36_perp);
console.log('  Ratio beats both residuals? ' + (ratioBeatsBoth ? 'YES — ratio state is fundamental' : 'NO'));

console.log('');

// ============================================================
// TEST 3: STRUCTURE CONTROLS
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: STRUCTURE CONTROLS');
console.log('  Does log(MHI/L3.6) survive controlling for size,');
console.log('  surface density, and baryonic compactness?');
console.log('══════════════════════════════════════════════════════════════\n');

const validBC = highV.filter(g => isFinite(g.baryonCompact));

const structControls = [
  { name: 'Rdisk (size)', ctrl: [highV.map(g => g.logRdisk)] },
  { name: 'Sigma0 (surface density)', ctrl: [highV.map(g => g.logSigma0)] },
  { name: 'BaryonCompact', ctrl: [validBC.map(g => g.baryonCompact)], subset: validBC },
  { name: 'Rdisk + Sigma0', ctrl: [highV.map(g => g.logRdisk), highV.map(g => g.logSigma0)] },
  { name: 'ALL structure', ctrl: [highV.map(g => g.logRdisk), highV.map(g => g.logSigma0)] },
];

for (const sc of structControls) {
  const sub = sc.subset || highV;
  const pr = partialR(sub.map(g => g.logRatio), sub.map(g => g.logOMD), sc.ctrl);
  const origR = pearsonR(sub.map(g => g.logRatio), sub.map(g => g.logOMD));
  const pct = (Math.abs(pr) / Math.abs(origR) * 100).toFixed(0);
  console.log('  partial r(logRatio, logOMD | ' + sc.name + ') = ' + fmt(pr) + '  (' + pct + '% of original r=' + fmt(origR) + ')');
}

const pr_sigma_ratio = partialR(highV.map(g => g.logSigma0), y, [highV.map(g => g.logRatio)]);
console.log('\n  Reverse: partial r(Sigma0, logOMD | logRatio) = ' + fmt(pr_sigma_ratio));
const pr_rdisk_ratio = partialR(highV.map(g => g.logRdisk), y, [highV.map(g => g.logRatio)]);
console.log('  Reverse: partial r(Rdisk, logOMD | logRatio) = ' + fmt(pr_rdisk_ratio));

console.log('');

// ============================================================
// TEST 4: TOTAL BARYONIC MASS CONTROL
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 4: TOTAL BARYONIC MASS CONTROL');
console.log('  Is the signal just because less massive = more gas-rich?');
console.log('══════════════════════════════════════════════════════════════\n');

const pr_ratio_mbar = partialR(highV.map(g => g.logRatio), y, [highV.map(g => g.logMbar)]);
const pr_mbar_ratio = partialR(highV.map(g => g.logMbar), y, [highV.map(g => g.logRatio)]);
const r_mbar = pearsonR(highV.map(g => g.logMbar), y);
const loo_mbar = looR2(highV.map(g => g.logMbar), y);

console.log('  log(Mbar) alone: r = ' + fmt(r_mbar) + ', LOO = ' + fmt(loo_mbar));
console.log('  partial r(logRatio, logOMD | logMbar) = ' + fmt(pr_ratio_mbar));
console.log('  partial r(logMbar, logOMD | logRatio) = ' + fmt(pr_mbar_ratio));

const pr_ratio_mbar_vflat = partialR(highV.map(g => g.logRatio), y,
  [highV.map(g => g.logMbar), highV.map(g => g.Vflat)]);
console.log('  partial r(logRatio, logOMD | logMbar + Vflat) = ' + fmt(pr_ratio_mbar_vflat));

const ratioSurv = Math.abs(pr_ratio_mbar) > 0.3;
console.log('\n  Ratio survives mass control? ' + (ratioSurv ? 'YES' : 'NO'));

console.log('');

// ============================================================
// TEST 5: SFE / EVOLUTIONARY STATE
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 5: STAR FORMATION EFFICIENCY / EVOLUTIONARY STATE');
console.log('  Is log(MHI/L3.6) just capturing evolutionary stage?');
console.log('══════════════════════════════════════════════════════════════\n');

const r_sfe = pearsonR(highV.map(g => g.logSFE), y);
const loo_sfe = looR2(highV.map(g => g.logSFE), y);
console.log('  log(SFE) alone: r = ' + fmt(r_sfe) + ', LOO = ' + fmt(loo_sfe));
console.log('  log(MHI/L3.6):  r = ' + fmt(r_ratio) + ', LOO = ' + fmt(looR2(highV.map(g => g.logRatio), y)));

const r_ratio_sfe = pearsonR(highV.map(g => g.logRatio), highV.map(g => g.logSFE));
console.log('  Correlation ratio vs SFE: r = ' + fmt(r_ratio_sfe));

const pr_ratio_sfe = partialR(highV.map(g => g.logRatio), y, [highV.map(g => g.logSFE)]);
const pr_sfe_ratio = partialR(highV.map(g => g.logSFE), y, [highV.map(g => g.logRatio)]);
console.log('  partial r(logRatio, logOMD | logSFE) = ' + fmt(pr_ratio_sfe));
console.log('  partial r(logSFE, logOMD | logRatio) = ' + fmt(pr_sfe_ratio));

const stellarFrac = highV.map(g => 1 - g.fgas);
const r_sfrac = pearsonR(stellarFrac, y);
const pr_ratio_sfrac = partialR(highV.map(g => g.logRatio), y, [stellarFrac]);
console.log('\n  Stellar fraction alone: r = ' + fmt(r_sfrac));
console.log('  partial r(logRatio, logOMD | stellarFrac) = ' + fmt(pr_ratio_sfrac));

console.log('');

// ============================================================
// TEST 6: HALO MECHANISM CROSS-TEST
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 6: HALO MECHANISM CROSS-TEST');
console.log('  Does ratio survive controlling halo? Does halo survive');
console.log('  controlling ratio? Both survive = multi-component.');
console.log('══════════════════════════════════════════════════════════════\n');

const validH = highV.filter(g => isFinite(g.innerDMfrac) && isFinite(g.logRho0) && isFinite(g.logRc));
const NH = validH.length;
console.log('  N with halo fits: ' + NH);

if (NH >= 20) {
  const haloVars = [
    { name: 'innerDMfrac', vals: validH.map(g => g.innerDMfrac) },
    { name: 'logRho0', vals: validH.map(g => g.logRho0) },
    { name: 'logRc', vals: validH.map(g => g.logRc) },
    { name: 'haloMassProxy', vals: validH.map(g => g.haloMassProxy) },
  ];

  const yH = validH.map(g => g.logOMD);
  const ratioH = validH.map(g => g.logRatio);
  const r_ratio_h = pearsonR(ratioH, yH);

  console.log('\n  Individual halo predictors of logOMD:');
  for (const hv of haloVars) {
    const r = pearsonR(hv.vals, yH);
    const loo = looR2(hv.vals, yH);
    console.log('    ' + hv.name.padEnd(16) + 'r=' + fmt(r) + '  LOO=' + fmt(loo));
  }

  console.log('\n  Ratio survives each halo control:');
  for (const hv of haloVars) {
    const pr = partialR(ratioH, yH, [hv.vals]);
    const pct = (Math.abs(pr) / Math.abs(r_ratio_h) * 100).toFixed(0);
    console.log('    partial r(ratio | ' + hv.name.padEnd(14) + ') = ' + fmt(pr) + '  (' + pct + '%)');
  }

  const allHaloCtrl = haloVars.map(h => h.vals);
  const pr_ratio_allH = partialR(ratioH, yH, allHaloCtrl);
  const pct_all = (Math.abs(pr_ratio_allH) / Math.abs(r_ratio_h) * 100).toFixed(0);
  console.log('    partial r(ratio | ALL halo)    = ' + fmt(pr_ratio_allH) + '  (' + pct_all + '%)');

  console.log('\n  Reverse: halo survives ratio control:');
  for (const hv of haloVars) {
    const pr = partialR(hv.vals, yH, [ratioH]);
    console.log('    partial r(' + hv.name.padEnd(14) + ' | ratio) = ' + fmt(pr));
  }

  const bothSurvive = Math.abs(pr_ratio_allH) > 0.3 &&
    haloVars.some(hv => Math.abs(partialR(hv.vals, yH, [ratioH])) > 0.2);
  const ratioOnly = Math.abs(pr_ratio_allH) > 0.3 &&
    haloVars.every(hv => Math.abs(partialR(hv.vals, yH, [ratioH])) < 0.15);

  console.log('\n  Both survive? ' + (bothSurvive ? 'YES — multi-component picture' : 'NO'));
  console.log('  Ratio alone survives? ' + (ratioOnly ? 'YES — ratio deeper than halo' : 'NO'));
}

console.log('');

// ============================================================
// TABLE 3: MECHANISM VERDICT
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  MECHANISM VERDICT TABLE');
console.log('══════════════════════════════════════════════════════════════\n');

const r_ratio_val = pearsonR(highV.map(g => g.logRatio), y);
const r_mhi_val = pearsonR(highV.map(g => g.logMHI), y);
const r_l36_val = pearsonR(highV.map(g => g.logL36), y);

const gasOnlyVerdict = Math.abs(r_mhi_val) >= Math.abs(r_ratio_val) - 0.02 ? 'PASS' : 'FAIL';
const starsOnlyVerdict = Math.abs(r_l36_val) >= Math.abs(r_ratio_val) - 0.02 ? 'PASS' : 'FAIL';

const ratioBoth2 = Math.abs(r_ratio_val) > Math.abs(r_mhi_val) + 0.02 &&
                   Math.abs(r_ratio_val) > Math.abs(r_l36_val) + 0.02;
const ratioStateVerdict = ratioBoth2 ? 'PASS' : (Math.abs(r_ratio_val) >= Math.abs(r_mhi_val) ? 'PARTIAL' : 'FAIL');

const pr_struct = partialR(highV.map(g => g.logRatio), y,
  [highV.map(g => g.logRdisk), highV.map(g => g.logSigma0)]);
const structOnlyVerdict = Math.abs(pr_struct) < 0.2 ? 'PASS' : 'FAIL';

const validHM = highV.filter(g => isFinite(g.innerDMfrac) && isFinite(g.logRho0));
let haloMedVerdict = 'N/A';
if (validHM.length >= 20) {
  const prH = partialR(validHM.map(g => g.logRatio), validHM.map(g => g.logOMD),
    [validHM.map(g => g.innerDMfrac), validHM.map(g => g.logRho0)]);
  haloMedVerdict = Math.abs(prH) < 0.2 ? 'PASS' : (Math.abs(prH) < 0.4 ? 'PARTIAL' : 'FAIL');
}

console.log('  Hypothesis             | Verdict  | Evidence');
console.log('  -----------------------|----------|------------------------------------------');
console.log('  Gas-only driver        | ' + gasOnlyVerdict.padEnd(8) + ' | logMHI r=' + fmt(r_mhi_val) + ' vs ratio r=' + fmt(r_ratio_val));
console.log('  Stars-only driver      | ' + starsOnlyVerdict.padEnd(8) + ' | logL36 r=' + fmt(r_l36_val) + ' vs ratio r=' + fmt(r_ratio_val));
console.log('  Ratio-state (balance)  | ' + ratioStateVerdict.padEnd(8) + ' | ratio beats both components');
console.log('  Structure-only         | ' + structOnlyVerdict.padEnd(8) + ' | partial r(ratio|struct) = ' + fmt(pr_struct));
console.log('  Halo-mediated          | ' + haloMedVerdict.padEnd(8) + ' | ratio survives halo controls');

console.log('');

// ============================================================
// FINAL SYNTHESIS
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  FINAL SYNTHESIS');
console.log('══════════════════════════════════════════════════════════════\n');

const canSay = ratioStateVerdict !== 'FAIL' &&
  structOnlyVerdict === 'FAIL' &&
  haloMedVerdict === 'FAIL' &&
  Math.abs(pr_ratio_mbar) > 0.3;

if (canSay) {
  console.log('  TARGET SENTENCE ACHIEVED:');
  console.log('  "log(MHI/L3.6) is not merely a proxy for gas mass alone,');
  console.log('   nor stellar mass alone, nor surface structure alone,');
  console.log('   but represents an independent relative physical state');
  console.log('   of the galaxy."');
} else {
  console.log('  TARGET SENTENCE STATUS:');
  console.log('  Ratio-state: ' + ratioStateVerdict);
  console.log('  Structure-only rejected: ' + structOnlyVerdict);
  console.log('  Halo-mediated rejected: ' + haloMedVerdict);
  console.log('  Survives mass control: ' + (Math.abs(pr_ratio_mbar) > 0.3));
  if (ratioStateVerdict === 'FAIL') {
    console.log('\n  CAUTION: Ratio does not clearly beat components.');
    console.log('  Consider: "compressed proxy of galaxy state" rather');
    console.log('  than "independent ratio state variable".');
  }
}

console.log('\n  KEY NUMBERS SUMMARY:');
console.log('    Ratio r=' + fmt(r_ratio_val) + ', logMHI r=' + fmt(r_mhi_val) + ', logL36 r=' + fmt(r_l36_val));
console.log('    Two-var LOO=' + fmt(loo2) + ' vs ratio LOO=' + fmt(looR2(highV.map(g => g.logRatio), y)));
console.log('    MHI_perp r=' + fmt(r_mhi_perp) + ', L36_perp r=' + fmt(r_l36_perp));
console.log('    Ratio survives Mbar: pr=' + fmt(pr_ratio_mbar));
console.log('    Ratio survives structure: pr=' + fmt(pr_struct));
console.log('    Ratio survives SFE: pr=' + fmt(pr_ratio_sfe));

const output = {
  phase: 120,
  title: 'Mechanism Map: Why does log(MHI/L3.6) predict outer support requirement?',
  N: N,
  table1_component: table1Results.map(t => ({ name: t.name, r: t.r, loo: t.loo, aic: t.aic, bic: t.bic })),
  table2_residual: {
    MHI_perp: { r: r_mhi_perp, loo: loo_mhi_perp },
    L36_perp: { r: r_l36_perp, loo: loo_l36_perp },
    ratio_surviving_structure: pr_struct,
    ratio_surviving_mass: pr_ratio_mbar,
    ratio_surviving_SFE: pr_ratio_sfe,
  },
  table3_verdict: {
    gas_only: gasOnlyVerdict,
    stars_only: starsOnlyVerdict,
    ratio_state: ratioStateVerdict,
    structure_only: structOnlyVerdict,
    halo_mediated: haloMedVerdict,
  },
  target_sentence_achieved: canSay,
};

const outPath = path.join(__dirname, '..', 'public', 'phase120-mechanism-map.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);

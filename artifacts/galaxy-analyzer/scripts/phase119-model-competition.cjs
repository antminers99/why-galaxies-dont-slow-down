#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 119: MODEL COMPETITION');
console.log('');
console.log('  Three competing physical interpretations of WHY');
console.log('  gas-to-stellar balance predicts outer support requirement:');
console.log('');
console.log('  MODEL A: HALO CONCENTRATION');
console.log('    Gas-rich galaxies have lower-concentration halos → more');
console.log('    extended DM → higher outer mass discrepancy.');
console.log('    Prediction: logOMD anticorrelates with halo concentration.');
console.log('');
console.log('  MODEL B: BARYON-HALO COUPLING / ADIABATIC CONTRACTION');
console.log('    Stellar-dominated galaxies compressed their halo more');
console.log('    (adiabatic contraction), reducing outer DM support.');
console.log('    Prediction: logOMD correlates with baryon concentration.');
console.log('');
console.log('  MODEL C: EVOLUTIONARY STATE / ASSEMBLY HISTORY');
console.log('    Gas-rich = less evolved. Less star formation means less');
console.log('    baryonic concentration, preserving original halo state.');
console.log('    Prediction: fgas is a proxy for evolutionary state, not');
console.log('    a causal driver. Star formation indicators should do');
console.log('    equally well as fgas.');
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
  return { beta: aug.map((row, i) => row[p] / row[i]) };
}
function partialR(x, y, controls) {
  const n = x.length;
  if (n < controls.length + 5) return NaN;
  const X_ctrl = x.map((_, i) => [1, ...controls.map(c => c[i])]);
  const fitX = olsRegress(X_ctrl, x);
  const fitY = olsRegress(X_ctrl, y);
  if (!fitX || !fitY) return NaN;
  const resX = x.map((v, i) => v - X_ctrl[i].reduce((s, c, j) => s + c * fitX.beta[j], 0));
  const resY = y.map((v, i) => v - X_ctrl[i].reduce((s, c, j) => s + c * fitY.beta[j], 0));
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

function fitIsoHalo(dmPoints) {
  let bestRho0 = NaN, bestRc = NaN, bestChiSq = Infinity;
  for (let logRho0 = 5; logRho0 <= 9; logRho0 += 0.1) {
    for (let logRc = -1; logRc <= 2; logRc += 0.1) {
      const rho0 = Math.pow(10, logRho0);
      const rc = Math.pow(10, logRc);
      let chiSq = 0;
      for (const p of dmPoints) {
        const x = p.r / rc;
        const vdm_model_sq = 4 * Math.PI * 4.3009e-6 * rho0 * rc * rc * (1 - Math.atan(x) / x);
        chiSq += (p.vdm_sq - Math.max(vdm_model_sq, 0)) ** 2 / (Math.max(p.vdm_sq, 1) ** 2);
      }
      if (chiSq < bestChiSq) { bestChiSq = chiSq; bestRho0 = rho0; bestRc = rc; }
    }
  }
  return { rho0: bestRho0, rc: bestRc };
}

const galaxies = [];

for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat <= 0) continue;
  const sorted = rcPoints.filter(p => p.rad > 0 && p.vobs > 0).sort((a, b) => a.rad - b.r);
  if (sorted.length < 8) continue;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logMHI_L36 = Math.log10(t1.MHI) - Math.log10(t1.L36);
  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
  const logL36 = Math.log10(t1.L36);
  if (!isFinite(fgas) || !isFinite(logMHI_L36)) continue;

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

  const Mstar = t1.L36 * UPSILON_DISK * 1e9;
  const Mbar = Mstar + t1.MHI * 1.33 * 1e9;
  const Vflat = t1.Vflat;
  const Mdyn = 2.326e5 * Vflat ** 2 * sorted[sorted.length - 1].rad;
  const baryonFrac = Mbar / Math.max(Mdyn, 1);

  const barConc = t1.Rdisk > 0 ? Math.log10(sorted[sorted.length - 1].rad / t1.Rdisk) : NaN;
  const SFE = t1.L36 / (t1.MHI + t1.L36 * UPSILON_DISK);

  let innerDMfrac = NaN;
  const innerPts = sorted.filter(p => p.rad <= t1.Rdisk);
  if (innerPts.length > 0) {
    const innerVdm = innerPts.map(p => {
      const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                     UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                     p.vgas * Math.abs(p.vgas);
      return Math.max(p.vobs ** 2 - vBarSq, 0) / (p.vobs ** 2);
    });
    innerDMfrac = mean(innerVdm);
  }

  galaxies.push({
    name, fgas, logMHI_L36, logSigma0, logL36, Vflat, logOMD,
    logRho0, logRc, baryonFrac, barConc, SFE, innerDMfrac,
    logMbar: Math.log10(Mbar), logMdyn: Math.log10(Mdyn),
  });
}

const highV = galaxies.filter(g => g.Vflat >= 70);
console.log('  High-Vflat sample: N=' + highV.length + '\n');

// ============================================================
// MODEL A: HALO CONCENTRATION
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  MODEL A: HALO CONCENTRATION');
console.log('  Prediction: gas-rich → lower concentration → higher logOMD');
console.log('  Test: logRho0 (proxy for concentration) correlates with');
console.log('  both fgas and logOMD.');
console.log('══════════════════════════════════════════════════════════════\n');

const validA = highV.filter(g => isFinite(g.logRho0) && isFinite(g.logRc));
console.log('  N with halo fits = ' + validA.length);

if (validA.length >= 15) {
  const r_fgas_rho = pearsonR(validA.map(g => g.fgas), validA.map(g => g.logRho0));
  const r_fgas_rc = pearsonR(validA.map(g => g.fgas), validA.map(g => g.logRc));
  const r_rho_omd = pearsonR(validA.map(g => g.logRho0), validA.map(g => g.logOMD));
  const r_rc_omd = pearsonR(validA.map(g => g.logRc), validA.map(g => g.logOMD));

  console.log('  fgas → logRho0: r = ' + r_fgas_rho.toFixed(3));
  console.log('  fgas → logRc: r = ' + r_fgas_rc.toFixed(3));
  console.log('  logRho0 → logOMD: r = ' + r_rho_omd.toFixed(3));
  console.log('  logRc → logOMD: r = ' + r_rc_omd.toFixed(3));

  const pr_fgas_rho = partialR(
    validA.map(g => g.fgas), validA.map(g => g.logOMD),
    [validA.map(g => g.logRho0)]);
  console.log('  partial r(fgas, logOMD | logRho0) = ' + pr_fgas_rho.toFixed(3));

  const modelA_support = Math.abs(r_fgas_rho) > 0.3 && Math.abs(r_rho_omd) > 0.3;
  console.log('\n  Model A supported? ' + (modelA_support ? 'YES — halo concentration mediates' : 'WEAK/NO'));
}
console.log('');

// ============================================================
// MODEL B: BARYON-HALO COUPLING / ADIABATIC CONTRACTION
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  MODEL B: BARYON-HALO COUPLING');
console.log('  Prediction: higher baryon concentration → more adiabatic');
console.log('  contraction → lower outer DM → lower logOMD');
console.log('  Test: baryon fraction and concentration correlate with logOMD');
console.log('══════════════════════════════════════════════════════════════\n');

const validB = highV.filter(g => isFinite(g.barConc) && isFinite(g.baryonFrac));

if (validB.length >= 15) {
  const r_barFrac_omd = pearsonR(validB.map(g => g.baryonFrac), validB.map(g => g.logOMD));
  const r_barConc_omd = pearsonR(validB.map(g => g.barConc), validB.map(g => g.logOMD));
  const r_fgas_barFrac = pearsonR(validB.map(g => g.fgas), validB.map(g => g.baryonFrac));
  const r_fgas_barConc = pearsonR(validB.map(g => g.fgas), validB.map(g => g.barConc));

  console.log('  N = ' + validB.length);
  console.log('  baryonFrac → logOMD: r = ' + r_barFrac_omd.toFixed(3));
  console.log('  baryon concentration → logOMD: r = ' + r_barConc_omd.toFixed(3));
  console.log('  fgas → baryonFrac: r = ' + r_fgas_barFrac.toFixed(3));
  console.log('  fgas → baryon concentration: r = ' + r_fgas_barConc.toFixed(3));

  const pr_fgas_barConc = partialR(
    validB.map(g => g.fgas), validB.map(g => g.logOMD),
    [validB.map(g => g.barConc), validB.map(g => g.baryonFrac)]);
  console.log('  partial r(fgas, logOMD | barConc + barFrac) = ' + pr_fgas_barConc.toFixed(3));

  const modelB_support = Math.abs(r_barConc_omd) > 0.3;
  console.log('\n  Model B supported? ' + (modelB_support ? 'YES — baryon structure mediates' : 'WEAK/NO'));
}
console.log('');

// ============================================================
// MODEL C: EVOLUTIONARY STATE
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  MODEL C: EVOLUTIONARY STATE / ASSEMBLY HISTORY');
console.log('  Prediction: fgas = proxy for evolutionary state.');
console.log('  Test: SFE (star formation efficiency L3.6/Mbar) should');
console.log('  predict logOMD equally well or better than fgas.');
console.log('  If fgas >> SFE → fgas is more fundamental than evolution.');
console.log('  If SFE ≈ fgas → both are proxies for same thing.');
console.log('══════════════════════════════════════════════════════════════\n');

const validC = highV.filter(g => isFinite(g.SFE) && g.SFE > 0);

if (validC.length >= 15) {
  const logSFE = validC.map(g => Math.log10(g.SFE));
  const r_fgas_omd = pearsonR(validC.map(g => g.fgas), validC.map(g => g.logOMD));
  const r_sfe_omd = pearsonR(logSFE, validC.map(g => g.logOMD));
  const loo_fgas = looR2(validC.map(g => g.fgas), validC.map(g => g.logOMD));
  const loo_sfe = looR2(logSFE, validC.map(g => g.logOMD));

  console.log('  N = ' + validC.length);
  console.log('  fgas → logOMD: r = ' + r_fgas_omd.toFixed(3) + ', LOO = ' + loo_fgas.toFixed(3));
  console.log('  log(SFE) → logOMD: r = ' + r_sfe_omd.toFixed(3) + ', LOO = ' + loo_sfe.toFixed(3));

  const r_fgas_sfe = pearsonR(validC.map(g => g.fgas), logSFE);
  console.log('  fgas vs SFE correlation: r = ' + r_fgas_sfe.toFixed(3));

  const pr_fgas_sfe = partialR(validC.map(g => g.fgas), validC.map(g => g.logOMD), [logSFE]);
  const pr_sfe_fgas = partialR(logSFE, validC.map(g => g.logOMD), [validC.map(g => g.fgas)]);

  console.log('  partial r(fgas, logOMD | SFE) = ' + pr_fgas_sfe.toFixed(3));
  console.log('  partial r(SFE, logOMD | fgas) = ' + pr_sfe_fgas.toFixed(3));

  const modelC_support = Math.abs(r_sfe_omd) >= Math.abs(r_fgas_omd) - 0.05;
  console.log('\n  Model C supported? ' + (modelC_support ? 'YES — SFE and fgas are equally good (shared underlying state)' : 'NO — fgas outperforms SFE, suggesting gas ratio is more fundamental than evolutionary state'));
}
console.log('');

// ============================================================
// MODEL COMPARISON: WHICH WINS?
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  MODEL COMPARISON DEATH MATCH');
console.log('  Combined model: logOMD ~ fgas + logRho0 + barConc + SFE');
console.log('  Which variables survive in the combined model?');
console.log('══════════════════════════════════════════════════════════════\n');

const validAll = highV.filter(g =>
  isFinite(g.logRho0) && isFinite(g.barConc) && isFinite(g.SFE) && g.SFE > 0);

if (validAll.length >= 20) {
  const y = validAll.map(g => g.logOMD);
  const preds = [
    { name: 'fgas', vals: validAll.map(g => g.fgas) },
    { name: 'logRho0', vals: validAll.map(g => g.logRho0) },
    { name: 'barConc', vals: validAll.map(g => g.barConc) },
    { name: 'log(SFE)', vals: validAll.map(g => Math.log10(g.SFE)) },
    { name: 'log(MHI/L36)', vals: validAll.map(g => g.logMHI_L36) },
    { name: 'innerDMfrac', vals: validAll.filter(g => isFinite(g.innerDMfrac)).length === validAll.length ? validAll.map(g => g.innerDMfrac) : null },
  ].filter(p => p.vals !== null);

  console.log('  N = ' + validAll.length + '\n');
  console.log('  Individual predictions of logOMD:');
  for (const pred of preds) {
    const r = pearsonR(pred.vals, y);
    const loo = looR2(pred.vals, y);
    console.log('    ' + pred.name + ': r=' + r.toFixed(3) + ', LOO=' + loo.toFixed(3));
  }

  console.log('\n  Partial correlations (each controlling ALL others):');
  for (let i = 0; i < preds.length; i++) {
    const controls = preds.filter((_, j) => j !== i).map(p => p.vals);
    const pr = partialR(preds[i].vals, y, controls);
    console.log('    partial r(' + preds[i].name + ' | all others) = ' + pr.toFixed(3));
  }
}

// ============================================================
// OVERALL VERDICT
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  OVERALL VERDICT: WHICH MODEL BEST EXPLAINS THE SIGNAL?');
console.log('══════════════════════════════════════════════════════════════\n');

console.log('  The answer to "WHY does gas-to-stellar balance predict');
console.log('  outer support requirement?" depends on which model wins.');
console.log('  See individual model results above for the answer.\n');

const output = {
  phase: 119,
  title: 'Model Competition: Which physical mechanism explains the signal?',
  nHighV: highV.length,
};

const outPath = path.join(__dirname, '..', 'public', 'phase119-model-competition.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('  Results saved to: ' + outPath);

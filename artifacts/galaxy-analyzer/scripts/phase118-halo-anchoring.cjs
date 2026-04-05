#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 118: HALO ANCHORING');
console.log('');
console.log('  CRITICAL QUESTION: Does gas-to-stellar balance connect to');
console.log('  DARK MATTER HALO properties? If yes → bridge to physics.');
console.log('');
console.log('  From Vobs and Vbar we derive Vdm² = Vobs² - Vbar².');
console.log('  From Vdm(r) we fit pseudo-isothermal halo profiles to get:');
console.log('    - rho0 (central DM density)');
console.log('    - rc (core radius)');
console.log('    - rho0 × rc (surface density, Donato+ 2009)');
console.log('    - inner DM fraction');
console.log('    - halo compactness');
console.log('  Then test: does gas-to-stellar balance predict these?');
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
        const residSq = (p.vdm_sq - Math.max(vdm_model_sq, 0)) ** 2;
        chiSq += residSq / (Math.max(p.vdm_sq, 1) ** 2);
      }

      if (chiSq < bestChiSq) {
        bestChiSq = chiSq;
        bestRho0 = rho0;
        bestRc = rc;
      }
    }
  }
  return { rho0: bestRho0, rc: bestRc, chiSq: bestChiSq };
}

const galaxies = [];

for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat <= 0) continue;
  const sorted = rcPoints.filter(p => p.rad > 0 && p.vobs > 0).sort((a, b) => a.rad - b.rad);
  if (sorted.length < 8) continue;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logMHI_L36 = Math.log10(t1.MHI) - Math.log10(t1.L36);
  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
  if (!isFinite(fgas) || !isFinite(logMHI_L36)) continue;

  const dmPoints = [];
  let innerDMfrac_num = 0, innerDMfrac_den = 0;

  for (const p of sorted) {
    const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                   UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                   p.vgas * Math.abs(p.vgas);
    const vdm_sq = p.vobs ** 2 - Math.max(vBarSq, 0);

    dmPoints.push({ r: p.rad, vdm_sq, vobs_sq: p.vobs ** 2, vbar_sq: vBarSq });

    if (p.rad <= t1.Rdisk) {
      innerDMfrac_num += Math.max(vdm_sq, 0);
      innerDMfrac_den += p.vobs ** 2;
    }
  }

  if (dmPoints.length < 5) continue;

  const positiveDM = dmPoints.filter(p => p.vdm_sq > 0);
  if (positiveDM.length < 4) continue;

  const haloFit = fitIsoHalo(positiveDM);
  if (!isFinite(haloFit.rho0) || !isFinite(haloFit.rc)) continue;

  const logRho0 = Math.log10(haloFit.rho0);
  const logRc = Math.log10(haloFit.rc);
  const logRho0Rc = Math.log10(haloFit.rho0 * haloFit.rc);
  const innerDMfrac = innerDMfrac_den > 0 ? innerDMfrac_num / innerDMfrac_den : NaN;

  const rLast = sorted[sorted.length - 1].rad;
  const haloCompactness = haloFit.rc > 0 ? rLast / haloFit.rc : NaN;

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

  galaxies.push({
    name, fgas, logMHI_L36, logSigma0, Vflat: t1.Vflat, logOMD,
    logRho0, logRc, logRho0Rc, innerDMfrac,
    haloCompactness: isFinite(haloCompactness) ? haloCompactness : NaN,
  });
}

const highV = galaxies.filter(g => g.Vflat >= 70);
console.log('  Galaxies with halo fits: N=' + galaxies.length);
console.log('  High-Vflat: N=' + highV.length + '\n');

// ============================================================
// TEST 1: Gas state vs halo properties
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: GAS STATE vs HALO PROPERTIES');
console.log('  (Pseudo-isothermal halo fits)');
console.log('══════════════════════════════════════════════════════════════\n');

const haloPreds = [
  { name: 'fgas', get: g => g.fgas },
  { name: 'log(MHI/L36)', get: g => g.logMHI_L36 },
  { name: 'logSigma0', get: g => g.logSigma0 },
];

const haloTargets = [
  { name: 'logRho0 (central DM density)', get: g => g.logRho0 },
  { name: 'logRc (core radius)', get: g => g.logRc },
  { name: 'logRho0*Rc (surface density)', get: g => g.logRho0Rc },
  { name: 'innerDMfrac (DM frac at r<Rd)', get: g => g.innerDMfrac },
];

const test1Results = {};
for (const ht of haloTargets) {
  const valid = highV.filter(g => isFinite(ht.get(g)));
  console.log('  ' + ht.name + ' (N=' + valid.length + ')');
  test1Results[ht.name] = {};

  for (const pred of haloPreds) {
    const x = valid.map(pred.get);
    const y = valid.map(ht.get);
    const r = pearsonR(x, y);
    const loo = valid.length >= 15 ? looR2(x, y) : NaN;
    const marker = Math.abs(r) > 0.4 ? ' <<<' : Math.abs(r) > 0.3 ? ' <<' : '';
    console.log('    ' + pred.name + ': r=' + r.toFixed(3) + (isFinite(loo) ? ', LOO=' + loo.toFixed(3) : '') + marker);
    test1Results[ht.name][pred.name] = { r: parseFloat(r.toFixed(3)) };
  }
  console.log('');
}

// ============================================================
// TEST 2: Halo properties vs logOMD
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: HALO PROPERTIES vs logOMD');
console.log('  Which halo property best predicts outer mass discrepancy?');
console.log('══════════════════════════════════════════════════════════════\n');

const test2Results = {};
for (const ht of haloTargets) {
  const valid = highV.filter(g => isFinite(ht.get(g)));
  const r = pearsonR(valid.map(ht.get), valid.map(g => g.logOMD));
  const loo = valid.length >= 15 ? looR2(valid.map(ht.get), valid.map(g => g.logOMD)) : NaN;
  const marker = Math.abs(r) > 0.4 ? ' <<<' : Math.abs(r) > 0.3 ? ' <<' : '';
  console.log('  ' + ht.name + ' → logOMD: r=' + r.toFixed(3) + (isFinite(loo) ? ', LOO=' + loo.toFixed(3) : '') + marker);
  test2Results[ht.name] = { r: parseFloat(r.toFixed(3)), loo: isFinite(loo) ? parseFloat(loo.toFixed(3)) : null };
}
console.log('');

// ============================================================
// TEST 3: MEDIATION — does fgas predict logOMD THROUGH halo?
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: MEDIATION ANALYSIS');
console.log('  Path: fgas → halo property → logOMD');
console.log('  Does fgas predict logOMD THROUGH halo properties?');
console.log('══════════════════════════════════════════════════════════════\n');

const validAll = highV.filter(g =>
  isFinite(g.logRho0) && isFinite(g.logRc) && isFinite(g.innerDMfrac));

if (validAll.length >= 15) {
  const pr_fgas_direct = pearsonR(validAll.map(g => g.fgas), validAll.map(g => g.logOMD));
  console.log('  Direct: r(fgas, logOMD) = ' + pr_fgas_direct.toFixed(3));

  const pr_fgas_rho = partialR(
    validAll.map(g => g.fgas), validAll.map(g => g.logOMD),
    [validAll.map(g => g.logRho0)]);
  const pr_fgas_rc = partialR(
    validAll.map(g => g.fgas), validAll.map(g => g.logOMD),
    [validAll.map(g => g.logRc)]);
  const pr_fgas_dmfrac = partialR(
    validAll.map(g => g.fgas), validAll.map(g => g.logOMD),
    [validAll.map(g => g.innerDMfrac)]);
  const pr_fgas_allHalo = partialR(
    validAll.map(g => g.fgas), validAll.map(g => g.logOMD),
    [validAll.map(g => g.logRho0), validAll.map(g => g.logRc), validAll.map(g => g.innerDMfrac)]);

  console.log('  Partial r(fgas, logOMD | logRho0) = ' + pr_fgas_rho.toFixed(3));
  console.log('  Partial r(fgas, logOMD | logRc) = ' + pr_fgas_rc.toFixed(3));
  console.log('  Partial r(fgas, logOMD | innerDMfrac) = ' + pr_fgas_dmfrac.toFixed(3));
  console.log('  Partial r(fgas, logOMD | ALL halo) = ' + pr_fgas_allHalo.toFixed(3));
  console.log('');

  const mediated = Math.abs(pr_fgas_allHalo) < Math.abs(pr_fgas_direct) * 0.5;
  console.log('  fgas signal drops >50% after controlling ALL halo? ' + (mediated ? 'YES → MEDIATED through halo' : 'NO → DIRECT effect beyond halo'));
  console.log('  Remaining signal: ' + (Math.abs(pr_fgas_allHalo) / Math.abs(pr_fgas_direct) * 100).toFixed(0) + '% of original');
  console.log('');

  const pr_rho_fgas = partialR(
    validAll.map(g => g.logRho0), validAll.map(g => g.logOMD),
    [validAll.map(g => g.fgas)]);
  const pr_rc_fgas = partialR(
    validAll.map(g => g.logRc), validAll.map(g => g.logOMD),
    [validAll.map(g => g.fgas)]);
  const pr_dmfrac_fgas = partialR(
    validAll.map(g => g.innerDMfrac), validAll.map(g => g.logOMD),
    [validAll.map(g => g.fgas)]);

  console.log('  Reverse: do halo properties survive controlling fgas?');
  console.log('    partial r(logRho0, logOMD | fgas) = ' + pr_rho_fgas.toFixed(3));
  console.log('    partial r(logRc, logOMD | fgas) = ' + pr_rc_fgas.toFixed(3));
  console.log('    partial r(innerDMfrac, logOMD | fgas) = ' + pr_dmfrac_fgas.toFixed(3));
}

// ============================================================
// TEST 4: HALO PROPERTIES BY GAS TERCILE
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: HALO PROPERTIES BY GAS-FRACTION TERCILE');
console.log('══════════════════════════════════════════════════════════════\n');

const sortedByFgas = [...highV].filter(g => isFinite(g.logRho0) && isFinite(g.logRc))
  .sort((a, b) => a.fgas - b.fgas);
const t3 = Math.floor(sortedByFgas.length / 3);
const terciles = [
  { name: 'Low fgas (stellar-dominated)', g: sortedByFgas.slice(0, t3) },
  { name: 'Mid fgas', g: sortedByFgas.slice(t3, 2 * t3) },
  { name: 'High fgas (gas-dominated)', g: sortedByFgas.slice(2 * t3) },
];

for (const t of terciles) {
  console.log('  ' + t.name + ' (N=' + t.g.length + ')');
  console.log('    logRho0: ' + mean(t.g.map(g => g.logRho0)).toFixed(2) + ' +/- ' + sd(t.g.map(g => g.logRho0)).toFixed(2));
  console.log('    logRc:   ' + mean(t.g.map(g => g.logRc)).toFixed(2) + ' +/- ' + sd(t.g.map(g => g.logRc)).toFixed(2));
  console.log('    innerDM: ' + mean(t.g.filter(g => isFinite(g.innerDMfrac)).map(g => g.innerDMfrac)).toFixed(3));
  console.log('    logOMD:  ' + mean(t.g.map(g => g.logOMD)).toFixed(3));
  console.log('');
}

// ============================================================
// VERDICT
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  OVERALL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const strongHaloLinks = [];
for (const ht of haloTargets) {
  for (const pred of haloPreds) {
    const r = test1Results[ht.name]?.[pred.name]?.r || 0;
    if (Math.abs(r) > 0.3) strongHaloLinks.push(pred.name + ' → ' + ht.name + ' (r=' + r + ')');
  }
}

if (strongHaloLinks.length > 2) {
  console.log('  STRONG HALO ANCHORING: Gas-to-stellar balance connects to DM halo properties.');
} else if (strongHaloLinks.length > 0) {
  console.log('  PARTIAL HALO ANCHORING: Some halo connections exist.');
} else {
  console.log('  WEAK HALO ANCHORING: Gas state does not strongly predict halo properties.');
  console.log('  The signal may operate through a different physical channel.');
}
for (const link of strongHaloLinks) console.log('    ' + link);

const output = {
  phase: 118,
  title: 'Halo Anchoring: Does gas state connect to dark matter halo properties?',
  nHighV: highV.length,
  test1_gasVsHalo: test1Results,
  test2_haloVsLogOMD: test2Results,
  strongHaloLinks,
};

const outPath = path.join(__dirname, '..', 'public', 'phase118-halo-anchoring.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);

#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 115: BTFR RESIDUAL ANALYSIS');
console.log('  Ref: Verheijen 2001, McGaugh+ 2012, Lelli+ 2019');
console.log('');
console.log('  The Baryonic Tully-Fisher Relation (BTFR) links');
console.log('  Mbar (total baryonic mass) to Vflat.');
console.log('  Question: Does gas-to-stellar balance predict where');
console.log('  galaxies scatter off the BTFR?');
console.log('  If YES -> gas state drives BTFR scatter -> connects');
console.log('  our finding to the most fundamental scaling relation.');
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
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const rad = parseFloat(line.substring(19, 25).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, vdisk, vbul: vbul || 0, vgas });
}

const galaxies = [];

for (const [name, t1] of Object.entries(table1)) {
  if (!t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat <= 0) continue;

  const Mstar = t1.L36 * UPSILON_DISK * 1e9;
  const Mgas = t1.MHI * 1.33 * 1e9;
  const Mbar = Mstar + Mgas;

  const logMbar = Math.log10(Mbar);
  const logVflat = Math.log10(t1.Vflat);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logMHI_L36 = Math.log10(t1.MHI) - Math.log10(t1.L36);
  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
  const logL36 = Math.log10(t1.L36);

  if (!isFinite(logMbar) || !isFinite(logVflat) || !isFinite(fgas) || !isFinite(logMHI_L36)) continue;

  const rc = rcByGalaxy[name];
  let logOMD = NaN;
  if (rc && rc.length >= 8) {
    const sorted = rc.filter(p => p.rad > 0 && p.vobs > 0).sort((a, b) => a.rad - b.rad);
    const half = Math.floor(sorted.length / 2);
    const outerPts = sorted.slice(half);
    const outerMD = outerPts.map(p => {
      const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                     UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                     p.vgas * Math.abs(p.vgas);
      return (p.vobs ** 2) / Math.max(vBarSq, 0.01);
    }).filter(v => isFinite(v) && v > 0);
    logOMD = Math.log10(mean(outerMD));
  }

  galaxies.push({
    name, logMbar, logVflat, fgas, logMHI_L36, logSigma0, logL36,
    Vflat: t1.Vflat, logOMD,
  });
}

console.log('  Total galaxies with BTFR data: N=' + galaxies.length + '\n');

// ============================================================
// FIT THE BTFR: log(Mbar) = a + b * log(Vflat)
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  BTFR FIT: log(Mbar) = a + b * log(Vflat)');
console.log('══════════════════════════════════════════════════════════════\n');

const x_all = galaxies.map(g => g.logVflat);
const y_all = galaxies.map(g => g.logMbar);
const fitAll = olsRegress(x_all.map(v => [1, v]), y_all);

console.log('  Full sample (N=' + galaxies.length + '):');
console.log('  Intercept = ' + fitAll.beta[0].toFixed(3));
console.log('  Slope = ' + fitAll.beta[1].toFixed(3));
console.log('  (Literature: slope ~ 3.5-4.0)\n');

const btfrResiduals = galaxies.map((g, i) => g.logMbar - (fitAll.beta[0] + fitAll.beta[1] * g.logVflat));
const btfrRMS = Math.sqrt(mean(btfrResiduals.map(r => r ** 2)));
console.log('  BTFR RMS scatter: ' + btfrRMS.toFixed(4) + ' dex\n');

galaxies.forEach((g, i) => g.btfrResid = btfrResiduals[i]);

// ============================================================
// TEST 1: What predicts BTFR residuals?
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: WHAT PREDICTS BTFR RESIDUALS?');
console.log('  (Full sample, N=' + galaxies.length + ')');
console.log('══════════════════════════════════════════════════════════════\n');

const preds = [
  { name: 'fgas', get: g => g.fgas },
  { name: 'log(MHI/L36)', get: g => g.logMHI_L36 },
  { name: 'logSigma0', get: g => g.logSigma0 },
  { name: 'logL36', get: g => g.logL36 },
];

const test1Results = {};
for (const pred of preds) {
  const x = galaxies.map(pred.get);
  const y = btfrResiduals;
  const r = pearsonR(x, y);
  const loo = looR2(x, y);
  const marker = Math.abs(r) > 0.4 ? ' <<<' : Math.abs(r) > 0.3 ? ' <<' : '';
  console.log('  ' + pred.name + ': r=' + r.toFixed(3) + ', LOO=' + loo.toFixed(3) + marker);
  test1Results[pred.name] = { r: parseFloat(r.toFixed(3)), loo: parseFloat(loo.toFixed(3)) };
}
console.log('');

// ============================================================
// TEST 2: Same but high-Vflat only
// ============================================================

const highV = galaxies.filter(g => g.Vflat >= 70);
const x_hv = highV.map(g => g.logVflat);
const y_hv = highV.map(g => g.logMbar);
const fitHV = olsRegress(x_hv.map(v => [1, v]), y_hv);
const btfrResidHV = highV.map(g => g.logMbar - (fitHV.beta[0] + fitHV.beta[1] * g.logVflat));
highV.forEach((g, i) => g.btfrResidHV = btfrResidHV[i]);

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: HIGH-VFLAT BTFR RESIDUALS (N=' + highV.length + ')');
console.log('  BTFR slope = ' + fitHV.beta[1].toFixed(3) + ', RMS = ' + Math.sqrt(mean(btfrResidHV.map(r => r ** 2))).toFixed(4) + ' dex');
console.log('══════════════════════════════════════════════════════════════\n');

const test2Results = {};
for (const pred of preds) {
  const x = highV.map(pred.get);
  const y = btfrResidHV;
  const r = pearsonR(x, y);
  const loo = looR2(x, y);
  const marker = Math.abs(r) > 0.4 ? ' <<<' : Math.abs(r) > 0.3 ? ' <<' : '';
  console.log('  ' + pred.name + ': r=' + r.toFixed(3) + ', LOO=' + loo.toFixed(3) + marker);
  test2Results[pred.name] = { r: parseFloat(r.toFixed(3)), loo: parseFloat(loo.toFixed(3)) };
}
console.log('');

// ============================================================
// TEST 3: BTFR residual vs logOMD
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: BTFR RESIDUAL vs logOMD');
console.log('  Does BTFR offset connect to outer mass discrepancy?');
console.log('══════════════════════════════════════════════════════════════\n');

const validOMD = highV.filter(g => isFinite(g.logOMD));
if (validOMD.length >= 15) {
  const r_btfr_omd = pearsonR(validOMD.map(g => g.btfrResidHV), validOMD.map(g => g.logOMD));
  console.log('  r(BTFR residual, logOMD) = ' + r_btfr_omd.toFixed(3));
  console.log('  (N = ' + validOMD.length + ')\n');

  const pr_fgas_btfr = partialR(
    validOMD.map(g => g.fgas), validOMD.map(g => g.logOMD),
    [validOMD.map(g => g.btfrResidHV)]);
  const pr_btfr_fgas = partialR(
    validOMD.map(g => g.btfrResidHV), validOMD.map(g => g.logOMD),
    [validOMD.map(g => g.fgas)]);

  console.log('  Partial r(fgas, logOMD | BTFR resid) = ' + pr_fgas_btfr.toFixed(3));
  console.log('  Partial r(BTFR resid, logOMD | fgas) = ' + pr_btfr_fgas.toFixed(3));
  console.log('  -> fgas survives controlling BTFR residual? ' + (Math.abs(pr_fgas_btfr) > 0.3 ? 'YES' : 'NO'));
  console.log('  -> BTFR resid survives controlling fgas? ' + (Math.abs(pr_btfr_fgas) > 0.3 ? 'YES' : 'NO'));
  console.log('');
}

// ============================================================
// TEST 4: Tercile analysis — BTFR scatter by gas state
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 4: BTFR SCATTER BY GAS-FRACTION TERCILE');
console.log('══════════════════════════════════════════════════════════════\n');

const sortedByFgas = [...highV].sort((a, b) => a.fgas - b.fgas);
const t3 = Math.floor(sortedByFgas.length / 3);
const tercGroups = [
  { name: 'Low fgas', g: sortedByFgas.slice(0, t3) },
  { name: 'Mid fgas', g: sortedByFgas.slice(t3, 2 * t3) },
  { name: 'High fgas', g: sortedByFgas.slice(2 * t3) },
];

const tercileResults = [];
for (const t of tercGroups) {
  const resids = t.g.map(g => g.btfrResidHV);
  const rms = Math.sqrt(mean(resids.map(r => r ** 2)));
  const meanResid = mean(resids);
  console.log('  ' + t.name + ' (N=' + t.g.length + ', fgas=[' + Math.min(...t.g.map(g => g.fgas)).toFixed(3) + ', ' + Math.max(...t.g.map(g => g.fgas)).toFixed(3) + '])');
  console.log('    Mean BTFR residual: ' + meanResid.toFixed(4) + ' dex');
  console.log('    RMS scatter: ' + rms.toFixed(4) + ' dex');
  console.log('');
  tercileResults.push({ name: t.name, n: t.g.length, meanResid: parseFloat(meanResid.toFixed(4)), rms: parseFloat(rms.toFixed(4)) });
}

const meanDiff = tercileResults[2].meanResid - tercileResults[0].meanResid;
console.log('  High-Low fgas mean BTFR residual diff: ' + meanDiff.toFixed(4) + ' dex');
console.log('  Gas-rich galaxies sit ' + (meanDiff > 0 ? 'ABOVE' : 'BELOW') + ' gas-poor on BTFR\n');

// ============================================================
// TEST 5: BTFR residual vs gas state — circularity check
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 5: CIRCULARITY CHECK');
console.log('  fgas contains MHI; Mbar contains MHI. Is the BTFR');
console.log('  residual correlation just definitional overlap?');
console.log('══════════════════════════════════════════════════════════════\n');

const gasOnlyBTFR = highV.map(g => {
  const Mgas = Math.pow(10, g.logMHI_L36) * Math.pow(10, g.logL36) * 1.33 * 1e9;
  const logMgas = Math.log10(Mgas);
  return logMgas;
});

const r_gasMass_btfr = pearsonR(gasOnlyBTFR, btfrResidHV);
console.log('  r(log Mgas, BTFR residual) = ' + r_gasMass_btfr.toFixed(3));

const starMass = highV.map(g => g.logL36 + Math.log10(UPSILON_DISK * 1e9));
const r_starMass_btfr = pearsonR(starMass, btfrResidHV);
console.log('  r(log Mstar, BTFR residual) = ' + r_starMass_btfr.toFixed(3));
console.log('');

if (Math.abs(r_gasMass_btfr) > 0.5) {
  console.log('  WARNING: Gas mass directly correlates with BTFR residual.');
  console.log('  This is expected because Mbar = Mstar + Mgas.');
  console.log('  The fgas-BTFR correlation may be PARTIALLY definitional.\n');
}

const pr_fgas_btfr_stellar = partialR(
  highV.map(g => g.fgas), btfrResidHV,
  [starMass]);
console.log('  Partial r(fgas, BTFR resid | Mstar) = ' + pr_fgas_btfr_stellar.toFixed(3));
console.log('  -> fgas signal after removing stellar mass effect: ' + (Math.abs(pr_fgas_btfr_stellar) > 0.3 ? 'SURVIVES' : 'dies') + '\n');

// ============================================================
// TEST 6: Does gas-stellar balance reduce BTFR scatter?
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 6: CAN GAS STATE REDUCE BTFR SCATTER?');
console.log('  Fit: log(Mbar) = a + b*log(Vflat) + c*fgas');
console.log('══════════════════════════════════════════════════════════════\n');

const X_2d = highV.map(g => [1, g.logVflat, g.fgas]);
const y_2d = highV.map(g => g.logMbar);
const fit2d = olsRegress(X_2d, y_2d);

if (fit2d) {
  const resid2d = highV.map((g, i) => g.logMbar - (fit2d.beta[0] + fit2d.beta[1] * g.logVflat + fit2d.beta[2] * g.fgas));
  const rms_1d = Math.sqrt(mean(btfrResidHV.map(r => r ** 2)));
  const rms_2d = Math.sqrt(mean(resid2d.map(r => r ** 2)));
  const improvement = ((rms_1d - rms_2d) / rms_1d * 100);

  console.log('  BTFR alone: RMS = ' + rms_1d.toFixed(4) + ' dex');
  console.log('  BTFR + fgas: RMS = ' + rms_2d.toFixed(4) + ' dex');
  console.log('  Scatter reduction: ' + improvement.toFixed(1) + '%');
  console.log('  fgas coefficient: ' + fit2d.beta[2].toFixed(3));
  console.log('');

  const X_3d = highV.map(g => [1, g.logVflat, g.logMHI_L36]);
  const fit3d = olsRegress(X_3d, y_2d);
  if (fit3d) {
    const resid3d = highV.map((g, i) => g.logMbar - (fit3d.beta[0] + fit3d.beta[1] * g.logVflat + fit3d.beta[2] * g.logMHI_L36));
    const rms_3d = Math.sqrt(mean(resid3d.map(r => r ** 2)));
    const improvement3 = ((rms_1d - rms_3d) / rms_1d * 100);
    console.log('  BTFR + log(MHI/L36): RMS = ' + rms_3d.toFixed(4) + ' dex');
    console.log('  Scatter reduction: ' + improvement3.toFixed(1) + '%');
    console.log('  log(MHI/L36) coefficient: ' + fit3d.beta[2].toFixed(3));
  }
}

// ============================================================
// VERDICT
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  OVERALL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const fgasBTFR_r = test1Results['fgas']?.r || 0;
const strongBTFR = Math.abs(fgasBTFR_r) > 0.3;

if (strongBTFR) {
  console.log('  Gas-to-stellar balance DOES predict BTFR residuals.');
  console.log('  But beware definitional overlap: Mbar contains Mgas.');
} else {
  console.log('  Gas-to-stellar balance does NOT strongly predict BTFR residuals.');
  console.log('  The BTFR and our logOMD finding may be independent axes.');
}

const output = {
  phase: 115,
  title: 'BTFR Residual Analysis: Does gas state predict Baryonic Tully-Fisher scatter?',
  nAll: galaxies.length,
  nHighV: highV.length,
  btfrFit: { intercept: parseFloat(fitAll.beta[0].toFixed(3)), slope: parseFloat(fitAll.beta[1].toFixed(3)), rms: parseFloat(btfrRMS.toFixed(4)) },
  test1_fullSample: test1Results,
  test2_highV: test2Results,
  test4_terciles: tercileResults,
};

const outPath = path.join(__dirname, '..', 'public', 'phase115-btfr-residuals.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);

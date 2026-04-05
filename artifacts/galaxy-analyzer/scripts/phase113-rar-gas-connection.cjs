#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 113: RAR-GAS CONNECTION');
console.log('  Does gas-to-stellar balance predict WHERE galaxies sit');
console.log('  relative to the Radial Acceleration Relation?');
console.log('');
console.log('  RAR: gobs = f(gbar) point-by-point (McGaugh+ 2016, Lelli+ 2017)');
console.log('  Our result: log(MHI/L36) predicts outer mass discrepancy (galaxy-level)');
console.log('  NEW QUESTION: Do these connect? Does gas-to-stellar state');
console.log('  predict RAR residuals systematically?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const G = 4.3009e-6;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function median(a) { const s = [...a].sort((x, y) => x - y); const m = Math.floor(s.length / 2); return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function percentile(a, p) { const s = [...a].sort((x, y) => x - y); const i = p * (s.length - 1); const lo = Math.floor(i); return lo >= s.length - 1 ? s[s.length - 1] : s[lo] + (i - lo) * (s[lo + 1] - s[lo]); }

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
  return { beta };
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
  const D = parseFloat(line.substring(15, 21).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  table1[name] = { D, L36, Rdisk, MHI, Vflat };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const evobs = parseFloat(line.substring(32, 39).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, evobs: isFinite(evobs) ? evobs : 0, vgas, vdisk, vbul: vbul || 0 });
}

const galaxies = [];
const allRARpoints = [];

for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat <= 0) continue;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logMHI_L36 = Math.log10(t1.MHI) - Math.log10(t1.L36);
  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
  if (!isFinite(fgas) || !isFinite(logMHI_L36)) continue;

  const pts = [];
  const rarPts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;

    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * pt.vbul * Math.abs(pt.vbul) +
                   pt.vgas * Math.abs(pt.vgas);
    const vBar = Math.sqrt(Math.max(vBarSq, 0.01));

    const gobs = pt.vobs ** 2 / (pt.rad * 3.086e16);
    const gbar = vBarSq / (pt.rad * 3.086e16);

    if (gbar > 0 && gobs > 0 && isFinite(gbar) && isFinite(gobs)) {
      const log_gbar = Math.log10(gbar * 1e3);
      const log_gobs = Math.log10(gobs * 1e3);

      if (isFinite(log_gbar) && isFinite(log_gobs)) {
        rarPts.push({ log_gbar, log_gobs, rad: pt.rad });
        allRARpoints.push({ log_gbar, log_gobs, name, fgas, logMHI_L36, logSigma0, Vflat: t1.Vflat, rad: pt.rad });
      }
    }

    pts.push({ r: pt.rad, vobs: pt.vobs, vBar });
  }

  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);
  const outerMD = outerPts.map(p => (p.vobs ** 2) / (p.vBar ** 2)).filter(v => isFinite(v) && v > 0);
  const logOMD = Math.log10(mean(outerMD));
  if (!isFinite(logOMD)) continue;

  const rarResiduals = rarPts.map(p => {
    const rar_pred = p.log_gbar + Math.log10(1 / (1 - Math.exp(-Math.sqrt(Math.pow(10, p.log_gbar - Math.log10(1.2e-10))))));
    return p.log_gobs - (isFinite(rar_pred) ? rar_pred : p.log_gbar);
  }).filter(v => isFinite(v));

  const outerRarPts = rarPts.filter(p => p.rad >= sorted[half].r);
  const outerRarResiduals = outerRarPts.map(p => {
    const rar_pred = p.log_gbar + Math.log10(1 / (1 - Math.exp(-Math.sqrt(Math.pow(10, p.log_gbar - Math.log10(1.2e-10))))));
    return p.log_gobs - (isFinite(rar_pred) ? rar_pred : p.log_gbar);
  }).filter(v => isFinite(v));

  const innerRarPts = rarPts.filter(p => p.rad < sorted[half].r);
  const innerRarResiduals = innerRarPts.map(p => {
    const rar_pred = p.log_gbar + Math.log10(1 / (1 - Math.exp(-Math.sqrt(Math.pow(10, p.log_gbar - Math.log10(1.2e-10))))));
    return p.log_gobs - (isFinite(rar_pred) ? rar_pred : p.log_gbar);
  }).filter(v => isFinite(v));

  galaxies.push({
    name, logOMD, fgas, logMHI_L36, logSigma0,
    Vflat: t1.Vflat,
    meanRARresid: rarResiduals.length > 0 ? mean(rarResiduals) : NaN,
    medianRARresid: rarResiduals.length > 0 ? median(rarResiduals) : NaN,
    outerRARresid: outerRarResiduals.length > 0 ? mean(outerRarResiduals) : NaN,
    innerRARresid: innerRarResiduals.length > 0 ? mean(innerRarResiduals) : NaN,
    rarScatter: rarResiduals.length > 2 ? sd(rarResiduals) : NaN,
    nRARpts: rarPts.length,
  });
}

const highV = galaxies.filter(g => g.Vflat >= 70 && isFinite(g.meanRARresid));
console.log('  Total galaxies with RAR data: N=' + galaxies.filter(g => isFinite(g.meanRARresid)).length);
console.log('  High-Vflat (>= 70): N=' + highV.length);
console.log('  Total RAR data points: ' + allRARpoints.length + '\n');

// ============================================================
// TEST 1: Galaxy-level RAR residuals vs gas-to-stellar balance
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: GALAXY-LEVEL RAR RESIDUALS vs GAS STATE');
console.log('  Mean RAR residual per galaxy vs fgas, log(MHI/L36), Sigma0');
console.log('══════════════════════════════════════════════════════════════\n');

const predictors = [
  { name: 'fgas', arr: highV.map(g => g.fgas) },
  { name: 'log(MHI/L36)', arr: highV.map(g => g.logMHI_L36) },
  { name: 'logSigma0', arr: highV.map(g => g.logSigma0) },
];

const targets = [
  { name: 'mean RAR residual', arr: highV.map(g => g.meanRARresid) },
  { name: 'outer RAR residual', arr: highV.map(g => g.outerRARresid) },
  { name: 'inner RAR residual', arr: highV.map(g => g.innerRARresid) },
  { name: 'RAR scatter (sd)', arr: highV.map(g => g.rarScatter) },
];

const test1Results = {};
for (const tgt of targets) {
  const valid = tgt.arr.map((v, i) => isFinite(v) ? i : -1).filter(i => i >= 0);
  const yValid = valid.map(i => tgt.arr[i]);
  test1Results[tgt.name] = {};

  console.log('  Target: ' + tgt.name + ' (N=' + valid.length + ')');
  for (const pred of predictors) {
    const xValid = valid.map(i => pred.arr[i]);
    const r = pearsonR(xValid, yValid);
    const loo = valid.length >= 15 ? looR2(xValid, yValid) : NaN;
    console.log('    ' + pred.name + ': r=' + r.toFixed(3) + (isFinite(loo) ? ', LOO=' + loo.toFixed(3) : ''));
    test1Results[tgt.name][pred.name] = { r: parseFloat(r.toFixed(3)), loo: isFinite(loo) ? parseFloat(loo.toFixed(3)) : null };
  }
  console.log('');
}

// ============================================================
// TEST 2: OUTER vs INNER RAR residuals
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: OUTER vs INNER RAR RESIDUALS');
console.log('  Does gas-to-stellar balance predict OUTER RAR residuals');
console.log('  more than INNER ones? (connecting to our Phase 102+ result)');
console.log('══════════════════════════════════════════════════════════════\n');

const validOI = highV.filter(g => isFinite(g.outerRARresid) && isFinite(g.innerRARresid));

const r_outer_fgas = pearsonR(validOI.map(g => g.fgas), validOI.map(g => g.outerRARresid));
const r_inner_fgas = pearsonR(validOI.map(g => g.fgas), validOI.map(g => g.innerRARresid));
const r_outer_ratio = pearsonR(validOI.map(g => g.logMHI_L36), validOI.map(g => g.outerRARresid));
const r_inner_ratio = pearsonR(validOI.map(g => g.logMHI_L36), validOI.map(g => g.innerRARresid));

console.log('  N = ' + validOI.length + ' galaxies with both inner and outer RAR residuals\n');
console.log('  fgas vs outer RAR resid:  r=' + r_outer_fgas.toFixed(3));
console.log('  fgas vs inner RAR resid:  r=' + r_inner_fgas.toFixed(3));
console.log('  Outer/inner ratio:        ' + (Math.abs(r_outer_fgas) / Math.abs(r_inner_fgas)).toFixed(2) + 'x\n');
console.log('  log(MHI/L36) vs outer RAR resid:  r=' + r_outer_ratio.toFixed(3));
console.log('  log(MHI/L36) vs inner RAR resid:  r=' + r_inner_ratio.toFixed(3));
console.log('  Outer/inner ratio:                ' + (Math.abs(r_outer_ratio) / Math.abs(r_inner_ratio)).toFixed(2) + 'x\n');

const outerStronger = Math.abs(r_outer_fgas) > Math.abs(r_inner_fgas) + 0.05;
console.log('  Gas-to-stellar balance more predictive of OUTER RAR? ' + (outerStronger ? 'YES' : 'NO') + '\n');

// ============================================================
// TEST 3: RAR by gas-fraction terciles
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: RAR SHAPE BY GAS-FRACTION TERCILES');
console.log('  Do gas-rich vs gas-poor galaxies follow DIFFERENT RARs?');
console.log('══════════════════════════════════════════════════════════════\n');

const sortedByFgas = [...highV].sort((a, b) => a.fgas - b.fgas);
const tercileSize = Math.floor(sortedByFgas.length / 3);
const terciles = [
  { name: 'Low fgas (stellar-dominated)', galaxies: sortedByFgas.slice(0, tercileSize) },
  { name: 'Mid fgas', galaxies: sortedByFgas.slice(tercileSize, 2 * tercileSize) },
  { name: 'High fgas (gas-dominated)', galaxies: sortedByFgas.slice(2 * tercileSize) },
];

const tercileResults = [];
for (const t of terciles) {
  const names = new Set(t.galaxies.map(g => g.name));
  const tPts = allRARpoints.filter(p => names.has(p.name));

  const fgasRange = [Math.min(...t.galaxies.map(g => g.fgas)).toFixed(3), Math.max(...t.galaxies.map(g => g.fgas)).toFixed(3)];
  const meanResid = mean(t.galaxies.filter(g => isFinite(g.meanRARresid)).map(g => g.meanRARresid));
  const outerResid = mean(t.galaxies.filter(g => isFinite(g.outerRARresid)).map(g => g.outerRARresid));
  const rmsScatter = Math.sqrt(mean(tPts.map(p => {
    const rar_pred = p.log_gbar + Math.log10(1 / (1 - Math.exp(-Math.sqrt(Math.pow(10, p.log_gbar - Math.log10(1.2e-10))))));
    const resid = p.log_gobs - (isFinite(rar_pred) ? rar_pred : p.log_gbar);
    return resid ** 2;
  }).filter(v => isFinite(v))));

  console.log('  ' + t.name + ':');
  console.log('    N_gal=' + t.galaxies.length + ', N_pts=' + tPts.length + ', fgas=[' + fgasRange[0] + ', ' + fgasRange[1] + ']');
  console.log('    Mean RAR residual: ' + meanResid.toFixed(4));
  console.log('    Mean OUTER RAR residual: ' + outerResid.toFixed(4));
  console.log('    RMS scatter around RAR: ' + rmsScatter.toFixed(4));
  console.log('');

  tercileResults.push({
    name: t.name,
    nGal: t.galaxies.length,
    nPts: tPts.length,
    fgasRange,
    meanRARresid: parseFloat(meanResid.toFixed(4)),
    outerRARresid: parseFloat(outerResid.toFixed(4)),
    rmsScatter: parseFloat(rmsScatter.toFixed(4)),
  });
}

const residDiff = tercileResults[2].meanRARresid - tercileResults[0].meanRARresid;
const outerResidDiff = tercileResults[2].outerRARresid - tercileResults[0].outerRARresid;
console.log('  High-Low fgas RAR residual difference: ' + residDiff.toFixed(4));
console.log('  High-Low fgas OUTER RAR residual difference: ' + outerResidDiff.toFixed(4));
console.log('  Gas-rich galaxies sit ' + (residDiff > 0 ? 'ABOVE' : 'BELOW') + ' gas-poor on RAR\n');

// ============================================================
// TEST 4: RAR scatter vs gas-to-stellar balance
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 4: DOES GAS STATE PREDICT RAR SCATTER?');
console.log('  If gas-rich galaxies have more scatter around RAR,');
console.log('  it connects gas-to-stellar balance to RAR tightness.');
console.log('══════════════════════════════════════════════════════════════\n');

const validScatter = highV.filter(g => isFinite(g.rarScatter) && g.nRARpts >= 10);
if (validScatter.length >= 15) {
  const r_scatter_fgas = pearsonR(validScatter.map(g => g.fgas), validScatter.map(g => g.rarScatter));
  const r_scatter_ratio = pearsonR(validScatter.map(g => g.logMHI_L36), validScatter.map(g => g.rarScatter));
  const r_scatter_sig = pearsonR(validScatter.map(g => g.logSigma0), validScatter.map(g => g.rarScatter));

  console.log('  N = ' + validScatter.length + ' galaxies (>= 10 RAR points each)\n');
  console.log('  fgas vs RAR scatter:          r=' + r_scatter_fgas.toFixed(3));
  console.log('  log(MHI/L36) vs RAR scatter:  r=' + r_scatter_ratio.toFixed(3));
  console.log('  logSigma0 vs RAR scatter:     r=' + r_scatter_sig.toFixed(3));
  console.log('');
}

// ============================================================
// TEST 5: ACCELERATION REGIME SPLIT
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 5: ACCELERATION REGIME SPLIT');
console.log('  Does gas-to-stellar balance matter more in low-gbar regime?');
console.log('══════════════════════════════════════════════════════════════\n');

const highGbar = allRARpoints.filter(p => p.log_gbar > Math.log10(1.2e-10) && p.Vflat >= 70);
const lowGbar = allRARpoints.filter(p => p.log_gbar <= Math.log10(1.2e-10) && p.Vflat >= 70);

function rarResidual(p) {
  const rar_pred = p.log_gbar + Math.log10(1 / (1 - Math.exp(-Math.sqrt(Math.pow(10, p.log_gbar - Math.log10(1.2e-10))))));
  return p.log_gobs - (isFinite(rar_pred) ? rar_pred : p.log_gbar);
}

const highGbarResids = highGbar.map(p => ({ ...p, resid: rarResidual(p) })).filter(p => isFinite(p.resid));
const lowGbarResids = lowGbar.map(p => ({ ...p, resid: rarResidual(p) })).filter(p => isFinite(p.resid));

console.log('  High-gbar regime (gbar > a0): N_pts=' + highGbarResids.length);
console.log('    r(fgas, RAR resid) = ' + pearsonR(highGbarResids.map(p => p.fgas), highGbarResids.map(p => p.resid)).toFixed(3));
console.log('    r(log(MHI/L36), RAR resid) = ' + pearsonR(highGbarResids.map(p => p.logMHI_L36), highGbarResids.map(p => p.resid)).toFixed(3));
console.log('');

console.log('  Low-gbar regime (gbar <= a0): N_pts=' + lowGbarResids.length);
console.log('    r(fgas, RAR resid) = ' + pearsonR(lowGbarResids.map(p => p.fgas), lowGbarResids.map(p => p.resid)).toFixed(3));
console.log('    r(log(MHI/L36), RAR resid) = ' + pearsonR(lowGbarResids.map(p => p.logMHI_L36), lowGbarResids.map(p => p.resid)).toFixed(3));
console.log('');

// ============================================================
// TEST 6: logOMD vs RAR residuals — are they the same thing?
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 6: logOMD vs GALAXY-LEVEL RAR RESIDUALS');
console.log('  Are they measuring the same thing?');
console.log('══════════════════════════════════════════════════════════════\n');

const validBoth = highV.filter(g => isFinite(g.meanRARresid) && isFinite(g.outerRARresid));
const r_logOMD_rarResid = pearsonR(validBoth.map(g => g.logOMD), validBoth.map(g => g.meanRARresid));
const r_logOMD_outerRar = pearsonR(validBoth.map(g => g.logOMD), validBoth.map(g => g.outerRARresid));

console.log('  r(logOMD, mean RAR residual) = ' + r_logOMD_rarResid.toFixed(3));
console.log('  r(logOMD, outer RAR residual) = ' + r_logOMD_outerRar.toFixed(3));
console.log('  -> ' + (Math.abs(r_logOMD_outerRar) > 0.8 ? 'STRONGLY related — same underlying signal' : Math.abs(r_logOMD_outerRar) > 0.5 ? 'Moderately related — overlapping but not identical' : 'Weakly related — different information') + '\n');

// ============================================================
// VERDICT
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  OVERALL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const gasRARconnection = Math.abs(test1Results['mean RAR residual']['fgas'].r) > 0.3 ||
                         Math.abs(test1Results['outer RAR residual']['fgas'].r) > 0.3;
const outerSpecific = outerStronger;
const tercileSplit = Math.abs(residDiff) > 0.01;

let verdict;
if (gasRARconnection && outerSpecific && tercileSplit) {
  verdict = 'STRONG CONNECTION: Gas-to-stellar balance systematically predicts RAR residuals, especially in the outer regions. This connects our galaxy-level finding to the point-by-point RAR framework. Gas-rich galaxies deviate from the mean RAR differently than gas-poor ones.';
} else if (gasRARconnection && tercileSplit) {
  verdict = 'MODERATE CONNECTION: Gas-to-stellar balance correlates with RAR residuals. The signal exists but may not be outer-specific.';
} else if (tercileSplit) {
  verdict = 'WEAK CONNECTION: Gas-fraction terciles show different RAR residuals, but the correlation is weak at the galaxy level.';
} else {
  verdict = 'NO CONNECTION: Gas-to-stellar balance does not systematically predict RAR residuals. Our finding (logOMD vs gas state) and the RAR operate on independent axes.';
}

console.log('  ' + verdict);

const output = {
  phase: 113,
  title: 'RAR-Gas Connection: Does gas-to-stellar balance predict RAR residuals?',
  nHighV: highV.length,
  nRARpoints: allRARpoints.length,
  test1_galaxyLevel: test1Results,
  test2_outerVsInner: {
    fgas: { outerR: parseFloat(r_outer_fgas.toFixed(3)), innerR: parseFloat(r_inner_fgas.toFixed(3)) },
    logMHI_L36: { outerR: parseFloat(r_outer_ratio.toFixed(3)), innerR: parseFloat(r_inner_ratio.toFixed(3)) },
    outerStronger
  },
  test3_terciles: tercileResults,
  test6_logOMD_vs_RAR: {
    r_logOMD_meanRAR: parseFloat(r_logOMD_rarResid.toFixed(3)),
    r_logOMD_outerRAR: parseFloat(r_logOMD_outerRar.toFixed(3)),
  },
  verdict,
};

const outPath = path.join(__dirname, '..', 'public', 'phase113-rar-gas-connection.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);

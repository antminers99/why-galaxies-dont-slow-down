#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 114: DIVERSITY CONNECTION');
console.log('  Ref: Oman+ 2015, Santos-Santos+ 2020, Marasco+ 2020');
console.log('');
console.log('  Rotation curve SHAPE diversity is a major open problem.');
console.log('  Question: Can gas-to-stellar balance explain some of the');
console.log('  diversity in rotation curve shapes?');
console.log('');
console.log('  We measure shape metrics: rise steepness, flatness,');
console.log('  outer slope, Vmax/Vflat ratio, and test if gas state');
console.log('  predicts shape diversity.');
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
  if (rcPoints.length < 8) continue;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logMHI_L36 = Math.log10(t1.MHI) - Math.log10(t1.L36);
  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
  const logL36 = Math.log10(t1.L36);
  if (!isFinite(fgas) || !isFinite(logMHI_L36)) continue;

  const sorted = [...rcPoints].filter(p => p.rad > 0 && p.vobs > 0).sort((a, b) => a.rad - b.rad);
  if (sorted.length < 8) continue;

  const Vmax = Math.max(...sorted.map(p => p.vobs));
  const rMax = sorted.find(p => p.vobs === Vmax)?.rad || sorted[sorted.length - 1].rad;
  const rLast = sorted[sorted.length - 1].rad;
  const Vlast = sorted[sorted.length - 1].vobs;

  // Shape metric 1: Rise steepness — how fast V reaches Vmax (normalized to Rdisk)
  const v50idx = sorted.findIndex(p => p.vobs >= 0.5 * Vmax);
  const r50 = v50idx >= 0 ? sorted[v50idx].rad : sorted[0].rad;
  const riseSteepness = t1.Rdisk > 0 ? r50 / t1.Rdisk : NaN;

  // Shape metric 2: Outer slope — (Vlast - Vmax) / Vmax
  const outerSlope = (Vlast - Vmax) / Vmax;

  // Shape metric 3: Vmax / Vflat ratio (from Verheijen 2001)
  const VmaxVflat = Vmax / t1.Vflat;

  // Shape metric 4: Flatness — std of V / mean(V) in outer half
  const half = Math.floor(sorted.length / 2);
  const outerV = sorted.slice(half).map(p => p.vobs);
  const flatness = outerV.length > 2 ? sd(outerV) / mean(outerV) : NaN;

  // Shape metric 5: Inner concentration — V at 1Rd / Vmax
  const r1Rd = sorted.find(p => p.rad >= t1.Rdisk);
  const V1Rd = r1Rd ? r1Rd.vobs : sorted[Math.min(2, sorted.length - 1)].vobs;
  const concentration = V1Rd / Vmax;

  // Shape metric 6: Baryon dominance fraction — fraction of points where Vbar > 0.7 * Vobs
  const baryonDom = sorted.filter(p => {
    const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                   UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                   p.vgas * Math.abs(p.vgas);
    const vBar = Math.sqrt(Math.max(vBarSq, 0.01));
    return vBar > 0.7 * p.vobs;
  }).length / sorted.length;

  // logOMD (our standard target)
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
    name, fgas, logMHI_L36, logSigma0, logL36, Vflat: t1.Vflat, logOMD,
    riseSteepness: isFinite(riseSteepness) ? riseSteepness : NaN,
    outerSlope: isFinite(outerSlope) ? outerSlope : NaN,
    VmaxVflat: isFinite(VmaxVflat) ? VmaxVflat : NaN,
    flatness: isFinite(flatness) ? flatness : NaN,
    concentration: isFinite(concentration) ? concentration : NaN,
    baryonDom: isFinite(baryonDom) ? baryonDom : NaN,
  });
}

const highV = galaxies.filter(g => g.Vflat >= 70);
const all = galaxies;

console.log('  Total galaxies: N=' + all.length);
console.log('  High-Vflat (>=70): N=' + highV.length + '\n');

// ============================================================
// TEST 1: Gas state vs shape metrics (high-Vflat sample)
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: GAS STATE vs ROTATION CURVE SHAPE METRICS');
console.log('  (High-Vflat sample, N=' + highV.length + ')');
console.log('══════════════════════════════════════════════════════════════\n');

const shapeMetrics = [
  { name: 'outerSlope', label: 'Outer slope (Vlast-Vmax)/Vmax' },
  { name: 'VmaxVflat', label: 'Vmax/Vflat ratio (Verheijen 2001)' },
  { name: 'flatness', label: 'Flatness (sd/mean outer V)' },
  { name: 'concentration', label: 'Inner concentration (V@1Rd/Vmax)' },
  { name: 'riseSteepness', label: 'Rise steepness (r50/Rdisk)' },
  { name: 'baryonDom', label: 'Baryon dominance fraction' },
];

const preds = [
  { name: 'fgas', get: g => g.fgas },
  { name: 'log(MHI/L36)', get: g => g.logMHI_L36 },
  { name: 'logSigma0', get: g => g.logSigma0 },
  { name: 'logL36', get: g => g.logL36 },
];

const test1Results = {};

for (const sm of shapeMetrics) {
  const valid = highV.filter(g => isFinite(g[sm.name]));
  console.log('  ' + sm.label + ' (N=' + valid.length + ')');
  test1Results[sm.name] = {};

  for (const pred of preds) {
    const x = valid.map(pred.get);
    const y = valid.map(g => g[sm.name]);
    const r = pearsonR(x, y);
    const loo = valid.length >= 15 ? looR2(x, y) : NaN;
    const marker = Math.abs(r) > 0.4 ? ' <<<' : Math.abs(r) > 0.3 ? ' <<' : '';
    console.log('    ' + pred.name + ': r=' + r.toFixed(3) + (isFinite(loo) ? ', LOO=' + loo.toFixed(3) : '') + marker);
    test1Results[sm.name][pred.name] = { r: parseFloat(r.toFixed(3)), loo: isFinite(loo) ? parseFloat(loo.toFixed(3)) : null };
  }
  console.log('');
}

// ============================================================
// TEST 2: SAME TESTS ON FULL SAMPLE (including dwarfs)
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: GAS STATE vs SHAPE — FULL SAMPLE (N=' + all.length + ')');
console.log('  Including dwarfs: does the pattern change?');
console.log('══════════════════════════════════════════════════════════════\n');

const test2Results = {};
for (const sm of shapeMetrics) {
  const valid = all.filter(g => isFinite(g[sm.name]));
  console.log('  ' + sm.label + ' (N=' + valid.length + ')');
  test2Results[sm.name] = {};

  for (const pred of ['fgas', 'log(MHI/L36)']) {
    const predFn = pred === 'fgas' ? (g => g.fgas) : (g => g.logMHI_L36);
    const x = valid.map(predFn);
    const y = valid.map(g => g[sm.name]);
    const r = pearsonR(x, y);
    const marker = Math.abs(r) > 0.4 ? ' <<<' : Math.abs(r) > 0.3 ? ' <<' : '';
    console.log('    ' + pred + ': r=' + r.toFixed(3) + marker);
    test2Results[sm.name][pred] = { r: parseFloat(r.toFixed(3)) };
  }
  console.log('');
}

// ============================================================
// TEST 3: PARTIAL CORRELATIONS — shape metrics after controlling Vflat
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: PARTIAL CORRELATIONS (controlling Vflat)');
console.log('  Does gas state predict shape BEYOND mass/velocity?');
console.log('══════════════════════════════════════════════════════════════\n');

const test3Results = {};
for (const sm of shapeMetrics) {
  const valid = highV.filter(g => isFinite(g[sm.name]));
  if (valid.length < 15) continue;

  const y = valid.map(g => g[sm.name]);
  const vflat = valid.map(g => Math.log10(g.Vflat));

  const pr_fgas = partialR(valid.map(g => g.fgas), y, [vflat]);
  const pr_ratio = partialR(valid.map(g => g.logMHI_L36), y, [vflat]);
  const pr_sigma = partialR(valid.map(g => g.logSigma0), y, [vflat]);

  console.log('  ' + sm.label);
  console.log('    partial r(fgas | Vflat) = ' + pr_fgas.toFixed(3));
  console.log('    partial r(log(MHI/L36) | Vflat) = ' + pr_ratio.toFixed(3));
  console.log('    partial r(logSigma0 | Vflat) = ' + pr_sigma.toFixed(3));

  test3Results[sm.name] = {
    pr_fgas: parseFloat(pr_fgas.toFixed(3)),
    pr_logMHI_L36: parseFloat(pr_ratio.toFixed(3)),
    pr_logSigma0: parseFloat(pr_sigma.toFixed(3)),
  };
  console.log('');
}

// ============================================================
// TEST 4: DIVERSITY METRIC — total shape variance explained
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 4: SHAPE DIVERSITY INDEX');
console.log('  Combined shape diversity metric and how much gas explains');
console.log('══════════════════════════════════════════════════════════════\n');

const validAll = highV.filter(g =>
  isFinite(g.outerSlope) && isFinite(g.flatness) && isFinite(g.concentration) && isFinite(g.VmaxVflat));

if (validAll.length >= 15) {
  const metrics = ['outerSlope', 'flatness', 'concentration', 'VmaxVflat'];
  const normalized = metrics.map(m => {
    const vals = validAll.map(g => g[m]);
    const mu = mean(vals), s = sd(vals);
    return vals.map(v => s > 0 ? (v - mu) / s : 0);
  });

  const diversityIndex = validAll.map((_, i) =>
    Math.sqrt(metrics.reduce((s, _, j) => s + normalized[j][i] ** 2, 0) / metrics.length));

  const r_div_fgas = pearsonR(validAll.map(g => g.fgas), diversityIndex);
  const r_div_ratio = pearsonR(validAll.map(g => g.logMHI_L36), diversityIndex);
  const r_div_sigma = pearsonR(validAll.map(g => g.logSigma0), diversityIndex);
  const r_div_logOMD = pearsonR(validAll.map(g => g.logOMD), diversityIndex);

  console.log('  Combined shape diversity index (RMS of normalized shape metrics)');
  console.log('  N = ' + validAll.length + '\n');
  console.log('  r(fgas, diversity) = ' + r_div_fgas.toFixed(3));
  console.log('  r(log(MHI/L36), diversity) = ' + r_div_ratio.toFixed(3));
  console.log('  r(logSigma0, diversity) = ' + r_div_sigma.toFixed(3));
  console.log('  r(logOMD, diversity) = ' + r_div_logOMD.toFixed(3));
  console.log('');

  // Tercile comparison
  const sortedByFgas = [...validAll].sort((a, b) => a.fgas - b.fgas);
  const t3 = Math.floor(sortedByFgas.length / 3);
  const tercLabels = ['Low fgas', 'Mid fgas', 'High fgas'];
  const tercGroups = [sortedByFgas.slice(0, t3), sortedByFgas.slice(t3, 2 * t3), sortedByFgas.slice(2 * t3)];

  console.log('  Shape diversity by fgas tercile:');
  for (let t = 0; t < 3; t++) {
    const tIdx = tercGroups[t].map(g => validAll.indexOf(g));
    const tDiv = tIdx.map(i => diversityIndex[i]);
    console.log('    ' + tercLabels[t] + ' (N=' + tercGroups[t].length + '): mean diversity = ' + mean(tDiv).toFixed(3) + ', sd = ' + sd(tDiv).toFixed(3));
  }
  console.log('');
}

// ============================================================
// TEST 5: BARYON DOMINANCE FRACTION — key diversity link
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 5: BARYON DOMINANCE FRACTION');
console.log('  Fraction of RC points where Vbar > 0.7*Vobs');
console.log('  (connects to Santos-Santos+ 2020 "baryonic clues")');
console.log('══════════════════════════════════════════════════════════════\n');

const validBD = highV.filter(g => isFinite(g.baryonDom));
if (validBD.length >= 15) {
  const r_bd_fgas = pearsonR(validBD.map(g => g.fgas), validBD.map(g => g.baryonDom));
  const r_bd_ratio = pearsonR(validBD.map(g => g.logMHI_L36), validBD.map(g => g.baryonDom));
  const r_bd_logOMD = pearsonR(validBD.map(g => g.logOMD), validBD.map(g => g.baryonDom));
  const loo_bd_fgas = looR2(validBD.map(g => g.fgas), validBD.map(g => g.baryonDom));

  console.log('  N = ' + validBD.length);
  console.log('  r(fgas, baryonDom) = ' + r_bd_fgas.toFixed(3));
  console.log('  r(log(MHI/L36), baryonDom) = ' + r_bd_ratio.toFixed(3));
  console.log('  r(logOMD, baryonDom) = ' + r_bd_logOMD.toFixed(3));
  console.log('  LOO(fgas -> baryonDom) = ' + loo_bd_fgas.toFixed(3));
  console.log('');
}

// ============================================================
// VERDICT
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  OVERALL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const strongShapeLinks = [];
for (const sm of shapeMetrics) {
  for (const pred of preds) {
    const result = test1Results[sm.name]?.[pred.name];
    if (result && Math.abs(result.r) > 0.3) {
      strongShapeLinks.push(sm.label + ' <-> ' + pred.name + ' (r=' + result.r + ')');
    }
  }
}

if (strongShapeLinks.length > 3) {
  console.log('  STRONG: Gas-to-stellar balance predicts multiple shape metrics.');
} else if (strongShapeLinks.length > 0) {
  console.log('  MODERATE: Gas state predicts some shape metrics:');
} else {
  console.log('  WEAK: Gas state is not a strong predictor of individual shape metrics.');
}
for (const link of strongShapeLinks) {
  console.log('    ' + link);
}

const output = {
  phase: 114,
  title: 'Diversity Connection: Gas-to-stellar balance and rotation curve shape diversity',
  nHighV: highV.length,
  nAll: all.length,
  test1_highV_shapeVsGas: test1Results,
  test2_fullSample: test2Results,
  test3_partial: test3Results,
  strongShapeLinks,
};

const outPath = path.join(__dirname, '..', 'public', 'phase114-diversity-connection.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);

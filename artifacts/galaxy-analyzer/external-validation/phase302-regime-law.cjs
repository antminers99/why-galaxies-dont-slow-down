#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const pub = p => path.join(__dirname, '..', 'public', p);

console.log('======================================================================');
console.log('  PHASE 302: REGIME LAW');
console.log('  Does the hidden physics inside VfResid turn on with Vflat,');
console.log('  or does only its observability improve?');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(pub('stage-A-master-table.json'), 'utf8'));
const sparcTable = JSON.parse(fs.readFileSync(pub('sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(pub('sparc-results.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(pub('phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const extDataset = JSON.parse(fs.readFileSync(pub('phase200-external-dataset.json'), 'utf8'));
const salvageData = JSON.parse(fs.readFileSync(pub('phase300-sample-salvage.json'), 'utf8'));
const p301 = JSON.parse(fs.readFileSync(pub('phase301-vfresid-drivers.json'), 'utf8'));

const sparcMap = {}; sparcTable.forEach(s => { sparcMap[s.name] = s; });
const resMap = {}; sparcResults.perGalaxy.forEach(g => { if (g.rawName) resMap[g.rawName] = g; });
const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : 0; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function pearsonR(x, y) {
  const n = x.length; if (n < 3) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}
function solveLinear(A, b) {
  const n = A.length;
  const M = A.map((r, i) => [...r, b[i]]);
  for (let i = 0; i < n; i++) {
    let mx = i;
    for (let j = i + 1; j < n; j++) if (Math.abs(M[j][i]) > Math.abs(M[mx][i])) mx = j;
    [M[i], M[mx]] = [M[mx], M[i]];
    if (Math.abs(M[i][i]) < 1e-15) continue;
    for (let j = i + 1; j < n; j++) {
      const f = M[j][i] / M[i][i];
      for (let k = i; k <= n; k++) M[j][k] -= f * M[i][k];
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}
function olsFull(Y, X) {
  const n = Y.length, p = X[0].length + 1;
  const Xa = X.map(r => [1, ...r]);
  const XtX = Array.from({ length: p }, () => new Array(p).fill(0));
  const XtY = new Array(p).fill(0);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < p; j++) {
      XtY[j] += Xa[i][j] * Y[i];
      for (let l = 0; l < p; l++) XtX[j][l] += Xa[i][j] * Xa[i][l];
    }
  const beta = solveLinear(XtX, XtY);
  let sse = 0, sst = 0;
  const my = mean(Y);
  for (let i = 0; i < n; i++) {
    const pred = Xa[i].reduce((s, x, j) => s + x * beta[j], 0);
    sse += (Y[i] - pred) ** 2;
    sst += (Y[i] - my) ** 2;
  }
  return { beta, r2: sst > 0 ? 1 - sse / sst : 0, rmse: Math.sqrt(sse / n) };
}
function looR2(Y, X) {
  const n = Y.length;
  let press = 0;
  const my = mean(Y);
  const sst = Y.reduce((s, y) => s + (y - my) ** 2, 0);
  if (sst === 0) return 0;
  for (let i = 0; i < n; i++) {
    const Yt = Y.filter((_, j) => j !== i);
    const Xt = X.filter((_, j) => j !== i);
    try {
      const m = olsFull(Yt, Xt);
      const xi = [1, ...X[i]];
      const pred = xi.reduce((s, x, j) => s + x * m.beta[j], 0);
      press += (Y[i] - pred) ** 2;
    } catch (e) { press += (Y[i] - my) ** 2; }
  }
  return 1 - press / sst;
}
function gapPct(Y, X) {
  const n = Y.length;
  const my = mean(Y);
  const sst = Y.reduce((s, y) => s + (y - my) ** 2, 0);
  if (sst === 0) return 0;
  let press = 0;
  for (let i = 0; i < n; i++) {
    const Yt = Y.filter((_, j) => j !== i);
    const Xt = X.filter((_, j) => j !== i);
    try {
      const m = olsFull(Yt, Xt);
      const pred = [1, ...X[i]].reduce((s, x, j) => s + x * m.beta[j], 0);
      press += (Y[i] - pred) ** 2;
    } catch (e) { press += (Y[i] - my) ** 2; }
  }
  return (1 - press / sst) * 100;
}
function bootstrapSignStability(Y, X, nBoot) {
  const n = Y.length;
  const p = X[0].length + 1;
  const signs = Array.from({ length: p }, () => ({ pos: 0, neg: 0 }));
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(Math.random() * n));
    try {
      const m = olsFull(idx.map(i => Y[i]), idx.map(i => X[i]));
      for (let j = 0; j < p; j++) {
        if (m.beta[j] > 0) signs[j].pos++; else signs[j].neg++;
      }
    } catch (e) {}
  }
  return signs.map(s => Math.max(s.pos, s.neg) / (s.pos + s.neg) * 100);
}

const pubGals = stageA.galaxies.filter(g => pubNames.has(g.name));
const structTrainData = pubGals.map(g => {
  const s = sparcMap[g.name];
  return {
    logVflat: Math.log10(s.Vflat),
    logMbar: Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9),
    logL36: Math.log10(Math.max(s.L36, 0.01)),
    logRdisk: Math.log10(s.Rdisk),
    morphT: s.T ?? 5
  };
}).filter(g => isFinite(g.logVflat) && isFinite(g.logMbar));
const structModel = olsFull(structTrainData.map(g => g.logVflat), structTrainData.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]));

function buildGalaxy(g, s, r, td, source) {
  const logMbar = Math.log10(s.L36 * 0.5 * 1e9 + (g.logMHI ? Math.pow(10, g.logMHI) : s.MHI) * 1.33 * 1e9);
  const logVflat = Math.log10(g.Vflat || s.Vflat);
  const logL36 = Math.log10(Math.max(s.L36, 0.01));
  const logRdisk = Math.log10(s.Rdisk);
  const predLogVflat = [1, logMbar, logL36, logRdisk, s.T ?? 5].reduce((sum, x, j) => sum + x * structModel.beta[j], 0);
  const VfResid = logVflat - predLogVflat;
  const haloK = r && r.models && r.models.dark_halo_linear ? Math.log10(Math.max(r.models.dark_halo_linear.k, 1)) : null;
  const lhOuter = r && r.models && r.models.log_halo ? r.models.log_halo.outerImprovement : null;
  const mondImprove = r && r.models && r.models.mond ? r.models.mond.improvementVsNewton : null;
  return {
    name: g.name, VfResid, logA0: g.logA0,
    haloK, lhOuter, mondImprove,
    logMHI: g.logMHI ?? Math.log10(s.MHI),
    logMbar, logMhost: g.logMhost ?? (td ? td.logMhost : null),
    logMeanRun: g.logMeanRun, logSigma0: g.logSigma0,
    logL36, logRdisk, morphT: s.T ?? 5,
    envCode: g.envCode, Vflat: g.Vflat || s.Vflat, logVflat,
    Q: g.Q ?? s.Q, source
  };
}

const internal = pubGals.map(g => {
  const s = sparcMap[g.name]; const r = resMap[g.name]; const td = tdMap[g.name];
  return buildGalaxy(g, s, r, td, 'internal');
}).filter(g => isFinite(g.VfResid));

const extOriginal = extDataset.galaxies.map(g => {
  const s = sparcMap[g.name]; const r = resMap[g.name];
  return { ...g, haloK: g.haloK, lhOuter: g.lhOuter, mondImprove: r && r.models && r.models.mond ? r.models.mond.improvementVsNewton : null, source: 'external' };
});

const extSalvaged = (salvageData.salvagedGalaxies || []).map(g => {
  const s = sparcMap[g.name]; const r = resMap[g.name];
  return { ...g, haloK: g.haloK, lhOuter: g.lhOuter, mondImprove: r && r.models && r.models.mond ? r.models.mond.improvementVsNewton : null, source: 'salvaged' };
});

const extAll = [...extOriginal, ...extSalvaged];
const pooledOriginal = [...internal, ...extOriginal];
const pooledAll = [...internal, ...extAll];

console.log('  Internal: N=' + internal.length);
console.log('  External original: N=' + extOriginal.length);
console.log('  External salvaged: N=' + extSalvaged.length);
console.log('  External all: N=' + extAll.length);
console.log('  Pooled (int+ext_orig): N=' + pooledOriginal.length);
console.log('  Pooled (int+ext_all): N=' + pooledAll.length);

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: ACTIVATION SHAPE');
console.log('  Is the regime effect sharp, smooth, or two-population?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function activationProfile(data, label) {
  const sorted = [...data].sort((a, b) => a.Vflat - b.Vflat);
  const thresholds = [50, 70, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 200];
  console.log('  === ' + label + ' (N=' + data.length + ') ===');

  console.log('\n  A) Cumulative threshold (Vflat >= threshold):');
  console.log('  ' + 'Thresh'.padStart(6) + '   N' + ' r(VfR,a0)'.padStart(11) + ' VfR coeff'.padStart(10) +
    ' signStab'.padStart(9) + '  |VfR|_mean');
  console.log('  ' + '-'.repeat(65));

  const cumResults = [];
  for (const th of thresholds) {
    const sub = data.filter(g => g.Vflat >= th);
    if (sub.length < 5) continue;
    const r = pearsonR(sub.map(g => g.logA0), sub.map(g => g.VfResid));
    const model = olsFull(sub.map(g => g.logA0), sub.map(g => [g.VfResid]));
    const sign = bootstrapSignStability(sub.map(g => g.logA0), sub.map(g => [g.VfResid]), 500);
    const absVfR = mean(sub.map(g => Math.abs(g.VfResid)));
    cumResults.push({ threshold: th, n: sub.length, r: isNaN(r) ? 0 : r, coeff: model.beta[1], signStab: sign[1], absVfR });
    console.log('  ' + String(th).padStart(6) + String(sub.length).padStart(4) +
      (isNaN(r) ? '       N/A' : r.toFixed(3).padStart(11)) +
      model.beta[1].toFixed(3).padStart(10) +
      sign[1].toFixed(0).padStart(8) + '%' +
      absVfR.toFixed(4).padStart(11));
  }

  console.log('\n  B) Sliding windows (width=40 km/s):');
  console.log('  ' + 'Center'.padStart(6) + '   N' + ' r(VfR,a0)'.padStart(11) + ' |VfR|_mean'.padStart(11) + ' sd(logA0)'.padStart(10));
  console.log('  ' + '-'.repeat(50));

  const windowResults = [];
  for (let center = 60; center <= 260; center += 20) {
    const lo = center - 20, hi = center + 20;
    const sub = data.filter(g => g.Vflat >= lo && g.Vflat < hi);
    if (sub.length < 4) continue;
    const r = pearsonR(sub.map(g => g.logA0), sub.map(g => g.VfResid));
    const absVfR = mean(sub.map(g => Math.abs(g.VfResid)));
    const sdA0 = sd(sub.map(g => g.logA0));
    windowResults.push({ center, n: sub.length, r: isNaN(r) ? 0 : r, absVfR, sdA0 });
    console.log('  ' + String(center).padStart(6) + String(sub.length).padStart(4) +
      (isNaN(r) ? '       N/A' : r.toFixed(3).padStart(11)) +
      absVfR.toFixed(4).padStart(11) +
      sdA0.toFixed(4).padStart(10));
  }

  console.log('\n  C) Regime bins:');
  const bins = [
    { label: 'Very-low (<70)', filter: g => g.Vflat < 70 },
    { label: 'Low (70-120)', filter: g => g.Vflat >= 70 && g.Vflat < 120 },
    { label: 'Mid (120-160)', filter: g => g.Vflat >= 120 && g.Vflat < 160 },
    { label: 'High (160-200)', filter: g => g.Vflat >= 160 && g.Vflat < 200 },
    { label: 'Very-high (>=200)', filter: g => g.Vflat >= 200 },
  ];
  console.log('  ' + 'Bin'.padEnd(22) + '  N' + ' r(VfR,a0)'.padStart(11) + '  VfR_sd'.padStart(9) + '  a0_sd'.padStart(8));
  const binResults = [];
  for (const bin of bins) {
    const sub = data.filter(bin.filter);
    if (sub.length < 3) { console.log('  ' + bin.label.padEnd(22) + String(sub.length).padStart(3) + '  (too few)'); continue; }
    const r = pearsonR(sub.map(g => g.logA0), sub.map(g => g.VfResid));
    const vfSd = sd(sub.map(g => g.VfResid));
    const a0Sd = sd(sub.map(g => g.logA0));
    binResults.push({ label: bin.label, n: sub.length, r: isNaN(r) ? 0 : r, vfSd, a0Sd });
    console.log('  ' + bin.label.padEnd(22) + String(sub.length).padStart(3) +
      (isNaN(r) ? '       N/A' : r.toFixed(3).padStart(11)) +
      vfSd.toFixed(4).padStart(9) +
      a0Sd.toFixed(4).padStart(8));
  }

  return { cumResults, windowResults, binResults };
}

const actInt = activationProfile(internal, 'INTERNAL');
console.log('');
const actExtOrig = activationProfile(extOriginal, 'EXTERNAL ORIGINAL');
console.log('');
const actExtAll = activationProfile(extAll, 'EXTERNAL ALL (orig+salvaged)');
console.log('');
const actPoolOrig = activationProfile(pooledOriginal, 'POOLED (int+ext_orig)');

console.log('\n  D) ACTIVATION SHAPE ASSESSMENT:');
function assessShape(cumResults, label) {
  const filtered = cumResults.filter(c => c.n >= 5);
  if (filtered.length < 3) { console.log('  ' + label + ': insufficient data'); return 'INSUFFICIENT'; }
  const rVals = filtered.map(c => c.r);
  const isMonotonic = rVals.every((r, i) => i === 0 || r >= rVals[i - 1] - 0.05);
  const maxR = Math.max(...rVals);
  const minR = Math.min(...rVals.slice(0, Math.min(3, rVals.length)));
  const range = maxR - minR;

  const lowR = filtered.filter(c => c.threshold <= 100).map(c => c.r);
  const highR = filtered.filter(c => c.threshold >= 140).map(c => c.r);
  const avgLow = lowR.length ? mean(lowR) : 0;
  const avgHigh = highR.length ? mean(highR) : 0;
  const jump = avgHigh - avgLow;

  let shape = 'UNKNOWN';
  if (jump > 0.3 && range > 0.3) shape = 'SHARP_ACTIVATION';
  else if (jump > 0.15 && isMonotonic) shape = 'SMOOTH_RAMP';
  else if (jump > 0.1) shape = 'GRADUAL_STRENGTHENING';
  else if (range < 0.1) shape = 'FLAT_THROUGHOUT';
  else shape = 'NOISY';

  console.log('  ' + label + ': ' + shape +
    ' (low-V avg r=' + avgLow.toFixed(3) + ', high-V avg r=' + avgHigh.toFixed(3) +
    ', jump=' + jump.toFixed(3) + ', range=' + range.toFixed(3) + ')');
  return shape;
}

const shapeInt = assessShape(actInt.cumResults, 'Internal');
const shapeExtOrig = assessShape(actExtOrig.cumResults, 'External original');
const shapeExtAll = assessShape(actExtAll.cumResults, 'External all');
const shapePool = assessShape(actPoolOrig.cumResults, 'Pooled');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: DRIVER REGIME SPLIT');
console.log('  Does haloK strengthen with regime, or does the residual rise faster?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function driverRegimeSplit(data, label) {
  const regimes = [
    { label: 'Low (<120)', filter: g => g.Vflat < 120 },
    { label: 'High (>=120)', filter: g => g.Vflat >= 120 },
    { label: 'Very-high (>=180)', filter: g => g.Vflat >= 180 },
  ];
  const drivers = ['haloK', 'lhOuter', 'mondImprove', 'logMHI', 'logMbar', 'logMeanRun', 'logRdisk', 'logSBeff'];

  console.log('  === ' + label + ' ===');
  console.log('  ' + 'Driver'.padEnd(16) + regimes.map(r => (r.label + ' r').padStart(14)).join(''));
  console.log('  ' + '-'.repeat(60));

  const results = {};
  for (const d of drivers) {
    const row = [];
    for (const reg of regimes) {
      const sub = data.filter(reg.filter).filter(g => g[d] !== null && g[d] !== undefined && isFinite(g[d]));
      if (sub.length < 4) { row.push({ r: NaN, n: sub.length }); continue; }
      const r = pearsonR(sub.map(g => g.VfResid), sub.map(g => g[d]));
      row.push({ r, n: sub.length });
    }
    results[d] = row;
    console.log('  ' + d.padEnd(16) + row.map(v =>
      (isNaN(v.r) ? ('N/A(' + v.n + ')') : (v.r.toFixed(3) + '(' + v.n + ')')).padStart(14)
    ).join(''));
  }
  return results;
}

const drsInt = driverRegimeSplit(internal, 'INTERNAL DRIVERS by REGIME');
console.log('');
const drsExtOrig = driverRegimeSplit(extOriginal, 'EXTERNAL ORIGINAL DRIVERS by REGIME');
console.log('');
const drsPool = driverRegimeSplit(pooledOriginal, 'POOLED (int+ext_orig) DRIVERS by REGIME');

console.log('\n  HALO K REGIME PATTERN:');
function haloKRegimeDetail(data, label) {
  const regimes = [
    { label: '<80', lo: 0, hi: 80 },
    { label: '80-120', lo: 80, hi: 120 },
    { label: '120-160', lo: 120, hi: 160 },
    { label: '160-220', lo: 160, hi: 220 },
    { label: '>=220', lo: 220, hi: 9999 },
  ];
  console.log('  ' + label + ':');
  for (const reg of regimes) {
    const sub = data.filter(g => g.Vflat >= reg.lo && g.Vflat < reg.hi && g.haloK !== null && isFinite(g.haloK));
    if (sub.length < 3) { console.log('    ' + reg.label + ': N=' + sub.length + ' (too few)'); continue; }
    const rHK = pearsonR(sub.map(g => g.VfResid), sub.map(g => g.haloK));
    const meanHK = mean(sub.map(g => g.haloK));
    const sdVfR = sd(sub.map(g => g.VfResid));
    console.log('    ' + reg.label.padEnd(8) + ' N=' + String(sub.length).padStart(3) +
      ' r(haloK,VfR)=' + (isNaN(rHK) ? 'N/A' : rHK.toFixed(3)) +
      ' mean(haloK)=' + meanHK.toFixed(3) +
      ' sd(VfR)=' + sdVfR.toFixed(4));
  }
}

haloKRegimeDetail(internal, 'Internal');
haloKRegimeDetail(pooledOriginal, 'Pooled');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: RESIDUAL REGIME LAW');
console.log('  Does the unexplained VfResid also activate with Vflat?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function residualRegimeTest(data, features, label) {
  const valid = data.filter(g => features.every(f => g[f] !== null && g[f] !== undefined && isFinite(g[f])) && isFinite(g.logA0));
  if (valid.length < features.length + 5) {
    console.log('  ' + label + ': SKIP (N=' + valid.length + ')');
    return null;
  }

  const Y = valid.map(g => g.VfResid);
  const X = valid.map(g => features.map(f => g[f]));
  const model = olsFull(Y, X);
  const residuals = Y.map((y, i) => y - [1, ...X[i]].reduce((s, x, j) => s + x * model.beta[j], 0));

  valid.forEach((g, i) => { g._residVfR = residuals[i]; });

  console.log('  === ' + label + ' (N=' + valid.length + ', features: ' + features.join('+') + ') ===');
  console.log('  Model R² = ' + model.r2.toFixed(3));

  const regimes = [
    { label: 'Low (<120)', filter: g => g.Vflat < 120 },
    { label: 'High (>=120)', filter: g => g.Vflat >= 120 },
    { label: 'Very-high (>=180)', filter: g => g.Vflat >= 180 },
  ];

  console.log('\n  ' + 'Regime'.padEnd(22) + '  N' + ' r(resid,a0)'.padStart(13) + ' sd(resid)'.padStart(11) + '  |resid|_mean');
  console.log('  ' + '-'.repeat(65));

  const regimeResults = [];
  for (const reg of regimes) {
    const sub = valid.filter(reg.filter);
    if (sub.length < 3) { console.log('  ' + reg.label.padEnd(22) + String(sub.length).padStart(3) + '  (too few)'); continue; }
    const rRA0 = pearsonR(sub.map(g => g.logA0), sub.map(g => g._residVfR));
    const sdResid = sd(sub.map(g => g._residVfR));
    const absResid = mean(sub.map(g => Math.abs(g._residVfR)));
    regimeResults.push({ regime: reg.label, n: sub.length, rResidA0: isNaN(rRA0) ? 0 : rRA0, sdResid, absResid });
    console.log('  ' + reg.label.padEnd(22) + String(sub.length).padStart(3) +
      (isNaN(rRA0) ? '         N/A' : rRA0.toFixed(3).padStart(13)) +
      sdResid.toFixed(4).padStart(11) +
      absResid.toFixed(4).padStart(13));
  }

  const thresholds = [70, 90, 110, 130, 150, 170, 190];
  console.log('\n  Cumulative r(unexplained_VfResid, a0) >= threshold:');
  const cumResid = [];
  for (const th of thresholds) {
    const sub = valid.filter(g => g.Vflat >= th);
    if (sub.length < 4) continue;
    const rRA0 = pearsonR(sub.map(g => g.logA0), sub.map(g => g._residVfR));
    cumResid.push({ threshold: th, n: sub.length, r: isNaN(rRA0) ? 0 : rRA0 });
    console.log('    Vflat >= ' + th + ': N=' + sub.length + '  r=' + (isNaN(rRA0) ? 'N/A' : rRA0.toFixed(3)));
  }

  return { regimeResults, cumResid, modelR2: model.r2 };
}

const residIntHaloK = residualRegimeTest(internal, ['haloK'], 'Internal (haloK only)');
console.log('');
const residIntBest = residualRegimeTest(internal, ['haloK', 'lhOuter', 'logMeanRun'], 'Internal (haloK+lhOI+MR)');
console.log('');
const residExtHaloK = residualRegimeTest(extOriginal, ['haloK'], 'External orig (haloK only)');
console.log('');
const residPoolHaloK = residualRegimeTest(pooledOriginal, ['haloK'], 'Pooled (haloK only)');

console.log('\n  CRITICAL COMPARISON: Does unexplained VfResid activate faster than haloK?');
function compareActivation(data, label) {
  const valid = data.filter(g => g.haloK !== null && isFinite(g.haloK) && isFinite(g.logA0));
  if (valid.length < 10) return;

  const model = olsFull(valid.map(g => g.VfResid), valid.map(g => [g.haloK]));
  const residuals = valid.map((g, i) => g.VfResid - [1, g.haloK].reduce((s, x, j) => s + x * model.beta[j], 0));
  valid.forEach((g, i) => { g._hkResid = residuals[i]; });

  const regimes = [
    { label: '<100', filter: g => g.Vflat < 100 },
    { label: '100-150', filter: g => g.Vflat >= 100 && g.Vflat < 150 },
    { label: '150-200', filter: g => g.Vflat >= 150 && g.Vflat < 200 },
    { label: '>=200', filter: g => g.Vflat >= 200 },
  ];

  console.log('\n  ' + label + ' (N=' + valid.length + '):');
  console.log('  ' + 'Regime'.padEnd(12) + '  N' + ' r(haloK,a0)'.padStart(13) + ' r(resid,a0)'.padStart(13) + ' resid > haloK?');
  console.log('  ' + '-'.repeat(60));

  for (const reg of regimes) {
    const sub = valid.filter(reg.filter);
    if (sub.length < 3) { console.log('  ' + reg.label.padEnd(12) + String(sub.length).padStart(3) + '  (too few)'); continue; }
    const rHK = pearsonR(sub.map(g => g.logA0), sub.map(g => g.haloK));
    const rResid = pearsonR(sub.map(g => g.logA0), sub.map(g => g._hkResid));
    const faster = Math.abs(rResid) > Math.abs(rHK) ? 'YES' : 'no';
    console.log('  ' + reg.label.padEnd(12) + String(sub.length).padStart(3) +
      (isNaN(rHK) ? '         N/A' : rHK.toFixed(3).padStart(13)) +
      (isNaN(rResid) ? '         N/A' : rResid.toFixed(3).padStart(13)) +
      '  ' + faster);
  }
}

compareActivation(internal, 'Internal');
compareActivation(extOriginal, 'External original');
compareActivation(pooledOriginal, 'Pooled (int+ext_orig)');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROBUSTNESS: ORIGINAL vs AUGMENTED EXTERNAL');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const thresholds = [70, 100, 120, 150, 180];
console.log('  r(VfResid, logA0) at Vflat >= threshold:');
console.log('  ' + 'Thresh'.padStart(6) + ' Ext_orig'.padStart(12) + ' Ext_all'.padStart(12) + '  Consistent?');
console.log('  ' + '-'.repeat(45));
for (const th of thresholds) {
  const sOrig = extOriginal.filter(g => g.Vflat >= th);
  const sAll = extAll.filter(g => g.Vflat >= th);
  const rOrig = sOrig.length >= 4 ? pearsonR(sOrig.map(g => g.logA0), sOrig.map(g => g.VfResid)) : NaN;
  const rAll = sAll.length >= 4 ? pearsonR(sAll.map(g => g.logA0), sAll.map(g => g.VfResid)) : NaN;
  const consistent = (isNaN(rOrig) || isNaN(rAll)) ? '?' : (Math.abs(rOrig - rAll) < 0.1 ? 'YES' : 'DIFFERS');
  console.log('  ' + String(th).padStart(6) +
    (isNaN(rOrig) ? ('  N/A(' + sOrig.length + ')') : (rOrig.toFixed(3) + '(' + sOrig.length + ')').padStart(12)) +
    (isNaN(rAll) ? ('  N/A(' + sAll.length + ')') : (rAll.toFixed(3) + '(' + sAll.length + ')').padStart(12)) +
    '  ' + consistent);
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  SYNTHESIS: PHASE 302 VERDICT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const test1Pass = (shapeInt !== 'FLAT_THROUGHOUT' && shapeInt !== 'NOISY');
const test3IntResid = residIntHaloK && residIntHaloK.regimeResults;
const highVResidR = test3IntResid ? test3IntResid.find(r => r.regime.includes('High'))?.rResidA0 : 0;
const lowVResidR = test3IntResid ? test3IntResid.find(r => r.regime.includes('Low'))?.rResidA0 : 0;
const test3Pass = Math.abs(highVResidR) > Math.abs(lowVResidR) + 0.1;

let verdict = 'INCONCLUSIVE';
if (test1Pass && test3Pass) verdict = 'STRONG_REGIME_ACTIVATION';
else if (test1Pass) verdict = 'REGIME_STRENGTHENING';
else if (test3Pass) verdict = 'RESIDUAL_REGIME_ONLY';

console.log('  Test 1 (Activation shape):');
console.log('    Internal: ' + shapeInt);
console.log('    External original: ' + shapeExtOrig);
console.log('    External all: ' + shapeExtAll);
console.log('    Pooled: ' + shapePool);
console.log('    PASS: ' + test1Pass);

console.log('\n  Test 2 (Driver regime split):');
const haloKLow = drsInt['haloK'] ? drsInt['haloK'][0] : { r: NaN };
const haloKHigh = drsInt['haloK'] ? drsInt['haloK'][1] : { r: NaN };
console.log('    haloK: low-V r=' + (isNaN(haloKLow.r) ? 'N/A' : haloKLow.r.toFixed(3)) +
  ', high-V r=' + (isNaN(haloKHigh.r) ? 'N/A' : haloKHigh.r.toFixed(3)));
console.log('    Pattern: ' + (Math.abs(haloKHigh.r) > Math.abs(haloKLow.r) + 0.1 ? 'haloK STRENGTHENS with regime' : 'haloK relatively stable'));

console.log('\n  Test 3 (Residual regime law):');
console.log('    Unexplained VfResid r(a0): low-V=' + (lowVResidR !== undefined ? lowVResidR.toFixed(3) : 'N/A') +
  ', high-V=' + (highVResidR !== undefined ? highVResidR.toFixed(3) : 'N/A'));
console.log('    PASS: ' + test3Pass);
if (test3Pass) {
  console.log('    *** THE HIDDEN PHYSICS ACTIVATES WITH Vflat ***');
}

console.log('\n  OVERALL VERDICT: ' + verdict);

const output = {
  phase: '302',
  title: 'Regime Law — Does the hidden physics turn on with Vflat?',
  timestamp: new Date().toISOString(),
  verdict,
  test1_activationShape: {
    internal: { shape: shapeInt, cumulative: actInt.cumResults, windows: actInt.windowResults, bins: actInt.binResults },
    externalOriginal: { shape: shapeExtOrig, cumulative: actExtOrig.cumResults },
    externalAll: { shape: shapeExtAll, cumulative: actExtAll.cumResults },
    pooled: { shape: shapePool, cumulative: actPoolOrig.cumResults }
  },
  test2_driverRegimeSplit: {
    internal: Object.fromEntries(Object.entries(drsInt).map(([k, v]) => [k, v.map(x => ({ r: isNaN(x.r) ? null : +x.r.toFixed(3), n: x.n }))])),
    externalOriginal: Object.fromEntries(Object.entries(drsExtOrig).map(([k, v]) => [k, v.map(x => ({ r: isNaN(x.r) ? null : +x.r.toFixed(3), n: x.n }))])),
    pooled: Object.fromEntries(Object.entries(drsPool).map(([k, v]) => [k, v.map(x => ({ r: isNaN(x.r) ? null : +x.r.toFixed(3), n: x.n }))]))
  },
  test3_residualRegimeLaw: {
    internalHaloK: residIntHaloK,
    internalBest: residIntBest,
    externalHaloK: residExtHaloK,
    pooledHaloK: residPoolHaloK
  },
  robustness: { note: 'Original vs augmented external patterns compared above' },
  caveats: [
    'Sliding window N varies — small bins have high variance',
    'Internal N=45 makes regime splits very small (some bins N<10)',
    'Salvaged galaxies have lower data quality — augmented results should be compared against original',
    'Sharp vs smooth distinction depends on bin choice — intrinsic uncertainty exists',
    'Residual regime test assumes haloK model is the right baseline for subtraction'
  ]
};

const outPath = pub('phase302-regime-law.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nOutput: ' + outPath);
console.log('\n======================================================================');
console.log('  PHASE 302 COMPLETE');
console.log('  Next: Phase 303 — Physical interpretation');
console.log('======================================================================');

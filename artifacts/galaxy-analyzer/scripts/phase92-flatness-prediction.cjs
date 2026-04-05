#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 92: DOES THE HIGH-REGIME LAW PREDICT RC FLATNESS?');
console.log('  From a₀ law → rotation curve shape');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0(pts) {
  let lo = 0.5, hi = 5.0;
  for (let s = 0; s < 200; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gb = Math.pow(10, p.logGbar);
      c1 += (p.logGobs - Math.log10(mcgaughRAR(gb, Math.pow(10, m1)))) ** 2;
      c2 += (p.logGobs - Math.log10(mcgaughRAR(gb, Math.pow(10, m2)))) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  return { logA0: (lo + hi) / 2, a0: Math.pow(10, (lo + hi) / 2) };
}

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
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
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, rms: Math.sqrt(ss_res / n), residuals, yhat };
}

function looCV(X, y) {
  const n = y.length; let ss = 0;
  for (let i = 0; i < n; i++) {
    const Xt = X.filter((_, j) => j !== i);
    const yt = y.filter((_, j) => j !== i);
    const fit = olsRegress(Xt, yt);
    if (!fit) return NaN;
    const pred = X[i].reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss += (y[i] - pred) ** 2;
  }
  return Math.sqrt(ss / n);
}

function permTest(X, y, nPerm) {
  const realFit = olsRegress(X, y);
  if (!realFit) return { pValue: 1, realR2: 0 };
  let nBetter = 0;
  for (let p = 0; p < nPerm; p++) {
    const shuffled = [...y].sort(() => Math.random() - 0.5);
    const fit = olsRegress(X, shuffled);
    if (fit && fit.r2 >= realFit.r2) nBetter++;
  }
  return { pValue: nBetter / nPerm, realR2: realFit.r2 };
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { T, L36, Rdisk, MHI, Vflat, Q };
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
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul });
}

const p25 = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase25-group-membership.json'), 'utf-8'));
const envLookup = {};
for (const g of p25.galaxyAssignments) envLookup[g.name] = g;

const allGalaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    const gObs = pt.vobs * pt.vobs / pt.rad;
    const gBar = Math.abs(vBarSq) / pt.rad;
    if (gBar <= 0 || gObs <= 0 || !isFinite(gBar) || !isFinite(gObs)) continue;
    pts.push({ r: pt.rad, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), vobs: pt.vobs });
  }
  if (pts.length < 5) continue;

  const fit = fitA0(pts);
  if (!isFinite(fit.logA0) || fit.logA0 < 0.5 || fit.logA0 > 5.5) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);

  const half = Math.floor(sorted.length / 2);
  const outer = sorted.slice(half);
  const logR = outer.map(p => Math.log10(p.r));
  const logV = outer.map(p => Math.log10(p.vobs));
  const mrr = mean(logR), mvv = mean(logV);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mrr) * (logV[i] - mvv); sxx += (logR[i] - mrr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const vMax = Math.max(...sorted.map(p => p.vobs));
  const vLast3 = mean(sorted.slice(-3).map(p => p.vobs));
  const flatness = vMax > 0 ? vLast3 / vMax : 1;

  const vOuter = outer.map(p => p.vobs);
  const outerCV = sd(vOuter) / mean(vOuter);

  const rarResid = sorted.map(p => {
    const gb = Math.pow(10, p.logGbar);
    const pred = mcgaughRAR(gb, fit.a0);
    return p.logGobs - Math.log10(pred > 0 ? pred : 1e-20);
  });
  let currentSign = Math.sign(rarResid[0]);
  let runLen = 1;
  const residRun = [];
  for (let i = 1; i < rarResid.length; i++) {
    const s = Math.sign(rarResid[i]);
    if (s === currentSign) runLen++;
    else { residRun.push(runLen); runLen = 1; currentSign = s; }
  }
  residRun.push(runLen);
  const meanRunLen = residRun.reduce((s, v) => s + v, 0) / residRun.length;

  const ev = envLookup[name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;
  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI: Math.log10(t1.MHI),
    logMeanRun: Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1),
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    Vflat: t1.Vflat, Q: t1.Q, T: t1.T,
    envCode,
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    logL36: Math.log10(t1.L36),
    outerSlope, flatness, outerCV, fgas,
  });
}

const hi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);

console.log('  High-regime galaxies: N=' + hi.length);
console.log('  Flatness targets:');
console.log('    outerSlope: mean=' + mean(hi.map(g => g.outerSlope)).toFixed(4) + ', sd=' + sd(hi.map(g => g.outerSlope)).toFixed(4));
console.log('    flatness (vLast3/vMax): mean=' + mean(hi.map(g => g.flatness)).toFixed(4) + ', sd=' + sd(hi.map(g => g.flatness)).toFixed(4));
console.log('    outerCV: mean=' + mean(hi.map(g => g.outerCV)).toFixed(4) + ', sd=' + sd(hi.map(g => g.outerCV)).toFixed(4));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 1: CORRELATIONS — logA0 vs flatness metrics');
console.log('══════════════════════════════════════════════════════════════\n');

const targets = [
  { name: 'outerSlope', fn: g => g.outerSlope },
  { name: 'flatness', fn: g => g.flatness },
  { name: 'outerCV', fn: g => g.outerCV },
];

const predictors = [
  { name: 'logA0', fn: g => g.logA0 },
  { name: 'logMHI', fn: g => g.logMHI },
  { name: 'logMeanRun', fn: g => g.logMeanRun },
  { name: 'logSigma0', fn: g => g.logSigma0 },
  { name: 'logVflat', fn: g => g.logVflat },
  { name: 'envCode', fn: g => g.envCode },
  { name: 'fgas', fn: g => g.fgas },
  { name: 'T', fn: g => g.T },
];

console.log('  ' + 'Predictor'.padEnd(15));
for (const t of targets) process.stdout.write(t.name.padStart(12));
console.log();

const corrMatrix = {};
for (const p of predictors) {
  corrMatrix[p.name] = {};
  process.stdout.write('  ' + p.name.padEnd(15));
  for (const t of targets) {
    const r = pearsonR(hi.map(p.fn), hi.map(t.fn));
    corrMatrix[p.name][t.name] = +r.toFixed(3);
    process.stdout.write((isNaN(r) ? '---' : r.toFixed(3)).padStart(12));
  }
  console.log();
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: PREDICT outerSlope FROM GALAXY PROPERTIES');
console.log('  (same predictors as the a₀ law)');
console.log('══════════════════════════════════════════════════════════════\n');

function testTarget(data, targetFn, targetName) {
  const y = data.map(targetFn);
  const m0sd = sd(y);

  const models = [
    { name: 'logA0', make: () => data.map(g => [1, g.logA0]) },
    { name: 'logMHI', make: () => data.map(g => [1, g.logMHI]) },
    { name: 'logMeanRun', make: () => data.map(g => [1, g.logMeanRun]) },
    { name: 'logSigma0', make: () => data.map(g => [1, g.logSigma0]) },
    { name: 'logVflat', make: () => data.map(g => [1, g.logVflat]) },
    { name: 'fgas', make: () => data.map(g => [1, g.fgas]) },
    { name: 'M3 vars', make: () => data.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0]) },
    { name: 'M5 vars', make: () => data.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0, g.envCode, g.logVflat]) },
    { name: 'logA0 + logVflat', make: () => data.map(g => [1, g.logA0, g.logVflat]) },
    { name: 'logA0 + fgas', make: () => data.map(g => [1, g.logA0, g.fgas]) },
    { name: 'M5 + logA0', make: () => data.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0, g.envCode, g.logVflat, g.logA0]) },
    { name: 'fgas + logVflat', make: () => data.map(g => [1, g.fgas, g.logVflat]) },
    { name: 'fgas + logMHI', make: () => data.map(g => [1, g.fgas, g.logMHI]) },
    { name: 'fgas+Vflat+MHI', make: () => data.map(g => [1, g.fgas, g.logVflat, g.logMHI]) },
  ];

  console.log('  Target: ' + targetName + ' (N=' + data.length + ', M0 sd=' + m0sd.toFixed(4) + ')\n');
  console.log('  ' + 'Model'.padEnd(22) + 'R2'.padStart(8) + 'RMS'.padStart(9) +
    'LOO'.padStart(9) + 'Gap%'.padStart(8));

  const results = [];
  for (const m of models) {
    const X = m.make();
    const fit = olsRegress(X, y);
    if (!fit) { console.log('  ' + m.name.padEnd(22) + ' FAILED'); continue; }
    const loo = looCV(X, y);
    const gap = m0sd > 0 ? (m0sd - loo) / m0sd * 100 : 0;
    results.push({ name: m.name, r2: +fit.r2.toFixed(4), rms: +fit.rms.toFixed(4), loo: +loo.toFixed(4), gap: +gap.toFixed(1), beta: fit.beta.map(b => +b.toFixed(4)) });
    console.log('  ' + m.name.padEnd(22) + fit.r2.toFixed(4).padStart(8) +
      fit.rms.toFixed(4).padStart(9) + loo.toFixed(4).padStart(9) + gap.toFixed(1).padStart(8) + '%');
  }

  const best = results.filter(r => r.gap > 0).sort((a, b) => b.gap - a.gap);
  if (best.length > 0) {
    console.log('\n  Best: ' + best[0].name + ' (gap=' + best[0].gap + '%, R²=' + best[0].r2 + ')');
  } else {
    console.log('\n  No model beats M0 on LOO.');
  }
  return { results, m0sd: +m0sd.toFixed(4), best: best.length > 0 ? best[0] : null };
}

const slopeResults = testTarget(hi, g => g.outerSlope, 'outerSlope');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: PREDICT flatness (vLast3/vMax)');
console.log('══════════════════════════════════════════════════════════════\n');

const flatResults = testTarget(hi, g => g.flatness, 'flatness');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: PREDICT outerCV (variation coefficient)');
console.log('══════════════════════════════════════════════════════════════\n');

const cvResults = testTarget(hi, g => g.outerCV, 'outerCV');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: PERMUTATION TEST FOR BEST MODELS');
console.log('══════════════════════════════════════════════════════════════\n');

const permTargets = [
  { name: 'outerSlope', results: slopeResults, fn: g => g.outerSlope },
  { name: 'flatness', results: flatResults, fn: g => g.flatness },
  { name: 'outerCV', results: cvResults, fn: g => g.outerCV },
];

const permResults = {};
for (const pt of permTargets) {
  if (!pt.results.best) { permResults[pt.name] = null; continue; }
  const modelDef = [
    { name: 'logA0', make: () => hi.map(g => [1, g.logA0]) },
    { name: 'logMHI', make: () => hi.map(g => [1, g.logMHI]) },
    { name: 'logMeanRun', make: () => hi.map(g => [1, g.logMeanRun]) },
    { name: 'logSigma0', make: () => hi.map(g => [1, g.logSigma0]) },
    { name: 'logVflat', make: () => hi.map(g => [1, g.logVflat]) },
    { name: 'fgas', make: () => hi.map(g => [1, g.fgas]) },
    { name: 'M3 vars', make: () => hi.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0]) },
    { name: 'M5 vars', make: () => hi.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0, g.envCode, g.logVflat]) },
    { name: 'logA0 + logVflat', make: () => hi.map(g => [1, g.logA0, g.logVflat]) },
    { name: 'logA0 + fgas', make: () => hi.map(g => [1, g.logA0, g.fgas]) },
    { name: 'M5 + logA0', make: () => hi.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0, g.envCode, g.logVflat, g.logA0]) },
    { name: 'fgas + logVflat', make: () => hi.map(g => [1, g.fgas, g.logVflat]) },
    { name: 'fgas + logMHI', make: () => hi.map(g => [1, g.fgas, g.logMHI]) },
    { name: 'fgas+Vflat+MHI', make: () => hi.map(g => [1, g.fgas, g.logVflat, g.logMHI]) },
  ].find(m => m.name === pt.results.best.name);
  if (!modelDef) continue;
  const X = modelDef.make();
  const y = hi.map(pt.fn);
  const perm = permTest(X, y, 2000);
  permResults[pt.name] = { model: pt.results.best.name, r2: pt.results.best.r2, gap: pt.results.best.gap, pValue: perm.pValue };
  console.log('  ' + pt.name + ': model=' + pt.results.best.name + ', R²=' + pt.results.best.r2 +
    ', gap=' + pt.results.best.gap + '%, perm p=' + perm.pValue.toFixed(3));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 6: DOES logA0 DIRECTLY PREDICT FLATNESS?');
console.log('  (Is the a₀ law itself connected to RC shape?)');
console.log('══════════════════════════════════════════════════════════════\n');

const rA0slope = pearsonR(hi.map(g => g.logA0), hi.map(g => g.outerSlope));
const rA0flat = pearsonR(hi.map(g => g.logA0), hi.map(g => g.flatness));
const rA0cv = pearsonR(hi.map(g => g.logA0), hi.map(g => g.outerCV));

console.log('  r(logA0, outerSlope) = ' + rA0slope.toFixed(3));
console.log('  r(logA0, flatness)   = ' + rA0flat.toFixed(3));
console.log('  r(logA0, outerCV)    = ' + rA0cv.toFixed(3));

console.log('\n  Interpretation:');
if (Math.abs(rA0slope) > 0.3) {
  console.log('  ★ logA0 has a DIRECT connection to RC shape');
  console.log('    Higher/lower a₀ → different outer slope');
} else if (Math.abs(rA0slope) > 0.15) {
  console.log('  ★ logA0 has a WEAK connection to RC shape');
} else {
  console.log('  ★ logA0 has NO direct connection to RC shape');
  console.log('    a₀ variation and flatness are independent phenomena');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 7: WHAT PREDICTS WHETHER RC IS FLAT VS RISING?');
console.log('  (Binary classification: |slope| < 0.1 = flat)');
console.log('══════════════════════════════════════════════════════════════\n');

const hiFlat = hi.filter(g => Math.abs(g.outerSlope) <= 0.1);
const hiRising = hi.filter(g => g.outerSlope > 0.1);
const hiDecl = hi.filter(g => g.outerSlope < -0.1);

console.log('  RC shape in high regime:');
console.log('    Flat: N=' + hiFlat.length + ' (' + (hiFlat.length / hi.length * 100).toFixed(0) + '%)');
console.log('    Rising: N=' + hiRising.length + ' (' + (hiRising.length / hi.length * 100).toFixed(0) + '%)');
console.log('    Declining: N=' + hiDecl.length + ' (' + (hiDecl.length / hi.length * 100).toFixed(0) + '%)');

console.log('\n  Properties of flat vs rising (high regime):');
const compareVars = [
  { name: 'logA0', fn: g => g.logA0 },
  { name: 'logMHI', fn: g => g.logMHI },
  { name: 'fgas', fn: g => g.fgas },
  { name: 'logVflat', fn: g => g.logVflat },
  { name: 'logSigma0', fn: g => g.logSigma0 },
  { name: 'T', fn: g => g.T },
  { name: 'logMeanRun', fn: g => g.logMeanRun },
  { name: 'nPts', fn: g => g.nPts },
];

console.log('  ' + 'Variable'.padEnd(15) + 'Flat'.padStart(9) + 'Rising'.padStart(9) + 'Diff'.padStart(9));
for (const v of compareVars) {
  const fv = mean(hiFlat.map(v.fn)), rv = mean(hiRising.map(v.fn));
  console.log('  ' + v.name.padEnd(15) + fv.toFixed(3).padStart(9) + rv.toFixed(3).padStart(9) +
    (fv - rv).toFixed(3).padStart(9));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 8: RESIDUAL CONNECTION');
console.log('  After removing Vflat effect, does a₀ residual predict');
console.log('  flatness residual?');
console.log('══════════════════════════════════════════════════════════════\n');

const vflatX = hi.map(g => [1, g.logVflat]);
const a0y = hi.map(g => g.logA0);
const slopey = hi.map(g => g.outerSlope);

const a0fit = olsRegress(vflatX, a0y);
const slopefit = olsRegress(vflatX, slopey);

if (a0fit && slopefit) {
  const a0resid = a0fit.residuals;
  const sloperesid = slopefit.residuals;
  const rResid = pearsonR(a0resid, sloperesid);
  console.log('  r(logA0 residual, outerSlope residual | logVflat) = ' + rResid.toFixed(3));
  console.log('  (After removing the common Vflat dependence)');

  if (Math.abs(rResid) > 0.15) {
    console.log('  → a₀ and flatness share structure BEYOND their common Vflat dependence');
  } else {
    console.log('  → No residual connection: a₀ variation and flatness are mediated by');
    console.log('    different aspects of galaxy structure, connected only through Vflat');
  }
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  SYNTHESIS');
console.log('══════════════════════════════════════════════════════════════\n');

const a0PredSlope = slopeResults.results.find(r => r.name === 'logA0');
const bestSlopeModel = slopeResults.best;

console.log('  KEY FINDINGS:');
console.log();
console.log('  1. logA0 → outerSlope: r = ' + rA0slope.toFixed(3));
console.log('     ' + (Math.abs(rA0slope) > 0.3 ? 'STRONG' : Math.abs(rA0slope) > 0.15 ? 'WEAK' : 'NO') + ' direct connection');
console.log();
console.log('  2. Best flatness predictor: ' + (bestSlopeModel ? bestSlopeModel.name + ' (gap=' + bestSlopeModel.gap + '%)' : 'NONE'));
console.log();
console.log('  3. Can M3/M5 predict RC shape?');
const m3slope = slopeResults.results.find(r => r.name === 'M3 vars');
const m5slope = slopeResults.results.find(r => r.name === 'M5 vars');
if (m3slope) console.log('     M3 → outerSlope: R²=' + m3slope.r2 + ', LOO gap=' + m3slope.gap + '%');
if (m5slope) console.log('     M5 → outerSlope: R²=' + m5slope.r2 + ', LOO gap=' + m5slope.gap + '%');

const output = {
  phase: '92',
  title: 'Does the High-Regime Law Predict RC Flatness?',
  hiN: hi.length,
  correlations: corrMatrix,
  outerSlope: slopeResults,
  flatness: flatResults,
  outerCV: cvResults,
  permutation: permResults,
  directA0connection: { outerSlope: +rA0slope.toFixed(3), flatness: +rA0flat.toFixed(3), outerCV: +rA0cv.toFixed(3) },
  rcShapeDistribution: { flat: hiFlat.length, rising: hiRising.length, declining: hiDecl.length },
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase92-flatness-prediction.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase92-flatness-prediction.json');

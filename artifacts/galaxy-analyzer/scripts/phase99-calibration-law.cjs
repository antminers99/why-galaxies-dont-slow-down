#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 99: CALIBRATION LAW');
console.log('  Why does the slope differ between N45 and non-N45?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function median(a) { const s = [...a].sort((x, y) => x - y); return s.length % 2 ? s[Math.floor(s.length / 2)] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2; }
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
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, residuals, yhat };
}

function looR2(X_fn, y_fn, data) {
  const n = data.length;
  let ss_res = 0, ss_tot = 0;
  const yAll = data.map(y_fn);
  const yMean = mean(yAll);
  for (let i = 0; i < n; i++) {
    const train = data.filter((_, j) => j !== i);
    const fit = olsRegress(train.map(X_fn), train.map(y_fn));
    if (!fit) return NaN;
    const yPred = X_fn(data[i]).reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss_res += (yAll[i] - yPred) ** 2;
    ss_tot += (yAll[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
}

const n45csv = fs.readFileSync(path.join(__dirname, '..', 'public', 'replication', 'N45_final_dataset.csv'), 'utf-8').trim().split('\n');
const n45names = new Set();
const n45data = {};
for (let i = 1; i < n45csv.length; i++) {
  const cols = n45csv[i].split(',');
  n45names.add(cols[0]);
  n45data[cols[0]] = { D_Mpc: parseFloat(cols[10]) };
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const D = parseFloat(line.substring(15, 24).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { T, D, L36, Rdisk, MHI, Vflat, Q };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs });
}

const allGalaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0) continue;

  const pts = rcPoints.filter(p => p.rad > 0 && p.vobs > 0).sort((a, b) => a.rad - b.rad);
  if (pts.length < 8) continue;

  const half = Math.floor(pts.length / 2);
  const innerPts = pts.slice(0, half);
  const outerPts = pts.slice(half);

  const innerVmax = Math.max(...innerPts.map(p => p.vobs));

  const logR = outerPts.map(p => Math.log10(p.r || p.rad));
  const logV = outerPts.map(p => Math.log10(p.vobs));
  const mr = mean(logR), mv = mean(logV);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mr) * (logV[i] - mv); sxx += (logR[i] - mr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const rMax = pts[pts.length - 1].rad;
  const rMin = pts[0].rad;
  const radialCoverage = rMax / (t1.Rdisk > 0 ? t1.Rdisk : 1);
  const logRadialRange = Math.log10(rMax / (rMin > 0 ? rMin : 0.1));
  const innerFrac = half / pts.length;
  const innerRadialSpan = innerPts[innerPts.length - 1].rad - innerPts[0].rad;
  const outerRadialSpan = outerPts[outerPts.length - 1].rad - outerPts[0].rad;
  const angularRes = t1.D > 0 ? t1.Rdisk / t1.D * 206265 : NaN;

  const logRatio = Math.log10(t1.Vflat > 0 && innerVmax > 0 ? t1.Vflat / innerVmax : 1);

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);

  allGalaxies.push({
    name, nPts: pts.length, nInner: half, nOuter: pts.length - half,
    innerVmax, logInnerVmax: Math.log10(innerVmax > 0 ? innerVmax : 1),
    outerSlope, logRatio,
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    Vflat: t1.Vflat, fgas,
    logMHI: Math.log10(t1.MHI),
    T: t1.T, Q: t1.Q, D: t1.D,
    radialCoverage, logRadialRange,
    innerRadialSpan, outerRadialSpan,
    angularRes: isFinite(angularRes) ? angularRes : NaN,
    Rdisk: t1.Rdisk,
    isN45: n45names.has(name),
    sdOuterSlope: sd(outerPts.map(p => Math.log10(p.vobs))),
  });
}

const hi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);
const n45hi = hi.filter(g => g.isN45);
const nonN45hi = hi.filter(g => !g.isN45);

console.log('  High-regime: N=' + hi.length + ' (N45=' + n45hi.length + ', non-N45=' + nonN45hi.length + ')\n');

console.log('══════════════════════════════════════════════════════════════');
console.log('  PART A: SAMPLE DIFFERENCES — What differs between N45');
console.log('  and non-N45 that could affect the calibration?');
console.log('══════════════════════════════════════════════════════════════\n');

const props = [
  { name: 'nPts', fn: g => g.nPts },
  { name: 'Distance (Mpc)', fn: g => g.D },
  { name: 'logVflat', fn: g => g.logVflat },
  { name: 'logInnerVmax', fn: g => g.logInnerVmax },
  { name: 'logRatio', fn: g => g.logRatio },
  { name: 'outerSlope', fn: g => g.outerSlope },
  { name: 'Q', fn: g => g.Q },
  { name: 'radialCoverage (R/Rd)', fn: g => g.radialCoverage },
  { name: 'logRadialRange', fn: g => g.logRadialRange },
  { name: 'logSigma0', fn: g => g.logSigma0 },
  { name: 'fgas', fn: g => g.fgas },
  { name: 'nInner', fn: g => g.nInner },
  { name: 'angularRes (arcsec)', fn: g => g.angularRes },
  { name: 'Rdisk (kpc)', fn: g => g.Rdisk },
  { name: 'sd(outerSlope)', fn: g => g.sdOuterSlope },
];

console.log('  ' + 'Property'.padEnd(28) + 'N45 median'.padStart(12) + 'non-N45 med'.padStart(13) + '  Δ%'.padStart(8));
console.log('  ' + '-'.repeat(62));

const sampleDiffs = [];
for (const p of props) {
  const v1 = n45hi.map(p.fn).filter(v => isFinite(v));
  const v2 = nonN45hi.map(p.fn).filter(v => isFinite(v));
  if (v1.length < 3 || v2.length < 3) continue;
  const m1 = median(v1), m2 = median(v2);
  const avg = (Math.abs(m1) + Math.abs(m2)) / 2;
  const pctDiff = avg > 0 ? ((m1 - m2) / avg * 100).toFixed(1) : '—';
  console.log('  ' + p.name.padEnd(28) + m1.toFixed(3).padStart(12) + m2.toFixed(3).padStart(13) + pctDiff.toString().padStart(8));
  sampleDiffs.push({ name: p.name, n45med: +m1.toFixed(3), nonMed: +m2.toFixed(3), pctDiff: +pctDiff });
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART B: VARIANCE RANGE RESTRICTION');
console.log('  Does N45 have less variance in the key variables?');
console.log('══════════════════════════════════════════════════════════════\n');

const sdRatio45 = sd(n45hi.map(g => g.logRatio));
const sdRatioNon = sd(nonN45hi.map(g => g.logRatio));
const sdSlope45 = sd(n45hi.map(g => g.outerSlope));
const sdSlopeNon = sd(nonN45hi.map(g => g.outerSlope));
const rangeRatio45 = Math.max(...n45hi.map(g => g.logRatio)) - Math.min(...n45hi.map(g => g.logRatio));
const rangeRatioNon = Math.max(...nonN45hi.map(g => g.logRatio)) - Math.min(...nonN45hi.map(g => g.logRatio));

console.log('  sd(logRatio):    N45 = ' + sdRatio45.toFixed(3) + ', non-N45 = ' + sdRatioNon.toFixed(3) + ' (ratio ' + (sdRatio45 / sdRatioNon).toFixed(2) + ')');
console.log('  sd(outerSlope):  N45 = ' + sdSlope45.toFixed(3) + ', non-N45 = ' + sdSlopeNon.toFixed(3) + ' (ratio ' + (sdSlope45 / sdSlopeNon).toFixed(2) + ')');
console.log('  range(logRatio): N45 = ' + rangeRatio45.toFixed(3) + ', non-N45 = ' + rangeRatioNon.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART C: WHICH PROPERTY MODULATES THE SLOPE?');
console.log('  Interaction test: outerSlope ~ logRatio * moderator');
console.log('══════════════════════════════════════════════════════════════\n');

const moderators = [
  { name: 'nPts', fn: g => g.nPts },
  { name: 'logD', fn: g => Math.log10(g.D > 0 ? g.D : 1) },
  { name: 'Q', fn: g => g.Q },
  { name: 'radialCoverage', fn: g => g.radialCoverage },
  { name: 'logRadialRange', fn: g => g.logRadialRange },
  { name: 'nInner', fn: g => g.nInner },
  { name: 'logRdisk', fn: g => Math.log10(g.Rdisk > 0 ? g.Rdisk : 0.1) },
  { name: 'fgas', fn: g => g.fgas },
  { name: 'logSigma0', fn: g => g.logSigma0 },
];

const interResults = [];
for (const mod of moderators) {
  const valid = hi.filter(g => isFinite(mod.fn(g)));
  const X_main = valid.map(g => [1, g.logRatio]);
  const X_inter = valid.map(g => [1, g.logRatio, mod.fn(g), g.logRatio * mod.fn(g)]);
  const y = valid.map(g => g.outerSlope);

  const fitMain = olsRegress(X_main, y);
  const fitInter = olsRegress(X_inter, y);
  if (!fitMain || !fitInter) continue;

  const deltaR2 = fitInter.r2 - fitMain.r2;
  const interCoeff = fitInter.beta[3];

  console.log('  ' + mod.name + ':');
  console.log('    Interaction β = ' + interCoeff.toFixed(4) + ', ΔR² = ' + deltaR2.toFixed(3));

  interResults.push({ name: mod.name, interBeta: +interCoeff.toFixed(4), deltaR2: +deltaR2.toFixed(3) });
}

interResults.sort((a, b) => Math.abs(b.deltaR2) - Math.abs(a.deltaR2));
console.log('\n  Ranked by |ΔR²|:');
for (const ir of interResults) {
  console.log('    ' + ir.name + ': ΔR² = ' + ir.deltaR2 + ', β_inter = ' + ir.interBeta);
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART D: CORRECTED MODEL — include top moderator');
console.log('══════════════════════════════════════════════════════════════\n');

const topMod = moderators.find(m => m.name === interResults[0].name);
if (topMod) {
  const valid = hi.filter(g => isFinite(topMod.fn(g)));
  const X_corrected = g => [1, g.logRatio, topMod.fn(g), g.logRatio * topMod.fn(g)];
  const fitCorr = olsRegress(valid.map(X_corrected), valid.map(g => g.outerSlope));
  const looCorr = looR2(X_corrected, g => g.outerSlope, valid);

  console.log('  Corrected model: outerSlope ~ logRatio * ' + topMod.name);
  console.log('  R² = ' + fitCorr.r2.toFixed(3) + ', LOO R² = ' + looCorr.toFixed(3));
  console.log('  Coefficients: intercept=' + fitCorr.beta[0].toFixed(3) + ', logRatio=' + fitCorr.beta[1].toFixed(3) +
    ', ' + topMod.name + '=' + fitCorr.beta[2].toFixed(4) + ', interaction=' + fitCorr.beta[3].toFixed(4));

  const n45valid = valid.filter(g => g.isN45);
  const nonValid = valid.filter(g => !g.isN45);

  function trainTest(trainSet, testSet, label) {
    const fit = olsRegress(trainSet.map(X_corrected), trainSet.map(g => g.outerSlope));
    if (!fit) return null;
    const yTest = testSet.map(g => g.outerSlope);
    const yPred = testSet.map(g => X_corrected(g).reduce((s, v, j) => s + v * fit.beta[j], 0));
    const ss_res = yTest.reduce((s, v, i) => s + (v - yPred[i]) ** 2, 0);
    const ss_tot = yTest.reduce((s, v) => s + (v - mean(yTest)) ** 2, 0);
    const testR2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
    console.log('  ' + label + ': test R² = ' + testR2.toFixed(3) + ', r(pred,actual) = ' + pearsonR(yPred, yTest).toFixed(3));
    return testR2;
  }

  console.log('\n  Cross-validation with corrected model:');
  const corrNtoN45 = trainTest(nonValid, n45valid, 'Train non-N45 → test N45');
  const corrN45toN = trainTest(n45valid, nonValid, 'Train N45 → test non-N45');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART E: ALTERNATIVE — NORMALIZED RATIO');
console.log('  Use (InnerVmax/Vflat) directly as fraction [0,1+]');
console.log('══════════════════════════════════════════════════════════════\n');

const X_frac = g => [1, g.innerVmax / g.Vflat];
const fitFrac = olsRegress(hi.map(X_frac), hi.map(g => g.outerSlope));
const looFrac = looR2(X_frac, g => g.outerSlope, hi);
console.log('  outerSlope ~ InnerVmax/Vflat (linear fraction):');
console.log('  R² = ' + fitFrac.r2.toFixed(3) + ', LOO R² = ' + looFrac.toFixed(3));
console.log('  outerSlope = ' + fitFrac.beta[0].toFixed(3) + ' + ' + fitFrac.beta[1].toFixed(3) + ' * (InnerVmax/Vflat)');

const fitFracN45 = olsRegress(n45hi.map(X_frac), n45hi.map(g => g.outerSlope));
const fitFracNon = olsRegress(nonN45hi.map(X_frac), nonN45hi.map(g => g.outerSlope));
console.log('\n  Coefficient comparison (linear fraction):');
console.log('    N45:     a=' + fitFracN45.beta[0].toFixed(3) + ', b=' + fitFracN45.beta[1].toFixed(3));
console.log('    non-N45: a=' + fitFracNon.beta[0].toFixed(3) + ', b=' + fitFracNon.beta[1].toFixed(3));
const fracDiff = Math.abs(fitFracN45.beta[1] - fitFracNon.beta[1]) / ((Math.abs(fitFracN45.beta[1]) + Math.abs(fitFracNon.beta[1])) / 2) * 100;
console.log('    Slope difference: ' + fracDiff.toFixed(1) + '%');

function trainTestSimple(X_fn, trainSet, testSet) {
  const fit = olsRegress(trainSet.map(X_fn), trainSet.map(g => g.outerSlope));
  if (!fit) return NaN;
  const yTest = testSet.map(g => g.outerSlope);
  const yPred = testSet.map(g => X_fn(g).reduce((s, v, j) => s + v * fit.beta[j], 0));
  const ss_res = yTest.reduce((s, v, i) => s + (v - yPred[i]) ** 2, 0);
  const ss_tot = yTest.reduce((s, v) => s + (v - mean(yTest)) ** 2, 0);
  return ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
}

const fracNtoN45 = trainTestSimple(X_frac, nonN45hi, n45hi);
const fracN45toN = trainTestSimple(X_frac, n45hi, nonN45hi);
console.log('\n  Cross-validation (linear fraction):');
console.log('    Train non-N45 → test N45: R² = ' + fracNtoN45.toFixed(3));
console.log('    Train N45 → test non-N45: R² = ' + fracN45toN.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART F: ROBUST REGRESSION — median-based slope');
console.log('══════════════════════════════════════════════════════════════\n');

function theilSenSlope(x, y) {
  const n = x.length;
  const slopes = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (x[j] !== x[i]) slopes.push((y[j] - y[i]) / (x[j] - x[i]));
    }
  }
  const b = median(slopes);
  const a = median(y.map((yi, i) => yi - b * x[i]));
  return { a, b };
}

const tsAll = theilSenSlope(hi.map(g => g.logRatio), hi.map(g => g.outerSlope));
const tsN45 = theilSenSlope(n45hi.map(g => g.logRatio), n45hi.map(g => g.outerSlope));
const tsNon = theilSenSlope(nonN45hi.map(g => g.logRatio), nonN45hi.map(g => g.outerSlope));

console.log('  Theil-Sen slopes:');
console.log('    All:     b = ' + tsAll.b.toFixed(3) + ', a = ' + tsAll.a.toFixed(3));
console.log('    N45:     b = ' + tsN45.b.toFixed(3) + ', a = ' + tsN45.a.toFixed(3));
console.log('    non-N45: b = ' + tsNon.b.toFixed(3) + ', a = ' + tsNon.a.toFixed(3));
const tsDiff = Math.abs(tsN45.b - tsNon.b) / ((Math.abs(tsN45.b) + Math.abs(tsNon.b)) / 2) * 100;
console.log('    Theil-Sen slope difference: ' + tsDiff.toFixed(1) + '%');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const output = {
  phase: '99',
  title: 'Calibration Law',
  n45hi: n45hi.length, nonN45hi: nonN45hi.length,
  sampleDiffs,
  varianceRestriction: {
    sdRatio: { n45: +sdRatio45.toFixed(3), non: +sdRatioNon.toFixed(3), ratio: +(sdRatio45 / sdRatioNon).toFixed(2) },
    sdSlope: { n45: +sdSlope45.toFixed(3), non: +sdSlopeNon.toFixed(3), ratio: +(sdSlope45 / sdSlopeNon).toFixed(2) },
  },
  interactionRank: interResults.slice(0, 5),
  theilSen: {
    all: +tsAll.b.toFixed(3), n45: +tsN45.b.toFixed(3), non: +tsNon.b.toFixed(3),
    pctDiff: +tsDiff.toFixed(1),
  },
  linearFraction: {
    slopeDiff: +fracDiff.toFixed(1),
    crossVal: { nonToN45: +fracNtoN45.toFixed(3), n45ToNon: +fracN45toN.toFixed(3) },
  },
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase99-calibration-law.json'), JSON.stringify(output, null, 2));
console.log('  Saved: public/phase99-calibration-law.json');

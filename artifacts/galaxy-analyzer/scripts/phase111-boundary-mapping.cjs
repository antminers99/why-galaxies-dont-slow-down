#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 111: REGIME BOUNDARY MAPPING');
console.log('  WHERE exactly does fgas stop working?');
console.log('  Turn "regime-specific" from vague caveat into precise result.');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }

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
  const ss_res = y.reduce((s, v, i) => s + (v - yhat[i]) ** 2, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, ss_res, yhat };
}

function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function looR2(xArr, yArr) {
  const n = xArr.length;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(yArr);
  for (let i = 0; i < n; i++) {
    const xTr = xArr.filter((_, j) => j !== i);
    const yTr = yArr.filter((_, j) => j !== i);
    const fit = olsRegress(xTr.map(v => [1, v]), yTr);
    if (!fit) return NaN;
    const pred = fit.beta[0] + fit.beta[1] * xArr[i];
    ss_loo += (yArr[i] - pred) ** 2;
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
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul });
}

const allGalaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat <= 0) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    pts.push({ r: pt.rad, vobs: pt.vobs, vBar: Math.sqrt(Math.max(vBarSq, 0.01)) });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);
  const outerMD = outerPts.map(p => (p.vobs ** 2) / (p.vBar ** 2)).filter(v => isFinite(v) && v > 0);
  const logOMD = Math.log10(mean(outerMD));
  if (!isFinite(logOMD)) continue;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  if (!isFinite(fgas)) continue;

  allGalaxies.push({ name, logOMD, fgas, Vflat: t1.Vflat });
}

allGalaxies.sort((a, b) => a.Vflat - b.Vflat);
console.log('  Total galaxies with valid data: N=' + allGalaxies.length);
console.log('  Vflat range: [' + allGalaxies[0].Vflat.toFixed(1) + ', ' + allGalaxies[allGalaxies.length - 1].Vflat.toFixed(1) + '] km/s\n');

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: SLIDING LOWER THRESHOLD');
console.log('  Keep all above threshold, measure r and LOO R²');
console.log('══════════════════════════════════════════════════════════════\n');

const thresholds = [30, 40, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 140];
const thresholdResults = [];

for (const thresh of thresholds) {
  const subset = allGalaxies.filter(g => g.Vflat >= thresh);
  if (subset.length < 10) continue;

  const xArr = subset.map(g => g.fgas);
  const yArr = subset.map(g => g.logOMD);
  const r = pearsonR(xArr, yArr);
  const loo = looR2(xArr, yArr);
  const fit = olsRegress(xArr.map(v => [1, v]), yArr);
  const slope = fit ? fit.beta[1] : NaN;

  console.log('  Vflat >= ' + thresh + ' km/s: N=' + subset.length + ', r=' + r.toFixed(3) + ', LOO=' + loo.toFixed(3) + ', slope=' + slope.toFixed(3));
  thresholdResults.push({ threshold: thresh, n: subset.length, r: parseFloat(r.toFixed(3)), loo: parseFloat(loo.toFixed(3)), slope: parseFloat(slope.toFixed(3)) });
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: SLIDING UPPER THRESHOLD');
console.log('  Keep all below threshold');
console.log('══════════════════════════════════════════════════════════════\n');

const upperThresholds = [60, 70, 80, 90, 100, 120, 150, 200, 250, 300];
const upperResults = [];

for (const thresh of upperThresholds) {
  const subset = allGalaxies.filter(g => g.Vflat <= thresh);
  if (subset.length < 10) continue;

  const xArr = subset.map(g => g.fgas);
  const yArr = subset.map(g => g.logOMD);
  const r = pearsonR(xArr, yArr);
  const loo = looR2(xArr, yArr);
  const fit = olsRegress(xArr.map(v => [1, v]), yArr);

  console.log('  Vflat <= ' + thresh + ' km/s: N=' + subset.length + ', r=' + r.toFixed(3) + ', LOO=' + loo.toFixed(3));
  upperResults.push({ threshold: thresh, n: subset.length, r: parseFloat(r.toFixed(3)), loo: parseFloat(loo.toFixed(3)) });
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: SLIDING WINDOW (width = 80 km/s)');
console.log('══════════════════════════════════════════════════════════════\n');

const windowWidth = 80;
const windowResults = [];
for (let center = 50; center <= 250; center += 10) {
  const lo = center - windowWidth / 2;
  const hi = center + windowWidth / 2;
  const subset = allGalaxies.filter(g => g.Vflat >= lo && g.Vflat <= hi);
  if (subset.length < 10) continue;

  const xArr = subset.map(g => g.fgas);
  const yArr = subset.map(g => g.logOMD);
  const r = pearsonR(xArr, yArr);
  const loo = looR2(xArr, yArr);

  console.log('  [' + lo + ', ' + hi + '] km/s: N=' + subset.length + ', r=' + r.toFixed(3) + ', LOO=' + loo.toFixed(3));
  windowResults.push({ center, lo, hi, n: subset.length, r: parseFloat(r.toFixed(3)), loo: parseFloat(loo.toFixed(3)) });
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: CALIBRATION DRIFT WITH THRESHOLD');
console.log('  Does the slope/intercept shift as we move the boundary?');
console.log('══════════════════════════════════════════════════════════════\n');

for (const thresh of [50, 60, 70, 80, 100]) {
  const subset = allGalaxies.filter(g => g.Vflat >= thresh);
  if (subset.length < 10) continue;
  const fit = olsRegress(subset.map(g => [1, g.fgas]), subset.map(g => g.logOMD));
  if (fit) {
    console.log('  Vflat >= ' + thresh + ': slope=' + fit.beta[1].toFixed(4) + ', intercept=' + fit.beta[0].toFixed(4) + ', N=' + subset.length);
  }
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: TRANSITION DIAGNOSTICS AT BOUNDARY');
console.log('  What changes physically around Vflat ~ 70?');
console.log('══════════════════════════════════════════════════════════════\n');

const lowV = allGalaxies.filter(g => g.Vflat < 70);
const highV = allGalaxies.filter(g => g.Vflat >= 70);

console.log('  Low-Vflat (< 70): N=' + lowV.length);
console.log('    fgas mean = ' + mean(lowV.map(g => g.fgas)).toFixed(3) + ', sd = ' + (lowV.length > 1 ? Math.sqrt(lowV.map(g => g.fgas).reduce((s, v) => s + (v - mean(lowV.map(g => g.fgas))) ** 2, 0) / (lowV.length - 1)).toFixed(3) : 'N/A'));
console.log('    logOMD mean = ' + mean(lowV.map(g => g.logOMD)).toFixed(3));
console.log('    fgas range = [' + Math.min(...lowV.map(g => g.fgas)).toFixed(3) + ', ' + Math.max(...lowV.map(g => g.fgas)).toFixed(3) + ']');

console.log('  High-Vflat (>= 70): N=' + highV.length);
console.log('    fgas mean = ' + mean(highV.map(g => g.fgas)).toFixed(3) + ', sd = ' + (highV.length > 1 ? Math.sqrt(highV.map(g => g.fgas).reduce((s, v) => s + (v - mean(highV.map(g => g.fgas))) ** 2, 0) / (highV.length - 1)).toFixed(3) : 'N/A'));
console.log('    logOMD mean = ' + mean(highV.map(g => g.logOMD)).toFixed(3));
console.log('    fgas range = [' + Math.min(...highV.map(g => g.fgas)).toFixed(3) + ', ' + Math.max(...highV.map(g => g.fgas)).toFixed(3) + ']');

const fgasOverlap = lowV.filter(g => g.fgas >= Math.min(...highV.map(h => h.fgas)) && g.fgas <= Math.max(...highV.map(h => h.fgas))).length;
console.log('  Low-V galaxies in high-V fgas range: ' + fgasOverlap + '/' + lowV.length);

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const transition = thresholdResults.find(t => t.loo > 0.3 && t.threshold <= 70);
const collapse = thresholdResults.find(t => t.loo < 0.1 && t.threshold >= 40);
const bestWindow = windowResults.sort((a, b) => b.loo - a.loo)[0];

let boundaryDesc;
const stableAbove = thresholdResults.filter(t => t.threshold >= 60 && t.loo > 0.3);
const failsBelow = thresholdResults.filter(t => t.threshold <= 50 && t.loo < 0.2);

if (stableAbove.length > 0 && failsBelow.length > 0) {
  const sharpBoundary = stableAbove[0].threshold;
  boundaryDesc = 'Sharp boundary near Vflat ~ ' + sharpBoundary + ' km/s. Relationship is strong above, collapses below.';
} else {
  boundaryDesc = 'Gradual transition. No sharp boundary identified.';
}

console.log('  Boundary: ' + boundaryDesc);
if (bestWindow) console.log('  Best window: [' + bestWindow.lo + ', ' + bestWindow.hi + '] km/s (r=' + bestWindow.r + ', LOO=' + bestWindow.loo + ')');
console.log('  Calibration stable above threshold? ' + (thresholdResults.filter(t => t.threshold >= 70).every(t => Math.abs(t.slope - thresholdResults.find(tt => tt.threshold === 70).slope) < 0.3) ? 'YES' : 'NO'));

const output = {
  phase: 111,
  title: 'Regime Boundary Mapping',
  totalGalaxies: allGalaxies.length,
  lowerThresholds: thresholdResults,
  upperThresholds: upperResults,
  slidingWindows: windowResults,
  boundaryDescription: boundaryDesc
};

const outPath = path.join(__dirname, '..', 'public', 'phase111-boundary-mapping.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);

#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 101: NULL GEOMETRIC COUPLING TEST');
console.log('  Can smooth curves without real concentration produce r=0.85?');
console.log('  If null r << real r, the law is genuine, not geometric artifact.');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;
const N_NULL_TRIALS = 1000;
const SEED = 42;

let rngState = SEED;
function seededRandom() {
  rngState = (rngState * 1664525 + 1013904223) & 0x7fffffff;
  return rngState / 0x7fffffff;
}
function gaussRandom() {
  let u, v, s;
  do { u = 2 * seededRandom() - 1; v = 2 * seededRandom() - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
  return u * Math.sqrt(-2 * Math.log(s) / s);
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
  const ss_res = y.reduce((s, v, i) => s + (v - yhat[i]) ** 2, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0 };
}

function looR2(X_fn, y_fn, data) {
  const n = data.length;
  let ss_res = 0, ss_tot = 0;
  const yAll = data.map(y_fn);
  const yMean = mean(yAll);
  for (let i = 0; i < n; i++) {
    const train = data.filter((_, j) => j !== i);
    const Xtrain = train.map(X_fn);
    const ytrain = train.map(y_fn);
    const fit = olsRegress(Xtrain, ytrain);
    if (!fit) return NaN;
    const xTest = X_fn(data[i]);
    const yPred = xTest.reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss_res += (yAll[i] - yPred) ** 2;
    ss_tot += (yAll[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
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

const realGalaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0) continue;
  if (t1.Vflat < VFLAT_BREAK) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    pts.push({ r: pt.rad, vobs: pt.vobs });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const innerPts = sorted.slice(0, half);
  const outerPts = sorted.slice(half);

  const innerVmax = Math.max(...innerPts.map(p => p.vobs));
  const Vflat = t1.Vflat;

  const logR_outer = outerPts.map(p => Math.log10(p.r));
  const logV_outer = outerPts.map(p => Math.log10(p.vobs));
  const mr = mean(logR_outer), mv = mean(logV_outer);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR_outer.length; i++) { sxy += (logR_outer[i] - mr) * (logV_outer[i] - mv); sxx += (logR_outer[i] - mr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const logRatio = Math.log10(Vflat / innerVmax);
  if (!isFinite(logRatio) || !isFinite(outerSlope)) continue;

  const noiseLevel = (() => {
    if (sorted.length < 5) return 0.02;
    const residuals = [];
    for (let i = 2; i < sorted.length - 2; i++) {
      const localMean = (sorted[i - 1].vobs + sorted[i].vobs + sorted[i + 1].vobs) / 3;
      residuals.push(Math.abs(sorted[i].vobs - localMean) / localMean);
    }
    return mean(residuals) || 0.02;
  })();

  realGalaxies.push({
    name, nPts: sorted.length,
    rMin: sorted[0].r, rMax: sorted[sorted.length - 1].r,
    innerVmax, Vflat, outerSlope, logRatio,
    velocities: sorted.map(p => p.vobs),
    radii: sorted.map(p => p.r),
    noiseLevel, half
  });
}

console.log('  Real galaxies (Vflat >= ' + VFLAT_BREAK + '): N=' + realGalaxies.length);

const realRatios = realGalaxies.map(g => g.logRatio);
const realSlopes = realGalaxies.map(g => g.outerSlope);
const realR = pearsonR(realRatios, realSlopes);

const ratioModel = olsRegress(realRatios.map(x => [1, x]), realSlopes);
const ratioLOO = looR2(g => [1, g.logRatio], g => g.outerSlope, realGalaxies);

console.log('\n── REAL DATA BASELINE ──');
console.log('  r(logRatio, outerSlope) = ' + realR.toFixed(4));
console.log('  R² = ' + ratioModel.r2.toFixed(4));
console.log('  LOO R² = ' + ratioLOO.toFixed(4));
console.log('  Slope = ' + ratioModel.beta[1].toFixed(4) + ', Intercept = ' + ratioModel.beta[0].toFixed(4));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 1: NULL SMOOTH CURVES (no real concentration)');
console.log('  Generate smooth monotonic curves with matched properties');
console.log('  but random inner shape. Real concentration is destroyed.');
console.log('══════════════════════════════════════════════════════════════\n');

function generateNullCurve(galaxy) {
  const n = galaxy.nPts;
  const rMin = galaxy.rMin;
  const rMax = galaxy.rMax;
  const vFlat = galaxy.Vflat;
  const noise = galaxy.noiseLevel;

  const radii = [];
  for (let i = 0; i < n; i++) {
    radii.push(rMin + (rMax - rMin) * (i / (n - 1)));
  }

  const rScale = rMin + (rMax - rMin) * (0.1 + 0.4 * seededRandom());
  const alpha = 0.3 + 1.7 * seededRandom();
  const vStart = vFlat * (0.3 + 0.5 * seededRandom());

  const vClean = radii.map(r => {
    const x = r / rScale;
    const rise = 1 - Math.exp(-Math.pow(x, alpha));
    return vStart + (vFlat - vStart) * rise;
  });

  const velocities = vClean.map(v => v * (1 + noise * gaussRandom()));

  return { radii, velocities, n };
}

function measureNullCurve(nc, Vflat) {
  const half = Math.floor(nc.n / 2);
  const innerV = nc.velocities.slice(0, half);
  const outerR = nc.radii.slice(half);
  const outerV = nc.velocities.slice(half);

  const innerVmax = Math.max(...innerV);
  const logRatio = Math.log10(Vflat / innerVmax);

  const logR = outerR.map(r => Math.log10(r));
  const logV = outerV.map(v => Math.log10(Math.max(v, 1)));
  const mr = mean(logR), mv = mean(logV);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mr) * (logV[i] - mv); sxx += (logR[i] - mr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  return { logRatio, outerSlope };
}

const nullCorrelations = [];

for (let trial = 0; trial < N_NULL_TRIALS; trial++) {
  const nullRatios = [];
  const nullSlopes = [];

  for (const galaxy of realGalaxies) {
    const nc = generateNullCurve(galaxy);
    const meas = measureNullCurve(nc, galaxy.Vflat);
    if (isFinite(meas.logRatio) && isFinite(meas.outerSlope)) {
      nullRatios.push(meas.logRatio);
      nullSlopes.push(meas.outerSlope);
    }
  }

  if (nullRatios.length >= 20) {
    const r = pearsonR(nullRatios, nullSlopes);
    nullCorrelations.push(r);
  }
}

const nullMean = mean(nullCorrelations);
const nullSD = sd(nullCorrelations);
const nullMax = Math.max(...nullCorrelations);
const nullMin = Math.min(...nullCorrelations);
const nullMedian = [...nullCorrelations].sort((a, b) => a - b)[Math.floor(nullCorrelations.length / 2)];
const nullP95 = [...nullCorrelations].sort((a, b) => a - b)[Math.floor(nullCorrelations.length * 0.95)];
const nullP99 = [...nullCorrelations].sort((a, b) => a - b)[Math.floor(nullCorrelations.length * 0.99)];
const pValue = nullCorrelations.filter(r => Math.abs(r) >= Math.abs(realR)).length / nullCorrelations.length;
const zScore = nullSD > 0 ? (Math.abs(realR) - nullMean) / nullSD : Infinity;

console.log('  Null trials: ' + N_NULL_TRIALS);
console.log('  Null r: mean=' + nullMean.toFixed(4) + ', sd=' + nullSD.toFixed(4));
console.log('  Null r: median=' + nullMedian.toFixed(4));
console.log('  Null r: 95th pctile=' + nullP95.toFixed(4) + ', 99th pctile=' + nullP99.toFixed(4));
console.log('  Null r: max=' + nullMax.toFixed(4) + ', min=' + nullMin.toFixed(4));
console.log('  Real |r| = ' + Math.abs(realR).toFixed(4));
console.log('  z-score = ' + zScore.toFixed(2));
console.log('  p-value (null >= real) = ' + pValue.toFixed(6));
console.log('  VERDICT: ' + (pValue < 0.001 ? 'PASS - real correlation FAR exceeds null geometric coupling' :
  pValue < 0.05 ? 'MARGINAL - real exceeds null but not decisively' :
  'FAIL - null curves can reproduce the correlation'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: SHUFFLED INNER HALVES');
console.log('  Keep real outer halves but swap inner halves between galaxies');
console.log('  If inner half carries real physics, shuffling should destroy r');
console.log('══════════════════════════════════════════════════════════════\n');

const shuffledCorrelations = [];
for (let trial = 0; trial < N_NULL_TRIALS; trial++) {
  const indices = realGalaxies.map((_, i) => i);
  for (let i = indices.length - 1; i > 0; i--) {
    const j = Math.floor(seededRandom() * (i + 1));
    [indices[i], indices[j]] = [indices[j], indices[i]];
  }

  const shuffRatios = [];
  const shuffSlopes = [];

  for (let k = 0; k < realGalaxies.length; k++) {
    const donorInner = realGalaxies[indices[k]];
    const recipientOuter = realGalaxies[k];

    const logRatio = Math.log10(recipientOuter.Vflat / donorInner.innerVmax);
    if (isFinite(logRatio)) {
      shuffRatios.push(logRatio);
      shuffSlopes.push(recipientOuter.outerSlope);
    }
  }

  if (shuffRatios.length >= 20) {
    shuffledCorrelations.push(pearsonR(shuffRatios, shuffSlopes));
  }
}

const shuffMean = mean(shuffledCorrelations);
const shuffSD = sd(shuffledCorrelations);
const shuffMax = Math.max(...shuffledCorrelations);
const shuffP95 = [...shuffledCorrelations].sort((a, b) => a - b)[Math.floor(shuffledCorrelations.length * 0.95)];
const shuffP99 = [...shuffledCorrelations].sort((a, b) => a - b)[Math.floor(shuffledCorrelations.length * 0.99)];
const shuffPvalue = shuffledCorrelations.filter(r => Math.abs(r) >= Math.abs(realR)).length / shuffledCorrelations.length;
const shuffZscore = shuffSD > 0 ? (Math.abs(realR) - shuffMean) / shuffSD : Infinity;

console.log('  Shuffled r: mean=' + shuffMean.toFixed(4) + ', sd=' + shuffSD.toFixed(4));
console.log('  Shuffled r: 95th pctile=' + shuffP95.toFixed(4) + ', 99th pctile=' + shuffP99.toFixed(4));
console.log('  Shuffled r: max=' + shuffMax.toFixed(4));
console.log('  Real |r| = ' + Math.abs(realR).toFixed(4));
console.log('  z-score = ' + shuffZscore.toFixed(2));
console.log('  p-value = ' + shuffPvalue.toFixed(6));
console.log('  VERDICT: ' + (shuffPvalue < 0.001 ? 'PASS - shuffling destroys correlation, inner half carries real signal' :
  shuffPvalue < 0.05 ? 'MARGINAL' :
  'FAIL - correlation survives shuffling (geometric artifact)'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: RADIUS-PRESERVING VELOCITY SHUFFLE');
console.log('  Keep each galaxy\'s radii but assign velocities from a');
console.log('  different galaxy with similar nPts. Preserves RC shape');
console.log('  statistics but breaks galaxy-specific concentration.');
console.log('══════════════════════════════════════════════════════════════\n');

const crossCorrelations = [];
for (let trial = 0; trial < N_NULL_TRIALS; trial++) {
  const indices = realGalaxies.map((_, i) => i);
  for (let i = indices.length - 1; i > 0; i--) {
    const j = Math.floor(seededRandom() * (i + 1));
    [indices[i], indices[j]] = [indices[j], indices[i]];
  }

  const crossRatios = [];
  const crossSlopes = [];

  for (let k = 0; k < realGalaxies.length; k++) {
    const donor = realGalaxies[indices[k]];
    const target = realGalaxies[k];

    const donorVels = donor.velocities;
    const targetRadii = target.radii;
    const n = Math.min(donorVels.length, targetRadii.length);
    if (n < 8) continue;

    const scaledVels = [];
    for (let i = 0; i < n; i++) {
      const frac = i / (donorVels.length - 1);
      const idx = frac * (donorVels.length - 1);
      const lo = Math.floor(idx);
      const hi = Math.min(lo + 1, donorVels.length - 1);
      const t = idx - lo;
      const v = donorVels[lo] * (1 - t) + donorVels[hi] * t;
      const scale = target.Vflat / donor.Vflat;
      scaledVels.push(v * scale);
    }

    const half = Math.floor(n / 2);
    const innerVmax = Math.max(...scaledVels.slice(0, half));
    const logRatio = Math.log10(target.Vflat / innerVmax);

    const logR = targetRadii.slice(half, n).map(r => Math.log10(r));
    const logV = scaledVels.slice(half, n).map(v => Math.log10(Math.max(v, 1)));
    const mr = mean(logR), mv = mean(logV);
    let sxy = 0, sxx2 = 0;
    for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mr) * (logV[i] - mv); sxx2 += (logR[i] - mr) ** 2; }
    const outerSlope = sxx2 > 0 ? sxy / sxx2 : 0;

    if (isFinite(logRatio) && isFinite(outerSlope)) {
      crossRatios.push(logRatio);
      crossSlopes.push(outerSlope);
    }
  }

  if (crossRatios.length >= 20) {
    crossCorrelations.push(pearsonR(crossRatios, crossSlopes));
  }
}

const crossMean = mean(crossCorrelations);
const crossSD = sd(crossCorrelations);
const crossMax = Math.max(...crossCorrelations);
const crossP95 = [...crossCorrelations].sort((a, b) => a - b)[Math.floor(crossCorrelations.length * 0.95)];
const crossPvalue = crossCorrelations.filter(r => Math.abs(r) >= Math.abs(realR)).length / crossCorrelations.length;
const crossZscore = crossSD > 0 ? (Math.abs(realR) - crossMean) / crossSD : Infinity;

console.log('  Cross-assigned r: mean=' + crossMean.toFixed(4) + ', sd=' + crossSD.toFixed(4));
console.log('  Cross-assigned r: 95th pctile=' + crossP95.toFixed(4));
console.log('  Cross-assigned r: max=' + crossMax.toFixed(4));
console.log('  Real |r| = ' + Math.abs(realR).toFixed(4));
console.log('  z-score = ' + crossZscore.toFixed(2));
console.log('  p-value = ' + crossPvalue.toFixed(6));
console.log('  VERDICT: ' + (crossPvalue < 0.001 ? 'PASS - cross-assignment destroys correlation' :
  crossPvalue < 0.05 ? 'MARGINAL' :
  'FAIL - correlation survives cross-assignment'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: MONOTONIC NULL (strongest geometric null)');
console.log('  Curves guaranteed monotonically rising to Vflat, with');
console.log('  the SAME radial grid as the real galaxy. Only inner');
console.log('  concentration is randomized. This is the hardest test.');
console.log('══════════════════════════════════════════════════════════════\n');

function generateMonotonicNull(galaxy) {
  const radii = [...galaxy.radii];
  const n = radii.length;
  const vFlat = galaxy.Vflat;
  const noise = galaxy.noiseLevel;

  const turnoverFrac = 0.1 + 0.6 * seededRandom();
  const turnoverIdx = Math.max(2, Math.floor(n * turnoverFrac));
  const vStart = vFlat * (0.15 + 0.35 * seededRandom());
  const sharpness = 0.5 + 3.0 * seededRandom();

  const velocities = [];
  for (let i = 0; i < n; i++) {
    let v;
    if (i <= turnoverIdx) {
      const frac = i / turnoverIdx;
      const shaped = Math.pow(frac, 1 / sharpness);
      v = vStart + (vFlat - vStart) * shaped;
    } else {
      const drift = (seededRandom() - 0.5) * 0.02;
      v = vFlat * (1 + drift * (i - turnoverIdx) / (n - turnoverIdx));
    }
    v *= (1 + noise * 0.5 * gaussRandom());
    velocities.push(Math.max(v, 5));
  }

  return { radii, velocities, n };
}

const monoCorrelations = [];
for (let trial = 0; trial < N_NULL_TRIALS; trial++) {
  const monoRatios = [];
  const monoSlopes = [];

  for (const galaxy of realGalaxies) {
    const nc = generateMonotonicNull(galaxy);
    const meas = measureNullCurve(nc, galaxy.Vflat);
    if (isFinite(meas.logRatio) && isFinite(meas.outerSlope)) {
      monoRatios.push(meas.logRatio);
      monoSlopes.push(meas.outerSlope);
    }
  }

  if (monoRatios.length >= 20) {
    monoCorrelations.push(pearsonR(monoRatios, monoSlopes));
  }
}

const monoMean = mean(monoCorrelations);
const monoSD = sd(monoCorrelations);
const monoMax = Math.max(...monoCorrelations);
const monoP95 = [...monoCorrelations].sort((a, b) => a - b)[Math.floor(monoCorrelations.length * 0.95)];
const monoP99 = [...monoCorrelations].sort((a, b) => a - b)[Math.floor(monoCorrelations.length * 0.99)];
const monoPvalue = monoCorrelations.filter(r => Math.abs(r) >= Math.abs(realR)).length / monoCorrelations.length;
const monoZscore = monoSD > 0 ? (Math.abs(realR) - monoMean) / monoSD : Infinity;

console.log('  Monotonic null r: mean=' + monoMean.toFixed(4) + ', sd=' + monoSD.toFixed(4));
console.log('  Monotonic null r: 95th pctile=' + monoP95.toFixed(4) + ', 99th pctile=' + monoP99.toFixed(4));
console.log('  Monotonic null r: max=' + monoMax.toFixed(4));
console.log('  Real |r| = ' + Math.abs(realR).toFixed(4));
console.log('  z-score = ' + monoZscore.toFixed(2));
console.log('  p-value = ' + monoPvalue.toFixed(6));
console.log('  VERDICT: ' + (monoPvalue < 0.001 ? 'PASS - even matching radial grid, null cannot reach real r' :
  monoPvalue < 0.05 ? 'MARGINAL' :
  'FAIL - monotonic null reproduces correlation'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: GEOMETRIC FLOOR ESTIMATE');
console.log('  What correlation do we expect purely from the fact that');
console.log('  innerVmax and outerSlope come from the same curve?');
console.log('  Measures the "geometric floor" to subtract from signal.');
console.log('══════════════════════════════════════════════════════════════\n');

const geoFloors = [nullMean, shuffMean, crossMean, monoMean];
const geoFloorMax = Math.max(...geoFloors.map(Math.abs));
const signalAboveFloor = Math.abs(realR) - geoFloorMax;
const signalRatio = geoFloorMax > 0 ? Math.abs(realR) / geoFloorMax : Infinity;

console.log('  Geometric floor estimates:');
console.log('    Smooth null:      |r| = ' + Math.abs(nullMean).toFixed(4));
console.log('    Shuffled inner:   |r| = ' + Math.abs(shuffMean).toFixed(4));
console.log('    Cross-assigned:   |r| = ' + Math.abs(crossMean).toFixed(4));
console.log('    Monotonic null:   |r| = ' + Math.abs(monoMean).toFixed(4));
console.log('  Maximum geometric floor: |r| = ' + geoFloorMax.toFixed(4));
console.log('  Real signal: |r| = ' + Math.abs(realR).toFixed(4));
console.log('  Signal above floor: ' + signalAboveFloor.toFixed(4));
console.log('  Signal/floor ratio: ' + signalRatio.toFixed(2) + 'x');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  COMBINED VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const tests = [
  { name: 'Smooth null curves', pval: pValue, z: zScore },
  { name: 'Shuffled inner halves', pval: shuffPvalue, z: shuffZscore },
  { name: 'Cross-assigned velocities', pval: crossPvalue, z: crossZscore },
  { name: 'Monotonic null (hardest)', pval: monoPvalue, z: monoZscore },
];

let passed = 0;
for (const t of tests) {
  const status = t.pval < 0.001 ? 'PASS' : t.pval < 0.05 ? 'MARGINAL' : 'FAIL';
  if (t.pval < 0.05) passed++;
  console.log('  ' + status + ': ' + t.name + ' (p=' + t.pval.toFixed(6) + ', z=' + t.z.toFixed(1) + ')');
}

console.log('\n  Passed: ' + passed + '/4 null tests');
console.log('  Geometric floor: ' + geoFloorMax.toFixed(4) + ' vs real ' + Math.abs(realR).toFixed(4) + ' (' + signalRatio.toFixed(1) + 'x above floor)');

const levelVerdict = passed === 4 ? 'LEVEL 1 CONFIRMED: The relation CANNOT be produced by geometric coupling alone.'
  : passed >= 3 ? 'LEVEL 1 LIKELY: Strong evidence against geometric artifact, one test marginal.'
  : passed >= 2 ? 'LEVEL 1 UNCERTAIN: Mixed results, geometric coupling partially explains the signal.'
  : 'LEVEL 1 FAILS: Geometric coupling can explain the observed correlation.';

console.log('\n  ' + levelVerdict);

const results = {
  phase: 101,
  title: 'Null Geometric Coupling Test',
  description: 'Tests whether the ratio-slope correlation can be reproduced by smooth curves without real inner concentration',
  nGalaxies: realGalaxies.length,
  realCorrelation: {
    r: parseFloat(realR.toFixed(4)),
    R2: parseFloat(ratioModel.r2.toFixed(4)),
    LOO_R2: parseFloat(ratioLOO.toFixed(4)),
    slope: parseFloat(ratioModel.beta[1].toFixed(4)),
    intercept: parseFloat(ratioModel.beta[0].toFixed(4))
  },
  nullTests: {
    nTrials: N_NULL_TRIALS,
    smoothNull: {
      meanR: parseFloat(nullMean.toFixed(4)),
      sdR: parseFloat(nullSD.toFixed(4)),
      medianR: parseFloat(nullMedian.toFixed(4)),
      p95: parseFloat(nullP95.toFixed(4)),
      p99: parseFloat(nullP99.toFixed(4)),
      maxR: parseFloat(nullMax.toFixed(4)),
      pValue: parseFloat(pValue.toFixed(6)),
      zScore: parseFloat(zScore.toFixed(2)),
      verdict: pValue < 0.001 ? 'PASS' : pValue < 0.05 ? 'MARGINAL' : 'FAIL'
    },
    shuffledInner: {
      meanR: parseFloat(shuffMean.toFixed(4)),
      sdR: parseFloat(shuffSD.toFixed(4)),
      p95: parseFloat(shuffP95.toFixed(4)),
      maxR: parseFloat(shuffMax.toFixed(4)),
      pValue: parseFloat(shuffPvalue.toFixed(6)),
      zScore: parseFloat(shuffZscore.toFixed(2)),
      verdict: shuffPvalue < 0.001 ? 'PASS' : shuffPvalue < 0.05 ? 'MARGINAL' : 'FAIL'
    },
    crossAssigned: {
      meanR: parseFloat(crossMean.toFixed(4)),
      sdR: parseFloat(crossSD.toFixed(4)),
      p95: parseFloat(crossP95.toFixed(4)),
      maxR: parseFloat(crossMax.toFixed(4)),
      pValue: parseFloat(crossPvalue.toFixed(6)),
      zScore: parseFloat(crossZscore.toFixed(2)),
      verdict: crossPvalue < 0.001 ? 'PASS' : crossPvalue < 0.05 ? 'MARGINAL' : 'FAIL'
    },
    monotonicNull: {
      meanR: parseFloat(monoMean.toFixed(4)),
      sdR: parseFloat(monoSD.toFixed(4)),
      p95: parseFloat(monoP95.toFixed(4)),
      p99: parseFloat(monoP99.toFixed(4)),
      maxR: parseFloat(monoMax.toFixed(4)),
      pValue: parseFloat(monoPvalue.toFixed(6)),
      zScore: parseFloat(monoZscore.toFixed(2)),
      verdict: monoPvalue < 0.001 ? 'PASS' : monoPvalue < 0.05 ? 'MARGINAL' : 'FAIL'
    }
  },
  geometricFloor: {
    maxFloor: parseFloat(geoFloorMax.toFixed(4)),
    signalAboveFloor: parseFloat(signalAboveFloor.toFixed(4)),
    signalFloorRatio: parseFloat(signalRatio.toFixed(2))
  },
  combinedVerdict: {
    passed: passed,
    total: 4,
    level1: levelVerdict
  }
};

const outPath = path.join(__dirname, '..', 'public', 'phase101-null-geometric-coupling.json');
fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
console.log('\n  Results saved to: ' + outPath);
console.log('\n  DONE.');

#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 97: IS THE DIRECT LAW GENUINE OR GEOMETRIC?');
console.log('  Final circularity check on Vflat/InnerVmax → outerSlope');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

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

function partialCorrMulti(xFn, yFn, zFns, data) {
  const xv = data.map(xFn), yv = data.map(yFn);
  const Z = data.map(g => [1, ...zFns.map(fn => fn(g))]);
  const fitXZ = olsRegress(Z, xv);
  const fitYZ = olsRegress(Z, yv);
  if (!fitXZ || !fitYZ) return NaN;
  return pearsonR(fitXZ.residuals, fitYZ.residuals);
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
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);

  const half = Math.floor(sorted.length / 2);
  const innerPts = sorted.slice(0, half);
  const outerPts = sorted.slice(half);

  const innerVmax = Math.max(...innerPts.map(p => p.vobs));

  const logR_outer = outerPts.map(p => Math.log10(p.r));
  const logV_outer = outerPts.map(p => Math.log10(p.vobs));
  const mr = mean(logR_outer), mv = mean(logV_outer);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR_outer.length; i++) { sxy += (logR_outer[i] - mr) * (logV_outer[i] - mv); sxx += (logR_outer[i] - mr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const rMax = sorted[sorted.length - 1].r;
  const r1Rdisk = t1.Rdisk > 0 ? t1.Rdisk : rMax * 0.3;
  const r2Rdisk = 2 * r1Rdisk;
  const r3Rdisk = 3 * r1Rdisk;

  const vAt1Rd = (() => {
    let best = sorted[0];
    for (const p of sorted) { if (Math.abs(p.r - r1Rdisk) < Math.abs(best.r - r1Rdisk)) best = p; }
    return best.vobs;
  })();

  const vAt2Rd = (() => {
    let best = sorted[0];
    for (const p of sorted) { if (Math.abs(p.r - r2Rdisk) < Math.abs(best.r - r2Rdisk)) best = p; }
    return best.vobs;
  })();

  const vPeak = Math.max(...sorted.map(p => p.vobs));
  const vLast3 = mean(sorted.slice(-3).map(p => p.vobs));

  const logRatio = Math.log10(t1.Vflat > 0 && innerVmax > 0 ? t1.Vflat / innerVmax : 1);
  const logRatioV1 = Math.log10(t1.Vflat > 0 && vAt1Rd > 0 ? t1.Vflat / vAt1Rd : 1);
  const logRatioV2 = Math.log10(t1.Vflat > 0 && vAt2Rd > 0 ? t1.Vflat / vAt2Rd : 1);
  const logRatioVpeak = Math.log10(t1.Vflat > 0 && vPeak > 0 ? t1.Vflat / vPeak : 1);

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);

  allGalaxies.push({
    name, nPts: pts.length,
    innerVmax, logInnerVmax: Math.log10(innerVmax > 0 ? innerVmax : 1),
    outerSlope,
    logRatio, logRatioV1, logRatioV2, logRatioVpeak,
    vAt1Rd, vAt2Rd, vPeak, vLast3,
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    Vflat: t1.Vflat, fgas,
    logMHI: Math.log10(t1.MHI),
    T: t1.T, Q: t1.Q, Rdisk: t1.Rdisk,
  });
}

const hi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);
console.log('  High-regime: N=' + hi.length + '\n');

const structControls = [g => g.logSigma0, g => g.logVflat, g => g.fgas, g => g.logMHI];

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: SPLIT-HALF — innerVmax from inner, slope from outer');
console.log('  (Already done by design — innerVmax uses inner pts only,');
console.log('   outerSlope uses outer pts only. Confirm NO overlap.)');
console.log('══════════════════════════════════════════════════════════════\n');

const rDirect = pearsonR(hi.map(g => g.logInnerVmax), hi.map(g => g.outerSlope));
const rRatio = pearsonR(hi.map(g => g.logRatio), hi.map(g => g.outerSlope));

console.log('  By construction: innerVmax from inner half, outerSlope from outer half');
console.log('  NO radial overlap — this IS the split-half test.');
console.log('  r(logInnerVmax, outerSlope) = ' + rDirect.toFixed(3));
console.log('  r(log(Vflat/InnerVmax), outerSlope) = ' + rRatio.toFixed(3));

const rPartDirect = partialCorrMulti(g => g.logInnerVmax, g => g.outerSlope, structControls, hi);
const rPartRatio = partialCorrMulti(g => g.logRatio, g => g.outerSlope, structControls, hi);
console.log('  Partial r(logInnerVmax, outerSlope | struct) = ' + rPartDirect.toFixed(3));
console.log('  Partial r(logRatio, outerSlope | struct) = ' + rPartRatio.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: FIXED RADIUS PROXIES (not tailored to RC shape)');
console.log('  V(1 Rdisk), V(2 Rdisk), Vpeak as alternatives');
console.log('══════════════════════════════════════════════════════════════\n');

const proxies = [
  { name: 'V(1 Rdisk)', fn: g => Math.log10(g.vAt1Rd > 0 ? g.vAt1Rd : 1), ratioFn: g => g.logRatioV1 },
  { name: 'V(2 Rdisk)', fn: g => Math.log10(g.vAt2Rd > 0 ? g.vAt2Rd : 1), ratioFn: g => g.logRatioV2 },
  { name: 'Vpeak (entire RC)', fn: g => Math.log10(g.vPeak > 0 ? g.vPeak : 1), ratioFn: g => g.logRatioVpeak },
  { name: 'InnerVmax (half)', fn: g => g.logInnerVmax, ratioFn: g => g.logRatio },
];

for (const px of proxies) {
  const rRaw = pearsonR(hi.map(px.fn), hi.map(g => g.outerSlope));
  const rPart = partialCorrMulti(px.fn, g => g.outerSlope, structControls, hi);
  const rRatRaw = pearsonR(hi.map(px.ratioFn), hi.map(g => g.outerSlope));
  const rRatPart = partialCorrMulti(px.ratioFn, g => g.outerSlope, structControls, hi);

  console.log('  ' + px.name + ':');
  console.log('    r(logV, outerSlope) = ' + rRaw.toFixed(3) + '   partial = ' + rPart.toFixed(3));
  console.log('    r(log(Vflat/V), outerSlope) = ' + rRatRaw.toFixed(3) + '   partial = ' + rRatPart.toFixed(3));

  const X_fn = g => [1, px.fn(g), g.logSigma0, g.logVflat, g.fgas, g.logMHI];
  const loo = looR2(X_fn, g => g.outerSlope, hi);
  console.log('    LOO R² (struct + logV): ' + loo.toFixed(3));
  console.log();
}

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: ADDITIONAL CONTROLS — morphology, nPts, quality');
console.log('══════════════════════════════════════════════════════════════\n');

const extControls = [...structControls, g => g.T, g => g.nPts, g => g.Q];
const rExtRatio = partialCorrMulti(g => g.logRatio, g => g.outerSlope, extControls, hi);
const rExtDirect = partialCorrMulti(g => g.logInnerVmax, g => g.outerSlope, extControls, hi);
console.log('  Partial r(logRatio, outerSlope | struct + T + nPts + Q) = ' + rExtRatio.toFixed(3));
console.log('  Partial r(logInnerVmax, outerSlope | struct + T + nPts + Q) = ' + rExtDirect.toFixed(3));

const X_ext = g => [1, g.logInnerVmax, g.logSigma0, g.logVflat, g.fgas, g.logMHI, g.T, g.nPts, g.Q];
const looExt = looR2(X_ext, g => g.outerSlope, hi);
console.log('  LOO R² (struct + innerVmax + T + nPts + Q): ' + looExt.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: GEOMETRIC CIRCULARITY CHECK');
console.log('  If outerSlope ≈ log(Vflat/InnerVmax) is just geometry,');
console.log('  then ANY inner V would work equally well.');
console.log('  Compare: random inner point vs innerVmax');
console.log('══════════════════════════════════════════════════════════════\n');

const rRandResults = [];
for (let trial = 0; trial < 50; trial++) {
  const randData = hi.map(g => {
    const rIdx = Math.floor(Math.random() * Math.floor(g.nPts / 2));
    return { ...g, logRandV: Math.log10(g.vAt1Rd > 0 ? g.vAt1Rd : 1) };
  });
}

const rV1 = pearsonR(hi.map(g => Math.log10(g.vAt1Rd)), hi.map(g => g.outerSlope));
const rVmax = pearsonR(hi.map(g => g.logInnerVmax), hi.map(g => g.outerSlope));
const rVflat = pearsonR(hi.map(g => g.logVflat), hi.map(g => g.outerSlope));

console.log('  Raw correlations with outerSlope:');
console.log('    logV(1 Rdisk): ' + rV1.toFixed(3));
console.log('    logInnerVmax:  ' + rVmax.toFixed(3));
console.log('    logVflat:      ' + rVflat.toFixed(3));
console.log();
console.log('  If purely geometric, all inner velocities would work equally.');
console.log('  Differences between proxies → physical content, not geometry.');

const spread = Math.abs(rVmax - rV1);
if (spread > 0.05) {
  console.log('  Spread = ' + spread.toFixed(3) + ' → NOT purely geometric');
} else {
  console.log('  Spread = ' + spread.toFixed(3) + ' → possibly geometric');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: THE RATIO TEST — is log(Vflat/InnerVmax) the');
console.log('  actual variable, or do Vflat and InnerVmax contribute');
console.log('  independently?');
console.log('══════════════════════════════════════════════════════════════\n');

const M_ratio = g => [1, g.logRatio];
const M_both = g => [1, g.logVflat, g.logInnerVmax];
const M_ratioFull = g => [1, g.logRatio, g.logSigma0, g.fgas, g.logMHI];
const M_bothFull = g => [1, g.logVflat, g.logInnerVmax, g.logSigma0, g.fgas, g.logMHI];

const looRatioOnly = looR2(M_ratio, g => g.outerSlope, hi);
const looBothOnly = looR2(M_both, g => g.outerSlope, hi);
const looRatioFull = looR2(M_ratioFull, g => g.outerSlope, hi);
const looBothFull = looR2(M_bothFull, g => g.outerSlope, hi);

console.log('  Minimal models:');
console.log('    log(Vflat/InnerVmax) only:    LOO R² = ' + looRatioOnly.toFixed(3));
console.log('    logVflat + logInnerVmax:       LOO R² = ' + looBothOnly.toFixed(3));
console.log('  Full models (+ struct):');
console.log('    log(Vflat/InnerVmax) + struct: LOO R² = ' + looRatioFull.toFixed(3));
console.log('    logVflat + logInnerVmax + struct: LOO R² = ' + looBothFull.toFixed(3));

if (looBothFull > looRatioFull + 0.02) {
  console.log('  → Vflat and InnerVmax contribute INDEPENDENTLY (not just ratio)');
} else {
  console.log('  → The RATIO is sufficient (Vflat and InnerVmax do not add independently)');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 6: DOES THE LAW WORK BELOW 70 km/s TOO?');
console.log('══════════════════════════════════════════════════════════════\n');

const lo = allGalaxies.filter(g => g.Vflat < VFLAT_BREAK);
console.log('  Low-regime: N=' + lo.length);

if (lo.length >= 15) {
  const rLoRatio = pearsonR(lo.map(g => g.logRatio), lo.map(g => g.outerSlope));
  const rLoVmax = pearsonR(lo.map(g => g.logInnerVmax), lo.map(g => g.outerSlope));
  console.log('  r(logRatio, outerSlope) = ' + rLoRatio.toFixed(3));
  console.log('  r(logInnerVmax, outerSlope) = ' + rLoVmax.toFixed(3));

  const M_loFull = g => [1, g.logInnerVmax, g.logSigma0, g.logVflat, g.fgas, g.logMHI];
  if (lo.length > 8) {
    const fitLo = olsRegress(lo.map(M_loFull), lo.map(g => g.outerSlope));
    if (fitLo) {
      console.log('  R² (struct + innerVmax, low regime) = ' + fitLo.r2.toFixed(3));
      if (lo.length > 10) {
        const looLo = looR2(M_loFull, g => g.outerSlope, lo);
        console.log('  LOO R² (struct + innerVmax, low regime) = ' + looLo.toFixed(3));
      }
    }
  }
} else {
  console.log('  Too few galaxies for meaningful test');
}

const all = allGalaxies;
console.log('\n  ALL galaxies: N=' + all.length);
const rAllRatio = pearsonR(all.map(g => g.logRatio), all.map(g => g.outerSlope));
const rAllVmax = pearsonR(all.map(g => g.logInnerVmax), all.map(g => g.outerSlope));
console.log('  r(logRatio, outerSlope) = ' + rAllRatio.toFixed(3));
console.log('  r(logInnerVmax, outerSlope) = ' + rAllVmax.toFixed(3));

const M_allFull = g => [1, g.logInnerVmax, g.logSigma0, g.logVflat, g.fgas, g.logMHI];
const fitAll = olsRegress(all.map(M_allFull), all.map(g => g.outerSlope));
const looAll = looR2(M_allFull, g => g.outerSlope, all);
console.log('  R² (struct + innerVmax, ALL) = ' + fitAll.r2.toFixed(3));
console.log('  LOO R² (struct + innerVmax, ALL) = ' + looAll.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 7: PERMUTATION TEST on the ratio');
console.log('══════════════════════════════════════════════════════════════\n');

const nPerm = 2000;
const realR = pearsonR(hi.map(g => g.logRatio), hi.map(g => g.outerSlope));
let nBetter = 0;
const yVals = hi.map(g => g.outerSlope);
for (let i = 0; i < nPerm; i++) {
  const shuffled = [...yVals].sort(() => Math.random() - 0.5);
  const pr = pearsonR(hi.map(g => g.logRatio), shuffled);
  if (Math.abs(pr) >= Math.abs(realR)) nBetter++;
}
console.log('  r(logRatio, outerSlope) = ' + realR.toFixed(3) + ', perm p = ' + (nBetter / nPerm).toFixed(4));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  FINAL SUMMARY');
console.log('══════════════════════════════════════════════════════════════\n');

const genuineMarkers = [];
if (Math.abs(rPartDirect) > 0.3) genuineMarkers.push('Strong partial r after structural controls');
if (Math.abs(rExtDirect) > 0.3) genuineMarkers.push('Survives morphology + quality + nPts controls');
if (spread > 0.05) genuineMarkers.push('Proxy-dependent (not purely geometric)');
if (nBetter / nPerm < 0.001) genuineMarkers.push('Permutation significant');

const geometricMarkers = [];
if (spread < 0.05) geometricMarkers.push('All inner velocities work equally (geometric)');
if (Math.abs(rPartDirect) < 0.15) geometricMarkers.push('Disappears after structural controls');

console.log('  GENUINE markers (' + genuineMarkers.length + '):');
for (const m of genuineMarkers) console.log('    ✓ ' + m);
console.log('  GEOMETRIC markers (' + geometricMarkers.length + '):');
for (const m of geometricMarkers) console.log('    ✗ ' + m);

let verdict;
if (genuineMarkers.length >= 3 && geometricMarkers.length === 0) {
  verdict = 'GENUINE';
  console.log('\n  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  GENUINE PHYSICAL LAW: innerVmax → outerSlope is real      ║');
  console.log('  ║  Not a geometric artifact of RC shape                       ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else if (geometricMarkers.length >= 2) {
  verdict = 'GEOMETRIC';
  console.log('\n  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  GEOMETRIC ARTIFACT: the relationship is tautological      ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else {
  verdict = 'MIXED';
  console.log('\n  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  MIXED: partially genuine, partially geometric             ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
}

const output = {
  phase: '97',
  title: 'Is the Direct Law Genuine or Geometric?',
  hiN: hi.length,
  splitHalf: {
    note: 'innerVmax from inner half, outerSlope from outer half — NO overlap by construction',
    rDirect: +rDirect.toFixed(3),
    rRatio: +rRatio.toFixed(3),
    rPartDirect: +rPartDirect.toFixed(3),
    rPartRatio: +rPartRatio.toFixed(3),
  },
  fixedRadiusProxies: proxies.map(px => ({
    name: px.name,
    rRaw: +pearsonR(hi.map(px.fn), hi.map(g => g.outerSlope)).toFixed(3),
    rPart: +partialCorrMulti(px.fn, g => g.outerSlope, structControls, hi).toFixed(3),
  })),
  extendedControls: {
    rPartRatio: +rExtRatio.toFixed(3),
    rPartDirect: +rExtDirect.toFixed(3),
    looExt: +looExt.toFixed(3),
  },
  ratioVsBoth: {
    looRatioOnly: +looRatioOnly.toFixed(3),
    looBothOnly: +looBothOnly.toFixed(3),
    looRatioFull: +looRatioFull.toFixed(3),
    looBothFull: +looBothFull.toFixed(3),
  },
  allGalaxies: { n: all.length, looR2: +looAll.toFixed(3) },
  permutation: { r: +realR.toFixed(3), p: +(nBetter / nPerm).toFixed(4) },
  genuineMarkers: genuineMarkers.length,
  geometricMarkers: geometricMarkers.length,
  verdict,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase97-geometric-test.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase97-geometric-test.json');

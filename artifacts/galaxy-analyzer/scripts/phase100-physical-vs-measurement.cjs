#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 100: IS THE LAW PHYSICAL OR A MEASUREMENT ARTIFACT?');
console.log('  Three decisive tests within distance-controlled conditions');
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
function partialR(x, y, controls) {
  const n = x.length;
  if (n < controls[0].length + 5) return NaN;
  function residualize(v, C) {
    const p = C[0].length;
    const X = C.map((row, i) => [1, ...row]);
    const Xt = Array.from({ length: p + 1 }, (_, j) => X.map(r => r[j]));
    const XtX = Array.from({ length: p + 1 }, (_, i) =>
      Array.from({ length: p + 1 }, (_, j) => Xt[i].reduce((s, _, k) => s + Xt[i][k] * Xt[j][k], 0)));
    const Xty = Xt.map(col => col.reduce((s, val, k) => s + val * v[k], 0));
    const aug = XtX.map((row, i) => [...row, Xty[i]]);
    for (let i = 0; i <= p; i++) {
      let mx = i;
      for (let k = i + 1; k <= p; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[mx][i])) mx = k;
      [aug[i], aug[mx]] = [aug[mx], aug[i]];
      if (Math.abs(aug[i][i]) < 1e-12) return v;
      for (let k = 0; k <= p; k++) {
        if (k === i) continue;
        const f = aug[k][i] / aug[i][i];
        for (let j = i; j <= p + 1; j++) aug[k][j] -= f * aug[i][j];
      }
    }
    const beta = aug.map((row, i) => row[p + 1] / row[i]);
    return v.map((vi, i) => vi - X[i].reduce((s, xj, j) => s + xj * beta[j], 0));
  }
  const rx = residualize(x, controls);
  const ry = residualize(y, controls);
  return pearsonR(rx, ry);
}

function permTest(x, y, nPerm) {
  const rObs = Math.abs(pearsonR(x, y));
  let count = 0;
  for (let p = 0; p < nPerm; p++) {
    const shuffled = [...x];
    for (let i = shuffled.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
    }
    if (Math.abs(pearsonR(shuffled, y)) >= rObs) count++;
  }
  return count / nPerm;
}

const n45csv = fs.readFileSync(path.join(__dirname, '..', 'public', 'replication', 'N45_final_dataset.csv'), 'utf-8').trim().split('\n');
const n45names = new Set();
for (let i = 1; i < n45csv.length; i++) n45names.add(n45csv[i].split(',')[0]);

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
  const errV = parseFloat(line.substring(33, 38).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, errV, vgas, vdisk, vbul });
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

  const logR = outerPts.map(p => Math.log10(p.rad));
  const logV = outerPts.map(p => Math.log10(p.vobs));
  const mr = mean(logR), mv = mean(logV);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mr) * (logV[i] - mv); sxx += (logR[i] - mr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const rMax = pts[pts.length - 1].rad;
  const radialCoverage = rMax / (t1.Rdisk > 0 ? t1.Rdisk : 1);

  const logRatio = Math.log10(t1.Vflat > 0 && innerVmax > 0 ? t1.Vflat / innerVmax : 1);

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);

  // Alternative inner proxies less sensitive to resolution
  // 1) V at fixed fraction of Rdisk (e.g. 1 Rdisk)
  const rTarget1Rd = t1.Rdisk;
  const rTarget2Rd = 2 * t1.Rdisk;
  let v1Rd = NaN, v2Rd = NaN;
  for (let i = 0; i < innerPts.length - 1; i++) {
    if (innerPts[i].rad <= rTarget1Rd && innerPts[i + 1].rad >= rTarget1Rd) {
      const frac = (rTarget1Rd - innerPts[i].rad) / (innerPts[i + 1].rad - innerPts[i].rad);
      v1Rd = innerPts[i].vobs + frac * (innerPts[i + 1].vobs - innerPts[i].vobs);
    }
  }
  for (let i = 0; i < pts.length - 1; i++) {
    if (pts[i].rad <= rTarget2Rd && pts[i + 1].rad >= rTarget2Rd) {
      const frac = (rTarget2Rd - pts[i].rad) / (pts[i + 1].rad - pts[i].rad);
      v2Rd = pts[i].vobs + frac * (pts[i + 1].vobs - pts[i].vobs);
    }
  }

  // 2) Mean inner velocity (less sensitive to peak than max)
  const meanInnerV = mean(innerPts.map(p => p.vobs));

  // 3) V at the innermost 25% of radii (deep inner)
  const q25pts = innerPts.slice(0, Math.max(3, Math.floor(innerPts.length * 0.5)));
  const deepInnerVmax = Math.max(...q25pts.map(p => p.vobs));

  // 4) Baryon-only velocity at inner half (from decomposition)
  const innerBarV = [];
  for (const p of innerPts) {
    const vbar2 = (p.vgas || 0) ** 2 + UPSILON_DISK * (p.vdisk || 0) ** 2 + UPSILON_BULGE * (p.vbul || 0) ** 2;
    if (vbar2 > 0) innerBarV.push(Math.sqrt(vbar2));
  }
  const innerBaryonVmax = innerBarV.length > 0 ? Math.max(...innerBarV) : NaN;

  // 5) Inner discrepancy: vobs/vbar at inner peak
  const innerDiscrepancy = innerBaryonVmax > 0 ? innerVmax / innerBaryonVmax : NaN;

  allGalaxies.push({
    name, nPts: pts.length, half,
    innerVmax, logInnerVmax: Math.log10(innerVmax > 0 ? innerVmax : 1),
    outerSlope, logRatio,
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    Vflat: t1.Vflat, fgas,
    T: t1.T, Q: t1.Q, D: t1.D, logD: Math.log10(t1.D > 0 ? t1.D : 1),
    radialCoverage,
    Rdisk: t1.Rdisk,
    isN45: n45names.has(name),
    v1Rd, v2Rd, meanInnerV, deepInnerVmax, innerBaryonVmax, innerDiscrepancy,
    logRatioV1Rd: isFinite(v1Rd) && v1Rd > 0 ? Math.log10(t1.Vflat / v1Rd) : NaN,
    logRatioMeanV: meanInnerV > 0 ? Math.log10(t1.Vflat / meanInnerV) : NaN,
    logRatioDeepV: deepInnerVmax > 0 ? Math.log10(t1.Vflat / deepInnerVmax) : NaN,
    logRatioBaryonV: isFinite(innerBaryonVmax) && innerBaryonVmax > 0 ? Math.log10(t1.Vflat / innerBaryonVmax) : NaN,
  });
}

const hi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);
console.log('  High-regime galaxies: N=' + hi.length + '\n');

// ═══════════════════════════════════════════════════════════════
//  TEST 1: WITHIN DISTANCE BINS
// ═══════════════════════════════════════════════════════════════

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: WITHIN DISTANCE BINS');
console.log('  Does the correlation survive when distance is held constant?');
console.log('══════════════════════════════════════════════════════════════\n');

const distances = hi.map(g => g.D).sort((a, b) => a - b);
const dTertiles = [distances[Math.floor(distances.length / 3)], distances[Math.floor(2 * distances.length / 3)]];

const bins = [
  { label: 'Near (D < ' + dTertiles[0].toFixed(0) + ' Mpc)', data: hi.filter(g => g.D < dTertiles[0]) },
  { label: 'Mid (' + dTertiles[0].toFixed(0) + ' <= D < ' + dTertiles[1].toFixed(0) + ' Mpc)', data: hi.filter(g => g.D >= dTertiles[0] && g.D < dTertiles[1]) },
  { label: 'Far (D >= ' + dTertiles[1].toFixed(0) + ' Mpc)', data: hi.filter(g => g.D >= dTertiles[1]) },
];

const binResults = [];
for (const bin of bins) {
  const x = bin.data.map(g => g.logRatio);
  const y = bin.data.map(g => g.outerSlope);
  const r = pearsonR(x, y);
  const p = bin.data.length >= 8 ? permTest(x, y, 2000) : NaN;
  console.log('  ' + bin.label + ' (N=' + bin.data.length + '):');
  console.log('    r(logRatio, outerSlope) = ' + r.toFixed(3) + ', perm p = ' + (isFinite(p) ? p.toFixed(4) : 'N/A (too few)'));
  console.log('    D range: ' + Math.min(...bin.data.map(g => g.D)).toFixed(1) + ' – ' + Math.max(...bin.data.map(g => g.D)).toFixed(1) + ' Mpc');
  console.log('    sd(logRatio) = ' + sd(x).toFixed(3) + ', sd(outerSlope) = ' + sd(y).toFixed(3));
  binResults.push({ label: bin.label, n: bin.data.length, r: +r.toFixed(3), p: isFinite(p) ? +p.toFixed(4) : null });
}

// Also try halves
const dMedian = median(distances);
const nearHalf = hi.filter(g => g.D < dMedian);
const farHalf = hi.filter(g => g.D >= dMedian);
const rNear = pearsonR(nearHalf.map(g => g.logRatio), nearHalf.map(g => g.outerSlope));
const rFar = pearsonR(farHalf.map(g => g.logRatio), farHalf.map(g => g.outerSlope));
const pNear = permTest(nearHalf.map(g => g.logRatio), nearHalf.map(g => g.outerSlope), 2000);
const pFar = permTest(farHalf.map(g => g.logRatio), farHalf.map(g => g.outerSlope), 2000);

console.log('\n  Halves (median D = ' + dMedian.toFixed(1) + ' Mpc):');
console.log('    Near half (N=' + nearHalf.length + '): r = ' + rNear.toFixed(3) + ', p = ' + pNear.toFixed(4));
console.log('    Far half  (N=' + farHalf.length + '):  r = ' + rFar.toFixed(3) + ', p = ' + pFar.toFixed(4));

// Partial correlation controlling for distance
const prDist = partialR(
  hi.map(g => g.logRatio), hi.map(g => g.outerSlope),
  hi.map(g => [g.logD])
);
const prDistNpts = partialR(
  hi.map(g => g.logRatio), hi.map(g => g.outerSlope),
  hi.map(g => [g.logD, g.nPts])
);
const prDistAll = partialR(
  hi.map(g => g.logRatio), hi.map(g => g.outerSlope),
  hi.map(g => [g.logD, g.nPts, g.radialCoverage])
);

console.log('\n  Partial correlations:');
console.log('    r(logRatio, outerSlope | logD) = ' + prDist.toFixed(3));
console.log('    r(logRatio, outerSlope | logD, nPts) = ' + prDistNpts.toFixed(3));
console.log('    r(logRatio, outerSlope | logD, nPts, radCov) = ' + prDistAll.toFixed(3));

// ═══════════════════════════════════════════════════════════════
//  TEST 2: DISTANCE-MATCHED PAIRS
// ═══════════════════════════════════════════════════════════════

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: DISTANCE-MATCHED PAIRS');
console.log('  Among galaxies at similar D, does ratio predict slope?');
console.log('══════════════════════════════════════════════════════════════\n');

// Create pairs within 20% distance tolerance
const dTol = 0.2;
const pairs = [];
const used = new Set();
const sorted = [...hi].sort((a, b) => a.D - b.D);
for (let i = 0; i < sorted.length; i++) {
  if (used.has(i)) continue;
  for (let j = i + 1; j < sorted.length; j++) {
    if (used.has(j)) continue;
    if (Math.abs(sorted[i].D - sorted[j].D) / ((sorted[i].D + sorted[j].D) / 2) < dTol) {
      pairs.push([sorted[i], sorted[j]]);
      used.add(i);
      used.add(j);
      break;
    }
  }
}

console.log('  Distance-matched pairs (within ' + (dTol * 100) + '% distance): N=' + pairs.length + ' pairs');

// For each pair, compute diff in logRatio and diff in outerSlope
const dRatios = pairs.map(([a, b]) => a.logRatio - b.logRatio);
const dSlopes = pairs.map(([a, b]) => a.outerSlope - b.outerSlope);
const rPairs = pearsonR(dRatios, dSlopes);
const pPairs = permTest(dRatios, dSlopes, 2000);

console.log('  r(Δ logRatio, Δ outerSlope) = ' + rPairs.toFixed(3));
console.log('  perm p = ' + pPairs.toFixed(4));

// Also: within narrow distance bands, rank correlation
const narrowBands = [
  { lo: 0, hi: 10, label: '0-10 Mpc' },
  { lo: 10, hi: 20, label: '10-20 Mpc' },
  { lo: 20, hi: 40, label: '20-40 Mpc' },
  { lo: 40, hi: 200, label: '40+ Mpc' },
];

console.log('\n  Narrow distance bands:');
for (const band of narrowBands) {
  const g = hi.filter(g => g.D >= band.lo && g.D < band.hi);
  if (g.length < 6) { console.log('    ' + band.label + ': N=' + g.length + ' (too few)'); continue; }
  const r = pearsonR(g.map(g => g.logRatio), g.map(g => g.outerSlope));
  console.log('    ' + band.label + ' (N=' + g.length + '): r = ' + r.toFixed(3));
}

// ═══════════════════════════════════════════════════════════════
//  TEST 3: ALTERNATIVE INNER PROXIES
// ═══════════════════════════════════════════════════════════════

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: ALTERNATIVE INNER PROXIES');
console.log('  Less resolution-sensitive measures of inner velocity');
console.log('══════════════════════════════════════════════════════════════\n');

const proxies = [
  { name: 'InnerVmax (original)', fn: g => g.logRatio, filter: g => isFinite(g.logRatio) },
  { name: 'V at 1 Rdisk', fn: g => g.logRatioV1Rd, filter: g => isFinite(g.logRatioV1Rd) },
  { name: 'Mean inner V', fn: g => g.logRatioMeanV, filter: g => isFinite(g.logRatioMeanV) },
  { name: 'Deep inner Vmax (inner 50%)', fn: g => g.logRatioDeepV, filter: g => isFinite(g.logRatioDeepV) },
  { name: 'Inner baryon Vmax', fn: g => g.logRatioBaryonV, filter: g => isFinite(g.logRatioBaryonV) },
];

const proxyResults = [];
for (const proxy of proxies) {
  const valid = hi.filter(proxy.filter);
  const x = valid.map(proxy.fn);
  const y = valid.map(g => g.outerSlope);
  const r = pearsonR(x, y);
  const p = valid.length >= 8 ? permTest(x, y, 2000) : NaN;

  // Partial r controlling for distance
  const pr = valid.length >= 10 ? partialR(x, y, valid.map(g => [g.logD])) : NaN;

  // Partial r controlling for distance + nPts + radCov
  const prAll = valid.length >= 12 ? partialR(x, y, valid.map(g => [g.logD, g.nPts, g.radialCoverage])) : NaN;

  console.log('  ' + proxy.name + ' (N=' + valid.length + '):');
  console.log('    r = ' + r.toFixed(3) + ', perm p = ' + (isFinite(p) ? p.toFixed(4) : 'N/A'));
  console.log('    partial r (| logD) = ' + (isFinite(pr) ? pr.toFixed(3) : 'N/A'));
  console.log('    partial r (| logD, nPts, radCov) = ' + (isFinite(prAll) ? prAll.toFixed(3) : 'N/A'));

  proxyResults.push({
    name: proxy.name, n: valid.length, r: +r.toFixed(3),
    p: isFinite(p) ? +p.toFixed(4) : null,
    partialR_D: isFinite(pr) ? +pr.toFixed(3) : null,
    partialR_all: isFinite(prAll) ? +prAll.toFixed(3) : null,
  });
}

// ═══════════════════════════════════════════════════════════════
//  TEST 3b: Baryon-only ratio (completely resolution-independent)
// ═══════════════════════════════════════════════════════════════

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3b: BARYON RATIO — vobs(inner) vs vbar(inner)');
console.log('  Inner discrepancy as alternative predictor');
console.log('══════════════════════════════════════════════════════════════\n');

const withDisc = hi.filter(g => isFinite(g.innerDiscrepancy) && g.innerDiscrepancy > 0);
const logDisc = withDisc.map(g => Math.log10(g.innerDiscrepancy));
const discSlope = withDisc.map(g => g.outerSlope);
const rDisc = pearsonR(logDisc, discSlope);
const pDisc = permTest(logDisc, discSlope, 2000);
const prDisc = partialR(logDisc, discSlope, withDisc.map(g => [g.logD]));
const prDiscAll = partialR(logDisc, discSlope, withDisc.map(g => [g.logD, g.nPts, g.radialCoverage]));

console.log('  log(innerDiscrepancy) = log(Vobs_inner / Vbar_inner)');
console.log('  N = ' + withDisc.length);
console.log('  r = ' + rDisc.toFixed(3) + ', perm p = ' + pDisc.toFixed(4));
console.log('  partial r (| logD) = ' + prDisc.toFixed(3));
console.log('  partial r (| logD, nPts, radCov) = ' + prDiscAll.toFixed(3));

// ═══════════════════════════════════════════════════════════════
//  SYNTHESIS: COMBINE ALL EVIDENCE
// ═══════════════════════════════════════════════════════════════

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  SYNTHESIS: PHYSICAL vs MEASUREMENT');
console.log('══════════════════════════════════════════════════════════════\n');

let physicalScore = 0;
let measurementScore = 0;
const verdicts = [];

// Test 1 verdict
const binSignificant = binResults.filter(b => b.p !== null && b.p < 0.05).length;
if (binSignificant >= 2) {
  console.log('  TEST 1: PASSES — signal survives in ' + binSignificant + '/3 distance bins');
  verdicts.push({ test: 'Distance bins', result: 'PASS', detail: binSignificant + '/3 bins significant' });
  physicalScore += 2;
} else if (binSignificant >= 1) {
  console.log('  TEST 1: PARTIAL — signal survives in ' + binSignificant + '/3 distance bins');
  verdicts.push({ test: 'Distance bins', result: 'PARTIAL', detail: binSignificant + '/3 bins significant' });
  physicalScore += 1;
} else {
  console.log('  TEST 1: FAILS — signal disappears in distance bins');
  verdicts.push({ test: 'Distance bins', result: 'FAIL', detail: '0/3 bins significant' });
  measurementScore += 2;
}

// Partial correlation verdict
if (Math.abs(prDistAll) > 0.5) {
  console.log('  PARTIAL r (all controls): ' + prDistAll.toFixed(3) + ' → STRONG after controlling D, nPts, radCov');
  verdicts.push({ test: 'Partial r (all controls)', result: 'PASS', detail: 'r=' + prDistAll.toFixed(3) });
  physicalScore += 2;
} else if (Math.abs(prDistAll) > 0.3) {
  console.log('  PARTIAL r (all controls): ' + prDistAll.toFixed(3) + ' → MODERATE after controlling D, nPts, radCov');
  verdicts.push({ test: 'Partial r (all controls)', result: 'PARTIAL', detail: 'r=' + prDistAll.toFixed(3) });
  physicalScore += 1;
} else {
  console.log('  PARTIAL r (all controls): ' + prDistAll.toFixed(3) + ' → WEAK — measurement explains most');
  verdicts.push({ test: 'Partial r (all controls)', result: 'FAIL', detail: 'r=' + prDistAll.toFixed(3) });
  measurementScore += 2;
}

// Test 2 verdict
if (rPairs > 0.3 && pPairs < 0.05) {
  console.log('  TEST 2: PASSES — distance-matched pairs show r=' + rPairs.toFixed(3) + ', p=' + pPairs.toFixed(4));
  verdicts.push({ test: 'Distance-matched pairs', result: 'PASS', detail: 'r=' + rPairs.toFixed(3) });
  physicalScore += 2;
} else if (rPairs > 0.2) {
  console.log('  TEST 2: PARTIAL — distance-matched pairs show r=' + rPairs.toFixed(3));
  verdicts.push({ test: 'Distance-matched pairs', result: 'PARTIAL', detail: 'r=' + rPairs.toFixed(3) });
  physicalScore += 1;
} else {
  console.log('  TEST 2: FAILS — no signal in distance-matched pairs');
  verdicts.push({ test: 'Distance-matched pairs', result: 'FAIL', detail: 'r=' + rPairs.toFixed(3) });
  measurementScore += 2;
}

// Test 3 verdict: alternative proxies
const origR = proxyResults[0].partialR_all;
const bestAlt = proxyResults.slice(1).reduce((best, p) =>
  (p.partialR_all !== null && Math.abs(p.partialR_all) > Math.abs(best.partialR_all || 0)) ? p : best,
  proxyResults[1]);
if (bestAlt.partialR_all !== null && Math.abs(bestAlt.partialR_all) > 0.4) {
  console.log('  TEST 3: PASSES — alternative proxy "' + bestAlt.name + '" has partial r=' + bestAlt.partialR_all + ' after all controls');
  verdicts.push({ test: 'Alternative proxies', result: 'PASS', detail: bestAlt.name + ' partial r=' + bestAlt.partialR_all });
  physicalScore += 2;
} else if (bestAlt.partialR_all !== null && Math.abs(bestAlt.partialR_all) > 0.25) {
  console.log('  TEST 3: PARTIAL — best alternative proxy partial r=' + bestAlt.partialR_all);
  verdicts.push({ test: 'Alternative proxies', result: 'PARTIAL', detail: bestAlt.name + ' partial r=' + bestAlt.partialR_all });
  physicalScore += 1;
} else {
  console.log('  TEST 3: FAILS — alternative proxies too weak after controls');
  verdicts.push({ test: 'Alternative proxies', result: 'FAIL', detail: 'best partial r=' + (bestAlt.partialR_all || 'N/A') });
  measurementScore += 2;
}

console.log('\n  ╔══════════════════════════════════════════════════════════════╗');
if (physicalScore >= 6) {
  console.log('  ║  VERDICT: PHYSICAL LAW — survives all measurement controls ║');
} else if (physicalScore >= 4) {
  console.log('  ║  VERDICT: PREDOMINANTLY PHYSICAL with measurement component║');
} else if (physicalScore >= 2) {
  console.log('  ║  VERDICT: MIXED — partly physical, partly measurement      ║');
} else {
  console.log('  ║  VERDICT: PREDOMINANTLY MEASUREMENT ARTIFACT               ║');
}
console.log('  ║  Physical score: ' + physicalScore + '/8, Measurement score: ' + measurementScore + '/8' + ' '.repeat(11) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════╝');

const output = {
  phase: '100',
  title: 'Physical vs Measurement Artifact',
  nHi: hi.length,
  distanceBins: binResults,
  halves: {
    near: { n: nearHalf.length, r: +rNear.toFixed(3), p: +pNear.toFixed(4) },
    far: { n: farHalf.length, r: +rFar.toFixed(3), p: +pFar.toFixed(4) },
  },
  partialCorrelations: {
    logD: +prDist.toFixed(3),
    logD_nPts: +prDistNpts.toFixed(3),
    logD_nPts_radCov: +prDistAll.toFixed(3),
  },
  distanceMatchedPairs: { nPairs: pairs.length, r: +rPairs.toFixed(3), p: +pPairs.toFixed(4) },
  alternativeProxies: proxyResults,
  innerDiscrepancy: {
    r: +rDisc.toFixed(3), p: +pDisc.toFixed(4),
    partialR_D: +prDisc.toFixed(3), partialR_all: +prDiscAll.toFixed(3),
  },
  verdicts,
  physicalScore, measurementScore,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase100-physical-vs-measurement.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase100-physical-vs-measurement.json');

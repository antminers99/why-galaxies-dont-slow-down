#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 83: COMMON SUPPORT TEST');
console.log('  Match external galaxies to training regime, re-test frozen M3');
console.log('  Question: Does the law work on similar-type external galaxies?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

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
  const logA0 = (lo + hi) / 2, a0 = Math.pow(10, logA0);
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.logGbar);
    const pred = mcgaughRAR(gb, a0);
    ss += (p.logGobs - Math.log10(pred > 0 ? pred : 1e-20)) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length) };
}

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function median(a) { const s = [...a].sort((x, y) => x - y); return s.length % 2 ? s[(s.length - 1) / 2] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2; }

const csvPath = path.join(__dirname, '..', 'public', 'replication', 'N45_final_dataset.csv');
const csvLines = fs.readFileSync(csvPath, 'utf-8').trim().split('\n');
const headers = csvLines[0].split(',');
const trainCSV = csvLines.slice(1).map(l => {
  const vals = l.split(',');
  const obj = {};
  headers.forEach((h, i) => obj[h] = vals[i]);
  return obj;
});
const trainNames = new Set(trainCSV.map(d => d.galaxy_name));
const trainLogA0 = trainCSV.map(d => parseFloat(d.logA0));
const trainMeanA0 = mean(trainLogA0);

const M3_BETA = [5.182, -0.198, -0.155, 0.459];
console.log('  FROZEN M3: intercept=' + M3_BETA[0] + ', logMHI=' + M3_BETA[1] +
  ', logMhost=' + M3_BETA[2] + ', logMeanRun=' + M3_BETA[3]);

const M2P_BETA = [4.191, -0.136, 0.313];
console.log('  FROZEN M2\' (no Mhost): intercept=' + M2P_BETA[0] + ', logMHI=' + M2P_BETA[1] +
  ', logMeanRun=' + M2P_BETA[2]);
console.log();

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const D = parseFloat(line.substring(15, 21).trim());
  const fD = parseInt(line.substring(28, 29).trim());
  const inc = parseFloat(line.substring(30, 34).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { T, D, fD, inc, L36, Rdisk, MHI, Vflat, Q };
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

const sparcAll = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf-8'));
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
    pts.push({ r: pt.rad, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar) });
  }
  if (pts.length < 5) continue;

  const fit = fitA0(pts);
  if (!isFinite(fit.logA0) || fit.logA0 < 0.5 || fit.logA0 > 5.5) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
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
  const logMeanRun = Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1);

  const logMHI = Math.log10(t1.MHI);
  const logVflat = Math.log10(t1.Vflat > 0 ? t1.Vflat : 50);
  const logMhost = 4 * logVflat - 1.5;

  const ev = envLookup[name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI, logMeanRun, logMhost,
    Vflat: t1.Vflat, Q: t1.Q, D: t1.D, fD: t1.fD, T: t1.T,
    L36: t1.L36, Rdisk: t1.Rdisk, inc: t1.inc,
    envCode,
    inTraining: trainNames.has(name),
  });
}

const train = allGalaxies.filter(g => g.inTraining);
const external = allGalaxies.filter(g => !g.inTraining);

console.log('  Total galaxies computed: ' + allGalaxies.length);
console.log('  Training: ' + train.length + ', External: ' + external.length);
console.log();

console.log('======================================================================');
console.log('  STEP 1: Define training regime (common support bounds)');
console.log('======================================================================\n');

const matchVars = [
  { name: 'Vflat', key: 'Vflat' },
  { name: 'logMHI', key: 'logMHI' },
  { name: 'logMeanRun', key: 'logMeanRun' },
  { name: 'nPts', key: 'nPts' },
  { name: 'T (morphology)', key: 'T' },
];

const bounds = {};
for (const v of matchVars) {
  const tVals = train.map(g => g[v.key]).filter(x => isFinite(x));
  const lo = Math.min(...tVals);
  const hi = Math.max(...tVals);
  const p10 = [...tVals].sort((a, b) => a - b)[Math.floor(tVals.length * 0.10)];
  const p90 = [...tVals].sort((a, b) => a - b)[Math.floor(tVals.length * 0.90)];
  bounds[v.key] = { lo, hi, p10, p90, tMean: mean(tVals), tSD: sd(tVals) };
  console.log('  ' + v.name.padEnd(20) + ': range [' + lo.toFixed(2) + ', ' + hi.toFixed(2) + ']' +
    '  P10-P90 [' + p10.toFixed(2) + ', ' + p90.toFixed(2) + ']');
}
console.log();

console.log('======================================================================');
console.log('  STEP 2: Apply matching criteria (progressively strict)');
console.log('======================================================================\n');

function matchLevel(g, level) {
  if (level === 'loose') {
    return g.Vflat >= bounds.Vflat.lo &&
           g.logMHI >= bounds.logMHI.lo && g.logMHI <= bounds.logMHI.hi &&
           g.nPts >= 7;
  }
  if (level === 'moderate') {
    return g.Vflat >= bounds.Vflat.p10 && g.Vflat <= bounds.Vflat.p90 &&
           g.logMHI >= bounds.logMHI.p10 && g.logMHI <= bounds.logMHI.p90 &&
           g.logMeanRun >= bounds.logMeanRun.p10 &&
           g.nPts >= 10 &&
           g.T <= 8;
  }
  if (level === 'strict') {
    return g.Vflat >= bounds.Vflat.p10 && g.Vflat <= bounds.Vflat.p90 &&
           g.logMHI >= bounds.logMHI.p10 && g.logMHI <= bounds.logMHI.p90 &&
           g.logMeanRun >= bounds.logMeanRun.p10 && g.logMeanRun <= bounds.logMeanRun.p90 &&
           g.nPts >= 15 &&
           g.T <= 7 &&
           g.Q <= 2;
  }
}

const levels = ['loose', 'moderate', 'strict'];
const levelResults = {};

for (const level of levels) {
  const matched = external.filter(g => matchLevel(g, level));

  console.log('  ── ' + level.toUpperCase() + ' matching: N = ' + matched.length + ' ──');
  if (matched.length < 3) {
    console.log('    Too few galaxies for meaningful test.\n');
    levelResults[level] = { n: matched.length, verdict: 'INSUFFICIENT DATA' };
    continue;
  }

  console.log('    Matched galaxies: ' + matched.map(g => g.name).join(', '));
  console.log();

  const matchedVals = {};
  for (const v of matchVars) {
    const eVals = matched.map(g => g[v.key]).filter(x => isFinite(x));
    const tVals = train.map(g => g[v.key]).filter(x => isFinite(x));
    matchedVals[v.key] = {
      tMean: mean(tVals).toFixed(2), eMean: mean(eVals).toFixed(2),
      tMed: median(tVals).toFixed(2), eMed: median(eVals).toFixed(2),
    };
    console.log('    ' + v.name.padEnd(20) +
      ': train=' + mean(tVals).toFixed(2) + ' ± ' + sd(tVals).toFixed(2) +
      ', matched=' + mean(eVals).toFixed(2) + ' ± ' + sd(eVals).toFixed(2));
  }
  console.log();

  const extMeanA0 = mean(matched.map(g => g.logA0));

  let m0_ss = 0, m3_ss = 0, m2p_ss = 0;
  let m3_wins = 0, m2p_wins = 0;
  const perGalaxy = [];

  for (const g of matched) {
    const err_m0 = (g.logA0 - trainMeanA0) ** 2;
    const pred_m3 = M3_BETA[0] + M3_BETA[1] * g.logMHI + M3_BETA[2] * g.logMhost + M3_BETA[3] * g.logMeanRun;
    const err_m3 = (g.logA0 - pred_m3) ** 2;
    const pred_m2p = M2P_BETA[0] + M2P_BETA[1] * g.logMHI + M2P_BETA[2] * g.logMeanRun;
    const err_m2p = (g.logA0 - pred_m2p) ** 2;

    m0_ss += err_m0;
    m3_ss += err_m3;
    m2p_ss += err_m2p;
    if (err_m3 < err_m0) m3_wins++;
    if (err_m2p < err_m0) m2p_wins++;

    perGalaxy.push({
      name: g.name, logA0: +g.logA0.toFixed(3),
      pred_m3: +pred_m3.toFixed(3), err_m3: +Math.sqrt(err_m3).toFixed(3),
      pred_m2p: +pred_m2p.toFixed(3), err_m2p: +Math.sqrt(err_m2p).toFixed(3),
      err_m0: +Math.sqrt(err_m0).toFixed(3),
      m3_wins: err_m3 < err_m0,
    });
  }

  const n = matched.length;
  const rms_m0 = Math.sqrt(m0_ss / n);
  const rms_m3 = Math.sqrt(m3_ss / n);
  const rms_m2p = Math.sqrt(m2p_ss / n);
  const gap_m3 = ((rms_m0 - rms_m3) / rms_m0 * 100);
  const gap_m2p = ((rms_m0 - rms_m2p) / rms_m0 * 100);
  const winRate_m3 = m3_wins / n * 100;
  const winRate_m2p = m2p_wins / n * 100;

  console.log('    M0 (baseline)  RMS = ' + rms_m0.toFixed(3));
  console.log('    M3 (frozen)    RMS = ' + rms_m3.toFixed(3) + '  Gap% = ' + gap_m3.toFixed(1) + '%  Win rate = ' + winRate_m3.toFixed(1) + '%');
  console.log('    M2\' (no Mhost) RMS = ' + rms_m2p.toFixed(3) + '  Gap% = ' + gap_m2p.toFixed(1) + '%  Win rate = ' + winRate_m2p.toFixed(1) + '%');

  const m3_pass = rms_m3 < rms_m0 && winRate_m3 > 50;
  const m2p_pass = rms_m2p < rms_m0 && winRate_m2p > 50;
  const verdict = m3_pass ? 'M3 PASSES' : m2p_pass ? 'M2\' PASSES (M3 fails)' : 'BOTH FAIL';
  console.log('    *** VERDICT: ' + verdict + ' ***');
  console.log();

  const partialCorr = (matched_arr, key) => {
    const x = matched_arr.map(g => g[key]);
    const y = matched_arr.map(g => g.logA0);
    const mx = mean(x), my = mean(y);
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < x.length; i++) {
      sxy += (x[i] - mx) * (y[i] - my);
      sxx += (x[i] - mx) ** 2;
      syy += (y[i] - my) ** 2;
    }
    return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
  };

  const r_mhi = partialCorr(matched, 'logMHI');
  const r_run = partialCorr(matched, 'logMeanRun');
  console.log('    Correlation with logA0:');
  console.log('      logMHI:     r = ' + r_mhi.toFixed(3) + ' (training expects NEGATIVE)');
  console.log('      logMeanRun: r = ' + r_run.toFixed(3) + ' (training expects POSITIVE)');
  const signs_correct = r_mhi < 0 && r_run > 0;
  console.log('      Signs correct? ' + (signs_correct ? 'YES' : 'NO — REVERSED'));
  console.log();

  console.log('    Per-galaxy results:');
  console.log('    ' + 'Galaxy'.padEnd(14) + 'logA0'.padStart(8) + 'M3_pred'.padStart(8) + 'M3_err'.padStart(8) + 'M0_err'.padStart(8) + 'M3 wins?'.padStart(10));
  for (const pg of perGalaxy) {
    console.log('    ' + pg.name.padEnd(14) + pg.logA0.toFixed(3).padStart(8) + pg.pred_m3.toFixed(3).padStart(8) +
      pg.err_m3.toFixed(3).padStart(8) + pg.err_m0.toFixed(3).padStart(8) + (pg.m3_wins ? 'YES' : 'no').padStart(10));
  }
  console.log();

  levelResults[level] = {
    n, rms_m0: +rms_m0.toFixed(3), rms_m3: +rms_m3.toFixed(3), rms_m2p: +rms_m2p.toFixed(3),
    gap_m3: +gap_m3.toFixed(1), gap_m2p: +gap_m2p.toFixed(1),
    winRate_m3: +winRate_m3.toFixed(1), winRate_m2p: +winRate_m2p.toFixed(1),
    r_logMHI: +r_mhi.toFixed(3), r_logMeanRun: +r_run.toFixed(3),
    signsCorrect: signs_correct,
    verdict,
    matchedGalaxies: matched.map(g => g.name),
    perGalaxy,
    matchedStats: matchedVals,
  };
}

console.log('======================================================================');
console.log('  STEP 3: Propensity-style matched pairs');
console.log('  Match each training galaxy to nearest external by Mahalanobis');
console.log('======================================================================\n');

const matchKeys = ['Vflat', 'logMHI', 'logMeanRun', 'nPts'];
const trainMeans = {}, trainSDs = {};
for (const k of matchKeys) {
  const vals = train.map(g => g[k]);
  trainMeans[k] = mean(vals);
  trainSDs[k] = sd(vals);
}

function mahal(g1, g2) {
  let d2 = 0;
  for (const k of matchKeys) {
    const z1 = (g1[k] - trainMeans[k]) / (trainSDs[k] || 1);
    const z2 = (g2[k] - trainMeans[k]) / (trainSDs[k] || 1);
    d2 += (z1 - z2) ** 2;
  }
  return Math.sqrt(d2);
}

const pairs = [];
const usedExt = new Set();

for (const tg of train) {
  let bestDist = Infinity, bestMatch = null;
  for (const eg of external) {
    if (usedExt.has(eg.name)) continue;
    const d = mahal(tg, eg);
    if (d < bestDist) { bestDist = d; bestMatch = eg; }
  }
  if (bestMatch && bestDist < 3.0) {
    pairs.push({ train: tg, ext: bestMatch, distance: bestDist });
    usedExt.add(bestMatch.name);
  }
}

console.log('  Matched pairs found: ' + pairs.length + ' / ' + train.length);
console.log();

if (pairs.length >= 5) {
  console.log('  Pair details:');
  console.log('  ' + 'Train'.padEnd(14) + 'External'.padEnd(14) + 'dist'.padStart(6) +
    'tVflat'.padStart(8) + 'eVflat'.padStart(8) + 'tMHI'.padStart(8) + 'eMHI'.padStart(8));
  for (const p of pairs) {
    console.log('  ' + p.train.name.padEnd(14) + p.ext.name.padEnd(14) +
      p.distance.toFixed(2).padStart(6) +
      p.train.Vflat.toFixed(0).padStart(8) + p.ext.Vflat.toFixed(0).padStart(8) +
      p.train.logMHI.toFixed(2).padStart(8) + p.ext.logMHI.toFixed(2).padStart(8));
  }
  console.log();

  const matchedExt = pairs.map(p => p.ext);
  let m0_ss = 0, m3_ss = 0, m2p_ss = 0, m3_wins = 0;
  for (const g of matchedExt) {
    const err_m0 = (g.logA0 - trainMeanA0) ** 2;
    const pred_m3 = M3_BETA[0] + M3_BETA[1] * g.logMHI + M3_BETA[2] * g.logMhost + M3_BETA[3] * g.logMeanRun;
    const err_m3 = (g.logA0 - pred_m3) ** 2;
    const pred_m2p = M2P_BETA[0] + M2P_BETA[1] * g.logMHI + M2P_BETA[2] * g.logMeanRun;
    const err_m2p = (g.logA0 - pred_m2p) ** 2;
    m0_ss += err_m0; m3_ss += err_m3; m2p_ss += err_m2p;
    if (err_m3 < err_m0) m3_wins++;
  }
  const n = matchedExt.length;
  console.log('  Matched-pair M0 RMS = ' + Math.sqrt(m0_ss / n).toFixed(3));
  console.log('  Matched-pair M3 RMS = ' + Math.sqrt(m3_ss / n).toFixed(3) +
    '  Win rate = ' + (m3_wins / n * 100).toFixed(1) + '%');
  console.log('  Matched-pair M2\' RMS = ' + Math.sqrt(m2p_ss / n).toFixed(3));

  const r_mhi = (() => {
    const x = matchedExt.map(g => g.logMHI), y = matchedExt.map(g => g.logA0);
    const mx = mean(x), my = mean(y);
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < x.length; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
    return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
  })();
  const r_run = (() => {
    const x = matchedExt.map(g => g.logMeanRun), y = matchedExt.map(g => g.logA0);
    const mx = mean(x), my = mean(y);
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < x.length; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
    return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
  })();
  console.log('  Correlations: logMHI r=' + r_mhi.toFixed(3) + ', logMeanRun r=' + r_run.toFixed(3));
  console.log('  Signs correct? ' + (r_mhi < 0 && r_run > 0 ? 'YES' : 'NO'));
  console.log();
}

console.log('======================================================================');
console.log('  FINAL SUMMARY');
console.log('======================================================================\n');

for (const level of levels) {
  const r = levelResults[level];
  if (!r || r.verdict === 'INSUFFICIENT DATA') {
    console.log('  ' + level.toUpperCase().padEnd(12) + ': N=' + (r ? r.n : 0) + ' — INSUFFICIENT DATA');
  } else {
    console.log('  ' + level.toUpperCase().padEnd(12) + ': N=' + r.n +
      '  M3 RMS=' + r.rms_m3.toFixed(3) + ' vs M0=' + r.rms_m0.toFixed(3) +
      '  Win=' + r.winRate_m3.toFixed(0) + '%' +
      '  Signs=' + (r.signsCorrect ? 'OK' : 'REVERSED') +
      '  → ' + r.verdict);
  }
}
if (pairs.length >= 5) {
  const matchedExt = pairs.map(p => p.ext);
  const n = matchedExt.length;
  let m0_ss = 0, m3_ss = 0, m3_wins = 0;
  for (const g of matchedExt) {
    m0_ss += (g.logA0 - trainMeanA0) ** 2;
    const pred = M3_BETA[0] + M3_BETA[1] * g.logMHI + M3_BETA[2] * g.logMhost + M3_BETA[3] * g.logMeanRun;
    m3_ss += (g.logA0 - pred) ** 2;
    if ((g.logA0 - pred) ** 2 < (g.logA0 - trainMeanA0) ** 2) m3_wins++;
  }
  console.log('  PAIRED'.padEnd(12) + ': N=' + n +
    '  M3 RMS=' + Math.sqrt(m3_ss / n).toFixed(3) + ' vs M0=' + Math.sqrt(m0_ss / n).toFixed(3) +
    '  Win=' + (m3_wins / n * 100).toFixed(0) + '%');
}
console.log();

const anyPass = Object.values(levelResults).some(r => r.verdict && r.verdict.includes('PASS'));
if (anyPass) {
  console.log('  ╔════════════════════════════════════════════════════════════════╗');
  console.log('  ║  LAW WORKS WITHIN REGIME — sample-specific but real pattern  ║');
  console.log('  ╚════════════════════════════════════════════════════════════════╝');
} else {
  console.log('  ╔════════════════════════════════════════════════════════════════╗');
  console.log('  ║  LAW FAILS EVEN ON MATCHED GALAXIES — deeper problem         ║');
  console.log('  ╚════════════════════════════════════════════════════════════════╝');
}
console.log();

const output = {
  phase: '83',
  title: 'Common Support Test — Regime-Matched External Validation',
  nTraining: train.length,
  nExternal: external.length,
  matchingBounds: bounds,
  levelResults,
  nMatchedPairs: pairs.length,
  matchedPairGalaxies: pairs.map(p => ({ train: p.train.name, ext: p.ext.name, distance: +p.distance.toFixed(3) })),
  anyLevelPasses: anyPass,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase83-common-support.json'), JSON.stringify(output, null, 2));
console.log('  Saved: public/phase83-common-support.json');

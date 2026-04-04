#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 81b: DIAGNOSTIC EXTERNAL TESTS');
console.log('  Test 1: M2\' (logMHI + logMeanRun only — drop Mhost)');
console.log('  Test 2: M3_fair (Vflat-estimated Mhost on both train + ext)');
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
    ss += (p.logGobs - Math.log10(mcgaughRAR(gb, a0))) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length) };
}

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }

function olsFit(Y, X) {
  const n = Y.length, p = X[0].length;
  const Xa = X.map(row => [1, ...row]);
  const pp = p + 1;
  const XtX = Array.from({ length: pp }, () => new Array(pp).fill(0));
  const XtY = new Array(pp).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < pp; j++) {
      XtY[j] += Xa[i][j] * Y[i];
      for (let k = 0; k < pp; k++) XtX[j][k] += Xa[i][j] * Xa[i][k];
    }
  }
  const aug = XtX.map((row, i) => [...row, ...Array.from({ length: pp }, (_, j) => i === j ? 1 : 0)]);
  for (let col = 0; col < pp; col++) {
    let maxRow = col;
    for (let row = col + 1; row < pp; row++)
      if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row;
    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
    const pivot = aug[col][col];
    if (Math.abs(pivot) < 1e-12) return null;
    for (let j = 0; j < 2 * pp; j++) aug[col][j] /= pivot;
    for (let row = 0; row < pp; row++) {
      if (row === col) continue;
      const factor = aug[row][col];
      for (let j = 0; j < 2 * pp; j++) aug[row][j] -= factor * aug[col][j];
    }
  }
  const inv = aug.map(row => row.slice(pp));
  return inv.map(row => row.reduce((s, v, j) => s + v * XtY[j], 0));
}

function spearman(obs, pred) {
  const rank = arr => {
    const sorted = [...arr].map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
    const r = new Array(arr.length);
    sorted.forEach((s, rank) => r[s.i] = rank + 1);
    return r;
  };
  const ro = rank(obs), rp = rank(pred);
  const mo = mean(ro), mp = mean(rp);
  const num = ro.reduce((s, v, i) => s + (v - mo) * (rp[i] - mp), 0);
  const d1 = Math.sqrt(ro.reduce((s, v) => s + (v - mo) ** 2, 0));
  const d2 = Math.sqrt(rp.reduce((s, v) => s + (v - mp) ** 2, 0));
  return num / (d1 * d2 + 1e-12);
}

function pearson(obs, pred) {
  const mo = mean(obs), mp = mean(pred);
  const num = obs.reduce((s, v, i) => s + (v - mo) * (pred[i] - mp), 0);
  const d1 = Math.sqrt(obs.reduce((s, v) => s + (v - mo) ** 2, 0));
  const d2 = Math.sqrt(pred.reduce((s, v) => s + (v - mp) ** 2, 0));
  return num / (d1 * d2 + 1e-12);
}

const csvPath = path.join(__dirname, '..', 'public', 'replication', 'N45_final_dataset.csv');
const csvLines = fs.readFileSync(csvPath, 'utf-8').trim().split('\n');
const headers = csvLines[0].split(',');
const trainData = csvLines.slice(1).map(l => {
  const vals = l.split(',');
  const obj = {};
  headers.forEach((h, i) => obj[h] = vals[i]);
  return obj;
});
const trainNames = new Set(trainData.map(d => d.galaxy_name));
const trainY = trainData.map(d => parseFloat(d.logA0));
const trainMeanA0 = mean(trainY);

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const D = parseFloat(line.substring(15, 21).trim());
  const fD = parseInt(line.substring(28, 29).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { D, fD, L36, Rdisk, MHI, Vflat, Q };
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
  const logMhost_vflat = t1.Vflat > 0 ? 4.0 * Math.log10(t1.Vflat) + 2.4 : null;

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI, logMeanRun, logMhost_vflat,
    Vflat: t1.Vflat, Q: t1.Q, fD: t1.fD,
    inTraining: trainNames.has(name),
  });
}

const trainGals = allGalaxies.filter(g => g.inTraining);
const extGals = allGalaxies.filter(g => !g.inTraining && g.fD >= 2 && g.logA0 >= 2.0 && g.nPts >= 5);

console.log('  Total computed: ' + allGalaxies.length);
console.log('  Training: ' + trainGals.length + ', External (filtered): ' + extGals.length);
console.log();

console.log('======================================================================');
console.log('  TEST 1: M2\' = logMHI + logMeanRun (no Mhost at all)');
console.log('======================================================================\n');

const trainX_M2p = trainData.map(d => [parseFloat(d.logMHI), parseFloat(d.logMeanRun)]);
const betaM2p = olsFit(trainY, trainX_M2p);
console.log('  Frozen M2\' (fit on training N=45):');
console.log('    intercept  = ' + betaM2p[0].toFixed(4));
console.log('    logMHI     = ' + betaM2p[1].toFixed(4));
console.log('    logMeanRun = ' + betaM2p[2].toFixed(4));

let looSS = 0;
for (let i = 0; i < trainY.length; i++) {
  const Yr = trainY.filter((_, j) => j !== i);
  const Xr = trainX_M2p.filter((_, j) => j !== i);
  const b = olsFit(Yr, Xr);
  const pred = b[0] + b[1] * trainX_M2p[i][0] + b[2] * trainX_M2p[i][1];
  looSS += (trainY[i] - pred) ** 2;
}
const looRMS_M2p = Math.sqrt(looSS / trainY.length);
const sdTrain = sd(trainY);
const looGap_M2p = 100 * (1 - looRMS_M2p ** 2 / sdTrain ** 2);
console.log('    LOO RMS = ' + looRMS_M2p.toFixed(4) + ', LOO gap% = ' + looGap_M2p.toFixed(1) + '%');
console.log();

const extObs = extGals.map(g => g.logA0);
const extPredM0 = extGals.map(() => trainMeanA0);
const extPredM2p = extGals.map(g => betaM2p[0] + betaM2p[1] * g.logMHI + betaM2p[2] * g.logMeanRun);
const sdExt = sd(extObs);
const rmsM0 = Math.sqrt(extPredM0.reduce((s, p, i) => s + (extObs[i] - p) ** 2, 0) / extGals.length);
const rmsM2p = Math.sqrt(extPredM2p.reduce((s, p, i) => s + (extObs[i] - p) ** 2, 0) / extGals.length);
const biasM2p = extPredM2p.reduce((s, p, i) => s + (p - extObs[i]), 0) / extGals.length;
const gapM2p = 100 * (1 - rmsM2p ** 2 / sdExt ** 2);
const rsM2p = spearman(extObs, extPredM2p);
const rpM2p = pearson(extObs, extPredM2p);

console.log('  External N = ' + extGals.length + ', SD(logA0) = ' + sdExt.toFixed(4));
console.log('  ' + 'Model'.padEnd(10) + 'RMS'.padStart(8) + 'Bias'.padStart(8) + 'Gap%'.padStart(8) + 'rs'.padStart(8) + 'rp'.padStart(8));
console.log('  ' + '-'.repeat(50));
console.log('  ' + 'M0'.padEnd(10) + rmsM0.toFixed(4).padStart(8) + '—'.padStart(8) + '—'.padStart(8) + '—'.padStart(8) + '—'.padStart(8));
console.log('  ' + 'M2\''.padEnd(10) + rmsM2p.toFixed(4).padStart(8) + biasM2p.toFixed(4).padStart(8) + (gapM2p.toFixed(1) + '%').padStart(8) + rsM2p.toFixed(3).padStart(8) + rpM2p.toFixed(3).padStart(8));
console.log();

let m2pWins = 0;
console.log('  Per-galaxy:');
console.log('  ' + 'Galaxy'.padEnd(14) + 'nPts'.padStart(5) + 'Obs'.padStart(7) + 'M2p'.padStart(7) + 'Err'.padStart(8) + 'M0'.padStart(7) + 'M2p<M0?'.padStart(9));
console.log('  ' + '-'.repeat(57));
for (const g of extGals) {
  const o = g.logA0;
  const p = betaM2p[0] + betaM2p[1] * g.logMHI + betaM2p[2] * g.logMeanRun;
  const e = p - o;
  const wins = Math.abs(e) < Math.abs(trainMeanA0 - o);
  if (wins) m2pWins++;
  console.log('  ' + g.name.padEnd(14) + String(g.nPts).padStart(5) + o.toFixed(3).padStart(7) +
    p.toFixed(3).padStart(7) + ((e >= 0 ? '+' : '') + e.toFixed(3)).padStart(8) +
    trainMeanA0.toFixed(3).padStart(7) + (wins ? 'YES' : 'no').padStart(9));
}
console.log();
console.log('  M2\' wins: ' + m2pWins + '/' + extGals.length + ' (' + (m2pWins / extGals.length * 100).toFixed(1) + '%)');
console.log();

console.log('======================================================================');
console.log('  TEST 2: M3_fair = logMHI + logMhost_vflat + logMeanRun');
console.log('  (Refit on training using Vflat-estimated Mhost everywhere)');
console.log('======================================================================\n');

const trainVflat = trainData.map(d => parseFloat(d.Vflat_km_s));
const trainLogMhostVflat = trainVflat.map(v => v > 0 ? 4.0 * Math.log10(v) + 2.4 : 10.0);
const trainX_M3f = trainData.map((d, i) => [parseFloat(d.logMHI), trainLogMhostVflat[i], parseFloat(d.logMeanRun)]);
const betaM3f = olsFit(trainY, trainX_M3f);
console.log('  M3_fair (fit on training with Vflat Mhost):');
console.log('    intercept  = ' + betaM3f[0].toFixed(4));
console.log('    logMHI     = ' + betaM3f[1].toFixed(4));
console.log('    logMhost_vf= ' + betaM3f[2].toFixed(4));
console.log('    logMeanRun = ' + betaM3f[3].toFixed(4));
console.log();

const extWithMhost = extGals.filter(g => g.logMhost_vflat !== null);
console.log('  External with valid Vflat Mhost: ' + extWithMhost.length);

const extObs3f = extWithMhost.map(g => g.logA0);
const extPredM3f = extWithMhost.map(g =>
  betaM3f[0] + betaM3f[1] * g.logMHI + betaM3f[2] * g.logMhost_vflat + betaM3f[3] * g.logMeanRun
);
const extPredM0_3f = extWithMhost.map(() => trainMeanA0);

const sdExt3f = sd(extObs3f);
const rmsM0_3f = Math.sqrt(extPredM0_3f.reduce((s, p, i) => s + (extObs3f[i] - p) ** 2, 0) / extWithMhost.length);
const rmsM3f = Math.sqrt(extPredM3f.reduce((s, p, i) => s + (extObs3f[i] - p) ** 2, 0) / extWithMhost.length);
const biasM3f = extPredM3f.reduce((s, p, i) => s + (p - extObs3f[i]), 0) / extWithMhost.length;
const gapM3f = 100 * (1 - rmsM3f ** 2 / sdExt3f ** 2);
const rsM3f = spearman(extObs3f, extPredM3f);
const rpM3f = pearson(extObs3f, extPredM3f);

console.log('  ' + 'Model'.padEnd(10) + 'RMS'.padStart(8) + 'Bias'.padStart(8) + 'Gap%'.padStart(8) + 'rs'.padStart(8) + 'rp'.padStart(8));
console.log('  ' + '-'.repeat(50));
console.log('  ' + 'M0'.padEnd(10) + rmsM0_3f.toFixed(4).padStart(8));
console.log('  ' + 'M3_fair'.padEnd(10) + rmsM3f.toFixed(4).padStart(8) + biasM3f.toFixed(4).padStart(8) + (gapM3f.toFixed(1) + '%').padStart(8) + rsM3f.toFixed(3).padStart(8) + rpM3f.toFixed(3).padStart(8));
console.log();

let m3fWins = 0;
for (let i = 0; i < extWithMhost.length; i++) {
  if (Math.abs(extPredM3f[i] - extObs3f[i]) < Math.abs(trainMeanA0 - extObs3f[i])) m3fWins++;
}
console.log('  M3_fair wins: ' + m3fWins + '/' + extWithMhost.length + ' (' + (m3fWins / extWithMhost.length * 100).toFixed(1) + '%)');
console.log();

console.log('======================================================================');
console.log('  TEST 3: CDS-DERIVED M2\' ON TRAINING (sanity check)');
console.log('  Use CDS-computed logMHI+logMeanRun for training galaxies');
console.log('======================================================================\n');

const trainFromCDS = trainGals.filter(g => g.logA0 >= 2.0);
console.log('  Training galaxies with CDS data: ' + trainFromCDS.length);

const cdsTrainObs = trainFromCDS.map(g => g.logA0);
const cdsTrainPred = trainFromCDS.map(g => betaM2p[0] + betaM2p[1] * g.logMHI + betaM2p[2] * g.logMeanRun);
const cdsPredM0 = trainFromCDS.map(() => trainMeanA0);
const rmsM0_cds = Math.sqrt(cdsPredM0.reduce((s, p, i) => s + (cdsTrainObs[i] - p) ** 2, 0) / trainFromCDS.length);
const rmsM2p_cds = Math.sqrt(cdsTrainPred.reduce((s, p, i) => s + (cdsTrainObs[i] - p) ** 2, 0) / trainFromCDS.length);

console.log('  M0 RMS (CDS a0) = ' + rmsM0_cds.toFixed(4));
console.log('  M2\' RMS (CDS a0, CDS MHI+MR) = ' + rmsM2p_cds.toFixed(4));
console.log('  M2\' beats M0? ' + (rmsM2p_cds < rmsM0_cds ? 'YES' : 'NO'));
console.log();

console.log('======================================================================');
console.log('  TEST 4: Q=1 SUBSAMPLE for M2\'');
console.log('======================================================================\n');

const q1 = extGals.filter(g => g.Q === 1);
if (q1.length >= 3) {
  const q1Obs = q1.map(g => g.logA0);
  const q1Pred = q1.map(g => betaM2p[0] + betaM2p[1] * g.logMHI + betaM2p[2] * g.logMeanRun);
  const q1M0 = q1.map(() => trainMeanA0);
  const sdQ1 = sd(q1Obs);
  const rmsQ1_M0 = Math.sqrt(q1M0.reduce((s, p, i) => s + (q1Obs[i] - p) ** 2, 0) / q1.length);
  const rmsQ1_M2p = Math.sqrt(q1Pred.reduce((s, p, i) => s + (q1Obs[i] - p) ** 2, 0) / q1.length);
  const gapQ1 = 100 * (1 - rmsQ1_M2p ** 2 / sdQ1 ** 2);
  let q1Wins = 0;
  for (let i = 0; i < q1.length; i++) {
    if (Math.abs(q1Pred[i] - q1Obs[i]) < Math.abs(trainMeanA0 - q1Obs[i])) q1Wins++;
  }
  console.log('  Q=1 external: N = ' + q1.length);
  console.log('  M0 RMS = ' + rmsQ1_M0.toFixed(4) + ', M2\' RMS = ' + rmsQ1_M2p.toFixed(4));
  console.log('  Gap% = ' + gapQ1.toFixed(1) + '%');
  console.log('  M2\' beats M0? ' + (rmsQ1_M2p < rmsQ1_M0 ? 'YES' : 'NO'));
  console.log('  M2\' wins per-galaxy: ' + q1Wins + '/' + q1.length + ' (' + (q1Wins / q1.length * 100).toFixed(1) + '%)');
  console.log('  Spearman = ' + spearman(q1Obs, q1Pred).toFixed(3));
}
console.log();

console.log('======================================================================');
console.log('  AXIS SIGN CHECK: logMHI and logMeanRun partial correlations');
console.log('======================================================================\n');

const extA0 = extGals.map(g => g.logA0);
const extMHI = extGals.map(g => g.logMHI);
const extMR = extGals.map(g => g.logMeanRun);

function simpleCorr(x, y) {
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < x.length; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxy / Math.sqrt(sxx * syy + 1e-20);
}

function partialCorr2(x, y, z) {
  const rxy = simpleCorr(x, y);
  const rxz = simpleCorr(x, z);
  const ryz = simpleCorr(y, z);
  return (rxy - rxz * ryz) / Math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2) + 1e-20);
}

const rMHI_raw = simpleCorr(extMHI, extA0);
const rMR_raw = simpleCorr(extMR, extA0);
const rMHI_partial = partialCorr2(extMHI, extA0, extMR);
const rMR_partial = partialCorr2(extMR, extA0, extMHI);

console.log('  External set (N=' + extGals.length + '):');
console.log('  ' + 'Axis'.padEnd(14) + 'Raw r'.padStart(8) + 'Partial r'.padStart(11) + 'Training sign'.padStart(15));
console.log('  ' + '-'.repeat(48));
console.log('  ' + 'logMHI'.padEnd(14) + rMHI_raw.toFixed(3).padStart(8) + rMHI_partial.toFixed(3).padStart(11) + '    negative'.padStart(15));
console.log('  ' + 'logMeanRun'.padEnd(14) + rMR_raw.toFixed(3).padStart(8) + rMR_partial.toFixed(3).padStart(11) + '    positive'.padStart(15));
console.log();

const mhiSignMatch = rMHI_partial < 0;
const mrSignMatch = rMR_partial > 0;
console.log('  Sign agreement: ' + (mhiSignMatch ? '✓' : '✗') + ' logMHI (' + (mhiSignMatch ? 'matches' : 'WRONG') + ')');
console.log('  Sign agreement: ' + (mrSignMatch ? '✓' : '✗') + ' logMeanRun (' + (mrSignMatch ? 'matches' : 'WRONG') + ')');
console.log();

console.log('======================================================================');
console.log('  OVERALL SUMMARY');
console.log('======================================================================\n');

const m2pBeatsPer = m2pWins / extGals.length;
const m2pBeatsRMS = rmsM2p < rmsM0;

let verdict;
if (m2pBeatsRMS && m2pBeatsPer > 0.55 && gapM2p > 5 && mhiSignMatch && mrSignMatch) {
  verdict = 'PASS';
} else if (m2pBeatsRMS && m2pBeatsPer > 0.5 && mhiSignMatch) {
  verdict = 'PARTIAL-PASS';
} else if ((mhiSignMatch || mrSignMatch) && m2pBeatsPer >= 0.45) {
  verdict = 'MARGINAL';
} else {
  verdict = 'FAIL';
}

console.log('  M2\' (logMHI + logMeanRun, NO Mhost):');
console.log('    RMS M2\' < RMS M0? ' + (m2pBeatsRMS ? 'YES' : 'NO') + ' (' + rmsM2p.toFixed(3) + ' vs ' + rmsM0.toFixed(3) + ')');
console.log('    Per-galaxy win rate: ' + (m2pBeatsPer * 100).toFixed(1) + '%');
console.log('    Gap%: ' + gapM2p.toFixed(1) + '%');
console.log('    Spearman: ' + rsM2p.toFixed(3));
console.log('    Axis signs match: MHI ' + (mhiSignMatch ? 'YES' : 'NO') + ', MeanRun ' + (mrSignMatch ? 'YES' : 'NO'));
console.log('    VERDICT: ' + verdict);
console.log();

if (m3fWins !== undefined) {
  const m3fBeatsPer = m3fWins / extWithMhost.length;
  console.log('  M3_fair (Vflat Mhost everywhere):');
  console.log('    RMS M3_fair < RMS M0? ' + (rmsM3f < rmsM0_3f ? 'YES' : 'NO') + ' (' + rmsM3f.toFixed(3) + ' vs ' + rmsM0_3f.toFixed(3) + ')');
  console.log('    Per-galaxy win rate: ' + (m3fBeatsPer * 100).toFixed(1) + '%');
  console.log('    Spearman: ' + rsM3f.toFixed(3));
}
console.log();

const output = {
  phase: '81b',
  title: 'Diagnostic External Tests',
  m2prime: {
    description: 'logMHI + logMeanRun only (no Mhost)',
    frozenBeta: betaM2p,
    looGapPct: +looGap_M2p.toFixed(1),
    nExternal: extGals.length,
    rmsM0: +rmsM0.toFixed(4),
    rmsM2p: +rmsM2p.toFixed(4),
    biasM2p: +biasM2p.toFixed(4),
    gapPct: +gapM2p.toFixed(1),
    spearman: +rsM2p.toFixed(3),
    pearson: +rpM2p.toFixed(3),
    winRate: +(m2pBeatsPer * 100).toFixed(1),
  },
  m3fair: {
    description: 'logMHI + logMhost_vflat + logMeanRun (Vflat Mhost everywhere)',
    frozenBeta: betaM3f,
    nExternal: extWithMhost.length,
    rmsM0: +rmsM0_3f.toFixed(4),
    rmsM3f: +rmsM3f.toFixed(4),
    biasM3f: +biasM3f.toFixed(4),
    gapPct: +gapM3f.toFixed(1),
    spearman: +rsM3f.toFixed(3),
    pearson: +rpM3f.toFixed(3),
    winRate: +(m3fWins / extWithMhost.length * 100).toFixed(1),
  },
  axisSigns: {
    logMHI_partial: +rMHI_partial.toFixed(3),
    logMeanRun_partial: +rMR_partial.toFixed(3),
    mhiSignMatch,
    mrSignMatch,
  },
  verdict,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase81b-diagnostic-tests.json'), JSON.stringify(output, null, 2));
console.log('  Saved: public/phase81b-diagnostic-tests.json');

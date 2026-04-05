#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 84: SELECTION / COLLIDER AUDIT');
console.log('  What about the N=45 construction generates the structured law?');
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

function pearsonR(x, y) {
  const n = x.length;
  if (n < 3) return 0;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function looGapPct(Y, X) {
  const n = Y.length;
  const yMean = mean(Y);
  let m0_ss = 0, model_ss = 0;
  for (let i = 0; i < n; i++) {
    const yTrain = Y.filter((_, j) => j !== i);
    const xTrain = X.filter((_, j) => j !== i);
    const beta = olsFit(yTrain, xTrain);
    if (!beta) return -999;
    const pred = beta[0] + X[i].reduce((s, v, k) => s + v * beta[k + 1], 0);
    model_ss += (Y[i] - pred) ** 2;
    m0_ss += (Y[i] - mean(yTrain)) ** 2;
  }
  const perGalSS = Y.reduce((s, v) => s + (v - yMean) ** 2, 0);
  const m2Gap = Math.sqrt(perGalSS / n);
  const m0RMS = Math.sqrt(m0_ss / n);
  const modelRMS = Math.sqrt(model_ss / n);
  return (m0RMS - modelRMS) / (m0RMS - 0) * 100;
}

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

const train45 = allGalaxies.filter(g => g.inTraining);
const external = allGalaxies.filter(g => !g.inTraining);

console.log('  Total galaxies: ' + allGalaxies.length);
console.log('  Training N=45: ' + train45.length);
console.log('  External: ' + external.length);
console.log();

const trainOrigCSV = trainCSV.map(d => ({
  name: d.galaxy_name,
  logA0: parseFloat(d.logA0),
  logMHI: parseFloat(d.logMHI),
  logMhost: parseFloat(d.logMhost),
  logMeanRun: parseFloat(d.logMeanRun),
}));

function runM2pOnSample(sample) {
  if (sample.length < 6) return null;
  const Y = sample.map(g => g.logA0);
  const X = sample.map(g => [g.logMHI, g.logMeanRun]);
  const beta = olsFit(Y, X);
  if (!beta) return null;
  const yMean = mean(Y);
  let m0_ss = 0, model_ss = 0, loo_m0 = 0, loo_mod = 0;
  for (let i = 0; i < sample.length; i++) {
    const pred = beta[0] + beta[1] * X[i][0] + beta[2] * X[i][1];
    model_ss += (Y[i] - pred) ** 2;
    m0_ss += (Y[i] - yMean) ** 2;
  }
  for (let i = 0; i < sample.length; i++) {
    const yTr = Y.filter((_, j) => j !== i);
    const xTr = X.filter((_, j) => j !== i);
    const b = olsFit(yTr, xTr);
    if (!b) return null;
    const pred = b[0] + b[1] * X[i][0] + b[2] * X[i][1];
    loo_mod += (Y[i] - pred) ** 2;
    loo_m0 += (Y[i] - mean(yTr)) ** 2;
  }
  const rmsM0 = Math.sqrt(loo_m0 / sample.length);
  const rmsModel = Math.sqrt(loo_mod / sample.length);
  const gap = (rmsM0 - rmsModel) / rmsM0 * 100;
  const r_mhi = pearsonR(sample.map(g => g.logMHI), Y);
  const r_run = pearsonR(sample.map(g => g.logMeanRun), Y);
  return {
    n: sample.length,
    beta,
    rmsM0: +rmsM0.toFixed(3), rmsModel: +rmsModel.toFixed(3),
    gapPct: +gap.toFixed(1),
    r_logMHI: +r_mhi.toFixed(3), r_logMeanRun: +r_run.toFixed(3),
    signsCorrect: r_mhi < 0 && r_run > 0,
    betaLogMHI: +beta[1].toFixed(4), betaLogMeanRun: +beta[2].toFixed(4),
  };
}

console.log('======================================================================');
console.log('  TEST 1: Selection criteria — one-at-a-time relaxation');
console.log('  Start from all SPARC galaxies, apply each N=45 criterion');
console.log('======================================================================\n');

const selectionCriteria = [
  { name: 'Published distance (fD>=2)', filter: g => g.fD >= 2 },
  { name: 'Quality Q<=2', filter: g => g.Q <= 2 },
  { name: 'nPts >= 10', filter: g => g.nPts >= 10 },
  { name: 'Vflat >= 80 km/s', filter: g => g.Vflat >= 80 },
  { name: 'T <= 8 (not Sdm/Sm/Im)', filter: g => g.T <= 8 },
  { name: 'logA0 in [3.0, 4.2]', filter: g => g.logA0 >= 3.0 && g.logA0 <= 4.2 },
  { name: 'In groups/clusters', filter: g => g.envCode >= 1 },
];

console.log('  ── Baseline: ALL SPARC (N=' + allGalaxies.length + ') ──');
const baselineAll = runM2pOnSample(allGalaxies);
if (baselineAll) {
  console.log('    M2\' LOO gap = ' + baselineAll.gapPct + '%  r(MHI)=' + baselineAll.r_logMHI + '  r(MeanRun)=' + baselineAll.r_logMeanRun + '  Signs=' + (baselineAll.signsCorrect ? 'OK' : 'WRONG'));
}
console.log();

console.log('  ── Baseline: Training N=45 (using CDS-derived values) ──');
const baseline45 = runM2pOnSample(train45);
if (baseline45) {
  console.log('    M2\' LOO gap = ' + baseline45.gapPct + '%  r(MHI)=' + baseline45.r_logMHI + '  r(MeanRun)=' + baseline45.r_logMeanRun + '  Signs=' + (baseline45.signsCorrect ? 'OK' : 'WRONG'));
}
console.log();

const criteriaResults = [];

for (const c of selectionCriteria) {
  const subset = allGalaxies.filter(c.filter);
  const result = runM2pOnSample(subset);
  if (result) {
    console.log('  ── ' + c.name + ' (N=' + subset.length + ') ──');
    console.log('    M2\' LOO gap = ' + result.gapPct + '%  r(MHI)=' + result.r_logMHI + '  r(MeanRun)=' + result.r_logMeanRun + '  Signs=' + (result.signsCorrect ? 'OK' : 'WRONG'));
    console.log('    beta: MHI=' + result.betaLogMHI + '  MeanRun=' + result.betaLogMeanRun);
    criteriaResults.push({ criterion: c.name, ...result });
  }
  console.log();
}

console.log('======================================================================');
console.log('  TEST 2: Cumulative selection — build toward N=45 step by step');
console.log('======================================================================\n');

const cumulativeSteps = [
  { name: 'All SPARC', filters: [] },
  { name: '+ fD>=2', filters: [g => g.fD >= 2] },
  { name: '+ Q<=2', filters: [g => g.fD >= 2, g => g.Q <= 2] },
  { name: '+ nPts>=10', filters: [g => g.fD >= 2, g => g.Q <= 2, g => g.nPts >= 10] },
  { name: '+ Vflat>=80', filters: [g => g.fD >= 2, g => g.Q <= 2, g => g.nPts >= 10, g => g.Vflat >= 80] },
  { name: '+ T<=8', filters: [g => g.fD >= 2, g => g.Q <= 2, g => g.nPts >= 10, g => g.Vflat >= 80, g => g.T <= 8] },
  { name: '+ envCode>=1', filters: [g => g.fD >= 2, g => g.Q <= 2, g => g.nPts >= 10, g => g.Vflat >= 80, g => g.T <= 8, g => g.envCode >= 1] },
];

const cumulativeResults = [];

for (const step of cumulativeSteps) {
  let subset = [...allGalaxies];
  for (const f of step.filters) subset = subset.filter(f);
  const result = runM2pOnSample(subset);
  if (result) {
    const nIn45 = subset.filter(g => g.inTraining).length;
    console.log('  ' + step.name.padEnd(20) + '  N=' + String(subset.length).padStart(3) + '  (of which ' + nIn45 + ' in N=45)' +
      '  Gap=' + String(result.gapPct).padStart(6) + '%' +
      '  r(MHI)=' + String(result.r_logMHI).padStart(7) +
      '  r(Run)=' + String(result.r_logMeanRun).padStart(7) +
      '  Signs=' + (result.signsCorrect ? 'OK' : 'WRONG'));
    cumulativeResults.push({ step: step.name, nTotal: subset.length, nIn45, ...result });
  }
}
console.log();

console.log('======================================================================');
console.log('  TEST 3: "In-training" indicator — does sample membership absorb');
console.log('  the signal? If adding inTraining flag kills MHI/MeanRun, the');
console.log('  signal comes from the selection, not from physics.');
console.log('======================================================================\n');

const pooled = allGalaxies.filter(g => g.fD >= 2 && g.nPts >= 5);
console.log('  Pooled sample (fD>=2, nPts>=5): N = ' + pooled.length);

const Y_pool = pooled.map(g => g.logA0);
const X_noFlag = pooled.map(g => [g.logMHI, g.logMeanRun]);
const X_withFlag = pooled.map(g => [g.logMHI, g.logMeanRun, g.inTraining ? 1 : 0]);
const X_flagOnly = pooled.map(g => [g.inTraining ? 1 : 0]);
const X_flagInteract = pooled.map(g => [
  g.logMHI, g.logMeanRun, g.inTraining ? 1 : 0,
  g.logMHI * (g.inTraining ? 1 : 0),
  g.logMeanRun * (g.inTraining ? 1 : 0),
]);

const models = [
  { name: 'M2\' (MHI + MeanRun)', X: X_noFlag },
  { name: 'M2\' + inTraining flag', X: X_withFlag },
  { name: 'inTraining flag ONLY', X: X_flagOnly },
  { name: 'M2\' + flag + interactions', X: X_flagInteract },
];

for (const m of models) {
  const beta = olsFit(Y_pool, m.X);
  if (!beta) { console.log('  ' + m.name + ': OLS failed'); continue; }

  let loo_m0 = 0, loo_mod = 0;
  for (let i = 0; i < pooled.length; i++) {
    const yTr = Y_pool.filter((_, j) => j !== i);
    const xTr = m.X.filter((_, j) => j !== i);
    const b = olsFit(yTr, xTr);
    if (!b) { loo_m0 = -1; break; }
    const pred = b[0] + m.X[i].reduce((s, v, k) => s + v * b[k + 1], 0);
    loo_mod += (Y_pool[i] - pred) ** 2;
    loo_m0 += (Y_pool[i] - mean(yTr)) ** 2;
  }
  if (loo_m0 < 0) { console.log('  ' + m.name + ': LOO failed'); continue; }

  const rmsM0 = Math.sqrt(loo_m0 / pooled.length);
  const rmsModel = Math.sqrt(loo_mod / pooled.length);
  const gap = (rmsM0 - rmsModel) / rmsM0 * 100;

  console.log('  ' + m.name.padEnd(35) + '  Gap=' + gap.toFixed(1) + '%');
  console.log('    beta: ' + beta.map((b, i) => (i === 0 ? 'int' : 'b' + i) + '=' + b.toFixed(3)).join(', '));
}
console.log();

console.log('======================================================================');
console.log('  TEST 4: Correlation structure — does r(MHI, a0) exist in full');
console.log('  SPARC or only appear when selected?');
console.log('======================================================================\n');

const subsets = [
  { name: 'All SPARC', data: allGalaxies },
  { name: 'fD >= 2 only', data: allGalaxies.filter(g => g.fD >= 2) },
  { name: 'fD >= 2, Q <= 2', data: allGalaxies.filter(g => g.fD >= 2 && g.Q <= 2) },
  { name: 'Training N=45', data: train45 },
  { name: 'External (all)', data: external },
  { name: 'External fD>=2', data: external.filter(g => g.fD >= 2) },
  { name: 'fD==1 only (Hubble flow)', data: allGalaxies.filter(g => g.fD === 1) },
];

console.log('  ' + 'Subset'.padEnd(28) + 'N'.padStart(5) + 'r(MHI,a0)'.padStart(11) + 'r(Run,a0)'.padStart(11) + 'r(Vflat,a0)'.padStart(12) + 'sd(a0)'.padStart(8));
for (const s of subsets) {
  if (s.data.length < 5) continue;
  const r_mhi = pearsonR(s.data.map(g => g.logMHI), s.data.map(g => g.logA0));
  const r_run = pearsonR(s.data.map(g => g.logMeanRun), s.data.map(g => g.logA0));
  const r_vf = pearsonR(s.data.map(g => Math.log10(g.Vflat > 0 ? g.Vflat : 50)), s.data.map(g => g.logA0));
  const sdA0 = sd(s.data.map(g => g.logA0));
  console.log('  ' + s.name.padEnd(28) + String(s.data.length).padStart(5) +
    r_mhi.toFixed(3).padStart(11) + r_run.toFixed(3).padStart(11) + r_vf.toFixed(3).padStart(12) + sdA0.toFixed(3).padStart(8));
}
console.log();

console.log('======================================================================');
console.log('  TEST 5: Random subsample stability');
console.log('  Draw 1000 random N=45 from full SPARC, compute M2\' r-values');
console.log('  How often does a random sample show training-like correlations?');
console.log('======================================================================\n');

const N_PERM = 1000;
const trainR_MHI = pearsonR(train45.map(g => g.logMHI), train45.map(g => g.logA0));
const trainR_Run = pearsonR(train45.map(g => g.logMeanRun), train45.map(g => g.logA0));

let nStrongerMHI = 0, nStrongerRun = 0, nBothCorrectSign = 0;
const randR_MHI = [], randR_Run = [];

for (let p = 0; p < N_PERM; p++) {
  const shuffled = [...allGalaxies].sort(() => Math.random() - 0.5);
  const sample = shuffled.slice(0, 45);
  const r_mhi = pearsonR(sample.map(g => g.logMHI), sample.map(g => g.logA0));
  const r_run = pearsonR(sample.map(g => g.logMeanRun), sample.map(g => g.logA0));
  randR_MHI.push(r_mhi);
  randR_Run.push(r_run);
  if (r_mhi <= trainR_MHI) nStrongerMHI++;
  if (r_run >= trainR_Run) nStrongerRun++;
  if (r_mhi < 0 && r_run > 0) nBothCorrectSign++;
}

console.log('  Training r(MHI, a0) = ' + trainR_MHI.toFixed(3));
console.log('  Random N=45 r(MHI, a0): mean = ' + mean(randR_MHI).toFixed(3) + ' ± ' + sd(randR_MHI).toFixed(3));
console.log('  P(r <= training): ' + (nStrongerMHI / N_PERM).toFixed(3) + ' (' + nStrongerMHI + '/' + N_PERM + ')');
console.log();
console.log('  Training r(MeanRun, a0) = ' + trainR_Run.toFixed(3));
console.log('  Random N=45 r(MeanRun, a0): mean = ' + mean(randR_Run).toFixed(3) + ' ± ' + sd(randR_Run).toFixed(3));
console.log('  P(r >= training): ' + (nStrongerRun / N_PERM).toFixed(3) + ' (' + nStrongerRun + '/' + N_PERM + ')');
console.log();
console.log('  P(both signs correct in random N=45): ' + (nBothCorrectSign / N_PERM * 100).toFixed(1) + '%');
console.log();

console.log('======================================================================');
console.log('  TEST 6: Distance method as collider');
console.log('  Does conditioning on fD>=2 create artificial a0-MHI correlation?');
console.log('======================================================================\n');

const fdGroups = [
  { name: 'fD=1 (Hubble flow)', data: allGalaxies.filter(g => g.fD === 1) },
  { name: 'fD=2 (TRGB)', data: allGalaxies.filter(g => g.fD === 2) },
  { name: 'fD=3 (Cepheids)', data: allGalaxies.filter(g => g.fD === 3) },
  { name: 'fD=4 (Ursa Major)', data: allGalaxies.filter(g => g.fD === 4) },
  { name: 'fD>=2 (all published)', data: allGalaxies.filter(g => g.fD >= 2) },
];

console.log('  ' + 'Distance method'.padEnd(25) + 'N'.padStart(5) + 'r(MHI,a0)'.padStart(11) + 'r(Run,a0)'.padStart(11) + 'mean(a0)'.padStart(10) + 'sd(a0)'.padStart(8));
for (const fg of fdGroups) {
  if (fg.data.length < 5) {
    console.log('  ' + fg.name.padEnd(25) + String(fg.data.length).padStart(5) + '  (too few)');
    continue;
  }
  const r_mhi = pearsonR(fg.data.map(g => g.logMHI), fg.data.map(g => g.logA0));
  const r_run = pearsonR(fg.data.map(g => g.logMeanRun), fg.data.map(g => g.logA0));
  const ma0 = mean(fg.data.map(g => g.logA0));
  const sa0 = sd(fg.data.map(g => g.logA0));
  console.log('  ' + fg.name.padEnd(25) + String(fg.data.length).padStart(5) +
    r_mhi.toFixed(3).padStart(11) + r_run.toFixed(3).padStart(11) + ma0.toFixed(3).padStart(10) + sa0.toFixed(3).padStart(8));
}
console.log();

console.log('======================================================================');
console.log('  FINAL INTERPRETATION');
console.log('======================================================================\n');

const fullR_MHI = pearsonR(allGalaxies.map(g => g.logMHI), allGalaxies.map(g => g.logA0));
const onlyN45_R_MHI = trainR_MHI;
const extR_MHI = pearsonR(external.map(g => g.logMHI), external.map(g => g.logA0));

console.log('  r(logMHI, logA0) across subsets:');
console.log('    Full SPARC (N=' + allGalaxies.length + '): ' + fullR_MHI.toFixed(3));
console.log('    Training N=45:        ' + onlyN45_R_MHI.toFixed(3));
console.log('    External N=' + external.length + ':        ' + extR_MHI.toFixed(3));
console.log();

if (fullR_MHI < 0 && Math.abs(fullR_MHI) > 0.15) {
  console.log('  The MHI-a0 correlation EXISTS in full SPARC → selection AMPLIFIES it');
} else if (fullR_MHI < 0) {
  console.log('  Weak MHI-a0 correlation in full SPARC → selection may be amplifying a noise feature');
} else {
  console.log('  No MHI-a0 correlation in full SPARC → selection CREATES an artificial correlation');
}

if (extR_MHI > 0) {
  console.log('  Sign REVERSAL on external → classic collider/selection artifact signature');
}
console.log();

const output = {
  phase: '84',
  title: 'Selection / Collider Audit',
  nTotal: allGalaxies.length,
  nTraining: train45.length,
  nExternal: external.length,
  baselineAll: baselineAll,
  baseline45: baseline45,
  singleCriterionTests: criteriaResults,
  cumulativeSelection: cumulativeResults,
  correlationStructure: subsets.map(s => ({
    name: s.name, n: s.data.length,
    r_logMHI: +pearsonR(s.data.map(g => g.logMHI), s.data.map(g => g.logA0)).toFixed(3),
    r_logMeanRun: +pearsonR(s.data.map(g => g.logMeanRun), s.data.map(g => g.logA0)).toFixed(3),
  })).filter(s => s.n >= 5),
  randomSubsampleTest: {
    nPermutations: N_PERM,
    training_r_MHI: +trainR_MHI.toFixed(3),
    training_r_MeanRun: +trainR_Run.toFixed(3),
    pValueMHI: +(nStrongerMHI / N_PERM).toFixed(3),
    pValueMeanRun: +(nStrongerRun / N_PERM).toFixed(3),
    pctBothCorrectSign: +(nBothCorrectSign / N_PERM * 100).toFixed(1),
  },
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase84-selection-audit.json'), JSON.stringify(output, null, 2));
console.log('  Saved: public/phase84-selection-audit.json');

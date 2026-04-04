#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 81: CLEAN EXTERNAL TEST — SPARC-out N=27+');
console.log('  Using CDS/VizieR rotation curve data (Lelli+2016 Table 2)');
console.log('  Frozen M3 coefficients applied to galaxies OUTSIDE training N=45');
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
      const pred1 = mcgaughRAR(gb, Math.pow(10, m1));
      const pred2 = mcgaughRAR(gb, Math.pow(10, m2));
      c1 += (p.logGobs - Math.log10(pred1 > 0 ? pred1 : 1e-20)) ** 2;
      c2 += (p.logGobs - Math.log10(pred2 > 0 ? pred2 : 1e-20)) ** 2;
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

console.log('  Training set: N = ' + trainNames.size + ' galaxies');
console.log('  Training names sample: ' + [...trainNames].slice(0, 5).join(', ') + ' ...\n');

const trainY = trainData.map(d => parseFloat(d.logA0));
const trainMeanA0 = mean(trainY);
const M3_BETA = [5.182, -0.198, -0.155, 0.459];
console.log('  FROZEN M3: intercept=' + M3_BETA[0] + ', logMHI=' + M3_BETA[1] +
  ', logMhost=' + M3_BETA[2] + ', logMeanRun=' + M3_BETA[3]);
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STEP 1: Parse CDS Table 2 (rotation curves for 175 galaxies)');
console.log('----------------------------------------------------------------------\n');

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const D = parseFloat(line.substring(15, 21).trim());
  const eD = parseFloat(line.substring(22, 27).trim());
  const fD = parseInt(line.substring(28, 29).trim());
  const inc = parseFloat(line.substring(30, 34).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const SBdisk = parseFloat(line.substring(77, 85).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { T, D, eD, fD, inc, L36, Rdisk, SBdisk, MHI, Vflat, Q };
}
console.log('  Parsed Table 1: ' + Object.keys(table1).length + ' galaxies');

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const dist = parseFloat(line.substring(12, 18).trim());
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const evobs = parseFloat(line.substring(33, 38).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, evobs, vgas, vdisk, vbul });
}
console.log('  Parsed Table 2: ' + Object.keys(rcByGalaxy).length + ' galaxies with RC data');
console.log('  Total RC points: ' + table2Raw.length);
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STEP 2: Compute a0 + MeanRun for each galaxy from raw RC');
console.log('----------------------------------------------------------------------\n');

const sparcAll = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf-8'));
const p25 = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase25-group-membership.json'), 'utf-8'));
const envLookup = {};
for (const g of p25.galaxyAssignments) envLookup[g.name] = g;

const allGalaxies = [];

for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1) continue;
  if (!t1.MHI || t1.MHI <= 0) continue;
  if (!t1.L36 || t1.L36 <= 0) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    const gObs = pt.vobs * pt.vobs / pt.rad;
    const gBar = Math.abs(vBarSq) / pt.rad;
    if (gBar <= 0 || gObs <= 0 || !isFinite(gBar) || !isFinite(gObs)) continue;
    const logGobs = Math.log10(gObs);
    const logGbar = Math.log10(gBar);
    if (!isFinite(logGobs) || !isFinite(logGbar)) continue;
    pts.push({ r: pt.rad, logGobs, logGbar, vobs: pt.vobs });
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

  const sp = sparcAll.find(s => s.name === name);
  const ev = envLookup[name];

  const Rdisk = t1.Rdisk || 3.0;
  const Mbar_sun = t1.MHI * 1.33e9 + t1.L36 * 1e9 * 0.5;
  const Sigma0_bar = Mbar_sun / (Math.PI * (Rdisk * 1e3) ** 2);
  const logSigma0 = Math.log10(Sigma0_bar > 0 ? Sigma0_bar : 1e-3);

  const rArr = pts.map(p => p.r);
  const rMax = Math.max(...rArr);
  const rMid = rMax / 2;
  const innerPts = pts.filter(p => p.r <= rMid);
  const outerPts = pts.filter(p => p.r > rMid);

  function residualVariability(subset) {
    if (subset.length < 3) return 0;
    const xv = subset.map(p => p.logGbar);
    const yv = subset.map(p => p.logGobs);
    const mx = xv.reduce((s, v) => s + v, 0) / xv.length;
    const my = yv.reduce((s, v) => s + v, 0) / yv.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < xv.length; i++) { sxy += (xv[i] - mx) * (yv[i] - my); sxx += (xv[i] - mx) ** 2; }
    const b = sxx > 0 ? sxy / sxx : 0;
    const a = my - b * mx;
    let ss = 0;
    for (let i = 0; i < xv.length; i++) ss += (yv[i] - (a + b * xv[i])) ** 2;
    return Math.sqrt(ss / xv.length);
  }

  const rcWiggliness = residualVariability(innerPts) + residualVariability(outerPts);

  const logUpsilon_disk = Math.log10(UPSILON_DISK);

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0, rms: fit.rms,
    logMHI, logMeanRun, logSigma0, rcWiggliness,
    logUpsilon_disk,
    Vflat: t1.Vflat, Q: t1.Q, D: t1.D, fD: t1.fD, T: t1.T,
    inTraining: trainNames.has(name),
  });
}

console.log('  Total galaxies with valid a0 + MeanRun: ' + allGalaxies.length);
const inTrain = allGalaxies.filter(g => g.inTraining);
const external = allGalaxies.filter(g => !g.inTraining);
console.log('  In training set: ' + inTrain.length);
console.log('  External (not in training): ' + external.length);
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STEP 3: VALIDATE — Compare CDS-derived a0 with training a0');
console.log('----------------------------------------------------------------------\n');

let nMatch = 0, sumDiff = 0, sumDiffSq = 0;
for (const tg of inTrain) {
  const orig = trainData.find(d => d.galaxy_name === tg.name);
  if (orig) {
    const diff = tg.logA0 - parseFloat(orig.logA0);
    sumDiff += diff;
    sumDiffSq += diff ** 2;
    nMatch++;
    if (Math.abs(diff) > 0.3) {
      console.log('  WARNING: Large a0 diff for ' + tg.name + ': CDS=' + tg.logA0.toFixed(3) + ' train=' + parseFloat(orig.logA0).toFixed(3) + ' diff=' + diff.toFixed(3));
    }
  }
}
const meanDiff = sumDiff / nMatch;
const rmsDiff = Math.sqrt(sumDiffSq / nMatch);
console.log('  Matched training galaxies: ' + nMatch);
console.log('  a0 CDS vs original: mean diff = ' + meanDiff.toFixed(4) + ', RMS diff = ' + rmsDiff.toFixed(4) + ' dex');
if (rmsDiff > 0.2) {
  console.log('  *** WARNING: RMS diff > 0.2 dex — CDS pipeline may not match original ***');
}
console.log();

console.log('  MeanRun validation:');
let nMRmatch = 0, sumMRdiff = 0, sumMRdiffSq = 0;
for (const tg of inTrain) {
  const orig = trainData.find(d => d.galaxy_name === tg.name);
  if (orig) {
    const diff = tg.logMeanRun - parseFloat(orig.logMeanRun);
    sumMRdiff += diff;
    sumMRdiffSq += diff ** 2;
    nMRmatch++;
  }
}
const meanMRdiff = sumMRdiff / nMRmatch;
const rmsMRdiff = Math.sqrt(sumMRdiffSq / nMRmatch);
console.log('  logMeanRun CDS vs original: mean diff = ' + meanMRdiff.toFixed(4) + ', RMS diff = ' + rmsMRdiff.toFixed(4));
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STEP 4: APPLY FROZEN M3 TO EXTERNAL SAMPLE');
console.log('----------------------------------------------------------------------\n');

const extFiltered = external.filter(g => g.fD >= 2 && g.logA0 >= 2.0 && g.nPts >= 5);
console.log('  External galaxies after quality filter (fD>=2, logA0>=2, nPts>=5): ' + extFiltered.length);
console.log();

const spLookup = {};
for (const g of p25.galaxyAssignments) spLookup[g.name] = g;

for (const g of extFiltered) {
  const grp = spLookup[g.name];
  if (grp && grp.logMhost && isFinite(grp.logMhost)) {
    g.logMhost = grp.logMhost;
    g.mhostSource = 'group-catalog';
  } else {
    if (g.Vflat > 0) {
      g.logMhost = 4.0 * Math.log10(g.Vflat) + 2.4;
      g.mhostSource = 'Vflat-proxy';
    } else {
      g.logMhost = null;
      g.mhostSource = 'none';
    }
  }
}

const extValid = extFiltered.filter(g => g.logMhost !== null && isFinite(g.logMhost));
console.log('  External with valid Mhost: ' + extValid.length);
const nGroupCat = extValid.filter(g => g.mhostSource === 'group-catalog').length;
const nVflatProxy = extValid.filter(g => g.mhostSource === 'Vflat-proxy').length;
console.log('    from group catalog: ' + nGroupCat);
console.log('    from Vflat proxy: ' + nVflatProxy);
console.log();

if (extValid.length < 5) {
  console.log('  *** TOO FEW EXTERNAL GALAXIES — cannot run M3 test ***');
  process.exit(1);
}

const obs = extValid.map(g => g.logA0);
const sdObs = sd(obs);
const predM0 = extValid.map(() => trainMeanA0);
const predM3 = extValid.map(g =>
  M3_BETA[0] + M3_BETA[1] * g.logMHI + M3_BETA[2] * g.logMhost + M3_BETA[3] * g.logMeanRun
);

function rmsCalc(pred) { return Math.sqrt(pred.reduce((s, p, i) => s + (obs[i] - p) ** 2, 0) / obs.length); }
function biasCalc(pred) { return pred.reduce((s, p, i) => s + (p - obs[i]), 0) / obs.length; }
function pearson(pred) {
  const mo = mean(obs), mp = mean(pred);
  const num = obs.reduce((s, v, i) => s + (v - mo) * (pred[i] - mp), 0);
  const d1 = Math.sqrt(obs.reduce((s, v) => s + (v - mo) ** 2, 0));
  const d2 = Math.sqrt(pred.reduce((s, v) => s + (v - mp) ** 2, 0));
  return num / (d1 * d2 + 1e-12);
}
function spearman(pred) {
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

const rmsM0 = rmsCalc(predM0);
const rmsM3 = rmsCalc(predM3);
const biasM3 = biasCalc(predM3);
const gapM3 = 100 * (1 - rmsM3 ** 2 / sdObs ** 2);
const rsM3 = spearman(predM3);
const rpM3 = pearson(predM3);

console.log('  SD(observed logA0) = ' + sdObs.toFixed(4) + ' dex\n');
console.log('  ' + 'Model'.padEnd(10) + 'RMS'.padStart(8) + 'Bias'.padStart(8) + 'Gap%'.padStart(8) + 'Spearman'.padStart(10) + 'Pearson'.padStart(9));
console.log('  ' + '-'.repeat(53));
console.log('  ' + 'M0'.padEnd(10) + rmsM0.toFixed(4).padStart(8) + biasCalc(predM0).toFixed(4).padStart(8) + '0.0%'.padStart(8) + '—'.padStart(10) + '—'.padStart(9));
console.log('  ' + 'M3_frozen'.padEnd(10) + rmsM3.toFixed(4).padStart(8) + biasM3.toFixed(4).padStart(8) + (gapM3.toFixed(1) + '%').padStart(8) + rsM3.toFixed(3).padStart(10) + rpM3.toFixed(3).padStart(9));
console.log();

let m3wins = 0;
console.log('  Per-galaxy frozen M3 predictions:');
console.log('  ' + 'Galaxy'.padEnd(14) + 'nPts'.padStart(5) + 'Obs'.padStart(7) + 'PredM3'.padStart(8) + 'Err'.padStart(8) + 'PredM0'.padStart(8) + 'M3<M0?'.padStart(8) + 'MhSrc'.padStart(8));
console.log('  ' + '-'.repeat(66));
for (const g of extValid) {
  const o = g.logA0;
  const p3 = M3_BETA[0] + M3_BETA[1] * g.logMHI + M3_BETA[2] * g.logMhost + M3_BETA[3] * g.logMeanRun;
  const p0 = trainMeanA0;
  const e3 = p3 - o;
  const e0 = p0 - o;
  const wins = Math.abs(e3) < Math.abs(e0);
  if (wins) m3wins++;
  console.log('  ' + g.name.padEnd(14) + String(g.nPts).padStart(5) + o.toFixed(3).padStart(7) +
    p3.toFixed(3).padStart(8) + ((e3 >= 0 ? '+' : '') + e3.toFixed(3)).padStart(8) +
    p0.toFixed(3).padStart(8) + (wins ? 'YES' : 'no').padStart(8) + g.mhostSource.substring(0, 7).padStart(8));
}
console.log();
console.log('  M3 closer than M0: ' + m3wins + '/' + extValid.length + ' (' + (m3wins / extValid.length * 100).toFixed(1) + '%)');
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STEP 5: GROUP-CATALOG-ONLY SUBSAMPLE (cleanest Mhost)');
console.log('----------------------------------------------------------------------\n');

const gcOnly = extValid.filter(g => g.mhostSource === 'group-catalog');
if (gcOnly.length >= 3) {
  const obsGC = gcOnly.map(g => g.logA0);
  const predGC_M3 = gcOnly.map(g =>
    M3_BETA[0] + M3_BETA[1] * g.logMHI + M3_BETA[2] * g.logMhost + M3_BETA[3] * g.logMeanRun
  );
  const predGC_M0 = gcOnly.map(() => trainMeanA0);
  const rmsGC_M0 = Math.sqrt(predGC_M0.reduce((s, p, i) => s + (obsGC[i] - p) ** 2, 0) / gcOnly.length);
  const rmsGC_M3 = Math.sqrt(predGC_M3.reduce((s, p, i) => s + (obsGC[i] - p) ** 2, 0) / gcOnly.length);
  let gcWins = 0;
  for (let i = 0; i < gcOnly.length; i++) {
    if (Math.abs(predGC_M3[i] - obsGC[i]) < Math.abs(predGC_M0[i] - obsGC[i])) gcWins++;
  }
  console.log('  Group-catalog Mhost subsample: N = ' + gcOnly.length);
  console.log('  M0 RMS = ' + rmsGC_M0.toFixed(4) + ', M3 RMS = ' + rmsGC_M3.toFixed(4));
  console.log('  M3 beats M0? ' + (rmsGC_M3 < rmsGC_M0 ? 'YES' : 'NO'));
  console.log('  M3 closer per-galaxy: ' + gcWins + '/' + gcOnly.length);
} else {
  console.log('  Only ' + gcOnly.length + ' galaxies with group-catalog Mhost — too few for subsample');
}
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STEP 6: Q=1 SUBSAMPLE (highest quality RC)');
console.log('----------------------------------------------------------------------\n');

const q1 = extValid.filter(g => g.Q === 1);
if (q1.length >= 3) {
  const obsQ1 = q1.map(g => g.logA0);
  const predQ1_M3 = q1.map(g =>
    M3_BETA[0] + M3_BETA[1] * g.logMHI + M3_BETA[2] * g.logMhost + M3_BETA[3] * g.logMeanRun
  );
  const predQ1_M0 = q1.map(() => trainMeanA0);
  const sdQ1 = q1.length > 1 ? sd(obsQ1) : 0;
  const rmsQ1_M0 = Math.sqrt(predQ1_M0.reduce((s, p, i) => s + (obsQ1[i] - p) ** 2, 0) / q1.length);
  const rmsQ1_M3 = Math.sqrt(predQ1_M3.reduce((s, p, i) => s + (obsQ1[i] - p) ** 2, 0) / q1.length);
  const gapQ1 = sdQ1 > 0 ? 100 * (1 - rmsQ1_M3 ** 2 / sdQ1 ** 2) : 0;
  let q1wins = 0;
  for (let i = 0; i < q1.length; i++) {
    if (Math.abs(predQ1_M3[i] - obsQ1[i]) < Math.abs(predQ1_M0[i] - obsQ1[i])) q1wins++;
  }
  console.log('  Q=1 subsample: N = ' + q1.length);
  console.log('  SD(logA0) = ' + sdQ1.toFixed(4));
  console.log('  M0 RMS = ' + rmsQ1_M0.toFixed(4) + ', M3 RMS = ' + rmsQ1_M3.toFixed(4));
  console.log('  M3 gap% = ' + gapQ1.toFixed(1) + '%');
  console.log('  M3 beats M0? ' + (rmsQ1_M3 < rmsQ1_M0 ? 'YES' : 'NO'));
  console.log('  M3 closer per-galaxy: ' + q1wins + '/' + q1.length);
  console.log('  Spearman = ' + spearman(predQ1_M3).toFixed(3));
} else {
  console.log('  Only ' + q1.length + ' Q=1 galaxies — too few');
}
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STEP 7: AXIS-BY-AXIS PARTIAL CORRELATION ON EXTERNAL SET');
console.log('----------------------------------------------------------------------\n');

function partialCorr(x, y, controls) {
  if (controls.length === 0) {
    const mx = mean(x), my = mean(y);
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < x.length; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
    return sxy / Math.sqrt(sxx * syy + 1e-20);
  }
  const n = x.length;
  const X = Array.from({ length: n }, (_, i) => controls.map(c => c[i]));
  const regX = olsFit(x, X);
  const regY = olsFit(y, X);
  const resX = x.map((v, i) => v - (regX[0] + X[i].reduce((s, c, j) => s + regX[j + 1] * c, 0)));
  const resY = y.map((v, i) => v - (regY[0] + X[i].reduce((s, c, j) => s + regY[j + 1] * c, 0)));
  const mRX = mean(resX), mRY = mean(resY);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (resX[i] - mRX) * (resY[i] - mRY); sxx += (resX[i] - mRX) ** 2; syy += (resY[i] - mRY) ** 2; }
  return sxy / Math.sqrt(sxx * syy + 1e-20);
}

const extA0 = extValid.map(g => g.logA0);
const extMHI = extValid.map(g => g.logMHI);
const extMhost = extValid.map(g => g.logMhost);
const extMR = extValid.map(g => g.logMeanRun);

const rMHI_raw = partialCorr(extMHI, extA0, []);
const rMhost_raw = partialCorr(extMhost, extA0, []);
const rMR_raw = partialCorr(extMR, extA0, []);

const rMHI_ctrl = partialCorr(extMHI, extA0, [extMhost, extMR]);
const rMhost_ctrl = partialCorr(extMhost, extA0, [extMHI, extMR]);
const rMR_ctrl = partialCorr(extMR, extA0, [extMHI, extMhost]);

console.log('  ' + 'Axis'.padEnd(14) + 'Raw r'.padStart(8) + 'Partial r'.padStart(11) + ' Training sign');
console.log('  ' + '-'.repeat(45));
console.log('  ' + 'logMHI'.padEnd(14) + rMHI_raw.toFixed(3).padStart(8) + rMHI_ctrl.toFixed(3).padStart(11) + '    negative');
console.log('  ' + 'logMhost'.padEnd(14) + rMhost_raw.toFixed(3).padStart(8) + rMhost_ctrl.toFixed(3).padStart(11) + '    negative');
console.log('  ' + 'logMeanRun'.padEnd(14) + rMR_raw.toFixed(3).padStart(8) + rMR_ctrl.toFixed(3).padStart(11) + '    positive');
console.log();

const signsMatch = (rMHI_ctrl < 0 ? 1 : 0) + (rMhost_ctrl < 0 ? 1 : 0) + (rMR_ctrl > 0 ? 1 : 0);
console.log('  Sign agreement with training: ' + signsMatch + '/3');
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STEP 8: REFIT M3 ON EXTERNAL (for sign/structure comparison only)');
console.log('----------------------------------------------------------------------\n');

const extX = extValid.map(g => [g.logMHI, g.logMhost, g.logMeanRun]);
const betaExt = olsFit(obs, extX);
if (betaExt) {
  console.log('  M3 refit on external:');
  console.log('    intercept = ' + betaExt[0].toFixed(4) + '  (frozen: ' + M3_BETA[0] + ')');
  console.log('    logMHI    = ' + betaExt[1].toFixed(4) + '  (frozen: ' + M3_BETA[1] + ')');
  console.log('    logMhost  = ' + betaExt[2].toFixed(4) + '  (frozen: ' + M3_BETA[2] + ')');
  console.log('    logMeanRun= ' + betaExt[3].toFixed(4) + '  (frozen: ' + M3_BETA[3] + ')');

  const signMatch = [
    Math.sign(betaExt[1]) === Math.sign(M3_BETA[1]),
    Math.sign(betaExt[2]) === Math.sign(M3_BETA[2]),
    Math.sign(betaExt[3]) === Math.sign(M3_BETA[3]),
  ];
  console.log('    Sign agreement: ' + signMatch.filter(v => v).length + '/3 (' +
    ['logMHI', 'logMhost', 'logMeanRun'].filter((_, i) => signMatch[i]).join(', ') + ')');
}
console.log();

console.log('======================================================================');
console.log('  VERDICT');
console.log('======================================================================\n');

const m3BeatsPer = m3wins / extValid.length;
const m3BeatsRMS = rmsM3 < rmsM0;
const signsOK = signsMatch >= 2;

let verdict;
if (m3BeatsRMS && m3BeatsPer > 0.55 && gapM3 > 5 && signsOK) {
  verdict = 'PASS';
} else if (m3BeatsRMS && m3BeatsPer > 0.5) {
  verdict = 'PARTIAL-PASS';
} else if (signsOK && m3BeatsPer >= 0.45) {
  verdict = 'MARGINAL';
} else {
  verdict = 'FAIL';
}

console.log('  Frozen M3 external performance:');
console.log('    RMS M3 < RMS M0? ' + (m3BeatsRMS ? 'YES' : 'NO') + ' (' + rmsM3.toFixed(3) + ' vs ' + rmsM0.toFixed(3) + ')');
console.log('    Per-galaxy win rate: ' + (m3BeatsPer * 100).toFixed(1) + '%');
console.log('    Gap%: ' + gapM3.toFixed(1) + '%');
console.log('    Spearman: ' + rsM3.toFixed(3));
console.log('    Axis signs match training: ' + signsMatch + '/3');
console.log('    VERDICT: ' + verdict);
console.log();

if (verdict === 'PASS') {
  console.log('  The structured law M3 GENERALIZES to unseen SPARC galaxies.');
  console.log('  The multi-axis a0 variation is not an artifact of the training sample.');
} else if (verdict === 'PARTIAL-PASS') {
  console.log('  M3 shows partial generalization — directionally correct but noisy.');
  console.log('  This is consistent with: (a) the law being real but Mhost proxy adding noise,');
  console.log('  or (b) the law being weaker than measured in the training set.');
} else if (verdict === 'MARGINAL') {
  console.log('  M3 shows marginally correct structure (axis signs) but not enough');
  console.log('  predictive power on external data. May reflect measurement noise.');
} else {
  console.log('  M3 does NOT generalize to external galaxies.');
  console.log('  The multi-axis structure may be an artifact of the training sample.');
}
console.log();

const output = {
  phase: '81',
  title: 'Clean External Test — SPARC-out with CDS RC data',
  dataSource: 'CDS/VizieR J/AJ/152/157 (Lelli+2016)',
  nTotal: allGalaxies.length,
  nExternal: external.length,
  nExternalFiltered: extValid.length,
  nMhostGroupCatalog: nGroupCat,
  nMhostVflatProxy: nVflatProxy,
  a0ValidationVsTraining: {
    nMatched: nMatch,
    meanDiff: +meanDiff.toFixed(4),
    rmsDiff: +rmsDiff.toFixed(4),
  },
  meanRunValidation: {
    nMatched: nMRmatch,
    meanDiff: +meanMRdiff.toFixed(4),
    rmsDiff: +rmsMRdiff.toFixed(4),
  },
  frozenM3: {
    beta: M3_BETA,
    rmsM0: +rmsM0.toFixed(4),
    rmsM3: +rmsM3.toFixed(4),
    biasM3: +biasM3.toFixed(4),
    gapPct: +gapM3.toFixed(1),
    spearman: +rsM3.toFixed(3),
    pearson: +rpM3.toFixed(3),
    m3BeatsM0: m3BeatsRMS,
    perGalaxyWinRate: +(m3BeatsPer * 100).toFixed(1),
  },
  axisSignsExternal: {
    logMHI: +rMHI_ctrl.toFixed(3),
    logMhost: +rMhost_ctrl.toFixed(3),
    logMeanRun: +rMR_ctrl.toFixed(3),
    signsMatchTraining: signsMatch,
  },
  refitOnExternal: betaExt ? {
    intercept: +betaExt[0].toFixed(4),
    logMHI: +betaExt[1].toFixed(4),
    logMhost: +betaExt[2].toFixed(4),
    logMeanRun: +betaExt[3].toFixed(4),
  } : null,
  verdict,
  perGalaxy: extValid.map(g => ({
    name: g.name, nPts: g.nPts, logA0_obs: +g.logA0.toFixed(3),
    logMHI: +g.logMHI.toFixed(3), logMhost: +g.logMhost.toFixed(3),
    logMeanRun: +g.logMeanRun.toFixed(3),
    predM3: +(M3_BETA[0] + M3_BETA[1] * g.logMHI + M3_BETA[2] * g.logMhost + M3_BETA[3] * g.logMeanRun).toFixed(3),
    predM0: +trainMeanA0.toFixed(3),
    mhostSource: g.mhostSource,
    Q: g.Q,
  })),
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase81-external-sparc-test.json'), JSON.stringify(output, null, 2));
console.log('  Saved: public/phase81-external-sparc-test.json');

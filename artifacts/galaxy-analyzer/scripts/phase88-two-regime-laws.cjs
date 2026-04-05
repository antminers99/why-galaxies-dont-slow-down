#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 88: TWO-REGIME LAWS');
console.log('  Separate laws for Vflat < 70 and Vflat >= 70');
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

function pearsonR(x, y) {
  const n = x.length;
  if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function olsRegress(X, y) {
  const n = y.length, p = X[0].length;
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
  const r2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
  return { beta, r2, rms: Math.sqrt(ss_res / n), residuals, yhat, n, p: p };
}

function looCV(X, y) {
  const n = y.length;
  let ss = 0;
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

function kFoldCV(X, y, k) {
  const n = y.length;
  const indices = Array.from({ length: n }, (_, i) => i).sort(() => Math.random() - 0.5);
  const foldSize = Math.floor(n / k);
  let ss = 0, count = 0;
  for (let f = 0; f < k; f++) {
    const testIdx = new Set(indices.slice(f * foldSize, (f + 1) * foldSize));
    const Xt = X.filter((_, j) => !testIdx.has(j));
    const yt = y.filter((_, j) => !testIdx.has(j));
    const fit = olsRegress(Xt, yt);
    if (!fit) continue;
    for (const ti of testIdx) {
      const pred = X[ti].reduce((s, v, j) => s + v * fit.beta[j], 0);
      ss += (y[ti] - pred) ** 2;
      count++;
    }
  }
  return Math.sqrt(ss / count);
}

function permTest(X, y, nPerm) {
  const realFit = olsRegress(X, y);
  if (!realFit) return { pValue: 1 };
  const realR2 = realFit.r2;
  let nBetter = 0;
  for (let p = 0; p < nPerm; p++) {
    const shuffled = [...y].sort(() => Math.random() - 0.5);
    const fit = olsRegress(X, shuffled);
    if (fit && fit.r2 >= realR2) nBetter++;
  }
  return { pValue: nBetter / nPerm, realR2 };
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
  if (pts.length < 3) continue;
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
  const ev = envLookup[name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1));

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI: Math.log10(t1.MHI),
    logMeanRun,
    logSigma0,
    Vflat: t1.Vflat, Q: t1.Q, D: t1.D, fD: t1.fD, T: t1.T,
    L36: t1.L36, inc: t1.inc, envCode,
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    logL36: Math.log10(t1.L36),
    inTraining: trainNames.has(name),
  });
}

const loRegime = allGalaxies.filter(g => g.Vflat < VFLAT_BREAK);
const hiRegime = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);

console.log('  Total: ' + allGalaxies.length);
console.log('  Low-Vflat (<' + VFLAT_BREAK + '): N=' + loRegime.length);
console.log('  High-Vflat (>=' + VFLAT_BREAK + '): N=' + hiRegime.length);
console.log();

function testModels(data, regimeLabel) {
  console.log('══════════════════════════════════════════════════════════════');
  console.log('  ' + regimeLabel + ' (N=' + data.length + ')');
  console.log('══════════════════════════════════════════════════════════════\n');

  const y = data.map(g => g.logA0);
  const m0rms = sd(y);
  const m0mean = mean(y);
  console.log('  M0 (universal a0): mean=' + m0mean.toFixed(3) + ', sd=' + sd(y).toFixed(3) + '\n');

  const models = [
    { name: 'M0: constant', vars: [], make: () => data.map(() => [1]) },
    { name: 'Gas only: logMHI', vars: ['logMHI'], make: () => data.map(g => [1, g.logMHI]) },
    { name: 'Coherence only: logMeanRun', vars: ['logMeanRun'], make: () => data.map(g => [1, g.logMeanRun]) },
    { name: 'Density only: logSigma0', vars: ['logSigma0'], make: () => data.map(g => [1, g.logSigma0]) },
    { name: 'Env only: envCode', vars: ['envCode'], make: () => data.map(g => [1, g.envCode]) },
    { name: 'Vflat only: logVflat', vars: ['logVflat'], make: () => data.map(g => [1, g.logVflat]) },
    { name: 'Gas + Coherence', vars: ['logMHI', 'logMeanRun'], make: () => data.map(g => [1, g.logMHI, g.logMeanRun]) },
    { name: 'Gas + Density', vars: ['logMHI', 'logSigma0'], make: () => data.map(g => [1, g.logMHI, g.logSigma0]) },
    { name: 'Gas + Env', vars: ['logMHI', 'envCode'], make: () => data.map(g => [1, g.logMHI, g.envCode]) },
    { name: 'Gas + Coherence + Density', vars: ['logMHI', 'logMeanRun', 'logSigma0'], make: () => data.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0]) },
    { name: 'Gas + Coherence + Env', vars: ['logMHI', 'logMeanRun', 'envCode'], make: () => data.map(g => [1, g.logMHI, g.logMeanRun, g.envCode]) },
    { name: 'M3-like: Gas + Coherence + Density + Env', vars: ['logMHI', 'logMeanRun', 'logSigma0', 'envCode'], make: () => data.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0, g.envCode]) },
    { name: 'Full: Gas + Coh + Dens + Env + Vflat', vars: ['logMHI', 'logMeanRun', 'logSigma0', 'envCode', 'logVflat'], make: () => data.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0, g.envCode, g.logVflat]) },
  ];

  const results = [];

  console.log('  ' + 'Model'.padEnd(45) + 'R2'.padStart(7) + 'RMS'.padStart(8) +
    'LOO'.padStart(8) + '5-fold'.padStart(8) + 'LOO%gap'.padStart(9));

  for (const m of models) {
    const X = m.make();
    const fit = olsRegress(X, y);
    if (!fit) {
      console.log('  ' + m.name.padEnd(45) + ' FAILED');
      continue;
    }

    const loo = m.vars.length > 0 ? looCV(X, y) : m0rms;
    const folds = [];
    if (m.vars.length > 0 && data.length >= 15) {
      for (let rep = 0; rep < 20; rep++) folds.push(kFoldCV(X, y, 5));
    }
    const kfold5 = folds.length > 0 ? mean(folds) : m0rms;
    const looGap = m0rms > 0 ? ((m0rms - loo) / m0rms * 100) : 0;

    const result = {
      name: m.name, vars: m.vars,
      r2: +fit.r2.toFixed(4), rms: +fit.rms.toFixed(4),
      loo: +loo.toFixed(4), kfold5: +kfold5.toFixed(4),
      looGapPct: +looGap.toFixed(1),
      beta: fit.beta.map(b => +b.toFixed(4)),
    };
    results.push(result);

    console.log('  ' + m.name.padEnd(45) +
      fit.r2.toFixed(4).padStart(7) + fit.rms.toFixed(4).padStart(8) +
      loo.toFixed(4).padStart(8) + kfold5.toFixed(4).padStart(8) +
      looGap.toFixed(1).padStart(9) + '%');
  }

  const best = results.filter(r => r.vars.length > 0).sort((a, b) => a.loo - b.loo)[0];
  console.log('\n  Best by LOO: ' + best.name + ' (LOO=' + best.loo + ', gap=' + best.looGapPct + '%)');
  console.log('  Coefficients: ' + best.beta.join(', '));
  console.log();

  console.log('  ── Correlation matrix ──');
  const varFns = {
    logMHI: g => g.logMHI, logMeanRun: g => g.logMeanRun,
    logSigma0: g => g.logSigma0, envCode: g => g.envCode,
    logVflat: g => g.logVflat, logA0: g => g.logA0,
  };
  const varNames = Object.keys(varFns);
  console.log('  ' + ''.padEnd(12) + varNames.map(v => v.padStart(10)).join(''));
  for (const v1 of varNames) {
    let line = '  ' + v1.padEnd(12);
    for (const v2 of varNames) {
      const r = pearsonR(data.map(varFns[v1]), data.map(varFns[v2]));
      line += (isNaN(r) ? '   ---' : r.toFixed(3)).padStart(10);
    }
    console.log(line);
  }
  console.log();

  if (best.vars.length > 0 && data.length >= 20) {
    console.log('  ── Permutation test for best model (1000 perms) ──');
    const X = models.find(m => m.name === best.name).make();
    const perm = permTest(X, y, 1000);
    console.log('    p-value: ' + perm.pValue.toFixed(3));
    console.log('    R2: ' + perm.realR2.toFixed(4));
    best.permPvalue = perm.pValue;
  }

  return { results, best, m0rms: +m0rms.toFixed(4), m0mean: +m0mean.toFixed(3), n: data.length };
}

const loResults = testModels(loRegime, 'LOW-Vflat REGIME (Vflat < ' + VFLAT_BREAK + ')');
const hiResults = testModels(hiRegime, 'HIGH-Vflat REGIME (Vflat >= ' + VFLAT_BREAK + ')');

console.log('══════════════════════════════════════════════════════════════');
console.log('  88d: SIGN COMPARISON ACROSS REGIMES');
console.log('══════════════════════════════════════════════════════════════\n');

const signVars = ['logMHI', 'logMeanRun', 'logSigma0', 'envCode', 'logVflat'];
const signComparison = {};

console.log('  ' + 'Variable'.padEnd(15) + 'r(lo)'.padStart(8) + 'r(hi)'.padStart(8) +
  'sign_lo'.padStart(9) + 'sign_hi'.padStart(9) + 'SAME?'.padStart(8));

for (const v of signVars) {
  const fn = g => g[v];
  const rLo = pearsonR(loRegime.map(fn), loRegime.map(g => g.logA0));
  const rHi = pearsonR(hiRegime.map(fn), hiRegime.map(g => g.logA0));
  const signLo = isNaN(rLo) ? '?' : rLo > 0 ? '+' : '-';
  const signHi = isNaN(rHi) ? '?' : rHi > 0 ? '+' : '-';
  const same = signLo === signHi ? 'YES' : '** NO **';
  signComparison[v] = { rLo: +(rLo || 0).toFixed(3), rHi: +(rHi || 0).toFixed(3), signLo, signHi, same: signLo === signHi };

  console.log('  ' + v.padEnd(15) +
    (isNaN(rLo) ? '---' : rLo.toFixed(3)).padStart(8) +
    (isNaN(rHi) ? '---' : rHi.toFixed(3)).padStart(8) +
    signLo.padStart(9) + signHi.padStart(9) + same.padStart(8));
}

const nReversed = Object.values(signComparison).filter(s => !s.same).length;
console.log('\n  Variables with reversed signs: ' + nReversed + '/' + signVars.length);
if (nReversed === 1) {
  const which = Object.entries(signComparison).find(([_, s]) => !s.same);
  console.log('  Only ' + which[0] + ' reverses → GAS-REGIME TRANSITION');
} else if (nReversed >= 3) {
  console.log('  Multiple reversals → FULL STATE TRANSITION');
} else {
  console.log('  ' + nReversed + ' reversals → PARTIAL TRANSITION');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  88e: DOES THE TRANSITION AFFECT ALL AXES?');
console.log('  Testing each variable for broken behavior at Vflat=' + VFLAT_BREAK);
console.log('══════════════════════════════════════════════════════════════\n');

const axisTransition = {};
for (const v of signVars) {
  const fn = g => g[v];
  const X_single = allGalaxies.map(g => [1, fn(g)]);
  const indicator = allGalaxies.map(g => g.Vflat >= VFLAT_BREAK ? 1 : 0);
  const X_broken = allGalaxies.map((g, i) => [1, fn(g), indicator[i], fn(g) * indicator[i]]);
  const y = allGalaxies.map(g => g.logA0);

  const fitSingle = olsRegress(X_single, y);
  const fitBroken = olsRegress(X_broken, y);
  if (!fitSingle || !fitBroken) continue;

  const sseSingle = fitSingle.residuals.reduce((s, r) => s + r * r, 0);
  const sseBroken = fitBroken.residuals.reduce((s, r) => s + r * r, 0);
  const n = allGalaxies.length;
  const bicSingle = n * Math.log(sseSingle / n) + 2 * Math.log(n);
  const bicBroken = n * Math.log(sseBroken / n) + 4 * Math.log(n);
  const deltaBIC = bicBroken - bicSingle;

  const slopeLo = fitBroken.beta[1];
  const slopeHi = fitBroken.beta[1] + fitBroken.beta[3];

  axisTransition[v] = {
    slopeLo: +slopeLo.toFixed(4), slopeHi: +slopeHi.toFixed(4),
    deltaBIC: +deltaBIC.toFixed(2),
    reverses: Math.sign(slopeLo) !== Math.sign(slopeHi),
  };

  const status = deltaBIC < -6 ? 'STRONG' : deltaBIC < -2 ? 'MODERATE' : 'WEAK/NONE';
  const rev = Math.sign(slopeLo) !== Math.sign(slopeHi) ? 'REVERSES' : 'same sign';

  console.log('  ' + v.padEnd(15) + 'slope_lo=' + slopeLo.toFixed(4).padStart(8) +
    ', slope_hi=' + slopeHi.toFixed(4).padStart(8) +
    ', ΔBIC=' + deltaBIC.toFixed(2).padStart(7) +
    ' [' + status + '] ' + rev);
}

const nAxesReverse = Object.values(axisTransition).filter(a => a.reverses).length;
const nAxesStrong = Object.values(axisTransition).filter(a => a.deltaBIC < -6).length;
console.log('\n  Axes that reverse: ' + nAxesReverse + '/' + signVars.length);
console.log('  Axes with strong broken model: ' + nAxesStrong + '/' + signVars.length);

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  88f: CROSS-REGIME VALIDATION');
console.log('  Train on one regime, test on the other');
console.log('══════════════════════════════════════════════════════════════\n');

function crossRegimeTest(trainData, testData, trainLabel, testLabel) {
  console.log('  Train on ' + trainLabel + ' (N=' + trainData.length + '), test on ' + testLabel + ' (N=' + testData.length + ')');

  const bestTrainModel = trainLabel.includes('HIGH') ? hiResults.best : loResults.best;
  const vars = bestTrainModel.vars;
  if (vars.length === 0) { console.log('    No predictive model found\n'); return null; }

  const makeX = (data) => data.map(g => {
    const row = [1];
    for (const v of vars) row.push(g[v]);
    return row;
  });

  const Xtrain = makeX(trainData);
  const ytrain = trainData.map(g => g.logA0);
  const fit = olsRegress(Xtrain, ytrain);
  if (!fit) { console.log('    Fit failed\n'); return null; }

  const Xtest = makeX(testData);
  const ytest = testData.map(g => g.logA0);
  const preds = Xtest.map(row => row.reduce((s, v, j) => s + v * fit.beta[j], 0));
  const testResid = ytest.map((v, i) => v - preds[i]);
  const testRMS = Math.sqrt(testResid.reduce((s, r) => s + r * r, 0) / ytest.length);
  const m0RMS = sd(ytest);
  const gap = m0RMS > 0 ? (m0RMS - testRMS) / m0RMS * 100 : 0;

  console.log('    Model: ' + bestTrainModel.name);
  console.log('    Coefficients: ' + fit.beta.map(b => b.toFixed(4)).join(', '));
  console.log('    Train LOO RMS: ' + bestTrainModel.loo);
  console.log('    Test RMS: ' + testRMS.toFixed(4) + ' (M0 baseline: ' + m0RMS.toFixed(4) + ')');
  console.log('    Gap closed on test: ' + gap.toFixed(1) + '%');
  console.log('    ' + (gap > 0 ? 'TRANSFERS' : '** FAILS TO TRANSFER **'));
  console.log();

  return { model: bestTrainModel.name, trainLOO: bestTrainModel.loo, testRMS: +testRMS.toFixed(4), testM0: +m0RMS.toFixed(4), gap: +gap.toFixed(1) };
}

const hiToLo = crossRegimeTest(hiRegime, loRegime, 'HIGH-Vflat', 'LOW-Vflat');
const loToHi = crossRegimeTest(loRegime, hiRegime, 'LOW-Vflat', 'HIGH-Vflat');

console.log('══════════════════════════════════════════════════════════════');
console.log('  FINAL SUMMARY');
console.log('══════════════════════════════════════════════════════════════\n');

console.log('  Transition fixed at: Vflat = ' + VFLAT_BREAK + ' km/s\n');

console.log('  Low-Vflat regime:');
console.log('    N = ' + loResults.n);
console.log('    Best model = ' + loResults.best.name);
console.log('    LOO RMS = ' + loResults.best.loo);
console.log('    LOO gap = ' + loResults.best.looGapPct + '%');
console.log('    Coefficients = ' + loResults.best.beta.join(', '));
console.log();

console.log('  High-Vflat regime:');
console.log('    N = ' + hiResults.n);
console.log('    Best model = ' + hiResults.best.name);
console.log('    LOO RMS = ' + hiResults.best.loo);
console.log('    LOO gap = ' + hiResults.best.looGapPct + '%');
console.log('    Coefficients = ' + hiResults.best.beta.join(', '));
console.log();

console.log('  Sign comparison:');
for (const [v, s] of Object.entries(signComparison)) {
  console.log('    ' + v + ': lo=' + s.signLo + s.rLo + ' / hi=' + s.signHi + s.rHi +
    (s.same ? '' : ' ** REVERSED **'));
}
console.log();

const sameStructure = nReversed <= 1 && loResults.best.name === hiResults.best.name;
console.log('  Do both regimes share one law? ' + (sameStructure ? 'PARTIALLY' : 'NO'));
console.log('  Require two separate laws? ' + (sameStructure ? 'MAYBE — same structure, different coefficients' : 'YES'));
console.log();

let verdict;
if (nReversed >= 2 && loResults.best.looGapPct > 5 && hiResults.best.looGapPct > 5) {
  verdict = 'CONFIRMED-TWO-REGIME';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  CONFIRMED-TWO-REGIME: Different laws for different regimes ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else if (nReversed >= 1 && (loResults.best.looGapPct > 5 || hiResults.best.looGapPct > 5)) {
  verdict = 'PARTIAL';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  PARTIAL: One regime has structured law, other is weaker    ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else {
  verdict = 'FAIL';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  FAIL: No strong structured law in either regime            ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
}

const output = {
  phase: '88',
  title: 'Two-Regime Laws',
  vflatBreak: VFLAT_BREAK,
  lowRegime: loResults,
  highRegime: hiResults,
  signComparison,
  axisTransition,
  crossValidation: { hiToLo, loToHi },
  nSignReversals: nReversed,
  nAxesReverseInBroken: nAxesReverse,
  nAxesStrongBroken: nAxesStrong,
  verdict,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase88-two-regime-laws.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase88-two-regime-laws.json');

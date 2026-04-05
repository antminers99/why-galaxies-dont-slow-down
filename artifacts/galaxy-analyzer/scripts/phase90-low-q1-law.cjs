#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 90: LOW-REGIME Q=1 LAW TEST');
console.log('  Can we extract a structured law from Q=1 dwarfs?');
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
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, rms: Math.sqrt(ss_res / n), residuals, yhat };
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

function permTest(X, y, nPerm) {
  const realFit = olsRegress(X, y);
  if (!realFit) return { pValue: 1, realR2: 0 };
  let nBetter = 0;
  for (let p = 0; p < nPerm; p++) {
    const shuffled = [...y].sort(() => Math.random() - 0.5);
    const fit = olsRegress(X, shuffled);
    if (fit && fit.r2 >= realFit.r2) nBetter++;
  }
  return { pValue: nBetter / nPerm, realR2: realFit.r2 };
}

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
  const logGbarRange = Math.max(...pts.map(p => p.logGbar)) - Math.min(...pts.map(p => p.logGbar));

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI: Math.log10(t1.MHI),
    logMeanRun, logSigma0,
    Vflat: t1.Vflat, Q: t1.Q, T: t1.T,
    envCode, logGbarRange,
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    logL36: Math.log10(t1.L36),
    inc: t1.inc,
  });
}

const loQ1 = allGalaxies.filter(g => g.Vflat < VFLAT_BREAK && g.Q === 1);
const hiQ1 = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK && g.Q === 1);
const loAll = allGalaxies.filter(g => g.Vflat < VFLAT_BREAK);
const hiAll = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);

console.log('  Samples:');
console.log('    Lo Q=1: N=' + loQ1.length + ', sd(logA0)=' + sd(loQ1.map(g => g.logA0)).toFixed(3));
console.log('    Hi Q=1: N=' + hiQ1.length + ', sd(logA0)=' + sd(hiQ1.map(g => g.logA0)).toFixed(3));
console.log('    Lo all: N=' + loAll.length + ', sd(logA0)=' + sd(loAll.map(g => g.logA0)).toFixed(3));
console.log('    Hi all: N=' + hiAll.length + ', sd(logA0)=' + sd(hiAll.map(g => g.logA0)).toFixed(3));
console.log();

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: CORRELATION STRUCTURE IN LOW Q=1');
console.log('══════════════════════════════════════════════════════════════\n');

const corrVars = {
  logMHI: g => g.logMHI, logMeanRun: g => g.logMeanRun,
  logSigma0: g => g.logSigma0, logVflat: g => g.logVflat,
  logGbarRange: g => g.logGbarRange, nPts: g => g.nPts,
  inc: g => g.inc, T: g => g.T,
};

console.log('  r(var, logA0) in each subsample:\n');
console.log('  ' + 'Variable'.padEnd(15) + 'LoQ1'.padStart(8) + 'LoAll'.padStart(8) +
  'HiQ1'.padStart(8) + 'HiAll'.padStart(8));

const corrResults = {};
for (const [v, fn] of Object.entries(corrVars)) {
  const rLoQ1 = pearsonR(loQ1.map(fn), loQ1.map(g => g.logA0));
  const rLoAll = pearsonR(loAll.map(fn), loAll.map(g => g.logA0));
  const rHiQ1 = pearsonR(hiQ1.map(fn), hiQ1.map(g => g.logA0));
  const rHiAll = pearsonR(hiAll.map(fn), hiAll.map(g => g.logA0));
  corrResults[v] = { loQ1: +rLoQ1.toFixed(3), loAll: +rLoAll.toFixed(3), hiQ1: +rHiQ1.toFixed(3), hiAll: +rHiAll.toFixed(3) };
  console.log('  ' + v.padEnd(15) +
    (isNaN(rLoQ1) ? '---' : rLoQ1.toFixed(3)).padStart(8) +
    (isNaN(rLoAll) ? '---' : rLoAll.toFixed(3)).padStart(8) +
    (isNaN(rHiQ1) ? '---' : rHiQ1.toFixed(3)).padStart(8) +
    (isNaN(rHiAll) ? '---' : rHiAll.toFixed(3)).padStart(8));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: MODEL COMPETITION IN LOW Q=1 (N=' + loQ1.length + ')');
console.log('══════════════════════════════════════════════════════════════\n');

function testModelsOn(data, label) {
  const y = data.map(g => g.logA0);
  const m0sd = sd(y);

  const models = [
    { name: 'Gas: logMHI', make: () => data.map(g => [1, g.logMHI]) },
    { name: 'Coherence: logMeanRun', make: () => data.map(g => [1, g.logMeanRun]) },
    { name: 'Density: logSigma0', make: () => data.map(g => [1, g.logSigma0]) },
    { name: 'Vflat: logVflat', make: () => data.map(g => [1, g.logVflat]) },
    { name: 'gRange: logGbarRange', make: () => data.map(g => [1, g.logGbarRange]) },
    { name: 'nPts', make: () => data.map(g => [1, g.nPts]) },
    { name: 'Morph: T', make: () => data.map(g => [1, g.T]) },
    { name: 'inc', make: () => data.map(g => [1, g.inc]) },
    { name: 'Gas + Coh', make: () => data.map(g => [1, g.logMHI, g.logMeanRun]) },
    { name: 'Gas + Dens', make: () => data.map(g => [1, g.logMHI, g.logSigma0]) },
    { name: 'Gas + Vflat', make: () => data.map(g => [1, g.logMHI, g.logVflat]) },
    { name: 'Coh + Dens', make: () => data.map(g => [1, g.logMeanRun, g.logSigma0]) },
    { name: 'Gas + Coh + Dens', make: () => data.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0]) },
    { name: 'Gas + Coh + Vflat', make: () => data.map(g => [1, g.logMHI, g.logMeanRun, g.logVflat]) },
  ];

  const results = [];
  console.log('  ' + label + ' (N=' + data.length + ', M0 sd=' + m0sd.toFixed(3) + ')\n');
  console.log('  ' + 'Model'.padEnd(25) + 'R2'.padStart(7) + 'RMS'.padStart(8) +
    'LOO'.padStart(8) + 'Gap%'.padStart(7));

  for (const m of models) {
    const X = m.make();
    const fit = olsRegress(X, y);
    if (!fit) { console.log('  ' + m.name.padEnd(25) + ' FAILED'); continue; }
    const loo = looCV(X, y);
    const gap = m0sd > 0 ? (m0sd - loo) / m0sd * 100 : 0;
    results.push({ name: m.name, r2: +fit.r2.toFixed(4), rms: +fit.rms.toFixed(4), loo: +loo.toFixed(4), gap: +gap.toFixed(1), beta: fit.beta.map(b => +b.toFixed(4)) });
    console.log('  ' + m.name.padEnd(25) + fit.r2.toFixed(4).padStart(7) +
      fit.rms.toFixed(4).padStart(8) + loo.toFixed(4).padStart(8) + gap.toFixed(1).padStart(7) + '%');
  }

  const best = results.filter(r => r.gap > 0).sort((a, b) => b.gap - a.gap);
  if (best.length > 0) {
    console.log('\n  Best: ' + best[0].name + ' (LOO gap=' + best[0].gap + '%, coefficients=' + best[0].beta.join(', ') + ')');
  } else {
    console.log('\n  No model beats M0 on LOO.');
  }

  return { results, m0sd: +m0sd.toFixed(3), best: best.length > 0 ? best[0] : null };
}

const loQ1models = testModelsOn(loQ1, 'Low-Vflat Q=1');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: SAME MODELS ON HIGH Q=1 (for comparison)');
console.log('══════════════════════════════════════════════════════════════\n');

const hiQ1models = testModelsOn(hiQ1, 'High-Vflat Q=1');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: PERMUTATION TESTS FOR TOP MODELS IN LOW Q=1');
console.log('══════════════════════════════════════════════════════════════\n');

const topLoModels = loQ1models.results.filter(r => r.gap > 0).sort((a, b) => b.gap - a.gap).slice(0, 5);
const permResults = [];

for (const m of topLoModels) {
  const modelDef = [
    { name: 'Gas: logMHI', make: () => loQ1.map(g => [1, g.logMHI]) },
    { name: 'Coherence: logMeanRun', make: () => loQ1.map(g => [1, g.logMeanRun]) },
    { name: 'Density: logSigma0', make: () => loQ1.map(g => [1, g.logSigma0]) },
    { name: 'Vflat: logVflat', make: () => loQ1.map(g => [1, g.logVflat]) },
    { name: 'gRange: logGbarRange', make: () => loQ1.map(g => [1, g.logGbarRange]) },
    { name: 'nPts', make: () => loQ1.map(g => [1, g.nPts]) },
    { name: 'Morph: T', make: () => loQ1.map(g => [1, g.T]) },
    { name: 'inc', make: () => loQ1.map(g => [1, g.inc]) },
    { name: 'Gas + Coh', make: () => loQ1.map(g => [1, g.logMHI, g.logMeanRun]) },
    { name: 'Gas + Dens', make: () => loQ1.map(g => [1, g.logMHI, g.logSigma0]) },
    { name: 'Gas + Vflat', make: () => loQ1.map(g => [1, g.logMHI, g.logVflat]) },
    { name: 'Coh + Dens', make: () => loQ1.map(g => [1, g.logMeanRun, g.logSigma0]) },
    { name: 'Gas + Coh + Dens', make: () => loQ1.map(g => [1, g.logMHI, g.logMeanRun, g.logSigma0]) },
    { name: 'Gas + Coh + Vflat', make: () => loQ1.map(g => [1, g.logMHI, g.logMeanRun, g.logVflat]) },
  ].find(md => md.name === m.name);

  if (!modelDef) continue;
  const X = modelDef.make();
  const y = loQ1.map(g => g.logA0);
  const perm = permTest(X, y, 2000);
  permResults.push({ name: m.name, gap: m.gap, r2: m.r2, pValue: perm.pValue });
  console.log('  ' + m.name + ': R2=' + m.r2 + ', LOO gap=' + m.gap + '%, perm p=' + perm.pValue.toFixed(3));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: SIGN COMPARISON — Do Q=1 dwarfs share signs');
console.log('  with high-Vflat regime?');
console.log('══════════════════════════════════════════════════════════════\n');

const signCheck = {};
const checkVars = ['logMHI', 'logMeanRun', 'logSigma0', 'logVflat'];
console.log('  ' + 'Variable'.padEnd(15) + 'LoQ1'.padStart(8) + 'HiQ1'.padStart(8) + 'HiAll'.padStart(8) + 'Same?'.padStart(8));

for (const v of checkVars) {
  const fn = g => g[v];
  const rLoQ1 = pearsonR(loQ1.map(fn), loQ1.map(g => g.logA0));
  const rHiQ1 = pearsonR(hiQ1.map(fn), hiQ1.map(g => g.logA0));
  const rHiAll = pearsonR(hiAll.map(fn), hiAll.map(g => g.logA0));
  const sameAsHi = (Math.sign(rLoQ1) === Math.sign(rHiAll)) || isNaN(rLoQ1);
  signCheck[v] = { loQ1: +(rLoQ1 || 0).toFixed(3), hiQ1: +(rHiQ1 || 0).toFixed(3), hiAll: +(rHiAll || 0).toFixed(3), same: sameAsHi };
  console.log('  ' + v.padEnd(15) +
    (isNaN(rLoQ1) ? '---' : rLoQ1.toFixed(3)).padStart(8) +
    (isNaN(rHiQ1) ? '---' : rHiQ1.toFixed(3)).padStart(8) +
    rHiAll.toFixed(3).padStart(8) +
    (sameAsHi ? 'YES' : '** NO **').padStart(8));
}

const nSame = Object.values(signCheck).filter(s => s.same).length;
console.log('\n  Signs matching high regime: ' + nSame + '/' + checkVars.length);

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 6: CROSS-REGIME TRANSFER WITH Q=1 ONLY');
console.log('  Train on HiQ1, test on LoQ1 (and vice versa)');
console.log('══════════════════════════════════════════════════════════════\n');

function crossTest(trainData, testData, vars, trainLabel, testLabel) {
  const makeX = d => d.map(g => {
    const row = [1];
    for (const v of vars) row.push(g[v]);
    return row;
  });
  const Xt = makeX(trainData);
  const yt = trainData.map(g => g.logA0);
  const fit = olsRegress(Xt, yt);
  if (!fit) return null;
  const Xtest = makeX(testData);
  const ytest = testData.map(g => g.logA0);
  const preds = Xtest.map(r => r.reduce((s, v, j) => s + v * fit.beta[j], 0));
  const testRMS = Math.sqrt(preds.reduce((s, p, i) => s + (ytest[i] - p) ** 2, 0) / ytest.length);
  const m0RMS = sd(ytest);
  const gap = m0RMS > 0 ? (m0RMS - testRMS) / m0RMS * 100 : 0;
  console.log('  ' + trainLabel + ' → ' + testLabel + ': testRMS=' + testRMS.toFixed(3) +
    ', M0=' + m0RMS.toFixed(3) + ', gap=' + gap.toFixed(1) + '% ' + (gap > 0 ? 'TRANSFERS' : 'FAILS'));
  return { testRMS: +testRMS.toFixed(3), m0: +m0RMS.toFixed(3), gap: +gap.toFixed(1) };
}

const xferVars = ['logMHI', 'logMeanRun'];
console.log('  Using model: Gas + Coherence');
const hiToLo = crossTest(hiQ1, loQ1, xferVars, 'HiQ1', 'LoQ1');
const loToHi = crossTest(loQ1, hiQ1, xferVars, 'LoQ1', 'HiQ1');

const xferVars2 = ['logMHI'];
console.log('\n  Using model: Gas only');
const hiToLo2 = crossTest(hiQ1, loQ1, xferVars2, 'HiQ1', 'LoQ1');
const loToHi2 = crossTest(loQ1, hiQ1, xferVars2, 'LoQ1', 'HiQ1');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 7: ALL SPARC Q=1 — UNIFIED MODEL (no Vflat split)');
console.log('══════════════════════════════════════════════════════════════\n');

const allQ1 = allGalaxies.filter(g => g.Q === 1);
console.log('  All Q=1: N=' + allQ1.length + ', sd(logA0)=' + sd(allQ1.map(g => g.logA0)).toFixed(3));
const allQ1models = testModelsOn(allQ1, 'All SPARC Q=1');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 8: GALAXY-BY-GALAXY IN LOW Q=1');
console.log('══════════════════════════════════════════════════════════════\n');

console.log('  ' + 'Name'.padEnd(14) + 'Vflat'.padStart(6) + 'nPts'.padStart(6) +
  'logA0'.padStart(8) + 'logMHI'.padStart(8) + 'gRange'.padStart(8) + 'T'.padStart(4));
for (const g of loQ1.sort((a, b) => a.logA0 - b.logA0)) {
  console.log('  ' + g.name.padEnd(14) + g.Vflat.toFixed(0).padStart(6) +
    String(g.nPts).padStart(6) + g.logA0.toFixed(3).padStart(8) +
    g.logMHI.toFixed(2).padStart(8) + g.logGbarRange.toFixed(2).padStart(8) +
    String(g.T).padStart(4));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  FINAL SUMMARY');
console.log('══════════════════════════════════════════════════════════════\n');

const hasLaw = loQ1models.best && loQ1models.best.gap > 5;
const hasMarginalLaw = loQ1models.best && loQ1models.best.gap > 0;
const signsMatch = nSame >= 3;

console.log('  Low Q=1 (N=' + loQ1.length + '):');
console.log('    Best model: ' + (loQ1models.best ? loQ1models.best.name + ' (gap=' + loQ1models.best.gap + '%)' : 'NONE'));
console.log('    Signs match high regime: ' + nSame + '/' + checkVars.length);
console.log();

let verdict;
if (hasLaw && signsMatch) {
  verdict = 'SAME-LAW-BOTH-REGIMES';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  SAME LAW IN BOTH REGIMES (quality-masked in dwarfs)       ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else if (hasLaw && !signsMatch) {
  verdict = 'DIFFERENT-LAWS';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  DIFFERENT LAWS: Both regimes have structure but different  ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else if (hasMarginalLaw) {
  verdict = 'MARGINAL-SIGNAL';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  MARGINAL: Weak signal in low Q=1 (N too small?)           ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else {
  verdict = 'NO-LAW-IN-DWARFS';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  NO LAW IN DWARFS: Even Q=1 shows no structured a0        ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
}

const output = {
  phase: '90',
  title: 'Low-Regime Q=1 Law Test',
  samples: { loQ1: loQ1.length, hiQ1: hiQ1.length, loAll: loAll.length, hiAll: hiAll.length, allQ1: allQ1.length },
  correlations: corrResults,
  loQ1models,
  hiQ1models,
  permResults,
  signComparison: signCheck,
  nSignsMatch: nSame,
  crossRegime: { hiToLo, loToHi, hiToLo2, loToHi2 },
  allQ1models,
  verdict,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase90-low-q1-law.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase90-low-q1-law.json');

#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 125: ENVIRONMENTAL THRESHOLD LAW');
console.log('');
console.log('  Does M3 improve if we replace linear logMhost with a');
console.log('  threshold/broken-line representation at Mhost=11.5?');
console.log('');
console.log('  Models:');
console.log('    A) M3-linear:  MHI + Mhost + MeanRun');
console.log('    B) M3-broken:  MHI + Mhost_low + Mhost_high + MeanRun');
console.log('    C) M3-step:    MHI + I(Mhost>=11.5) + MeanRun');
console.log('    D) M3-hybrid:  MHI + I(Mhost>=11.5) + slope_above + MeanRun');
console.log('');
console.log('  Plus: source sensitivity (KT2017 vs non-KT2017)');
console.log('  Plus: does logSigma0 residual signal survive after threshold fix?');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'stage-A-master-table.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparc.forEach(s => { sparcMap[s.name] = s; });
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
const N = gals45.length;

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function pearsonR(x, y) {
  const n = x.length, mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function solveLinear(A, b) {
  const n = A.length;
  const M = A.map((r, i) => [...r, b[i]]);
  for (let i = 0; i < n; i++) {
    let mx = i;
    for (let j = i + 1; j < n; j++) if (Math.abs(M[j][i]) > Math.abs(M[mx][i])) mx = j;
    [M[i], M[mx]] = [M[mx], M[i]];
    if (Math.abs(M[i][i]) < 1e-15) continue;
    for (let j = i + 1; j < n; j++) {
      const f = M[j][i] / M[i][i];
      for (let k = i; k <= n; k++) M[j][k] -= f * M[i][k];
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}

function ols(Y, X) {
  const n = Y.length, p = X[0].length + 1;
  const Xa = X.map(r => [1, ...r]);
  const XtX = Array.from({ length: p }, () => new Array(p).fill(0));
  const XtY = new Array(p).fill(0);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < p; j++) {
      XtY[j] += Xa[i][j] * Y[i];
      for (let l = 0; l < p; l++) XtX[j][l] += Xa[i][j] * Xa[i][l];
    }
  const beta = solveLinear(XtX, XtY);
  const resid = Y.map((y, i) => y - Xa[i].reduce((s, x, j) => s + x * beta[j], 0));
  const rss = resid.reduce((s, r) => s + r * r, 0);
  const tss = Y.reduce((s, y) => s + (y - mean(Y)) ** 2, 0);
  const r2 = tss > 0 ? 1 - rss / tss : 0;
  const r2adj = 1 - (rss / (n - p)) / (tss / (n - 1));
  const aic = n * Math.log(rss / n) + 2 * p;
  const bic = n * Math.log(rss / n) + p * Math.log(n);
  return { beta, resid, rss, tss, r2, r2adj, aic, bic, n, k: p };
}

function looCV(Y, X) {
  const n = Y.length; let ss = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xt = [...X.slice(0, i), ...X.slice(i + 1)];
    const f = ols(Yt, Xt);
    const xi = [1, ...X[i]];
    ss += (Y[i] - xi.reduce((s, x, j) => s + x * f.beta[j], 0)) ** 2;
  }
  return Math.sqrt(ss / n);
}

function gapPct(rms, sdy) { return 100 * (1 - rms ** 2 / sdy ** 2); }

function mulberry32(a) {
  return function () {
    a |= 0; a = a + 0x6D2B79F5 | 0;
    var t = Math.imul(a ^ a >>> 15, 1 | a);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

function bootstrapSignStability(Y, X, nBoot, seed) {
  const rng = mulberry32(seed);
  const n = Y.length, p = X[0].length;
  const baseFit = ols(Y, X);
  const baseSigns = baseFit.beta.slice(1).map(b => Math.sign(b));
  const flips = new Array(p).fill(0);
  const coeffSamples = Array.from({ length: p }, () => []);
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(rng() * n));
    const Yb = idx.map(i => Y[i]); const Xb = idx.map(i => X[i]);
    const fit = ols(Yb, Xb);
    for (let j = 0; j < p; j++) {
      if (Math.sign(fit.beta[j + 1]) !== baseSigns[j]) flips[j]++;
      coeffSamples[j].push(fit.beta[j + 1]);
    }
  }
  return {
    flipRates: flips.map(f => f / nBoot),
    medians: coeffSamples.map(s => { s.sort((a, b) => a - b); return s[Math.floor(s.length / 2)]; }),
    ci025: coeffSamples.map(s => s[Math.floor(s.length * 0.025)]),
    ci975: coeffSamples.map(s => s[Math.floor(s.length * 0.975)])
  };
}

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logMHI = gals45.map(g => g.logMHI);
const logMR = gals45.map(g => g.logMeanRun);
const logSig0 = gals45.map(g => g.logSigma0);
const sources = gals45.map(g => tdMap[g.name].source);

const THRESHOLD = 11.5;

const mhostLow = logMhost.map(v => Math.min(v, THRESHOLD));
const mhostHigh = logMhost.map(v => Math.max(v - THRESHOLD, 0));
const mhostStep = logMhost.map(v => v >= THRESHOLD ? 1 : 0);
const mhostSlopeAbove = logMhost.map(v => v >= THRESHOLD ? (v - THRESHOLD) : 0);

const models = {
  A: {
    name: 'M3-linear (MHI + Mhost + MR)',
    X: gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]),
    varNames: ['logMHI', 'logMhost', 'logMR']
  },
  B: {
    name: 'M3-broken (MHI + Mhost_low + Mhost_high + MR)',
    X: gals45.map((_, i) => [logMHI[i], mhostLow[i], mhostHigh[i], logMR[i]]),
    varNames: ['logMHI', 'Mhost_low', 'Mhost_high', 'logMR']
  },
  C: {
    name: 'M3-step (MHI + I(Mhost>=11.5) + MR)',
    X: gals45.map((_, i) => [logMHI[i], mhostStep[i], logMR[i]]),
    varNames: ['logMHI', 'I(Mhost>=11.5)', 'logMR']
  },
  D: {
    name: 'M3-hybrid (MHI + I(Mhost>=11.5) + slope_above + MR)',
    X: gals45.map((_, i) => [logMHI[i], mhostStep[i], mhostSlopeAbove[i], logMR[i]]),
    varNames: ['logMHI', 'I(Mhost>=11.5)', 'slope_above', 'logMR']
  }
};

console.log('  N = ' + N + ', SD(logA0) = ' + sdY.toFixed(4) + ' dex');
console.log('  Threshold = ' + THRESHOLD);
console.log('  N below threshold: ' + mhostStep.filter(v => v === 0).length);
console.log('  N at/above threshold: ' + mhostStep.filter(v => v === 1).length);
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 1: HEAD-TO-HEAD COMPARISON');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const NBOOT = 2000;
const modelResults = {};

for (const [key, mod] of Object.entries(models)) {
  const fit = ols(Y, mod.X);
  const loo = looCV(Y, mod.X);
  const gap = gapPct(loo, sdY);
  const boot = bootstrapSignStability(Y, mod.X, NBOOT, 42);

  modelResults[key] = { fit, loo, gap, boot, name: mod.name, varNames: mod.varNames };

  console.log('  Model ' + key + ': ' + mod.name);
  console.log('    Coefficients:');
  console.log('      intercept = ' + fit.beta[0].toFixed(4));
  for (let j = 0; j < mod.varNames.length; j++) {
    console.log('      ' + mod.varNames[j].padEnd(18) + ' = ' + fit.beta[j + 1].toFixed(4) +
      '  [' + boot.ci025[j].toFixed(4) + ', ' + boot.ci975[j].toFixed(4) + ']' +
      '  flip=' + (boot.flipRates[j] * 100).toFixed(1) + '%');
  }
  console.log('    LOO gap%:  ' + gap.toFixed(1) + '%');
  console.log('    LOO RMS:   ' + loo.toFixed(4) + ' dex');
  console.log('    R2adj:     ' + fit.r2adj.toFixed(4));
  console.log('    AIC:       ' + fit.aic.toFixed(2));
  console.log('    BIC:       ' + fit.bic.toFixed(2));
  console.log('    Resid SD:  ' + sd(fit.resid).toFixed(4) + ' dex');
  console.log('');
}

console.log('  ┌───────┬────────────────────────────────────────────┬────────┬──────────┬──────────┬──────────┐');
console.log('  │ Model │ Description                                │ LOO gap│ AIC      │ BIC      │ Resid SD │');
console.log('  ├───────┼────────────────────────────────────────────┼────────┼──────────┼──────────┼──────────┤');
for (const [key, mr] of Object.entries(modelResults)) {
  console.log('  │   ' + key + '   │ ' + mr.name.padEnd(42) + ' │ ' +
    mr.gap.toFixed(1).padStart(5) + '% │ ' +
    mr.fit.aic.toFixed(2).padStart(8) + ' │ ' +
    mr.fit.bic.toFixed(2).padStart(8) + ' │ ' +
    sd(mr.fit.resid).toFixed(4).padStart(8) + ' │');
}
console.log('  └───────┴────────────────────────────────────────────┴────────┴──────────┴──────────┴──────────┘');
console.log('');

const bestKey = Object.entries(modelResults).sort((a, b) => b[1].gap - a[1].gap)[0][0];
const bestModel = modelResults[bestKey];
const baseModel = modelResults['A'];
console.log('  WINNER: Model ' + bestKey + ' (' + bestModel.name + ')');
console.log('  vs baseline (A): delta gap = ' + (bestModel.gap - baseModel.gap).toFixed(1) + 'pp');
console.log('  vs baseline (A): delta AIC = ' + (baseModel.fit.aic - bestModel.fit.aic).toFixed(2));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 2: THRESHOLD SCAN — Is 11.5 optimal?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const scanResults = [];
for (let t = 10.5; t <= 13.0; t += 0.1) {
  const tRound = +t.toFixed(1);
  const nBelow = logMhost.filter(v => v < tRound).length;
  const nAbove = logMhost.filter(v => v >= tRound).length;
  if (nBelow < 5 || nAbove < 5) continue;

  const low = logMhost.map(v => Math.min(v, tRound));
  const high = logMhost.map(v => Math.max(v - tRound, 0));
  const Xb = gals45.map((_, i) => [logMHI[i], low[i], high[i], logMR[i]]);
  const fit = ols(Y, Xb);
  const loo = looCV(Y, Xb);
  const gap = gapPct(loo, sdY);
  scanResults.push({ threshold: tRound, gap, aic: fit.aic, bic: fit.bic, nBelow, nAbove });
}

scanResults.sort((a, b) => b.gap - a.gap);
console.log('  ┌───────────┬─────────┬─────────┬──────────┬──────────┐');
console.log('  │ Threshold │ N_below │ N_above │ LOO gap% │ AIC      │');
console.log('  ├───────────┼─────────┼─────────┼──────────┼──────────┤');
for (const sr of scanResults.slice(0, 15)) {
  const marker = sr.threshold === THRESHOLD ? ' ◄' : '';
  console.log('  │ ' + sr.threshold.toFixed(1).padStart(7) + '   │ ' +
    String(sr.nBelow).padStart(7) + ' │ ' +
    String(sr.nAbove).padStart(7) + ' │ ' +
    sr.gap.toFixed(1).padStart(6) + '% │ ' +
    sr.aic.toFixed(2).padStart(8) + ' │' + marker);
}
console.log('  └───────────┴─────────┴─────────┴──────────┴──────────┘');
console.log('');

const optimalThresh = scanResults[0].threshold;
console.log('  Optimal threshold: ' + optimalThresh.toFixed(1));
console.log('  11.5 rank: #' + (scanResults.findIndex(s => s.threshold === THRESHOLD) + 1) + '/' + scanResults.length);
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 3: logSigma0 RESIDUAL TEST — Does threshold fix absorb it?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

for (const [key, mr] of Object.entries(modelResults)) {
  const rSig0 = pearsonR(mr.fit.resid, logSig0);
  console.log('  Model ' + key + ' residual vs logSigma0: r = ' + rSig0.toFixed(3) +
    (Math.abs(rSig0) > 0.25 ? ' ★★' : Math.abs(rSig0) > 0.15 ? ' ★' : ''));
}
console.log('');

const rSig0_A = pearsonR(modelResults.A.fit.resid, logSig0);
const rSig0_B = pearsonR(modelResults.B.fit.resid, logSig0);
const rSig0_D = pearsonR(modelResults.D.fit.resid, logSig0);

const sig0Absorbed = Math.abs(rSig0_B) < Math.abs(rSig0_A) * 0.6;
console.log('  logSigma0 absorbed by threshold fix?');
console.log('    r(A resid, Sig0) = ' + rSig0_A.toFixed(3));
console.log('    r(B resid, Sig0) = ' + rSig0_B.toFixed(3));
console.log('    Reduction: ' + ((1 - Math.abs(rSig0_B) / Math.abs(rSig0_A)) * 100).toFixed(0) + '%');
console.log('    Verdict: ' + (sig0Absorbed ? 'MOSTLY ABSORBED — Sigma0 was compensating for Mhost misfit' :
  'STILL PRESENT — Sigma0 carries independent information'));
console.log('');

const otherResidCorrs = [
  { name: 'logVflat', arr: gals45.map(g => Math.log10(sparcMap[g.name]?.Vflat || 150)) },
  { name: 'logL36', arr: gals45.map(g => Math.log10(Math.max(sparcMap[g.name]?.L36 || 1, 0.01))) },
  { name: 'morphT', arr: gals45.map(g => sparcMap[g.name]?.T ?? 5) },
  { name: 'logRdisk', arr: gals45.map(g => Math.log10(sparcMap[g.name]?.Rdisk || 3)) },
  { name: 'rcWig', arr: gals45.map(g => g.rcWiggliness) },
  { name: 'envCode', arr: gals45.map(g => g.envCode) }
];

console.log('  All residual correlations for best model (' + bestKey + '):\n');
const bestResid = bestModel.fit.resid;
const allResidCorrs = [{ name: 'logSigma0', r: pearsonR(bestResid, logSig0) }];
for (const v of otherResidCorrs) {
  allResidCorrs.push({ name: v.name, r: pearsonR(bestResid, v.arr) });
}
allResidCorrs.sort((a, b) => Math.abs(b.r) - Math.abs(a.r));
for (const rc of allResidCorrs) {
  const flag = Math.abs(rc.r) > 0.25 ? ' ★★' : Math.abs(rc.r) > 0.15 ? ' ★' : '';
  console.log('    r(resid, ' + rc.name.padEnd(12) + ') = ' + rc.r.toFixed(3) + flag);
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 4: SOURCE SENSITIVITY — Is the threshold robust to KT2017?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const isKT = sources.map(s => s === 'KT2017');
const nKT = isKT.filter(Boolean).length;
const nNonKT = N - nKT;
console.log('  KT2017: N=' + nKT + ', non-KT2017: N=' + nNonKT);
console.log('');

console.log('  4a: Full sample (N=' + N + ')\n');
for (const [key, mod] of Object.entries(models)) {
  console.log('    ' + key + ': gap=' + modelResults[key].gap.toFixed(1) + '%, AIC=' + modelResults[key].fit.aic.toFixed(2));
}
console.log('');

if (nKT >= 10) {
  console.log('  4b: KT2017-only (N=' + nKT + ')\n');
  const ktIdx = gals45.map((_, i) => i).filter(i => isKT[i]);
  const YKT = ktIdx.map(i => Y[i]);
  const sdYKT = sd(YKT);
  for (const [key, mod] of Object.entries(models)) {
    const XKT = ktIdx.map(i => mod.X[i]);
    if (XKT.length <= XKT[0].length + 2) {
      console.log('    ' + key + ': insufficient data');
      continue;
    }
    const fit = ols(YKT, XKT);
    const loo = looCV(YKT, XKT);
    const gap = gapPct(loo, sdYKT);
    console.log('    ' + key + ': gap=' + gap.toFixed(1) + '%, AIC=' + fit.aic.toFixed(2) + ', Resid SD=' + sd(fit.resid).toFixed(4));
  }
  console.log('');
}

if (nNonKT >= 8) {
  console.log('  4c: non-KT2017 only (N=' + nNonKT + ')\n');
  const nktIdx = gals45.map((_, i) => i).filter(i => !isKT[i]);
  const YNKT = nktIdx.map(i => Y[i]);
  const sdYNKT = sd(YNKT);
  for (const [key, mod] of Object.entries(models)) {
    const XNKT = nktIdx.map(i => mod.X[i]);
    if (XNKT.length <= XNKT[0].length + 2) {
      console.log('    ' + key + ': insufficient data');
      continue;
    }
    const fit = ols(YNKT, XNKT);
    const loo = looCV(YNKT, XNKT);
    const gap = gapPct(loo, sdYNKT);
    console.log('    ' + key + ': gap=' + gap.toFixed(1) + '%, AIC=' + fit.aic.toFixed(2) + ', Resid SD=' + sd(fit.resid).toFixed(4));
  }
  console.log('');
} else {
  console.log('  4c: non-KT2017 too small (N=' + nNonKT + ') for reliable model comparison\n');
}

console.log('  4d: Bootstrap source-weighted sensitivity (downsample KT2017)\n');

const rng = mulberry32(777);
const NBOOT_SRC = 1000;
const gapSamplesA = [], gapSamplesB = [];

for (let b = 0; b < NBOOT_SRC; b++) {
  const ktIdx = gals45.map((_, i) => i).filter(i => isKT[i]);
  const nktIdx = gals45.map((_, i) => i).filter(i => !isKT[i]);
  const ktSample = Array.from({ length: Math.min(nNonKT, nKT) }, () => ktIdx[Math.floor(rng() * ktIdx.length)]);
  const combined = [...nktIdx, ...ktSample];
  const Yb = combined.map(i => Y[i]);
  const sdYb = sd(Yb);
  if (sdYb < 0.01) continue;

  const XAb = combined.map(i => models.A.X[i]);
  const XBb = combined.map(i => models.B.X[i]);
  const looA = looCV(Yb, XAb);
  const looB = looCV(Yb, XBb);
  gapSamplesA.push(gapPct(looA, sdYb));
  gapSamplesB.push(gapPct(looB, sdYb));
}

gapSamplesA.sort((a, b) => a - b);
gapSamplesB.sort((a, b) => a - b);
const medA = gapSamplesA[Math.floor(gapSamplesA.length / 2)];
const medB = gapSamplesB[Math.floor(gapSamplesB.length / 2)];
const bWins = gapSamplesB.filter((v, i) => v > gapSamplesA[i]).length;

console.log('  Source-balanced bootstrap (downsample KT2017 to N=' + Math.min(nNonKT, nKT) + '):');
console.log('    Model A median gap: ' + medA.toFixed(1) + '% [' +
  gapSamplesA[Math.floor(gapSamplesA.length * 0.025)].toFixed(1) + ', ' +
  gapSamplesA[Math.floor(gapSamplesA.length * 0.975)].toFixed(1) + ']');
console.log('    Model B median gap: ' + medB.toFixed(1) + '% [' +
  gapSamplesB[Math.floor(gapSamplesB.length * 0.025)].toFixed(1) + ', ' +
  gapSamplesB[Math.floor(gapSamplesB.length * 0.975)].toFixed(1) + ']');
console.log('    B beats A in ' + (bWins / gapSamplesA.length * 100).toFixed(0) + '% of bootstrap samples');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PART 5: NESTED CV — Does the threshold model survive model selection?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const modKeys = ['A', 'B', 'C', 'D'];
let nestedSS = 0;
const selectionCounts = { A: 0, B: 0, C: 0, D: 0 };

for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];
  let bestGap = -Infinity, bestPred = 0, bestMod = '';

  for (const mk of modKeys) {
    const Xfull = models[mk].X;
    const Xtrain = [...Xfull.slice(0, i), ...Xfull.slice(i + 1)];

    let innerSS = 0;
    for (let j = 0; j < Ytrain.length; j++) {
      const Yinn = [...Ytrain.slice(0, j), ...Ytrain.slice(j + 1)];
      const Xinn = [...Xtrain.slice(0, j), ...Xtrain.slice(j + 1)];
      const f = ols(Yinn, Xinn);
      const xi = [1, ...Xtrain[j]];
      innerSS += (Ytrain[j] - xi.reduce((s, x, jj) => s + x * f.beta[jj], 0)) ** 2;
    }
    const innerGap = gapPct(Math.sqrt(innerSS / Ytrain.length), sd(Ytrain));

    if (innerGap > bestGap) {
      bestGap = innerGap;
      bestMod = mk;
      const f = ols(Ytrain, Xtrain);
      const xi = [1, ...Xfull[i]];
      bestPred = xi.reduce((s, x, j) => s + x * f.beta[j], 0);
    }
  }

  selectionCounts[bestMod]++;
  nestedSS += (Y[i] - bestPred) ** 2;
}

const nestedRms = Math.sqrt(nestedSS / N);
const nestedGap = gapPct(nestedRms, sdY);

console.log('  Nested CV (outer LOO, inner LOO for A/B/C/D selection):');
console.log('    Nested gap% = ' + nestedGap.toFixed(1) + '%');
console.log('    Model selection frequency:');
for (const mk of modKeys) {
  console.log('      ' + mk + ': ' + selectionCounts[mk] + '/' + N + ' (' + (selectionCounts[mk] / N * 100).toFixed(0) + '%)');
}
console.log('    Best standard model gap: ' + Math.max(...Object.values(modelResults).map(m => m.gap)).toFixed(1) + '%');
console.log('    Optimism: ' + (Math.max(...Object.values(modelResults).map(m => m.gap)) - nestedGap).toFixed(1) + 'pp');
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 125: SYNTHESIS');
console.log('══════════════════════════════════════════════════════════════════════\n');

const thresholdHelps = bestModel.gap > baseModel.gap + 2;
const sig0StillPresent = Math.abs(pearsonR(bestModel.fit.resid, logSig0)) > 0.2;
const sourceRobust = bWins / gapSamplesA.length > 0.55;
const nestedSurvives = nestedGap > 20;

console.log('  KEY FINDINGS:\n');
console.log('  1. Does threshold representation improve M3?');
console.log('     Winner: Model ' + bestKey + ' (gap ' + bestModel.gap.toFixed(1) + '% vs A=' + baseModel.gap.toFixed(1) + '%)');
console.log('     Delta: ' + (bestModel.gap - baseModel.gap).toFixed(1) + 'pp');
console.log('     Verdict: ' + (thresholdHelps ? 'YES — threshold form significantly improves M3' : 'MARGINAL or NO — threshold does not clearly help'));
console.log('');
console.log('  2. Does logSigma0 survive after threshold correction?');
console.log('     r(best resid, Sig0) = ' + pearsonR(bestModel.fit.resid, logSig0).toFixed(3));
console.log('     Verdict: ' + (sig0StillPresent ? 'YES — Sigma0 is independently real → Phase 126 warranted' : 'NO — Sigma0 was compensating for Mhost misfit'));
console.log('');
console.log('  3. Is the threshold robust to source composition?');
console.log('     B beats A in ' + (bWins / gapSamplesA.length * 100).toFixed(0) + '% of source-balanced bootstrap');
console.log('     Verdict: ' + (sourceRobust ? 'YES — threshold holds under source rebalancing' : 'UNCERTAIN — may depend on KT2017 dominance'));
console.log('');
console.log('  4. Nested CV survival:');
console.log('     gap = ' + nestedGap.toFixed(1) + '%');
console.log('     Verdict: ' + (nestedSurvives ? 'PASS' : 'MARGINAL'));
console.log('');
console.log('  OPTIMAL THRESHOLD: ' + optimalThresh.toFixed(1) + ' (11.5 is #' + (scanResults.findIndex(s => s.threshold === THRESHOLD) + 1) + ')');
console.log('');

let nextStep;
if (thresholdHelps && sig0StillPresent) {
  nextStep = 'Phase 126: M4 candidate with threshold-corrected Mhost + logSigma0';
} else if (thresholdHelps && !sig0StillPresent) {
  nextStep = 'Sigma0 absorbed by threshold fix. Best law = threshold-env + MHI + MR.';
} else if (!thresholdHelps && sig0StillPresent) {
  nextStep = 'Threshold does not help. Phase 126: M4 candidate with linear Mhost + Sigma0.';
} else {
  nextStep = 'Neither threshold nor Sigma0 adds value. M3-linear is the final form.';
}
console.log('  RECOMMENDED NEXT STEP: ' + nextStep);
console.log('');

const output = {
  phase: '125',
  title: 'Environmental Threshold Law',
  threshold: THRESHOLD,
  optimal_threshold: optimalThresh,
  models: {},
  threshold_scan: scanResults.slice(0, 15).map(s => ({
    threshold: s.threshold, gap: +s.gap.toFixed(1), aic: +s.aic.toFixed(2),
    n_below: s.nBelow, n_above: s.nAbove
  })),
  sigma0_test: {
    r_A: +rSig0_A.toFixed(3),
    r_B: +rSig0_B.toFixed(3),
    r_D: +rSig0_D.toFixed(3),
    r_best: +pearsonR(bestModel.fit.resid, logSig0).toFixed(3),
    absorbed: sig0Absorbed,
    still_present: sig0StillPresent
  },
  source_sensitivity: {
    kt2017_n: nKT,
    non_kt_n: nNonKT,
    b_beats_a_pct: +(bWins / gapSamplesA.length * 100).toFixed(0),
    robust: sourceRobust
  },
  nested_cv: {
    gap: +nestedGap.toFixed(1),
    selection: selectionCounts,
    survives: nestedSurvives
  },
  synthesis: {
    threshold_helps: thresholdHelps,
    sigma0_independent: sig0StillPresent,
    source_robust: sourceRobust,
    recommended_next: nextStep
  }
};

for (const [key, mr] of Object.entries(modelResults)) {
  output.models[key] = {
    name: mr.name,
    gap: +mr.gap.toFixed(1),
    loo_rms: +mr.loo.toFixed(4),
    r2adj: +mr.fit.r2adj.toFixed(4),
    aic: +mr.fit.aic.toFixed(2),
    bic: +mr.fit.bic.toFixed(2),
    resid_sd: +sd(mr.fit.resid).toFixed(4),
    coefficients: Object.fromEntries(
      [['intercept', +mr.fit.beta[0].toFixed(4)],
      ...mr.varNames.map((v, j) => [v, +mr.fit.beta[j + 1].toFixed(4)])]
    ),
    bootstrap_flip_rates: Object.fromEntries(mr.varNames.map((v, j) => [v, +(mr.boot.flipRates[j] * 100).toFixed(1)]))
  };
}

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase125-environmental-threshold.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase125-environmental-threshold.json');

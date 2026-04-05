#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 123: M3 STRESS TEST + REGIME TEST');
console.log('');
console.log('  Question: Is M3 a genuine state law, or only a compressed');
console.log('  statistical bundle?');
console.log('');
console.log('  M3: log(a0) = c0 + c1*logMHI + c2*logMhost + c3*logMeanRun');
console.log('');
console.log('  Three tests:');
console.log('    A) Death match vs all competing models (1/2/3-var)');
console.log('    B) Vflat regime threshold scan');
console.log('    C) Physical axis decomposition');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'stage-A-master-table.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparc.forEach(s => { sparcMap[s.name] = s; });
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));

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
  const r2adj = 1 - (rss / (n - p)) / (tss / (n - 1));
  const aic = n * Math.log(rss / n) + 2 * p;
  const bic = n * Math.log(rss / n) + p * Math.log(n);
  return { beta, resid, rss, tss, r2adj, aic, bic, n, k: p };
}

function looCV(Y, X) {
  const n = Y.length;
  let ss = 0;
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

function shuffle(arr, rng) {
  const a = [...arr];
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
}

function bootstrapSignStability(Y, X, nBoot, seed) {
  const rng = mulberry32(seed);
  const n = Y.length, p = X[0].length;
  const baseFit = ols(Y, X);
  const baseSigns = baseFit.beta.slice(1).map(b => Math.sign(b));
  const flips = new Array(p).fill(0);
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(rng() * n));
    const Yb = idx.map(i => Y[i]);
    const Xb = idx.map(i => X[i]);
    const fit = ols(Yb, Xb);
    for (let j = 0; j < p; j++) {
      if (Math.sign(fit.beta[j + 1]) !== baseSigns[j]) flips[j]++;
    }
  }
  return flips.map(f => f / nBoot);
}

const upsilonMap = {
  'NGC0024':0.50,'NGC0289':0.47,'NGC0891':0.61,'NGC1003':0.40,
  'NGC1090':0.45,'NGC1705':0.26,'NGC2403':0.45,'NGC2683':0.52,
  'NGC2841':0.74,'NGC2903':0.57,'NGC2915':0.22,'NGC3198':0.47,
  'NGC3521':0.60,'NGC3726':0.33,'NGC3741':0.18,'NGC3769':0.37,
  'NGC3893':0.44,'NGC4013':0.50,'NGC4100':0.49,'NGC4138':0.79,
  'NGC4157':0.47,'NGC4217':0.55,'NGC4559':0.22,'NGC5005':0.53,
  'NGC5033':0.53,'NGC5055':0.56,'NGC5371':0.50,'NGC5907':0.48,
  'NGC6015':0.47,'NGC6503':0.52,'NGC6674':0.55,'NGC7331':0.58,
  'NGC7814':0.71,
  'UGC01281':0.25,'UGC02953':0.55,'UGC03205':0.55,'UGC03546':0.60,
  'UGC03580':0.55,'UGC05721':0.30,'UGC06786':0.55,'UGC06787':0.55,
  'UGC06973':0.50,'UGC08490':0.25,'UGC08699':0.55,'UGC09133':0.50,
  'F571-8':0.30
};

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const morphT = gals45.map(g => sparcMap[g.name]?.T ?? 5);
const logUps = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logMHI = gals45.map(g => g.logMHI);
const logSig0 = gals45.map(g => g.logSigma0);
const logMR = gals45.map(g => g.logMeanRun);
const rcWig = gals45.map(g => g.rcWiggliness);
const Vflat = gals45.map(g => sparcMap[g.name]?.Vflat || g.Vflat_km_s || 150);
const logRdisk = gals45.map(g => Math.log10(sparcMap[g.name]?.Rdisk || 3));
const L36 = gals45.map(g => sparcMap[g.name]?.L36 || 1);
const logL36 = L36.map(v => Math.log10(v > 0 ? v : 1));
const logVflat = Vflat.map(v => Math.log10(v));
const logConc = gals45.map((g, i) => logMHI[i] - logL36[i]);

const Xconf = gals45.map((g, i) => [logMHI[i], logSig0[i], morphT[i]]);
const fConf = ols(logUps, Xconf);
const upsPerp = fConf.resid;

const allVarArrays = {
  logMHI, logMhost, logMR, logSig0, logVflat, logL36, logConc, logRdisk
};
const varNames = Object.keys(allVarArrays);

console.log('  N = ' + gals45.length + ', SD(logA0) = ' + sdY.toFixed(4) + ' dex\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST A: DEATH MATCH — M3 vs ALL competing models');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const results = [];

for (let i = 0; i < varNames.length; i++) {
  const X = gals45.map((_, k) => [allVarArrays[varNames[i]][k]]);
  const loo = looCV(Y, X);
  const fit = ols(Y, X);
  results.push({
    name: varNames[i], nvars: 1,
    gap: gapPct(loo, sdY), r2adj: fit.r2adj,
    aic: fit.aic, bic: fit.bic, loo_rms: loo
  });
}

for (let i = 0; i < varNames.length; i++) {
  for (let j = i + 1; j < varNames.length; j++) {
    const X = gals45.map((_, k) => [allVarArrays[varNames[i]][k], allVarArrays[varNames[j]][k]]);
    const loo = looCV(Y, X);
    const fit = ols(Y, X);
    results.push({
      name: varNames[i] + '+' + varNames[j], nvars: 2,
      gap: gapPct(loo, sdY), r2adj: fit.r2adj,
      aic: fit.aic, bic: fit.bic, loo_rms: loo
    });
  }
}

for (let i = 0; i < varNames.length; i++) {
  for (let j = i + 1; j < varNames.length; j++) {
    for (let l = j + 1; l < varNames.length; l++) {
      const X = gals45.map((_, k) => [allVarArrays[varNames[i]][k], allVarArrays[varNames[j]][k], allVarArrays[varNames[l]][k]]);
      const loo = looCV(Y, X);
      const fit = ols(Y, X);
      results.push({
        name: varNames[i] + '+' + varNames[j] + '+' + varNames[l], nvars: 3,
        gap: gapPct(loo, sdY), r2adj: fit.r2adj,
        aic: fit.aic, bic: fit.bic, loo_rms: loo
      });
    }
  }
}

results.sort((a, b) => b.gap - a.gap);

console.log('  ┌─────┬──────────────────────────────────────────────┬────────┬────────┬──────────┐');
console.log('  │ Rank│ Model                                        │ k      │ LOO gap│ R²adj    │');
console.log('  ├─────┼──────────────────────────────────────────────┼────────┼────────┼──────────┤');
const m3Rank = results.findIndex(r => r.name === 'logMHI+logMhost+logMR') + 1;
for (let i = 0; i < Math.min(25, results.length); i++) {
  const r = results[i];
  const isM3 = r.name === 'logMHI+logMhost+logMR';
  const marker = isM3 ? ' ◄◄◄' : '';
  console.log('  │ ' + String(i + 1).padStart(3) + ' │ ' +
    (r.name + marker).padEnd(44) + ' │ ' +
    String(r.nvars).padStart(6) + ' │ ' +
    r.gap.toFixed(1).padStart(5) + '% │ ' +
    r.r2adj.toFixed(4).padStart(8) + ' │');
}
console.log('  └─────┴──────────────────────────────────────────────┴────────┴────────┴──────────┘');
console.log('');
console.log('  M3 rank: #' + m3Rank + ' out of ' + results.length + ' models tested');

const m3Result = results.find(r => r.name === 'logMHI+logMhost+logMR');
const sameOrBetter3var = results.filter(r => r.nvars === 3 && r.gap >= m3Result.gap).length;
const total3var = results.filter(r => r.nvars === 3).length;
console.log('  Among 3-var models: M3 is #' + sameOrBetter3var + ' of ' + total3var);
console.log('');

const bestByK = {};
for (const r of results) {
  if (!bestByK[r.nvars] || r.gap > bestByK[r.nvars].gap) bestByK[r.nvars] = r;
}
console.log('  Best by complexity:');
for (const k of [1, 2, 3]) {
  const b = bestByK[k];
  if (b) console.log('    k=' + k + ': ' + b.name + ' (gap=' + b.gap.toFixed(1) + '%)');
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST A2: BOOTSTRAP SIGN STABILITY — M3 vs top competitors');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const NBOOT = 2000;
const m3X = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]);
const m3FlipRates = bootstrapSignStability(Y, m3X, NBOOT, 42);
console.log('  M3 bootstrap sign-flip rates (N=' + NBOOT + '):');
console.log('    logMHI:    ' + (m3FlipRates[0] * 100).toFixed(1) + '%');
console.log('    logMhost:  ' + (m3FlipRates[1] * 100).toFixed(1) + '%');
console.log('    logMeanRun: ' + (m3FlipRates[2] * 100).toFixed(1) + '%');
const m3AllStable = m3FlipRates.every(f => f < 0.10);
console.log('    All < 10%? ' + (m3AllStable ? 'YES ✓' : 'NO ✗'));
console.log('');

const top3 = results.filter(r => r.nvars === 3 && r.name !== 'logMHI+logMhost+logMR').slice(0, 3);
for (const comp of top3) {
  const parts = comp.name.split('+');
  const compX = gals45.map((_, i) => parts.map(p => allVarArrays[p][i]));
  const compFlips = bootstrapSignStability(Y, compX, NBOOT, 42);
  console.log('  ' + comp.name + ' sign-flip rates:');
  parts.forEach((p, j) => console.log('    ' + p + ': ' + (compFlips[j] * 100).toFixed(1) + '%'));
  const allStable = compFlips.every(f => f < 0.10);
  console.log('    All < 10%? ' + (allStable ? 'YES ✓' : 'NO ✗'));
  console.log('');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST A3: PERMUTATION TEST — Is M3 significantly better than chance?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const NPERM = 5000;
const m3RealGap = m3Result.gap;
let nBeat = 0;
const rng = mulberry32(123);
for (let p = 0; p < NPERM; p++) {
  const Yp = shuffle(Y, rng);
  const loo = looCV(Yp, m3X);
  const gap = gapPct(loo, sd(Yp));
  if (gap >= m3RealGap) nBeat++;
}
const permP = (nBeat + 1) / (NPERM + 1);
console.log('  M3 real LOO gap% = ' + m3RealGap.toFixed(1) + '%');
console.log('  Permutations beating real: ' + nBeat + '/' + NPERM);
console.log('  Permutation p = ' + permP.toFixed(4));
console.log('  Verdict: ' + (permP < 0.01 ? 'SIGNIFICANT (p < 0.01) ✓' : permP < 0.05 ? 'SIGNIFICANT (p < 0.05) ✓' : 'NOT SIGNIFICANT ✗'));
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST B: VFLAT REGIME THRESHOLD SCAN');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const thresholds = [0, 60, 70, 80, 90, 100, 120, 150];
const threshResults = [];

console.log('  ┌──────────┬─────┬──────────┬──────────┬────────────────────────────────────┐');
console.log('  │ Vflat ≥  │ N   │ LOO gap% │ LOO RMS  │ M3 coefficients [MHI, Mhost, MR]   │');
console.log('  ├──────────┼─────┼──────────┼──────────┼────────────────────────────────────┤');

for (const thresh of thresholds) {
  const mask = gals45.map((g, i) => Vflat[i] >= thresh);
  const Ys = Y.filter((_, i) => mask[i]);
  const Xs = m3X.filter((_, i) => mask[i]);
  const n = Ys.length;

  if (n < 8) {
    console.log('  │ ' + String(thresh).padStart(6) + '   │ ' + String(n).padStart(3) + ' │    —     │    —     │    (insufficient data)              │');
    threshResults.push({ thresh, n, gap: null, rms: null, coeffs: null });
    continue;
  }

  const loo = looCV(Ys, Xs);
  const gap = gapPct(loo, sd(Ys));
  const fit = ols(Ys, Xs);
  const coeffs = fit.beta.slice(1);
  const coeffStr = coeffs.map(c => (c >= 0 ? '+' : '') + c.toFixed(3)).join(', ');

  console.log('  │ ' + String(thresh).padStart(6) + '   │ ' + String(n).padStart(3) + ' │ ' +
    gap.toFixed(1).padStart(6) + '%  │ ' + loo.toFixed(4).padStart(8) + ' │ [' + coeffStr + ']' +
    ' '.repeat(Math.max(0, 30 - coeffStr.length)) + '│');

  threshResults.push({ thresh, n, gap, rms: loo, coeffs, signs: coeffs.map(c => Math.sign(c)) });
}
console.log('  └──────────┴─────┴──────────┴──────────┴────────────────────────────────────┘');
console.log('');

const fullGap = threshResults.find(t => t.thresh === 0)?.gap || 0;
const signStable = threshResults.filter(t => t.coeffs && t.signs.join(',') === '-1,-1,1');
console.log('  Sign stability: ' + signStable.length + '/' + threshResults.filter(t => t.coeffs).length +
  ' thresholds maintain correct signs [-,-,+]');

const peakThresh = threshResults.filter(t => t.gap !== null).sort((a, b) => b.gap - a.gap)[0];
console.log('  Peak performance: Vflat >= ' + peakThresh.thresh + ' (N=' + peakThresh.n + ', gap=' + peakThresh.gap.toFixed(1) + '%)');
console.log('');

const strengthening = threshResults.filter(t => t.gap !== null && t.thresh >= 70 && t.gap > fullGap);
if (strengthening.length > 0) {
  console.log('  REGIME SIGNAL: M3 strengthens above Vflat thresholds:');
  strengthening.forEach(t => console.log('    Vflat >= ' + t.thresh + ': gap=' + t.gap.toFixed(1) + '% (vs full=' + fullGap.toFixed(1) + '%)'));
} else {
  console.log('  No clear regime strengthening detected.');
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST C: PHYSICAL AXIS DECOMPOSITION');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  C1: Individual axis contributions (1-var LOO)\n');

const axes = [
  { name: 'logMHI', arr: logMHI, phys: 'Gas reservoir / evolutionary state' },
  { name: 'logMhost', arr: logMhost, phys: 'Environmental depth / tidal processing' },
  { name: 'logMeanRun', arr: logMR, phys: 'Dynamical coherence / inner structure' }
];

for (const ax of axes) {
  const r = pearsonR(ax.arr, Y);
  const X1 = gals45.map((_, i) => [ax.arr[i]]);
  const loo = looCV(Y, X1);
  const gap = gapPct(loo, sdY);
  console.log('  ' + ax.name);
  console.log('    r(logA0) = ' + r.toFixed(3));
  console.log('    1-var LOO gap% = ' + gap.toFixed(1) + '%');
  console.log('    Physics: ' + ax.phys);
  console.log('');
}

console.log('  C2: Axis independence (pairwise correlations among M3 predictors)\n');

for (let i = 0; i < axes.length; i++) {
  for (let j = i + 1; j < axes.length; j++) {
    const r = pearsonR(axes[i].arr, axes[j].arr);
    console.log('    r(' + axes[i].name + ', ' + axes[j].name + ') = ' + r.toFixed(3));
  }
}
console.log('');

console.log('  C3: Incremental contribution (drop-one analysis)\n');

const m3Fit = ols(Y, m3X);
const m3LooRms = looCV(Y, m3X);
const m3Gap = gapPct(m3LooRms, sdY);
console.log('  Full M3: gap=' + m3Gap.toFixed(1) + '%, R2adj=' + m3Fit.r2adj.toFixed(4));
console.log('');

for (let drop = 0; drop < 3; drop++) {
  const kept = [0, 1, 2].filter(j => j !== drop);
  const Xdrop = gals45.map((_, i) => kept.map(j => m3X[i][j]));
  const loo = looCV(Y, Xdrop);
  const gap = gapPct(loo, sdY);
  const fit = ols(Y, Xdrop);
  const delta = m3Gap - gap;
  console.log('  Drop ' + axes[drop].name + ':');
  console.log('    gap=' + gap.toFixed(1) + '% (delta=' + (delta > 0 ? '-' : '+') + Math.abs(delta).toFixed(1) + 'pp)');
  console.log('    R2adj=' + fit.r2adj.toFixed(4));
  console.log('');
}

console.log('  C4: logMeanRun deep dive — what does it encode?\n');

const meanRunCorrs = [
  { name: 'logMHI', arr: logMHI },
  { name: 'logMhost', arr: logMhost },
  { name: 'logVflat', arr: logVflat },
  { name: 'logSigma0', arr: logSig0 },
  { name: 'logL36', arr: logL36 },
  { name: 'logRdisk', arr: logRdisk },
  { name: 'logConc(MHI/L36)', arr: logConc },
  { name: 'morphT', arr: morphT.map(t => t) }
];

console.log('  Correlations with logMeanRun:');
for (const v of meanRunCorrs) {
  const r = pearsonR(logMR, v.arr);
  console.log('    r(' + v.name + ') = ' + r.toFixed(3) + (Math.abs(r) > 0.4 ? ' ***' : Math.abs(r) > 0.3 ? ' **' : ''));
}
console.log('');

const logMR_residM2 = (() => {
  const X2 = gals45.map((_, i) => [logMHI[i], logMhost[i]]);
  const fitMR = ols(logMR, X2);
  return fitMR.resid;
})();

const r_MRresid_a0 = pearsonR(logMR_residM2, Y);
console.log('  logMeanRun residual (after removing MHI+Mhost dependence):');
console.log('    r(residMR, logA0) = ' + r_MRresid_a0.toFixed(3));
console.log('    This is the INDEPENDENT information MeanRun adds beyond mass axes.');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST D: NESTED CV — Does M3 survive model selection inside CV?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

let m3wins = 0, bestWins = 0;
let nestedSS = 0;
const N = gals45.length;
const competitors3var = results.filter(r => r.nvars === 3).slice(0, 5).map(r => r.name);

for (let i = 0; i < N; i++) {
  const Ytrain = [...Y.slice(0, i), ...Y.slice(i + 1)];

  let bestGap = -Infinity, bestPred = 0;

  for (const modName of competitors3var) {
    const parts = modName.split('+');
    const Xfull = gals45.map((_, k) => parts.map(p => allVarArrays[p][k]));
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
      const f = ols(Ytrain, Xtrain);
      const xi = [1, ...Xfull[i]];
      bestPred = xi.reduce((s, x, j) => s + x * f.beta[j], 0);
      if (modName === 'logMHI+logMhost+logMR') m3wins++;
    }
  }

  nestedSS += (Y[i] - bestPred) ** 2;
}

const nestedRms = Math.sqrt(nestedSS / N);
const nestedGap = gapPct(nestedRms, sdY);
console.log('  Nested CV (outer LOO, inner LOO for model selection):');
console.log('    Nested gap% = ' + nestedGap.toFixed(1) + '%');
console.log('    M3 selected in ' + m3wins + '/' + N + ' outer folds (' + (m3wins / N * 100).toFixed(0) + '%)');
console.log('    Standard M3 gap% = ' + m3Gap.toFixed(1) + '%');
console.log('    Optimism = ' + (m3Gap - nestedGap).toFixed(1) + ' pp');
console.log('    Verdict: ' + (nestedGap > 25 ? 'PASS' : 'MARGINAL') + ' (nested gap > 25% threshold)');
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 123: FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const verdictA = m3Rank <= 5 && m3AllStable && permP < 0.01;
const verdictB = signStable.length >= 4;
const verdictC = Math.abs(r_MRresid_a0) > 0.2;
const verdictD = nestedGap > 25;

console.log('  A) Death match rank: #' + m3Rank + '/' + results.length +
  (m3Rank <= 5 ? ' ✓' : ' ✗') +
  ' | Signs stable: ' + (m3AllStable ? '✓' : '✗') +
  ' | Permutation: p=' + permP.toFixed(4) + (permP < 0.01 ? ' ✓' : ' ✗'));
console.log('  B) Regime stability: ' + signStable.length + '/' + threshResults.filter(t => t.coeffs).length +
  ' thresholds sign-stable' + (verdictB ? ' ✓' : ' ✗') +
  ' | Peak at Vflat>=' + peakThresh.thresh);
console.log('  C) MeanRun independence: r(resid,a0)=' + r_MRresid_a0.toFixed(3) +
  (verdictC ? ' ✓' : ' ✗') +
  ' — carries independent dynamical info');
console.log('  D) Nested CV: gap=' + nestedGap.toFixed(1) + '%' + (verdictD ? ' ✓' : ' ✗'));
console.log('');

const overall = [verdictA, verdictB, verdictC, verdictD].filter(Boolean).length;
let verdict;
if (overall === 4) verdict = 'STRONG PASS — M3 behaves like a genuine state law';
else if (overall >= 3) verdict = 'PASS — M3 is more than a statistical bundle';
else if (overall >= 2) verdict = 'MARGINAL — M3 shows promise but needs more evidence';
else verdict = 'FAIL — M3 may be a compressed statistical artifact';

console.log('  OVERALL: ' + overall + '/4 tests passed');
console.log('  VERDICT: ' + verdict);
console.log('');

const output = {
  phase: '123',
  title: 'M3 Stress Test + Regime Test',
  question: 'Is M3 a genuine state law, or only a compressed statistical bundle?',
  m3: {
    formula: 'logA0 = ' + m3Fit.beta[0].toFixed(4) + ' + ' +
      m3Fit.beta[1].toFixed(4) + '*logMHI + ' +
      m3Fit.beta[2].toFixed(4) + '*logMhost + ' +
      m3Fit.beta[3].toFixed(4) + '*logMeanRun',
    gap_pct: +m3Gap.toFixed(1),
    loo_rms: +m3LooRms.toFixed(4),
    r2adj: +m3Fit.r2adj.toFixed(4),
    coefficients: {
      intercept: +m3Fit.beta[0].toFixed(4),
      logMHI: +m3Fit.beta[1].toFixed(4),
      logMhost: +m3Fit.beta[2].toFixed(4),
      logMeanRun: +m3Fit.beta[3].toFixed(4)
    }
  },
  testA: {
    title: 'Death Match',
    total_models: results.length,
    m3_rank: m3Rank,
    top10: results.slice(0, 10).map(r => ({ name: r.name, k: r.nvars, gap: +r.gap.toFixed(1) })),
    bootstrap_sign_flips: {
      logMHI: +(m3FlipRates[0] * 100).toFixed(1),
      logMhost: +(m3FlipRates[1] * 100).toFixed(1),
      logMeanRun: +(m3FlipRates[2] * 100).toFixed(1),
      all_stable: m3AllStable
    },
    permutation: { p: +permP.toFixed(4), n_perms: NPERM, significant: permP < 0.01 }
  },
  testB: {
    title: 'Vflat Regime Scan',
    thresholds: threshResults.map(t => ({
      vflat_min: t.thresh, n: t.n,
      gap: t.gap !== null ? +t.gap.toFixed(1) : null,
      rms: t.rms !== null ? +t.rms.toFixed(4) : null,
      signs_correct: t.signs ? t.signs.join(',') === '-1,-1,1' : null
    })),
    sign_stable_count: signStable.length,
    peak: { vflat_min: peakThresh.thresh, gap: +peakThresh.gap.toFixed(1) }
  },
  testC: {
    title: 'Physical Axis Decomposition',
    individual: axes.map(ax => ({
      name: ax.name,
      r_logA0: +pearsonR(ax.arr, Y).toFixed(3),
      loo_gap_1var: +(gapPct(looCV(Y, gals45.map((_, i) => [ax.arr[i]])), sdY)).toFixed(1),
      physics: ax.phys
    })),
    meanRun_independence: +r_MRresid_a0.toFixed(3),
    axis_correlations: {
      MHI_Mhost: +pearsonR(logMHI, logMhost).toFixed(3),
      MHI_MR: +pearsonR(logMHI, logMR).toFixed(3),
      Mhost_MR: +pearsonR(logMhost, logMR).toFixed(3)
    }
  },
  testD: {
    title: 'Nested CV',
    nested_gap: +nestedGap.toFixed(1),
    standard_gap: +m3Gap.toFixed(1),
    optimism: +(m3Gap - nestedGap).toFixed(1),
    m3_selected_pct: +(m3wins / N * 100).toFixed(0)
  },
  verdict: {
    tests_passed: overall,
    tests_total: 4,
    result: verdict,
    details: { A: verdictA, B: verdictB, C: verdictC, D: verdictD }
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase123-m3-stress-test.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase123-m3-stress-test.json');

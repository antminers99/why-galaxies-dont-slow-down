#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 132A: HALO DRIVER DEATH MATCH');
console.log('');
console.log('  Every available halo/DM proxy vs VfResid in a direct comparison.');
console.log('  Which halo observable best explains the coupling residual?');
console.log('  Can any halo variable approach VfResid predictive power?');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'stage-A-master-table.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const sparcTable = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-results.json'), 'utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparcTable.forEach(s => { sparcMap[s.name] = s; });
const resMap = {}; sparcResults.perGalaxy.forEach(g => { if (g.rawName) resMap[g.rawName] = g; });
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
  return { beta, resid, rss, tss, r2: tss > 0 ? 1 - rss / tss : 0 };
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
function nestedCV(Y, Xbase, Xfull) {
  const n = Y.length; let wins = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xbt = [...Xbase.slice(0, i), ...Xbase.slice(i + 1)];
    const Xft = [...Xfull.slice(0, i), ...Xfull.slice(i + 1)];
    const fb = ols(Yt, Xbt); const ff = ols(Yt, Xft);
    const pb = [1, ...Xbase[i]].reduce((s, x, j) => s + x * fb.beta[j], 0);
    const pf = [1, ...Xfull[i]].reduce((s, x, j) => s + x * ff.beta[j], 0);
    if (Math.abs(Y[i] - pf) <= Math.abs(Y[i] - pb)) wins++;
  }
  return wins;
}
function gapPct(rms, sdy) { return 100 * (1 - rms ** 2 / sdy ** 2); }
function partialR(x, y, controls) {
  const rx = ols(x, controls).resid;
  const ry = ols(y, controls).resid;
  return pearsonR(rx, ry);
}
function bootstrap(Y, X, nBoot) {
  const n = Y.length;
  const refBeta = ols(Y, X).beta;
  const signs = refBeta.slice(1).map(b => Math.sign(b));
  let flips = 0;
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(Math.random() * n));
    const Yb = idx.map(i => Y[i]); const Xb = idx.map(i => X[i]);
    try {
      const fb = ols(Yb, Xb);
      const bSigns = fb.beta.slice(1).map(b => Math.sign(b));
      if (bSigns.some((s, i) => s !== signs[i])) flips++;
    } catch (e) { flips++; }
  }
  return flips / nBoot;
}

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const logMHI = gals45.map(g => g.logMHI);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logMR = gals45.map(g => g.logMeanRun);
const logVflat = gals45.map(g => Math.log10(sparcMap[g.name].Vflat));
const logL36 = gals45.map(g => Math.log10(Math.max(sparcMap[g.name].L36, 0.01)));
const logRdisk = gals45.map(g => Math.log10(sparcMap[g.name].Rdisk));
const morphT = gals45.map(g => sparcMap[g.name].T);
const logMbar = gals45.map(g => {
  const s = sparcMap[g.name];
  return Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
});
const envCode = gals45.map(g => g.envCode);
const logSig0 = gals45.map(g => g.logSigma0);

const structX = gals45.map((_, i) => [logMbar[i], logL36[i], logRdisk[i], morphT[i]]);
const VfResid = ols(logVflat, structX).resid;
console.log('  VfResid: SD=' + sd(VfResid).toFixed(4) + ' r(logA0)=' + pearsonR(VfResid, Y).toFixed(3) + '\n');

const core3X = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]);
const coreStructX = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i], logMbar[i], logL36[i], logRdisk[i], morphT[i]]);

function safeLog(v) { return v > 0 ? Math.log10(v) : Math.log10(Math.max(Math.abs(v), 1)); }

const haloProxies = [];

gals45.forEach((g, i) => {
  const r = resMap[g.name];
  const s = sparcMap[g.name];
  const dhl = r.models.dark_halo_linear;
  const dhf = r.models.dark_halo_flat;
  const lh = r.models.log_halo;
  const mo = r.models.mond;
  const mgh = r.models.modified_gravity_halo;
  const tr = r.models.transition;
  const newt = r.models.newtonian;

  const Vflat2 = s.Vflat * s.Vflat;
  const Mstar = s.L36 * 0.5 * 1e9;
  const Mgas = Math.pow(10, g.logMHI) * 1.33 * 1e9;
  const Mbar_val = Mstar + Mgas;
  const Vbar2_approx = 6.674e-11 * Mbar_val * 1.989e30 / (s.Rdisk * 3.086e19) / 1e6;
  const fDM_flat = Math.max(0, 1 - Math.min(1, Vbar2_approx / Vflat2));
  const VDM2 = Vflat2 - Vbar2_approx;

  const row = {
    haloK_linear: safeLog(Math.max(dhl.k, 1)),
    haloK_flat: safeLog(Math.max(dhf.k, 1)),
    haloK_log: safeLog(Math.max(lh.k, 1)),
    dhl_improve: dhl.improvementVsNewton,
    dhf_improve: dhf.improvementVsNewton,
    lh_improve: lh.improvementVsNewton,
    mond_a0: safeLog(Math.max(mo.a, 0.01)),
    mond_improve: mo.improvementVsNewton,
    mond_mse: Math.log10(Math.max(mo.mse, 0.1)),
    mgh_k: safeLog(Math.max(mgh.k, 1)),
    mgh_a: safeLog(Math.max(mgh.a, 0.01)),
    transition_k: safeLog(Math.max(tr.k, 0.1)),
    transition_a: safeLog(Math.max(tr.a, 0.01)),
    newt_mse: Math.log10(Math.max(newt.mse, 0.1)),
    dhl_mse: Math.log10(Math.max(dhl.mse, 0.1)),
    lh_mse: Math.log10(Math.max(lh.mse, 0.1)),
    dhl_innerImprove: dhl.innerImprovement,
    dhl_outerImprove: dhl.outerImprovement,
    lh_innerImprove: lh.innerImprovement,
    lh_outerImprove: lh.outerImprovement,
    fDM_flat: fDM_flat,
    logVDM2: VDM2 > 0 ? Math.log10(VDM2) : 0,
    logFdm: fDM_flat > 0 ? Math.log10(fDM_flat) : -2,
    halo_baryon_ratio: VDM2 > 0 && Vbar2_approx > 0 ? Math.log10(VDM2 / Vbar2_approx) : 0,
    dhl_mse_ratio: Math.log10(Math.max(dhl.mse / Math.max(newt.mse, 0.1), 0.001)),
    lh_mse_ratio: Math.log10(Math.max(lh.mse / Math.max(newt.mse, 0.1), 0.001)),
    mond_mse_ratio: Math.log10(Math.max(mo.mse / Math.max(newt.mse, 0.1), 0.001)),
    dhl_outer_inner: dhl.outerImprovement - dhl.innerImprovement,
    lh_outer_inner: lh.outerImprovement - lh.innerImprovement,
    logMaxV_Vflat: Math.log10(r.maxV / s.Vflat),
    logExtent: Math.log10(r.maxR / s.Rdisk),
    logMhost: tdMap[g.name].logMhost,
    envCode: g.envCode,
    logMeanRun: g.logMeanRun,
    concentration: Math.log10(s.Reff / s.Rdisk),
    logSBeff: Math.log10(Math.max(s.SBeff, 0.01)),
  };
  haloProxies.push(row);
});

const proxyNames = Object.keys(haloProxies[0]);

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 1: RAW CORRELATIONS — r(VfResid, proxy)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const corrs = proxyNames.map(name => {
  const arr = haloProxies.map(h => h[name]);
  const r = pearsonR(VfResid, arr);
  return { name, r, rAbs: Math.abs(r) };
}).sort((a, b) => b.rAbs - a.rAbs);

for (const c of corrs) {
  const flag = c.rAbs > 0.5 ? '★★★' : c.rAbs > 0.3 ? '★★ ' : c.rAbs > 0.15 ? '★  ' : '   ';
  console.log('  ' + flag + ' r=' + c.r.toFixed(3).padStart(7) + '  ' + c.name);
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 2: PARTIAL CORRELATIONS — After removing Core + Structure');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const top15 = corrs.slice(0, 15);
const partCorrs = top15.map(c => {
  const arr = haloProxies.map(h => h[c.name]);
  const rPartial = partialR(VfResid, arr, coreStructX);
  return { name: c.name, rRaw: c.r, rPartial, rPartialAbs: Math.abs(rPartial) };
}).sort((a, b) => b.rPartialAbs - a.rPartialAbs);

console.log('  ' + 'Proxy'.padEnd(22) + '  r(raw)  r(partial|core+struct)');
console.log('  ' + '─'.repeat(55));
for (const pc of partCorrs) {
  const flag = pc.rPartialAbs > 0.3 ? '★★' : pc.rPartialAbs > 0.15 ? '★ ' : '  ';
  console.log('  ' + flag + ' ' + pc.name.padEnd(22) + pc.rRaw.toFixed(3).padStart(7) + '  ' + pc.rPartial.toFixed(3).padStart(7));
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 3: PREDICTIVE DEATH MATCH — Core + candidate vs Core + VfResid');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const coreGap = gapPct(looCV(Y, core3X), sdY);
const vfResidX = gals45.map((_, i) => [...core3X[i], VfResid[i]]);
const vfResidGap = gapPct(looCV(Y, vfResidX), sdY);
const vflatX = gals45.map((_, i) => [...core3X[i], logVflat[i]]);
const vflatGap = gapPct(looCV(Y, vflatX), sdY);

console.log('  BENCHMARKS:');
console.log('    Core only:         ' + coreGap.toFixed(1) + '% gap');
console.log('    Core + Vflat:      ' + vflatGap.toFixed(1) + '% gap');
console.log('    Core + VfResid:    ' + vfResidGap.toFixed(1) + '% gap\n');

const candidateResults = [];
const bestProxies = corrs.filter(c => c.rAbs > 0.15).slice(0, 20);

console.log('  ' + 'Candidate'.padEnd(24) + '  gap%   delta  nested  flip%  r(VfR)');
console.log('  ' + '─'.repeat(70));

for (const bp of bestProxies) {
  const arr = haloProxies.map(h => h[bp.name]);
  const candX = gals45.map((_, i) => [...core3X[i], arr[i]]);
  const gap = gapPct(looCV(Y, candX), sdY);
  const nested = nestedCV(Y, core3X, candX);
  const flip = bootstrap(Y, candX, 500);
  const delta = gap - coreGap;

  candidateResults.push({ name: bp.name, gap, delta, nested, flip, rVfResid: bp.r });

  const flag = gap >= vfResidGap - 3 ? '◆' : gap >= vflatGap - 3 ? '●' : gap >= coreGap + 3 ? '○' : ' ';
  console.log('  ' + flag + ' ' + bp.name.padEnd(24) + gap.toFixed(1).padStart(5) + '%  ' +
    (delta > 0 ? '+' : '') + delta.toFixed(1).padStart(5) + 'pp  ' + nested + '/' + N + '   ' +
    (flip * 100).toFixed(1).padStart(5) + '%  ' + bp.r.toFixed(3));
}
console.log('');
console.log('  Legend: ◆ = within 3pp of VfResid, ● = within 3pp of Vflat, ○ = above core+3pp\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 4: BEST HALO BUNDLES vs VfResid');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const topCandNames = candidateResults.sort((a, b) => b.gap - a.gap).slice(0, 6).map(c => c.name);
console.log('  Top 6 candidates: ' + topCandNames.join(', ') + '\n');

const bundleTests = [];

for (let i = 0; i < Math.min(topCandNames.length, 5); i++) {
  for (let j = i + 1; j < Math.min(topCandNames.length, 6); j++) {
    const arr1 = haloProxies.map(h => h[topCandNames[i]]);
    const arr2 = haloProxies.map(h => h[topCandNames[j]]);
    const rBetween = Math.abs(pearsonR(arr1, arr2));
    if (rBetween > 0.85) continue;
    const bX = gals45.map((_, k) => [...core3X[k], arr1[k], arr2[k]]);
    const gap = gapPct(looCV(Y, bX), sdY);
    const nested = nestedCV(Y, core3X, bX);
    bundleTests.push({ name: topCandNames[i] + ' + ' + topCandNames[j], gap, nested, r12: rBetween });
  }
}

bundleTests.sort((a, b) => b.gap - a.gap);
console.log('  ' + 'Bundle (Core + ...)'.padEnd(50) + '  gap%   nested  r12');
console.log('  ' + '─'.repeat(75));
console.log('  ' + 'VfResid (target)'.padEnd(50) + vfResidGap.toFixed(1).padStart(6) + '%   ' + nestedCV(Y, core3X, vfResidX) + '/' + N);

for (const bt of bundleTests.slice(0, 10)) {
  console.log('  ' + bt.name.padEnd(50) + bt.gap.toFixed(1).padStart(6) + '%   ' + bt.nested + '/' + N + '    ' + bt.r12.toFixed(2));
}
console.log('');

if (topCandNames.length >= 3) {
  const tripleTests = [];
  for (let i = 0; i < Math.min(topCandNames.length, 4); i++) {
    for (let j = i + 1; j < Math.min(topCandNames.length, 5); j++) {
      for (let l = j + 1; l < Math.min(topCandNames.length, 6); l++) {
        const arr1 = haloProxies.map(h => h[topCandNames[i]]);
        const arr2 = haloProxies.map(h => h[topCandNames[j]]);
        const arr3 = haloProxies.map(h => h[topCandNames[l]]);
        const tX = gals45.map((_, k) => [...core3X[k], arr1[k], arr2[k], arr3[k]]);
        const gap = gapPct(looCV(Y, tX), sdY);
        tripleTests.push({ name: topCandNames[i] + ' + ' + topCandNames[j] + ' + ' + topCandNames[l], gap });
      }
    }
  }
  tripleTests.sort((a, b) => b.gap - a.gap);
  console.log('  Best triple bundles:');
  for (const tt of tripleTests.slice(0, 5)) {
    console.log('    Core + ' + tt.name.padEnd(55) + tt.gap.toFixed(1) + '%');
  }
  console.log('');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 5: WHAT DOES VfResid KNOW THAT THE BEST HALO PROXY DOESN\'T?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const bestSingle = candidateResults.sort((a, b) => b.gap - a.gap)[0];
const bestArr = haloProxies.map(h => h[bestSingle.name]);
const rVfResidBest = pearsonR(VfResid, bestArr);

console.log('  Best single halo proxy: ' + bestSingle.name + ' (gap=' + bestSingle.gap.toFixed(1) + '%)');
console.log('  r(VfResid, ' + bestSingle.name + ') = ' + rVfResidBest.toFixed(3) + '\n');

const vfResidAfterBest = ols(VfResid, gals45.map((_, i) => [bestArr[i]])).resid;
console.log('  VfResid after removing ' + bestSingle.name + ':');
console.log('    Variance remaining: ' + (sd(vfResidAfterBest) ** 2 / sd(VfResid) ** 2 * 100).toFixed(1) + '%');
console.log('    r(residual, logA0) = ' + pearsonR(vfResidAfterBest, Y).toFixed(3));

const coreVfResidAfterBest = gals45.map((_, i) => [...core3X[i], vfResidAfterBest[i]]);
const gapAfter = gapPct(looCV(Y, coreVfResidAfterBest), sdY);
console.log('    Core + residual gap = ' + gapAfter.toFixed(1) + '% (still +' + (gapAfter - coreGap).toFixed(1) + 'pp above core)');
console.log('');

const bestAfterVfResid = ols(bestArr, gals45.map((_, i) => [VfResid[i]])).resid;
console.log('  ' + bestSingle.name + ' after removing VfResid:');
console.log('    Variance remaining: ' + (sd(bestAfterVfResid) ** 2 / sd(bestArr) ** 2 * 100).toFixed(1) + '%');
console.log('    r(residual, logA0) = ' + pearsonR(bestAfterVfResid, Y).toFixed(3));

const coreBestAfterVfResid = gals45.map((_, i) => [...core3X[i], bestAfterVfResid[i]]);
const gapBestAfter = gapPct(looCV(Y, coreBestAfterVfResid), sdY);
console.log('    Core + residual gap = ' + gapBestAfter.toFixed(1) + '% (' + (gapBestAfter - coreGap > 0 ? '+' : '') + (gapBestAfter - coreGap).toFixed(1) + 'pp above core)');
console.log('');

console.log('  ASYMMETRY TEST:');
const vfResidSurvives = gapAfter - coreGap > 3;
const bestSurvives = gapBestAfter - coreGap > 3;
if (vfResidSurvives && !bestSurvives) {
  console.log('    VfResid DOMINATES — survives after best halo, but best halo dies after VfResid');
} else if (!vfResidSurvives && bestSurvives) {
  console.log('    ' + bestSingle.name + ' DOMINATES — survives after VfResid, but VfResid dies');
} else if (vfResidSurvives && bestSurvives) {
  console.log('    MUTUAL SURVIVAL — both carry independent information');
} else {
  console.log('    MUTUAL ABSORPTION — both largely captured by the other');
}
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 132A: VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const bestHaloGap = candidateResults.sort((a, b) => b.gap - a.gap)[0].gap;
const deficit = vfResidGap - bestHaloGap;
const bestBundleGap = bundleTests.length > 0 ? bundleTests[0].gap : bestHaloGap;
const bundleDeficit = vfResidGap - bestBundleGap;

let verdict;
if (deficit > 10) {
  verdict = 'VFRESID_UNTOUCHABLE';
} else if (deficit > 5) {
  verdict = 'VFRESID_DOMINANT';
} else if (deficit > 2) {
  verdict = 'VFRESID_BETTER';
} else if (deficit > -2) {
  verdict = 'TIE';
} else {
  verdict = 'HALO_PROXY_WINS';
}

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + verdict.padEnd(55) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════════╝\n');

console.log('  SCORECARD:');
console.log('    VfResid gap:           ' + vfResidGap.toFixed(1) + '%');
console.log('    Best halo single:      ' + bestSingle.name + ' = ' + bestHaloGap.toFixed(1) + '% (deficit ' + deficit.toFixed(1) + 'pp)');
console.log('    Best halo bundle:      ' + (bundleTests.length > 0 ? bundleTests[0].name : 'N/A') + ' = ' + bestBundleGap.toFixed(1) + '% (deficit ' + bundleDeficit.toFixed(1) + 'pp)');
console.log('    Model C (Vflat):       ' + vflatGap.toFixed(1) + '%');
console.log('    Core baseline:         ' + coreGap.toFixed(1) + '%');
console.log('');

if (verdict === 'VFRESID_UNTOUCHABLE' || verdict === 'VFRESID_DOMINANT') {
  console.log('  INTERPRETATION:');
  console.log('    No single halo proxy approaches VfResid predictive power.');
  console.log('    VfResid compresses halo-coupling information more efficiently');
  console.log('    than any individual observable or small bundle.');
  console.log('    The best halo proxy (' + bestSingle.name + ') captures part of');
  console.log('    the coupling physics but misses the compressed multi-scale');
  console.log('    dynamical integration that Vflat naturally performs.');
}
console.log('');

const output = {
  phase: '132A',
  title: 'Halo Driver Death Match',
  verdict,
  vfresid_gap: +vfResidGap.toFixed(1),
  vflat_gap: +vflatGap.toFixed(1),
  core_gap: +coreGap.toFixed(1),
  best_halo_single: { name: bestSingle.name, gap: +bestHaloGap.toFixed(1), deficit: +deficit.toFixed(1), r_vfresid: +rVfResidBest.toFixed(3) },
  best_halo_bundle: bundleTests.length > 0 ? { name: bundleTests[0].name, gap: +bestBundleGap.toFixed(1) } : null,
  ranked_correlations: corrs.slice(0, 20).map(c => ({ name: c.name, r: +c.r.toFixed(3) })),
  partial_correlations: partCorrs.map(pc => ({ name: pc.name, rRaw: +pc.rRaw.toFixed(3), rPartial: +pc.rPartial.toFixed(3) })),
  candidate_results: candidateResults.sort((a, b) => b.gap - a.gap).slice(0, 15).map(c => ({
    name: c.name, gap: +c.gap.toFixed(1), nested: c.nested + '/' + N, flip: +(c.flip * 100).toFixed(1)
  })),
  asymmetry: {
    vfresid_survives_after_best: vfResidSurvives,
    best_survives_after_vfresid: bestSurvives,
    vfresid_remaining_pct: +(sd(vfResidAfterBest) ** 2 / sd(VfResid) ** 2 * 100).toFixed(1),
    best_remaining_pct: +(sd(bestAfterVfResid) ** 2 / sd(bestArr) ** 2 * 100).toFixed(1)
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase132a-halo-death-match.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase132a-halo-death-match.json');

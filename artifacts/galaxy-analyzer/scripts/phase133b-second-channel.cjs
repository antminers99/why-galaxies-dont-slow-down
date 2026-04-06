#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 133B: SECOND HALO CHANNEL — lh_outerImprove');
console.log('');
console.log('  Is lh_outerImprove a genuine 5th axis beyond VfResid?');
console.log('  Full nested CV + bootstrap + fold-internal validation.');
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
    try {
      const f = ols(Yt, Xt);
      const xi = [1, ...X[i]];
      ss += (Y[i] - xi.reduce((s, x, j) => s + x * f.beta[j], 0)) ** 2;
    } catch (e) { ss += (Y[i] - mean(Y)) ** 2; }
  }
  return Math.sqrt(ss / n);
}
function nestedCV(Y, Xbase, Xfull) {
  const n = Y.length; let wins = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xbt = [...Xbase.slice(0, i), ...Xbase.slice(i + 1)];
    const Xft = [...Xfull.slice(0, i), ...Xfull.slice(i + 1)];
    try {
      const fb = ols(Yt, Xbt); const ff = ols(Yt, Xft);
      const pb = [1, ...Xbase[i]].reduce((s, x, j) => s + x * fb.beta[j], 0);
      const pf = [1, ...Xfull[i]].reduce((s, x, j) => s + x * ff.beta[j], 0);
      if (Math.abs(Y[i] - pf) <= Math.abs(Y[i] - pb)) wins++;
    } catch (e) { }
  }
  return wins;
}
function gapPct(rms, sdy) { return sdy > 0 ? 100 * (1 - rms ** 2 / sdy ** 2) : 0; }

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const logMHI = gals45.map(g => g.logMHI);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logMR = gals45.map(g => g.logMeanRun);
const logVflat = gals45.map(g => Math.log10(sparcMap[g.name].Vflat));
const logMbar = gals45.map(g => {
  const s = sparcMap[g.name];
  return Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
});
const logL36 = gals45.map(g => Math.log10(Math.max(sparcMap[g.name].L36, 0.01)));
const logRdisk = gals45.map(g => Math.log10(sparcMap[g.name].Rdisk));
const morphT = gals45.map(g => sparcMap[g.name].T);

const structX = gals45.map((_, i) => [logMbar[i], logL36[i], logRdisk[i], morphT[i]]);
const VfResid = ols(logVflat, structX).resid;

const lhOuter = gals45.map(g => resMap[g.name].models.log_halo.outerImprovement);
const lhImprove = gals45.map(g => resMap[g.name].models.log_halo.improvementVsNewton);
const lhInner = gals45.map(g => resMap[g.name].models.log_halo.innerImprovement);
const dhlMse = gals45.map(g => Math.log10(Math.max(resMap[g.name].models.dark_halo_linear.mse, 0.1)));
const mondImprove = gals45.map(g => resMap[g.name].models.mond.improvementVsNewton);

const core3X = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]);
const coreVfResidX = gals45.map((_, i) => [...core3X[i], VfResid[i]]);

const coreGap = gapPct(looCV(Y, core3X), sdY);
const vfResidGap = gapPct(looCV(Y, coreVfResidX), sdY);

console.log('  BASELINES: Core=' + coreGap.toFixed(1) + '%, Core+VfResid=' + vfResidGap.toFixed(1) + '%\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: 5th AXIS CANDIDATES — Core + VfResid + candidate');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const candidates5th = [
  { name: 'lh_outerImprove', arr: lhOuter },
  { name: 'lh_improve', arr: lhImprove },
  { name: 'lh_innerImprove', arr: lhInner },
  { name: 'dhl_mse', arr: dhlMse },
  { name: 'mond_improve', arr: mondImprove },
  { name: 'logSig0', arr: gals45.map(g => g.logSigma0) },
  { name: 'logVflat', arr: logVflat },
  { name: 'concentration', arr: gals45.map(g => Math.log10(sparcMap[g.name].Reff / sparcMap[g.name].Rdisk)) },
  { name: 'morphT', arr: morphT },
  { name: 'envCode', arr: gals45.map(g => g.envCode) }
];

console.log('  ' + 'Candidate'.padEnd(20) + '  gap%   delta  nested  flip%  r(cand,VfR)');
console.log('  ' + '─'.repeat(70));

const cand5Results = [];
for (const c of candidates5th) {
  const fullX = gals45.map((_, i) => [...coreVfResidX[i], c.arr[i]]);
  const gap5 = gapPct(looCV(Y, fullX), sdY);
  const nested = nestedCV(Y, coreVfResidX, fullX);
  const rCandVfR = pearsonR(c.arr, VfResid);

  let flipRate = 0;
  const nBoot = 500;
  const refBeta = ols(Y, fullX).beta;
  const refSigns = refBeta.slice(1).map(b => Math.sign(b));
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: N }, () => Math.floor(Math.random() * N));
    try {
      const fb = ols(idx.map(i => Y[i]), idx.map(i => fullX[i]));
      if (fb.beta.slice(1).some((b, j) => Math.sign(b) !== refSigns[j])) flipRate++;
    } catch (e) { flipRate++; }
  }
  flipRate /= nBoot;

  const delta = gap5 - vfResidGap;
  cand5Results.push({
    name: c.name, gap5: +gap5.toFixed(1), delta: +delta.toFixed(1),
    nested: nested + '/' + N, flipRate: +(flipRate * 100).toFixed(1),
    rCandVfR: +rCandVfR.toFixed(3)
  });

  const flag = delta > 3 && nested >= 25 && flipRate < 0.15 ? '◆' :
    delta > 2 && nested >= 23 ? '●' : delta > 0 ? '○' : ' ';
  console.log('  ' + flag + ' ' + c.name.padEnd(20) + gap5.toFixed(1).padStart(5) + '%  ' +
    (delta > 0 ? '+' : '') + delta.toFixed(1).padStart(5) + 'pp  ' + nested + '/' + N + '   ' +
    (flipRate * 100).toFixed(1).padStart(5) + '%  ' + rCandVfR.toFixed(3));
}
console.log('\n  ◆ = strong 5th axis  ● = moderate  ○ = positive but weak\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: FOLD-INTERNAL VALIDATION — lh_outerImprove');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

let foldInternalGap4 = 0, foldInternalGap5 = 0;
{
  let ss4 = 0, ss5 = 0;
  for (let i = 0; i < N; i++) {
    const trainIdx = [...Array(i), ...Array(N - i - 1).fill(0).map((_, j) => j < i ? j : j + 1 + i)].filter((_, j) => j < N - 1);
    const trIdx = gals45.map((_, j) => j).filter(j => j !== i);

    const logVflat_tr = trIdx.map(j => logVflat[j]);
    const structX_tr = trIdx.map(j => structX[j]);
    const structModel_tr = ols(logVflat_tr, structX_tr);

    const VfResid_tr = logVflat_tr.map((v, j) => v - [1, ...structX_tr[j]].reduce((s, x, k) => s + x * structModel_tr.beta[k], 0));
    const VfResid_i = logVflat[i] - [1, ...structX[i]].reduce((s, x, k) => s + x * structModel_tr.beta[k], 0);

    const Y_tr = trIdx.map(j => Y[j]);
    const core_tr = trIdx.map(j => core3X[j]);
    const coreVfR_tr = trIdx.map((_, j) => [...core_tr[j], VfResid_tr[j]]);
    const coreVfR5_tr = trIdx.map((_, j) => [...core_tr[j], VfResid_tr[j], lhOuter[trIdx[j]]]);

    const f4 = ols(Y_tr, coreVfR_tr);
    const f5 = ols(Y_tr, coreVfR5_tr);

    const pred4 = [1, ...core3X[i], VfResid_i].reduce((s, x, j) => s + x * f4.beta[j], 0);
    const pred5 = [1, ...core3X[i], VfResid_i, lhOuter[i]].reduce((s, x, j) => s + x * f5.beta[j], 0);

    ss4 += (Y[i] - pred4) ** 2;
    ss5 += (Y[i] - pred5) ** 2;
  }
  foldInternalGap4 = gapPct(Math.sqrt(ss4 / N), sdY);
  foldInternalGap5 = gapPct(Math.sqrt(ss5 / N), sdY);
}

console.log('  Fold-internal (VfResid computed inside each fold):');
console.log('    Core + VfResid (4-axis):               ' + foldInternalGap4.toFixed(1) + '% gap');
console.log('    Core + VfResid + lh_outerImprove (5-axis): ' + foldInternalGap5.toFixed(1) + '% gap');
console.log('    Delta: ' + (foldInternalGap5 - foldInternalGap4 > 0 ? '+' : '') + (foldInternalGap5 - foldInternalGap4).toFixed(1) + 'pp');
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: WHAT IS lh_outerImprove PHYSICALLY?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const physCorrelates = [
  { name: 'logVflat', arr: logVflat },
  { name: 'logMbar', arr: logMbar },
  { name: 'logMHI', arr: logMHI.map(v => v) },
  { name: 'logRdisk', arr: logRdisk },
  { name: 'morphT', arr: morphT },
  { name: 'logSBeff', arr: gals45.map(g => Math.log10(Math.max(sparcMap[g.name].SBeff, 0.01))) },
  { name: 'logExtent', arr: gals45.map(g => Math.log10(resMap[g.name].maxR / sparcMap[g.name].Rdisk)) },
  { name: 'pointCount', arr: gals45.map(g => resMap[g.name].pointCount) },
  { name: 'inc', arr: gals45.map(g => sparcMap[g.name].inc) },
  { name: 'D', arr: gals45.map(g => sparcMap[g.name].D) },
  { name: 'Q', arr: gals45.map(g => sparcMap[g.name].Q) },
  { name: 'haloK_linear', arr: gals45.map(g => Math.log10(Math.max(resMap[g.name].models.dark_halo_linear.k, 1))) },
  { name: 'VfResid', arr: VfResid },
  { name: 'envCode', arr: gals45.map(g => g.envCode) },
  { name: 'logMhost', arr: logMhost },
  { name: 'logMeanRun', arr: gals45.map(g => g.logMeanRun) }
];

console.log('  ' + 'Correlate'.padEnd(18) + '  r(lhOI)  r_partial(|core+VfR)');
console.log('  ' + '─'.repeat(50));

const lhOuterResid = ols(lhOuter, coreVfResidX.map(x => x.slice(0))).resid;

for (const pc of physCorrelates) {
  const rRaw = pearsonR(lhOuter, pc.arr);
  const pcResid = ols(pc.arr, coreVfResidX.map(x => x.slice(0))).resid;
  const rPartial = pearsonR(lhOuterResid, pcResid);
  const flag = Math.abs(rPartial) > 0.3 ? '★★' : Math.abs(rPartial) > 0.15 ? '★ ' : '  ';
  console.log('  ' + flag + ' ' + pc.name.padEnd(18) + rRaw.toFixed(3).padStart(7) + '  ' + rPartial.toFixed(3).padStart(7));
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: SYSTEMATICS CHECK — Is lh_outerImprove clean?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const systematics = [
  { name: 'inc', arr: gals45.map(g => sparcMap[g.name].inc) },
  { name: 'D', arr: gals45.map(g => sparcMap[g.name].D) },
  { name: 'Q', arr: gals45.map(g => sparcMap[g.name].Q) },
  { name: 'eVflat/Vflat', arr: gals45.map(g => sparcMap[g.name].eVflat / sparcMap[g.name].Vflat) },
  { name: 'pointCount', arr: gals45.map(g => resMap[g.name].pointCount) }
];

let systemClean = true;
for (const s of systematics) {
  const r = pearsonR(lhOuter, s.arr);
  const rPartial = pearsonR(lhOuterResid, ols(s.arr, coreVfResidX.map(x => x.slice(0))).resid);
  const flag = Math.abs(rPartial) > 0.3 ? '⚠️' : '✓';
  if (Math.abs(rPartial) > 0.3) systemClean = false;
  console.log('  ' + flag + ' ' + s.name.padEnd(18) + 'r=' + r.toFixed(3) + '  r_partial=' + rPartial.toFixed(3));
}
console.log('\n  Systematics verdict: ' + (systemClean ? 'CLEAN' : 'CONTAMINATED — proceed with caution') + '\n');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 133B: VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const bestCand = cand5Results.sort((a, b) => b.gap5 - a.gap5)[0];
const foldDelta = foldInternalGap5 - foldInternalGap4;
const isGenuine5th = bestCand.delta > 2 && foldDelta > 1 && bestCand.flipRate < 15;
const verdict = isGenuine5th ? 'GENUINE_5TH_AXIS' :
  bestCand.delta > 2 ? 'SUGGESTIVE_NOT_STABLE' :
    bestCand.delta > 0 ? 'MARGINAL' : 'NO_5TH_AXIS';

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + verdict.padEnd(55) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════════╝\n');

console.log('  Best 5th axis candidate: ' + bestCand.name);
console.log('    LOO gap (5-axis): ' + bestCand.gap5 + '% (delta +' + bestCand.delta + 'pp)');
console.log('    Nested CV: ' + bestCand.nested);
console.log('    Flip rate: ' + bestCand.flipRate + '%');
console.log('    Fold-internal delta: ' + (foldDelta > 0 ? '+' : '') + foldDelta.toFixed(1) + 'pp');
console.log('    Systematics: ' + (systemClean ? 'CLEAN' : 'FLAGGED'));
console.log('');

const output = {
  phase: '133B',
  title: 'Second Halo Channel — 5th Axis Test',
  verdict,
  candidates: cand5Results,
  foldInternal: {
    gap4: +foldInternalGap4.toFixed(1),
    gap5: +foldInternalGap5.toFixed(1),
    delta: +foldDelta.toFixed(1)
  },
  systematicsClean: systemClean
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase133b-second-channel.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase133b-second-channel.json');

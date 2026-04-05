#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 132C: EXTERNAL VALIDATION / SOURCE-ROBUSTNESS TEST');
console.log('');
console.log('  Does the coupling signal survive outside our specific N=45 subset?');
console.log('  Tests: leave-source-out, quality splits, domain splits, bootstrap.');
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
function gapPct(rms, sdy) { return sdy > 0 ? 100 * (1 - rms ** 2 / sdy ** 2) : 0; }

const Y_all = gals45.map(g => g.logA0);
const sdY_all = sd(Y_all);
const logMHI_all = gals45.map(g => g.logMHI);
const logMhost_all = gals45.map(g => tdMap[g.name].logMhost);
const logMR_all = gals45.map(g => g.logMeanRun);
const logVflat_all = gals45.map(g => Math.log10(sparcMap[g.name].Vflat));
const logMbar_all = gals45.map(g => {
  const s = sparcMap[g.name];
  return Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
});
const logL36_all = gals45.map(g => Math.log10(Math.max(sparcMap[g.name].L36, 0.01)));
const logRdisk_all = gals45.map(g => Math.log10(sparcMap[g.name].Rdisk));
const morphT_all = gals45.map(g => sparcMap[g.name].T);

function computeVfResid(indices) {
  const logVflat_sub = indices.map(i => logVflat_all[i]);
  const structX_sub = indices.map(i => [logMbar_all[i], logL36_all[i], logRdisk_all[i], morphT_all[i]]);
  return ols(logVflat_sub, structX_sub).resid;
}

function evaluateSubset(indices, label) {
  const n = indices.length;
  if (n < 10) return null;
  const Y = indices.map(i => Y_all[i]);
  const sdY = sd(Y);
  if (sdY < 0.01) return null;

  const logMHI = indices.map(i => logMHI_all[i]);
  const logMhost = indices.map(i => logMhost_all[i]);
  const logMR = indices.map(i => logMR_all[i]);
  const logVflat = indices.map(i => logVflat_all[i]);
  const VfResid = computeVfResid(indices);

  const core3X = indices.map((_, j) => [logMHI[j], logMhost[j], logMR[j]]);
  const coreVflatX = indices.map((_, j) => [...core3X[j], logVflat[j]]);
  const coreVfResidX = indices.map((_, j) => [...core3X[j], VfResid[j]]);

  const coreGap = gapPct(looCV(Y, core3X), sdY);
  const vflatGap = gapPct(looCV(Y, coreVflatX), sdY);
  const vfResidGap = gapPct(looCV(Y, coreVfResidX), sdY);

  const coreBeta = ols(Y, core3X).beta;
  const vflatBeta = ols(Y, coreVflatX).beta;

  const rVfResidA0 = pearsonR(VfResid, Y);

  return {
    label, n, sdY: +sdY.toFixed(4),
    coreGap: +coreGap.toFixed(1),
    vflatGap: +vflatGap.toFixed(1),
    vfResidGap: +vfResidGap.toFixed(1),
    vfResidDelta: +(vfResidGap - coreGap).toFixed(1),
    rVfResidA0: +rVfResidA0.toFixed(3),
    coreSigns: coreBeta.slice(1).map(b => Math.sign(b) > 0 ? '+' : '-').join(''),
    vflatSign: Math.sign(vflatBeta[4]) > 0 ? '+' : '-'
  };
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: LEAVE-SOURCE-OUT — Drop each source, measure stability');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const sources = {};
gals45.forEach((g, i) => {
  const src = g.source || 'unknown';
  if (!sources[src]) sources[src] = [];
  sources[src].push(i);
});

console.log('  Sources:');
for (const [src, idx] of Object.entries(sources)) {
  console.log('    ' + src + ': N=' + idx.length);
}
console.log('');

const fullResult = evaluateSubset(gals45.map((_, i) => i), 'Full N=45');
console.log('  ' + 'Split'.padEnd(28) + '  N   core%  Vflat%  VfRes%  delta  signs  r(VfR,a0)');
console.log('  ' + '─'.repeat(85));
console.log('  ' + 'FULL SAMPLE'.padEnd(28) + fullResult.n.toString().padStart(3) + '   ' +
  fullResult.coreGap.toFixed(1).padStart(5) + '  ' + fullResult.vflatGap.toFixed(1).padStart(6) + '  ' +
  fullResult.vfResidGap.toFixed(1).padStart(6) + '  ' + ('+' + fullResult.vfResidDelta.toFixed(1)).padStart(5) + '  ' +
  fullResult.coreSigns + fullResult.vflatSign + '    ' + fullResult.rVfResidA0.toFixed(3));

const lsoResults = [];
for (const [src, idx] of Object.entries(sources)) {
  if (idx.length < 3) continue;
  const remaining = gals45.map((_, i) => i).filter(i => !idx.includes(i));
  const r = evaluateSubset(remaining, 'drop-' + src);
  if (r) {
    lsoResults.push(r);
    console.log('  ' + r.label.padEnd(28) + r.n.toString().padStart(3) + '   ' +
      r.coreGap.toFixed(1).padStart(5) + '  ' + r.vflatGap.toFixed(1).padStart(6) + '  ' +
      r.vfResidGap.toFixed(1).padStart(6) + '  ' + (r.vfResidDelta > 0 ? '+' : '') + r.vfResidDelta.toFixed(1).padStart(5) + '  ' +
      r.coreSigns + r.vflatSign + '    ' + r.rVfResidA0.toFixed(3));
  }
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: QUALITY-BASED SPLITS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const qualityScore = gals45.map((g, i) => {
  const s = sparcMap[g.name];
  const r = resMap[g.name];
  return {
    idx: i,
    Q: s.Q,
    eVflat: s.eVflat / s.Vflat,
    eD: s.eD / s.D,
    inc: s.inc,
    pointCount: r.pointCount
  };
});

const highQ = qualityScore.filter(q => q.Q <= 2 && q.inc >= 30).map(q => q.idx);
const lowQ = qualityScore.filter(q => q.Q > 2 || q.inc < 30).map(q => q.idx);

const manyPts = qualityScore.sort((a, b) => b.pointCount - a.pointCount).slice(0, Math.floor(N / 2)).map(q => q.idx);
const fewPts = qualityScore.sort((a, b) => a.pointCount - b.pointCount).slice(0, Math.floor(N / 2)).map(q => q.idx);

const lowErrV = qualityScore.sort((a, b) => a.eVflat - b.eVflat).slice(0, Math.floor(N / 2)).map(q => q.idx);
const highErrV = qualityScore.sort((a, b) => b.eVflat - a.eVflat).slice(0, Math.floor(N / 2)).map(q => q.idx);

const splits = [
  ['High-Q (Q<=2, inc>=30)', highQ],
  ['Low-Q (Q>2 or inc<30)', lowQ],
  ['Many points (top half)', manyPts],
  ['Few points (bottom half)', fewPts],
  ['Low eVflat (top half)', lowErrV],
  ['High eVflat (bottom half)', highErrV]
];

console.log('  ' + 'Split'.padEnd(28) + '  N   core%  Vflat%  VfRes%  delta  signs  r(VfR,a0)');
console.log('  ' + '─'.repeat(85));

const splitResults = [];
for (const [label, idx] of splits) {
  const r = evaluateSubset(idx, label);
  if (r) {
    splitResults.push(r);
    console.log('  ' + r.label.padEnd(28) + r.n.toString().padStart(3) + '   ' +
      r.coreGap.toFixed(1).padStart(5) + '  ' + r.vflatGap.toFixed(1).padStart(6) + '  ' +
      r.vfResidGap.toFixed(1).padStart(6) + '  ' + (r.vfResidDelta > 0 ? '+' : '') + r.vfResidDelta.toFixed(1).padStart(5) + '  ' +
      r.coreSigns + r.vflatSign + '    ' + r.rVfResidA0.toFixed(3));
  }
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: DOMAIN SPLITS — Physical property divisions');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const medVflat = [...logVflat_all].sort((a, b) => a - b)[Math.floor(N / 2)];
const medMbar = [...logMbar_all].sort((a, b) => a - b)[Math.floor(N / 2)];
const medMHI = [...logMHI_all].sort((a, b) => a - b)[Math.floor(N / 2)];
const medT = [...morphT_all].sort((a, b) => a - b)[Math.floor(N / 2)];

const domainSplits = [
  ['High Vflat', gals45.map((_, i) => i).filter(i => logVflat_all[i] >= medVflat)],
  ['Low Vflat', gals45.map((_, i) => i).filter(i => logVflat_all[i] < medVflat)],
  ['High Mbar', gals45.map((_, i) => i).filter(i => logMbar_all[i] >= medMbar)],
  ['Low Mbar', gals45.map((_, i) => i).filter(i => logMbar_all[i] < medMbar)],
  ['High MHI', gals45.map((_, i) => i).filter(i => logMHI_all[i] >= medMHI)],
  ['Low MHI', gals45.map((_, i) => i).filter(i => logMHI_all[i] < medMHI)],
  ['Early type (T<=' + medT + ')', gals45.map((_, i) => i).filter(i => morphT_all[i] <= medT)],
  ['Late type (T>' + medT + ')', gals45.map((_, i) => i).filter(i => morphT_all[i] > medT)],
  ['Field (env>=3)', gals45.map((_, i) => i).filter(i => gals45[i].envCode >= 3)],
  ['Group (env<3)', gals45.map((_, i) => i).filter(i => gals45[i].envCode < 3)]
];

console.log('  ' + 'Split'.padEnd(28) + '  N   core%  Vflat%  VfRes%  delta  signs  r(VfR,a0)');
console.log('  ' + '─'.repeat(85));

const domainResults = [];
for (const [label, idx] of domainSplits) {
  const r = evaluateSubset(idx, label);
  if (r) {
    domainResults.push(r);
    console.log('  ' + r.label.padEnd(28) + r.n.toString().padStart(3) + '   ' +
      r.coreGap.toFixed(1).padStart(5) + '  ' + r.vflatGap.toFixed(1).padStart(6) + '  ' +
      r.vfResidGap.toFixed(1).padStart(6) + '  ' + (r.vfResidDelta > 0 ? '+' : '') + r.vfResidDelta.toFixed(1).padStart(5) + '  ' +
      r.coreSigns + r.vflatSign + '    ' + r.rVfResidA0.toFixed(3));
  }
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: SOURCE-BALANCED BOOTSTRAP — 1000 draws');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const nBoot = 1000;
const bootDeltas = [];
const bootSignStable = [];
const bootVfResidR = [];

for (let b = 0; b < nBoot; b++) {
  const idx = Array.from({ length: N }, () => Math.floor(Math.random() * N));

  const Y = idx.map(i => Y_all[i]);
  const sdY = sd(Y);
  if (sdY < 0.01) continue;

  const logMHI = idx.map(i => logMHI_all[i]);
  const logMhost = idx.map(i => logMhost_all[i]);
  const logMR = idx.map(i => logMR_all[i]);
  const logVflat = idx.map(i => logVflat_all[i]);
  const logMbar = idx.map(i => logMbar_all[i]);
  const logL36 = idx.map(i => logL36_all[i]);
  const logRdisk = idx.map(i => logRdisk_all[i]);
  const morphT = idx.map(i => morphT_all[i]);

  try {
    const structXb = idx.map((_, j) => [logMbar[j], logL36[j], logRdisk[j], morphT[j]]);
    const VfResid = ols(logVflat, structXb).resid;

    const core3X = idx.map((_, j) => [logMHI[j], logMhost[j], logMR[j]]);
    const coreVfResidX = idx.map((_, j) => [...core3X[j], VfResid[j]]);

    const coreBeta = ols(Y, core3X);
    const vfResBeta = ols(Y, coreVfResidX);

    const coreGap = gapPct(looCV(Y, core3X), sdY);
    const vfResidGap = gapPct(looCV(Y, coreVfResidX), sdY);

    const delta = vfResidGap - coreGap;
    bootDeltas.push(delta);

    const signs = coreBeta.beta.slice(1).map(b => Math.sign(b));
    const refSigns = [+1, -1, +1];
    const signMatch = signs.every((s, j) => s === refSigns[j]);
    bootSignStable.push(signMatch ? 1 : 0);

    bootVfResidR.push(pearsonR(VfResid, Y));
  } catch (e) { }
}

bootDeltas.sort((a, b) => a - b);
bootVfResidR.sort((a, b) => a - b);

const p5 = bootDeltas[Math.floor(bootDeltas.length * 0.05)];
const p25 = bootDeltas[Math.floor(bootDeltas.length * 0.25)];
const p50 = bootDeltas[Math.floor(bootDeltas.length * 0.5)];
const p75 = bootDeltas[Math.floor(bootDeltas.length * 0.75)];
const p95 = bootDeltas[Math.floor(bootDeltas.length * 0.95)];

console.log('  VfResid delta (VfResid gap - Core gap) distribution:');
console.log('    5th:  ' + p5.toFixed(1) + 'pp');
console.log('    25th: ' + p25.toFixed(1) + 'pp');
console.log('    50th: ' + p50.toFixed(1) + 'pp');
console.log('    75th: ' + p75.toFixed(1) + 'pp');
console.log('    95th: ' + p95.toFixed(1) + 'pp');
console.log('    P(delta > 0): ' + (bootDeltas.filter(d => d > 0).length / bootDeltas.length * 100).toFixed(1) + '%');
console.log('    P(delta > 5): ' + (bootDeltas.filter(d => d > 5).length / bootDeltas.length * 100).toFixed(1) + '%');
console.log('    P(delta > 10): ' + (bootDeltas.filter(d => d > 10).length / bootDeltas.length * 100).toFixed(1) + '%\n');

console.log('  Core sign stability:');
console.log('    Expected signs: MHI(+) Mhost(-) MR(+)');
console.log('    Stable: ' + (bootSignStable.reduce((s, v) => s + v, 0) / bootSignStable.length * 100).toFixed(1) + '%\n');

const r5 = bootVfResidR[Math.floor(bootVfResidR.length * 0.05)];
const r50 = bootVfResidR[Math.floor(bootVfResidR.length * 0.5)];
const r95 = bootVfResidR[Math.floor(bootVfResidR.length * 0.95)];
console.log('  r(VfResid, logA0) bootstrap:');
console.log('    5th: ' + r5.toFixed(3) + '  50th: ' + r50.toFixed(3) + '  95th: ' + r95.toFixed(3));
console.log('    P(r > 0.5): ' + (bootVfResidR.filter(r => r > 0.5).length / bootVfResidR.length * 100).toFixed(1) + '%\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 5: CROSS-DOMAIN TRANSFER — Train on one half, predict other');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function crossDomainTest(trainIdx, testIdx, label) {
  const n_train = trainIdx.length, n_test = testIdx.length;
  if (n_train < 8 || n_test < 8) return null;

  const Y_tr = trainIdx.map(i => Y_all[i]);
  const Y_te = testIdx.map(i => Y_all[i]);
  const sdY_te = sd(Y_te);
  if (sdY_te < 0.01) return null;

  const logVflat_tr = trainIdx.map(i => logVflat_all[i]);
  const structX_tr = trainIdx.map(i => [logMbar_all[i], logL36_all[i], logRdisk_all[i], morphT_all[i]]);
  const structModel = ols(logVflat_tr, structX_tr);

  const logVflat_te = testIdx.map(i => logVflat_all[i]);
  const structX_te = testIdx.map(i => [logMbar_all[i], logL36_all[i], logRdisk_all[i], morphT_all[i]]);
  const VfResid_te = logVflat_te.map((v, j) => v - [1, ...structX_te[j]].reduce((s, x, k) => s + x * structModel.beta[k], 0));

  const core_tr = trainIdx.map(i => [logMHI_all[i], logMhost_all[i], logMR_all[i]]);
  const core_te = testIdx.map(i => [logMHI_all[i], logMhost_all[i], logMR_all[i]]);

  const coreVfResid_tr = trainIdx.map((_, j) => {
    const logVflat_j = logVflat_tr[j];
    const struct_j = structX_tr[j];
    const pred_j = [1, ...struct_j].reduce((s, x, k) => s + x * structModel.beta[k], 0);
    const vfr_j = logVflat_j - pred_j;
    return [...core_tr[j], vfr_j];
  });

  const coreModel = ols(Y_tr, core_tr);
  const coreVfResidModel = ols(Y_tr, coreVfResid_tr);

  const corePred_te = core_te.map(x => [1, ...x].reduce((s, v, k) => s + v * coreModel.beta[k], 0));
  const coreVfResidPred_te = testIdx.map((_, j) => {
    const xFull = [1, ...core_te[j], VfResid_te[j]];
    return xFull.reduce((s, v, k) => s + v * coreVfResidModel.beta[k], 0);
  });

  const coreRMSE = Math.sqrt(Y_te.reduce((s, y, j) => s + (y - corePred_te[j]) ** 2, 0) / n_test);
  const vfResidRMSE = Math.sqrt(Y_te.reduce((s, y, j) => s + (y - coreVfResidPred_te[j]) ** 2, 0) / n_test);

  const coreGap = gapPct(coreRMSE, sdY_te);
  const vfResidGap = gapPct(vfResidRMSE, sdY_te);

  return { label, n_train, n_test, coreGap: +coreGap.toFixed(1), vfResidGap: +vfResidGap.toFixed(1), delta: +(vfResidGap - coreGap).toFixed(1) };
}

const halfN = Math.floor(N / 2);
const vflatSorted = gals45.map((_, i) => ({ i, v: logVflat_all[i] })).sort((a, b) => a.v - b.v);
const lowVflatIdx = vflatSorted.slice(0, halfN).map(x => x.i);
const highVflatIdx = vflatSorted.slice(halfN).map(x => x.i);

const mbarSorted = gals45.map((_, i) => ({ i, v: logMbar_all[i] })).sort((a, b) => a.v - b.v);
const lowMbarIdx = mbarSorted.slice(0, halfN).map(x => x.i);
const highMbarIdx = mbarSorted.slice(halfN).map(x => x.i);

const oddIdx = gals45.map((_, i) => i).filter(i => i % 2 === 1);
const evenIdx = gals45.map((_, i) => i).filter(i => i % 2 === 0);

const transferTests = [
  crossDomainTest(lowVflatIdx, highVflatIdx, 'Train: low-Vflat → Test: high-Vflat'),
  crossDomainTest(highVflatIdx, lowVflatIdx, 'Train: high-Vflat → Test: low-Vflat'),
  crossDomainTest(lowMbarIdx, highMbarIdx, 'Train: low-Mbar → Test: high-Mbar'),
  crossDomainTest(highMbarIdx, lowMbarIdx, 'Train: high-Mbar → Test: low-Mbar'),
  crossDomainTest(oddIdx, evenIdx, 'Train: odd-index → Test: even-index'),
  crossDomainTest(evenIdx, oddIdx, 'Train: even-index → Test: odd-index')
];

console.log('  ' + 'Transfer'.padEnd(42) + 'Ntr  Nte  core%  VfRes%  delta');
console.log('  ' + '─'.repeat(75));
const transferResults = [];
for (const t of transferTests) {
  if (t) {
    transferResults.push(t);
    console.log('  ' + t.label.padEnd(42) + t.n_train.toString().padStart(3) + '  ' +
      t.n_test.toString().padStart(3) + '  ' + t.coreGap.toFixed(1).padStart(5) + '  ' +
      t.vfResidGap.toFixed(1).padStart(6) + '  ' + (t.delta > 0 ? '+' : '') + t.delta.toFixed(1).padStart(5));
  }
}
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 132C: VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const allDeltas = [...lsoResults, ...splitResults, ...domainResults].map(r => r.vfResidDelta);
const positiveDeltas = allDeltas.filter(d => d > 0).length;
const totalDeltas = allDeltas.length;
const signConsistency = positiveDeltas / totalDeltas;

const allSigns = [...lsoResults, ...splitResults, ...domainResults].map(r => r.coreSigns);
const canonicalSign = fullResult.coreSigns;
const signStable = allSigns.filter(s => s === canonicalSign).length / allSigns.length;

const transferPositive = transferResults.filter(t => t.delta > 0).length;
const bootP = bootDeltas.filter(d => d > 0).length / bootDeltas.length;

let verdict;
if (signConsistency > 0.8 && bootP > 0.9 && transferPositive >= transferResults.length * 0.6) {
  verdict = 'ROBUST_CONFIRMED';
} else if (signConsistency > 0.6 && bootP > 0.8) {
  verdict = 'MOSTLY_ROBUST';
} else if (signConsistency > 0.4) {
  verdict = 'PARTIALLY_ROBUST';
} else {
  verdict = 'FRAGILE';
}

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + verdict.padEnd(55) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════════╝\n');

console.log('  STABILITY SCORES:');
console.log('    VfResid delta > 0 in ' + positiveDeltas + '/' + totalDeltas + ' splits (' + (signConsistency * 100).toFixed(0) + '%)');
console.log('    Core sign stable in ' + (signStable * 100).toFixed(0) + '% of splits');
console.log('    Bootstrap P(delta>0): ' + (bootP * 100).toFixed(1) + '%');
console.log('    Transfer tests positive: ' + transferPositive + '/' + transferResults.length);
console.log('');

const output = {
  phase: '132C',
  title: 'External Validation / Source-Robustness Test',
  verdict,
  fullSample: fullResult,
  leaveSourceOut: lsoResults,
  qualitySplits: splitResults,
  domainSplits: domainResults,
  bootstrap: {
    nBoot: bootDeltas.length,
    p5: +p5.toFixed(1), p25: +p25.toFixed(1), p50: +p50.toFixed(1),
    p75: +p75.toFixed(1), p95: +p95.toFixed(1),
    pDeltaGt0: +(bootP * 100).toFixed(1),
    signStability: +(bootSignStable.reduce((s, v) => s + v, 0) / bootSignStable.length * 100).toFixed(1)
  },
  transferTests: transferResults,
  stability: {
    deltaPositivePct: +(signConsistency * 100).toFixed(0),
    coreSignStablePct: +(signStable * 100).toFixed(0),
    transferPositive: transferPositive + '/' + transferResults.length
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase132c-external-robustness.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase132c-external-robustness.json');

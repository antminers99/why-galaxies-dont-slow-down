#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const pub = p => path.join(__dirname, '..', 'public', p);

console.log('======================================================================');
console.log('  PHASE 200: EXTERNAL VALIDATION DATA ASSEMBLY');
console.log('  Build complete dataset for SPARC galaxies outside N=56 stageA');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(pub('stage-A-master-table.json'), 'utf8'));
const sparcTable = JSON.parse(fs.readFileSync(pub('sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(pub('sparc-results.json'), 'utf8'));
const tsData = JSON.parse(fs.readFileSync(pub('transition-scale.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(pub('phase58a2-tidal-expansion.json'), 'utf8')).tidalData;

const stageAnames = new Set(stageA.galaxies.map(g => g.name));
const sparcMap = {}; sparcTable.forEach(s => { sparcMap[s.name] = s; });
const resMap = {}; sparcResults.perGalaxy.forEach(g => { if (g.rawName) resMap[g.rawName] = g; });
const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : 0; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function pearsonR(x, y) {
  const n = x.length, mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0(pts) {
  let lo = 2.0, hi = 5.0;
  for (let s = 0; s < 150; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gb = Math.pow(10, p.log_g_bar);
      const pred1 = mcgaughRAR(gb, Math.pow(10, m1));
      const pred2 = mcgaughRAR(gb, Math.pow(10, m2));
      if (pred1 > 0) c1 += (p.log_g_obs - Math.log10(pred1)) ** 2;
      if (pred2 > 0) c2 += (p.log_g_obs - Math.log10(pred2)) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  const logA0 = (lo + hi) / 2, a0 = Math.pow(10, logA0);
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gb, a0);
    if (pred > 0) ss += (p.log_g_obs - Math.log10(pred)) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length), n: pts.length };
}

function computeLogMeanRun(pts, a0) {
  const sorted = [...pts].sort((a, b) => a.log_g_bar - b.log_g_bar);
  const rarResid = sorted.map(p => {
    const gb = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gb, a0);
    return p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10);
  });
  let currentSign = Math.sign(rarResid[0]);
  let runLen = 1;
  const runs = [];
  for (let i = 1; i < rarResid.length; i++) {
    const s = Math.sign(rarResid[i]);
    if (s === currentSign) { runLen++; }
    else { runs.push(runLen); runLen = 1; currentSign = s; }
  }
  runs.push(runLen);
  const meanRunLen = runs.reduce((s, v) => s + v, 0) / runs.length;
  return Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1);
}

function estimateHaloMass(Vflat, L36) {
  if (Vflat && Vflat > 0) return 11.5 + 3.0 * (Math.log10(Vflat) - 2.2);
  if (L36 && L36 > 0) return 11.5 + 1.0 * (Math.log10(L36) - 1.5);
  return 11.5;
}

const byGal = {};
tsData.plotPoints.forEach(p => {
  if (!byGal[p.g]) byGal[p.g] = [];
  byGal[p.g].push({ log_g_bar: p.x, log_g_obs: p.x + p.y });
});

const nameMap = {};
sparcTable.forEach(s => { nameMap[s.name.substring(0, 12)] = s.name; });

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 1: REBUILD TRAINING STRUCTURAL MODEL (N=45)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const pubGals = stageA.galaxies.filter(g => pubNames.has(g.name));

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
  return { beta };
}

const trainData = pubGals.map(g => {
  const s = sparcMap[g.name];
  return {
    logVflat: Math.log10(s.Vflat),
    logMbar: Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9),
    logL36: Math.log10(Math.max(s.L36, 0.01)),
    logRdisk: Math.log10(s.Rdisk),
    morphT: s.T
  };
}).filter(g => isFinite(g.logVflat) && isFinite(g.logMbar));

const structModel = ols(
  trainData.map(g => g.logVflat),
  trainData.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT])
);

console.log('  Structural model (logVflat ~ logMbar + logL36 + logRdisk + T):');
console.log('  beta = [' + structModel.beta.map(b => b.toFixed(4)).join(', ') + ']');
console.log('  (intercept, logMbar, logL36, logRdisk, morphT)\n');

const trainLogA0 = pubGals.map(g => g.logA0);
const trainMeanLogA0 = mean(trainLogA0);
const trainSDLogA0 = sd(trainLogA0);
console.log('  Training logA0: mean=' + trainMeanLogA0.toFixed(4) + ' sd=' + trainSDLogA0.toFixed(4));

const trainLogMbar = trainData.map(g => g.logMbar);
const trainLogVflat = trainData.map(g => g.logVflat);
console.log('  Training logMbar range: [' + Math.min(...trainLogMbar).toFixed(2) + ', ' + Math.max(...trainLogMbar).toFixed(2) + ']');
console.log('  Training logVflat range: [' + Math.min(...trainLogVflat).toFixed(3) + ', ' + Math.max(...trainLogVflat).toFixed(3) + ']\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 2: ASSEMBLE EXTERNAL GALAXY DATA');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const externalGalaxies = [];
const rejections = { noRAR: 0, tooFewPts: 0, noSPARC: 0, noVflat: 0, boundaryA0: 0, fitFail: 0 };

for (const [shortName, rarPts] of Object.entries(byGal)) {
  const fullName = nameMap[shortName];
  if (!fullName || stageAnames.has(fullName)) continue;

  const s = sparcMap[fullName];
  if (!s) { rejections.noSPARC++; continue; }
  if (!s.Vflat || s.Vflat <= 0) { rejections.noVflat++; continue; }
  if (!s.MHI || s.MHI <= 0 || !s.L36 || s.L36 <= 0) { rejections.noSPARC++; continue; }
  if (rarPts.length < 5) { rejections.tooFewPts++; continue; }

  const fit = fitA0(rarPts);
  if (fit.logA0 <= 2.3 || fit.logA0 >= 4.8) { rejections.boundaryA0++; continue; }
  if (!isFinite(fit.logA0)) { rejections.fitFail++; continue; }

  const logMHI = Math.log10(s.MHI);
  const Mbar = s.L36 * 0.5 * 1e9 + s.MHI * 1.33 * 1e9;
  const logMbar = Math.log10(Mbar);
  const logL36 = Math.log10(Math.max(s.L36, 0.01));
  const logRdisk = Math.log10(s.Rdisk);
  const logVflat = Math.log10(s.Vflat);

  const predLogVflat = [1, logMbar, logL36, logRdisk, s.T]
    .reduce((sum, x, j) => sum + x * structModel.beta[j], 0);
  const VfResid = logVflat - predLogVflat;

  const logMeanRun = computeLogMeanRun(rarPts, fit.a0);

  const logMhost = estimateHaloMass(s.Vflat, s.L36);

  const r = resMap[fullName];
  const lhOuter = r && r.models && r.models.log_halo ? r.models.log_halo.outerImprovement : null;
  const haloK = r && r.models && r.models.dark_halo_linear ? Math.log10(Math.max(r.models.dark_halo_linear.k, 1)) : null;

  const extrapolating = logMbar > Math.max(...trainLogMbar) * 1.01 ||
    logVflat > Math.max(...trainLogVflat) * 1.01 ||
    logMbar < Math.min(...trainLogMbar) * 0.99 ||
    logVflat < Math.min(...trainLogVflat) * 0.99;

  externalGalaxies.push({
    name: fullName,
    logA0: +fit.logA0.toFixed(4),
    logMHI: +logMHI.toFixed(4),
    logMbar: +logMbar.toFixed(4),
    logMhost: +logMhost.toFixed(4),
    logMeanRun: +logMeanRun.toFixed(4),
    VfResid: +VfResid.toFixed(4),
    lhOuter: lhOuter !== null ? +lhOuter.toFixed(2) : null,
    haloK: haloK !== null ? +haloK.toFixed(4) : null,
    Vflat: s.Vflat,
    logVflat: +logVflat.toFixed(4),
    logL36: +logL36.toFixed(4),
    logRdisk: +logRdisk.toFixed(4),
    morphT: s.T,
    Q: s.Q,
    inc: s.inc,
    nRARpts: rarPts.length,
    fitRMS: +fit.rms.toFixed(4),
    extrapolating
  });
}

console.log('  Assembled: ' + externalGalaxies.length + ' external galaxies');
console.log('  Rejections:');
for (const [reason, count] of Object.entries(rejections)) {
  if (count > 0) console.log('    ' + reason + ': ' + count);
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 3: SAMPLE OVERVIEW');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const ext = externalGalaxies;
const highV = ext.filter(g => g.Vflat >= 120);
const veryHighV = ext.filter(g => g.Vflat >= 180);
const goodQ = ext.filter(g => g.Q <= 2);
const extrapolators = ext.filter(g => g.extrapolating);

console.log('  Total assembled: ' + ext.length);
console.log('  By quality: Q=1: ' + ext.filter(g => g.Q === 1).length +
  '  Q=2: ' + ext.filter(g => g.Q === 2).length +
  '  Q=3: ' + ext.filter(g => g.Q === 3).length);
console.log('  By Vflat: <120: ' + ext.filter(g => g.Vflat < 120).length +
  '  120-180: ' + ext.filter(g => g.Vflat >= 120 && g.Vflat < 180).length +
  '  >=180: ' + veryHighV.length);
console.log('  Extrapolating: ' + extrapolators.length);
console.log('  RAR points: min=' + Math.min(...ext.map(g => g.nRARpts)) +
  ' median=' + ext.map(g => g.nRARpts).sort((a, b) => a - b)[Math.floor(ext.length / 2)] +
  ' max=' + Math.max(...ext.map(g => g.nRARpts)));

console.log('\n  logA0 distribution:');
console.log('    External: mean=' + mean(ext.map(g => g.logA0)).toFixed(3) +
  ' sd=' + sd(ext.map(g => g.logA0)).toFixed(3) +
  ' range=[' + Math.min(...ext.map(g => g.logA0)).toFixed(3) + ', ' + Math.max(...ext.map(g => g.logA0)).toFixed(3) + ']');
console.log('    Training: mean=' + trainMeanLogA0.toFixed(3) + ' sd=' + trainSDLogA0.toFixed(3));

console.log('\n  Per-galaxy summary:');
console.log('  ' + 'Galaxy'.padEnd(16) + 'Vf'.padStart(5) + ' Q' + ' nPts'.padStart(5) +
  ' logA0'.padStart(7) + ' VfRes'.padStart(7) + ' lhOI'.padStart(6) + '  Extrap');
console.log('  ' + '─'.repeat(60));
ext.sort((a, b) => b.Vflat - a.Vflat);
ext.forEach(g => {
  console.log('  ' + g.name.padEnd(16) +
    g.Vflat.toFixed(0).padStart(5) + ' ' + g.Q +
    String(g.nRARpts).padStart(5) +
    g.logA0.toFixed(3).padStart(7) +
    g.VfResid.toFixed(3).padStart(7) +
    (g.lhOuter !== null ? g.lhOuter.toFixed(1).padStart(6) : '   n/a') +
    (g.extrapolating ? '  YES' : ''));
});

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 4: INITIAL CORRELATIONS (PREVIEW ONLY)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const logA0s = ext.map(g => g.logA0);
const vars = [
  { name: 'VfResid', vals: ext.map(g => g.VfResid) },
  { name: 'logMHI', vals: ext.map(g => g.logMHI) },
  { name: 'logMbar', vals: ext.map(g => g.logMbar) },
  { name: 'logMhost', vals: ext.map(g => g.logMhost) },
  { name: 'logMeanRun', vals: ext.map(g => g.logMeanRun) },
  { name: 'logVflat', vals: ext.map(g => g.logVflat) },
  { name: 'haloK', vals: ext.map(g => g.haloK ?? 0) },
];

console.log('  Correlations with logA0 (full external sample, N=' + ext.length + '):');
vars.forEach(v => {
  const valid = v.vals.filter((_, i) => isFinite(v.vals[i]));
  const r = pearsonR(logA0s.filter((_, i) => isFinite(v.vals[i])), valid);
  console.log('    ' + v.name.padEnd(15) + 'r = ' + r.toFixed(3));
});

if (highV.length >= 5) {
  console.log('\n  Correlations with logA0 (Vflat>=120, N=' + highV.length + '):');
  const hvA0 = highV.map(g => g.logA0);
  [
    { name: 'VfResid', vals: highV.map(g => g.VfResid) },
    { name: 'logMHI', vals: highV.map(g => g.logMHI) },
    { name: 'logMeanRun', vals: highV.map(g => g.logMeanRun) },
    { name: 'haloK', vals: highV.map(g => g.haloK ?? 0) },
  ].forEach(v => {
    const r = pearsonR(hvA0, v.vals);
    console.log('    ' + v.name.padEnd(15) + 'r = ' + r.toFixed(3));
  });
}

if (veryHighV.length >= 4) {
  console.log('\n  Correlations with logA0 (Vflat>=180, N=' + veryHighV.length + '):');
  const vhA0 = veryHighV.map(g => g.logA0);
  [
    { name: 'VfResid', vals: veryHighV.map(g => g.VfResid) },
    { name: 'logMHI', vals: veryHighV.map(g => g.logMHI) },
    { name: 'logMeanRun', vals: veryHighV.map(g => g.logMeanRun) },
  ].forEach(v => {
    const r = pearsonR(vhA0, v.vals);
    console.log('    ' + v.name.padEnd(15) + 'r = ' + r.toFixed(3));
  });
}

console.log('\n  COMPARE training r(VfResid, logA0) = 0.801 on N=10 holdout');
console.log('  This preview is NOT the formal transfer test (Phase 201).\n');

const output = {
  phase: '200',
  title: 'External Validation Data Assembly',
  timestamp: new Date().toISOString(),
  trainingReference: {
    N: pubGals.length,
    meanLogA0: +trainMeanLogA0.toFixed(4),
    sdLogA0: +trainSDLogA0.toFixed(4),
    structModelBeta: structModel.beta.map(b => +b.toFixed(6)),
    logMbarRange: [+Math.min(...trainLogMbar).toFixed(2), +Math.max(...trainLogMbar).toFixed(2)],
    logVflatRange: [+Math.min(...trainLogVflat).toFixed(3), +Math.max(...trainLogVflat).toFixed(3)]
  },
  externalSample: {
    N: ext.length,
    byQuality: { Q1: ext.filter(g => g.Q === 1).length, Q2: ext.filter(g => g.Q === 2).length, Q3: ext.filter(g => g.Q === 3).length },
    byVflat: {
      below120: ext.filter(g => g.Vflat < 120).length,
      '120to180': ext.filter(g => g.Vflat >= 120 && g.Vflat < 180).length,
      above180: veryHighV.length
    },
    nExtrapolating: extrapolators.length,
    logA0Stats: {
      mean: +mean(ext.map(g => g.logA0)).toFixed(4),
      sd: +sd(ext.map(g => g.logA0)).toFixed(4)
    }
  },
  galaxies: ext.sort((a, b) => b.Vflat - a.Vflat),
  rejections,
  caveats: [
    'logA0 computed from transition-scale RAR plotPoints (fewer points than phase56 localProfile)',
    'logMhost estimated from Vflat-based formula (no group membership data for external galaxies)',
    'logMeanRun ordered by log_g_bar (not radius), may differ slightly from phase56 method',
    'VfResid uses N=45 training structural model applied out-of-sample',
    'Some galaxies have very few RAR points (min=5), a0 fit may be imprecise'
  ]
};

const outPath = path.join(__dirname, '..', 'public', 'phase200-external-dataset.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('Output: ' + outPath);
console.log('\n======================================================================');
console.log('  PHASE 200 COMPLETE');
console.log('  Next: Phase 201 — Blind prediction using frozen coefficients');
console.log('======================================================================');

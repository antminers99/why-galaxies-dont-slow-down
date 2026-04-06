#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const pub = p => path.join(__dirname, '..', 'public', p);

console.log('======================================================================');
console.log('  PHASE 300: SAMPLE SALVAGE — RECOVER UNUSED SPARC GALAXIES');
console.log('  Goal: Augment external sample from 60 unused SPARC galaxies');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(pub('stage-A-master-table.json'), 'utf8'));
const sparcTable = JSON.parse(fs.readFileSync(pub('sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(pub('sparc-results.json'), 'utf8'));
const tsData = JSON.parse(fs.readFileSync(pub('transition-scale.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(pub('phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const extDataset = JSON.parse(fs.readFileSync(pub('phase200-external-dataset.json'), 'utf8'));

const stageAnames = new Set(stageA.galaxies.map(g => g.name));
const extNames = new Set(extDataset.galaxies.map(g => g.name));
const allUsed = new Set([...stageAnames, ...extNames]);

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
  return { beta: solveLinear(XtX, XtY) };
}

const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const pubGals = stageA.galaxies.filter(g => pubNames.has(g.name));

const trainData = pubGals.map(g => {
  const s = sparcMap[g.name];
  return {
    logVflat: Math.log10(s.Vflat),
    logMbar: Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9),
    logL36: Math.log10(Math.max(s.L36, 0.01)),
    logRdisk: Math.log10(s.Rdisk),
    morphT: s.T ?? 5
  };
}).filter(g => isFinite(g.logVflat) && isFinite(g.logMbar));

const structModel = ols(
  trainData.map(g => g.logVflat),
  trainData.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT])
);

console.log('  Structural model beta: [' + structModel.beta.map(b => b.toFixed(4)).join(', ') + ']');

const trainLogMbar = trainData.map(g => g.logMbar);
const trainLogVflat = trainData.map(g => g.logVflat);

const byGal = {};
tsData.plotPoints.forEach(p => {
  if (!byGal[p.g]) byGal[p.g] = [];
  byGal[p.g].push({ log_g_bar: p.x, log_g_obs: p.x + p.y });
});

const nameMap = {};
sparcTable.forEach(s => { nameMap[s.name.substring(0, 12)] = s.name; });

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STRATEGY 1: RECOVER GALAXIES WITHOUT Vflat');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const noVflatGals = sparcTable.filter(g => !allUsed.has(g.name) && (!g.Vflat || g.Vflat <= 0));
console.log('  Galaxies without Vflat: ' + noVflatGals.length);

let vflatEstimated = 0;
const recoveredFromVflat = [];

for (const s of noVflatGals) {
  const r = resMap[s.name];
  if (!r || !r.maxV || r.maxV <= 0) continue;
  if (r.pointCount < 5) continue;

  const outerV = r.maxV;
  if (outerV < 20) continue;

  const shortName = Object.keys(byGal).find(k => nameMap[k] === s.name);
  if (!shortName || !byGal[shortName]) continue;
  const rarPts = byGal[shortName];
  if (rarPts.length < 3) continue;

  if (!s.MHI || s.MHI <= 0 || !s.L36 || s.L36 <= 0) continue;

  const fit = fitA0(rarPts);
  if (fit.logA0 <= 2.3 || fit.logA0 >= 4.8 || !isFinite(fit.logA0)) continue;

  const logMHI = Math.log10(s.MHI);
  const Mbar = s.L36 * 0.5 * 1e9 + s.MHI * 1.33 * 1e9;
  const logMbar = Math.log10(Mbar);
  const logL36 = Math.log10(Math.max(s.L36, 0.01));
  const logRdisk = Math.log10(s.Rdisk);
  const logVflat = Math.log10(outerV);

  const predLogVflat = [1, logMbar, logL36, logRdisk, s.T ?? 5]
    .reduce((sum, x, j) => sum + x * structModel.beta[j], 0);
  const VfResid = logVflat - predLogVflat;

  const logMeanRun = computeLogMeanRun(rarPts, fit.a0);
  const logMhost = 11.5 + 3.0 * (logVflat - Math.log10(100));

  const lhOuter = r.models && r.models.log_halo ? r.models.log_halo.outerImprovement : null;
  const haloK = r.models && r.models.dark_halo_linear ? Math.log10(Math.max(r.models.dark_halo_linear.k, 1)) : null;

  recoveredFromVflat.push({
    name: s.name,
    logA0: +fit.logA0.toFixed(4),
    logMHI: +logMHI.toFixed(4),
    logMbar: +logMbar.toFixed(4),
    logMhost: +logMhost.toFixed(4),
    logMeanRun: +logMeanRun.toFixed(4),
    VfResid: +VfResid.toFixed(4),
    lhOuter: lhOuter !== null ? +lhOuter.toFixed(2) : null,
    haloK: haloK !== null ? +haloK.toFixed(4) : null,
    Vflat: outerV,
    VflatSource: 'maxV_from_RC',
    logVflat: +logVflat.toFixed(4),
    logL36: +logL36.toFixed(4),
    logRdisk: +logRdisk.toFixed(4),
    morphT: s.T ?? 5,
    Q: s.Q,
    inc: s.inc,
    nRARpts: rarPts.length,
    fitRMS: +fit.rms.toFixed(4),
    salvageSource: 'noVflat_maxV'
  });
  vflatEstimated++;
}

console.log('  Recovered using maxV from rotation curve: ' + vflatEstimated);
if (recoveredFromVflat.length > 0) {
  console.log('  Recovered galaxies:');
  recoveredFromVflat.sort((a, b) => b.Vflat - a.Vflat).forEach(g =>
    console.log('    ' + g.name.padEnd(16) + ' V=' + g.Vflat.toFixed(1) + ' Q=' + g.Q + ' nPts=' + g.nRARpts)
  );
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STRATEGY 2: LOWER POINT THRESHOLD (5 → 3)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const recoveredFromPts = [];
const alreadyRecovered = new Set(recoveredFromVflat.map(g => g.name));

for (const [shortName, rarPts] of Object.entries(byGal)) {
  const fullName = nameMap[shortName];
  if (!fullName || allUsed.has(fullName) || alreadyRecovered.has(fullName)) continue;

  const s = sparcMap[fullName];
  if (!s) continue;
  if (!s.Vflat || s.Vflat <= 0) continue;
  if (!s.MHI || s.MHI <= 0 || !s.L36 || s.L36 <= 0) continue;
  if (rarPts.length < 3 || rarPts.length >= 5) continue;

  const fit = fitA0(rarPts);
  if (fit.logA0 <= 2.3 || fit.logA0 >= 4.8 || !isFinite(fit.logA0)) continue;

  const logMHI = Math.log10(s.MHI);
  const Mbar = s.L36 * 0.5 * 1e9 + s.MHI * 1.33 * 1e9;
  const logMbar = Math.log10(Mbar);
  const logL36 = Math.log10(Math.max(s.L36, 0.01));
  const logRdisk = Math.log10(s.Rdisk);
  const logVflat = Math.log10(s.Vflat);

  const predLogVflat = [1, logMbar, logL36, logRdisk, s.T ?? 5]
    .reduce((sum, x, j) => sum + x * structModel.beta[j], 0);
  const VfResid = logVflat - predLogVflat;
  const logMeanRun = computeLogMeanRun(rarPts, fit.a0);
  const logMhost = 11.5 + 3.0 * (logVflat - Math.log10(100));

  const r = resMap[fullName];
  const lhOuter = r && r.models && r.models.log_halo ? r.models.log_halo.outerImprovement : null;
  const haloK = r && r.models && r.models.dark_halo_linear ? Math.log10(Math.max(r.models.dark_halo_linear.k, 1)) : null;

  recoveredFromPts.push({
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
    VflatSource: 'sparc_catalog',
    logVflat: +logVflat.toFixed(4),
    logL36: +logL36.toFixed(4),
    logRdisk: +logRdisk.toFixed(4),
    morphT: s.T ?? 5,
    Q: s.Q,
    inc: s.inc,
    nRARpts: rarPts.length,
    fitRMS: +fit.rms.toFixed(4),
    salvageSource: 'lowered_threshold'
  });
}

console.log('  Recovered with 3-4 RAR points: ' + recoveredFromPts.length);
if (recoveredFromPts.length > 0) {
  recoveredFromPts.sort((a, b) => b.Vflat - a.Vflat).forEach(g =>
    console.log('    ' + g.name.padEnd(16) + ' V=' + g.Vflat.toFixed(1) + ' Q=' + g.Q + ' nPts=' + g.nRARpts)
  );
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STRATEGY 3: ADD 3 KNOWN HIGH-V GALAXIES');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const knownHighV = ['NGC3953', 'NGC3949', 'NGC4051'];
const recoveredHighV = [];
const allRecoveredNames = new Set([...alreadyRecovered, ...recoveredFromPts.map(g => g.name)]);

for (const gname of knownHighV) {
  if (allRecoveredNames.has(gname)) {
    console.log('  ' + gname + ': already recovered in previous strategy');
    continue;
  }

  const s = sparcMap[gname];
  if (!s) { console.log('  ' + gname + ': not in SPARC table'); continue; }

  const shortName = Object.keys(byGal).find(k => nameMap[k] === gname);
  if (!shortName || !byGal[shortName]) { console.log('  ' + gname + ': no RAR data'); continue; }

  const rarPts = byGal[shortName];
  if (rarPts.length < 3) { console.log('  ' + gname + ': only ' + rarPts.length + ' RAR points'); continue; }
  if (!s.MHI || s.MHI <= 0 || !s.L36 || s.L36 <= 0) { console.log('  ' + gname + ': missing MHI/L36'); continue; }

  const fit = fitA0(rarPts);
  if (fit.logA0 <= 2.3 || fit.logA0 >= 4.8 || !isFinite(fit.logA0)) {
    console.log('  ' + gname + ': a0 fit boundary (' + fit.logA0.toFixed(3) + ')');
    continue;
  }

  const Vflat = s.Vflat || 0;
  const r = resMap[gname];
  let usedV = Vflat;
  let vSource = 'sparc_catalog';
  if (usedV <= 0 && r && r.maxV > 0) {
    usedV = r.maxV;
    vSource = 'maxV_from_RC';
  }
  if (usedV <= 0) { console.log('  ' + gname + ': no velocity estimate'); continue; }

  const logMHI = Math.log10(s.MHI);
  const Mbar = s.L36 * 0.5 * 1e9 + s.MHI * 1.33 * 1e9;
  const logMbar = Math.log10(Mbar);
  const logL36 = Math.log10(Math.max(s.L36, 0.01));
  const logRdisk = Math.log10(s.Rdisk);
  const logVflat = Math.log10(usedV);

  const predLogVflat = [1, logMbar, logL36, logRdisk, s.T ?? 5]
    .reduce((sum, x, j) => sum + x * structModel.beta[j], 0);
  const VfResid = logVflat - predLogVflat;
  const logMeanRun = computeLogMeanRun(rarPts, fit.a0);
  const logMhost = 11.5 + 3.0 * (logVflat - Math.log10(100));

  const lhOuter = r && r.models && r.models.log_halo ? r.models.log_halo.outerImprovement : null;
  const haloK = r && r.models && r.models.dark_halo_linear ? Math.log10(Math.max(r.models.dark_halo_linear.k, 1)) : null;

  recoveredHighV.push({
    name: gname,
    logA0: +fit.logA0.toFixed(4),
    logMHI: +logMHI.toFixed(4),
    logMbar: +logMbar.toFixed(4),
    logMhost: +logMhost.toFixed(4),
    logMeanRun: +logMeanRun.toFixed(4),
    VfResid: +VfResid.toFixed(4),
    lhOuter: lhOuter !== null ? +lhOuter.toFixed(2) : null,
    haloK: haloK !== null ? +haloK.toFixed(4) : null,
    Vflat: usedV,
    VflatSource: vSource,
    logVflat: +logVflat.toFixed(4),
    logL36: +logL36.toFixed(4),
    logRdisk: +logRdisk.toFixed(4),
    morphT: s.T ?? 5,
    Q: s.Q,
    inc: s.inc,
    nRARpts: rarPts.length,
    fitRMS: +fit.rms.toFixed(4),
    salvageSource: 'known_highV'
  });
  console.log('  ' + gname + ': RECOVERED V=' + usedV.toFixed(1) + ' Q=' + s.Q + ' nPts=' + rarPts.length + ' (' + vSource + ')');
}

const allSalvaged = [...recoveredFromVflat, ...recoveredFromPts, ...recoveredHighV];

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  SALVAGE SUMMARY');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  Total salvaged: ' + allSalvaged.length);
console.log('    From noVflat (using maxV): ' + recoveredFromVflat.length);
console.log('    From lowered threshold (3-4 pts): ' + recoveredFromPts.length);
console.log('    Known high-V targets: ' + recoveredHighV.length);
console.log('');

const salvHighV = allSalvaged.filter(g => g.Vflat >= 120);
const salvVeryHighV = allSalvaged.filter(g => g.Vflat >= 180);
console.log('  Salvaged by Vflat:');
console.log('    Vflat >= 180: ' + salvVeryHighV.length);
console.log('    Vflat 120-180: ' + allSalvaged.filter(g => g.Vflat >= 120 && g.Vflat < 180).length);
console.log('    Vflat < 120: ' + allSalvaged.filter(g => g.Vflat < 120).length);

const augmented = [...extDataset.galaxies, ...allSalvaged];
console.log('\n  AUGMENTED EXTERNAL SAMPLE:');
console.log('    Previous: N=' + extDataset.galaxies.length);
console.log('    Salvaged: N=' + allSalvaged.length);
console.log('    Total: N=' + augmented.length);
console.log('    Vflat >= 120: ' + augmented.filter(g => g.Vflat >= 120).length +
  ' (was ' + extDataset.galaxies.filter(g => g.Vflat >= 120).length + ')');
console.log('    Vflat >= 180: ' + augmented.filter(g => g.Vflat >= 180).length +
  ' (was ' + extDataset.galaxies.filter(g => g.Vflat >= 180).length + ')');

if (allSalvaged.length > 0) {
  console.log('\n  All salvaged galaxies:');
  console.log('  ' + 'Galaxy'.padEnd(16) + 'Vf'.padStart(6) + ' Q' + ' nPts'.padStart(5) +
    ' logA0'.padStart(7) + ' VfRes'.padStart(7) + ' Source');
  console.log('  ' + '-'.repeat(60));
  allSalvaged.sort((a, b) => b.Vflat - a.Vflat).forEach(g => {
    console.log('  ' + g.name.padEnd(16) +
      g.Vflat.toFixed(1).padStart(6) + ' ' + g.Q +
      String(g.nRARpts).padStart(5) +
      g.logA0.toFixed(3).padStart(7) +
      g.VfResid.toFixed(3).padStart(7) +
      '  ' + g.salvageSource);
  });
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  QUICK VALIDATION: CORRELATIONS ON AUGMENTED SAMPLE');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const augLogA0 = augmented.map(g => g.logA0);
const augVfResid = augmented.map(g => g.VfResid);
console.log('  Full augmented (N=' + augmented.length + '):');
console.log('    r(VfResid, logA0) = ' + pearsonR(augLogA0, augVfResid).toFixed(3));

const augHV = augmented.filter(g => g.Vflat >= 120);
if (augHV.length >= 5) {
  console.log('  High-V augmented (N=' + augHV.length + '):');
  console.log('    r(VfResid, logA0) = ' + pearsonR(augHV.map(g => g.logA0), augHV.map(g => g.VfResid)).toFixed(3));
}

const augVHV = augmented.filter(g => g.Vflat >= 180);
if (augVHV.length >= 4) {
  console.log('  Very-high-V augmented (N=' + augVHV.length + '):');
  console.log('    r(VfResid, logA0) = ' + pearsonR(augVHV.map(g => g.logA0), augVHV.map(g => g.VfResid)).toFixed(3));
}

const prevLogA0 = extDataset.galaxies.map(g => g.logA0);
const prevVfResid = extDataset.galaxies.map(g => g.VfResid);
console.log('\n  Compare with previous N=59:');
console.log('    r(VfResid, logA0) = ' + pearsonR(prevLogA0, prevVfResid).toFixed(3));
const prevHV = extDataset.galaxies.filter(g => g.Vflat >= 120);
if (prevHV.length >= 5) {
  console.log('  Previous high-V (N=' + prevHV.length + '):');
  console.log('    r(VfResid, logA0) = ' + pearsonR(prevHV.map(g => g.logA0), prevHV.map(g => g.VfResid)).toFixed(3));
}

const output = {
  phase: '300',
  title: 'Sample Salvage — Recover Unused SPARC Galaxies',
  timestamp: new Date().toISOString(),
  strategies: {
    noVflat_maxV: { attempted: noVflatGals.length, recovered: recoveredFromVflat.length },
    lowered_threshold: { attempted: '3-4 RAR pts', recovered: recoveredFromPts.length },
    known_highV: { targets: knownHighV, recovered: recoveredHighV.length }
  },
  totalSalvaged: allSalvaged.length,
  salvageByVflat: {
    above180: salvVeryHighV.length,
    '120to180': allSalvaged.filter(g => g.Vflat >= 120 && g.Vflat < 180).length,
    below120: allSalvaged.filter(g => g.Vflat < 120).length
  },
  augmentedSample: {
    previous: extDataset.galaxies.length,
    salvaged: allSalvaged.length,
    total: augmented.length,
    highV_total: augmented.filter(g => g.Vflat >= 120).length,
    veryHighV_total: augmented.filter(g => g.Vflat >= 180).length
  },
  salvagedGalaxies: allSalvaged,
  augmentedGalaxies: augmented.sort((a, b) => b.Vflat - a.Vflat),
  caveats: [
    'noVflat galaxies use maxV from rotation curve as Vflat proxy — NOT catalog Vflat',
    'Lowered-threshold galaxies have only 3-4 RAR points — logA0 less precise',
    'Salvaged galaxies are lower-quality by construction — treat as bonus, not primary',
    'VfResid still uses frozen N=45 structural model',
    'logMhost still uses crude Vflat-based formula'
  ]
};

const outPath = pub('phase300-sample-salvage.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nOutput: ' + outPath);

console.log('\n======================================================================');
console.log('  PHASE 300 COMPLETE — Sample salvage finished');
console.log('  Next: Phase 301 — Model VfResid itself');
console.log('======================================================================');

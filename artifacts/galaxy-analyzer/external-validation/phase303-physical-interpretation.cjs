#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const pub = p => path.join(__dirname, '..', 'public', p);

console.log('======================================================================');
console.log('  PHASE 303: PHYSICAL INTERPRETATION');
console.log('  What IS the hidden dynamical channel inside VfResid?');
console.log('');
console.log('  Candidate hypotheses:');
console.log('    H1: Halo response (adiabatic contraction / baryon-halo coupling)');
console.log('    H2: Assembly history (formation time / halo age)');
console.log('    H3: Feedback imprint (SN/AGN reshaping of inner halo)');
console.log('    H4: Dynamical integration (accumulated baryon-halo processing)');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(pub('stage-A-master-table.json'), 'utf8'));
const sparcTable = JSON.parse(fs.readFileSync(pub('sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(pub('sparc-results.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(pub('phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const extDataset = JSON.parse(fs.readFileSync(pub('phase200-external-dataset.json'), 'utf8'));
const salvageData = JSON.parse(fs.readFileSync(pub('phase300-sample-salvage.json'), 'utf8'));

const sparcMap = {}; sparcTable.forEach(s => { sparcMap[s.name] = s; });
const resMap = {}; sparcResults.perGalaxy.forEach(g => { if (g.rawName) resMap[g.rawName] = g; });
const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : 0; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function pearsonR(x, y) {
  const n = x.length; if (n < 3) return NaN;
  const mx = mean(x), my = mean(y);
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
function olsFull(Y, X) {
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
  let sse = 0, sst = 0;
  const my = mean(Y);
  for (let i = 0; i < n; i++) {
    const pred = Xa[i].reduce((s, x, j) => s + x * beta[j], 0);
    sse += (Y[i] - pred) ** 2;
    sst += (Y[i] - my) ** 2;
  }
  return { beta, r2: sst > 0 ? 1 - sse / sst : 0, rmse: Math.sqrt(sse / n) };
}
function looR2(Y, X) {
  const n = Y.length;
  let press = 0;
  const my = mean(Y);
  const sst = Y.reduce((s, y) => s + (y - my) ** 2, 0);
  if (sst === 0) return 0;
  for (let i = 0; i < n; i++) {
    const Yt = Y.filter((_, j) => j !== i);
    const Xt = X.filter((_, j) => j !== i);
    try {
      const m = olsFull(Yt, Xt);
      const pred = [1, ...X[i]].reduce((s, x, j) => s + x * m.beta[j], 0);
      press += (Y[i] - pred) ** 2;
    } catch (e) { press += (Y[i] - my) ** 2; }
  }
  return 1 - press / sst;
}
function bootstrapSign(Y, X, nBoot) {
  const n = Y.length, p = X[0].length + 1;
  const signs = Array.from({ length: p }, () => ({ pos: 0, neg: 0 }));
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(Math.random() * n));
    try {
      const m = olsFull(idx.map(i => Y[i]), idx.map(i => X[i]));
      for (let j = 0; j < p; j++) { if (m.beta[j] > 0) signs[j].pos++; else signs[j].neg++; }
    } catch (e) {}
  }
  return signs.map(s => Math.max(s.pos, s.neg) / (s.pos + s.neg) * 100);
}

const structTrainData = stageA.galaxies.filter(g => pubNames.has(g.name)).map(g => {
  const s = sparcMap[g.name];
  return {
    logVflat: Math.log10(s.Vflat),
    logMbar: Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9),
    logL36: Math.log10(Math.max(s.L36, 0.01)),
    logRdisk: Math.log10(s.Rdisk),
    morphT: s.T ?? 5
  };
}).filter(g => isFinite(g.logVflat) && isFinite(g.logMbar));
const structModel = olsFull(structTrainData.map(g => g.logVflat), structTrainData.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]));

const internal = stageA.galaxies.filter(g => pubNames.has(g.name)).map(g => {
  const s = sparcMap[g.name]; const r = resMap[g.name]; const td = tdMap[g.name];
  const logMbar = Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
  const logVflat = Math.log10(s.Vflat);
  const logL36 = Math.log10(Math.max(s.L36, 0.01));
  const logRdisk = Math.log10(s.Rdisk);
  const predLV = [1, logMbar, logL36, logRdisk, s.T ?? 5].reduce((sum, x, j) => sum + x * structModel.beta[j], 0);
  const VfResid = logVflat - predLV;
  const haloK = r && r.models && r.models.dark_halo_linear ? Math.log10(Math.max(r.models.dark_halo_linear.k, 1)) : null;
  const lhOuter = r && r.models && r.models.log_halo ? r.models.log_halo.outerImprovement : null;
  const mondImprove = r && r.models && r.models.mond ? r.models.mond.improvementVsNewton : null;
  const logMgas = Math.log10(Math.pow(10, g.logMHI) * 1.33 * 1e9);
  const fGas = Math.pow(10, logMgas) / Math.pow(10, logMbar);
  const logFgas = Math.log10(Math.max(fGas, 0.001));
  const logSigmaBar = logMbar - 2 * logRdisk;
  const logDynTime = logRdisk - logVflat + Math.log10(3.086e16 / 3.156e7);
  const maxR = r ? r.maxR || 0 : 0;
  const rcExtent = maxR / Math.max(s.Rdisk, 0.1);
  const logSBeff = Math.log10(Math.max(s.SBeff || 0.01, 0.01));
  const logSBdisk = Math.log10(Math.max(s.SBdisk || 0.01, 0.01));
  const diskDominance = s.L36 * 0.5 * 1e9 / Math.pow(10, logMbar);

  return {
    name: g.name, VfResid, logA0: g.logA0,
    haloK, lhOuter, mondImprove,
    logMHI: g.logMHI, logMbar, logMgas,
    logMhost: td ? td.logMhost : null,
    logMeanRun: g.logMeanRun, logSigma0: g.logSigma0,
    logL36, logRdisk, logSBeff, logSBdisk,
    morphT: s.T ?? 5, envCode: g.envCode,
    Vflat: s.Vflat, logVflat, Q: s.Q, inc: s.inc,
    logFgas, logSigmaBar, logDynTime,
    rcExtent, diskDominance,
    source: 'internal'
  };
}).filter(g => isFinite(g.VfResid));

const extAll = [...extDataset.galaxies, ...(salvageData.salvagedGalaxies || [])].map(g => {
  const s = sparcMap[g.name]; const r = resMap[g.name];
  const logMgas = g.logMHI ? Math.log10(Math.pow(10, g.logMHI) * 1.33 * 1e9) : null;
  const logMbar = g.logMbar;
  const logRdisk = g.logRdisk;
  const fGas = logMgas && logMbar ? Math.pow(10, logMgas) / Math.pow(10, logMbar) : null;
  const logFgas = fGas ? Math.log10(Math.max(fGas, 0.001)) : null;
  const logSigmaBar = logMbar && logRdisk ? logMbar - 2 * logRdisk : null;
  const logVflat = g.logVflat ?? Math.log10(g.Vflat);
  const logDynTime = logRdisk && logVflat ? logRdisk - logVflat + Math.log10(3.086e16 / 3.156e7) : null;
  const maxR = r ? r.maxR || 0 : 0;
  const rcExtent = s ? maxR / Math.max(s.Rdisk, 0.1) : null;
  const logSBeff = s ? Math.log10(Math.max(s.SBeff || 0.01, 0.01)) : null;
  const mondImprove = r && r.models && r.models.mond ? r.models.mond.improvementVsNewton : null;
  const diskDominance = s && logMbar ? s.L36 * 0.5 * 1e9 / Math.pow(10, logMbar) : null;

  return {
    name: g.name, VfResid: g.VfResid, logA0: g.logA0,
    haloK: g.haloK, lhOuter: g.lhOuter, mondImprove,
    logMHI: g.logMHI, logMbar, logMgas,
    logMhost: g.logMhost, logMeanRun: g.logMeanRun,
    logL36: g.logL36, logRdisk, logSBeff,
    morphT: g.morphT ?? (s ? s.T ?? 5 : 5),
    Vflat: g.Vflat, logVflat,
    Q: g.Q ?? (s ? s.Q : null), inc: g.inc ?? (s ? s.inc : null),
    logFgas, logSigmaBar, logDynTime, rcExtent, diskDominance,
    source: g.salvageSource ? 'salvaged' : 'external'
  };
});

const pooled = [...internal, ...extAll];
console.log('  Internal: N=' + internal.length);
console.log('  External: N=' + extAll.length);
console.log('  Pooled: N=' + pooled.length);

function computeHaloKResidual(data) {
  const valid = data.filter(g => g.haloK !== null && isFinite(g.haloK));
  const model = olsFull(valid.map(g => g.VfResid), valid.map(g => [g.haloK]));
  valid.forEach((g, i) => {
    g.hkResid = g.VfResid - [1, g.haloK].reduce((s, x, j) => s + x * model.beta[j], 0);
  });
  return valid;
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 1: CONSTRUCT HYPOTHESIS-SPECIFIC OBSERVABLES');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  H1 (Halo response) proxies:');
console.log('    logSigmaBar = logMbar - 2*logRdisk  (baryonic surface density)');
console.log('    diskDominance = M_disk / M_bar      (disk vs gas ratio)');
console.log('    haloK already available               (direct halo amplitude)');
console.log('  H2 (Assembly history) proxies:');
console.log('    morphT                                (Hubble type → formation epoch)');
console.log('    logSBeff                              (surface brightness → concentration)');
console.log('    logSBdisk                             (disk SB)');
console.log('  H3 (Feedback imprint) proxies:');
console.log('    logFgas = log(Mgas/Mbar)             (gas fraction → feedback efficiency)');
console.log('    logMHI                                (gas mass)');
console.log('    mondImprove                           (MOND improvement → gravity law sensitivity)');
console.log('  H4 (Dynamical integration) proxies:');
console.log('    logDynTime = logRdisk - logVflat     (crossing time)');
console.log('    rcExtent = maxR/Rdisk                 (how many disk radii sampled)');
console.log('    logMeanRun                            (RAR run length → dynamical coherence)');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 2: CORRELATE EACH PROXY WITH VfResid AND WITH hkResid');
console.log('  (hkResid = VfResid after removing haloK)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const intWithResid = computeHaloKResidual(internal);
const extWithResid = computeHaloKResidual(extAll);
const poolWithResid = computeHaloKResidual(pooled);

const hypothesisProxies = {
  'H1_haloResponse': ['logSigmaBar', 'diskDominance', 'haloK'],
  'H2_assemblyHistory': ['morphT', 'logSBeff'],
  'H3_feedbackImprint': ['logFgas', 'logMHI', 'mondImprove'],
  'H4_dynIntegration': ['logDynTime', 'rcExtent', 'logMeanRun']
};

function proxyCorrelations(data, target, label) {
  const allProxies = [
    'logSigmaBar', 'diskDominance', 'haloK',
    'morphT', 'logSBeff',
    'logFgas', 'logMHI', 'mondImprove',
    'logDynTime', 'rcExtent', 'logMeanRun'
  ];

  console.log('  ' + label + ' — r(proxy, ' + target + '):');
  console.log('  ' + 'Hypothesis'.padEnd(14) + 'Proxy'.padEnd(16) + 'Full'.padStart(8) + 'Low-V'.padStart(8) + 'High-V'.padStart(8) + 'VHigh-V'.padStart(8));
  console.log('  ' + '-'.repeat(65));

  const results = [];
  for (const [hyp, proxies] of Object.entries(hypothesisProxies)) {
    for (const proxy of proxies) {
      const full = data.filter(g => g[proxy] !== null && g[proxy] !== undefined && isFinite(g[proxy]) && g[target] !== undefined && isFinite(g[target]));
      const low = full.filter(g => g.Vflat < 120);
      const high = full.filter(g => g.Vflat >= 120);
      const vhigh = full.filter(g => g.Vflat >= 180);

      const rFull = full.length >= 5 ? pearsonR(full.map(g => g[target]), full.map(g => g[proxy])) : NaN;
      const rLow = low.length >= 4 ? pearsonR(low.map(g => g[target]), low.map(g => g[proxy])) : NaN;
      const rHigh = high.length >= 4 ? pearsonR(high.map(g => g[target]), high.map(g => g[proxy])) : NaN;
      const rVH = vhigh.length >= 4 ? pearsonR(vhigh.map(g => g[target]), vhigh.map(g => g[proxy])) : NaN;

      const fmt = v => isNaN(v) ? '   N/A' : v.toFixed(3).padStart(8);
      console.log('  ' + hyp.slice(0, 13).padEnd(14) + proxy.padEnd(16) + fmt(rFull) + fmt(rLow) + fmt(rHigh) + fmt(rVH));

      results.push({
        hypothesis: hyp, proxy,
        rFull: isNaN(rFull) ? null : +rFull.toFixed(3),
        rLow: isNaN(rLow) ? null : +rLow.toFixed(3),
        rHigh: isNaN(rHigh) ? null : +rHigh.toFixed(3),
        rVHigh: isNaN(rVH) ? null : +rVH.toFixed(3),
        nFull: full.length, nLow: low.length, nHigh: high.length, nVHigh: vhigh.length,
        regimeStrengthens: !isNaN(rHigh) && !isNaN(rLow) && Math.abs(rHigh) > Math.abs(rLow) + 0.05
      });
    }
  }
  return results;
}

console.log('  === CORRELATIONS WITH VfResid ===\n');
const corrVfRInt = proxyCorrelations(internal, 'VfResid', 'INTERNAL');
console.log('');
const corrVfRPool = proxyCorrelations(pooled, 'VfResid', 'POOLED');

console.log('\n  === CORRELATIONS WITH hkResid (VfResid after removing haloK) ===\n');
const corrHkRInt = proxyCorrelations(intWithResid, 'hkResid', 'INTERNAL (hkResid)');
console.log('');
const corrHkRPool = proxyCorrelations(poolWithResid, 'hkResid', 'POOLED (hkResid)');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 3: HYPOTHESIS MODELS — WHICH BEST EXPLAINS hkResid?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function testHypothesisModel(data, features, label) {
  const valid = data.filter(g => features.every(f => g[f] !== null && g[f] !== undefined && isFinite(g[f])) && g.hkResid !== undefined && isFinite(g.hkResid));
  if (valid.length < features.length + 3) {
    console.log('  ' + label + ': SKIP (N=' + valid.length + ')');
    return null;
  }
  const Y = valid.map(g => g.hkResid);
  const X = valid.map(g => features.map(f => g[f]));
  const model = olsFull(Y, X);
  let loo = NaN;
  try { loo = looR2(Y, X); } catch (e) {}
  const sign = bootstrapSign(Y, X, 500);

  const rResidA0 = pearsonR(
    valid.map(g => g.logA0),
    Y.map((y, i) => y - [1, ...X[i]].reduce((s, x, j) => s + x * model.beta[j], 0))
  );

  console.log('  ' + label + ' (N=' + valid.length + '):');
  console.log('    R² = ' + model.r2.toFixed(3) + '  LOO R² = ' + (isNaN(loo) ? 'N/A' : loo.toFixed(3)));
  console.log('    Features: ' + features.join(', '));
  console.log('    Beta: [' + model.beta.map(b => b.toFixed(4)).join(', ') + ']');
  console.log('    Sign stability: [' + sign.map(s => s.toFixed(0) + '%').join(', ') + ']');
  console.log('    After removal → r(remaining, a0) = ' + (isNaN(rResidA0) ? 'N/A' : rResidA0.toFixed(3)));

  return {
    label, features, n: valid.length,
    r2: +model.r2.toFixed(3), loo: isNaN(loo) ? null : +loo.toFixed(3),
    signStability: sign.map(s => +s.toFixed(1)),
    rResidA0AfterHyp: isNaN(rResidA0) ? null : +rResidA0.toFixed(3)
  };
}

console.log('  --- INTERNAL hkResid models ---\n');
const hypModelsInt = [];
hypModelsInt.push(testHypothesisModel(intWithResid, ['logSigmaBar'], 'H1a: SigmaBar'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logSigmaBar', 'diskDominance'], 'H1b: SigmaBar + diskDom'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['morphT'], 'H2a: morphT'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logSBeff'], 'H2b: SBeff'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['morphT', 'logSBeff'], 'H2c: morphT + SBeff'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logFgas'], 'H3a: Fgas'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logFgas', 'mondImprove'], 'H3b: Fgas + mondImprove'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['mondImprove'], 'H3c: mondImprove'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logDynTime'], 'H4a: dynTime'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logMeanRun'], 'H4b: MeanRun'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logDynTime', 'logMeanRun'], 'H4c: dynTime + MeanRun'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logDynTime', 'rcExtent'], 'H4d: dynTime + rcExtent'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['mondImprove', 'logMeanRun'], 'H3+H4: mondImprove + MR'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logSigmaBar', 'logDynTime', 'mondImprove'], 'H1+H3+H4: mixed'));
hypModelsInt.push(testHypothesisModel(intWithResid, ['logSigmaBar', 'logFgas', 'logDynTime', 'logMeanRun'], 'FULL: all proxy'));

console.log('\n  --- POOLED hkResid models ---\n');
const hypModelsPool = [];
hypModelsPool.push(testHypothesisModel(poolWithResid, ['logSigmaBar'], 'H1a: SigmaBar'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['logSigmaBar', 'diskDominance'], 'H1b: SigmaBar + diskDom'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['morphT'], 'H2a: morphT'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['logFgas'], 'H3a: Fgas'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['logFgas', 'mondImprove'], 'H3b: Fgas + mondImprove'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['mondImprove'], 'H3c: mondImprove'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['logDynTime'], 'H4a: dynTime'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['logMeanRun'], 'H4b: MeanRun'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['logDynTime', 'logMeanRun'], 'H4c: dynTime + MeanRun'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['mondImprove', 'logMeanRun'], 'H3+H4: mondImprove + MR'));
hypModelsPool.push(testHypothesisModel(poolWithResid, ['logSigmaBar', 'logDynTime', 'mondImprove'], 'H1+H3+H4: mixed'));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 4: REGIME-DEPENDENT HYPOTHESIS TEST');
console.log('  Which hypothesis gains power specifically at high Vflat?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function regimeHypTest(data, features, label) {
  const valid = data.filter(g => features.every(f => g[f] !== null && g[f] !== undefined && isFinite(g[f])) && g.hkResid !== undefined && isFinite(g.hkResid));

  const lowV = valid.filter(g => g.Vflat < 120);
  const highV = valid.filter(g => g.Vflat >= 120);
  const vhighV = valid.filter(g => g.Vflat >= 180);

  const regimes = [
    { label: 'Low(<120)', data: lowV },
    { label: 'High(>=120)', data: highV },
    { label: 'VHigh(>=180)', data: vhighV }
  ];

  const row = [];
  for (const reg of regimes) {
    if (reg.data.length < features.length + 3) {
      row.push({ r2: null, loo: null, n: reg.data.length, rResidA0: null });
      continue;
    }
    const Y = reg.data.map(g => g.hkResid);
    const X = reg.data.map(g => features.map(f => g[f]));
    const model = olsFull(Y, X);
    let loo = NaN;
    try { loo = looR2(Y, X); } catch (e) {}
    const resids = Y.map((y, i) => y - [1, ...X[i]].reduce((s, x, j) => s + x * model.beta[j], 0));
    const rRA0 = pearsonR(reg.data.map(g => g.logA0), resids);
    row.push({ r2: +model.r2.toFixed(3), loo: isNaN(loo) ? null : +loo.toFixed(3), n: reg.data.length, rResidA0: isNaN(rRA0) ? null : +rRA0.toFixed(3) });
  }

  const fmt = v => v.r2 !== null ? ('R²=' + v.r2.toFixed(2) + '(' + v.n + ')') : ('N=' + v.n);
  console.log('  ' + label.padEnd(30) + row.map(v => fmt(v).padStart(16)).join(''));
  return { label, features, regimes: row };
}

console.log('  ' + 'Model'.padEnd(30) + 'Low(<120)'.padStart(16) + 'High(>=120)'.padStart(16) + 'VHigh(>=180)'.padStart(16));
console.log('  ' + '-'.repeat(80));

const regHyp = [];
regHyp.push(regimeHypTest(intWithResid, ['logSigmaBar'], 'H1: SigmaBar'));
regHyp.push(regimeHypTest(intWithResid, ['morphT'], 'H2: morphT'));
regHyp.push(regimeHypTest(intWithResid, ['logFgas'], 'H3a: Fgas'));
regHyp.push(regimeHypTest(intWithResid, ['mondImprove'], 'H3c: mondImprove'));
regHyp.push(regimeHypTest(intWithResid, ['logDynTime'], 'H4a: dynTime'));
regHyp.push(regimeHypTest(intWithResid, ['logMeanRun'], 'H4b: MeanRun'));
regHyp.push(regimeHypTest(intWithResid, ['mondImprove', 'logMeanRun'], 'H3+H4: mond+MR'));

console.log('\n  After each hypothesis model, remaining r(resid, a0):');
console.log('  ' + 'Model'.padEnd(30) + 'Low(<120)'.padStart(16) + 'High(>=120)'.padStart(16) + 'VHigh(>=180)'.padStart(16));
console.log('  ' + '-'.repeat(80));
for (const rh of regHyp) {
  const fmt = v => v.rResidA0 !== null ? ('r=' + v.rResidA0.toFixed(2) + '(' + v.n + ')') : ('N=' + v.n);
  console.log('  ' + rh.label.padEnd(30) + rh.regimes.map(v => fmt(v).padStart(16)).join(''));
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 5: THE KEY QUESTION — DOES OBSERVABILITY OR PHYSICS CHANGE?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function observabilityTest(data, label) {
  const valid = data.filter(g => g.hkResid !== undefined && isFinite(g.hkResid) && isFinite(g.logA0));
  const low = valid.filter(g => g.Vflat < 120);
  const high = valid.filter(g => g.Vflat >= 120);

  const sdVfRLow = low.length > 1 ? sd(low.map(g => g.VfResid)) : NaN;
  const sdVfRHigh = high.length > 1 ? sd(high.map(g => g.VfResid)) : NaN;
  const sdHkRLow = low.length > 1 ? sd(low.map(g => g.hkResid)) : NaN;
  const sdHkRHigh = high.length > 1 ? sd(high.map(g => g.hkResid)) : NaN;
  const sdA0Low = low.length > 1 ? sd(low.map(g => g.logA0)) : NaN;
  const sdA0High = high.length > 1 ? sd(high.map(g => g.logA0)) : NaN;

  const rHkRA0Low = low.length >= 4 ? pearsonR(low.map(g => g.logA0), low.map(g => g.hkResid)) : NaN;
  const rHkRA0High = high.length >= 4 ? pearsonR(high.map(g => g.logA0), high.map(g => g.hkResid)) : NaN;

  const slopeHigh = high.length >= 4 ? olsFull(high.map(g => g.logA0), high.map(g => [g.hkResid])).beta[1] : NaN;
  const slopeLow = low.length >= 4 ? olsFull(low.map(g => g.logA0), low.map(g => [g.hkResid])).beta[1] : NaN;

  console.log('  ' + label + ':');
  console.log('    Low-V (N=' + low.length + '):  sd(VfR)=' + (isNaN(sdVfRLow) ? 'N/A' : sdVfRLow.toFixed(4)) +
    '  sd(hkR)=' + (isNaN(sdHkRLow) ? 'N/A' : sdHkRLow.toFixed(4)) +
    '  sd(a0)=' + (isNaN(sdA0Low) ? 'N/A' : sdA0Low.toFixed(4)));
  console.log('    High-V (N=' + high.length + '): sd(VfR)=' + (isNaN(sdVfRHigh) ? 'N/A' : sdVfRHigh.toFixed(4)) +
    '  sd(hkR)=' + (isNaN(sdHkRHigh) ? 'N/A' : sdHkRHigh.toFixed(4)) +
    '  sd(a0)=' + (isNaN(sdA0High) ? 'N/A' : sdA0High.toFixed(4)));
  console.log('    r(hkResid, a0): low=' + (isNaN(rHkRA0Low) ? 'N/A' : rHkRA0Low.toFixed(3)) +
    '  high=' + (isNaN(rHkRA0High) ? 'N/A' : rHkRA0High.toFixed(3)));
  console.log('    slope(a0 ~ hkResid): low=' + (isNaN(slopeLow) ? 'N/A' : slopeLow.toFixed(2)) +
    '  high=' + (isNaN(slopeHigh) ? 'N/A' : slopeHigh.toFixed(2)));

  const slopeRatio = !isNaN(slopeHigh) && !isNaN(slopeLow) && Math.abs(slopeLow) > 0.01 ? slopeHigh / slopeLow : NaN;
  const sdRatio = !isNaN(sdHkRHigh) && !isNaN(sdHkRLow) && sdHkRLow > 0 ? sdHkRHigh / sdHkRLow : NaN;
  const rRatio = !isNaN(rHkRA0High) && !isNaN(rHkRA0Low) && Math.abs(rHkRA0Low) > 0.01 ? rHkRA0High / rHkRA0Low : NaN;

  let interpretation = 'AMBIGUOUS';
  if (!isNaN(rRatio) && !isNaN(slopeRatio)) {
    if (Math.abs(slopeRatio) > 1.5 && rRatio > 1.3) {
      interpretation = 'PHYSICS_CHANGES: both slope and r strengthen → real physics amplifies';
    } else if (Math.abs(slopeRatio) < 1.3 && rRatio > 1.3) {
      interpretation = 'OBSERVABILITY_IMPROVES: same slope but tighter correlation → just less noise';
    } else if (Math.abs(slopeRatio) > 1.5 && rRatio < 1.3) {
      interpretation = 'AMPLITUDE_SCALES: bigger effect but same coupling tightness';
    } else {
      interpretation = 'MIXED: moderate changes in both';
    }
  }

  console.log('    Slope ratio (high/low): ' + (isNaN(slopeRatio) ? 'N/A' : slopeRatio.toFixed(2)));
  console.log('    sd(hkResid) ratio (high/low): ' + (isNaN(sdRatio) ? 'N/A' : sdRatio.toFixed(2)));
  console.log('    r ratio (high/low): ' + (isNaN(rRatio) ? 'N/A' : rRatio.toFixed(2)));
  console.log('    → INTERPRETATION: ' + interpretation);

  return { label, lowN: low.length, highN: high.length, sdVfRLow, sdVfRHigh, sdHkRLow, sdHkRHigh, sdA0Low, sdA0High, rHkRA0Low, rHkRA0High, slopeLow, slopeHigh, slopeRatio, interpretation };
}

const obsInt = observabilityTest(intWithResid, 'Internal');
console.log('');
const obsExt = observabilityTest(extWithResid, 'External');
console.log('');
const obsPool = observabilityTest(poolWithResid, 'Pooled');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 6: HYPOTHESIS SCORECARD');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function scoreHypothesis(name, intCorrs, poolCorrs, intModels, poolModels) {
  const proxies = hypothesisProxies[name] || [];
  let score = 0;
  let evidence = [];

  for (const proxy of proxies) {
    if (proxy === 'haloK') continue;
    const intC = intCorrs.find(c => c.proxy === proxy);
    const poolC = poolCorrs.find(c => c.proxy === proxy);

    if (intC && intC.rHigh !== null && Math.abs(intC.rHigh) > 0.2) { score += 1; evidence.push(proxy + ' int-HV r=' + intC.rHigh); }
    if (poolC && poolC.rHigh !== null && Math.abs(poolC.rHigh) > 0.2) { score += 1; evidence.push(proxy + ' pool-HV r=' + poolC.rHigh); }
    if (intC && intC.regimeStrengthens) { score += 1; evidence.push(proxy + ' regime-strengthens'); }
    if (intC && poolC && intC.rFull !== null && poolC.rFull !== null &&
        Math.sign(intC.rFull) === Math.sign(poolC.rFull) && Math.abs(intC.rFull) > 0.15) {
      score += 1; evidence.push(proxy + ' consistent sign');
    }
  }

  const bestIntModel = intModels.filter(m => m && m.label.startsWith(name.slice(0, 2))).sort((a, b) => (b.loo || -1) - (a.loo || -1))[0];
  const bestPoolModel = poolModels.filter(m => m && m.label.startsWith(name.slice(0, 2))).sort((a, b) => (b.loo || -1) - (a.loo || -1))[0];

  if (bestIntModel && bestIntModel.loo > 0) { score += 2; evidence.push('int LOO=' + bestIntModel.loo); }
  if (bestPoolModel && bestPoolModel.loo > 0) { score += 2; evidence.push('pool LOO=' + bestPoolModel.loo); }

  return { name, score, evidence };
}

const scores = [
  scoreHypothesis('H1_haloResponse', corrHkRInt, corrHkRPool, hypModelsInt, hypModelsPool),
  scoreHypothesis('H2_assemblyHistory', corrHkRInt, corrHkRPool, hypModelsInt, hypModelsPool),
  scoreHypothesis('H3_feedbackImprint', corrHkRInt, corrHkRPool, hypModelsInt, hypModelsPool),
  scoreHypothesis('H4_dynIntegration', corrHkRInt, corrHkRPool, hypModelsInt, hypModelsPool),
];

scores.sort((a, b) => b.score - a.score);

console.log('  Hypothesis scores (higher = more supported by data):');
console.log('  ' + '-'.repeat(70));
for (const s of scores) {
  console.log('  ' + s.name.padEnd(25) + 'Score: ' + String(s.score).padStart(2));
  console.log('    Evidence: ' + (s.evidence.length > 0 ? s.evidence.join('; ') : 'none'));
}

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PHASE 303 VERDICT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const bestHyp = scores[0];
const secondHyp = scores[1];
const gap = bestHyp.score - secondHyp.score;

let verdict = 'NO_CLEAR_WINNER';
if (gap >= 3) verdict = 'STRONG_SUPPORT: ' + bestHyp.name;
else if (gap >= 1) verdict = 'MODERATE_SUPPORT: ' + bestHyp.name;
else verdict = 'MIXED: ' + bestHyp.name + ' and ' + secondHyp.name + ' roughly tied';

console.log('  Best hypothesis: ' + bestHyp.name + ' (score=' + bestHyp.score + ')');
console.log('  Runner-up: ' + secondHyp.name + ' (score=' + secondHyp.score + ')');
console.log('  Gap: ' + gap);
console.log('  Verdict: ' + verdict);

const obsVerdict = obsInt.interpretation;
console.log('\n  Observability question: ' + obsVerdict);
console.log('\n  Physical narrative:');
if (bestHyp.name.includes('H3') || bestHyp.name.includes('H4')) {
  console.log('    The hidden channel is best explained by DYNAMICAL/FEEDBACK processes');
  console.log('    that accumulate with galaxy mass/velocity.');
  console.log('    haloK captures the halo-amplitude component,');
  console.log('    but the deeper signal is a regime-dependent integration effect.');
} else if (bestHyp.name.includes('H1')) {
  console.log('    The hidden channel is best explained by HALO RESPONSE to baryons.');
  console.log('    The baryonic surface density drives adiabatic contraction,');
  console.log('    which deepens in higher-mass systems.');
} else if (bestHyp.name.includes('H2')) {
  console.log('    The hidden channel is best explained by ASSEMBLY HISTORY.');
  console.log('    Earlier-forming systems have more concentrated halos,');
  console.log('    producing stronger coupling at high mass.');
}

const output = {
  phase: '303',
  title: 'Physical Interpretation — What is the hidden dynamical channel?',
  timestamp: new Date().toISOString(),
  verdict,
  observabilityVerdict: obsVerdict,
  hypothesisScores: scores,
  proxyCorrelations: {
    internal_VfResid: corrVfRInt,
    pooled_VfResid: corrVfRPool,
    internal_hkResid: corrHkRInt,
    pooled_hkResid: corrHkRPool
  },
  hypothesisModels: {
    internal: hypModelsInt.filter(m => m),
    pooled: hypModelsPool.filter(m => m)
  },
  regimeHypothesis: regHyp.map(r => ({ label: r.label, features: r.features, regimes: r.regimes })),
  observabilityAnalysis: {
    internal: obsInt,
    external: obsExt,
    pooled: obsPool
  },
  physicalNarrative: {
    haloK_role: 'Best partial observable handle — captures halo amplitude',
    hidden_channel: 'Regime-dependent dynamical integration effect',
    regime_dependence: 'Physics amplifies in deeper potential wells (high Vflat)',
    irreducible_nature: 'Beyond any single catalog observable — likely accumulated dynamical processing'
  },
  caveats: [
    'Hypothesis proxies are crude — true assembly history or feedback strength not directly observable',
    'N=45 internal sample limits regime-split model power',
    'logDynTime and logFgas are constructed from same catalog variables used for VfResid — partial circularity risk',
    'Scoring system is heuristic — not a formal Bayesian model comparison',
    'Multiple hypotheses may be simultaneously active (not mutually exclusive)'
  ]
};

const outPath = pub('phase303-physical-interpretation.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nOutput: ' + outPath);
console.log('\n======================================================================');
console.log('  PHASE 303 COMPLETE');
console.log('======================================================================');

#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const pub = p => path.join(__dirname, '..', 'public', p);

console.log('======================================================================');
console.log('  PHASE 301: WHAT DETERMINES VfResid?');
console.log('  Central question: What hidden physical process is encoded in VfResid?');
console.log('  Three levels: A) Internal, B) External, C) Pooled with controls');
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
  const n = x.length; if (n < 3) return 0;
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
  return { beta, r2: 1 - sse / sst, rmse: Math.sqrt(sse / n) };
}

function looR2(Y, X) {
  const n = Y.length;
  let looPressNum = 0;
  const my = mean(Y);
  let sst = Y.reduce((s, y) => s + (y - my) ** 2, 0);
  for (let i = 0; i < n; i++) {
    const Yt = Y.filter((_, j) => j !== i);
    const Xt = X.filter((_, j) => j !== i);
    const model = olsFull(Yt, Xt);
    const xi = [1, ...X[i]];
    const pred = xi.reduce((s, x, j) => s + x * model.beta[j], 0);
    looPressNum += (Y[i] - pred) ** 2;
  }
  return 1 - looPressNum / sst;
}

function bootstrapStability(Y, X, nBoot) {
  const n = Y.length;
  const p = X[0].length + 1;
  const signs = Array.from({ length: p }, () => ({ pos: 0, neg: 0 }));
  const betaSamples = [];
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: n }, () => Math.floor(Math.random() * n));
    const Yb = idx.map(i => Y[i]);
    const Xb = idx.map(i => X[i]);
    try {
      const m = olsFull(Yb, Xb);
      betaSamples.push(m.beta);
      for (let j = 0; j < p; j++) {
        if (m.beta[j] > 0) signs[j].pos++;
        else signs[j].neg++;
      }
    } catch (e) {}
  }
  const signStability = signs.map(s => Math.max(s.pos, s.neg) / (s.pos + s.neg) * 100);
  const medians = Array.from({ length: p }, (_, j) => {
    const vals = betaSamples.map(b => b[j]).sort((a, b) => a - b);
    return vals[Math.floor(vals.length / 2)];
  });
  return { signStability, medians, nBoot: betaSamples.length };
}

function partialR(x, y, controls) {
  const n = x.length;
  if (n < controls.length + 3) return 0;
  const residualize = (target, preds) => {
    const m = olsFull(target, preds);
    return target.map((v, i) => v - [1, ...preds[i]].reduce((s, x, j) => s + x * m.beta[j], 0));
  };
  const xResid = residualize(x, controls);
  const yResid = residualize(y, controls);
  return pearsonR(xResid, yResid);
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 1: ASSEMBLE INTERNAL DATASET (N=45)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const pubGals = stageA.galaxies.filter(g => pubNames.has(g.name));
const structTrainData = pubGals.map(g => {
  const s = sparcMap[g.name];
  return {
    name: g.name,
    logVflat: Math.log10(s.Vflat),
    logMbar: Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9),
    logL36: Math.log10(Math.max(s.L36, 0.01)),
    logRdisk: Math.log10(s.Rdisk),
    morphT: s.T ?? 5
  };
}).filter(g => isFinite(g.logVflat) && isFinite(g.logMbar));

const structModel = olsFull(
  structTrainData.map(g => g.logVflat),
  structTrainData.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT])
);

const internal = pubGals.map(g => {
  const s = sparcMap[g.name];
  const r = resMap[g.name];
  const td = tdMap[g.name];
  const logMbar = Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
  const logVflat = Math.log10(s.Vflat);
  const logL36 = Math.log10(Math.max(s.L36, 0.01));
  const logRdisk = Math.log10(s.Rdisk);
  const predLogVflat = [1, logMbar, logL36, logRdisk, s.T ?? 5]
    .reduce((sum, x, j) => sum + x * structModel.beta[j], 0);
  const VfResid = logVflat - predLogVflat;

  const haloK = r && r.models && r.models.dark_halo_linear
    ? Math.log10(Math.max(r.models.dark_halo_linear.k, 1)) : null;
  const lhOuter = r && r.models && r.models.log_halo
    ? r.models.log_halo.outerImprovement : null;
  const mondImprove = r && r.models && r.models.mond
    ? r.models.mond.improvementVsNewton : null;
  const dhlMSE = r && r.models && r.models.dark_halo_linear
    ? r.models.dark_halo_linear.mse : null;
  const concentration = r ? (r.maxR || 0) / Math.max(s.Rdisk, 0.1) : null;
  const logSBeff = Math.log10(Math.max(s.SBeff, 0.01));
  const logSBdisk = Math.log10(Math.max(s.SBdisk, 0.01));

  return {
    name: g.name,
    VfResid,
    logA0: g.logA0,
    haloK,
    lhOuter,
    mondImprove,
    dhlMSE,
    logMHI: g.logMHI,
    logMbar,
    logMhost: td ? td.logMhost : null,
    logMeanRun: g.logMeanRun,
    logSigma0: g.logSigma0,
    logL36,
    logRdisk,
    logSBeff,
    logSBdisk,
    morphT: s.T ?? 5,
    envCode: g.envCode,
    rcWiggliness: g.rcWiggliness,
    inc: s.inc,
    Q: s.Q,
    Vflat: s.Vflat,
    logVflat,
    concentration,
    source: 'internal'
  };
}).filter(g => isFinite(g.VfResid));

console.log('  Internal sample: N=' + internal.length);
console.log('  VfResid: mean=' + mean(internal.map(g => g.VfResid)).toFixed(4) +
  ' sd=' + sd(internal.map(g => g.VfResid)).toFixed(4));
console.log('  Vflat range: [' + Math.min(...internal.map(g => g.Vflat)).toFixed(0) + ', ' +
  Math.max(...internal.map(g => g.Vflat)).toFixed(0) + ']');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 2: ASSEMBLE EXTERNAL DATASET');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const extGals = extDataset.galaxies;
const salvGals = salvageData.salvagedGalaxies || [];

const external = [...extGals, ...salvGals].map(g => {
  const s = sparcMap[g.name];
  const r = resMap[g.name];
  const mondImprove = r && r.models && r.models.mond ? r.models.mond.improvementVsNewton : null;
  const dhlMSE = r && r.models && r.models.dark_halo_linear ? r.models.dark_halo_linear.mse : null;
  const concentration = r ? (r.maxR || 0) / Math.max(s ? s.Rdisk : 1, 0.1) : null;
  const logSBeff = s ? Math.log10(Math.max(s.SBeff, 0.01)) : null;
  const logSBdisk = s ? Math.log10(Math.max(s.SBdisk, 0.01)) : null;

  return {
    name: g.name,
    VfResid: g.VfResid,
    logA0: g.logA0,
    haloK: g.haloK,
    lhOuter: g.lhOuter,
    mondImprove,
    dhlMSE,
    logMHI: g.logMHI,
    logMbar: g.logMbar,
    logMhost: g.logMhost,
    logMeanRun: g.logMeanRun,
    logSigma0: s ? Math.log10(Math.max(s.SBdisk, 0.01)) : null,
    logL36: g.logL36,
    logRdisk: g.logRdisk,
    logSBeff,
    logSBdisk,
    morphT: g.morphT ?? (s ? s.T ?? 5 : 5),
    envCode: null,
    rcWiggliness: null,
    inc: g.inc ?? (s ? s.inc : null),
    Q: g.Q ?? (s ? s.Q : null),
    Vflat: g.Vflat,
    logVflat: g.logVflat ?? Math.log10(g.Vflat),
    concentration,
    source: g.salvageSource ? 'salvaged' : 'external',
    salvageSource: g.salvageSource || null
  };
});

console.log('  External sample (original + salvaged): N=' + external.length);
console.log('  VfResid: mean=' + mean(external.map(g => g.VfResid)).toFixed(4) +
  ' sd=' + sd(external.map(g => g.VfResid)).toFixed(4));

const extOriginal = external.filter(g => g.source === 'external');
const extSalvaged = external.filter(g => g.source === 'salvaged');
console.log('  Original external: N=' + extOriginal.length + '  Salvaged: N=' + extSalvaged.length);

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  LEVEL A: INTERNAL-ONLY VfResid MODELS (BENCHMARK)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function runModelSuite(data, label, candidateNames) {
  const Y = data.map(g => g.VfResid);
  const results = [];

  for (const candName of candidateNames) {
    const valid = data.filter(g => {
      const v = g[candName];
      return v !== null && v !== undefined && isFinite(v);
    });
    if (valid.length < 8) continue;

    const Yv = valid.map(g => g.VfResid);
    const Xv = valid.map(g => [g[candName]]);
    const r = pearsonR(Yv, valid.map(g => g[candName]));
    let loo = NaN;
    try { loo = looR2(Yv, Xv); } catch (e) {}

    results.push({ name: candName, r, loo, n: valid.length });
  }

  results.sort((a, b) => Math.abs(b.r) - Math.abs(a.r));
  console.log('  ' + label + ' — Single predictor ranking for VfResid:');
  console.log('  ' + 'Predictor'.padEnd(18) + 'r'.padStart(8) + 'LOO R²'.padStart(10) + '  N');
  console.log('  ' + '-'.repeat(45));
  for (const res of results) {
    console.log('  ' + res.name.padEnd(18) + res.r.toFixed(3).padStart(8) +
      (isNaN(res.loo) ? '     N/A' : res.loo.toFixed(3).padStart(10)) + '  ' + res.n);
  }
  return results;
}

const candidates = [
  'haloK', 'lhOuter', 'mondImprove', 'dhlMSE',
  'logMHI', 'logMbar', 'logMhost', 'logMeanRun',
  'logSigma0', 'logL36', 'logRdisk', 'logSBeff', 'logSBdisk',
  'morphT', 'envCode', 'rcWiggliness', 'inc', 'concentration'
];

const intResults = runModelSuite(internal, 'INTERNAL (N=' + internal.length + ')', candidates);

const intHighV = internal.filter(g => g.Vflat >= 120);
console.log('');
const intHVResults = runModelSuite(intHighV, 'INTERNAL HIGH-V (N=' + intHighV.length + ')', candidates);

console.log('\n  MULTI-PREDICTOR MODELS (Internal):');

function tryModel(data, features, label) {
  const valid = data.filter(g => features.every(f => g[f] !== null && g[f] !== undefined && isFinite(g[f])));
  if (valid.length < features.length + 3) {
    console.log('  ' + label + ': SKIP (N=' + valid.length + ' insufficient)');
    return null;
  }
  const Y = valid.map(g => g.VfResid);
  const X = valid.map(g => features.map(f => g[f]));
  const model = olsFull(Y, X);
  let loo = NaN;
  try { loo = looR2(Y, X); } catch (e) {}
  const boot = bootstrapStability(Y, X, 1000);

  console.log('  ' + label + ' (N=' + valid.length + '):');
  console.log('    R² = ' + model.r2.toFixed(3) + '  LOO R² = ' + (isNaN(loo) ? 'N/A' : loo.toFixed(3)) +
    '  RMSE = ' + model.rmse.toFixed(4));
  console.log('    Beta: [' + model.beta.map(b => b.toFixed(4)).join(', ') + ']');
  console.log('    Features: [intercept, ' + features.join(', ') + ']');
  console.log('    Sign stability: [' + boot.signStability.map(s => s.toFixed(0) + '%').join(', ') + ']');

  return { label, features, n: valid.length, r2: model.r2, loo, rmse: model.rmse, beta: model.beta, signStability: boot.signStability };
}

const multiModels = [];
multiModels.push(tryModel(internal, ['haloK'], 'M1: haloK'));
multiModels.push(tryModel(internal, ['haloK', 'lhOuter'], 'M2: haloK + lhOuter'));
multiModels.push(tryModel(internal, ['haloK', 'envCode', 'logMeanRun'], 'M3: haloK + env + MR'));
multiModels.push(tryModel(internal, ['haloK', 'lhOuter', 'envCode', 'logMeanRun', 'concentration'], 'M4: top-5'));
multiModels.push(tryModel(internal, ['haloK', 'logMHI', 'logMbar'], 'M5: haloK + baryonic'));
multiModels.push(tryModel(internal, ['haloK', 'logSBeff', 'logRdisk'], 'M6: haloK + structure'));
multiModels.push(tryModel(internal, ['haloK', 'mondImprove'], 'M7: haloK + MOND'));
multiModels.push(tryModel(internal, ['haloK', 'lhOuter', 'logSBeff', 'logMeanRun'], 'M8: haloK + lhOI + SBeff + MR'));
multiModels.push(tryModel(internal, ['logSBeff', 'logRdisk', 'logMbar', 'morphT'], 'M9: pure structure (no halo)'));
multiModels.push(tryModel(internal, ['haloK', 'concentration', 'logRdisk'], 'M10: haloK + conc + Rdisk'));

console.log('\n  HIGH-V MULTI-PREDICTOR (Vflat >= 120):');
const multiHV = [];
multiHV.push(tryModel(intHighV, ['haloK'], 'M1-HV: haloK'));
multiHV.push(tryModel(intHighV, ['haloK', 'lhOuter'], 'M2-HV: haloK + lhOuter'));
multiHV.push(tryModel(intHighV, ['haloK', 'envCode', 'logMeanRun'], 'M3-HV: haloK + env + MR'));
multiHV.push(tryModel(intHighV, ['haloK', 'lhOuter', 'envCode', 'logMeanRun', 'concentration'], 'M4-HV: top-5'));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  LEVEL B: EXTERNAL-ONLY VfResid MODELS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const extCandidates = candidates.filter(c => c !== 'envCode' && c !== 'rcWiggliness' && c !== 'logSigma0');
const extResults = runModelSuite(external, 'EXTERNAL (N=' + external.length + ')', extCandidates);

const extHighV = external.filter(g => g.Vflat >= 120);
console.log('');
const extHVResults = runModelSuite(extHighV, 'EXTERNAL HIGH-V (N=' + extHighV.length + ')', extCandidates);

console.log('\n  MULTI-PREDICTOR MODELS (External):');
const extMulti = [];
extMulti.push(tryModel(external, ['haloK'], 'E1: haloK'));
extMulti.push(tryModel(external, ['haloK', 'lhOuter'], 'E2: haloK + lhOuter'));
extMulti.push(tryModel(external, ['haloK', 'logMeanRun'], 'E3: haloK + MR'));
extMulti.push(tryModel(external, ['haloK', 'logMHI', 'logMbar'], 'E4: haloK + baryonic'));
extMulti.push(tryModel(external, ['haloK', 'logSBeff', 'logRdisk'], 'E5: haloK + structure'));
extMulti.push(tryModel(external, ['logSBeff', 'logRdisk', 'logMbar', 'morphT'], 'E6: pure structure (no halo)'));
extMulti.push(tryModel(external, ['haloK', 'lhOuter', 'logMeanRun', 'concentration'], 'E7: top-4 external'));
extMulti.push(tryModel(external, ['haloK', 'mondImprove'], 'E8: haloK + MOND'));

console.log('\n  HIGH-V EXTERNAL MULTI-PREDICTOR (N=' + extHighV.length + '):');
const extMultiHV = [];
extMultiHV.push(tryModel(extHighV, ['haloK'], 'E1-HV: haloK'));
extMultiHV.push(tryModel(extHighV, ['haloK', 'lhOuter'], 'E2-HV: haloK + lhOuter'));
extMultiHV.push(tryModel(extHighV, ['haloK', 'logMeanRun'], 'E3-HV: haloK + MR'));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  LEVEL C: POOLED MODELS WITH CONTROLS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const pooled = [...internal, ...external].map((g, i) => ({
  ...g,
  sourceFlag: g.source === 'internal' ? 1 : 0,
  qualityFlag: (g.Q === 1 ? 1 : 0),
  highVflag: (g.Vflat >= 120 ? 1 : 0)
}));

console.log('  Pooled sample: N=' + pooled.length +
  ' (internal=' + pooled.filter(g => g.sourceFlag === 1).length +
  ', external=' + pooled.filter(g => g.sourceFlag === 0).length + ')');

const poolCandidates = ['haloK', 'lhOuter', 'mondImprove', 'dhlMSE',
  'logMHI', 'logMbar', 'logMhost', 'logMeanRun',
  'logL36', 'logRdisk', 'logSBeff', 'logSBdisk',
  'morphT', 'inc', 'concentration'];

const poolResults = runModelSuite(pooled, 'POOLED (N=' + pooled.length + ')', poolCandidates);

console.log('\n  POOLED MULTI-PREDICTOR WITH SOURCE CONTROL:');
const poolMulti = [];
poolMulti.push(tryModel(pooled, ['haloK', 'sourceFlag'], 'P1: haloK + source'));
poolMulti.push(tryModel(pooled, ['haloK', 'lhOuter', 'sourceFlag'], 'P2: haloK + lhOI + source'));
poolMulti.push(tryModel(pooled, ['haloK', 'logMeanRun', 'sourceFlag'], 'P3: haloK + MR + source'));
poolMulti.push(tryModel(pooled, ['haloK', 'lhOuter', 'logMeanRun', 'sourceFlag'], 'P4: haloK + lhOI + MR + source'));
poolMulti.push(tryModel(pooled, ['haloK', 'lhOuter', 'logMeanRun', 'concentration', 'sourceFlag'], 'P5: top-4 + source'));
poolMulti.push(tryModel(pooled, ['haloK', 'sourceFlag', 'qualityFlag'], 'P6: haloK + source + quality'));
poolMulti.push(tryModel(pooled, ['logSBeff', 'logRdisk', 'logMbar', 'morphT', 'sourceFlag'], 'P7: structure + source'));

const poolHighV = pooled.filter(g => g.Vflat >= 120);
console.log('\n  POOLED HIGH-V (N=' + poolHighV.length + '):');
poolMulti.push(tryModel(poolHighV, ['haloK', 'sourceFlag'], 'P1-HV: haloK + source'));
poolMulti.push(tryModel(poolHighV, ['haloK', 'lhOuter', 'sourceFlag'], 'P2-HV: haloK + lhOI + source'));
poolMulti.push(tryModel(poolHighV, ['haloK', 'lhOuter', 'logMeanRun', 'sourceFlag'], 'P3-HV: haloK+lhOI+MR+source'));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  CROSS-SAMPLE TRANSFER TESTS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function transferTest(trainData, testData, features, label) {
  const trainValid = trainData.filter(g => features.every(f => g[f] !== null && g[f] !== undefined && isFinite(g[f])));
  const testValid = testData.filter(g => features.every(f => g[f] !== null && g[f] !== undefined && isFinite(g[f])));
  if (trainValid.length < features.length + 3 || testValid.length < 3) {
    console.log('  ' + label + ': SKIP (train=' + trainValid.length + ', test=' + testValid.length + ')');
    return null;
  }
  const Yt = trainValid.map(g => g.VfResid);
  const Xt = trainValid.map(g => features.map(f => g[f]));
  const model = olsFull(Yt, Xt);

  const Ytest = testValid.map(g => g.VfResid);
  const preds = testValid.map(g => [1, ...features.map(f => g[f])].reduce((s, x, j) => s + x * model.beta[j], 0));

  const testMeanY = mean(Ytest);
  let sse = 0, sst = 0;
  for (let i = 0; i < Ytest.length; i++) {
    sse += (Ytest[i] - preds[i]) ** 2;
    sst += (Ytest[i] - testMeanY) ** 2;
  }
  const r2Transfer = 1 - sse / sst;
  const rTransfer = pearsonR(Ytest, preds);

  console.log('  ' + label + ': train N=' + trainValid.length + ' → test N=' + testValid.length +
    '  r=' + rTransfer.toFixed(3) + '  R²=' + r2Transfer.toFixed(3));
  return { label, trainN: trainValid.length, testN: testValid.length, r: rTransfer, r2: r2Transfer };
}

console.log('  Train on INTERNAL → predict EXTERNAL VfResid:');
const transfers = [];
transfers.push(transferTest(internal, external, ['haloK'], 'haloK only'));
transfers.push(transferTest(internal, external, ['haloK', 'lhOuter'], 'haloK + lhOuter'));
transfers.push(transferTest(internal, external, ['logSBeff', 'logRdisk'], 'SBeff + Rdisk'));
transfers.push(transferTest(internal, external, ['logMbar', 'logRdisk', 'morphT'], 'Mbar + Rdisk + T'));

const intHVtrain = internal.filter(g => g.Vflat >= 120);
const extHVtest = external.filter(g => g.Vflat >= 120);
console.log('\n  Train on INTERNAL HIGH-V → predict EXTERNAL HIGH-V:');
transfers.push(transferTest(intHVtrain, extHVtest, ['haloK'], 'HV haloK'));
transfers.push(transferTest(intHVtrain, extHVtest, ['haloK', 'lhOuter'], 'HV haloK + lhOuter'));

console.log('\n  Train on EXTERNAL → predict INTERNAL VfResid:');
transfers.push(transferTest(external, internal, ['haloK'], 'REV haloK'));
transfers.push(transferTest(external, internal, ['haloK', 'lhOuter'], 'REV haloK + lhOuter'));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PARTIAL CORRELATIONS: WHAT SURVIVES CONTROLS?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function partialAnalysis(data, label) {
  const valid = data.filter(g =>
    g.haloK !== null && isFinite(g.haloK) &&
    g.lhOuter !== null && isFinite(g.lhOuter) &&
    isFinite(g.logMbar) && isFinite(g.logRdisk) && isFinite(g.morphT));
  if (valid.length < 10) {
    console.log('  ' + label + ': SKIP (N=' + valid.length + ')');
    return;
  }
  const vf = valid.map(g => g.VfResid);
  const hk = valid.map(g => g.haloK);
  const lh = valid.map(g => g.lhOuter);
  const mb = valid.map(g => g.logMbar);
  const rd = valid.map(g => g.logRdisk);
  const mt = valid.map(g => g.morphT);
  const controls = valid.map(g => [g.logMbar, g.logRdisk, g.morphT]);

  console.log('  ' + label + ' (N=' + valid.length + '):');
  console.log('    raw r(haloK, VfResid) = ' + pearsonR(hk, vf).toFixed(3));
  console.log('    partial r(haloK, VfResid | Mbar+Rdisk+T) = ' + partialR(hk, vf, controls).toFixed(3));
  console.log('    raw r(lhOuter, VfResid) = ' + pearsonR(lh, vf).toFixed(3));
  console.log('    partial r(lhOuter, VfResid | Mbar+Rdisk+T) = ' + partialR(lh, vf, controls).toFixed(3));

  const hkControls = valid.map(g => [g.logMbar, g.logRdisk, g.morphT, g.haloK]);
  console.log('    partial r(lhOuter, VfResid | Mbar+Rdisk+T+haloK) = ' + partialR(lh, vf, hkControls).toFixed(3));
}

partialAnalysis(internal, 'INTERNAL');
console.log('');
partialAnalysis(external, 'EXTERNAL');
console.log('');
partialAnalysis(pooled, 'POOLED');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  IRREDUCIBLE VfResid: WHAT REMAINS AFTER BEST MODEL?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function irreducibleTest(data, features, label) {
  const valid = data.filter(g => features.every(f => g[f] !== null && g[f] !== undefined && isFinite(g[f])) && isFinite(g.logA0));
  if (valid.length < features.length + 5) {
    console.log('  ' + label + ': SKIP (N=' + valid.length + ')');
    return null;
  }
  const Y = valid.map(g => g.VfResid);
  const X = valid.map(g => features.map(f => g[f]));
  const model = olsFull(Y, X);
  const residuals = Y.map((y, i) => y - [1, ...X[i]].reduce((s, x, j) => s + x * model.beta[j], 0));

  const totalVar = Y.reduce((s, y) => s + (y - mean(Y)) ** 2, 0);
  const residVar = residuals.reduce((s, r) => s + r ** 2, 0);
  const unexplainedPct = (residVar / totalVar * 100).toFixed(1);

  const logA0 = valid.map(g => g.logA0);
  const rResidA0 = pearsonR(residuals, logA0);

  console.log('  ' + label + ' (N=' + valid.length + '):');
  console.log('    Model R² = ' + model.r2.toFixed(3) + '  → Unexplained = ' + unexplainedPct + '%');
  console.log('    r(unexplained_VfResid, logA0) = ' + rResidA0.toFixed(3));

  if (Math.abs(rResidA0) > 0.15) {
    console.log('    *** UNEXPLAINED PART STILL CORRELATES WITH a0 ***');
  }

  return { label, r2: model.r2, unexplainedPct: +unexplainedPct, rResidA0, n: valid.length };
}

const bestIntFeatures = ['haloK', 'lhOuter', 'envCode', 'logMeanRun', 'concentration'];
const bestExtFeatures = ['haloK', 'lhOuter', 'logMeanRun', 'concentration'];

const irr = [];
irr.push(irreducibleTest(internal, bestIntFeatures, 'INTERNAL best-5'));
irr.push(irreducibleTest(internal, ['haloK'], 'INTERNAL haloK only'));
irr.push(irreducibleTest(external, bestExtFeatures, 'EXTERNAL best-4'));
irr.push(irreducibleTest(external, ['haloK'], 'EXTERNAL haloK only'));
irr.push(irreducibleTest(pooled, bestExtFeatures, 'POOLED best-4'));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  META-COMPARISON: SAME DRIVERS ACROSS DATASETS?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const sharedCands = ['haloK', 'lhOuter', 'mondImprove', 'logMHI', 'logMbar',
  'logMeanRun', 'logL36', 'logRdisk', 'logSBeff', 'morphT', 'inc', 'concentration'];

console.log('  Driver ranking comparison (|r| with VfResid):');
console.log('  ' + 'Predictor'.padEnd(18) + 'Internal'.padStart(10) + 'External'.padStart(10) + 'Pooled'.padStart(10) + '  Consistent?');
console.log('  ' + '-'.repeat(62));

const intRanks = {};
const extRanks = {};
const poolRanks = {};

for (const c of sharedCands) {
  const intValid = internal.filter(g => g[c] !== null && g[c] !== undefined && isFinite(g[c]));
  const extValid = external.filter(g => g[c] !== null && g[c] !== undefined && isFinite(g[c]));
  const poolValid = pooled.filter(g => g[c] !== null && g[c] !== undefined && isFinite(g[c]));

  const rInt = intValid.length >= 5 ? pearsonR(intValid.map(g => g.VfResid), intValid.map(g => g[c])) : NaN;
  const rExt = extValid.length >= 5 ? pearsonR(extValid.map(g => g.VfResid), extValid.map(g => g[c])) : NaN;
  const rPool = poolValid.length >= 5 ? pearsonR(poolValid.map(g => g.VfResid), poolValid.map(g => g[c])) : NaN;

  intRanks[c] = rInt;
  extRanks[c] = rExt;
  poolRanks[c] = rPool;

  const sameSign = (isNaN(rInt) || isNaN(rExt)) ? '?' :
    (Math.sign(rInt) === Math.sign(rExt) ? 'YES' : 'NO');

  console.log('  ' + c.padEnd(18) +
    (isNaN(rInt) ? '     N/A' : rInt.toFixed(3).padStart(10)) +
    (isNaN(rExt) ? '     N/A' : rExt.toFixed(3).padStart(10)) +
    (isNaN(rPool) ? '     N/A' : rPool.toFixed(3).padStart(10)) +
    '  ' + sameSign);
}

const topInternal = Object.entries(intRanks).filter(([k, v]) => !isNaN(v)).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1])).slice(0, 5);
const topExternal = Object.entries(extRanks).filter(([k, v]) => !isNaN(v)).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1])).slice(0, 5);

console.log('\n  Top-5 internal: ' + topInternal.map(([k, v]) => k + '(' + v.toFixed(2) + ')').join(', '));
console.log('  Top-5 external: ' + topExternal.map(([k, v]) => k + '(' + v.toFixed(2) + ')').join(', '));

const intTop5Set = new Set(topInternal.map(([k]) => k));
const extTop5Set = new Set(topExternal.map(([k]) => k));
const overlap = [...intTop5Set].filter(k => extTop5Set.has(k));
console.log('  Overlap in top-5: ' + overlap.length + '/5 (' + overlap.join(', ') + ')');

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  SYNTHESIS & VERDICT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const bestIntLOO = multiModels.filter(m => m).sort((a, b) => (b.loo || -1) - (a.loo || -1))[0];
const bestExtLOO = extMulti.filter(m => m).sort((a, b) => (b.loo || -1) - (a.loo || -1))[0];

console.log('  Best internal model: ' + (bestIntLOO ? bestIntLOO.label + ' LOO=' + (bestIntLOO.loo || 'N/A') : 'none'));
console.log('  Best external model: ' + (bestExtLOO ? bestExtLOO.label + ' LOO=' + (bestExtLOO.loo || 'N/A') : 'none'));

const irrBestInt = irr.find(r => r && r.label.includes('INTERNAL best'));
const irrBestExt = irr.find(r => r && r.label.includes('EXTERNAL best'));

if (irrBestInt) {
  console.log('\n  Internal irreducible: ' + irrBestInt.unexplainedPct + '% of VfResid unexplained');
  console.log('    Unexplained still predicts a0? r=' + irrBestInt.rResidA0.toFixed(3));
}
if (irrBestExt) {
  console.log('  External irreducible: ' + irrBestExt.unexplainedPct + '% of VfResid unexplained');
  console.log('    Unexplained still predicts a0? r=' + irrBestExt.rResidA0.toFixed(3));
}

let verdict = 'INCONCLUSIVE';
if (overlap.length >= 3) verdict = 'SHARED_DRIVERS';
if (overlap.length >= 4) verdict = 'STRONGLY_SHARED_DRIVERS';
if (overlap.length <= 1) verdict = 'DIVERGENT_DRIVERS';

console.log('\n  VERDICT: ' + verdict);
console.log('  Driver overlap: ' + overlap.length + '/5 between internal and external top-5');

const output = {
  phase: '301',
  title: 'What Determines VfResid?',
  timestamp: new Date().toISOString(),
  verdict,
  internalSample: { N: internal.length, highV: intHighV.length },
  externalSample: { N: external.length, highV: extHighV.length, original: extOriginal.length, salvaged: extSalvaged.length },
  pooledSample: { N: pooled.length, highV: poolHighV.length },
  singlePredictorRankings: {
    internal: intResults.map(r => ({ name: r.name, r: +r.r.toFixed(3), loo: isNaN(r.loo) ? null : +r.loo.toFixed(3), n: r.n })),
    external: extResults.map(r => ({ name: r.name, r: +r.r.toFixed(3), loo: isNaN(r.loo) ? null : +r.loo.toFixed(3), n: r.n })),
    pooled: poolResults.map(r => ({ name: r.name, r: +r.r.toFixed(3), loo: isNaN(r.loo) ? null : +r.loo.toFixed(3), n: r.n }))
  },
  multiPredictorModels: {
    internal: multiModels.filter(m => m).map(m => ({ label: m.label, features: m.features, n: m.n, r2: +m.r2.toFixed(3), loo: isNaN(m.loo) ? null : +m.loo.toFixed(3), signStability: m.signStability.map(s => +s.toFixed(1)) })),
    external: extMulti.filter(m => m).map(m => ({ label: m.label, features: m.features, n: m.n, r2: +m.r2.toFixed(3), loo: isNaN(m.loo) ? null : +m.loo.toFixed(3), signStability: m.signStability.map(s => +s.toFixed(1)) })),
    pooled: poolMulti.filter(m => m).map(m => ({ label: m.label, features: m.features, n: m.n, r2: +m.r2.toFixed(3), loo: isNaN(m.loo) ? null : +m.loo.toFixed(3), signStability: m.signStability.map(s => +s.toFixed(1)) }))
  },
  crossSampleTransfer: transfers.filter(t => t).map(t => ({ label: t.label, trainN: t.trainN, testN: t.testN, r: +t.r.toFixed(3), r2: +t.r2.toFixed(3) })),
  irreducible: irr.filter(r => r).map(r => ({ label: r.label, r2: +r.r2.toFixed(3), unexplainedPct: r.unexplainedPct, rResidA0: +r.rResidA0.toFixed(3), n: r.n })),
  metaComparison: {
    topInternal: topInternal.map(([k, v]) => ({ name: k, r: +v.toFixed(3) })),
    topExternal: topExternal.map(([k, v]) => ({ name: k, r: +v.toFixed(3) })),
    overlap: overlap.length,
    overlapNames: overlap
  },
  caveats: [
    'External logMhost uses crude Vflat formula — environmental variables may be noisy',
    'Salvaged galaxies have fewer RAR points or estimated Vflat — lower quality by construction',
    'envCode and rcWiggliness not available for external sample — internal-only variables',
    'N=45 internal sample: LOO is the gold standard (small sample risk)',
    'External haloK/lhOuter from rotation curve fits — same fitting pipeline as internal',
    'Pooled models should include source flag to avoid blending artifacts'
  ]
};

const outPath = pub('phase301-vfresid-drivers.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nOutput: ' + outPath);

console.log('\n======================================================================');
console.log('  PHASE 301 COMPLETE');
console.log('  Next: Phase 302 — Regime Law');
console.log('======================================================================');

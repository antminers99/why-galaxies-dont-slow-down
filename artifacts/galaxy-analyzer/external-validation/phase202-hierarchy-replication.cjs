#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const pub = p => path.join(__dirname, '..', 'public', p);

console.log('======================================================================');
console.log('  PHASE 202: EXTERNAL HIERARCHY REPLICATION');
console.log('  Does the full hierarchical structure reproduce outside N=45?');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(pub('stage-A-master-table.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(pub('phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const sparcTable = JSON.parse(fs.readFileSync(pub('sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(pub('sparc-results.json'), 'utf8'));
const extDataset = JSON.parse(fs.readFileSync(pub('phase200-external-dataset.json'), 'utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparcTable.forEach(s => { sparcMap[s.name] = s; });
const resMap = {}; sparcResults.perGalaxy.forEach(g => { if (g.rawName) resMap[g.rawName] = g; });

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : 0; }
function sd(a) { const m = mean(a); return a.length > 1 ? Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)) : 0; }
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
  const rmse = Math.sqrt(resid.reduce((s, e) => s + e * e, 0) / n);
  return { beta, resid, rmse };
}

function transferGap(rmse, baselineRmse) {
  return baselineRmse > 0 ? 100 * (1 - rmse ** 2 / baselineRmse ** 2) : 0;
}

const pubGals = stageA.galaxies.filter(g => pubNames.has(g.name));

const structModel = ols(
  pubGals.map(g => Math.log10(sparcMap[g.name].Vflat)),
  pubGals.map(g => {
    const s = sparcMap[g.name];
    return [
      Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9),
      Math.log10(Math.max(s.L36, 0.01)),
      Math.log10(s.Rdisk),
      s.T
    ];
  })
);

function getTrainGal(g) {
  const s = sparcMap[g.name]; const r = resMap[g.name]; const td = tdMap[g.name];
  const logMbar = Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
  const logVflat = Math.log10(s.Vflat);
  const predLV = [1, logMbar, Math.log10(Math.max(s.L36, 0.01)), Math.log10(s.Rdisk), s.T]
    .reduce((sum, x, j) => sum + x * structModel.beta[j], 0);
  return {
    name: g.name, logA0: g.logA0,
    logMHI: g.logMHI, logMhost: td.logMhost, logMR: g.logMeanRun,
    VfResid: logVflat - predLV,
    lhOuter: r.models.log_halo.outerImprovement,
    Vflat: s.Vflat
  };
}

const train = pubGals.map(getTrainGal);
const Y_tr = train.map(g => g.logA0);
const meanY_tr = mean(Y_tr);

const coreModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR]));
const vfResidModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR, g.VfResid]));
const fiveAxisModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR, g.VfResid, g.lhOuter]));

const ext = extDataset.galaxies;
const ext5 = ext.filter(g => g.lhOuter !== null);

function predictModel(model, features) {
  return [1, ...features].reduce((s, x, j) => s + x * model.beta[j], 0);
}

function runTransfer(name, model, featureFn, galaxies, baseline) {
  const preds = galaxies.map(g => predictModel(model, featureFn(g)));
  const errs = galaxies.map((g, i) => g.logA0 - preds[i]);
  const rmse = Math.sqrt(errs.reduce((s, e) => s + e * e, 0) / errs.length);
  const gap = transferGap(rmse, baseline);
  const r = pearsonR(galaxies.map(g => g.logA0), preds);
  return { name, rmse, gap, r, preds, errs };
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: NESTED HIERARCHY — EACH TIER TESTED');
console.log('  Internal reference: Core 44.1% → +VfResid 61.1% → +lhOuter 65.4%');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const vfOnlyModel = ols(Y_tr, train.map(g => [g.VfResid]));
const lhOnlyModel = ols(Y_tr, train.map(g => [g.lhOuter]));
const vfLhModel = ols(Y_tr, train.map(g => [g.VfResid, g.lhOuter]));
const coreLhModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR, g.lhOuter]));

function testRegime(label, gals, gals5) {
  if (gals.length < 4) { console.log('  ' + label + ': insufficient N=' + gals.length); return null; }
  const naive = Math.sqrt(gals.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / gals.length);
  const naive5 = gals5.length >= 4 ? Math.sqrt(gals5.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / gals5.length) : naive;

  const core = runTransfer('Core', coreModel, g => [g.logMHI, g.logMhost, g.logMeanRun], gals, naive);
  const vfr = runTransfer('Core+VfResid', vfResidModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid], gals, naive);
  const five = gals5.length >= 4 ? runTransfer('5-axis', fiveAxisModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid, g.lhOuter], gals5, naive5) : null;

  const vfOnly = runTransfer('VfResid-only', vfOnlyModel, g => [g.VfResid], gals, naive);
  const coreLh = runTransfer('Core+lhOuter', coreLhModel, g => [g.logMHI, g.logMhost, g.logMeanRun, g.lhOuter], gals5, naive5);

  console.log('  ' + label + ' (N=' + gals.length + ')');
  console.log('  Model'.padEnd(25) + 'RMSE'.padStart(8) + 'Gap%'.padStart(8) + 'Delta'.padStart(8) + 'r'.padStart(8));
  console.log('  ' + '─'.repeat(55));
  console.log('  ' + 'Naive'.padEnd(25) + naive.toFixed(4).padStart(8));
  console.log('  ' + 'VfResid-only'.padEnd(25) + vfOnly.rmse.toFixed(4).padStart(8) + (vfOnly.gap.toFixed(1) + '%').padStart(8) + ''.padStart(8) + vfOnly.r.toFixed(3).padStart(8));
  console.log('  ' + 'Core (3-axis)'.padEnd(25) + core.rmse.toFixed(4).padStart(8) + (core.gap.toFixed(1) + '%').padStart(8) + ''.padStart(8) + core.r.toFixed(3).padStart(8));
  console.log('  ' + 'Core+lhOuter'.padEnd(25) + coreLh.rmse.toFixed(4).padStart(8) + (coreLh.gap.toFixed(1) + '%').padStart(8) + ((coreLh.gap - core.gap).toFixed(1) + 'pp').padStart(8) + coreLh.r.toFixed(3).padStart(8));
  console.log('  ' + 'Core+VfResid'.padEnd(25) + vfr.rmse.toFixed(4).padStart(8) + (vfr.gap.toFixed(1) + '%').padStart(8) + ((vfr.gap - core.gap).toFixed(1) + 'pp').padStart(8) + vfr.r.toFixed(3).padStart(8));
  if (five) console.log('  ' + '5-axis (full)'.padEnd(25) + five.rmse.toFixed(4).padStart(8) + (five.gap.toFixed(1) + '%').padStart(8) + ((five.gap - vfr.gap).toFixed(1) + 'pp').padStart(8) + five.r.toFixed(3).padStart(8));

  const vfrDominates = vfr.gap > core.gap && vfr.gap > coreLh.gap;
  const coreFails = core.gap < 5;
  const lhAdds = five ? five.gap > vfr.gap : false;

  console.log('');
  console.log('  Hierarchy check:');
  console.log('    Core alone fails (<5% gap):      ' + (coreFails ? 'YES' : 'NO') + ' (gap=' + core.gap.toFixed(1) + '%)');
  console.log('    VfResid dominates over Core+lh:   ' + (vfrDominates ? 'YES' : 'NO'));
  console.log('    VfResid-only beats Core:          ' + (vfOnly.gap > core.gap ? 'YES' : 'NO'));
  if (five) console.log('    lhOuter adds above VfResid:      ' + (lhAdds ? 'YES (+' + (five.gap - vfr.gap).toFixed(1) + 'pp)' : 'NO'));
  console.log('');

  return {
    label, N: gals.length, naive, coreFails, vfrDominates, lhAdds,
    core: { gap: +core.gap.toFixed(1), r: +core.r.toFixed(3) },
    vfOnly: { gap: +vfOnly.gap.toFixed(1), r: +vfOnly.r.toFixed(3) },
    coreLh: { gap: +coreLh.gap.toFixed(1), r: +coreLh.r.toFixed(3) },
    vfResid: { gap: +vfr.gap.toFixed(1), r: +vfr.r.toFixed(3) },
    fiveAxis: five ? { gap: +five.gap.toFixed(1), r: +five.r.toFixed(3) } : null
  };
}

const fullResult = testRegime('FULL SAMPLE', ext, ext5);
const highV = ext.filter(g => g.Vflat >= 120);
const highV5 = highV.filter(g => g.lhOuter !== null);
const hvResult = testRegime('HIGH-VFLAT (>=120)', highV, highV5);

const vhV = ext.filter(g => g.Vflat >= 180);
const vhV5 = vhV.filter(g => g.lhOuter !== null);
const vhResult = testRegime('VERY-HIGH-VFLAT (>=180)', vhV, vhV5);

const lowV = ext.filter(g => g.Vflat < 120);
const lowV5 = lowV.filter(g => g.lhOuter !== null);
const lvResult = testRegime('LOW-VFLAT (<120)', lowV, lowV5);

const q1hv = ext.filter(g => g.Q <= 1 && g.Vflat >= 120);
const q1hv5 = q1hv.filter(g => g.lhOuter !== null);
const q1hvResult = testRegime('Q=1 + VFLAT>=120', q1hv, q1hv5);

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: CHANNEL DOMINANCE — VfResid vs ALL ALTERNATIVES');
console.log('  Internal: VfResid dominates every halo proxy by >=7pp (Phase 132A)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const altChannels = [
  { name: 'VfResid', fn: g => g.VfResid },
  { name: 'logVflat', fn: g => g.logVflat },
  { name: 'logMHI', fn: g => g.logMHI },
  { name: 'logMbar', fn: g => g.logMbar },
  { name: 'logMhost', fn: g => g.logMhost },
  { name: 'logMeanRun', fn: g => g.logMeanRun },
];

const channelResults = {};

function testChannels(label, gals) {
  if (gals.length < 4) return;
  const naive = Math.sqrt(gals.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / gals.length);

  console.log('  ' + label + ' (N=' + gals.length + ')');
  console.log('  ' + 'Channel'.padEnd(15) + 'r(ch,a0)'.padStart(10) + '  Core+ch gap'.padStart(14));
  console.log('  ' + '─'.repeat(42));

  const results = [];
  for (const ch of altChannels) {
    const vals = gals.map(ch.fn);
    const a0s = gals.map(g => g.logA0);
    const r = pearsonR(vals, a0s);

    const chModel = ols(Y_tr, train.map(g => [g.logMHI, g.logMhost, g.logMR, ch.fn({
      VfResid: g.VfResid, logVflat: Math.log10(sparcMap[g.name].Vflat),
      logMHI: g.logMHI, logMbar: Math.log10(sparcMap[g.name].L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9),
      logMhost: tdMap[g.name].logMhost, logMeanRun: g.logMR
    })]));

    const preds = gals.map(g => predictModel(chModel, [g.logMHI, g.logMhost, g.logMeanRun, ch.fn(g)]));
    const rmse = Math.sqrt(gals.reduce((s, g, i) => s + (g.logA0 - preds[i]) ** 2, 0) / gals.length);
    const gap = transferGap(rmse, naive);

    const marker = ch.name === 'VfResid' ? ' <<<' : '';
    console.log('  ' + ch.name.padEnd(15) + r.toFixed(3).padStart(10) + (gap.toFixed(1) + '%').padStart(14) + marker);
    results.push({ name: ch.name, r: +r.toFixed(3), gap: +gap.toFixed(1) });
  }

  const vfrGap = results.find(r => r.name === 'VfResid').gap;
  const bestAlt = results.filter(r => r.name !== 'VfResid').sort((a, b) => b.gap - a.gap)[0];
  console.log('  VfResid margin over best alternative (' + bestAlt.name + '): ' + (vfrGap - bestAlt.gap).toFixed(1) + 'pp');
  console.log('  VfResid dominates: ' + (vfrGap > bestAlt.gap ? 'YES' : 'NO'));
  console.log('');
  channelResults[label] = { results, vfrGap, bestAlt: bestAlt.name, bestAltGap: bestAlt.gap, margin: +(vfrGap - bestAlt.gap).toFixed(1) };
}

testChannels('Full sample', ext);
testChannels('High-Vflat', highV);
testChannels('Very-high-Vflat', vhV);

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: BOOTSTRAP HIERARCHY STABILITY (10000 resamples)');
console.log('  Internal reference: P(VfResid delta>0) = 99.9%');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const BOOT = 10000;

function bootstrapHierarchy(label, gals, gals5) {
  if (gals.length < 6) { console.log('  ' + label + ': insufficient N=' + gals.length); return null; }

  let vfrBeatsCore = 0;
  let lhAddsBeyondVfr = 0;
  let vfrBeatsCoreAndLh = 0;
  let coreFails = 0;
  let vfrPositiveCoeff = 0;
  const vfrDeltas = [];
  const lhDeltas = [];

  for (let b = 0; b < BOOT; b++) {
    const idx = Array.from({ length: gals.length }, () => Math.floor(Math.random() * gals.length));
    const sample = idx.map(i => gals[i]);
    const naive = Math.sqrt(sample.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / sample.length);
    if (naive < 0.01) continue;

    const coreP = sample.map(g => predictModel(coreModel, [g.logMHI, g.logMhost, g.logMeanRun]));
    const vfrP = sample.map(g => predictModel(vfResidModel, [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid]));
    const coreRMSE = Math.sqrt(sample.reduce((s, g, i) => s + (g.logA0 - coreP[i]) ** 2, 0) / sample.length);
    const vfrRMSE = Math.sqrt(sample.reduce((s, g, i) => s + (g.logA0 - vfrP[i]) ** 2, 0) / sample.length);
    const coreGap = transferGap(coreRMSE, naive);
    const vfrGap = transferGap(vfrRMSE, naive);

    if (vfrGap > coreGap) vfrBeatsCore++;
    if (coreGap < 5) coreFails++;
    vfrDeltas.push(vfrGap - coreGap);

    if (gals5.length >= 4) {
      const s5 = idx.map(i => gals5[i % gals5.length]);
      const naive5 = Math.sqrt(s5.reduce((s, g) => s + (g.logA0 - meanY_tr) ** 2, 0) / s5.length);
      if (naive5 > 0.01) {
        const fiveP = s5.map(g => predictModel(fiveAxisModel, [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid, g.lhOuter]));
        const fiveRMSE = Math.sqrt(s5.reduce((s, g, i) => s + (g.logA0 - fiveP[i]) ** 2, 0) / s5.length);
        const fiveGap = transferGap(fiveRMSE, naive5);
        const vfr5P = s5.map(g => predictModel(vfResidModel, [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid]));
        const vfr5RMSE = Math.sqrt(s5.reduce((s, g, i) => s + (g.logA0 - vfr5P[i]) ** 2, 0) / s5.length);
        const vfr5Gap = transferGap(vfr5RMSE, naive5);
        if (fiveGap > vfr5Gap) lhAddsBeyondVfr++;
        lhDeltas.push(fiveGap - vfr5Gap);

        const coreLhP = s5.map(g => predictModel(coreLhModel, [g.logMHI, g.logMhost, g.logMeanRun, g.lhOuter]));
        const coreLhRMSE = Math.sqrt(s5.reduce((s, g, i) => s + (g.logA0 - coreLhP[i]) ** 2, 0) / s5.length);
        const coreLhGap = transferGap(coreLhRMSE, naive5);
        if (vfrGap > coreLhGap) vfrBeatsCoreAndLh++;
      }
    }
  }

  const pVfrBeats = (vfrBeatsCore / BOOT * 100);
  const pCoreFails = (coreFails / BOOT * 100);
  const pLhAdds = lhDeltas.length > 0 ? (lhAddsBeyondVfr / lhDeltas.length * 100) : 0;
  const pVfrDominatesLh = gals5.length >= 4 ? (vfrBeatsCoreAndLh / BOOT * 100) : 0;
  const meanVfrDelta = mean(vfrDeltas);
  const meanLhDelta = lhDeltas.length > 0 ? mean(lhDeltas) : 0;
  const vfrDelta5th = vfrDeltas.sort((a, b) => a - b)[Math.floor(BOOT * 0.05)];
  const vfrDelta95th = vfrDeltas[Math.floor(BOOT * 0.95)];

  console.log('  ' + label + ' (N=' + gals.length + ', ' + BOOT + ' bootstrap resamples):');
  console.log('    P(VfResid > Core):           ' + pVfrBeats.toFixed(1) + '%');
  console.log('    P(Core gap < 5%):            ' + pCoreFails.toFixed(1) + '%');
  console.log('    Mean VfResid delta:          ' + meanVfrDelta.toFixed(1) + 'pp');
  console.log('    VfResid delta 90% CI:        [' + vfrDelta5th.toFixed(1) + ', ' + vfrDelta95th.toFixed(1) + ']pp');
  if (lhDeltas.length > 0) {
    console.log('    P(lhOuter adds > VfResid):   ' + pLhAdds.toFixed(1) + '%');
    console.log('    Mean lhOuter delta:          ' + meanLhDelta.toFixed(1) + 'pp');
  }
  if (gals5.length >= 4) {
    console.log('    P(VfResid > Core+lhOuter):   ' + pVfrDominatesLh.toFixed(1) + '%');
  }
  console.log('');

  return {
    label, N: gals.length,
    pVfrBeatsCore: +pVfrBeats.toFixed(1),
    pCoreFails: +pCoreFails.toFixed(1),
    meanVfrDelta: +meanVfrDelta.toFixed(1),
    vfrDeltaCI: [+vfrDelta5th.toFixed(1), +vfrDelta95th.toFixed(1)],
    pLhAdds: +pLhAdds.toFixed(1),
    meanLhDelta: +meanLhDelta.toFixed(1),
    pVfrDominatesLh: +pVfrDominatesLh.toFixed(1)
  };
}

const bootFull = bootstrapHierarchy('Full sample', ext, ext5);
const bootHV = bootstrapHierarchy('High-Vflat', highV, highV5);
const bootVH = bootstrapHierarchy('Very-high-Vflat', vhV, vhV5);
const bootQ1HV = bootstrapHierarchy('Q=1+Vflat>=120', q1hv, q1hv5);

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: PARTIAL CORRELATIONS — VfResid AFTER CORE REMOVAL');
console.log('  Internal: VfResid partial r = 0.60+ after controlling for Core');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function partialCorrelation(label, gals) {
  if (gals.length < 6) return null;
  const coreResidY = gals.map(g => {
    const pred = predictModel(coreModel, [g.logMHI, g.logMhost, g.logMeanRun]);
    return g.logA0 - pred;
  });

  const coreResidVfr = gals.map(g => {
    const coreVfr = ols(
      train.map(t => t.VfResid),
      train.map(t => [t.logMHI, t.logMhost, t.logMR])
    );
    return g.VfResid - predictModel(coreVfr, [g.logMHI, g.logMhost, g.logMeanRun]);
  });

  const r_partial = pearsonR(coreResidY, coreResidVfr);

  const rVfResid = pearsonR(gals.map(g => g.VfResid), gals.map(g => g.logA0));
  const rLogMHI = pearsonR(coreResidY, gals.map(g => g.logMHI));

  const lhGals = gals.filter(g => g.lhOuter !== null);
  let rLhPartial = null;
  if (lhGals.length >= 6) {
    const coreVfrResidY = lhGals.map(g => {
      const pred = predictModel(vfResidModel, [g.logMHI, g.logMhost, g.logMeanRun, g.VfResid]);
      return g.logA0 - pred;
    });
    rLhPartial = pearsonR(coreVfrResidY, lhGals.map(g => g.lhOuter));
  }

  console.log('  ' + label + ' (N=' + gals.length + '):');
  console.log('    r(VfResid, a0) raw:                 ' + rVfResid.toFixed(3));
  console.log('    r(VfResid, a0 | Core):              ' + r_partial.toFixed(3));
  console.log('    r(logMHI, a0 | Core) check:         ' + rLogMHI.toFixed(3) + ' (should be ~0 if Core removes MHI)');
  if (rLhPartial !== null) {
    console.log('    r(lhOuter, a0 | Core+VfResid):      ' + rLhPartial.toFixed(3));
  }
  console.log('');

  return {
    label, N: gals.length,
    rVfResidRaw: +rVfResid.toFixed(3),
    rVfResidPartial: +r_partial.toFixed(3),
    rLogMHIPartial: +rLogMHI.toFixed(3),
    rLhPartial: rLhPartial !== null ? +rLhPartial.toFixed(3) : null
  };
}

const partFull = partialCorrelation('Full sample', ext);
const partHV = partialCorrelation('High-Vflat', highV);
const partVH = partialCorrelation('Very-high-Vflat', vhV);
const partQ1HV = partialCorrelation('Q=1+Vflat>=120', q1hv);

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 5: INTERNAL vs EXTERNAL COMPARISON TABLE');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  Property'.padEnd(35) + 'Internal (N=45)'.padStart(16) + 'External (N=59)'.padStart(18));
console.log('  ' + '─'.repeat(65));
console.log('  ' + 'Core LOO gap / transfer gap'.padEnd(35) + '44.1%'.padStart(16) + (fullResult.core.gap + '%').padStart(18));
console.log('  ' + 'Core+VfResid LOO / transfer'.padEnd(35) + '61.1%'.padStart(16) + (fullResult.vfResid.gap + '%').padStart(18));
console.log('  ' + '5-axis LOO / transfer'.padEnd(35) + '65.4%'.padStart(16) + (fullResult.fiveAxis ? fullResult.fiveAxis.gap + '%' : 'n/a').padStart(18));
console.log('  ' + 'Core alone fails'.padEnd(35) + 'YES (-11.5% N=10)'.padStart(16) + (fullResult.coreFails ? 'YES' : 'NO').padStart(18));
console.log('  ' + 'VfResid dominates'.padEnd(35) + 'YES'.padStart(16) + (fullResult.vfrDominates ? 'YES' : 'NO').padStart(18));
console.log('  ' + 'r(VfResid, a0) full'.padEnd(35) + '0.690'.padStart(16) + partFull.rVfResidRaw.toFixed(3).padStart(18));
console.log('  ' + 'r(VfResid, a0) high-V'.padEnd(35) + '0.801 (N=10)'.padStart(16) + partHV.rVfResidRaw.toFixed(3).padStart(18));
console.log('  ' + 'P(VfResid > Core) bootstrap'.padEnd(35) + '99.9%'.padStart(16) + (bootFull.pVfrBeatsCore + '%').padStart(18));

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  VERDICT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const hierChecks = {
  coreFails_full: fullResult.coreFails,
  coreFails_hv: hvResult.coreFails,
  vfrDominates_full: fullResult.vfrDominates,
  vfrDominates_hv: hvResult.vfrDominates,
  bootVfrBeats_full: bootFull.pVfrBeatsCore > 90,
  bootVfrBeats_hv: bootHV.pVfrBeatsCore > 90,
  partialSignificant: partHV.rVfResidPartial > 0.3,
  regimePreserved: hvResult.vfResid.gap > fullResult.vfResid.gap,
};

const passed = Object.values(hierChecks).filter(Boolean).length;
const total = Object.values(hierChecks).length;

console.log('  Hierarchy replication checks:');
for (const [k, v] of Object.entries(hierChecks)) {
  console.log('    ' + k.padEnd(30) + (v ? 'PASS' : 'FAIL'));
}
console.log('  Passed: ' + passed + '/' + total);

let verdict;
if (passed >= 7) verdict = 'STRONG_HIERARCHY_REPLICATION';
else if (passed >= 5) verdict = 'MODERATE_HIERARCHY_REPLICATION';
else if (passed >= 3) verdict = 'WEAK_HIERARCHY_REPLICATION';
else verdict = 'HIERARCHY_NOT_REPLICATED';

console.log('\n  VERDICT: ' + verdict);

const lhVerdict = fullResult.fiveAxis && fullResult.fiveAxis.gap > fullResult.vfResid.gap
  ? 'LH_OUTER_ADDS' : 'LH_OUTER_NO_CLEAR_ADDITION';
const hvLhVerdict = hvResult.fiveAxis && hvResult.fiveAxis.gap > hvResult.vfResid.gap
  ? 'LH_OUTER_ADDS_HV' : 'LH_OUTER_NO_CLEAR_ADDITION_HV';
console.log('  lhOuter verdict (full): ' + lhVerdict);
console.log('  lhOuter verdict (high-V): ' + hvLhVerdict);

const output = {
  phase: '202',
  title: 'External Hierarchy Replication',
  verdict,
  hierarchyChecks: hierChecks,
  passed, total,
  regimes: {
    full: fullResult,
    highVflat: hvResult,
    veryHighVflat: vhResult,
    lowVflat: lvResult,
    q1HighVflat: q1hvResult
  },
  channelDominance: channelResults,
  bootstrap: {
    full: bootFull,
    highVflat: bootHV,
    veryHighVflat: bootVH,
    q1HighVflat: bootQ1HV
  },
  partialCorrelations: {
    full: partFull,
    highVflat: partHV,
    veryHighVflat: partVH,
    q1HighVflat: partQ1HV
  },
  lhOuterVerdict: lhVerdict,
  lhOuterVerdictHV: hvLhVerdict,
  internalReference: {
    coreLOO: 44.1, vfResidLOO: 61.1, fiveAxisLOO: 65.4,
    coreTransfer: -11.5, vfResidTransfer: 56.9, fiveAxisTransfer: 66.0,
    rVfResidHoldout: 0.801, bootPVfrBeats: 99.9
  }
};

fs.writeFileSync(pub('phase202-hierarchy-replication.json'), JSON.stringify(output, null, 2));
console.log('\n  Output: public/phase202-hierarchy-replication.json');
console.log('\n======================================================================');
console.log('  PHASE 202 COMPLETE');
console.log('======================================================================');

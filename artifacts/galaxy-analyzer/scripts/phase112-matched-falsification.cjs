#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 112: MATCHED FALSIFICATION OF GAS-DOMINANCE HYPOTHESIS');
console.log('  Does the signal survive when we control for confounders');
console.log('  via matching? Which variable is truly fundamental?');
console.log('  Testing: fgas, log(MHI/L36), logSigma0, logBaryonCompact');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

let rng = 2026;
function rand() { rng = (rng * 1664525 + 1013904223) & 0x7fffffff; return rng / 0x7fffffff; }
function shuffle(arr) {
  const a = [...arr];
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(rand() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
}

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function percentile(a, p) { const s = [...a].sort((x, y) => x - y); const i = p * (s.length - 1); const lo = Math.floor(i); return lo >= s.length - 1 ? s[s.length - 1] : s[lo] + (i - lo) * (s[lo + 1] - s[lo]); }

function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function olsRegress(X, y) {
  const n = y.length, p = X[0].length;
  if (n <= p + 1) return null;
  const Xt = Array.from({ length: p }, (_, i) => X.map(row => row[i]));
  const XtX = Array.from({ length: p }, (_, i) =>
    Array.from({ length: p }, (_, j) =>
      Xt[i].reduce((s, _, k) => s + Xt[i][k] * Xt[j][k], 0)));
  const Xty = Xt.map(col => col.reduce((s, v, k) => s + v * y[k], 0));
  const aug = XtX.map((row, i) => [...row, Xty[i]]);
  for (let i = 0; i < p; i++) {
    let maxRow = i;
    for (let k = i + 1; k < p; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) maxRow = k;
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];
    if (Math.abs(aug[i][i]) < 1e-12) return null;
    for (let k = 0; k < p; k++) {
      if (k === i) continue;
      const f = aug[k][i] / aug[i][i];
      for (let j = i; j <= p; j++) aug[k][j] -= f * aug[i][j];
    }
  }
  const beta = aug.map((row, i) => row[p] / row[i]);
  return { beta };
}

function looR2(xArr, yArr) {
  const n = xArr.length;
  if (n < 8) return NaN;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(yArr);
  for (let i = 0; i < n; i++) {
    const xTr = xArr.filter((_, j) => j !== i);
    const yTr = yArr.filter((_, j) => j !== i);
    const fit = olsRegress(xTr.map(v => [1, v]), yTr);
    if (!fit) return NaN;
    ss_loo += (yArr[i] - (fit.beta[0] + fit.beta[1] * xArr[i])) ** 2;
    ss_tot += (yArr[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_loo / ss_tot : 0;
}

function residualize(target, predictors) {
  const X = predictors[0].map((_, i) => [1, ...predictors.map(p => p[i])]);
  const fit = olsRegress(X, target);
  return fit ? target.map((v, i) => v - X[i].reduce((s, xv, j) => s + xv * fit.beta[j], 0)) : target;
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  table1[name] = { L36, Rdisk, MHI, Vflat };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul });
}

const galaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat < 70) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    pts.push({ r: pt.rad, vobs: pt.vobs, vBar: Math.sqrt(Math.max(vBarSq, 0.01)) });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);
  const outerMD = outerPts.map(p => (p.vobs ** 2) / (p.vBar ** 2)).filter(v => isFinite(v) && v > 0);
  const logOMD = Math.log10(mean(outerMD));
  if (!isFinite(logOMD)) continue;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1);
  const baryonCompact = t1.L36 * UPSILON_DISK / (t1.Rdisk > 0 ? t1.Rdisk : 1);

  if (!isFinite(fgas) || !isFinite(Sigma0) || Sigma0 <= 0) continue;

  galaxies.push({
    name, logOMD, fgas,
    logMHI_L36: Math.log10(t1.MHI) - Math.log10(t1.L36),
    logSigma0: Math.log10(Sigma0),
    logL36: Math.log10(t1.L36),
    logRdisk: Math.log10(t1.Rdisk > 0 ? t1.Rdisk : 0.1),
    logMHI: Math.log10(t1.MHI),
    logBaryonCompact: Math.log10(baryonCompact > 0 ? baryonCompact : 1e-3),
    Vflat: t1.Vflat,
  });
}

const N = galaxies.length;
console.log('  N = ' + N + ' high-Vflat galaxies\n');

const predictors = {
  fgas: galaxies.map(g => g.fgas),
  'log(MHI/L36)': galaxies.map(g => g.logMHI_L36),
  logSigma0: galaxies.map(g => g.logSigma0),
  logBaryonCompact: galaxies.map(g => g.logBaryonCompact),
};
const yArr = galaxies.map(g => g.logOMD);

console.log('  BASELINE (unmatched):');
for (const [pName, pArr] of Object.entries(predictors)) {
  const r = pearsonR(pArr, yArr);
  const loo = looR2(pArr, yArr);
  console.log('    ' + pName + ': r=' + r.toFixed(3) + ', LOO=' + loo.toFixed(3));
}

// ============================================================
// TEST 1: NEAREST-NEIGHBOR MATCHING ON CONFOUNDERS
// For each predictor X, match galaxies on the OTHER variables,
// then measure within-matched-pair residual correlation
// ============================================================

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 1: NEAREST-NEIGHBOR MATCHING');
console.log('  Match on confounders, then test residual signal');
console.log('══════════════════════════════════════════════════════════════\n');

function nearestNeighborMatch(galaxies, matchVars, testVar, target) {
  const matchSDs = matchVars.map(v => sd(galaxies.map(g => g[v])));
  const usedIndices = new Set();
  const pairs = [];

  for (let i = 0; i < galaxies.length; i++) {
    let bestJ = -1, bestDist = Infinity;
    for (let j = 0; j < galaxies.length; j++) {
      if (i === j || usedIndices.has(i) || usedIndices.has(j)) continue;
      let dist = 0;
      for (let k = 0; k < matchVars.length; k++) {
        const d = (galaxies[i][matchVars[k]] - galaxies[j][matchVars[k]]) / (matchSDs[k] || 1);
        dist += d * d;
      }
      dist = Math.sqrt(dist);
      if (dist < bestDist) { bestDist = dist; bestJ = j; }
    }
    if (bestJ >= 0 && bestDist < 1.5) {
      usedIndices.add(i);
      usedIndices.add(bestJ);
      pairs.push([i, bestJ]);
    }
  }

  if (pairs.length < 8) return { r: NaN, nPairs: pairs.length, diffs: [] };

  const testDiffs = pairs.map(([i, j]) => galaxies[i][testVar] - galaxies[j][testVar]);
  const targetDiffs = pairs.map(([i, j]) => galaxies[i][target] - galaxies[j][target]);

  const r = pearsonR(testDiffs, targetDiffs);
  return { r, nPairs: pairs.length, matchQuality: mean(pairs.map(([i, j]) => {
    let d = 0;
    for (let k = 0; k < matchVars.length; k++) {
      const dd = (galaxies[i][matchVars[k]] - galaxies[j][matchVars[k]]) / (matchSDs[k] || 1);
      d += dd * dd;
    }
    return Math.sqrt(d);
  })) };
}

const matchingTests = [
  { name: 'Match on logSigma0 + logL36', matchOn: ['logSigma0', 'logL36'] },
  { name: 'Match on logSigma0 + logL36 + Vflat', matchOn: ['logSigma0', 'logL36', 'Vflat'] },
  { name: 'Match on logBaryonCompact + logRdisk', matchOn: ['logBaryonCompact', 'logRdisk'] },
  { name: 'Match on logSigma0 + logBaryonCompact', matchOn: ['logSigma0', 'logBaryonCompact'] },
  { name: 'Match on logL36 + logRdisk + Vflat', matchOn: ['logL36', 'logRdisk', 'Vflat'] },
];

const matchResults = {};

for (const test of matchingTests) {
  console.log('  ' + test.name + ':');
  const testVars = ['fgas', 'logMHI_L36', 'logSigma0', 'logBaryonCompact'];
  const results = {};

  for (const tv of testVars) {
    if (test.matchOn.includes(tv)) {
      console.log('    ' + tv + ': (matched on this — skip)');
      results[tv] = { r: NaN, nPairs: 0, note: 'matched on this variable' };
      continue;
    }
    const res = nearestNeighborMatch(galaxies, test.matchOn, tv, 'logOMD');
    console.log('    ' + tv + ': r=' + (isNaN(res.r) ? 'N/A' : res.r.toFixed(3)) + ' (' + res.nPairs + ' pairs, match quality=' + (res.matchQuality ? res.matchQuality.toFixed(3) : 'N/A') + ')');
    results[tv] = { r: isNaN(res.r) ? null : parseFloat(res.r.toFixed(3)), nPairs: res.nPairs };
  }
  console.log('');
  matchResults[test.name] = results;
}

// ============================================================
// TEST 2: WITHIN-BIN ANALYSIS
// Divide into bins on confounder, compute correlation within
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: WITHIN-BIN ANALYSIS');
console.log('  Divide into terciles on confounder, measure r within each');
console.log('══════════════════════════════════════════════════════════════\n');

function withinBinAnalysis(galaxies, binVar, testVars, nBins) {
  const sorted = [...galaxies].sort((a, b) => a[binVar] - b[binVar]);
  const binSize = Math.floor(sorted.length / nBins);
  const bins = [];
  for (let b = 0; b < nBins; b++) {
    const start = b * binSize;
    const end = b === nBins - 1 ? sorted.length : (b + 1) * binSize;
    bins.push(sorted.slice(start, end));
  }

  const results = {};
  for (const tv of testVars) {
    const binRs = [];
    for (let b = 0; b < bins.length; b++) {
      const bin = bins[b];
      const x = bin.map(g => g[tv]);
      const y = bin.map(g => g.logOMD);
      const r = bin.length >= 8 ? pearsonR(x, y) : NaN;
      binRs.push({ binIdx: b, n: bin.length, r: isNaN(r) ? null : parseFloat(r.toFixed(3)),
        binRange: [bin[0][binVar].toFixed(2), bin[bin.length - 1][binVar].toFixed(2)] });
    }
    results[tv] = binRs;
  }
  return results;
}

const binConfounders = ['logSigma0', 'logL36', 'Vflat'];
const testVarNames = ['fgas', 'logMHI_L36', 'logSigma0', 'logBaryonCompact'];
const binResults = {};

for (const binVar of binConfounders) {
  console.log('  Bins on ' + binVar + ' (terciles):');
  const res = withinBinAnalysis(galaxies, binVar, testVarNames, 3);
  binResults[binVar] = res;

  for (const tv of testVarNames) {
    const line = res[tv].map(b => 'r=' + (b.r !== null ? b.r.toFixed(3) : 'N/A') + '(N=' + b.n + ')').join(' | ');
    console.log('    ' + tv + ': ' + line);
  }
  console.log('');
}

// ============================================================
// TEST 3: WITHIN-BIN SHUFFLE (PERMUTATION)
// Shuffle predictor within bins; if real r >> null r, signal is genuine
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: WITHIN-BIN SHUFFLE (PERMUTATION)');
console.log('  Shuffle predictor within bins of confounder.');
console.log('  If real signal >> shuffled, it is NOT confounded.');
console.log('══════════════════════════════════════════════════════════════\n');

const N_PERM = 1000;

function withinBinShuffle(galaxies, binVar, testVar, nBins, nPerm) {
  const sorted = [...galaxies].sort((a, b) => a[binVar] - b[binVar]);
  const binSize = Math.floor(sorted.length / nBins);

  const realX = sorted.map(g => g[testVar]);
  const realY = sorted.map(g => g.logOMD);
  const realR = Math.abs(pearsonR(realX, realY));

  let countGE = 0;
  const nullRs = [];
  for (let p = 0; p < nPerm; p++) {
    const shuffledX = new Array(sorted.length);
    for (let b = 0; b < nBins; b++) {
      const start = b * binSize;
      const end = b === nBins - 1 ? sorted.length : (b + 1) * binSize;
      const binValues = [];
      for (let i = start; i < end; i++) binValues.push(sorted[i][testVar]);
      const shuffled = shuffle(binValues);
      for (let i = 0; i < shuffled.length; i++) shuffledX[start + i] = shuffled[i];
    }
    const nullR = Math.abs(pearsonR(shuffledX, realY));
    nullRs.push(nullR);
    if (nullR >= realR) countGE++;
  }

  return {
    realR: parseFloat(realR.toFixed(3)),
    nullMean: parseFloat(mean(nullRs).toFixed(3)),
    null95: parseFloat(percentile(nullRs, 0.95).toFixed(3)),
    pValue: parseFloat(((countGE + 1) / (nPerm + 1)).toFixed(4)),
  };
}

const shuffleConfounders = ['logSigma0', 'logL36', 'Vflat', 'logBaryonCompact'];
const shuffleResults = {};

for (const binVar of shuffleConfounders) {
  console.log('  Shuffle within bins of ' + binVar + ' (4 bins, ' + N_PERM + ' permutations):');
  shuffleResults[binVar] = {};

  for (const tv of testVarNames) {
    if (tv === binVar) continue;
    const res = withinBinShuffle(galaxies, binVar, tv, 4, N_PERM);
    const survived = res.pValue < 0.05 ? 'SURVIVES' : 'FAILS';
    console.log('    ' + tv + ': |r|=' + res.realR + ', null mean=' + res.nullMean + ', null 95%=' + res.null95 + ', p=' + res.pValue + ' ' + survived);
    shuffleResults[binVar][tv] = res;
  }
  console.log('');
}

// ============================================================
// TEST 4: DOUBLE MATCHING — Sigma0 + luminosity simultaneously
// The strongest possible confounder control
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 4: DOUBLE MATCHING — logSigma0 + logL36');
console.log('  The strongest confounder control.');
console.log('  Match pairs on BOTH Sigma0 and luminosity,');
console.log('  then ask: does the predictor difference still');
console.log('  predict the logOMD difference?');
console.log('══════════════════════════════════════════════════════════════\n');

const doubleMatchVars = ['logSigma0', 'logL36'];

for (const tv of ['fgas', 'logMHI_L36']) {
  const res = nearestNeighborMatch(galaxies, doubleMatchVars, tv, 'logOMD');
  console.log('  ' + tv + ' after matching on Sigma0+L36:');
  console.log('    r=' + (isNaN(res.r) ? 'N/A' : res.r.toFixed(3)) + ', nPairs=' + res.nPairs + ', match quality=' + (res.matchQuality ? res.matchQuality.toFixed(3) : 'N/A'));

  if (res.nPairs >= 8 && !isNaN(res.r)) {
    let nullCount = 0;
    for (let p = 0; p < N_PERM; p++) {
      const pairs = [];
      const usedI = new Set();
      for (let i = 0; i < galaxies.length; i++) {
        if (usedI.has(i)) continue;
        let bestJ = -1, bestDist = Infinity;
        const matchSDs = doubleMatchVars.map(v => sd(galaxies.map(g => g[v])));
        for (let j = i + 1; j < galaxies.length; j++) {
          if (usedI.has(j)) continue;
          let dist = 0;
          for (let k = 0; k < doubleMatchVars.length; k++) {
            const d = (galaxies[i][doubleMatchVars[k]] - galaxies[j][doubleMatchVars[k]]) / (matchSDs[k] || 1);
            dist += d * d;
          }
          dist = Math.sqrt(dist);
          if (dist < bestDist) { bestDist = dist; bestJ = j; }
        }
        if (bestJ >= 0 && bestDist < 1.5) {
          usedI.add(i); usedI.add(bestJ);
          pairs.push([i, bestJ]);
        }
      }
      if (pairs.length < 8) continue;
      const testDiffs = pairs.map(([i, j]) => galaxies[i][tv] - galaxies[j][tv]);
      const targetDiffs = shuffle(pairs.map(([i, j]) => galaxies[i].logOMD - galaxies[j].logOMD));
      const nullR = Math.abs(pearsonR(testDiffs, targetDiffs));
      if (nullR >= Math.abs(res.r)) nullCount++;
    }
    const pVal = (nullCount + 1) / (N_PERM + 1);
    console.log('    Permutation p-value: ' + pVal.toFixed(4) + (pVal < 0.05 ? ' SIGNIFICANT' : ' NOT significant'));
  }
  console.log('');
}

// ============================================================
// TEST 5: PARTIAL CORRELATION CASCADE — definitive hierarchy
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 5: PARTIAL CORRELATION CASCADE');
console.log('  Residualize each variable against progressively more controls');
console.log('══════════════════════════════════════════════════════════════\n');

const controlSets = [
  { name: '| logSigma0', controls: ['logSigma0'] },
  { name: '| logSigma0 + logL36', controls: ['logSigma0', 'logL36'] },
  { name: '| logSigma0 + logL36 + Vflat', controls: ['logSigma0', 'logL36', 'Vflat'] },
  { name: '| logSigma0 + logBaryonCompact', controls: ['logSigma0', 'logBaryonCompact'] },
  { name: '| ALL (Sig0+L36+Rdisk+Vflat+Compact)', controls: ['logSigma0', 'logL36', 'logRdisk', 'Vflat', 'logBaryonCompact'] },
];

const cascadeResults = {};
const varMap = {
  fgas: galaxies.map(g => g.fgas),
  'log(MHI/L36)': galaxies.map(g => g.logMHI_L36),
  logSigma0: galaxies.map(g => g.logSigma0),
  logBaryonCompact: galaxies.map(g => g.logBaryonCompact),
  logL36: galaxies.map(g => g.logL36),
  logRdisk: galaxies.map(g => g.logRdisk),
  Vflat: galaxies.map(g => g.Vflat),
};

for (const cs of controlSets) {
  console.log('  Controls: ' + cs.name);
  cascadeResults[cs.name] = {};
  const controlArrs = cs.controls.map(c => varMap[c]);
  const yPerp = residualize(yArr, controlArrs);

  for (const tv of ['fgas', 'log(MHI/L36)', 'logSigma0', 'logBaryonCompact']) {
    if (cs.controls.includes(tv)) {
      console.log('    ' + tv + ': (controlled — skip)');
      cascadeResults[cs.name][tv] = { r: null, note: 'controlled' };
      continue;
    }
    const xPerp = residualize(varMap[tv], controlArrs);
    const r = pearsonR(xPerp, yPerp);
    console.log('    ' + tv + ': partial r = ' + r.toFixed(3));
    cascadeResults[cs.name][tv] = { r: parseFloat(r.toFixed(3)) };
  }
  console.log('');
}

// ============================================================
// VERDICT
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  OVERALL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const key1 = matchResults['Match on logSigma0 + logL36'];
const fgasMatch = key1 && key1.fgas && key1.fgas.r !== null ? key1.fgas.r : NaN;
const ratioMatch = key1 && key1.logMHI_L36 && key1.logMHI_L36.r !== null ? key1.logMHI_L36.r : NaN;
const sig0Match = key1 && key1.logSigma0 ? key1.logSigma0.r : NaN;

const cascadeAll = cascadeResults['| ALL (Sig0+L36+Rdisk+Vflat+Compact)'];
const fgasCascade = cascadeAll && cascadeAll.fgas ? cascadeAll.fgas.r : NaN;
const ratioCascade = cascadeAll && cascadeAll['log(MHI/L36)'] ? cascadeAll['log(MHI/L36)'].r : NaN;

console.log('  MATCHING SUMMARY (Sigma0 + L36 matched):');
console.log('    fgas:           r=' + (isNaN(fgasMatch) ? 'N/A' : fgasMatch.toFixed(3)));
console.log('    log(MHI/L36):   r=' + (isNaN(ratioMatch) ? 'N/A' : ratioMatch.toFixed(3)));
console.log('    logSigma0:      (matched on)');
console.log('');

console.log('  CASCADE SUMMARY (all controls):');
console.log('    fgas:           partial r=' + (isNaN(fgasCascade) ? 'N/A' : fgasCascade.toFixed(3)));
console.log('    log(MHI/L36):   partial r=' + (isNaN(ratioCascade) ? 'N/A' : ratioCascade.toFixed(3)));
console.log('');

const shuffleSig0_fgas = shuffleResults.logSigma0 && shuffleResults.logSigma0.fgas ? shuffleResults.logSigma0.fgas.pValue : 1;
const shuffleSig0_ratio = shuffleResults.logSigma0 && shuffleResults.logSigma0['logMHI_L36'] ? shuffleResults.logSigma0['logMHI_L36'].pValue : 1;
const shuffleL36_fgas = shuffleResults.logL36 && shuffleResults.logL36.fgas ? shuffleResults.logL36.fgas.pValue : 1;

console.log('  PERMUTATION SUMMARY (shuffle within Sigma0 bins):');
console.log('    fgas:           p=' + shuffleSig0_fgas.toFixed(4));
console.log('    log(MHI/L36):   p=' + shuffleSig0_ratio.toFixed(4));
console.log('');

let winner, verdict;
const fgasSurvived = !isNaN(fgasMatch) && Math.abs(fgasMatch) > 0.2;
const ratioSurvived = !isNaN(ratioMatch) && Math.abs(ratioMatch) > 0.2;
const fgasPermSurvived = shuffleSig0_fgas < 0.05;
const ratioPermSurvived = shuffleSig0_ratio < 0.05;

if (fgasSurvived && ratioSurvived) {
  if (Math.abs(ratioMatch) > Math.abs(fgasMatch) + 0.05) {
    winner = 'log(MHI/L36)';
    verdict = 'STRONG PASS: Both fgas and log(MHI/L36) survive matched falsification. log(MHI/L36) is the stronger survivor — the gas-to-stellar RATIO is the more fundamental variable.';
  } else if (Math.abs(fgasMatch) > Math.abs(ratioMatch) + 0.05) {
    winner = 'fgas';
    verdict = 'STRONG PASS: Both survive, but fgas is the stronger survivor after matching.';
  } else {
    winner = 'tied (fgas ~ log(MHI/L36))';
    verdict = 'STRONG PASS: Both fgas and log(MHI/L36) survive equally. They carry essentially the same independent signal — the gas-to-stellar balance.';
  }
} else if (fgasSurvived && !ratioSurvived) {
  winner = 'fgas';
  verdict = 'MODERATE PASS: fgas survives but log(MHI/L36) does not. fgas form is preferred.';
} else if (!fgasSurvived && ratioSurvived) {
  winner = 'log(MHI/L36)';
  verdict = 'MODERATE PASS: log(MHI/L36) survives but fgas does not. Ratio form is more fundamental.';
} else {
  winner = 'none';
  verdict = 'FAIL: Neither fgas nor log(MHI/L36) survive matched falsification. Signal may be fully confounded.';
}

console.log('  WINNER: ' + winner);
console.log('  VERDICT: ' + verdict);

const output = {
  phase: 112,
  title: 'Matched Falsification of Gas-Dominance Hypothesis',
  n: N,
  baseline: Object.fromEntries(Object.entries(predictors).map(([k, v]) => [k, { r: parseFloat(pearsonR(v, yArr).toFixed(3)), loo: parseFloat(looR2(v, yArr).toFixed(3)) }])),
  matchResults,
  binResults,
  shuffleResults,
  cascadeResults,
  winner,
  verdict
};

const outPath = path.join(__dirname, '..', 'public', 'phase112-matched-falsification.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);

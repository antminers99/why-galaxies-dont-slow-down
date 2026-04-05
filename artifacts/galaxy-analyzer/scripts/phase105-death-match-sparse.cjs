#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 105: DEATH MATCH — SPARSE MODEL COMPETITION');
console.log('  Which single variable is the TRUE backbone predictor?');
console.log('  7 models, 6 metrics each. No mercy.');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }

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
  const yhat = X.map(row => row.reduce((s, v, j) => s + v * beta[j], 0));
  const ss_res = y.reduce((s, v, i) => s + (v - yhat[i]) ** 2, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, ss_res, yhat, n, p };
}

function looR2(data, varFns) {
  const n = data.length;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(data.map(g => g.logOMD));
  for (let i = 0; i < n; i++) {
    const train = data.filter((_, j) => j !== i);
    const test = data[i];
    const X = train.map(g => [1, ...varFns.map(fn => fn(g))]);
    const y = train.map(g => g.logOMD);
    const fit = olsRegress(X, y);
    if (!fit) return NaN;
    const pred = fit.beta[0] + varFns.reduce((s, fn, j) => s + fit.beta[j + 1] * fn(test), 0);
    ss_loo += (test.logOMD - pred) ** 2;
    ss_tot += (test.logOMD - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_loo / ss_tot : 0;
}

function kFoldR2(data, varFns, k, seed) {
  let rng = seed;
  function rand() { rng = (rng * 1664525 + 1013904223) & 0x7fffffff; return rng / 0x7fffffff; }
  const indices = data.map((_, i) => i);
  for (let i = indices.length - 1; i > 0; i--) {
    const j = Math.floor(rand() * (i + 1));
    [indices[i], indices[j]] = [indices[j], indices[i]];
  }
  const folds = Array.from({ length: k }, () => []);
  indices.forEach((idx, i) => folds[i % k].push(idx));

  let ss_pred = 0, ss_tot = 0;
  const yMean = mean(data.map(g => g.logOMD));

  for (let f = 0; f < k; f++) {
    const testIdx = new Set(folds[f]);
    const train = data.filter((_, i) => !testIdx.has(i));
    const test = data.filter((_, i) => testIdx.has(i));
    const X = train.map(g => [1, ...varFns.map(fn => fn(g))]);
    const y = train.map(g => g.logOMD);
    const fit = olsRegress(X, y);
    if (!fit) return NaN;
    for (const g of test) {
      const pred = fit.beta[0] + varFns.reduce((s, fn, j) => s + fit.beta[j + 1] * fn(g), 0);
      ss_pred += (g.logOMD - pred) ** 2;
      ss_tot += (g.logOMD - yMean) ** 2;
    }
  }
  return ss_tot > 0 ? 1 - ss_pred / ss_tot : 0;
}

function aicBic(n, k, ss_res) {
  const logLik = -n / 2 * Math.log(2 * Math.PI * ss_res / n) - n / 2;
  const aic = -2 * logLik + 2 * (k + 1);
  const bic = -2 * logLik + Math.log(n) * (k + 1);
  return { aic, bic };
}

function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function partialCorr(x, y, controls) {
  function residualize(target, predictors) {
    const X = predictors[0].map((_, i) => [1, ...predictors.map(p => p[i])]);
    const fit = olsRegress(X, target);
    return fit ? fit.yhat.map((v, i) => target[i] - v) : target;
  }
  const xResid = residualize(x, controls);
  const yResid = residualize(y, controls);
  return pearsonR(xResid, yResid);
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

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const baryonCompact = t1.L36 * UPSILON_DISK / (t1.Rdisk > 0 ? t1.Rdisk : 1);

  if (!isFinite(fgas) || !isFinite(Sigma0) || Sigma0 <= 0) continue;

  galaxies.push({
    name, logOMD, fgas,
    logSigma0: Math.log10(Sigma0),
    logBaryonCompact: Math.log10(baryonCompact > 0 ? baryonCompact : 1e-3),
    logL36: Math.log10(t1.L36),
    Vflat: t1.Vflat
  });
}

console.log('  N = ' + galaxies.length + ' high-Vflat galaxies\n');

const models = [
  { name: 'fgas only', vars: [g => g.fgas], k: 1 },
  { name: 'logSigma0 only', vars: [g => g.logSigma0], k: 1 },
  { name: 'logBaryonCompact only', vars: [g => g.logBaryonCompact], k: 1 },
  { name: 'logL36 only', vars: [g => g.logL36], k: 1 },
  { name: 'fgas + logSigma0', vars: [g => g.fgas, g => g.logSigma0], k: 2 },
  { name: 'fgas + logBaryonCompact', vars: [g => g.fgas, g => g.logBaryonCompact], k: 2 },
  { name: 'fgas + logL36', vars: [g => g.fgas, g => g.logL36], k: 2 },
];

console.log('══════════════════════════════════════════════════════════════');
console.log('  PART 1: MODEL COMPARISON TABLE');
console.log('══════════════════════════════════════════════════════════════\n');

const y = galaxies.map(g => g.logOMD);
const results = [];

const nReps5 = 50;
const nReps10 = 50;

for (const model of models) {
  const X = galaxies.map(g => [1, ...model.vars.map(fn => fn(g))]);
  const fit = olsRegress(X, y);
  const r2 = fit.r2;

  const loo = looR2(galaxies, model.vars);

  let cv5s = [];
  for (let s = 0; s < nReps5; s++) cv5s.push(kFoldR2(galaxies, model.vars, 5, 42 + s * 7));
  const cv5 = mean(cv5s);

  let cv10s = [];
  for (let s = 0; s < nReps10; s++) cv10s.push(kFoldR2(galaxies, model.vars, 10, 99 + s * 13));
  const cv10 = mean(cv10s);

  const { aic, bic } = aicBic(galaxies.length, model.k, fit.ss_res);

  const res = { name: model.name, k: model.k, r2, loo, cv5, cv10, aic, bic };
  results.push(res);

  console.log('  ' + model.name);
  console.log('    R²=' + r2.toFixed(4) + '  LOO=' + loo.toFixed(4) + '  5-fold=' + cv5.toFixed(4) + '  10-fold=' + cv10.toFixed(4));
  console.log('    AIC=' + aic.toFixed(1) + '  BIC=' + bic.toFixed(1));
  console.log('    coefficients: [' + fit.beta.map(b => b.toFixed(4)).join(', ') + ']');
  console.log('');
}

console.log('══════════════════════════════════════════════════════════════');
console.log('  PART 2: PARTIAL CORRELATIONS');
console.log('  Each predictor controlling for each other');
console.log('══════════════════════════════════════════════════════════════\n');

const varMap = {
  fgas: galaxies.map(g => g.fgas),
  logSigma0: galaxies.map(g => g.logSigma0),
  logBaryonCompact: galaxies.map(g => g.logBaryonCompact),
  logL36: galaxies.map(g => g.logL36)
};
const yArr = galaxies.map(g => g.logOMD);

const partials = {};
for (const [xName, xArr] of Object.entries(varMap)) {
  partials[xName] = {};
  for (const [cName, cArr] of Object.entries(varMap)) {
    if (cName === xName) continue;
    const pr = partialCorr(xArr, yArr, [cArr]);
    partials[xName][cName] = pr;
  }
  const allOthers = Object.entries(varMap).filter(([n]) => n !== xName).map(([, a]) => a);
  partials[xName]['all_others'] = partialCorr(xArr, yArr, allOthers);
}

console.log('  r(X, logOMD | control):');
console.log('  ' + ''.padEnd(22) + 'raw'.padEnd(8) + '|Sig0'.padEnd(8) + '|Comp'.padEnd(8) + '|L36'.padEnd(8) + '|ALL'.padEnd(8));
for (const [xName, xArr] of Object.entries(varMap)) {
  const raw = pearsonR(xArr, yArr);
  const line = '  ' + xName.padEnd(22) +
    raw.toFixed(3).padEnd(8) +
    (partials[xName].logSigma0 || 0).toFixed(3).padEnd(8) +
    (partials[xName].logBaryonCompact || 0).toFixed(3).padEnd(8) +
    (partials[xName].logL36 || 0).toFixed(3).padEnd(8) +
    (partials[xName].all_others || 0).toFixed(3).padEnd(8);
  console.log(line);
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART 3: RESIDUAL IMPROVEMENT OVER fgas');
console.log('  Does adding X to fgas improve prediction?');
console.log('══════════════════════════════════════════════════════════════\n');

const fgasOnly = results.find(r => r.name === 'fgas only');

for (const r of results) {
  if (r.name === 'fgas only') continue;
  const deltaLOO = r.loo - fgasOnly.loo;
  const deltaAIC = r.aic - fgasOnly.aic;
  const deltaBIC = r.bic - fgasOnly.bic;
  const worth = deltaLOO > 0.02 && deltaAIC < -2;
  console.log('  ' + r.name + ':');
  console.log('    delta LOO = ' + (deltaLOO > 0 ? '+' : '') + deltaLOO.toFixed(4));
  console.log('    delta AIC = ' + (deltaAIC > 0 ? '+' : '') + deltaAIC.toFixed(1));
  console.log('    delta BIC = ' + (deltaBIC > 0 ? '+' : '') + deltaBIC.toFixed(1));
  console.log('    Worth adding? ' + (worth ? 'YES' : 'NO') + '\n');
}

console.log('══════════════════════════════════════════════════════════════');
console.log('  PART 4: HEAD-TO-HEAD KNOCKOUT (70/30 splits)');
console.log('══════════════════════════════════════════════════════════════\n');

const nSplits = 500;
let rng = 42;
function rand() { rng = (rng * 1664525 + 1013904223) & 0x7fffffff; return rng / 0x7fffffff; }

const singleModels = models.filter(m => m.k === 1);
const wins = {};
singleModels.forEach(m => wins[m.name] = 0);

for (let iter = 0; iter < nSplits; iter++) {
  const shuffled = [...galaxies];
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(rand() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  const cut = Math.floor(shuffled.length * 0.7);
  const train = shuffled.slice(0, cut);
  const test = shuffled.slice(cut);

  let bestName = '';
  let bestR2 = -Infinity;
  for (const m of singleModels) {
    const Xtr = train.map(g => [1, ...m.vars.map(fn => fn(g))]);
    const ytr = train.map(g => g.logOMD);
    const fit = olsRegress(Xtr, ytr);
    if (!fit) continue;
    const preds = test.map(g => fit.beta[0] + m.vars.reduce((s, fn, j) => s + fit.beta[j + 1] * fn(g), 0));
    const actuals = test.map(g => g.logOMD);
    const ss_r = actuals.reduce((s, v, i) => s + (v - preds[i]) ** 2, 0);
    const ss_t = actuals.reduce((s, v) => s + (v - mean(actuals)) ** 2, 0);
    const exR2 = ss_t > 0 ? 1 - ss_r / ss_t : 0;
    if (exR2 > bestR2) { bestR2 = exR2; bestName = m.name; }
  }
  if (bestName) wins[bestName]++;
}

console.log('  Winner frequency (' + nSplits + ' random 70/30 splits):');
for (const [name, w] of Object.entries(wins)) {
  console.log('    ' + name.padEnd(25) + ': ' + w + ' wins (' + (w / nSplits * 100).toFixed(1) + '%)');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART 5: PAIRWISE CORRELATION MATRIX');
console.log('══════════════════════════════════════════════════════════════\n');

const varNames = Object.keys(varMap);
console.log('  ' + ''.padEnd(22) + varNames.map(n => n.substring(0, 7).padEnd(8)).join(''));
for (const n1 of varNames) {
  let line = '  ' + n1.padEnd(22);
  for (const n2 of varNames) {
    line += pearsonR(varMap[n1], varMap[n2]).toFixed(3).padEnd(8);
  }
  console.log(line);
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const bestSingle = results.filter(r => r.k === 1).sort((a, b) => b.loo - a.loo)[0];
const bestPair = results.filter(r => r.k === 2).sort((a, b) => b.loo - a.loo)[0];
const bestOverall = results.sort((a, b) => {
  const aScore = a.loo - a.k * 0.001;
  const bScore = b.loo - b.k * 0.001;
  return bScore - aScore;
})[0];

const pairImprovement = bestPair.loo - bestSingle.loo;
const pairWorth = pairImprovement > 0.02 && bestPair.aic < bestSingle.aic;

console.log('  Best single: ' + bestSingle.name + ' (LOO=' + bestSingle.loo.toFixed(4) + ')');
console.log('  Best pair:   ' + bestPair.name + ' (LOO=' + bestPair.loo.toFixed(4) + ')');
console.log('  Pair improvement over best single: ' + (pairImprovement > 0 ? '+' : '') + pairImprovement.toFixed(4));
console.log('  Pair worth adding (LOO>+0.02 AND AIC lower)? ' + (pairWorth ? 'YES' : 'NO'));
console.log('');

const fgasPartialAfterAll = partials.fgas.all_others;
const strongBackbone = bestSingle.name === 'fgas only' && Math.abs(fgasPartialAfterAll) > 0.3;

let verdict;
if (strongBackbone && !pairWorth) {
  verdict = 'STRONG BACKBONE: fgas is the irreducible backbone predictor. No second variable adds meaningful improvement.';
} else if (strongBackbone && pairWorth) {
  verdict = 'BACKBONE + SUPPLEMENT: fgas is backbone but a second variable adds genuine (small) improvement.';
} else if (!strongBackbone && bestSingle.name !== 'fgas only') {
  verdict = 'DETHRONED: ' + bestSingle.name + ' beats fgas. The backbone predictor is not gas fraction.';
} else {
  verdict = 'AMBIGUOUS: fgas wins nominally but partial correlation drops below 0.3 after controlling competitors.';
}

console.log('  VERDICT: ' + verdict);

const output = {
  phase: 105,
  title: 'Death Match: Sparse Model Competition',
  n: galaxies.length,
  models: results.map(r => ({
    name: r.name, k: r.k, r2: parseFloat(r.r2.toFixed(4)),
    loo: parseFloat(r.loo.toFixed(4)), cv5: parseFloat(r.cv5.toFixed(4)),
    cv10: parseFloat(r.cv10.toFixed(4)), aic: parseFloat(r.aic.toFixed(1)),
    bic: parseFloat(r.bic.toFixed(1))
  })),
  partialCorrelations: partials,
  headToHead: wins,
  correlationMatrix: Object.fromEntries(varNames.map(n1 =>
    [n1, Object.fromEntries(varNames.map(n2 => [n2, parseFloat(pearsonR(varMap[n1], varMap[n2]).toFixed(3))]))])),
  verdict
};

const outPath = path.join(__dirname, '..', 'public', 'phase105-death-match-sparse.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);
